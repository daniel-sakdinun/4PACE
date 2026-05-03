import json
import numpy as np
import networkx as nx
import cvxpy as cp
import warnings

# Suppress CVXPY warnings regarding inaccurate solutions to keep console output clean
warnings.filterwarnings("ignore", message="Solution may be inaccurate.*")

from pandas.core.frame import DataFrame
from fourpace.psys import Grid
from fourpace.model import SynchronousMachine, AsynchronousMachine, Load, Shunt, Transformer


class NumpyEncoder(json.JSONEncoder):
    """
    Custom JSON Encoder to seamlessly serialize NumPy data types and Python complex numbers 
    into standard JSON format for exporting simulation results.
    """
    def default(self, obj):
        if type(obj).__module__ == np.__name__:
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            else:
                return obj.item() 
        
        if isinstance(obj, complex):
            return {"real": obj.real, "imag": obj.imag}
            
        return super(NumpyEncoder, self).default(obj)


def MPOPF(grid: Grid, profile_df: DataFrame = None, relax: str = 'SOCP', solver: str = 'SCS'):
    """
    Multi-Period Optimal Power Flow (MPOPF).
    
    Optimizes the dispatch schedule (Active and Reactive power) for generators, smart inverters, 
    and energy storage systems (BESS) over a given time horizon. Minimizes total generation cost 
    and BESS degradation while strictly adhering to AC power flow physics via convex relaxations.
    
    Parameters:
    -----------
    grid : Grid
        The power system grid object containing topology and components.
    profile_df : DataFrame, optional
        Time-series data for loads and renewable availability. If None, uses the grid's attached profile.
    relax : str
        The convex relaxation method to use for AC power flow ('SOCP' or 'SDP'). Default is 'SOCP'.
    solver : str
        The CVXPY backend solver to use (e.g., 'SCS', 'MOSEK', 'CLARABEL'). Default is 'SCS'.
    """
    active_profile: DataFrame = profile_df if profile_df is not None else grid.load_profile
    if active_profile is None:
        raise ValueError("❌ No load profile provided and no profile attached to Grid!")
        
    T = len(active_profile)
    num_bus = len(grid.buses)
    edges = list(grid.edges(data=True))
    num_branch = len(edges)
    nodes_list = list(grid.nodes)

    # ====================================================================
    # 1. Classify Devices & Map to Buses
    # ====================================================================
    machines, inverters, batteries = [], [], []
    machine_bus, inverter_bus, battery_bus = [], [], []

    for i, bus in enumerate(grid.buses):
        for comp in bus.components:
            c_name = type(comp).__name__
            if c_name == 'SynchronousMachine':
                machines.append(comp)
                machine_bus.append(i)
            elif c_name == 'Inverter':
                inverters.append(comp)
                inverter_bus.append(i)
            elif c_name == 'Battery':
                batteries.append(comp)
                battery_bus.append(i)

    num_gen, num_inv, num_bat = len(machines), len(inverters), len(batteries)

    # ====================================================================
    # 2. Extract Load Profile & Solar Availability Parameters
    # ====================================================================
    P_load_pu = np.zeros((num_bus, T))
    Q_load_pu = np.zeros((num_bus, T))
    Q_shunt_pu = np.zeros((num_bus, T))
    P_solar_avail_pu = np.zeros((num_inv, T)) if num_inv > 0 else None

    for t in range(T):
        grid.apply_profile(active_profile.iloc[t].to_dict())
        for i, bus in enumerate(grid.buses):
            for comp in bus.components:
                c_name = type(comp).__name__
                if c_name == 'Load':
                    P_load_pu[i, t] += abs(comp.P) / grid.Sbase
                    Q_load_pu[i, t] += abs(comp.Q) / grid.Sbase
                elif c_name == 'AsynchronousMachine':
                    comp.update_pq_from_slip(1.0, grid.Sbase)
                    P_load_pu[i, t] += abs(comp.P) / grid.Sbase
                    Q_load_pu[i, t] += abs(comp.Q) / grid.Sbase
                elif c_name == 'Shunt':
                    Q_shunt_pu[i, t] += comp.Q_nom / grid.Sbase
                elif c_name == 'Inverter':
                    inv_idx = inverters.index(comp)
                    P_solar_avail_pu[inv_idx, t] = comp.P / grid.Sbase

    grid.build_ybus()
    G_bus = grid.Ybus.real
    B_bus = grid.Ybus.imag

    # ====================================================================
    # 3. SHARED VARIABLES (Cross-temporal links for the optimization horizon)
    # ====================================================================
    Pg = cp.Variable((num_gen, T)) if num_gen > 0 else None
    Qg = cp.Variable((num_gen, T)) if num_gen > 0 else None
    P_inv = cp.Variable((num_inv, T)) if num_inv > 0 else None
    Q_inv = cp.Variable((num_inv, T)) if num_inv > 0 else None
    P_ch = cp.Variable((num_bat, T)) if num_bat > 0 else None
    P_dis = cp.Variable((num_bat, T)) if num_bat > 0 else None
    SoC = cp.Variable((num_bat, T)) if num_bat > 0 else None

    constraints = []
    slack_idx = next(idx for idx, b in enumerate(grid.buses) if b.type == 'Slack')

    # ====================================================================
    # 4. SWITCHABLE RELAXATION CORE (AC Power Flow Constraints)
    # ====================================================================
    if relax.upper() == 'SOCP':
        # --- Second-Order Cone Programming (Branch Flow Model) ---
        v = cp.Variable((num_bus, T))       # Squared voltage magnitude
        P_ij = cp.Variable((num_branch, T)) # Active power flow
        Q_ij = cp.Variable((num_branch, T)) # Reactive power flow
        l_ij = cp.Variable((num_branch, T)) # Squared current magnitude

        for t in range(T):
            for i in range(num_bus):
                constraints.append(v[i, t] == 1.0 if i == slack_idx else v[i, t] >= 0.85**2)
                if i != slack_idx: 
                    constraints.append(v[i, t] <= 1.15**2)

            P_inj = {i: 0.0 for i in range(num_bus)}
            Q_inj = {i: 0.0 for i in range(num_bus)}

            for k, (u, target_v, data) in enumerate(edges):
                i, j = nodes_list.index(u), nodes_list.index(target_v)
                branch = data.get('obj')
                R, X, Z2 = branch.R, branch.X, branch.R**2 + branch.X**2
                tau = branch.tap_ratio if type(branch).__name__ == 'Transformer' else 1.0
                v_i_eff = v[i, t] / (tau**2)

                # Voltage Drop Equation
                constraints.append(v[j, t] == v_i_eff - 2 * (R * P_ij[k, t] + X * Q_ij[k, t]) + Z2 * l_ij[k, t])
                
                # Conic Relaxation for Apparent Power Limit
                constraints.append(cp.SOC(v_i_eff + l_ij[k, t], cp.vstack([2 * P_ij[k, t], 2 * Q_ij[k, t], v_i_eff - l_ij[k, t]])))

                # Thermal Line Limits
                if getattr(branch, 'S_max', None):
                    constraints.append(cp.norm(cp.vstack([P_ij[k, t], Q_ij[k, t]])) <= branch.S_max / grid.Sbase)

                # Branch Power Balance Accumulation
                P_inj[i] += P_ij[k, t]
                Q_inj[i] += Q_ij[k, t]
                P_inj[j] -= (P_ij[k, t] - R * l_ij[k, t])
                Q_inj[j] -= (Q_ij[k, t] - X * l_ij[k, t])

            # Nodal Power Balance Constraints
            for i in range(num_bus):
                p_g = sum([Pg[k, t] for k, b_idx in enumerate(machine_bus) if b_idx == i]) if num_gen > 0 else 0.0
                q_g = sum([Qg[k, t] for k, b_idx in enumerate(machine_bus) if b_idx == i]) if num_gen > 0 else 0.0
                p_i = sum([P_inv[k, t] for k, b_idx in enumerate(inverter_bus) if b_idx == i]) if num_inv > 0 else 0.0
                q_i = sum([Q_inv[k, t] for k, b_idx in enumerate(inverter_bus) if b_idx == i]) if num_inv > 0 else 0.0
                p_d = sum([P_dis[k, t] for k, b_idx in enumerate(battery_bus) if b_idx == i]) if num_bat > 0 else 0.0
                p_c = sum([P_ch[k, t] for k, b_idx in enumerate(battery_bus) if b_idx == i]) if num_bat > 0 else 0.0
                
                constraints.append(p_g + p_i + p_d - p_c - P_load_pu[i, t] == P_inj[i])
                constraints.append(q_g + q_i - Q_load_pu[i, t] + (Q_shunt_pu[i, t] * v[i, t]) == Q_inj[i])

    elif relax.upper() == 'SDP':
        # --- Semidefinite Programming (Bus Injection Model) ---
        W = [cp.Variable((num_bus, num_bus), hermitian=True) for _ in range(T)]
        for t in range(T):
            constraints.append(W[t] >> 0) # W must be Positive Semidefinite
            WR, WI = cp.real(W[t]), cp.imag(W[t])
            
            for i in range(num_bus):
                constraints.append(WR[i, i] == 1.0 if i == slack_idx else WR[i, i] >= 0.85**2)
                if i != slack_idx: 
                    constraints.append(WR[i, i] <= 1.15**2)

                # Nodal Injection Equations derived from Y-bus
                P_calc = G_bus[i, :] @ WR[i, :] + B_bus[i, :] @ WI[i, :]
                Q_calc = G_bus[i, :] @ WI[i, :] - B_bus[i, :] @ WR[i, :]

                p_g = sum([Pg[k, t] for k, b_idx in enumerate(machine_bus) if b_idx == i]) if num_gen > 0 else 0.0
                q_g = sum([Qg[k, t] for k, b_idx in enumerate(machine_bus) if b_idx == i]) if num_gen > 0 else 0.0
                p_i = sum([P_inv[k, t] for k, b_idx in enumerate(inverter_bus) if b_idx == i]) if num_inv > 0 else 0.0
                q_i = sum([Q_inv[k, t] for k, b_idx in enumerate(inverter_bus) if b_idx == i]) if num_inv > 0 else 0.0
                p_d = sum([P_dis[k, t] for k, b_idx in enumerate(battery_bus) if b_idx == i]) if num_bat > 0 else 0.0
                p_c = sum([P_ch[k, t] for k, b_idx in enumerate(battery_bus) if b_idx == i]) if num_bat > 0 else 0.0
                
                constraints.append(P_calc == p_g + p_i + p_d - p_c - P_load_pu[i, t])
                constraints.append(Q_calc == q_g + q_i - Q_load_pu[i, t] + (Q_shunt_pu[i, t] * WR[i, i]))
    
    # ====================================================================
    # 5. COMMON DEVICE CONSTRAINTS (Time-Coupling & Limits)
    # ====================================================================
    for t in range(T):
        
        # 5.1 Generator Limits (Active and Reactive Power Bounds)
        if num_gen > 0:
            for k, m in enumerate(machines):
                constraints.extend([
                    Pg[k, t] >= (m.Pmin if m.Pmin != float('-inf') else 0.0) / grid.Sbase,
                    Pg[k, t] <= (m.Pmax if m.Pmax != float('inf') else 9999.0) / grid.Sbase,
                    Qg[k, t] >= (m.Qmin if m.Qmin != float('-inf') else -9999.0) / grid.Sbase,
                    Qg[k, t] <= (m.Qmax if m.Qmax != float('inf') else 9999.0) / grid.Sbase
                ])

        # 5.2 Smart Inverter Limits (Active Curtailment & Apparent Power Cone)
        if num_inv > 0:
            for k, inv in enumerate(inverters):
                constraints.extend([
                    P_inv[k, t] >= 0,
                    P_inv[k, t] <= P_solar_avail_pu[k, t],  # Allow solar curtailment if needed
                    cp.norm(cp.vstack([P_inv[k, t], Q_inv[k, t]])) <= inv.S_max / grid.Sbase 
                ])

        # 5.3 Battery Limits & Time-Coupling Dynamics (State of Charge)
        if num_bat > 0:
            for k, bat in enumerate(batteries):
                p_max_pu = bat.P_max / grid.Sbase
                e_max_pu = bat.E_max / grid.Sbase
                
                constraints.extend([
                    P_ch[k, t] >= 0, P_ch[k, t] <= p_max_pu,
                    P_dis[k, t] >= 0, P_dis[k, t] <= p_max_pu,
                    SoC[k, t] >= 0.1, SoC[k, t] <= 1.0  # Maintain at least 10% reserve capacity (DoD limit)
                ])

                # SoC Inter-temporal constraint (Energy balance across time bridging t and t-1)
                delta_soc = (P_ch[k, t] * bat.eta - P_dis[k, t] / bat.eta) / e_max_pu
                if t == 0:
                    constraints.append(SoC[k, t] == bat.init_soc + delta_soc)
                else:
                    constraints.append(SoC[k, t] == SoC[k, t-1] + delta_soc)

    # ====================================================================
    # 6. OBJECTIVE FUNCTION & SOLVE
    # ====================================================================
    cost = 0
    for t in range(T):
        
        # 6.1 Generator Operating Cost (Quadratic formulation)
        if num_gen > 0:
            for k, m in enumerate(machines):
                Pg_mw = Pg[k, t] * grid.Sbase
                cost += m.a + (m.b * Pg_mw) + (m.c * cp.square(Pg_mw))
        
        # 6.2 Renewable Incentives & ESS Degradation
        if num_inv > 0:
            # Reward AI slightly for utilizing solar energy (Negative cost acts as incentive)
            for k, inv in enumerate(inverters): 
                cost -= 0.01 * (P_inv[k, t] * grid.Sbase)
                
        if num_bat > 0:
            # Apply minor degradation cost to prevent unnecessary micro-charge/discharge cycles
            for k, bat in enumerate(batteries): 
                cost += 0.1 * (P_ch[k, t] + P_dis[k, t]) * grid.Sbase

    prob = cp.Problem(cp.Minimize(cost), constraints)
    print(f"\n⏳ Solving 24-Hr Master Plan with MPOPF ({relax.upper()} Relaxation)...")
    
    try:
        prob.solve(solver=getattr(cp, solver.upper()), verbose=False)
    except Exception as e:
        print(f"❌ CVXPY Error: {e}")

    # ====================================================================
    # 7. PARSE AND STORE RESULTS
    # ====================================================================
    if prob.status in ["optimal", "optimal_inaccurate"]:
        print(f"✅ MPOPF ({relax.upper()}) Converged! Master Plan 24h Total Cost: ${prob.value:.2f}")
        
        if num_gen > 0:
            for k, m in enumerate(machines):
                m.P_series = np.array(Pg.value[k, :]) * grid.Sbase
                m.Q_series = np.array(Qg.value[k, :]) * grid.Sbase
        if num_inv > 0:
            for k, inv in enumerate(inverters):
                inv.P_series = np.array(P_inv.value[k, :]) * grid.Sbase
                inv.Q_series = np.array(Q_inv.value[k, :]) * grid.Sbase
        if num_bat > 0:
            for k, bat in enumerate(batteries):
                bat.P_ch_series = np.array(P_ch.value[k, :]) * grid.Sbase
                bat.P_dis_series = np.array(P_dis.value[k, :]) * grid.Sbase
                bat.SoC_series = np.array(SoC.value[k, :])
    else:
        raise Exception(f"❌ MPOPF Infeasible! AI could not find a safe 24h plan (Status: {prob.status})")


def NR(grid: Grid, tol: float = 1e-6, max_iter: int = 100):
    """
    Non-Linear AC Power Flow Solver using the Newton-Raphson (NR) Method.
    
    Computes precise steady-state bus voltages (magnitude and phase angle), branch power flows, 
    and system losses. Includes dynamic PV-to-PQ bus switching to enforce generator reactive 
    power limits (Q-limits), and On-Load Tap Changer (OLTC) auto-adjustments.
    
    Parameters:
    -----------
    grid : Grid
        The power system grid object.
    tol : float
        Convergence tolerance for maximum power mismatch (p.u.). Default is 1e-6.
    max_iter : int
        Maximum number of iterations allowed before aborting. Default is 100.
    """
    print("\n🚀 Starting Newton-Raphson Power Flow...")
    
    Ybus = grid.build_ybus()
    num_bus = len(grid.buses)
    
    iteration = 0
    converged = False
    
    while iteration < max_iter:
        # =======================================================
        # 1. Update Voltage Dependence Dynamics
        # =======================================================
        for bus in grid.buses:
            for comp in bus.components:
                # Update induction motor slip based on current voltage estimation
                if isinstance(comp, AsynchronousMachine):
                    grid.update_motor_slip(comp, bus.V)
                    comp.update_pq_from_slip(bus.V, grid.Sbase)
                # Update capacitor bank injection (scales with V^2)
                elif isinstance(comp, Shunt):
                    comp.update_voltage_dependence(bus.V)

        # =======================================================
        # 2. Classify Dynamic Bus Types (Slack, PV, PQ)
        # =======================================================
        slack, pv, pq = [], [], []
        for i, bus in enumerate(grid.buses):
            if bus.type == 'Slack': 
                slack.append(i)
            elif bus.type == 'PV': 
                pv.append(i)
            else: 
                pq.append(i)
        non_slack = pv + pq

        # =======================================================
        # 3. Prepare State Variables (V, theta) and Specified P, Q
        # =======================================================
        V = np.array([b.V for b in grid.buses], dtype=float)
        theta = np.array([b.theta for b in grid.buses], dtype=float)
        
        P_spec = np.zeros(num_bus)
        Q_spec = np.zeros(num_bus)
        
        for i, bus in enumerate(grid.buses):
            for comp in bus.components:
                c_type = type(comp).__name__
                if c_type == 'SynchronousMachine':
                    P_spec[i] += comp.P / grid.Sbase
                    Q_spec[i] += comp.Q / grid.Sbase
                elif c_type == 'Load':
                    P_spec[i] -= abs(comp.P) / grid.Sbase
                    Q_spec[i] -= abs(comp.Q) / grid.Sbase
                elif c_type in ['Inverter', 'Battery']:
                    P_spec[i] += comp.P / grid.Sbase
                    Q_spec[i] += comp.Q / grid.Sbase
                elif c_type == 'AsynchronousMachine':
                    P_spec[i] -= abs(comp.P) / grid.Sbase
                    Q_spec[i] -= abs(comp.Q) / grid.Sbase
                elif c_type == 'Shunt':
                    Q_spec[i] += comp.Q / grid.Sbase

        # =======================================================
        # 4. Calculate Computed P_calc, Q_calc from Y-bus
        # =======================================================
        Vc = V * np.exp(1j * theta)
        I = Ybus @ Vc
        S_calc = Vc * np.conj(I)
        
        P_calc = S_calc.real
        Q_calc = S_calc.imag

        # =======================================================
        # 5. Generator Q-Limit Checks (PV to PQ Bus Switching)
        # =======================================================
        limit_hit_this_iter = False
        for i in pv:
            bus = grid.buses[i]
            Q_load_pu = sum([abs(c.Q) for c in bus.components if isinstance(c, (Load, AsynchronousMachine))]) / grid.Sbase
            Q_gen_req_pu = Q_calc[i] + Q_load_pu

            machines = [c for c in bus.components if isinstance(c, SynchronousMachine)]
            if machines:
                Qmax_total_pu = sum(getattr(m, 'Qmax', 9999.0) for m in machines) / grid.Sbase
                Qmin_total_pu = sum(getattr(m, 'Qmin', -9999.0) for m in machines) / grid.Sbase
                
                if Q_gen_req_pu > Qmax_total_pu:
                    print(f"⚠️ Iteration {iteration}: Bus {bus.name} exceeds Qmax ({Q_gen_req_pu*grid.Sbase:.2f} > {Qmax_total_pu*grid.Sbase:.2f}) -> Converted to PQ Bus!")
                    bus.type = 'PQ'
                    for m in machines: 
                        m.Q = getattr(m, 'Qmax', 9999.0)
                    limit_hit_this_iter = True
                        
                elif Q_gen_req_pu < Qmin_total_pu:
                    print(f"⚠️ Iteration {iteration}: Bus {bus.name} exceeds Qmin ({Q_gen_req_pu*grid.Sbase:.2f} < {Qmin_total_pu*grid.Sbase:.2f}) -> Converted to PQ Bus!")
                    bus.type = 'PQ'
                    for m in machines: 
                        m.Q = getattr(m, 'Qmin', -9999.0)
                    limit_hit_this_iter = True

        # If limits were hit, re-evaluate mismatches in the next iteration before calculating Jacobian
        if limit_hit_this_iter:
            iteration += 1
            continue

        # =======================================================
        # 6. Compute Mismatches and Check Convergence
        # =======================================================
        dP = P_spec - P_calc
        dQ = Q_spec - Q_calc
        
        dP_ns = dP[non_slack]
        dQ_pq = dQ[pq]
        mismatch = np.concatenate((dP_ns, dQ_pq))
        
        max_error = np.max(np.abs(mismatch)) if len(mismatch) > 0 else 0
        if max_error < tol:
            converged = True
            break  # Escaped NR Loop successfully!
            
        # =======================================================
        # 7. Construct and Solve Jacobian Matrix [J]
        # =======================================================
        diagVc = np.diag(Vc)
        diagI_conj = np.diag(np.conj(I))
        diagVc_conj = np.diag(np.conj(Vc))
        diagV_phase = np.diag(Vc / V)
        diagV_phase_conj = np.diag(np.conj(Vc) / V)
        
        dS_dTheta = 1j * diagVc @ diagI_conj - 1j * diagVc @ np.conj(Ybus) @ diagVc_conj
        dS_dV = diagV_phase @ diagI_conj + diagVc @ np.conj(Ybus) @ diagV_phase_conj
        
        H = dS_dTheta.real
        N = dS_dV.real
        M = dS_dTheta.imag
        L = dS_dV.imag
        
        J11 = H[np.ix_(non_slack, non_slack)]
        J12 = N[np.ix_(non_slack, pq)]
        J21 = M[np.ix_(pq, non_slack)]
        J22 = L[np.ix_(pq, pq)]
        
        J_top = np.hstack((J11, J12)) if J12.size else J11
        J_bottom = np.hstack((J21, J22)) if J12.size else J21
        J = np.vstack((J_top, J_bottom))
        
        try:
            dx = np.linalg.solve(J, mismatch)
        except np.linalg.LinAlgError:
            print("❌ Singular Jacobian Matrix! The grid might be physically collapsed (Voltage Collapse).")
            break
            
        dTheta = dx[:len(non_slack)]
        dV = dx[len(non_slack):]
        
        for idx, i in enumerate(non_slack): 
            grid.buses[i].theta += dTheta[idx]
        for idx, i in enumerate(pq): 
            grid.buses[i].V += dV[idx]

        # =======================================================
        # 8. On-Load Tap Changer (OLTC) Auto-Tap Updates
        # =======================================================
        tap_changed = False
        for u, v, data in grid.edges(data=True):
            branch = data.get('obj')
            if isinstance(branch, Transformer) and getattr(branch, 'auto_tap', False) and getattr(branch, 'controlled_bus', None):
                target_bus = grid.bus(branch.controlled_bus)
                deadband = 0.015  # Voltage tolerance band before stepping tap
                
                if target_bus.V < branch.target_V - deadband and branch.tap_ratio > branch.tap_min:
                    branch.tap_ratio -= branch.tap_step
                    branch.tap_ratio = max(branch.tap_ratio, branch.tap_min)
                    tap_changed = True
                    print(f"🔄 Iteration {iteration}: {branch.name} Step DOWN Tap -> {branch.tap_ratio:.4f}")
                    
                elif target_bus.V > branch.target_V + deadband and branch.tap_ratio < branch.tap_max:
                    branch.tap_ratio += branch.tap_step
                    branch.tap_ratio = min(branch.tap_ratio, branch.tap_max)
                    tap_changed = True
                    print(f"🔄 Iteration {iteration}: {branch.name} Step UP Tap -> {branch.tap_ratio:.4f}")

        # If tap changes, physical Y-bus matrix is altered and must be rebuilt
        if tap_changed:
            Ybus = grid.build_ybus()

        iteration += 1
        
    # =======================================================
    # 9. Post-Processing (Calculate Slack & PV actual generation)
    # =======================================================
    V_final = np.array([b.V for b in grid.buses])
    theta_final = np.array([b.theta for b in grid.buses])
    Vc_final = V_final * np.exp(1j * theta_final)
    S_final = Vc_final * np.conj(Ybus @ Vc_final)
    
    for i, bus in enumerate(grid.buses):
        P_load_total = sum([abs(c.P) / grid.Sbase for c in bus.components if type(c).__name__ in ['Load', 'AsynchronousMachine']])
        Q_load_total = sum([abs(c.Q) / grid.Sbase for c in bus.components if type(c).__name__ in ['Load', 'AsynchronousMachine']])
        
        P_inv_total = sum([c.P / grid.Sbase for c in bus.components if type(c).__name__ in ['Inverter', 'Battery']])
        Q_inv_total = sum([c.Q / grid.Sbase for c in bus.components if type(c).__name__ in ['Inverter', 'Battery']])
        
        Q_shunt_total = sum([c.Q / grid.Sbase for c in bus.components if type(c).__name__ == 'Shunt'])
        
        if bus.type == 'Slack':
            P_slack_gen = S_final.real[i] + P_load_total - P_inv_total
            Q_slack_gen = S_final.imag[i] + Q_load_total - Q_inv_total - Q_shunt_total
            for comp in bus.components:
                if type(comp).__name__ == 'SynchronousMachine':
                    comp.P = float(P_slack_gen * grid.Sbase)
                    comp.Q = float(Q_slack_gen * grid.Sbase)
                    
        elif bus.type == 'PV':
            Q_pv_gen = S_final.imag[i] + Q_load_total - Q_inv_total - Q_shunt_total
            for comp in bus.components:
                if type(comp).__name__ == 'SynchronousMachine':
                    comp.Q = float(Q_pv_gen * grid.Sbase)

    if converged:
        print(f"✅ Newton-Raphson Converged seamlessly in {iteration+1} iterations!")
    else:
        print(f"⚠️ Newton-Raphson Failed after {max_iter} iterations (Max Error: {max_error:.6f})")
        raise Exception(f"Voltage Collapse! (NR Max Error: {max_error:.2f})")


def plan(grid: Grid, profile_df: DataFrame = None, path: str = "settings.json", relax: str = 'SOCP', solver: str = 'SCS', tol: float = 1e-6, max_iter: int = 100):
    """
    Executes a complete Time-Series Simulation and Master Plan.
    
    This function wraps the MPOPF (to find the optimal economic dispatch schedules) and the 
    Newton-Raphson solver (to validate the exact AC physics at each time step). The final results, 
    including voltages, angles, and line flows, are exported to a JSON file.
    
    Parameters:
    -----------
    grid : Grid
        The power system grid object.
    profile_df : DataFrame, optional
        Time-series data for the planning horizon.
    path : str
        File path to export the resulting settings and states. Default is "settings.json".
    relax : str
        Relaxation method for the MPOPF ('SOCP' or 'SDP').
    solver : str
        Optimization backend solver.
    tol : float
        Tolerance for the NR solver.
    max_iter : int
        Maximum iterations for the NR solver.
        
    Returns:
    --------
    str
        The JSON string containing the comprehensive simulation history.
    """
    active_profile = profile_df if profile_df is not None else grid.load_profile
    T = len(active_profile)
    
    print(f"\n🔮 Starting Modern Grid Simulation (MPOPF) for {T} steps...")
    
    # 1. Generate optimal generation schedules using OPF
    MPOPF(grid, profile_df=active_profile, relax=relax, solver=solver)

    history = []
        
    # 2. Replay the schedules through the precise Non-Linear NR engine to get exact voltages
    for t in range(T):
        print(f"\n{'='*15} 🕒 Step {t} {'='*15}")
        grid.apply_profile(active_profile.iloc[t].to_dict())
            
        # Apply OPF schedules to hardware components
        for bus in grid.buses:
            for comp in bus.components:
                c_name = type(comp).__name__
                if c_name == 'SynchronousMachine':
                    comp.P = float(comp.P_series[t])
                    comp.Q = float(comp.Q_series[t])
                elif c_name == 'Inverter':
                    comp.P = float(comp.P_series[t])
                    comp.Q = float(comp.Q_series[t])
                elif c_name == 'Battery':
                    comp.P = float(comp.P_dis_series[t] - comp.P_ch_series[t])  # Net active power injection
                    comp.SoC = float(comp.SoC_series[t])
            
        # Validate exact AC Physics for the step
        NR(grid, max_iter=max_iter, tol=tol)
        total_cost = sum([bus.total_cost() for bus in grid.buses])
            
        # Compile Step Dashboard
        step_result = {
            'step': t,
            'system': {'total_cost_per_hr': round(total_cost, 2)},
            'buses': {}, 'components': {}, 'branches': {}
        }
            
        for bus in grid.buses:
            step_result['buses'][bus.name] = {
                'V_pu': round(bus.V, 4), 'theta_rad': round(bus.theta, 4),
                'P_total_MW': round(bus.P * grid.Sbase, 2), 'Q_total_MVAr': round(bus.Q * grid.Sbase, 2)
            }
            for comp in bus.components:
                comp_data = {
                    'bus': bus.name, 'type': type(comp).__name__,
                    'P_MW': round(comp.P * grid.Sbase, 2) if hasattr(comp, 'P') else 0.0,
                    'Q_MVAr': round(comp.Q * grid.Sbase, 2) if hasattr(comp, 'Q') else 0.0
                }
                if type(comp).__name__ == 'Battery':
                    comp_data['SoC'] = round(comp.SoC, 4)
                step_result['components'][comp.name] = comp_data
            
        for u, v, data in grid.edges(data=True):
            branch = data.get('obj')
            if branch:
                from_bus, to_bus = grid.bus(u), grid.bus(v)
                V_i = from_bus.V * np.exp(1j * from_bus.theta)
                V_j = to_bus.V * np.exp(1j * to_bus.theta)
                
                # Flow calculation depending on Tap Ratio (Transformer) or strict impedance (Line)
                if type(branch).__name__ == 'Transformer':
                    t_tap = branch.tap_ratio * np.exp(1j * branch.phase_shift)
                    I_ij = (V_i / t_tap - V_j) * branch.Y
                else:
                    I_ij = (V_i - V_j) * branch.Y
                    
                S_flow_MVA = abs(V_i * np.conj(I_ij)) * grid.Sbase
                branch_data = {
                    'from_bus': u, 'to_bus': v, 'type': type(branch).__name__,
                    'flow_MVA': round(S_flow_MVA, 2)
                }
                if getattr(branch, 'S_max', None):
                    loading_pct = (S_flow_MVA / branch.S_max) * 100
                    branch_data['loading_percent'] = round(loading_pct, 2)
                    branch_data['is_overload'] = loading_pct > 100
                step_result['branches'][branch.name] = branch_data
                    
        history.append(step_result)
    
    grid.loading_status()
    
    # Export cleanly formatted JSON Data
    result = json.dumps(history, indent=2, cls=NumpyEncoder)
    with open(path, "w") as f:
        f.write(result)
            
    print("\n✅ MPOPF Time-Series Simulation Completed!")
    return result


def CEP(grid: Grid, profile_df: DataFrame = None, relax: str = 'SOCP', solver: str = 'SCS'):
    """
    Capacity Expansion Planning (CEP) Engine.
    
    Co-optimizes long-term investment sizing (Capital Expenditure - CapEx) against short-term 
    operational schedules (Operational Expenditure - OpEx). Determines the absolute most cost-effective 
    capacity (MW / MWh) of candidate Renewable Energy and BESS assets required to satisfy load profiles 
    while respecting all AC grid physics constraints.
    
    Parameters:
    -----------
    grid : Grid
        The power system grid object containing designated 'candidate' components.
    profile_df : DataFrame, optional
        Time-series load profile.
    relax : str
        Relaxation formulation ('SOCP' or 'SDP').
    solver : str
        Optimization solver.
    """
    print("\n🏗️ Starting Capacity Expansion Planning (CEP)...")
    
    active_profile: DataFrame = profile_df if profile_df is not None else grid.load_profile
    if active_profile is None:
        raise ValueError("❌ No load profile provided and no profile attached to Grid!")
    
    T = len(active_profile)
    num_bus = len(grid.buses)
    edges = list(grid.edges(data=True))
    num_branch = len(edges)
    nodes_list = list(grid.nodes)

    machines, inverters, batteries = [], [], []
    machine_bus, inverter_bus, battery_bus = [], [], []

    for i, bus in enumerate(grid.buses):
        for comp in bus.components:
            c_name = type(comp).__name__
            if c_name == 'SynchronousMachine': 
                machines.append(comp); machine_bus.append(i)
            elif c_name == 'Inverter': 
                inverters.append(comp); inverter_bus.append(i)
            elif c_name == 'Battery': 
                batteries.append(comp); battery_bus.append(i)

    num_gen, num_inv, num_bat = len(machines), len(inverters), len(batteries)

    P_load_pu = np.zeros((num_bus, T))
    Q_load_pu = np.zeros((num_bus, T))
    Q_shunt_pu = np.zeros((num_bus, T))
    P_solar_avail_pu = np.zeros((num_inv, T)) if num_inv > 0 else None

    for t in range(T):
        grid.apply_profile(active_profile.iloc[t].to_dict())
        for i, bus in enumerate(grid.buses):
            for comp in bus.components:
                c_name = type(comp).__name__
                if c_name == 'Load':
                    P_load_pu[i, t] += abs(comp.P) / grid.Sbase
                    Q_load_pu[i, t] += abs(comp.Q) / grid.Sbase
                elif c_name == 'AsynchronousMachine':
                    comp.update_pq_from_slip(1.0, grid.Sbase)
                    P_load_pu[i, t] += abs(comp.P) / grid.Sbase
                    Q_load_pu[i, t] += abs(comp.Q) / grid.Sbase
                elif c_name == 'Shunt':
                    Q_shunt_pu[i, t] += comp.Q_nom / grid.Sbase
                elif c_name == 'Inverter':
                    inv_idx = inverters.index(comp)
                    # For candidate inverters, extract the raw multiplier directly from the profile
                    if getattr(comp, 'is_candidate', False):
                        multiplier = active_profile.iloc[t].to_dict().get(comp.name, 0.0)
                        P_solar_avail_pu[inv_idx, t] = multiplier
                    else:
                        P_solar_avail_pu[inv_idx, t] = comp.P / grid.Sbase

    grid.build_ybus()
    G_bus = grid.Ybus.real
    B_bus = grid.Ybus.imag
    
    Pg = cp.Variable((num_gen, T)) if num_gen > 0 else None
    Qg = cp.Variable((num_gen, T)) if num_gen > 0 else None
    P_inv = cp.Variable((num_inv, T)) if num_inv > 0 else None
    Q_inv = cp.Variable((num_inv, T)) if num_inv > 0 else None
    P_ch = cp.Variable((num_bat, T)) if num_bat > 0 else None
    P_dis = cp.Variable((num_bat, T)) if num_bat > 0 else None
    
    # Use E_stored (MWh equivalent) instead of SoC to avoid nonlinear division in optimization constraints
    E_stored = cp.Variable((num_bat, T)) if num_bat > 0 else None

    constraints = []
    slack_idx = next(idx for idx, b in enumerate(grid.buses) if b.type == 'Slack')

    # Apply Standard Power Flow Relaxations (Similar to MPOPF)
    if relax.upper() == 'SOCP':
        v = cp.Variable((num_bus, T))
        P_ij = cp.Variable((num_branch, T))
        Q_ij = cp.Variable((num_branch, T))
        l_ij = cp.Variable((num_branch, T))

        for t in range(T):
            for i in range(num_bus):
                constraints.append(v[i, t] == 1.0 if i == slack_idx else v[i, t] >= 0.85**2)
                if i != slack_idx: constraints.append(v[i, t] <= 1.15**2)

            P_inj = {i: 0.0 for i in range(num_bus)}
            Q_inj = {i: 0.0 for i in range(num_bus)}

            for k, (u, target_v, data) in enumerate(edges):
                i, j = nodes_list.index(u), nodes_list.index(target_v)
                branch = data.get('obj')
                R, X, Z2 = branch.R, branch.X, branch.R**2 + branch.X**2
                tau = branch.tap_ratio if type(branch).__name__ == 'Transformer' else 1.0
                v_i_eff = v[i, t] / (tau**2)

                constraints.append(v[j, t] == v_i_eff - 2 * (R * P_ij[k, t] + X * Q_ij[k, t]) + Z2 * l_ij[k, t])
                constraints.append(cp.SOC(v_i_eff + l_ij[k, t], cp.vstack([2 * P_ij[k, t], 2 * Q_ij[k, t], v_i_eff - l_ij[k, t]])))

                if getattr(branch, 'S_max', None):
                    constraints.append(cp.norm(cp.vstack([P_ij[k, t], Q_ij[k, t]])) <= branch.S_max / grid.Sbase)

                P_inj[i] += P_ij[k, t]; Q_inj[i] += Q_ij[k, t]
                P_inj[j] -= (P_ij[k, t] - R * l_ij[k, t]); Q_inj[j] -= (Q_ij[k, t] - X * l_ij[k, t])

            for i in range(num_bus):
                p_g = sum([Pg[k, t] for k, b_idx in enumerate(machine_bus) if b_idx == i]) if num_gen > 0 else 0.0
                q_g = sum([Qg[k, t] for k, b_idx in enumerate(machine_bus) if b_idx == i]) if num_gen > 0 else 0.0
                p_i = sum([P_inv[k, t] for k, b_idx in enumerate(inverter_bus) if b_idx == i]) if num_inv > 0 else 0.0
                q_i = sum([Q_inv[k, t] for k, b_idx in enumerate(inverter_bus) if b_idx == i]) if num_inv > 0 else 0.0
                p_d = sum([P_dis[k, t] for k, b_idx in enumerate(battery_bus) if b_idx == i]) if num_bat > 0 else 0.0
                p_c = sum([P_ch[k, t] for k, b_idx in enumerate(battery_bus) if b_idx == i]) if num_bat > 0 else 0.0
                constraints.append(p_g + p_i + p_d - p_c - P_load_pu[i, t] == P_inj[i])
                constraints.append(q_g + q_i - Q_load_pu[i, t] + (Q_shunt_pu[i, t] * v[i, t]) == Q_inj[i])

    elif relax.upper() == 'SDP':
        W = [cp.Variable((num_bus, num_bus), hermitian=True) for _ in range(T)]
        for t in range(T):
            constraints.append(W[t] >> 0)
            WR, WI = cp.real(W[t]), cp.imag(W[t])
            for i in range(num_bus):
                constraints.append(WR[i, i] == 1.0 if i == slack_idx else WR[i, i] >= 0.85**2)
                if i != slack_idx: constraints.append(WR[i, i] <= 1.15**2)

                P_calc = G_bus[i, :] @ WR[i, :] + B_bus[i, :] @ WI[i, :]
                Q_calc = G_bus[i, :] @ WI[i, :] - B_bus[i, :] @ WR[i, :]

                p_g = sum([Pg[k, t] for k, b_idx in enumerate(machine_bus) if b_idx == i]) if num_gen > 0 else 0.0
                q_g = sum([Qg[k, t] for k, b_idx in enumerate(machine_bus) if b_idx == i]) if num_gen > 0 else 0.0
                p_i = sum([P_inv[k, t] for k, b_idx in enumerate(inverter_bus) if b_idx == i]) if num_inv > 0 else 0.0
                q_i = sum([Q_inv[k, t] for k, b_idx in enumerate(inverter_bus) if b_idx == i]) if num_inv > 0 else 0.0
                p_d = sum([P_dis[k, t] for k, b_idx in enumerate(battery_bus) if b_idx == i]) if num_bat > 0 else 0.0
                p_c = sum([P_ch[k, t] for k, b_idx in enumerate(battery_bus) if b_idx == i]) if num_bat > 0 else 0.0
                constraints.append(P_calc == p_g + p_i + p_d - p_c - P_load_pu[i, t])
                constraints.append(Q_calc == q_g + q_i - Q_load_pu[i, t] + (Q_shunt_pu[i, t] * WR[i, i]))

    # ====================================================================
    # 5. CEP SPECIFIC: INVESTMENT VARIABLES & CONSTRAINTS (CAPEX formulation)
    # ====================================================================
    total_capex = 0.0
    bat_p_max_capacity = []
    bat_e_max_capacity = []
    inv_s_max_capacity = [] 

    # --- Consider Solar Investments ---
    if num_inv > 0:
        for k, inv in enumerate(inverters):
            if getattr(inv, 'is_candidate', False):
                built_s_mw = cp.Variable(nonneg=True) # AI decides how many MW to build
                constraints.append(built_s_mw <= inv.max_build_mw)
                
                inv_s_max_capacity.append(built_s_mw / grid.Sbase)
                total_capex += (built_s_mw * inv.capex_per_mw) * getattr(inv, 'daily_capex_factor', 1.0)
            else:
                inv_s_max_capacity.append(inv.S_max / grid.Sbase)

    # --- Consider Battery Investments ---
    if num_bat > 0:
        for k, bat in enumerate(batteries):
            if getattr(bat, 'is_candidate', False):
                built_p_mw = cp.Variable(nonneg=True)  # Power Capacity Investment
                built_e_mwh = cp.Variable(nonneg=True) # Energy Storage Investment
                
                constraints.extend([
                    built_p_mw <= bat.max_build_mw,
                    built_e_mwh <= bat.max_build_mwh
                ])
                
                bat_p_max_capacity.append(built_p_mw / grid.Sbase)
                bat_e_max_capacity.append(built_e_mwh / grid.Sbase)
                total_capex += (built_p_mw * bat.capex_per_mw + built_e_mwh * bat.capex_per_mwh) * getattr(bat, 'daily_capex_factor', 1.0)
            else:
                bat_p_max_capacity.append(bat.P_max / grid.Sbase)
                bat_e_max_capacity.append(bat.E_max / grid.Sbase)

    # ====================================================================
    # 6. TIME-COUPLING WITH DYNAMIC CAPACITY
    # ====================================================================
    for t in range(T):
        if num_gen > 0:
            for k, m in enumerate(machines):
                constraints.extend([
                    Pg[k, t] >= (m.Pmin if m.Pmin != float('-inf') else 0.0) / grid.Sbase,
                    Pg[k, t] <= (m.Pmax if m.Pmax != float('inf') else 9999.0) / grid.Sbase,
                    Qg[k, t] >= (m.Qmin if m.Qmin != float('-inf') else -9999.0) / grid.Sbase,
                    Qg[k, t] <= (m.Qmax if m.Qmax != float('inf') else 9999.0) / grid.Sbase
                ])

        if num_inv > 0:
            for k, inv in enumerate(inverters):
                # Available generation bounds scale linearly with built capacity
                if getattr(inv, 'is_candidate', False):
                    p_avail = P_solar_avail_pu[k, t] * inv_s_max_capacity[k]
                else:
                    p_avail = P_solar_avail_pu[k, t]

                constraints.extend([
                    P_inv[k, t] >= 0,
                    P_inv[k, t] <= p_avail,
                    cp.norm(cp.vstack([P_inv[k, t], Q_inv[k, t]])) <= inv_s_max_capacity[k]
                ])

        if num_bat > 0:
            for k, bat in enumerate(batteries):
                p_max_pu = bat_p_max_capacity[k]
                e_max_pu = bat_e_max_capacity[k]
                
                constraints.extend([
                    P_ch[k, t] >= 0, P_ch[k, t] <= p_max_pu,
                    P_dis[k, t] >= 0, P_dis[k, t] <= p_max_pu,
                    E_stored[k, t] >= 0.1 * e_max_pu,
                    E_stored[k, t] <= 1.0 * e_max_pu
                ])

                delta_e = (P_ch[k, t] * bat.eta - P_dis[k, t] / bat.eta)
                if t == 0:
                    init_e = bat.init_soc * e_max_pu
                    constraints.append(E_stored[k, t] == init_e + delta_e)
                else:
                    constraints.append(E_stored[k, t] == E_stored[k, t-1] + delta_e)

    # ====================================================================
    # 7. OBJECTIVE: MINIMIZE (CAPEX + OPEX)
    # ====================================================================
    total_opex = 0
    for t in range(T):
        if num_gen > 0:
            for k, m in enumerate(machines):
                Pg_mw = Pg[k, t] * grid.Sbase
                total_opex += m.a + (m.b * Pg_mw) + (m.c * cp.square(Pg_mw))
    
    total_cost = total_capex + total_opex

    prob = cp.Problem(cp.Minimize(total_cost), constraints)
    print(f"⏳ Optimizing Sizing and Operation ({relax.upper()})...")
    
    try:
        prob.solve(solver=getattr(cp, solver.upper()), verbose=False)
    except Exception as e:
        print(f"❌ CVXPY Error: {e}")

    # ====================================================================
    # 8. POST-PROCESS AND ASSIGN BUILT CAPACITIES
    # ====================================================================
    if prob.status in ["optimal", "optimal_inaccurate"]:
        print(f"✅ CEP Converged! Total Cost (CapEx + OpEx): ${prob.value:.2f}")
        
        if num_inv > 0:
            for k, inv in enumerate(inverters):
                if getattr(inv, 'is_candidate', False):
                    inv.built_S_max = float(inv_s_max_capacity[k].value * grid.Sbase)
                    inv.S_max = inv.built_S_max
                    print(f"   ☀️ Investment Decision -> {inv.name}: S_max = {inv.built_S_max:.2f} MW")

        if num_bat > 0:
            for k, bat in enumerate(batteries):
                if getattr(bat, 'is_candidate', False):
                    bat.built_P_max = float(bat_p_max_capacity[k].value * grid.Sbase)
                    bat.built_E_max = float(bat_e_max_capacity[k].value * grid.Sbase)
                    bat.P_max = bat.built_P_max
                    bat.E_max = bat.built_E_max
                    print(f"   🔋 Investment Decision -> {bat.name}: P_max = {bat.built_P_max:.2f} MW, E_max = {bat.built_E_max:.2f} MWh")
    else:
        raise Exception(f"❌ CEP Infeasible! Status: {prob.status}")


def N1_Screening(grid, peak_hour: int = None, relax: str = 'SOCP', solver: str = 'SCS'):
    """
    N-1 Contingency Screening utilizing Optimal Power Flow.
    
    Systematically disconnects one transmission branch at a time to test grid resilience 
    during peak loading. Identifies weak physical links that cause system collapse or unsolvable 
    overloads.
    
    Parameters:
    -----------
    grid : Grid
        The power system grid object.
    peak_hour : int, optional
        Specific hour index to test. If None, auto-detects the hour with maximum active power load.
    relax : str
        Relaxation method ('SOCP' or 'SDP').
    solver : str
        Solver backend (e.g., 'SCS').
    """
    if peak_hour is None:
        peak_hour = grid.get_peak_load_hour()
    
    print(f"\n🌪️ Starting N-1 Contingency Screening (Snapshot: Hour {peak_hour})...")
    single_hour_df = grid.load_profile.iloc[[peak_hour]].reset_index(drop=True)
    
    # ========================================================
    # 1. Verify standard system health (N-0) before any outages
    # ========================================================
    print("🏥 Testing Base Case (N-0) health...")
    try:
        MPOPF(grid, profile_df=single_hour_df, relax=relax, solver=solver)
        print("✅ Base Case is healthy! Proceeding to cut lines...\n")
    except Exception as e:
        print(f"\n💀 FATAL: The system is already collapsing BEFORE cutting any lines!")
        print(f"🚨 Actual Error: {e}")
        print("💡 Hint: Verify Battery E_max is non-zero, or check for undervoltage limits below 0.85 p.u.")
        return 
    
    # ========================================================
    # 2. Systematically disconnect lines to test N-1 scenarios
    # ========================================================
    original_edges = list(grid.edges(data=True))
    num_branches = len(original_edges)
    print(f"🔍 Testing {num_branches} branch contingencies...")
    report = {}
    
    for u, v, data in original_edges:
        branch = data['obj']
        branch_name = branch.name
        grid.remove_edge(u, v)
        
        try:
            MPOPF(grid, single_hour_df, relax=relax, solver=solver)
            report[branch_name] = "✅ Survived"
            
        except Exception as e:
            report[branch_name] = f"❌ CRITICAL: {str(e)}"
            
        finally:
            grid.add_edge(u, v, **data) # Always restore the branch before next iteration
            
    print("\n========================================")
    print("📊 N-1 Contingency Report:")
    print("========================================")
    for name, status in report.items():
        print(f"Outage of {name:<8} -> {status}")


def SCOPF(grid, peak_hour: int = None, relax: str = 'SOCP', solver: str = 'SCS'):
    """
    Security-Constrained Optimal Power Flow (SCOPF) with ESS Rescue.
    
    Formulates a massive "Multiverse" optimization problem that simultaneously calculates 
    the base case (N-0) and all single-branch outages (N-1). Ensures that synchronous 
    generators have a single, fixed setpoint across all universes (preventive), while allowing 
    fast-acting Inverters and Batteries to react dynamically (corrective) to survive outages.
    
    Parameters:
    -----------
    grid : Grid
        The power system grid object.
    peak_hour : int, optional
        Hour index representing peak loading. Auto-detected if None.
    relax : str
        Relaxation formulation ('SOCP' or 'SDP').
    solver : str
        Solver backend (e.g., 'CLARABEL', 'SCS').
        
    Returns:
    --------
    dict
        The 'Rescue Plan' mapping each contingency to optimal Battery/Inverter setpoints.
    """
    if peak_hour is None:
        peak_hour = grid.get_peak_load_hour()
        
    print(f"\n🛡️ Starting Preventive SCOPF (Multiverse: Hour {peak_hour}) with ESS Rescue...")
    
    # 1. Inject 1-hour peak profile into the grid snapshot
    single_hour_df = grid.load_profile.iloc[[peak_hour]].reset_index(drop=True)
    original_profile = grid.load_profile 
    grid.attach_profile(single_hour_df)  
    
    grid.apply_profile(0)
    
    edges = list(grid.edges(data=True))
    nodes_list = list(grid.nodes)
    num_bus = len(grid.buses)
    num_branch = len(edges)
    K = num_branch + 1 # Total universes: Base Case (k=0) + N-1 Outages (k=1 to K-1)

    # 2. Group all active power devices
    machines, inverters, batteries = [], [], []
    machine_bus, inverter_bus, battery_bus = [], [], []
    
    for i, bus in enumerate(grid.buses):
        for comp in bus.components:
            c_type = type(comp).__name__
            if c_type == 'SynchronousMachine': machines.append(comp); machine_bus.append(i)
            elif c_type == 'Inverter': inverters.append(comp); inverter_bus.append(i)
            elif c_type == 'Battery': batteries.append(comp); battery_bus.append(i)
                
    num_gen, num_inv, num_bat = len(machines), len(inverters), len(batteries)
    slack_idx = next(idx for idx, b in enumerate(grid.buses) if b.type == 'Slack')

    # ====================================================================
    # 3. MULTIVERSE VARIABLES (Simulating K networks in parallel)
    # ====================================================================
    Pg = cp.Variable((num_gen, K)) if num_gen > 0 else None
    Qg = cp.Variable((num_gen, K)) if num_gen > 0 else None
    P_inv = cp.Variable((num_inv, K)) if num_inv > 0 else None
    Q_inv = cp.Variable((num_inv, K)) if num_inv > 0 else None
    P_ch = cp.Variable((num_bat, K)) if num_bat > 0 else None
    P_dis = cp.Variable((num_bat, K)) if num_bat > 0 else None
    
    v = cp.Variable((num_bus, K))
    P_ij = cp.Variable((num_branch, K))
    Q_ij = cp.Variable((num_branch, K))
    l_ij = cp.Variable((num_branch, K))

    constraints = []

    for k in range(K):
        # 3.1 Voltage Limits per universe
        for i in range(num_bus):
            constraints.append(v[i, k] == 1.0 if i == slack_idx else v[i, k] >= 0.85**2)
            if i != slack_idx: constraints.append(v[i, k] <= 1.15**2)

        P_inj = {i: 0.0 for i in range(num_bus)}; Q_inj = {i: 0.0 for i in range(num_bus)}

        # 3.2 Branch Flow Formulation (SOCP)
        for b_idx, (u, target_v, data) in enumerate(edges):
            # If we are in universe k > 0, we slice the branch corresponding to (k-1)
            if k > 0 and b_idx == (k - 1):
                constraints.extend([P_ij[b_idx, k] == 0, Q_ij[b_idx, k] == 0, l_ij[b_idx, k] == 0])
                continue 

            i, j = nodes_list.index(u), nodes_list.index(target_v)
            branch = data.get('obj')
            R, X, Z2 = branch.R, branch.X, branch.R**2 + branch.X**2
            tau = branch.tap_ratio if type(branch).__name__ == 'Transformer' else 1.0
            v_i_eff = v[i, k] / (tau**2)

            constraints.append(v[j, k] == v_i_eff - 2*(R*P_ij[b_idx, k] + X*Q_ij[b_idx, k]) + Z2*l_ij[b_idx, k])
            constraints.append(cp.SOC(v_i_eff + l_ij[b_idx, k], cp.vstack([2*P_ij[b_idx, k], 2*Q_ij[b_idx, k], v_i_eff - l_ij[b_idx, k]])))

            if getattr(branch, 'S_max', None):
                constraints.append(cp.norm(cp.vstack([P_ij[b_idx, k], Q_ij[b_idx, k]])) <= branch.S_max / grid.Sbase)

            P_inj[i] += P_ij[b_idx, k]; Q_inj[i] += Q_ij[b_idx, k]
            P_inj[j] -= (P_ij[b_idx, k] - R*l_ij[b_idx, k]); Q_inj[j] -= (Q_ij[b_idx, k] - X*l_ij[b_idx, k])

        # 3.3 Nodal Power Balance per universe
        for i in range(num_bus):
            bus = grid.buses[i]
            P_L = sum([abs(c.P)/grid.Sbase for c in bus.components if isinstance(c, Load)])
            Q_L = sum([abs(c.Q)/grid.Sbase for c in bus.components if isinstance(c, Load)])
            Q_sh = sum([c.Q_nom/grid.Sbase for c in bus.components if isinstance(c, Shunt)])

            p_g = sum([Pg[m, k] for m, b_idx in enumerate(machine_bus) if b_idx == i]) if num_gen > 0 else 0.0
            q_g = sum([Qg[m, k] for m, b_idx in enumerate(machine_bus) if b_idx == i]) if num_gen > 0 else 0.0
            p_i = sum([P_inv[m, k] for m, b_idx in enumerate(inverter_bus) if b_idx == i]) if num_inv > 0 else 0.0
            q_i = sum([Q_inv[m, k] for m, b_idx in enumerate(inverter_bus) if b_idx == i]) if num_inv > 0 else 0.0
            p_d = sum([P_dis[m, k] for m, b_idx in enumerate(battery_bus) if b_idx == i]) if num_bat > 0 else 0.0
            p_c = sum([P_ch[m, k] for m, b_idx in enumerate(battery_bus) if b_idx == i]) if num_bat > 0 else 0.0

            constraints.append(p_g + p_i + p_d - p_c - P_L == P_inj[i])
            constraints.append(q_g + q_i - Q_L + (Q_sh * v[i, k]) == Q_inj[i])

        # 3.4 Synchronous Generators (Locked to Base Case)
        # Generators cannot react instantly to line drops; their P setpoint must remain identical across all scenarios.
        if num_gen > 0:
            for m, gen in enumerate(machines):
                constraints.extend([
                    Pg[m, k] >= (gen.Pmin / grid.Sbase), Pg[m, k] <= (gen.Pmax / grid.Sbase),
                    Qg[m, k] >= (gen.Qmin / grid.Sbase), Qg[m, k] <= (gen.Qmax / grid.Sbase)
                ])
                if k > 0 and machine_bus[m] != slack_idx:
                    constraints.append(Pg[m, k] == Pg[m, 0]) # 🔗 Force lock to Base Case

        # 3.5 Solar Inverters (Flexible Real-Time Curtailment)
        if num_inv > 0:
            for m, inv in enumerate(inverters):
                solar_avail = single_hour_df.iloc[0].to_dict().get(inv.name, inv.P / grid.Sbase)
                constraints.extend([
                    P_inv[m, k] >= 0, P_inv[m, k] <= solar_avail,
                    cp.norm(cp.vstack([P_inv[m, k], Q_inv[m, k]])) <= inv.S_max / grid.Sbase
                ])

        # 3.6 Batteries (Corrective Rescue Action)
        # Fast-acting ESS can instantly alter dispatch to rescue the grid during N-1 contingencies
        if num_bat > 0:
            for m, bat in enumerate(batteries):
                p_max_pu = bat.P_max / grid.Sbase
                if p_max_pu > 0:
                    constraints.extend([
                        P_ch[m, k] >= 0, P_ch[m, k] <= p_max_pu,
                        P_dis[m, k] >= 0, P_dis[m, k] <= p_max_pu
                    ])
                else:
                    constraints.extend([P_ch[m, k] == 0, P_dis[m, k] == 0])

    # ====================================================================
    # 4. OBJECTIVE: Minimize Generation Cost of the Base Case (k=0) Only
    # ====================================================================
    cost = 0
    if num_gen > 0:
        for m, gen in enumerate(machines):
            Pg_mw = Pg[m, 0] * grid.Sbase
            cost += gen.a + (gen.b * Pg_mw) + (gen.c * cp.square(Pg_mw))

    prob = cp.Problem(cp.Minimize(cost), constraints)
    print(f"⏳ Solving {K} Simultaneous Networks ({relax.upper()})...")
    
    try:
        prob.solve(solver=getattr(cp, solver.upper()), verbose=False)
    except Exception as e:
        print(f"❌ CVXPY Error: {e}")

    # ====================================================================
    # 5. PARSE SCOPF RESULTS AND GENERATE RESCUE PLAN
    # ====================================================================
    if prob.status in ["optimal", "optimal_inaccurate"]:
        print(f"✅ SCOPF Converged! Preventive Security Total Cost: ${prob.value:.2f}")
        
        # Restore safe P, Q setpoints into Hardware Objects
        if num_gen > 0:
            for m, gen in enumerate(machines):
                gen.P = float(Pg.value[m, 0] * grid.Sbase)
                gen.Q = float(Qg.value[m, 0] * grid.Sbase)
                print(f"   🏭 {gen.name:<5} Setpoint -> P = {gen.P:7.2f} MW (Safe for all N-1!)")
        
        if num_bat > 0:
             print("   🔋 Battery Base-Case Dispatch:")
             for m, bat in enumerate(batteries):
                 bat.P = float(P_dis.value[m, 0] * grid.Sbase - P_ch.value[m, 0] * grid.Sbase)
                 
                 if abs(bat.P) < 1e-4:
                     state = "Idle"
                 elif bat.P > 0:
                     state = "Discharging"
                 else:
                     state = "Charging"
                     
                 print(f"      {bat.name:<5} -> {abs(bat.P):5.2f} MW ({state})")

        if num_inv > 0:
            print("   ☀️ Inverter Base-Case Dispatch:")
            for m, inv in enumerate(inverters):
                inv.P = float(P_inv.value[m, 0] * grid.Sbase)
                inv.Q = float(Q_inv.value[m, 0] * grid.Sbase)
                print(f"      {inv.name:<5} -> P = {inv.P:7.2f} MW | Q = {inv.Q:7.2f} MVAr")

        rescue_plan = {}
        # Extract N-1 scenario emergency actions (k=1 to K-1)
        for k in range(1, K): 
            b_idx = k - 1
            branch_name = edges[b_idx][2]['obj'].name
            rescue_plan[branch_name] = {}
            
            if num_bat > 0:
                for m, bat in enumerate(batteries):
                    p_net = float(P_dis.value[m, k] * grid.Sbase - P_ch.value[m, k] * grid.Sbase)
                    rescue_plan[branch_name][bat.name] = {'P': p_net}
            
            if num_inv > 0:
                for m, inv in enumerate(inverters):
                    p_inv_val = float(P_inv.value[m, k] * grid.Sbase)
                    q_inv_val = float(Q_inv.value[m, k] * grid.Sbase)
                    rescue_plan[branch_name][inv.name] = {'P': p_inv_val, 'Q': q_inv_val}
                    
        grid.attach_profile(original_profile)
        return rescue_plan
    else:
        print(f"❌ SCOPF Infeasible! AI cannot find a safe state (Status: {prob.status})")
        
    # Restore original time-series profile on failure
    grid.attach_profile(original_profile)


def Validate_N1(grid: Grid, rescue_plan: dict = None, tol: float = 1e-6, max_iter: int = 100):
    """
    Physical Verification Auditor for N-1 Constraints.
    
    Acts as a strict physics validation engine. It drops branches sequentially, applies the 
    AI-generated emergency rescue plans (from SCOPF), and runs the strict non-linear Newton-Raphson 
    AC load flow to ensure that line flows and voltage deviations are practically sound and do not 
    trigger subsequent overloads.
    
    Parameters:
    -----------
    grid : Grid
        The power system grid object to validate.
    rescue_plan : dict, optional
        A dictionary mapping branch outages to specific Battery/Inverter corrective actions.
    tol : float
        Tolerance for the NR solver.
    max_iter : int
        Maximum iterations for the NR solver.
    """
    print("\n🔎 ========================================")
    print("🛡️ Starting N-1 Physics Validation (NR)...")
    print("========================================")
    
    report = {}
    edges = list(grid.edges(data=True))
    
    # Verify N-0 baseline first
    try:
        NR(grid, tol=tol, max_iter=max_iter)
        status = grid.check_overload()
        if any(pct > 100 for pct in status['branches'].values()):
            print("❌ Base Case FAILED: Overload detected in N-0!")
            return
        print("✅ Base Case (N-0) Passed Physics Test.")
    except Exception as e:
        print(f"❌ Base Case FAILED: NR did not converge! ({e})")
        return
    
    # 1.5 Create a system state checkpoint (Record the pristine N-0 state)
    base_case_state = {}
    for bus in grid.buses:
        b_data = {'V': bus.V, 'theta': bus.theta, 'type': bus.type, 'comps': {}}
        for comp in bus.components:
            if type(comp).__name__ == 'SynchronousMachine':
                b_data['comps'][comp.name] = {'Q': comp.Q}
            elif type(comp).__name__ in ['Battery', 'Inverter']:
                b_data['comps'][comp.name] = {'P': comp.P, 'Q': comp.Q}
        base_case_state[bus.name] = b_data
        
    base_case_taps = {}
    for u, v, data in edges:
        branch = data['obj']
        if type(branch).__name__ == 'Transformer':
            base_case_taps[branch.name] = branch.tap_ratio

    # 2. Begin cutting lines iteratively (N-1 testing)
    for u, v, data in edges:
        branch = data['obj']
        b_name = branch.name
        
        # ✂️ Simulate physical branch severing
        grid.remove_edge(u, v)
        
        # Check for islanding (disconnected graph topology preventing convergence)
        if not nx.is_connected(grid):
            report[b_name] = "⏩ SKIPPED (Radial Line / Islanding)"
            grid.add_edge(u, v, **data)
            continue
        
        # 2.1 Apply emergency rescue plan (Activate Batteries and Inverters based on SCOPF logic)
        if rescue_plan and b_name in rescue_plan:
            for bus in grid.buses:
                for comp in bus.components:
                    if type(comp).__name__ == 'Battery' and comp.name in rescue_plan[b_name]:
                        plan_data = rescue_plan[b_name][comp.name]
                        comp.P = plan_data['P'] if isinstance(plan_data, dict) else plan_data
                        
                    elif type(comp).__name__ == 'Inverter' and comp.name in rescue_plan[b_name]:
                        plan_data = rescue_plan[b_name][comp.name]
                        if isinstance(plan_data, dict):
                            comp.P = plan_data.get('P', comp.P)
                            comp.Q = plan_data.get('Q', comp.Q)
        
        grid.build_ybus()  # Rebuild Admittance Matrix without the cut line
        
        # 🧪 Physics Test via exact AC Newton-Raphson Load Flow
        try:
            NR(grid, tol=tol, max_iter=max_iter) 
            status = grid.check_overload()
            overloaded = [name for name, pct in status['branches'].items() if pct > 100]
            
            if overloaded:
                report[b_name] = f"❌ FAILED: Overloaded {overloaded}"
            else:
                report[b_name] = "✅ PASSED (No Overload)"
                
        except Exception as e:
            report[b_name] = f"❌ FAILED: {type(e).__name__} - {e}"
            
        finally:
            # 🩹 Reconnect the severed line
            grid.add_edge(u, v, **data)
            
            # 2.2 Load Save Game! (Restore everything completely to the pristine baseline N-0 state)
            for bus in grid.buses:
                b_data = base_case_state[bus.name]
                bus.V = b_data['V']
                bus.theta = b_data['theta']
                bus.type = b_data['type']
                
                for comp in bus.components:
                    if comp.name in b_data['comps']:
                        c_data = b_data['comps'][comp.name]
                        if 'P' in c_data: comp.P = c_data['P']
                        if 'Q' in c_data: comp.Q = c_data['Q']
                        
            for u_edge, v_edge, t_data in edges:
                t_branch = t_data['obj']
                if type(t_branch).__name__ == 'Transformer' and t_branch.name in base_case_taps:
                    t_branch.tap_ratio = base_case_taps[t_branch.name]
            
            grid.build_ybus()  # Reconstruct Ybus from the intact state
            
    # Final Reporting Output
    print("\n📊 Physics Validation Report:")
    for name, res in report.items():
        print(f"   Drop {name:<8} -> {res}")