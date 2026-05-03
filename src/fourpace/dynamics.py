import numpy as np
import pandas as pd

from fourpace.psys import Grid
from fourpace.model import SynchronousMachine, TransmissionLine
from fourpace.facts import ShuntFACTS, SeriesFACTS

# =====================================================================
# 1. State Mapping & Dynamic Allocation
# =====================================================================
def get_state_indices(machines, shunt_facts, series_facts):
    """
    Dynamically maps and allocates contiguous memory arrays for every generator, controller (AVR, GOV, PSS), 
    and FACTS device in the system. This allows the RK4 solver to process all states efficiently in a single 
    flattened vector.
    
    Parameters
    -----------
    machines : list
        List of tuples containing (bus_index, SynchronousMachine_object).
    shunt_facts : list
        List of tuples containing (bus_index, ShuntFACTS_object).
    series_facts : list
        List of tuples containing (branch_dict, SeriesFACTS_object).
        
    Returns
    --------
    tuple
        (indices_dict, total_number_of_states)
    """
    indices = {'machines': {}, 'shunt_facts': {}, 'series_facts': {}}
    current_idx = 0
    
    # 1. Map Generator States (Core + Controllers)
    for i, (bus_idx, gen) in enumerate(machines):
        start_idx = current_idx
        core_len = 3 # Basic machine states: delta, omega, Eq'
        avr_len = gen.avr.n_states if gen.avr else 0
        gov_len = gen.gov.n_states if gen.gov else 0
        pss_len = gen.pss.n_states if gen.pss else 0
        total_len = core_len + avr_len + gov_len + pss_len
        
        indices['machines'][i] = {
            'core': (start_idx, start_idx + 3),
            'avr':  (start_idx + 3, start_idx + 3 + avr_len) if avr_len > 0 else None,
            'gov':  (start_idx + 3 + avr_len, start_idx + 3 + avr_len + gov_len) if gov_len > 0 else None,
            'pss':  (start_idx + 3 + avr_len + gov_len, start_idx + total_len) if pss_len > 0 else None,
        }
        current_idx += total_len
        
    # 2. Map Shunt FACTS States (e.g., SVC, STATCOM)
    for i, (bus_idx, fact) in enumerate(shunt_facts):
        length = fact.n_states
        indices['shunt_facts'][i] = (current_idx, current_idx + length)
        current_idx += length
        
    # 3. Map Series FACTS States (e.g., TCSC)
    for i, (branch, fact) in enumerate(series_facts):
        length = fact.n_states
        indices['series_facts'][i] = (current_idx, current_idx + length)
        current_idx += length
        
    return indices, current_idx

# =====================================================================
# 2. The Universal Algebraic ODE Engine
# =====================================================================
def ode_engine(t: float, state: np.ndarray, Y_base: np.ndarray, 
                   machines: list, shunt_facts: list, series_facts: list, 
                   indices: dict, omega_s: float,
                   V_ref_dict: dict, P_ref_dict: dict, E_fd_0_dict: dict, 
                   algebraic_vars: dict) -> np.ndarray:
    """
    The core Differential-Algebraic Equation (DAE) solver engine. 
    It evaluates the derivatives (dx/dt) of all state variables at time 't'.
    
    This function utilizes a Norton Equivalent Network approach: 
    1. Generator voltage sources are converted to Norton current injections.
    2. Shunt and Series FACTS devices dynamically modify the Admittance Matrix (Y_eval).
    3. The algebraic network equation [I] = [Y][V] is solved to find bus voltages.
    4. These new voltages drive the differential equations of the machines and controllers.
    """
    n_bus = Y_base.shape[0]
    d_state = np.zeros_like(state)
    
    Y_eval = Y_base.copy()
    I_inj = np.zeros(n_bus, dtype=complex)
    last_V_bus = algebraic_vars['V_bus']
    
    # --- ASSEMBLE THE DYNAMIC NETWORK (ALGEBRAIC PART) ---
    
    # A. Inject Series FACTS (e.g., TCSC) directly into Y-bus
    for i, (branch, tcsc) in enumerate(series_facts):
        idx = indices['series_facts'][i]
        local_state = state[idx[0]:idx[1]]
        X_tcsc = tcsc.get_X_series(local_state) # Active reactance compensation
        
        frm, to = branch['from_idx'], branch['to_idx']
        base_R, base_X = branch['R'], branch['X']
        
        # 🛡️ SAFETY CLAMP: Prevent Series Resonance Explosion (Limit to 80% compensation)
        max_cap_limit = -0.8 * base_X
        if X_tcsc < max_cap_limit:
            X_tcsc = max_cap_limit
            
        y_old = 1 / (base_R + 1j * base_X)
        y_new = 1 / (base_R + 1j * (base_X + X_tcsc))
        delta_y = y_new - y_old
        
        Y_eval[frm, frm] += delta_y
        Y_eval[to, to] += delta_y
        Y_eval[frm, to] -= delta_y
        Y_eval[to, frm] -= delta_y
        branch['y_new'] = y_new # Store for power calculation later

    # B. Inject Shunt FACTS (e.g., SVC, STATCOM) into Y-bus diagonals
    for i, (bus_idx, fact) in enumerate(shunt_facts):
        idx = indices['shunt_facts'][i]
        local_state = state[idx[0]:idx[1]]
        
        if type(fact).__name__ == 'CSVGN1': # Static Var Compensator (SVC)
            B_svc = fact.get_susceptance(local_state)
            Y_eval[bus_idx, bus_idx] += 1j * B_svc
            
        elif type(fact).__name__ == 'STATCOM1': # STATCOM
            I_q = fact.get_Iq(local_state)
            # 🛡️ ALGEBRAIC LOOP FIX: Treat the current injection as a dynamic voltage-dependent admittance
            # Admittance Y = I / V = -j * (I_q / |V|)
            v_mag = max(abs(last_V_bus[bus_idx]), 0.05) # Prevent division by zero during faults
            Y_eval[bus_idx, bus_idx] += -1j * (I_q / v_mag)

    # C. Inject Generators (Norton Equivalent Current Sources)
    E_internal_list = []
    for i, (bus_idx, gen) in enumerate(machines):
        core_idx = indices['machines'][i]['core']
        delta, omega, Eq_prime = state[core_idx[0]], state[core_idx[0]+1], state[core_idx[0]+2]
        
        # Generator internal voltage phasor
        E_internal = Eq_prime * np.exp(1j * delta)
        E_internal_list.append(E_internal)
        
        # Generator internal admittance
        y_gen = 1 / (1j * gen.Xd_prime)
        
        # Add to network
        Y_eval[bus_idx, bus_idx] += y_gen
        I_inj[bus_idx] += E_internal * y_gen

    # --- SOLVE THE NETWORK ALGEBRAIC EQUATIONS ---
    # Find all bus voltages using [V] = [Y]^-1 * [I]
    V_bus = np.linalg.solve(Y_eval, I_inj)
    algebraic_vars['V_bus'] = V_bus 
    
    # --- EVALUATE DIFFERENTIAL EQUATIONS (dx/dt) ---
    
    # A. Series FACTS Derivatives
    for i, (branch, tcsc) in enumerate(series_facts):
        idx = indices['series_facts'][i]
        local_state = state[idx[0]:idx[1]]
        frm, to = branch['from_idx'], branch['to_idx']
        
        # 🛡️ TRANSIENT BYPASS LOGIC: If grid voltage collapses (e.g., severe fault nearby),
        # the TCSC's Metal-Oxide Varistor (MOV) conducts to bypass and protect the series capacitor.
        if abs(V_bus[frm]) < 0.8 or abs(V_bus[to]) < 0.8:
            d_state[idx[0]:idx[1]] = 0.0  # Freeze the internal PI controller integral states
        else:
            # Calculate actual line power to drive the TCSC Power Oscillation Damping (POD) loop
            I_line = (V_bus[frm] - V_bus[to]) * branch['y_new']
            P_line = (V_bus[frm] * np.conj(I_line)).real
            d_state[idx[0]:idx[1]] = tcsc.get_derivatives(local_state, P_line)
        
    # B. Shunt FACTS Derivatives
    for i, (bus_idx, fact) in enumerate(shunt_facts):
        idx = indices['shunt_facts'][i]
        local_state = state[idx[0]:idx[1]]
        V_mag = abs(V_bus[bus_idx]) # Voltage magnitude drives the SVC/STATCOM voltage regulator
        d_state[idx[0]:idx[1]] = fact.get_derivatives(local_state, V_mag)

    # C. Generator & Controller Derivatives (Swing Equation & Flux Decay)
    for i, (bus_idx, gen) in enumerate(machines):
        idx = indices['machines'][i]
        core_idx = idx['core']
        delta, omega, Eq_prime = state[core_idx[0]], state[core_idx[0]+1], state[core_idx[0]+2]
        E_internal = E_internal_list[i]
        
        y_gen = 1 / (1j * gen.Xd_prime)
        I_gen = (E_internal - V_bus[bus_idx]) * y_gen
        P_e = (E_internal * np.conj(I_gen)).real # Electrical Power Output
        V_t_mag = abs(V_bus[bus_idx])            # Terminal Voltage Magnitude
        
        # Transform current to d-q synchronous reference frame
        angle_I, angle_d = np.angle(I_gen), delta - (np.pi / 2)
        I_d = abs(I_gen) * np.cos(angle_I - angle_d) # d-axis current
        
        omega_pu = omega / omega_s
        V_pss = 0.0
        
        # 1. Power System Stabilizer (PSS) 
        if gen.pss:
            local_pss = state[idx['pss'][0]:idx['pss'][1]]
            V_pss = gen.pss.get_Vpss(local_pss, omega_pu, P_e)
            d_state[idx['pss'][0]:idx['pss'][1]] = gen.pss.get_derivatives(local_pss, omega_pu, P_e)

        E_fd = E_fd_0_dict[i] # Default field voltage if no AVR
        
        # 2. Automatic Voltage Regulator (AVR)
        if gen.avr:
            local_avr = state[idx['avr'][0]:idx['avr'][1]]
            E_fd = gen.avr.get_Efd(local_avr)
            d_state[idx['avr'][0]:idx['avr'][1]] = gen.avr.get_derivatives(local_avr, V_ref_dict[i], V_t_mag, V_pss)

        P_m = P_ref_dict[i] # Default mechanical power if no Governor
        
        # 3. Turbine Governor (GOV)
        if gen.gov:
            local_gov = state[idx['gov'][0]:idx['gov'][1]]
            P_m = gen.gov.get_Pm(local_gov)
            d_state[idx['gov'][0]:idx['gov'][1]] = gen.gov.get_derivatives(local_gov, omega_pu, P_ref_dict[i])

        # 4. Core Generator Physics (Swing Equation & Flux dynamics)
        d_state[core_idx[0]] = omega # d(delta)/dt = omega
        d_state[core_idx[0]+1] = (omega_s / (2 * gen.H)) * (P_m - P_e) # d(omega)/dt
        d_state[core_idx[0]+2] = (E_fd - Eq_prime - I_d * (gen.Xd - gen.Xd_prime)) / gen.Td0_prime # d(Eq')/dt

    return d_state

def rk4_step(state: np.ndarray, t: float, dt: float, *args) -> np.ndarray:
    """Standard 4th-Order Runge-Kutta numerical integration step."""
    k1 = ode_engine(t, state, *args)
    k2 = ode_engine(t + dt/2, state + (dt/2)*k1, *args)
    k3 = ode_engine(t + dt/2, state + (dt/2)*k2, *args)
    k4 = ode_engine(t + dt, state + dt*k3, *args)
    return state + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

# =====================================================================
# 3. Dynamic Orchestrator
# =====================================================================
def analyze_transient(grid: Grid, fault_bus: str, t_clear: float, t_end: float = 5.0, dt: float = 0.01, verbose: bool = True) -> pd.DataFrame:
    """
    Executes a time-domain transient stability simulation.
    
    A 3-Phase fault is applied at `fault_bus` at t=0 and is cleared at `t_clear`.
    Tracks the rotor angle (delta) of all synchronous machines to determine if synchronism is maintained.
    
    Parameters
    -----------
    grid : Grid
        The initialized and solved (steady-state load flow) power system grid.
    fault_bus : str
        Name of the bus to apply the solid 3-phase fault.
    t_clear : float
        Time in seconds when the fault is cleared.
    t_end : float
        Total duration of the simulation in seconds.
    dt : float
        Integration time step (s).
        
    Returns
    --------
    pd.DataFrame
        Time-series data containing simulation time and rotor angles.
    """
    if verbose: print(f"\n🎢 INITIATING DYNAMIC SIMULATION (FACTS Norton Engine)")
    
    omega_s = 2 * np.pi * 50.0
    nodes_list = list(grid.nodes)
    
    # 1. Harvest & Map Devices from the Network
    machines, shunt_facts = [], []
    for bus in grid.buses:
        bus_idx = nodes_list.index(bus.name)
        for comp in bus.components:
            if isinstance(comp, SynchronousMachine): machines.append((bus_idx, comp))
            elif isinstance(comp, ShuntFACTS): shunt_facts.append((bus_idx, comp))
            
    series_facts = []
    for fact in grid.series_facts:
        for u, v, data in grid.edges(data=True):
            if data.get('obj') and data['obj'].name == fact.branch_name:
                branch = {'from_idx': nodes_list.index(u), 'to_idx': nodes_list.index(v), 
                          'R': data['obj'].R, 'X': data['obj'].X}
                series_facts.append((branch, fact))
                break

    indices, total_states = get_state_indices(machines, shunt_facts, series_facts)
    initial_state = np.zeros(total_states)
    
    # 2. Base Y-Bus Generation
    # We borrow the highly-accurate Steady-State Ybus (which includes Transformer Taps and Line shunts)
    Y_bare = grid.build_ybus().copy() 
    
    # Convert constant power loads to constant admittance (Z-model) for dynamic stability evaluation
    for bus in grid.buses:
        i = nodes_list.index(bus.name)
        for comp in bus.components:
            if type(comp).__name__ in ['Load', 'Shunt', 'AsynchronousMachine']:
                S_load = (comp.P + 1j * comp.Q) / grid.Sbase
                # Y = P - jQ / |V|^2
                Y_bare[i, i] += np.conj(S_load) / (bus.V**2)

    # Create the faulted Y-bus (solid short to ground at fault_bus)
    Y_bare_fault = Y_bare.copy()
    fault_idx = nodes_list.index(fault_bus)
    Y_bare_fault[fault_idx, fault_idx] += 1e6 
    
    # 3. Initialization
    V_ref_dict, P_ref_dict, E_fd_0_dict = {}, {}, {}
    algebraic_vars = {'V_bus': np.array([bus.V * np.exp(1j * bus.theta) for bus in grid.buses])}
    
    # Initialize Generators to match the initial load flow state
    for i, (bus_idx, gen) in enumerate(machines):
        bus = grid.bus(nodes_list[bus_idx])
        P_pu, Q_pu = gen.P / grid.Sbase, gen.Q / grid.Sbase
        V_term = bus.V * np.exp(1j * bus.theta)
        
        V_ref_dict[i], P_ref_dict[i] = bus.V, P_pu
        S_gen = P_pu + 1j * Q_pu
        I_gen = np.conj(S_gen / V_term)
        
        E_int = V_term + (1j * gen.Xd_prime) * I_gen
        delta_0, Eq_prime_0 = np.angle(E_int), abs(E_int)
        
        angle_I, angle_d = np.angle(I_gen), delta_0 - (np.pi / 2)
        I_d_init = abs(I_gen) * np.cos(angle_I - angle_d)
        E_fd_0_dict[i] = Eq_prime_0 + I_d_init * (gen.Xd - gen.Xd_prime)
        
        idx = indices['machines'][i]
        initial_state[idx['core'][0]:idx['core'][0]+3] = [delta_0, 0.0, Eq_prime_0]
        
        # Initialize sub-controllers
        if gen.avr: initial_state[idx['avr'][0]:idx['avr'][1]] = gen.avr.initialize(V_ref_dict[i], bus.V, E_fd_0_dict[i])
        if gen.gov: initial_state[idx['gov'][0]:idx['gov'][1]] = gen.gov.initialize(P_ref_dict[i], 0.0)
        if gen.pss: initial_state[idx['pss'][0]:idx['pss'][1]] = gen.pss.initialize()

    # Initialize Shunt FACTS
    for i, (bus_idx, fact) in enumerate(shunt_facts):
        idx = indices['shunt_facts'][i]
        initial_state[idx[0]:idx[1]] = fact.initialize(grid.bus(nodes_list[bus_idx]).V)
        
    # 🛡️ STEADY-STATE FIX: Auto-Initialize Series FACTS to match exact Load Flow
    for i, (branch, fact) in enumerate(series_facts):
        idx = indices['series_facts'][i]
        
        frm, to = branch['from_idx'], branch['to_idx']
        V_frm = algebraic_vars['V_bus'][frm]
        V_to = algebraic_vars['V_bus'][to]
        y_line = 1 / (branch['R'] + 1j * branch['X'])
        
        I_line = (V_frm - V_to) * y_line
        P_initial = (V_frm * np.conj(I_line)).real
        
        # Override the user's P_ref with the actual steady-state load flow measurement
        if hasattr(fact, 'P_ref'):
            fact.P_ref = P_initial
            
        initial_state[idx[0]:idx[1]] = fact.initialize()

    # 4. Numerical Integration Loop
    time_steps = np.arange(0, t_end, dt)
    history = {'Time_s': time_steps}
    for i, (_, gen) in enumerate(machines): history[f"{gen.name}_Delta_deg"] = []
    
    current_state = initial_state.copy()
    
    for t in time_steps:
        # Record rotor angles
        for i, (_, gen) in enumerate(machines):
            history[f"{gen.name}_Delta_deg"].append(np.rad2deg(current_state[indices['machines'][i]['core'][0]]))
            
        # Switch Y-bus based on fault status
        Y_current = Y_bare_fault if (0 <= t < t_clear) else Y_bare
        
        # Execute RK4 step
        args = (Y_current, machines, shunt_facts, series_facts, indices, omega_s, V_ref_dict, P_ref_dict, E_fd_0_dict, algebraic_vars)
        current_state = rk4_step(current_state, t, dt, *args)
        
    if verbose: print("✅ FACTS Transient Simulation Complete!")
    return pd.DataFrame(history)

# =====================================================================
# 4. Critical Clearing Time (CCT) Finder
# =====================================================================
def find_cct(grid: Grid, fault_bus: str, t_min: float = 0.05, t_max: float = 0.50, tol: float = 0.005, path: str | None = None) -> dict:
    """
    Automated Binary Search Algorithm to locate the Critical Clearing Time (CCT).
    
    The CCT is the maximum time a fault can remain on the system before a generator 
    loses synchronism (rotor angle deviation > 180 degrees).
    
    Parameters
    -----------
    grid : Grid
        The power system grid object.
    fault_bus : str
        The bus to test the 3-phase fault on.
    t_min : float
        Lower bound of the search window (seconds). Default 0.05.
    t_max : float
        Upper bound of the search window (seconds). Default 0.50.
    tol : float
        Search precision tolerance. Search stops when high-low < tol. Default 0.005.
    path : str, optional
        Filepath to save the CSV trajectory of the final stable run.
        
    Returns
    --------
    dict
        A dictionary containing the final results of the binary search:
        - 'CCT_s' (float): The calculated Critical Clearing Time.
        - 'Simulation' (pandas.DataFrame or None): The time-series trajectory of the most critical stable scenario.
    """
    print(f"\n🔍 INITIATING CCT BINARY SEARCH (Norton FACTS Engine)")
    
    low, high, cct, best_df = t_min, t_max, 0.0, None
    
    while (high - low) > tol:
        t_test = (low + high) / 2.0
        print(f"   ➤ Testing t_clear = {t_test:.4f} s...", end="")
        try:
            # Simulate a 3-second horizon to check for stability
            df = analyze_transient(grid, fault_bus, t_clear=t_test, t_end=3.0, verbose=False) 
        except Exception as e:
            print(f" ❌ ERROR: {e}"); break
            
        # Stability Criterion: Max rotor angle difference among all machines must not exceed 180 degrees
        is_stable, angle_cols = True, [c for c in df.columns if 'Delta_deg' in c]
        for idx, row in df.iterrows():
            angles = row[angle_cols].values
            if (np.max(angles) - np.min(angles)) > 180.0:
                is_stable = False; break
                
        # Update Binary Search bounds
        if is_stable: 
            print(" ✅ STABLE")
            low, cct, best_df = t_test, t_test, df
        else: 
            print(" 💥 UNSTABLE")
            high = t_test

    print(f"\n🎯 BINARY SEARCH COMPLETE | 🏆 Ultimate CCT: {cct:.3f} seconds")
    
    if path and best_df is not None: 
        best_df.to_csv(path, index=False)
    
    return {'CCT_s': cct, 'Simulation': best_df}