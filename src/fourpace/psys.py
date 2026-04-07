import yaml
from pathlib import Path
import numpy as np
import scipy.optimize as opt
import networkx as nx

from fourpace.model import BusComponent, SynchronousMachine, AsynchronousMachine, Load, Shunt, Inverter, BranchComponent, TransmissionLine, Transformer

class Bus(nx.Graph):
    def __init__(self, name: str, Vbase: float, type: str = 'PQ'):
        super().__init__()
        self.name = name
        self.Vbase: float | None = Vbase
        
        self.type = type
        self.V:float = 1.0
        self.theta:float = 0.0
        self.add_node(name, obj=self)
    
    @property
    def P(self) -> float:
        total_p = 0.0
        for comp in self.components:
            if hasattr(comp, 'P'):
                total_p += comp.P
        return total_p

    @property
    def Q(self) -> float:
        total_q = 0.0
        for comp in self.components:
            if hasattr(comp, 'Q'):
                total_q += comp.Q
        return total_q
    
    @property
    def S(self) -> complex:
        return self.P + 1j*self.Q;
    
    def component(self, name: str) -> BusComponent:
        if name in self.nodes:
            return self.nodes[name]['obj']
        
        raise KeyError(f"❌ '{name}' not found in this bus.")
    
    @property
    def components(self) -> list:
        comp_list = []
        for _, data in self.nodes(data=True):
            obj = data.get('obj')
            if obj is not None and obj is not self:
                comp_list.append(obj)
        return comp_list
    
    def get(self):
        return [self.V, self.theta, self.P, self.Q]
    
    def add_component(self, component: BusComponent):
        self.add_node(component.name, obj=component)
        self.add_edge(self.name, component.name)
    
    def add_components(self, components: list[BusComponent]):
        for component in components:
            self.add_component(component)
        
    def total_cost(self) -> float:
        components = [component for _, component in self.nodes(data='obj')]
        components.remove(self)
        sum = 0
        for com in components:
            sum += com.cost()
        return sum

class Grid(nx.Graph):
    def __init__(self, Sbase: float):
        super().__init__()
        self.Ybus: np.ndarray | None = None
        self.Sbase: float = Sbase
    
    @classmethod
    def load(cls, filepath: str) -> 'Grid':
        path = Path(filepath)
        with open(path, 'r', encoding='utf-8') as f:
            if path.suffix in ['.yaml', '.yml']:
                data = yaml.safe_load(f)
            else:
                raise ValueError("❌ Unsupported format! Please use .yaml or .json")

        grid = cls(Sbase=data.get('Sbase', 100.0))

        comp_classes = {
            'SynchronousMachine': SynchronousMachine,
            'AsynchronousMachine': AsynchronousMachine,
            'Load': Load,
            'Shunt': Shunt,
            'Inverter': Inverter,
            'TransmissionLine': TransmissionLine,
            'Transformer': Transformer
        }

        for b_data in data.get('buses', []):
            bus = Bus(name=b_data['name'], Vbase=b_data['Vbase'], type=b_data.get('bus_type', 'PQ'))
            grid.add_bus(bus)
            
            for comp_data in b_data.get('components', []):
                comp_type = comp_data.pop('type')
                if comp_type in comp_classes:
                    comp_obj = comp_classes[comp_type](**comp_data) 
                    grid.bus(bus.name).add_component(comp_obj)

        for branch_data in data.get('branches', []):
            branch_type = branch_data.pop('type')
            from_bus = branch_data.pop('from_bus')
            to_bus = branch_data.pop('to_bus')
            
            if branch_type in comp_classes:
                branch_obj = comp_classes[branch_type](**branch_data)
                grid.connect(from_bus, to_bus, branch_obj)

        print(f"✅ Successfully loaded grid from {path.name}")
        return grid
    
    def add_bus(self, bus:Bus):
        self.add_node(bus.name, bus=bus)
    
    def add_busses(self, busses:list[Bus]):
        for bus in busses:
            self.add_bus(bus)
        
    def connect(self, from_bus:str, to_bus:str, branch: BranchComponent):
        self.add_edge(from_bus, to_bus, obj=branch)
    
    def bus(self, name: str) -> Bus:
        if name in self.nodes:
            return self.nodes[name]['bus']
        
        raise KeyError(f"❌  Bus  '{name}' not found in this grid.")
    
    @property
    def buses(self) -> list[Bus]:
        return [bus for _, bus in self.nodes(data='bus')]
    
    def build_ybus(self) -> np.ndarray:
        nodes = list(self.nodes)
        n = len(nodes)
        Y = np.zeros((n, n), dtype=complex)
        
        for u, v, data in self.edges(data=True):
            i = nodes.index(u)
            j = nodes.index(v)
            obj = data.get('obj')
            
            if obj is None: continue
            
            y = obj.Y
            
            if isinstance(obj, TransmissionLine):
                b_sh = obj.B_shunt
                Y[i, i] += y + (1j * b_sh / 2)
                Y[j, j] += y + (1j * b_sh / 2)
                Y[i, j] -= y
                Y[j, i] -= y
                
            elif isinstance(obj, Transformer):
                a = obj.tap_ratio
                alpha = obj.phase_shift
                
                tap_complex = a * np.exp(1j * alpha)
                
                Y[i, i] += y / (a**2)
                Y[j, j] += y
                Y[i, j] -= y / np.conj(tap_complex)
                Y[j, i] -= y / tap_complex
                
        self.Ybus:np.ndarray = Y
        return Y
    
    def result(self):
        print("\n=== 📊 POWER FLOW RESULTS ===")
        for i, bus in enumerate(self.buses):
            P_pu, Q_pu = self.calculate_PQ(i)
            
            P_actual:float = P_pu * self.Sbase
            Q_actual:float = Q_pu * self.Sbase
            
            print(f"Bus {bus.name} | V = {bus.V:.4f} pu | phase = {np.rad2deg(bus.theta):.2f}° | P = {P_actual:7.2f} MW | Q = {Q_actual:7.2f} MVAr")
    
    def calculate_PQ(self, i) -> tuple[float, float]:
        bus_i:Bus = self.buses[i]
        Vi = bus_i.V
        theta_i = bus_i.theta

        G:np.ndarray = self.Ybus.real
        B:np.ndarray = self.Ybus.imag

        P_i:float = 0.0
        Q_i:float = 0.0

        for j in range(len(self.buses)):
            bus_j:Bus = self.buses[j]
            Vj:float = bus_j.V
            theta_j:float = bus_j.theta
            delta_ij:float = theta_i - theta_j

            P_i += Vi * Vj * (G[i, j] * np.cos(delta_ij) + B[i, j] * np.sin(delta_ij))
            Q_i += Vi * Vj * (G[i, j] * np.sin(delta_ij) - B[i, j] * np.cos(delta_ij))
        
        return P_i, Q_i
    
    def update_motor_slip(self, motor, V_pu):
        def torque_mismatch(s):
            if s <= 0: return 1e6
            
            V_phase = (V_pu * motor.V_rated) / np.sqrt(3)
            ws = 2 * np.pi * motor.freq / (motor.poles / 2)
            
            R_total = motor.Rs + (motor.Rr / s)
            X_total = motor.Xs + motor.Xr
            Z_sq = R_total**2 + X_total**2
            
            Te = (3 * V_phase**2 * (motor.Rr / s)) / (ws * Z_sq)
            
            if motor.load_type == 'constant_torque':
                Tm = (motor.P_rated * 1000) / ws 
            else:
                Tm = ((motor.P_rated * 1000) / ws) * ((1 - s)**2)
                
            return Te - Tm

        new_s = opt.fsolve(torque_mismatch, x0=motor.s)[0]
        motor.s = max(0.0001, min(new_s, 0.99))
    
    def check_overload(self):
        print("\n=== 🔍 TRANSMISSION LINE LOADING CHECK ===")
        overload_found = False
        
        for u, v, data in self.edges(data=True):
            branch = data.get('obj')
            if branch is None: continue
            
            from_bus = self.bus(u)
            to_bus = self.bus(v)
            
            V_i = from_bus.V * np.exp(1j * from_bus.theta)
            V_j = to_bus.V * np.exp(1j * to_bus.theta)
            
            I_ij = (V_i - V_j) * branch.Y
            S_ij = V_i * np.conj(I_ij)
            
            S_mag = abs(S_ij) * self.Sbase
            
            S_max = getattr(branch, 'S_max', None)
            
            if S_max is not None:
                loading_pct = (S_mag / S_max) * 100
                if loading_pct > 100:
                    print(f"🚨 OVERLOAD: {branch.name} ({from_bus.name} -> {to_bus.name})")
                    print(f"   Flow: {S_mag:.2f} MVA | Limit: {S_max} MVA | Loading: {loading_pct:.1f}%")
                    overload_found = True
            elif S_max == 0:
                print(f"🚨 ERROR: Smax = 0 at {branch.name} ({from_bus.name} -> {to_bus.name})")
                    
        if not overload_found:
            print("✅ All branch flows are within limits.")
    
    def solve(self, tol=1e-6, max_iter=100):
        """
        Fast Decoupled Power Flow (FDPF) Solver with Dynamic Q-Limit Enforcement
        """
        print("\n🚀 Starting Fast Decoupled Power Flow (FDPF)...")
        
        Ybus = self.build_ybus()
        B_mat = Ybus.imag
        G_mat = Ybus.real
        
        bus_idx = {bus.name: i for i, bus in enumerate(self.buses)}
        
        iteration = 0
        converged = False

        while iteration < max_iter:
            for bus in self.buses:
                for comp in bus.components:
                    if isinstance(comp, (AsynchronousMachine, Shunt)):
                        if isinstance(comp, AsynchronousMachine):
                            self.update_motor_slip(comp, bus.V)
                            comp.update_pq_from_slip(bus.V, self.Sbase)
                        elif isinstance(comp, Shunt):
                            comp.update_voltage_dependence(bus.V)
            
            pv_buses = [b for b in self.buses if b.type == 'PV']
            pq_buses = [b for b in self.buses if b.type == 'PQ']
            
            pv_pq_idx = [bus_idx[b.name] for b in pv_buses + pq_buses]
            pq_idx = [bus_idx[b.name] for b in pq_buses]

            B_prime = -B_mat[np.ix_(pv_pq_idx, pv_pq_idx)]
            B_double_prime = -B_mat[np.ix_(pq_idx, pq_idx)]

            dP_dict = {}
            dQ_dict = {}
            
            limit_hit_this_iter = False
            
            for i, bus in enumerate(self.buses):
                P_calc = 0.0
                Q_calc = 0.0
                for j, bus_k in enumerate(self.buses):
                    theta_ij = bus.theta - bus_k.theta
                    P_calc += bus.V * bus_k.V * (G_mat[i,j]*np.cos(theta_ij) + B_mat[i,j]*np.sin(theta_ij))
                    Q_calc += bus.V * bus_k.V * (G_mat[i,j]*np.sin(theta_ij) - B_mat[i,j]*np.cos(theta_ij))

                P_spec = sum(c.P for c in bus.components if hasattr(c, 'P')) / self.Sbase
                
                if bus.type in ['PV', 'PQ']:
                    dP_dict[bus.name] = P_spec - P_calc
                
                if bus.type == 'PQ':
                    Q_spec = sum(c.Q for c in bus.components if hasattr(c, 'Q')) / self.Sbase
                    dQ_dict[bus.name] = Q_spec - Q_calc

                if bus.type == 'PV':
                    Q_load = sum(c.Q for c in bus.components if isinstance(c, Load)) / self.Sbase
                    Q_gen_req = Q_calc - Q_load
                    
                    machines = [c for c in bus.components if isinstance(c, SynchronousMachine)]
                    if machines:
                        Qmax_total = sum(getattr(m, 'Qmax', 9999.0) for m in machines)
                        Qmin_total = sum(getattr(m, 'Qmin', -9999.0) for m in machines)
                        
                        Q_gen_req_mvar = Q_gen_req * self.Sbase
                        
                        if Q_gen_req_mvar > Qmax_total:
                            print(f"⚠️ Iteration {iteration}: Bus {bus.name} exceed Qmax ({Q_gen_req_mvar:.2f} > {Qmax_total}) -> Change to PQ Bus!")
                            bus.type = 'PQ'
                            for m in machines: m.Q = getattr(m, 'Qmax', 9999.0)
                            limit_hit_this_iter = True
                            
                        elif Q_gen_req_mvar < Qmin_total:
                            print(f"⚠️ Iteration {iteration}: Bus {bus.name} exceed Qmin ({Q_gen_req_mvar:.2f} < {Qmin_total}) -> Change to PQ Bus!")
                            bus.type = 'PQ'
                            for m in machines: m.Q = getattr(m, 'Qmin', -9999.0)
                            limit_hit_this_iter = True

            if limit_hit_this_iter:
                iteration += 1
                continue

            max_dP = max(abs(val) for val in dP_dict.values()) if dP_dict else 0
            max_dQ = max(abs(val) for val in dQ_dict.values()) if dQ_dict else 0
            
            if max_dP < tol and max_dQ < tol:
                print(f"✅ FDPF Converged in {iteration} iterations!")
                converged = True
                break

            if dP_dict:
                dP_vec = np.array([dP_dict[b.name] / b.V for b in pv_buses + pq_buses])
                dTheta = np.linalg.solve(B_prime, dP_vec)
                for idx, b_name in enumerate([b.name for b in pv_buses + pq_buses]):
                    b_obj = next(b for b in self.buses if b.name == b_name)
                    b_obj.theta += dTheta[idx]

            if dQ_dict:
                dQ_vec = np.array([dQ_dict[b.name] / b.V for b in pq_buses])
                dV = np.linalg.solve(B_double_prime, dQ_vec)
                for idx, b_name in enumerate([b.name for b in pq_buses]):
                    b_obj = next(b for b in self.buses if b.name == b_name)
                    b_obj.V += dV[idx]

            tap_changed = False
            for u, v, data in self.edges(data=True):
                branch = data.get('obj')
                if isinstance(branch, Transformer) and getattr(branch, 'auto_tap', False) and getattr(branch, 'controlled_bus', None):
                    target_bus = self.bus(branch.controlled_bus)
                    
                    deadband = 0.015 
                    
                    if target_bus.V < branch.target_V - deadband and branch.tap_ratio > branch.tap_min:
                        branch.tap_ratio -= branch.tap_step
                        branch.tap_ratio = max(branch.tap_ratio, branch.tap_min)
                        tap_changed = True
                        print(f"🔄 Iteration {iteration}: {branch.name} Step DOWN Tap -> {branch.tap_ratio:.4f} (Target V: {target_bus.V:.4f})")
                        
                    elif target_bus.V > branch.target_V + deadband and branch.tap_ratio < branch.tap_max:
                        branch.tap_ratio += branch.tap_step
                        branch.tap_ratio = min(branch.tap_ratio, branch.tap_max)
                        tap_changed = True
                        print(f"🔄 Iteration {iteration}: {branch.name} Step UP Tap -> {branch.tap_ratio:.4f} (Target V: {target_bus.V:.4f})")

            if tap_changed:
                Ybus = self.build_ybus()
                B_mat = Ybus.imag
                G_mat = Ybus.real
                iteration += 1
                continue
            
            iteration += 1

        if not converged:
            print(f"❌ FDPF Did Not Converge after {max_iter} iterations.")
    
    def eco_dispatch(self):
        machines = []
        machine_bus_idx = []
        P_load_bus = np.zeros(len(self.buses))
        Q_load_bus = np.zeros(len(self.buses))
        
        for i, bus in enumerate(self.buses):
            for _, data in bus.nodes(data=True):
                obj = data.get('obj')
                if type(obj).__name__ == 'SynchronousMachine':
                    machines.append(obj)
                    machine_bus_idx.append(i)
                elif type(obj).__name__ == 'Load':
                    P_load_bus[i] += abs(obj.P) 
                    Q_load_bus[i] += abs(obj.Q)
                    
        num_gen = len(machines)
        num_bus = len(self.buses)
        
        if self.Ybus is None:
            self.build_ybus()
            
        G = self.Ybus.real
        B = self.Ybus.imag
        
        idx_Pg = slice(0, num_gen)
        idx_Qg = slice(num_gen, 2*num_gen)
        idx_V = slice(2*num_gen, 2*num_gen + num_bus)
        idx_theta = slice(2*num_gen + num_bus, 2*num_gen + 2*num_bus)
        
        def objective(x):
            Pg = x[idx_Pg]
            cost = 0.0
            for i, m in enumerate(machines):
                if Pg[i] > 0:
                    cost += m.a + (m.b * Pg[i]) + (m.c * (Pg[i]**2))
            return cost

        def power_balance(x):
            Pg = x[idx_Pg]
            Qg = x[idx_Qg]
            V = x[idx_V]
            theta = x[idx_theta]
            
            P_gen_bus = np.zeros(num_bus)
            Q_gen_bus = np.zeros(num_bus)
            for i, bus_idx in enumerate(machine_bus_idx):
                P_gen_bus[bus_idx] += Pg[i]
                Q_gen_bus[bus_idx] += Qg[i]
            
            P_total_load_MVA = P_load_bus.copy() 
            Q_total_load_MVA = Q_load_bus.copy()
            
            for i, bus in enumerate(self.buses):
                for comp in bus.components:
                    if isinstance(comp, AsynchronousMachine):
                        comp.update_pq_from_slip(V[i], self.Sbase)
            
                        P_total_load_MVA[i] += abs(comp.P * self.Sbase)
                        Q_total_load_MVA[i] += abs(comp.Q * self.Sbase)
                    elif isinstance(comp, Shunt):
                        comp.update_voltage_dependence(V[i])
                        Q_gen_bus[i] += comp.Q
                
            P_inj_pu = (P_gen_bus - P_total_load_MVA) / self.Sbase
            Q_inj_pu = (Q_gen_bus - Q_total_load_MVA) / self.Sbase
            
            P_calc_pu = np.zeros(num_bus)
            Q_calc_pu = np.zeros(num_bus)
            
            for i in range(num_bus):
                for j in range(num_bus):
                    delta_ij = theta[i] - theta[j]
                    P_calc_pu[i] += V[i] * V[j] * (G[i, j] * np.cos(delta_ij) + B[i, j] * np.sin(delta_ij))
                    Q_calc_pu[i] += V[i] * V[j] * (G[i, j] * np.sin(delta_ij) - B[i, j] * np.cos(delta_ij))
                    
            mismatch_P = P_inj_pu - P_calc_pu
            mismatch_Q = Q_inj_pu - Q_calc_pu
            
            return np.concatenate((mismatch_P, mismatch_Q))
        
        def line_limits(x):
            V = x[idx_V]
            theta = x[idx_theta]
            V_complex = V * np.exp(1j * theta)
            
            margins = []
            nodes_list = list(self.nodes)
            
            for u, v, data in self.edges(data=True):
                obj = data.get('obj')
                if obj and getattr(obj, 'S_max', None) is not None:
                    i = nodes_list.index(u)
                    j = nodes_list.index(v)
                    
                    Vi = V_complex[i]
                    Vj = V_complex[j]
                    
                    if type(obj).__name__ == 'Transformer':
                        t = obj.tap_ratio * np.exp(1j * obj.phase_shift)
                        I_ij = (Vi/t - Vj) * obj.Y
                    else:
                        I_ij = (Vi - Vj) * obj.Y
                        
                    # Apparent Power (S = V * I*)
                    S_flow_pu = abs(Vi * np.conj(I_ij))
                    S_flow_MVA = S_flow_pu * self.Sbase
                    
                    margins.append(obj.S_max - S_flow_MVA)
            
            if not margins:
                return [1.0]
            return np.array(margins)

        bounds = []
        
        for m in machines:
            p_min = m.Pmin if m.Pmin != float('-inf') else 0.0
            p_max = m.Pmax if m.Pmax != float('inf') else 9999.0 
            bounds.append((p_min, p_max))
            
        for m in machines:
            q_min = m.Qmin if m.Qmin != float('-inf') else -9999.0
            q_max = m.Qmax if m.Qmax != float('inf') else 9999.0
            bounds.append((q_min, q_max))
            
        for bus in self.buses:
            bounds.append((0.95, 1.05))
            
        slack_idx = next(i for i, b in enumerate(self.buses) if b.type == 'Slack')
        for i in range(num_bus):
            if i == slack_idx:
                bounds.append((0.0, 0.0))
            else:
                bounds.append((-np.pi, np.pi))

        Pg0 = [m.Pmax / 2 if m.Pmax != float('inf') else 50.0 for m in machines]
        Qg0 = [0.0 for m in machines]
        V0 = [1.0 for b in self.buses]
        theta0 = [0.0 for b in self.buses]
        x0 = np.concatenate((Pg0, Qg0, V0, theta0))

        constraints = [
            {'type': 'eq', 'fun': power_balance},
            {'type': 'ineq', 'fun': line_limits}
        ]
        
        print("\n⏳ Running AC Optimal Power Flow...")
        result = opt.minimize(objective, x0, bounds=bounds, constraints=constraints, 
                              method='SLSQP', options={'maxiter': 500, 'ftol': 1e-6})
        
        if result.success:
            print("✅ ACOPF Converged Successfully!")
            Pg_opt = result.x[idx_Pg]
            Qg_opt = result.x[idx_Qg]
            V_opt = result.x[idx_V]
            theta_opt = result.x[idx_theta]
            
            for i, m in enumerate(machines):
                m.P = Pg_opt[i]
                m.Q = Qg_opt[i]
            for i, bus in enumerate(self.buses):
                bus.V = V_opt[i]
                bus.theta = theta_opt[i]
                
            print(f"Total Optimal Cost: ${result.fun:.2f}/hr")
            self.result()
        else:
            print("❌ ACOPF Failed to Converge!")
            print(result.message)