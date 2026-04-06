import yaml
from pathlib import Path
import numpy as np
import scipy.optimize as opt
import networkx as nx

from model import BusComponent, SynchronousMachine, Load, BranchComponent, TransmissionLine, Transformer

class Bus(nx.Graph):
    def __init__(self, name: str, Vbase: float, bus_type: str = 'PQ'):
        super().__init__()
        self.name = name
        self.Vbase: float | None = Vbase
        
        self.type = bus_type
        self.v:float = 1.0
        self.phi:float = 0.0
        self.add_node(name, obj=self)
    
    @property
    def P(self) -> float:
        total_p = 0.0
        for _, data in self.nodes(data=True):
            obj = data.get('obj')
            if obj is not self and hasattr(obj, 'P'):
                total_p += obj.P
        return total_p

    @property
    def Q(self) -> float:
        total_q = 0.0
        for _, data in self.nodes(data=True):
            obj = data.get('obj')
            if obj is not self and hasattr(obj, 'Q'):
                total_q += obj.Q
        return total_q
    
    @property
    def S(self) -> complex:
        return self.P + 1j*self.Q;
    
    def flat_start(self):
        self.V = 1.0
        self.phi = 0.0
        
    def get(self):
        return [self.v, self.phi, self.P, self.Q]
    
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
            'Load': Load,
            'TransmissionLine': TransmissionLine,
            'Transformer': Transformer
        }

        for b_data in data.get('buses', []):
            bus = Bus(name=b_data['name'], Vbase=b_data['Vbase'], bus_type=b_data.get('bus_type', 'PQ'))
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
    def busses(self) -> list[Bus]:
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
                
        self.Ybus = Y
        return Y
    
    def result(self):
        print("\n=== 📊 POWER FLOW RESULTS ===")
        for i, bus in enumerate(self.busses):
            P_pu, Q_pu = self.calculate_PQ(i)
            
            P_actual = P_pu * self.Sbase
            Q_actual = Q_pu * self.Sbase
            
            print(f"Bus {bus.name} | V = {bus.v:.4f} pu | phi = {bus.phi:8.4f} rad | P = {P_actual:7.2f} MW | Q = {Q_actual:7.2f} MVAr")
    
    def flat_start(self):
        for bus in self.busses:
            bus.flat_start()
    
    def update_busses(self, correction: np.ndarray, damping = .5):
        pv_pq_idx = [i for i, b in enumerate(self.busses) if b.type in ['PV', 'PQ']]
        pq_idx = [i for i, b in enumerate(self.busses) if b.type == 'PQ']
        
        n1 = len(pv_pq_idx)
        
        for row_i, i in enumerate(pv_pq_idx):
            self.busses[i].phi += damping * float(correction[row_i])
            
        for row_i, i in enumerate(pq_idx):
            self.busses[i].v += damping * float(correction[n1 + row_i])
    
    def calculate_PQ(self, i):
        bus_i = self.busses[i]
        Vi = bus_i.v
        phi_i = bus_i.phi
        
        G = self.Ybus.real
        B = self.Ybus.imag
        
        P_i = 0
        Q_i = 0
        
        for j in range(len(self.busses)):
            bus_j = self.busses[j]
            Vj = bus_j.v
            phi_j = bus_j.phi
            delta_ij = phi_i - phi_j

            P_i += Vi * Vj * (G[i, j] * np.cos(delta_ij) + B[i, j] * np.sin(delta_ij))
            Q_i += Vi * Vj * (G[i, j] * np.sin(delta_ij) - B[i, j] * np.cos(delta_ij))

        return P_i, Q_i

    def jacobian(self):
        n = len(self.busses)
        G = self.Ybus.real
        B = self.Ybus.imag
        
        pv_pq_idx = [i for i, b in enumerate(self.busses) if b.type in ['PV', 'PQ']]
        pq_idx = [i for i, b in enumerate(self.busses) if b.type == 'PQ']
        
        n1 = len(pv_pq_idx)
        n2 = len(pq_idx)
        J = np.zeros((n1 + n2, n1 + n2))

        for row_i, i in enumerate(pv_pq_idx):
            Vi = self.busses[i].v
            phi_i = self.busses[i].phi
            
            for col_j, j in enumerate(pv_pq_idx):
                Vj = self.busses[j].v
                phi_j = self.busses[j].phi
                delta_ij = phi_i - phi_j
                
                if i == j: # Diagonal elements
                    # dPi/dphi_i = -Qi - (Vi^2 * Bii)
                    _, Qi_calc = self.calculate_PQ(i)
                    J[row_i, col_j] = -Qi_calc - (Vi**2 * B[i, i])
                else: # Off-diagonal
                    # dPi/dphi_j = Vi*Vj*(Gij*sin(d_ij) - Bij*cos(d_ij))
                    J[row_i, col_j] = Vi * Vj * (G[i, j] * np.sin(delta_ij) - B[i, j] * np.cos(delta_ij))

        for row_i, i in enumerate(pv_pq_idx):
            Vi = self.busses[i].v
            phi_i = self.busses[i].phi
            
            for col_j, j in enumerate(pq_idx):
                Vj = self.busses[j].v
                phi_j = self.busses[j].phi
                delta_ij = phi_i - phi_j
                
                if i == j: # Diagonal
                    # dPi/dVi = (Pi/Vi) + (Gii * Vi)
                    Pi_calc, _ = self.calculate_PQ(i)
                    J[row_i, n1 + col_j] = (Pi_calc / Vi) + (G[i, i] * Vi)
                else: # Off-diagonal
                    # dPi/dVj = Vi*(Gij*cos(d_ij) + Bij*sin(d_ij))
                    J[row_i, n1 + col_j] = Vi * (G[i, j] * np.cos(delta_ij) + B[i, j] * np.sin(delta_ij))

        for row_i, i in enumerate(pq_idx):
            Vi = self.busses[i].v
            phi_i = self.busses[i].phi
            Pi_calc, Qi_calc = self.calculate_PQ(i)
            actual_row = n1 + row_i
            
            # dQ/dphi
            for col_j, j in enumerate(pv_pq_idx):
                Vj = self.busses[j].v
                phi_j = self.busses[j].phi
                delta_ij = phi_i - phi_j
                
                if i == j:
                    # dQi/dphi_i = Pi - (Vi^2 * Gii)
                    Pi_calc, _ = self.calculate_PQ(i)
                    J[actual_row, col_j] = Pi_calc - (Vi**2 * G[i, i])
                else:
                    # dQi/dphi_j = -Vi*Vj*(Gij*cos(d_ij) + Bij*sin(d_ij))
                    J[actual_row, col_j] = -Vi * Vj * (G[i, j] * np.cos(delta_ij) + B[i, j] * np.sin(delta_ij))
            
            # dQ/dV
            for col_j, j in enumerate(pq_idx):
                Vj = self.busses[j].v
                phi_j = self.busses[j].phi
                delta_ij = phi_i - phi_j
                
                if i == j:
                    J[actual_row, n1 + col_j] = (Qi_calc / Vi) - (B[i, i] * Vi)
                else:
                    delta_ij = phi_i - phi_j
                    J[actual_row, n1 + col_j] = Vi * (G[i, j] * np.sin(delta_ij) - B[i, j] * np.cos(delta_ij))
                    
        return J
    
    def mismatch(self):
        delta_P = []
        delta_Q = []
        
        for i, bus in enumerate(self.busses):
            if bus.type == 'Slack':
                continue
                
            P_calc, Q_calc = self.calculate_PQ(i)
            
            delta_P.append((bus.P/self.Sbase) - P_calc)
            
            if bus.type == 'PQ':
                delta_Q.append((bus.Q/self.Sbase) - Q_calc)
        return np.array(delta_P + delta_Q).flatten()
    
    def solve(self, tolerance=1e-6, max_iteration=100):
        self.flat_start()
        self.build_ybus()
        
        for i in range(max_iteration):
            m = self.mismatch()
            if np.max(np.abs(m)) < tolerance:
                print(f"\n✅ Power Flow Converged in {i} iterations!")
                break
            
            dx = np.linalg.solve(self.jacobian(), m)
            self.update_busses(dx)
        else:
            print(f"\n⚠️ Warning: Power Flow did not converge within max iterations ({max_iteration}).")
        
        self.result()
    
    def eco_dispatch(self):
        machines = []
        machine_bus_idx = []
        P_load_bus = np.zeros(len(self.busses))
        Q_load_bus = np.zeros(len(self.busses))
        
        for i, bus in enumerate(self.busses):
            for _, data in bus.nodes(data=True):
                obj = data.get('obj')
                if type(obj).__name__ == 'SynchronousMachine':
                    machines.append(obj)
                    machine_bus_idx.append(i)
                elif type(obj).__name__ == 'Load':
                    P_load_bus[i] += abs(obj.P) 
                    Q_load_bus[i] += abs(obj.Q)
                    
        num_gen = len(machines)
        num_bus = len(self.busses)
        
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
                
            P_inj_pu = (P_gen_bus - P_load_bus) / self.Sbase
            Q_inj_pu = (Q_gen_bus - Q_load_bus) / self.Sbase
            
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
            
        for bus in self.busses:
            bounds.append((0.95, 1.05))
            
        slack_idx = next(i for i, b in enumerate(self.busses) if b.type == 'Slack')
        for i in range(num_bus):
            if i == slack_idx:
                bounds.append((0.0, 0.0))
            else:
                bounds.append((-np.pi, np.pi))

        Pg0 = [m.Pmax / 2 if m.Pmax != float('inf') else 50.0 for m in machines]
        Qg0 = [0.0 for m in machines]
        V0 = [1.0 for b in self.busses]
        theta0 = [0.0 for b in self.busses]
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
            for i, bus in enumerate(self.busses):
                bus.v = V_opt[i]
                bus.phi = theta_opt[i]
                
            print(f"Total Optimal Cost: ${result.fun:.2f}/hr")
            self.result()
        else:
            print("❌ ACOPF Failed to Converge!")
            print(result.message)