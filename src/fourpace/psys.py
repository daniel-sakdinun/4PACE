import yaml
import numpy as np
import networkx as nx
import scipy.optimize as opt

from pathlib import Path
from pandas.core.frame import DataFrame
from fourpace.model import (
    BusComponent, SynchronousMachine, AsynchronousMachine, Load, 
    Shunt, Inverter, Battery, BranchComponent, TransmissionLine, Transformer
)
from fourpace.control import SEXS, TGOV1, PSS1A
from fourpace.facts import CSVGN1, STATCOM1, TCSC1

class Bus(nx.Graph):
    """
    Representation of a Bus (Node) in a power system grid.
    Inherits from nx.Graph, allowing the bus to act as a subgraph containing connected components.
    """
    def __init__(self, name: str, Vbase: float, type: str = 'PQ'):
        """
        Parameters:
        -----------
        name : str
            The unique identifier/name of the bus.
        Vbase : float
            The base voltage of the bus in kV.
        type : str
            The bus type for Load Flow calculations ('PQ', 'PV', or 'Slack'). Default is 'PQ'.
        """
        super().__init__()
        self.name = name
        self.Vbase: float | None = Vbase
        self.type = type
        self.V: float = 1.0     # Voltage magnitude (pu)
        self.theta: float = 0.0 # Voltage phase angle (rad)
        self.add_node(name, obj=self)
    
    @property
    def P(self) -> float:
        """Calculates the net active power (P) injected into this bus by all connected components."""
        total_p = 0.0
        for comp in self.components:
            if hasattr(comp, 'P'):
                total_p += comp.P
        return total_p

    @property
    def Q(self) -> float:
        """Calculates the net reactive power (Q) injected into this bus by all connected components."""
        total_q = 0.0
        for comp in self.components:
            if hasattr(comp, 'Q'):
                total_q += comp.Q
        return total_q
    
    @property
    def S(self) -> complex:
        """Calculates the net apparent power (S) injected into this bus (P + jQ)."""
        return self.P + 1j * self.Q
    
    def component(self, name: str) -> BusComponent:
        """Retrieves a specific component connected to this bus by its name."""
        if name in self.nodes:
            return self.nodes[name]['obj']
        raise KeyError(f"❌ '{name}' not found in this bus.")
    
    @property
    def components(self) -> list:
        """Returns a list of all components (e.g., loads, generators) attached to this bus."""
        comp_list = []
        for _, data in self.nodes(data=True):
            obj = data.get('obj')
            if obj is not None and obj is not self:
                comp_list.append(obj)
        return comp_list
    
    def get(self):
        """Returns the current state variables of the bus: [V, theta, P, Q]."""
        return [self.V, self.theta, self.P, self.Q]
    
    def add_component(self, component: BusComponent):
        """Connects a single component to this bus."""
        self.add_node(component.name, obj=component)
        self.add_edge(self.name, component.name)
    
    def add_components(self, components: list[BusComponent]):
        """Connects multiple components to this bus simultaneously."""
        for component in components:
            self.add_component(component)
        
    def total_cost(self) -> float:
        """Calculates the total generation cost of all generators connected to this bus."""
        components = [component for _, component in self.nodes(data='obj')]
        components.remove(self)
        total = 0.0
        for com in components:
            total += com.cost()
        return total


class Grid(nx.Graph):
    """
    The main class representing the macro-level Power Grid network.
    Handles system topology, Y-bus matrix construction, and basic power flow operations.
    """
    def __init__(self, Sbase: float):
        """
        Parameters:
        -----------
        Sbase : float
            The system base apparent power in MVA.
        """
        super().__init__()
        self.Sbase: float = Sbase
        self.Ybus: np.ndarray | None = None
        self.load_profile: DataFrame = None
        self.T: int = 1
        
        self.series_facts = [] # Stores Series FACTS devices (e.g., TCSC) installed on branches.
    
    @classmethod
    def load(cls, filepath: str) -> 'Grid':
        """
        Loads the grid topology and parameters from a configuration file (YAML or JSON).
        
        Parameters:
        -----------
        filepath : str
            The path to the configuration file (e.g., 'config.yaml').
            
        Returns:
        --------
        Grid
            A fully initialized Grid object ready for simulation.
        """
        path = Path(filepath)
        with open(path, 'r', encoding='utf-8') as f:
            if path.suffix in ['.yaml', '.yml']:
                data = yaml.safe_load(f)
            else:
                raise ValueError("❌ Unsupported format! Please use .yaml or .json")

        grid = cls(Sbase=data.get('Sbase', 100.0))

        # Include Shunt FACTS in the standard bus component registry
        comp_classes = {
            'SynchronousMachine': SynchronousMachine,
            'AsynchronousMachine': AsynchronousMachine,
            'Load': Load,
            'Shunt': Shunt,
            'Inverter': Inverter,
            'Battery': Battery,
            'TransmissionLine': TransmissionLine,
            'Transformer': Transformer,
            'CSVGN1': CSVGN1,     # SVC
            'STATCOM1': STATCOM1  # STATCOM
        }

        ctrl_classes = {
            'SEXS': SEXS,
            'TGOV1': TGOV1,
            'PSS1A': PSS1A
        }
        
        # 1. Load Buses and Shunt Components (Including SVC & STATCOM)
        for b_data in data.get('buses', []):
            bus = Bus(name=b_data['name'], Vbase=b_data['Vbase'], type=b_data.get('bus_type', 'PQ'))
            grid.add_bus(bus)
            
            for comp_data in b_data.get('components', []):
                comp_type = comp_data.pop('type')
                
                if comp_type == 'SynchronousMachine':
                    avr_data = comp_data.pop('avr', None)
                    gov_data = comp_data.pop('gov', None)
                    pss_data = comp_data.pop('pss', None)
                    
                    avr_obj, gov_obj, pss_obj = None, None, None
                    gen_name = comp_data.get('name', 'UnknownGen')
                    
                    if avr_data:
                        ctrl_type = avr_data.pop('type')
                        if ctrl_type in ctrl_classes:
                            avr_obj = ctrl_classes[ctrl_type](name=f"AVR_{gen_name}", **avr_data)
                    
                    if gov_data:
                        ctrl_type = gov_data.pop('type')
                        if ctrl_type in ctrl_classes:
                            gov_obj = ctrl_classes[ctrl_type](name=f"GOV_{gen_name}", **gov_data)
                            
                    if pss_data:
                        ctrl_type = pss_data.pop('type')
                        if ctrl_type in ctrl_classes:
                            pss_obj = ctrl_classes[ctrl_type](name=f"PSS_{gen_name}", **pss_data)
                            
                    comp_data['avr'] = avr_obj
                    comp_data['gov'] = gov_obj
                    comp_data['pss'] = pss_obj

                if comp_type in comp_classes:
                    comp_obj = comp_classes[comp_type](**comp_data) 
                    grid.bus(bus.name).add_component(comp_obj)

        # 2. Load Branches
        for branch_data in data.get('branches', []):
            branch_type = branch_data.pop('type')
            from_bus = branch_data.pop('from_bus')
            to_bus = branch_data.pop('to_bus')
            
            if branch_type in comp_classes:
                branch_obj = comp_classes[branch_type](**branch_data)
                grid.connect(from_bus, to_bus, branch_obj)

        # 3. ADDED: Load Series FACTS Devices (TCSC)
        series_classes = {
            'TCSC1': TCSC1
        }
        
        for facts_data in data.get('series_facts', []):
            facts_type = facts_data.pop('type')
            if facts_type in series_classes:
                facts_obj = series_classes[facts_type](**facts_data)
                grid.series_facts.append(facts_obj)

        print(f"✅ Successfully loaded grid from {path.name}")
        return grid
    
    def add_bus(self, bus: Bus):
        """Adds a single bus to the grid."""
        self.add_node(bus.name, bus=bus)
    
    def add_busses(self, busses: list[Bus]):
        """Adds a list of buses to the grid."""
        for bus in busses:
            self.add_bus(bus)
        
    def connect(self, from_bus: str, to_bus: str, branch: BranchComponent):
        """
        Connects two buses via a branch component (e.g., Transmission Line or Transformer).
        """
        self.add_edge(from_bus, to_bus, obj=branch)
    
    def bus(self, name: str) -> Bus:
        """Retrieves a Bus object by its name."""
        if name in self.nodes:
            return self.nodes[name]['bus']
        raise KeyError(f"❌ Bus '{name}' not found in this grid.")
    
    @property
    def buses(self) -> list[Bus]:
        """Returns a list of all Bus objects in the system."""
        return [bus for _, bus in self.nodes(data='bus')]
    
    def build_ybus(self) -> np.ndarray:
        """
        Constructs the Steady-State Admittance Matrix (Y-bus).
        Used for standard Load Flow calculations.
        """
        nodes = list(self.nodes)
        n = len(nodes)
        Y = np.zeros((n, n), dtype=complex)
        
        for u, v, data in self.edges(data=True):
            i = nodes.index(u)
            j = nodes.index(v)
            obj = data.get('obj')
            
            if obj is None: 
                continue
            
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
    
    def build_ybus_pos(self):
        """
        Constructs the Positive-Sequence Y-bus matrix (Y1).
        Incorporates generator subtransient reactances (Xd").
        Used for fault calculation and transient stability analysis.
        """
        num_bus = len(self.buses)
        nodes_list = list(self.nodes)
        Ybus_fault = np.zeros((num_bus, num_bus), dtype=complex)

        for u, v, data in self.edges(data=True):
            i = nodes_list.index(u)
            j = nodes_list.index(v)
            branch = data.get('obj')
            
            Y = branch.Y
            Ybus_fault[i, j] -= Y
            Ybus_fault[j, i] -= Y
            Ybus_fault[i, i] += Y
            Ybus_fault[j, j] += Y
            
            if isinstance(branch, TransmissionLine) and getattr(branch, 'B_shunt', 0) > 0:
                Y_shunt = 1j * (branch.B_shunt / 2)
                Ybus_fault[i, i] += Y_shunt
                Ybus_fault[j, j] += Y_shunt
                
            elif isinstance(branch, Transformer) and branch.tap_ratio != 1.0:
                a = branch.tap_ratio
                Ybus_fault[i, i] += Y * ((1 - a) / (a**2))
                Ybus_fault[j, j] += Y * ((a - 1) / a)

        for i, bus in enumerate(self.buses):
            for comp in bus.components:
                if isinstance(comp, SynchronousMachine):
                    Xd_subtransient = getattr(comp, 'Xd_sub', 0.2) 
                    Y_gen = 1 / (1j * Xd_subtransient)
                    Ybus_fault[i, i] += Y_gen

        return Ybus_fault


    def build_ybus_zero(self):
        """
        Constructs the Zero-Sequence Y-bus Matrix (Y0).
        *** CRITICAL: This structure is highly dependent on transformer grounding/connection types. ***
        Used for asymmetrical fault analysis (e.g., SLG, DLG).
        """
        num_bus = len(self.buses)
        nodes_list = list(self.nodes)
        Ybus_zero = np.zeros((num_bus, num_bus), dtype=complex)

        for u, v, data in self.edges(data=True):
            i = nodes_list.index(u)
            j = nodes_list.index(v)
            branch = data.get('obj')
            
            if isinstance(branch, TransmissionLine):
                Y0 = getattr(branch, 'Y0', branch.Y) 
                Ybus_zero[i, j] -= Y0
                Ybus_zero[j, i] -= Y0
                Ybus_zero[i, i] += Y0
                Ybus_zero[j, j] += Y0
                
                if getattr(branch, 'B0_shunt', 0) > 0:
                    Y0_shunt = 1j * (branch.B0_shunt / 2)
                    Ybus_zero[i, i] += Y0_shunt
                    Ybus_zero[j, j] += Y0_shunt

            elif isinstance(branch, Transformer):
                conn_attr = getattr(branch, 'connection_type', 'yg-yg')
                conn = conn_attr.lower().split('-')
                
                if len(conn) != 2:
                    pri_type, sec_type = 'yg', 'yg'
                else:
                    pri_type, sec_type = conn[0], conn[1]
                    
                Y0 = getattr(branch, 'Y0', branch.Y)
                
                if pri_type == 'yg' and sec_type == 'yg':
                    Ybus_zero[i, j] -= Y0
                    Ybus_zero[j, i] -= Y0
                    Ybus_zero[i, i] += Y0
                    Ybus_zero[j, j] += Y0
                elif pri_type == 'delta' and sec_type == 'yg':
                    Ybus_zero[j, j] += Y0
                elif pri_type == 'yg' and sec_type == 'delta':
                    Ybus_zero[i, i] += Y0
                
        for i, bus in enumerate(self.buses):
            for comp in bus.components:
                if isinstance(comp, SynchronousMachine):
                    X0 = getattr(comp, 'X0', 0.05)
                    if X0 > 0:
                        Y0_gen = 1 / (1j * X0)
                        Ybus_zero[i, i] += Y0_gen

        # Artificial shunt to prevent singular matrices in grids with floating neutrals
        Ybus_zero += np.eye(num_bus) * 1e-6j
        
        return Ybus_zero
    
    def kron_reduction(self, Y_matrix: np.ndarray, keep_nodes: list[str]) -> np.ndarray:
        """
        Reduces the dimensions of the Y-bus matrix using Kron Reduction (Network Equivalencing).
        Eliminates non-essential nodes while preserving the electrical characteristics 
        between the specified 'keep_nodes'. Essential for Transient Stability Analysis.
        
        Parameters:
        -----------
        Y_matrix : np.ndarray
            The original admittance matrix.
        keep_nodes : list[str]
            A list of bus names to retain (e.g., buses with generators).
            
        Returns:
        --------
        np.ndarray
            The reduced admittance matrix.
        """
        nodes_list = list(self.nodes)
        n = Y_matrix.shape[0]
        
        # Find the indices of the nodes to keep and the nodes to eliminate
        keep_indices = [nodes_list.index(name) for name in keep_nodes if name in nodes_list]
        eliminate_indices = [i for i in range(n) if i not in keep_indices]
        
        Y_new = Y_matrix.copy()
        
        # Loop to eliminate nodes one by one (from highest index to lowest to prevent index shifting)
        for k in sorted(eliminate_indices, reverse=True):
            Y_kk = Y_new[k, k]
            if abs(Y_kk) < 1e-6:
                continue # Skip unconnected nodes (Floating nodes)
                
            Y_ik = Y_new[:, k:k+1]
            Y_kj = Y_new[k:k+1, :]
            
            # Kron Reduction equation
            Y_new = Y_new - (Y_ik @ Y_kj) / Y_kk
            
            # Delete row and column k
            Y_new = np.delete(Y_new, k, axis=0)
            Y_new = np.delete(Y_new, k, axis=1)
            
        return Y_new
    
    def result(self):
        """Prints the load flow calculation results (Voltage, Angle, Active & Reactive Power) to the console."""
        print("\n=== 📊 POWER FLOW RESULTS ===")
        for i, bus in enumerate(self.buses):
            P_pu, Q_pu = self.calculate_PQ(i)
            P_actual = P_pu * self.Sbase
            Q_actual = Q_pu * self.Sbase
            print(f"Bus {bus.name} | V = {bus.V:.4f} pu | phase = {np.rad2deg(bus.theta):.2f}° | P = {P_actual:7.2f} MW | Q = {Q_actual:7.2f} MVAr")
    
    def loading_status(self):
        """Checks and prints the loading status (overload warnings) for all branches and generators."""
        grid_status = self.check_overload()
        print("\n📊 Grid Loading Status:")
        print("--- Branches ---")
        for name, pct in grid_status['branches'].items():
            status = "🔥 OVERLOAD!" if pct > 100 else "✅ OK"
            print(f"{name}: {pct}% {status}")

        print("\n--- Generators ---")
        for name, pct in grid_status['generators'].items():
            status = "🔥 OVERLOAD!" if pct > 100 else "✅ OK"
            print(f"{name}: {pct}% {status}")
    
    def attach_profile(self, df: DataFrame):
        """
        Attaches a time-series profile (e.g., hourly load or solar generation profile) to the grid.
        
        Parameters:
        -----------
        df : pandas.DataFrame
            The time-series data profile.
        """
        self.load_profile = df
        self.T = len(df)
        print(f"📅 Attached Load/Generation Profile with {self.T} time steps.")
    
    def apply_profile(self, step_or_dict):
        """
        Updates the P and Q values of grid loads based on the specified time step in the attached profile.
        
        Parameters:
        -----------
        step_or_dict : int or dict
            The hour index (int) to fetch from the DataFrame, or a dictionary specifying load multipliers directly.
        """
        if isinstance(step_or_dict, int):
            if self.load_profile is None:
                raise ValueError("❌ No profile attached to Grid! Use grid.attach_profile(df) first.")
            data_dict = self.load_profile.iloc[step_or_dict].to_dict()
        else:
            data_dict = step_or_dict
        
        for bus in self.buses:
            for comp in bus.components:
                if type(comp).__name__ == 'Load':
                    if comp.name in data_dict:
                        multiplier = data_dict[comp.name]
                        comp.P = comp.P_nom * multiplier
                        comp.Q = comp.Q_nom * multiplier
        
    def get_peak_load_hour(self) -> int:
        """Finds and returns the hour (index) at which the system experiences the maximum total active power load."""
        if self.load_profile is None:
            raise ValueError("❌ No profile attached to Grid! Use grid.attach_profile(df) first.")

        max_load = -1.0
        peak_hour = 0

        for t in range(self.T):
            current_total_load = 0.0
            data_dict = self.load_profile.iloc[t].to_dict()

            for bus in self.buses:
                for comp in bus.components:
                    comp_type = type(comp).__name__
                    if comp_type == 'Load':
                        multiplier = data_dict.get(comp.name, 1.0)
                        current_total_load += abs(comp.P_nom * multiplier)
                    elif comp_type == 'AsynchronousMachine':
                        multiplier = data_dict.get(comp.name, 1.0)
                        current_total_load += abs(comp.P_rated_base * multiplier)

            if current_total_load > max_load:
                max_load = current_total_load
                peak_hour = t

        print(f"📈 Peak Load Detected at Hour {peak_hour} (Total System Load: {max_load:.2f} MW)")
        return peak_hour
    
    def calculate_PQ(self, i: int) -> tuple[float, float]:
        """
        Calculates the active (P) and reactive (Q) power injected at bus 'i' using power flow equations and the Y-bus matrix.
        
        Parameters:
        -----------
        i : int
            The index of the bus in the `self.buses` list.
            
        Returns:
        --------
        tuple[float, float]
            (P_i, Q_i) values in per-unit (p.u.).
        """
        bus_i = self.buses[i]
        Vi = bus_i.V
        theta_i = bus_i.theta

        G = self.Ybus.real
        B = self.Ybus.imag

        P_i = 0.0
        Q_i = 0.0

        for j in range(len(self.buses)):
            bus_j = self.buses[j]
            Vj = bus_j.V
            theta_j = bus_j.theta
            delta_ij = theta_i - theta_j

            P_i += Vi * Vj * (G[i, j] * np.cos(delta_ij) + B[i, j] * np.sin(delta_ij))
            Q_i += Vi * Vj * (G[i, j] * np.sin(delta_ij) - B[i, j] * np.cos(delta_ij))
        
        return P_i, Q_i
    
    def update_motor_slip(self, motor: AsynchronousMachine, V_pu: float):
        """
        Updates the slip of an induction motor to find the equilibrium point where 
        electromagnetic torque (Te) equals mechanical torque (Tm).
        
        Parameters:
        -----------
        motor : AsynchronousMachine
            The induction motor object.
        V_pu : float
            The current terminal voltage (p.u.).
        """
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

        # Root finding to determine the new equilibrium slip
        new_s = opt.fsolve(torque_mismatch, x0=motor.s)[0]
        motor.s = max(0.0001, min(new_s, 0.99))
    
    def check_overload(self) -> dict:
        """
        Inspects the loading status of all branches and generators in the system.
        
        Returns:
        --------
        dict
            A dictionary containing loading percentages categorized by 'branches' and 'generators'.
        """
        loadings = {
            'branches': {},
            'generators': {}
        }
        
        # 1. Branches Check
        for u, v, data in self.edges(data=True):
            branch = data.get('obj')
            if branch is None or getattr(branch, 'S_max', None) is None:
                continue
                
            from_bus = self.bus(u)
            to_bus = self.bus(v)
            
            V_i = from_bus.V * np.exp(1j * from_bus.theta)
            V_j = to_bus.V * np.exp(1j * to_bus.theta)
            
            if type(branch).__name__ == 'Transformer':
                t = branch.tap_ratio * np.exp(1j * branch.phase_shift)
                I_ij = (V_i / t - V_j) * branch.Y
            else:
                I_ij = (V_i - V_j) * branch.Y
                
            S_flow_MVA = abs(V_i * np.conj(I_ij)) * self.Sbase
            loading_pct = (S_flow_MVA / branch.S_max) * 100
            loadings['branches'][branch.name] = round(loading_pct, 2)

        # 2. Generator Check
        for bus in self.buses:
            for comp in bus.components:
                if type(comp).__name__ == 'SynchronousMachine':
                    S_gen_MVA = np.sqrt(comp.P**2 + comp.Q**2)
                    s_rated = getattr(comp, 'S_rated', None)
                    
                    if s_rated is None:
                        p_max = getattr(comp, 'Pmax', 9999.0)
                        q_max = getattr(comp, 'Qmax', 9999.0)
                        
                        if p_max == float('inf') or q_max == float('inf') or p_max == 9999.0:
                            s_rated = 9999.0
                        else:
                            s_rated = np.sqrt(p_max**2 + q_max**2)
                            
                    if s_rated > 0 and s_rated != 9999.0:
                        loading_pct = (S_gen_MVA / s_rated) * 100
                        loadings['generators'][comp.name] = round(loading_pct, 2)
                    else:
                        loadings['generators'][comp.name] = 0.0
                        
        return loadings