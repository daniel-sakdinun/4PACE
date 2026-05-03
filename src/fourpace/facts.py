from abc import ABC, abstractmethod
import numpy as np
from fourpace.model import BusComponent, BranchComponent

# =====================================================================
# Base Class for Shunt FACTS Devices
# =====================================================================
class ShuntFACTS(BusComponent, ABC):
    """
    Abstract base class for all Shunt Flexible AC Transmission System (FACTS) devices.
    Designed for shunt-connected dynamic compensators like SVCs and STATCOMs that inject or absorb reactive power to regulate bus voltage.
    """
    def __init__(self, name: str, P: float = 0.0, Q: float = 0.0):
        # FACTS devices generally don't consume real power (ignoring tiny losses)
        # We initialize R=0, X=0 because their admittance is dynamic.
        super().__init__(name, P=0.0, Q=0.0, R=0.0, X=0.0)
        self.n_states = 0
        
    @abstractmethod
    def initialize(self, V_initial: float) -> np.ndarray:
        """
        Calculates the initial steady-state values for the device's internal states based on the initial bus voltage magnitude.
        """
        pass

    @abstractmethod
    def get_derivatives(self, local_state: np.ndarray, V_bus_mag: float) -> np.ndarray:
        """
        Computes the time derivatives (dx/dt) for the internal state variables given the current bus voltage magnitude.
        """
        pass
        
    @abstractmethod
    def get_susceptance(self, local_state: np.ndarray) -> float:
        """ Returns the equivalent dynamic susceptance $B$ (p.u.) provided by the device to the network. """
        pass

    # Required by BusComponent, but usually 0 for FACTS
    def cost(self) -> float: return 0.0
    def incremental_cost(self) -> float: return 0.0

# =====================================================================
# Base Class for Series FACTS Devices
# =====================================================================

class SeriesFACTS(ABC):
    """
    Base class for dynamic series compensators (TCSC, SSSC)
    These do not sit on a bus; they are attached to a specific Branch/Line.
    """
    def __init__(self, name: str, branch_name: str):
        self.name = name
        self.branch_name = branch_name  # The TransmissionLine this sits on
        self.n_states = 0
        
    @abstractmethod
    def initialize(self, *args, **kwargs) -> np.ndarray:
        pass

    @abstractmethod
    def get_derivatives(self, local_state: np.ndarray, P_line_flow: float) -> np.ndarray:
        pass

# =====================================================================
# Concrete IEEE Standard SVC Model (CSVGN1)
# =====================================================================
class CSVGN1(ShuntFACTS):
    """
    Standard Static Var Compensator (SVC) Model.
    Model with 1 State Variable: [B_svc]
    - B_svc: The dynamic susceptance of the device.
    """
    def __init__(self, name: str, V_ref: float = 1.0, 
                 K_svc: float = 100.0, T_1: float = 0.05, 
                 B_max: float = 1.0, B_min: float = -1.0):
        super().__init__(name)
        self.n_states = 1
        
        self.V_ref = V_ref
        self.K_svc = K_svc    # Control Gain
        self.T_1 = T_1        # Time constant (firing delay)
        
        # B_max > 0 means capacitive (injecting Q)
        # B_min < 0 means inductive (absorbing Q)
        self.B_max = B_max
        self.B_min = B_min
        
    def initialize(self, V_initial: float) -> np.ndarray:
        """ 
        Calculate initial steady-state B_svc required to maintain V_initial.
        Usually, in power flow, the SVC might be at B=0 if voltage is fine.
        """
        # For simplicity, assume it starts at 0 or whatever is needed.
        # A rigorous power flow would dictate this, but we'll assume 0 for now.
        B_init = 0.0
        return np.array([B_init])

    def get_derivatives(self, local_state: np.ndarray, V_bus_mag: float) -> np.ndarray:
        B_svc = local_state[0]
        
        # Error signal (V_ref - V_bus)
        error = self.V_ref - V_bus_mag
        
        # First order delay block
        dB = (self.K_svc * error - B_svc) / self.T_1
        
        # Anti-windup Limiter
        if B_svc >= self.B_max and dB > 0: dB = 0.0
        if B_svc <= self.B_min and dB < 0: dB = 0.0
            
        return np.array([dB])
        
    def get_susceptance(self, local_state: np.ndarray) -> float:
        return local_state[0]

# =====================================================================
# Concrete STATCOM Model (VSC-based Shunt)
# =====================================================================
class STATCOM1(ShuntFACTS):
    """
    Generic Static Synchronous Compensator (STATCOM) Model.
    Model with 1 State Variable: [I_q]
    - I_q: Reactive current injection (Voltage Source Converter acting as a current source)
    """
    def __init__(self, name: str, V_ref: float = 1.0, 
                 K_r: float = 50.0, T_r: float = 0.05, 
                 Iq_max: float = 1.0, Iq_min: float = -1.0):
        super().__init__(name)
        self.n_states = 1
        
        self.V_ref = V_ref
        self.K_r = K_r        # Regulator Gain
        self.T_r = T_r        # Converter Delay (very fast, e.g., 0.02 - 0.05s)
        
        # Iq_max > 0 means capacitive (injecting reactive current)
        # Iq_min < 0 means inductive (absorbing reactive current)
        self.Iq_max = Iq_max
        self.Iq_min = Iq_min
        
    def initialize(self, V_initial: float) -> np.ndarray:
        """ STATCOM starts with 0 reactive current injection unless power flow dictates otherwise. """
        return np.array([0.0])

    def get_derivatives(self, local_state: np.ndarray, V_bus_mag: float) -> np.ndarray:
        I_q = local_state[0]
        
        # Voltage error
        error = self.V_ref - V_bus_mag
        
        # First order PI/Delay block for current command
        dI_q = (self.K_r * error - I_q) / self.T_r
        
        # Hard limits on converter current rating (Anti-windup)
        if I_q >= self.Iq_max and dI_q > 0: dI_q = 0.0
        if I_q <= self.Iq_min and dI_q < 0: dI_q = 0.0
            
        return np.array([dI_q])
        
    def get_susceptance(self, local_state: np.ndarray) -> float:
        """ STATCOMs are modeled as current sources, NOT susceptances. """
        return 0.0 
        
    def get_Iq(self, local_state: np.ndarray) -> float:
        """ Returns the dynamic reactive current injection. """
        return local_state[0]

# =====================================================================
# Concrete TCSC Model (Thyristor Controlled Series Capacitor)
# =====================================================================
class TCSC1(SeriesFACTS):
    """
    Generic Thyristor Controlled Series Capacitor (TCSC).
    Model with 1 State Variable: [X_tcsc]
    - X_tcsc: The dynamic series reactance added to the line.
    """
    def __init__(self, name: str, branch_name: str, P_ref: float, 
                 K_p: float = 1.0, T_p: float = 0.05, 
                 X_max: float = 0.0, X_min: float = -0.5):
        super().__init__(name, branch_name)
        self.n_states = 1
        
        self.P_ref = P_ref    # The desired active power flow through the line
        self.K_p = K_p        # Power controller gain
        self.T_p = T_p        # Thyristor firing delay
        
        # X_min is usually negative (capacitive, reducing total line reactance)
        # X_max is usually 0 or slightly positive (inductive)
        self.X_max = X_max
        self.X_min = X_min
        
    def initialize(self, *args, **kwargs) -> np.ndarray:
        """ Starts with 0 compensation. """
        return np.array([0.0])

    def get_derivatives(self, local_state: np.ndarray, P_line_flow: float) -> np.ndarray:
        X_tcsc = local_state[0]
        
        # Error between desired line flow and actual line flow
        error = self.P_ref - P_line_flow
        
        # Calculate derivative for the required reactance
        # Note: If P_flow < P_ref, we want MORE power, so we need LESS reactance (more negative X).
        # Therefore, we subtract the error.
        dX_tcsc = (-self.K_p * error - X_tcsc) / self.T_p
        
        # Anti-windup Limiter
        if X_tcsc >= self.X_max and dX_tcsc > 0: dX_tcsc = 0.0
        if X_tcsc <= self.X_min and dX_tcsc < 0: dX_tcsc = 0.0
            
        return np.array([dX_tcsc])
        
    def get_X_series(self, local_state: np.ndarray) -> float:
        """ Returns the dynamic series reactance. """
        return local_state[0]