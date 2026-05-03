from abc import ABC, abstractmethod
import numpy as np

# =====================================================================
# Base Classes
# =====================================================================
class AVR(ABC):
    """ Base template for all types of Automatic Voltage Regulators (AVR) """
    def __init__(self, name: str):
        self.name: str = name
        self.n_states: int = 0  
        
    @abstractmethod
    def initialize(self, V_ref: float, V_t: float, E_fd0: float, *args, **kwargs) -> np.ndarray:
        pass

    @abstractmethod
    def get_derivatives(self, local_state: np.ndarray, V_ref: float, V_t: float, V_pss: float, *args, **kwargs) -> np.ndarray:
        pass
        
    @abstractmethod
    def get_Efd(self, local_state: np.ndarray) -> float:
        pass

class Governor(ABC):
    """ Base template for all types of Turbine-Governors """
    def __init__(self, name: str):
        self.name: str = name
        self.n_states: int = 0  
        
    @abstractmethod
    def initialize(self, P_m0: float, omega_pu: float, *args, **kwargs) -> np.ndarray:
        pass

    @abstractmethod
    def get_derivatives(self, local_state: np.ndarray, omega_pu: float, P_ref: float, *args, **kwargs) -> np.ndarray:
        pass
        
    @abstractmethod
    def get_Pm(self, local_state: np.ndarray) -> float:
        pass

class PSS(ABC):
    """ Base template for all types of Power System Stabilizers (PSS) """
    def __init__(self, name: str):
        self.name: str = name
        self.n_states: int = 0  
        
    @abstractmethod
    def initialize(self, *args, **kwargs) -> np.ndarray:
        pass

    @abstractmethod
    def get_derivatives(self, local_state: np.ndarray, omega_pu: float, P_e: float, *args, **kwargs) -> np.ndarray:
        pass
        
    @abstractmethod
    def get_Vpss(self, local_state: np.ndarray, omega_pu: float, P_e: float) -> float:
        pass

# =====================================================================
# Concrete Standard Models
# =====================================================================

class SEXS(AVR):
    """
    Simplified Excitation System (IEEE Type 1 Equivalent)
    Model with 1 State Variable: [E_fd]
    """
    def __init__(self, name: str, Ka: float = 200.0, Ta: float = 0.02, 
                 Efd_min: float = -5.0, Efd_max: float = 5.0):
        super().__init__(name)
        self.n_states = 1
        
        self.Ka = Ka
        self.Ta = Ta
        self.Efd_min = Efd_min
        self.Efd_max = Efd_max
        
    def initialize(self, V_ref: float, V_t: float, E_fd0: float, *args, **kwargs) -> np.ndarray:
        return np.array([E_fd0])

    def get_derivatives(self, local_state: np.ndarray, V_ref: float, V_t: float, V_pss: float, *args, **kwargs) -> np.ndarray:
        E_fd = local_state[0]
        error = V_ref - V_t + V_pss
        deriv = (self.Ka * error - E_fd) / self.Ta
        
        # Non-windup limiter
        if E_fd >= self.Efd_max and deriv > 0:
            deriv = 0.0
        elif E_fd <= self.Efd_min and deriv < 0:
            deriv = 0.0
            
        return np.array([deriv])
        
    def get_Efd(self, local_state: np.ndarray) -> float:
        return local_state[0]


class TGOV1(Governor):
    """
    IEEE Standard Steam Turbine-Governor Model (TGOV1)
    Model with 2 State Variables: [P_gv, x_2]
    - P_gv: Valve position
    - x_2: Internal turbine steam chest state
    """
    def __init__(self, name: str, R: float = 0.05, T1: float = 0.5, 
                 T2: float = 1.0, T3: float = 3.0, 
                 Vmax: float = 1.0, Vmin: float = 0.0, Dt: float = 0.0):
        super().__init__(name)
        self.n_states = 2
        
        self.R = R       # Speed droop (pu)
        self.T1 = T1     # Governor time constant (s)
        self.T2 = T2     # Numerator time constant (s)
        self.T3 = T3     # Denominator/Turbine time constant (s)
        self.Vmax = Vmax # Maximum valve position (pu)
        self.Vmin = Vmin # Minimum valve position (pu)
        self.Dt = Dt     # Turbine damping coefficient
        
    def initialize(self, P_m0: float, omega_pu: float, *args, **kwargs) -> np.ndarray:
        """ At steady-state, valve position and internal state equal the initial mechanical power. """
        return np.array([P_m0, P_m0])

    def get_derivatives(self, local_state: np.ndarray, omega_pu: float, P_ref_setpoint: float, *args, **kwargs) -> np.ndarray:
        P_gv = local_state[0]
        x_2 = local_state[1]
        
        delta_w = omega_pu 
        
        # 1. Valve logic: Clip the input request FIRST
        P_in = P_ref_setpoint - (delta_w / self.R)
        P_in = np.clip(P_in, self.Vmin, self.Vmax)
        
        # Calculate valve derivative
        dP_gv = (P_in - P_gv) / self.T1
        
        # Hard limits on state to prevent integration drift
        if P_gv >= self.Vmax and dP_gv > 0: dP_gv = 0.0
        if P_gv <= self.Vmin and dP_gv < 0: dP_gv = 0.0
            
        # 2. Steam Chest Derivative
        # Prevent numerical instability if T3 is extremely small
        if self.T3 < 0.001:
            dx_2 = 0.0 # Acts instantly, handled in get_Pm
        else:
            dx_2 = (P_gv - x_2) / self.T3
        
        return np.array([dP_gv, dx_2])

    def get_Pm(self, local_state: np.ndarray) -> float:
        P_gv = local_state[0]
        x_2 = local_state[1]
        
        # If T3 is near zero, bypass the chest delay completely
        if self.T3 < 0.001:
            return P_gv
            
        P_m = x_2 + (self.T2 / self.T3) * (P_gv - x_2)
        
        # Mechanical power cannot be negative (turbines don't motor)
        return max(0.0, P_m)

class PSS1A(PSS):
    """
    IEEE Standard Power System Stabilizer (PSS1A)
    Model with 3 State Variables: [x_w, x_1, x_2]
    - x_w: Washout filter state
    - x_1: First lead-lag state
    - x_2: Second lead-lag state
    """
    def __init__(self, name: str, K_pss: float = 10.0, T_w: float = 10.0, 
                 T1: float = 0.05, T2: float = 0.02, 
                 T3: float = 0.05, T4: float = 0.02, 
                 V_max: float = 0.1, V_min: float = -0.1):
        super().__init__(name)
        self.n_states = 3
        
        self.K_pss = K_pss
        self.T_w = T_w
        self.T1, self.T2 = T1, T2
        self.T3, self.T4 = T3, T4
        self.V_max, self.V_min = V_max, V_min
        
    def initialize(self, *args, **kwargs) -> np.ndarray:
        """ At steady state, speed deviation is 0, so all internal states are 0. """
        return np.zeros(3)

    def get_derivatives(self, local_state: np.ndarray, omega_pu: float, P_e: float, *args, **kwargs) -> np.ndarray:
        x_w, x_1, x_2 = local_state
        delta_w = omega_pu  # Input is speed deviation
        
        # 1. Washout Filter Derivative
        dx_w = (delta_w - x_w) / self.T_w
        v_1 = self.K_pss * (delta_w - x_w)  # Washout algebraic output
        
        # 2. First Lead-Lag Derivative
        dx_1 = (v_1 - x_1) / self.T2
        v_2 = x_1 + (self.T1 / self.T2) * (v_1 - x_1)  # Lead-lag algebraic output
        
        # 3. Second Lead-Lag Derivative
        dx_2 = (v_2 - x_2) / self.T4
        
        return np.array([dx_w, dx_1, dx_2])
        
    def get_Vpss(self, local_state: np.ndarray, omega_pu: float, P_e: float) -> float:
        """ Calculates the final clamped V_pss signal """
        x_w, x_1, x_2 = local_state
        delta_w = omega_pu
        
        # Reconstruct the algebraic chain to get the final output
        v_1 = self.K_pss * (delta_w - x_w)
        v_2 = x_1 + (self.T1 / self.T2) * (v_1 - x_1)
        v_3 = x_2 + (self.T3 / self.T4) * (v_2 - x_2)
        
        # Apply output limits
        V_pss = np.clip(v_3, self.V_min, self.V_max)
        return float(V_pss)