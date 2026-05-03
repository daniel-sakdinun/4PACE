from abc import ABC, abstractmethod
import numpy as np
from typing import Optional

from fourpace.control import AVR, Governor, PSS

def get_daily_capex_factor(interest_rate: float, lifetime_years: int) -> float:
    """
    Calculates the Capital Recovery Factor (CRF) converted to a daily rate.
    Used for amortizing Capital Expenditure (CapEx) in microgrid investment planning.
    
    Parameters:
    -----------
    interest_rate : float
        The annual discount or interest rate (e.g., 0.05 for 5%).
    lifetime_years : int
        The expected operational lifetime of the asset in years.
        
    Returns:
    --------
    float
        The daily fraction of the total capital cost.
    """
    if lifetime_years <= 0 or interest_rate <= 0:
        return 0.0
    crf = (interest_rate * (1 + interest_rate)**lifetime_years) / ((1 + interest_rate)**lifetime_years - 1)
    return crf / 365.0

class BusComponent(ABC):
    """
    Abstract base class for all elements connected to a single bus (node) in the power system.
    (e.g., Generators, Loads, Shunts, Inverters).
    """
    def __init__(self, name: str, P: float = 0.0, Q: float = 0.0, R: float = 0.0, X: float = 0.0):
        self.name = name
        self.P = P  # Active power injection (p.u.)
        self.Q = Q  # Reactive power injection (p.u.)
        self.Z = complex(R, X)  # Internal impedance (p.u.)

    @property
    def S(self) -> complex:
        """Returns the apparent power (S = P + jQ) of the component."""
        return complex(self.P, self.Q)

    @abstractmethod
    def cost(self) -> float:
        """Calculates the total operating cost of the component."""
        pass
    
    @abstractmethod
    def incremental_cost(self) -> float:
        """Calculates the marginal/incremental cost (dC/dP) for economic dispatch."""
        pass

class BranchComponent(ABC):
    """
    Abstract base class for all elements connecting two buses in the power system.
    (e.g., Transmission Lines, Transformers, Series FACTS).
    """
    def __init__(self, name: str, R: float, X: float, S_max: float | None = None):
        self.name = name
        self.R = R
        self.X = X
        self.S_max = S_max  # Thermal rating / Apparent power limit (MVA)
        self.Z = complex(R, X)
        self.Y = 1 / self.Z if self.Z != 0 else 0j

class SynchronousMachine(BusComponent):
    """
    Represents a synchronous generator, motor, or synchronous condenser.
    Supports steady-state power flow, economic dispatch, and transient dynamic analysis.
    """
    def __init__(self, name: str, P: float = 0.0, Q: float = 0.0,
                 a: float = 0.0, b: float = 0.0, c: float= 0.0,
                 Pmin: float = 0.0, Pmax: Optional[float] = None, 
                 Qmin: Optional[float] = None, Qmax: Optional[float] = None,
                 S_rated: Optional[float] = None, pf: float = 0.85, mode: str = 'generator',
                 R: float = 0.0, X: float = 0.0,
                 # --- Dynamics & Fault Sequence Parameters ---
                 H: float = 5.0, Xd: float = 1.0, Xd_prime: float = 0.3, Xd_sub: float = 0.2,       
                 X2: float = 0.2, X0: float = 0.05, Td0_prime: float = 5.0,    
                 # --- Controllers (Composition) ---
                 avr: Optional[AVR] = None,
                 gov: Optional[Governor] = None,
                 pss: Optional[PSS] = None,
                 # --- CEP Parameters ---
                 is_candidate: bool = False, capex_per_mw: float = 0.0, 
                 max_build_mw: float = 0.0, lifetime_years: int = 20, interest_rate: float = 0.05):
        super().__init__(name, P, Q, R, X)
        self.a, self.b, self.c = a, b, c  # Cost curve coefficients: C(P) = a + bP + cP^2
        self.mode: str = mode
        
        # Dynamics Parameters for Transient Stability
        self.H: float = H                     # Inertia constant (MW.s/MVA)
        self.Xd: float = Xd                   # Synchronous reactance (p.u.)
        self.Xd_prime: float = Xd_prime       # Transient reactance (p.u.)
        self.Xd_sub: float = Xd_sub           # Subtransient reactance (p.u.) (Used in fault calcs)
        self.X2: float = X2                   # Negative sequence reactance (p.u.)
        self.X0: float = X0                   # Zero sequence reactance (p.u.)
        self.Td0_prime: float = Td0_prime     # d-axis transient open-circuit time constant (s)
        
        # Dynamic Controllers
        self.avr: AVR = avr
        self.gov: Governor = gov
        self.pss: PSS = pss
        
        # Power Limits Calculation based on mode and power factor
        if S_rated is not None:
            max_p = S_rated * pf
            max_q = S_rated * np.sin(np.acos(pf))
            if mode == 'generator': self.Pmax, self.Pmin = max_p, Pmin if Pmin is not None else 0.0
            elif mode == 'motor': self.Pmax, self.Pmin = Pmax if Pmax is not None else 0.0, -max_p
            elif mode == 'condenser': self.Pmax, self.Pmin = 0.0, 0.0
            elif mode == 'pumped_storage': self.Pmax, self.Pmin = max_p, -max_p
            self.Qmax, self.Qmin = max_q, Qmin if Qmin is not None else (-max_q * 0.3)
        else:
            self.Pmax, self.Pmin = Pmax if Pmax is not None else float('inf'), Pmin if Pmin is not None else (0.0 if mode == 'generator' else float('-inf'))
            self.Qmax, self.Qmin = Qmax if Qmax is not None else float('inf'), Qmin if Qmin is not None else float('-inf')

        # Capacity Expansion Planning (CEP) Config
        self.is_candidate, self.capex_per_mw, self.max_build_mw = is_candidate, capex_per_mw, max_build_mw
        self.daily_capex_factor = get_daily_capex_factor(interest_rate, lifetime_years)
        self.built_P_max = 0.0
    
    def cost(self) -> float:
        """Calculates the quadratic generation cost: C(P) = a + b*P + c*P^2."""
        if self.P <= 0: return 0.0
        return self.a + self.b*self.P + self.c*(self.P**2)
    
    def incremental_cost(self) -> float:
        """Calculates the marginal cost: dC/dP = b + 2*c*P."""
        if self.P <= 0: return 0.0
        return self.b + (2 * self.c * self.P)

class AsynchronousMachine(BusComponent):
    """
    Represents a three-phase induction motor.
    Power consumption dynamically updates based on voltage, slip, and torque-speed characteristics.
    """
    def __init__(self, name: str, P_rated: float, V_rated: float,
                 Rs: float, Xs: float, Rr: float, Xr: float, Xm: float,
                 poles: int = 4, freq: float = 50.0,
                 s: float = 0.02, load_type: str = 'constant_torque'):
        super().__init__(name, P=0.0, Q=0.0)
        self.Rs, self.Xs, self.Rr, self.Xr, self.Xm = Rs, Xs, Rr, Xr, Xm # Equivalent circuit parameters
        self.P_rated, self.P_rated_base = P_rated, P_rated
        self.V_rated, self.poles, self.freq = V_rated, poles, freq
        self.s, self.load_type = s, load_type # Operating slip and mechanical load model
    
    def cost(self) -> float: return 0.0
    def incremental_cost(self) -> float: return 0.0

    def update_pq_from_slip(self, V_mag_pu: float, Sbase:float):
        """Updates the active (P) and reactive (Q) power draw based on the current slip and terminal voltage."""
        V_phase = (V_mag_pu * self.V_rated) / np.sqrt(3)
        Z_rotor = (self.Rr / self.s) + 1j*self.Xr
        Z_parallel = (1j*self.Xm * Z_rotor) / (1j*self.Xm + Z_rotor)
        Z_total = (self.Rs + 1j*self.Xs) + Z_parallel
        I_s = V_phase / Z_total
        S_motor = 3 * V_phase * np.conj(I_s)
        
        # Power is drawn (consumed), hence the negative sign for standard injection models
        self.P = -S_motor.real / Sbase
        self.Q = -S_motor.imag / Sbase

class Load(BusComponent):
    """
    Represents a static electrical load.
    Supports ZIP modeling (Constant Impedance 'Z', Constant Current 'I', Constant Power 'P').
    """
    def __init__(self, name: str, model: str ='P', P: float = 0, Q: float = 0, R: float = 0, X: float = 0):
        super().__init__(name, -abs(P), -Q, R, X) # Loads consume power (negative injection)
        self.model = model
        self.P_nom, self.Q_nom = self.P, self.Q # Nominal power at 1.0 p.u. voltage
        
    def cost(self) -> float: return 0.0
    def incremental_cost(self) -> float: return 0.0
        
    def update_voltage_dependence(self, V_mag: float, V_nom: float = 1.0):
        """Modifies actual power consumption based on the terminal voltage and selected load model."""
        if self.model == 'Z':
            self.P = self.P_nom * (V_mag / V_nom)**2
            self.Q = self.Q_nom * (V_mag / V_nom)**2
        elif self.model == 'I':
            self.P = self.P_nom * (V_mag / V_nom)
            self.Q = self.Q_nom * (V_mag / V_nom)
        elif self.model == 'P':
            self.P, self.Q = self.P_nom, self.Q_nom

class Shunt(BusComponent):
    """
    Represents a passive shunt compensator (Capacitor bank or Reactor).
    Reactive power injection varies with the square of the voltage.
    """
    def __init__(self, name: str, Q_nom: float = 0.0, V_nom: float = 1.0):
        super().__init__(name, P=0.0, Q=Q_nom)
        self.Q_nom, self.Q_nom_base, self.V_nom = Q_nom, Q_nom, V_nom

    def update_voltage_dependence(self, V_mag: float):
        """Updates reactive power injection based on actual bus voltage."""
        self.Q = self.Q_nom * (V_mag / self.V_nom)**2

    def cost(self) -> float: return 0.0
    def incremental_cost(self) -> float: return 0.0

class Inverter(BusComponent):
    """
    Represents an Inverter-Based Resource (IBR) such as a Solar PV farm or a BESS inverter.
    Models grid-following behavior and implements fault current limitations during short-circuits.
    """
    def __init__(self, name: str, S_max: float, P: float = 0.0, Q: float = 0.0,
                 control_mode: str = 'grid_following', source_type: str = 'solar',
                 R: float = 0.0, X: float = 0.0,
                 # --- Fault Sequence Parameters ---
                 I_fault_limit_pu: float = 1.2, # Inverters act as limited current sources during faults
                 # --- CEP Parameters ---
                 is_candidate: bool = False, capex_per_mw: float = 0.0, 
                 max_build_mw: float = 0.0, lifetime_years: int = 15, interest_rate: float = 0.05):
        super().__init__(name, P, Q, R, X)
        self.S_max, self.S_max_base = S_max, S_max
        self.control_mode, self.source_type = control_mode, source_type
        
        self.I_fault_limit_pu = I_fault_limit_pu # Max current injection capability during severe voltage sags
        
        # Capability Curve Limits
        self.Qmax, self.Qmin = S_max, -S_max
        if source_type in ['solar', 'wind']:
            self.Pmax, self.Pmin = S_max, 0.0
        elif source_type == 'bess':
            self.Pmax, self.Pmin = S_max, -S_max

        # Capacity Expansion Planning (CEP) Config
        self.is_candidate = is_candidate
        self.capex_per_mw = capex_per_mw
        self.max_build_mw = max_build_mw
        self.daily_capex_factor = get_daily_capex_factor(interest_rate, lifetime_years)
        self.built_S_max = 0.0 

    def cost(self) -> float: return 0.0
    def incremental_cost(self) -> float: return 0.0
    
class Battery(BusComponent):
    """
    Represents a Battery Energy Storage System (BESS).
    Manages State of Charge (SoC) across time-series simulations and optimization routines.
    """
    def __init__(self, name: str, P_max: float, E_max: float, init_soc: float = 0.5, eta: float = 0.95,
                 # --- Fault Sequence Parameters ---
                 I_fault_limit_pu: float = 1.5,
                 # --- CEP Parameters ---
                 is_candidate: bool = False, capex_per_mw: float = 0.0, capex_per_mwh: float = 0.0, 
                 max_build_mw: float = 0.0, max_build_mwh: float = 0.0, 
                 lifetime_years: int = 10, interest_rate: float = 0.05):
        super().__init__(name)
        self.P_max, self.E_max = P_max, E_max           # Power (MW) and Energy (MWh) capacities
        self.init_soc, self.eta = init_soc, eta         # Initial State of Charge and Round-trip efficiency
        self.I_fault_limit_pu = I_fault_limit_pu
        
        self.P_ch_series, self.P_dis_series, self.SoC_series = [], [], []
        self.P, self.Q, self.SoC = 0.0, 0.0, init_soc

        # Capacity Expansion Planning (CEP) Config (Co-optimizing Power & Energy)
        self.is_candidate = is_candidate
        self.capex_per_mw = capex_per_mw
        self.capex_per_mwh = capex_per_mwh
        self.max_build_mw = max_build_mw
        self.max_build_mwh = max_build_mwh
        self.daily_capex_factor = get_daily_capex_factor(interest_rate, lifetime_years)
        self.built_P_max = 0.0
        self.built_E_max = 0.0

    def cost(self) -> float: return 0.0
    def incremental_cost(self) -> float: return 0.0

class TransmissionLine(BranchComponent):
    """
    Represents a three-phase AC transmission line using the Pi-equivalent circuit model.
    Includes zero-sequence impedances crucial for asymmetrical fault and open-conductor analysis.
    """
    def __init__(self, name: str, R: float, X: float, B_shunt:float = 0.0,
                 S_max: float | None = None, length_km: float = 1.0,
                 # --- Fault Sequence Parameters ---
                 R0: float | None = None, X0: float | None = None, B0_shunt: float | None = None,
                 # --- CEP Parameters ---
                 is_candidate: bool = False, capex_per_mva: float = 0.0, 
                 max_build_mva: float = 0.0, lifetime_years: int = 40, interest_rate: float = 0.05):
        super().__init__(name, R, X, S_max)
        self.B_shunt = B_shunt      # Positive sequence total line charging susceptance
        self.length_km = length_km
        
        # Fault Config (Zero Sequence)
        # If not provided, roughly estimate Zero Sequence impedance as 3x Positive Sequence
        self.R0 = R0 if R0 is not None else 3 * R
        self.X0 = X0 if X0 is not None else 3 * X
        self.B0_shunt = B0_shunt if B0_shunt is not None else 0.0
        self.Z0 = complex(self.R0, self.X0)
        self.Y0 = 1 / self.Z0 if self.Z0 != 0 else 0j
        
        # Capacity Expansion Planning (CEP) Config
        self.is_candidate = is_candidate
        self.capex_per_mva = capex_per_mva
        self.max_build_mva = max_build_mva
        self.daily_capex_factor = get_daily_capex_factor(interest_rate, lifetime_years)
        self.built_S_max = 0.0

class Transformer(BranchComponent):
    """
    Represents a two-winding power transformer.
    Supports tap changing (voltage regulation), phase shifting, and specific grounding configurations 
    (connection types) which dictate the zero-sequence path during unbalanced faults.
    """
    def __init__(self, name: str, R: float, X: float,
                 tap_ratio: float = 1.0, phase_shift: float = 0.0, S_max: float | None = None,
                 auto_tap: bool = False, controlled_bus: str | None = None,
                 target_V: float = 1.0, tap_step: float = 0.0125,
                 tap_min: float = 0.90, tap_max: float = 1.10,
                 # --- Fault Sequence Parameters ---
                 R0: float | None = None, X0: float | None = None, 
                 connection_type: str = 'Yg-Yg', # Critical for determining the Zero-Sequence network topology
                 # --- CEP Parameters ---
                 is_candidate: bool = False, capex_per_mva: float = 0.0, 
                 max_build_mva: float = 0.0, lifetime_years: int = 30, interest_rate: float = 0.05):
        super().__init__(name, R, X, S_max)
        self.tap_ratio, self.phase_shift = tap_ratio, phase_shift
        self.auto_tap, self.controlled_bus = auto_tap, controlled_bus
        self.target_V, self.tap_step = target_V, tap_step
        self.tap_min, self.tap_max = tap_min, tap_max
        
        # Fault Config (Zero Sequence & Grounding)
        # Options typically include: "Yg-Yg", "Delta-Yg", "Yg-Delta", "Delta-Delta"
        self.connection_type = connection_type
        # For transformers, Z0 is typically equal to Z1 (Positive Sequence) unless specified otherwise
        self.R0 = R0 if R0 is not None else R
        self.X0 = X0 if X0 is not None else X
        self.Z0 = complex(self.R0, self.X0)
        self.Y0 = 1 / self.Z0 if self.Z0 != 0 else 0j
        
        # Capacity Expansion Planning (CEP) Config
        self.is_candidate = is_candidate
        self.capex_per_mva = capex_per_mva
        self.max_build_mva = max_build_mva
        self.daily_capex_factor = get_daily_capex_factor(interest_rate, lifetime_years)
        self.built_S_max = 0.0