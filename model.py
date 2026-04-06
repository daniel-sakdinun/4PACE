from abc import ABC, abstractmethod
import numpy as np

class BusComponent(ABC):
    def __init__(self, name: str, P: float = 0.0, Q: float = 0.0, R: float = 0.0, X: float = 0.0):
        self.name = name
        self.P = P
        self.Q = Q
        
        self.Z = complex(R, X)

    @property
    def S(self) -> complex:
        """Complex Power (S)"""
        return complex(self.P, self.Q)

    @abstractmethod
    def cost(self) -> float:
        pass
    
    @abstractmethod
    def incremental_cost(self) -> float:
        pass

class BranchComponent(ABC):
    def __init__(self, name: str, R: float, X: float, S_max: float | None = None):
        self.name = name
        self.R = R
        self.X = X
        self.S_max = S_max
        
        self.Z = complex(R, X)
        self.Y = 1 / self.Z if self.Z != 0 else 0j

class SynchronousMachine(BusComponent):
    def __init__(self, name: str, P: float = 0.0, Q: float = 0.0,
                 a: float = 0.0, b: float = 0.0, c: float= 0.0,
                 Pmin: float = 0.0, Pmax: float | None = None, 
                 Qmin: float | None = None, Qmax: float | None = None,
                 S_rated: float | None = None, pf: float = 0.85, mode: str = 'generator',
                 R: float = 0.0, X: float = 0.0):
        super().__init__(name, P, Q, R, X)
        self.a, self.b, self.c = a, b, c
        self.mode:str = mode
        
        if S_rated is not None:
            max_p = S_rated * pf
            max_q = S_rated * np.sin(np.acos(pf))
            
            if mode == 'generator':
                self.Pmax = max_p
                self.Pmin = Pmin if Pmin is not None else 0.0
            elif mode == 'motor':
                self.Pmax = Pmax if Pmax is not None else 0.0
                self.Pmin = -max_p
            elif mode == 'condenser':
                self.Pmax = 0.0
                self.Pmin = 0.0
            elif mode == 'pumped_storage': # เป็นได้ทั้งสองอย่าง
                self.Pmax = max_p
                self.Pmin = -max_p
                
            self.Qmax = max_q
            self.Qmin = Qmin if Qmin is not None else (-max_q * 0.3)
        else:
            self.Pmax = Pmax if Pmax is not None else float('inf')
            self.Pmin = Pmin if Pmin is not None else (0.0 if mode == 'generator' else float('-inf'))
            self.Qmax = Qmax if Qmax is not None else float('inf')
            self.Qmin = Qmin if Qmin is not None else float('-inf')
    
    def cost(self) -> float:
        if self.P <= 0:
            return 0.0
        return self.a + self.b*self.P + self.c*(self.P**2)
    
    def incremental_cost(self) -> float:
        if self.P <= 0:
            return 0.0
        return self.b + (2 * self.c * self.P)

class Load(BusComponent):
    def __init__(self, name: str, model: str ='P',
                 P: float = 0, Q: float = 0,
                 R: float = 0, X: float = 0):
        """
        model: 'Z' (constant Z, P ∝ V^2)
               'I' (constant I, P ∝ V), 
               'P' (constant PQ), 
        """
        
        super().__init__(name, -abs(P), -Q, R, X)
        self.model = model
        
        self.P_nom = self.P
        self.Q_nom = self.Q
        
    def cost(self) -> float:
        return 0.0

    def incremental_cost(self) -> float:
        return 0.0
        
    def update_voltage_dependence(self, V_mag: float, V_nom: float = 1.0):
        if self.model == 'Z':
            self.P = self.P_nom * (V_mag / V_nom)**2
            self.Q = self.Q_nom * (V_mag / V_nom)**2
        elif self.model == 'I':
            self.P = self.P_nom * (V_mag / V_nom)
            self.Q = self.Q_nom * (V_mag / V_nom)
        elif self.model == 'P':
            self.P = self.P_nom
            self.Q = self.Q_nom

class Inverter(BusComponent):
    def __init__(self, name: str, S_max: float, 
                 P: float = 0.0, Q: float = 0.0,
                 control_mode: str = 'grid_following',
                 source_type: str = 'solar',
                 R: float = 0.0, X: float = 0.0):
        super().__init__(name, P, Q, R, X)
        self.S_max = S_max
        self.control_mode = control_mode
        self.source_type = source_type
        
        self.Qmax = S_max
        self.Qmin = -S_max
        
        if source_type in ['solar', 'wind']:
            self.Pmax = S_max
            self.Pmin = 0.0
        elif source_type == 'bess':
            self.Pmax = S_max
            self.Pmin = -S_max

    def cost(self) -> float:
        return 0.0
    
    def incremental_cost(self) -> float:
        return 0.0

class TransmissionLine(BranchComponent):
    def __init__(self, name: str,
                 R: float, X: float, B_shunt:float = 0.0,
                 S_max: float | None = None, length_km: float = 1.0):
        super().__init__(name, R, X, S_max)
        self.B_shunt = B_shunt
        self.length_km = length_km

class Transformer(BranchComponent):
    def __init__(self, name: str,
                 R: float, X: float,
                 tap_ratio: float = 1.0, phase_shift: float = 0.0,
                 S_max: float | None = None):
        super().__init__(name, R, X, S_max)
        self.tap_ratio = tap_ratio
        self.phase_shift = phase_shift