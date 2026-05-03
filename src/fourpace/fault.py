import numpy as np
import pandas as pd

from fourpace.psys import Grid

def ThreePhaseFault(grid: Grid, fault_bus: str, Z_fault: complex = 0j, max_iter: int = 15, tol: float = 1e-4, verbose: int = 0) -> dict:
    """ 
    Iterative Solver for Symmetrical Three-Phase (3PH) Fault. 
    
    Calculates fault currents and voltage profiles during a balanced three-phase short circuit.
    This advanced iterative solver actively monitors bus voltages and allows Inverter-Based 
    Resources (IBRs) to dynamically limit and inject reactive fault current.
    
    Parameters:
    -----------
    grid : Grid
        The power system grid object.
    fault_bus : str
        The name of the bus where the fault occurs.
    Z_fault : complex
        Fault impedance. Default is 0j (bolted fault/dead short).
    max_iter : int
        Maximum number of iterations for the IBR current injection loop.
    tol : float
        Convergence tolerance for the voltage profile differences.
    verbose : int
        Verbosity level for logging.
        
    Returns:
    --------
    dict
        A dictionary containing total fault current (p.u.), Thevenin impedance magnitude, 
        and the full system voltage profile (p.u.).
    """
    if verbose > 0: 
        print(f"\n⚡ Initiating Iterative 3-Phase Fault Analysis (IBR-Enabled) at Bus '{fault_bus}'...")
        
    nodes_list = list(grid.nodes)
    k = nodes_list.index(fault_bus)
    num_bus = len(grid.buses)
    
    # 1. Build Positive Sequence Y-bus and inject the Fault Admittance
    Ybus_1 = grid.build_ybus_pos()
    Y_fault = Ybus_1.copy()
    
    y_f = 1 / Z_fault if Z_fault != 0 else 1e6  # Use a large admittance for a dead short
    Y_fault[k, k] += y_f
    
    # Pre-invert the faulted Y-bus matrix to accelerate the iterative loop
    try:
        Z_fault_mat = np.linalg.inv(Y_fault)
    except np.linalg.LinAlgError:
        Y_fault += np.eye(num_bus) * 1e-6j # Add artificial shunt to prevent singularity
        Z_fault_mat = np.linalg.inv(Y_fault)

    Z_th = Z_fault_mat[k, k] 

    # 2. Prepare Initial Norton Current Injections
    I_sync = np.zeros(num_bus, dtype=complex)
    ibrs = []
    
    V_pre = 1.0 + 0j # Assume a flat start 1.0 pu pre-fault voltage profile
    
    for i, bus in enumerate(grid.buses):
        for comp in bus.components:
            c_name = type(comp).__name__
            # A. Synchronous Generators -> Modeled as a Constant Norton Equivalent Injection
            if c_name == 'SynchronousMachine':
                Xd_sub = getattr(comp, 'Xd_sub', 0.2)
                I_sync[i] += V_pre / (1j * Xd_sub)
                
            # B. Gather IBRs (Solar & Battery) -> Registered for Dynamic Iterative Injection
            elif c_name in ['Inverter', 'Battery']:
                limit = getattr(comp, 'I_fault_limit_pu', 1.2)
                S_rated = getattr(comp, 'S_max', getattr(comp, 'P_max', 0.0)) / grid.Sbase
                I_max = limit * S_rated 
                ibrs.append({'bus_idx': i, 'I_max': I_max, 'name': comp.name, 'I_inj': 0j})

    # 3. The Iterative Voltage & Current Solver Loop
    V_profile = np.ones(num_bus, dtype=complex) * V_pre
    converged = False
    
    for iteration in range(max_iter):
        I_total_inj = I_sync.copy()
        
        # Recalculate dynamic IBR current injection based on the latest estimated voltage
        for ibr in ibrs:
            idx = ibr['bus_idx']
            v_mag = abs(V_profile[idx])
            v_angle = np.angle(V_profile[idx])
            
            # IBR Fault Ride-Through Logic: 
            # If voltage drops below 0.9 p.u., inject purely reactive current to support the grid
            if v_mag < 0.9:
                I_inj_mag = ibr['I_max']
                ibr['I_inj'] = I_inj_mag * np.exp(1j * (v_angle + np.pi/2))
            else:
                ibr['I_inj'] = 0j # Normal operation (simplified; ignores pre-fault load current)
                
            I_total_inj[idx] += ibr['I_inj']

        # Solve for new grid voltages: [V] = [Z_fault] * [I_total]
        V_new = Z_fault_mat @ I_total_inj
        
        # Check Convergence Threshold
        v_diff = np.max(np.abs(V_new - V_profile))
        V_profile = V_new
        
        if v_diff < tol:
            converged = True
            if verbose > 1:
                print(f"   🔄 Iterative Fault Solver converged in {iteration+1} iterations.")
            break

    if not converged and verbose > 0:
        print("   ⚠️ Iterative Fault Solver reached max iterations without strict convergence.")

    # 4. Extract and Format Results
    # Total fault current is the current flowing from the faulted bus to ground
    I_fault = V_profile[k] * y_f

    return {
        'I_fault_pu': abs(I_fault), 
        'Z_th_mag': abs(Z_th), 
        'V_profile': [abs(v) for v in V_profile]
    }


def LineToGroundFault(grid: Grid, fault_bus: str, Z_fault: complex = 0j, max_iter: int = 15, tol: float = 1e-4, verbose: int = 0) -> dict:
    """ 
    Iterative Solver for Single Line-to-Ground (SLG) Fault. 
    
    Simulates an asymmetrical SLG fault by connecting the Positive, Negative, and Zero sequence 
    networks in series. Smart IBRs are constrained to inject balanced reactive current exclusively 
    into the Positive Sequence network.
    """
    if verbose > 0: 
        print(f"\n⚡ Initiating Iterative SLG Fault Analysis (IBR-Enabled) at Bus '{fault_bus}'...")
        
    nodes_list = list(grid.nodes)
    k = nodes_list.index(fault_bus)
    num_bus = len(grid.buses)
    
    # 1. Build Sequence Y-bus Matrices
    Ybus_1 = grid.build_ybus_pos()
    Ybus_0 = grid.build_ybus_zero()
    
    # Pre-invert to derive Z-bus matrices (crucial for fetching Thevenin equivalents)
    try:
        Zbus_1 = np.linalg.inv(Ybus_1)
    except np.linalg.LinAlgError:
        Ybus_1 += np.eye(num_bus) * 1e-6j
        Zbus_1 = np.linalg.inv(Ybus_1)
        
    try:
        Zbus_0 = np.linalg.inv(Ybus_0)
    except np.linalg.LinAlgError:
        Ybus_0 += np.eye(num_bus) * 1e-6j
        Zbus_0 = np.linalg.inv(Ybus_0)
        
    Z1_kk = Z2_kk = Zbus_1[k, k]  # Assumption: Negative Sequence (Z2) closely mirrors Positive Sequence (Z1)
    Z0_kk = Zbus_0[k, k]
    
    # 2. Prepare Norton Current Injections (Applied to Positive Sequence ONLY)
    I_sync = np.zeros(num_bus, dtype=complex)
    ibrs = []
    V_pre = 1.0 + 0j 
    
    for i, bus in enumerate(grid.buses):
        for comp in bus.components:
            c_name = type(comp).__name__
            if c_name == 'SynchronousMachine':
                Xd_sub = getattr(comp, 'Xd_sub', 0.2)
                I_sync[i] += V_pre / (1j * Xd_sub)
                
            elif c_name in ['Inverter', 'Battery']:
                limit = getattr(comp, 'I_fault_limit_pu', 1.2)
                S_rated = getattr(comp, 'S_max', getattr(comp, 'P_max', 0.0)) / grid.Sbase
                I_max = limit * S_rated 
                ibrs.append({'bus_idx': i, 'I_max': I_max, 'name': comp.name, 'I_inj': 0j})

    # 3. Iterative Solver Loop
    V1_profile = np.ones(num_bus, dtype=complex) * V_pre
    V2_new = np.zeros(num_bus, dtype=complex)
    V0_new = np.zeros(num_bus, dtype=complex)
    converged = False
    
    for iteration in range(max_iter):
        I_1_inj = I_sync.copy()
        
        # Update IBR current injection dynamically based on the Positive Sequence Voltage
        for ibr in ibrs:
            idx = ibr['bus_idx']
            v1_mag = abs(V1_profile[idx])
            v1_angle = np.angle(V1_profile[idx])
            
            if v1_mag < 0.9:
                ibr['I_inj'] = ibr['I_max'] * np.exp(1j * (v1_angle + np.pi/2))
            else:
                ibr['I_inj'] = 0j
                
            I_1_inj[idx] += ibr['I_inj']

        # Calculate pre-fault (Thevenin) Positive Sequence Voltage due to all network injections
        E_th = Zbus_1 @ I_1_inj
        E_th_k = E_th[k]
        
        # Solve the series-connected sequence networks at the faulted bus
        I_f1 = E_th_k / (Z1_kk + Z2_kk + Z0_kk + (3 * Z_fault))
        
        # Re-calculate specific sequence voltages at all system buses
        V1_new = E_th - (Zbus_1[:, k] * I_f1)
        V2_new = - (Zbus_1[:, k] * I_f1) # Because Zbus_2 = Zbus_1 and I_f2 = I_f1
        V0_new = - (Zbus_0[:, k] * I_f1) # Because I_f0 = I_f1
        
        # Convergence is strictly evaluated based on the Positive Sequence Voltage stability
        v_diff = np.max(np.abs(V1_new - V1_profile))
        V1_profile = V1_new
        
        if v_diff < tol:
            converged = True
            if verbose > 1:
                print(f"   🔄 Iterative SLG Solver converged in {iteration+1} iterations.")
            break

    if not converged and verbose > 0:
        print("   ⚠️ Iterative SLG Solver reached max iterations without strict convergence.")

    # 4. Compile Physical Results
    I_fault_total = 3 * I_f1
    Va_profile = V1_profile + V2_new + V0_new # Reconstruct Total Phase A voltage

    return {
        'I_fault_pu': abs(I_fault_total), 
        'Z_th_mag': abs(Z1_kk + Z2_kk + Z0_kk) / 3, 
        'V_profile': [abs(v) for v in Va_profile]
    }


def LineToLineFault(grid: Grid, fault_bus: str, Z_fault: complex = 0j, max_iter: int = 15, tol: float = 1e-4, verbose: int = 0) -> dict:
    """ 
    Iterative Solver for Line-to-Line (LL) Fault. 
    
    Simulates a short circuit between two phases. Connects the Positive and Negative sequence 
    networks in parallel. The Zero Sequence network is mathematically excluded from this analysis.
    """
    if verbose > 0: 
        print(f"\n⚡ Initiating Iterative LL Fault Analysis (IBR-Enabled) at Bus '{fault_bus}'...")
        
    nodes_list = list(grid.nodes)
    k = nodes_list.index(fault_bus)
    num_bus = len(grid.buses)
    
    # 1. Build Sequence Y-bus Matrices (Zero Sequence is completely ignored in LL faults)
    Ybus_1 = grid.build_ybus_pos()
    
    # Pre-invert to get Z-bus matrix
    try:
        Zbus_1 = np.linalg.inv(Ybus_1)
    except np.linalg.LinAlgError:
        Ybus_1 += np.eye(num_bus) * 1e-6j
        Zbus_1 = np.linalg.inv(Ybus_1)
        
    Z1_kk = Z2_kk = Zbus_1[k, k]  
    
    # 2. Prepare Norton Current Injections (Positive Sequence ONLY)
    I_sync = np.zeros(num_bus, dtype=complex)
    ibrs = []
    V_pre = 1.0 + 0j 
    
    for i, bus in enumerate(grid.buses):
        for comp in bus.components:
            c_name = type(comp).__name__
            if c_name == 'SynchronousMachine':
                Xd_sub = getattr(comp, 'Xd_sub', 0.2)
                I_sync[i] += V_pre / (1j * Xd_sub)
                
            elif c_name in ['Inverter', 'Battery']:
                limit = getattr(comp, 'I_fault_limit_pu', 1.2)
                S_rated = getattr(comp, 'S_max', getattr(comp, 'P_max', 0.0)) / grid.Sbase
                I_max = limit * S_rated 
                ibrs.append({'bus_idx': i, 'I_max': I_max, 'name': comp.name, 'I_inj': 0j})

    # 3. Iterative Solver Loop
    V1_profile = np.ones(num_bus, dtype=complex) * V_pre
    V2_new = np.zeros(num_bus, dtype=complex)
    converged = False
    
    for iteration in range(max_iter):
        I_1_inj = I_sync.copy()
        
        # Update IBR reactive current injection limits
        for ibr in ibrs:
            idx = ibr['bus_idx']
            v1_mag = abs(V1_profile[idx])
            v1_angle = np.angle(V1_profile[idx])
            
            if v1_mag < 0.9:
                ibr['I_inj'] = ibr['I_max'] * np.exp(1j * (v1_angle + np.pi/2))
            else:
                ibr['I_inj'] = 0j
                
            I_1_inj[idx] += ibr['I_inj']

        E_th = Zbus_1 @ I_1_inj
        E_th_k = E_th[k]
        
        # Solve the parallel-connected sequence networks at the fault bus
        I_f1 = E_th_k / (Z1_kk + Z2_kk + Z_fault)
        
        V1_new = E_th - (Zbus_1[:, k] * I_f1)
        
        # In an LL fault, I_f2 = -I_f1
        # Therefore, V_2 = - Z_2 * I_f2 = - Z_1 * (-I_f1) = Z_1 * I_f1
        V2_new = Zbus_1[:, k] * I_f1 
        
        v_diff = np.max(np.abs(V1_new - V1_profile))
        V1_profile = V1_new
        
        if v_diff < tol:
            converged = True
            if verbose > 1:
                print(f"   🔄 Iterative LL Solver converged in {iteration+1} iterations.")
            break

    if not converged and verbose > 0:
        print("   ⚠️ Iterative LL Solver reached max iterations without strict convergence.")

    # 4. Compile Results
    # Total physical fault current between the two shorted lines is sqrt(3) * I_f1
    I_fault_total = np.sqrt(3) * I_f1 
    
    # Phase A is the unfaulted phase: V_a = V_1 + V_2 + V_0 (where V_0 is exactly 0)
    Va_profile = V1_profile + V2_new 

    return {
        'I_fault_pu': abs(I_fault_total), 
        'Z_th_mag': abs(Z1_kk + Z2_kk), 
        'V_profile': [abs(v) for v in Va_profile]
    }


def DoubleLineToGroundFault(grid: Grid, fault_bus: str, Z_fault: complex = 0j, max_iter: int = 15, tol: float = 1e-4, verbose: int = 0) -> dict:
    """ 
    Iterative Solver for Double Line-to-Ground (DLG) Fault. 
    
    Simulates a severe asymmetrical fault by connecting all three sequence networks (Positive, 
    Negative, and Zero) in a parallel configuration.
    """
    if verbose > 0: 
        print(f"\n⚡ Initiating Iterative DLG Fault Analysis (IBR-Enabled) at Bus '{fault_bus}'...")
        
    nodes_list = list(grid.nodes)
    k = nodes_list.index(fault_bus)
    num_bus = len(grid.buses)
    
    # 1. Build Sequence Y-bus Matrices
    Ybus_1 = grid.build_ybus_pos()
    Ybus_0 = grid.build_ybus_zero()
    
    try:
        Zbus_1 = np.linalg.inv(Ybus_1)
    except np.linalg.LinAlgError:
        Ybus_1 += np.eye(num_bus) * 1e-6j
        Zbus_1 = np.linalg.inv(Ybus_1)
        
    try:
        Zbus_0 = np.linalg.inv(Ybus_0)
    except np.linalg.LinAlgError:
        Ybus_0 += np.eye(num_bus) * 1e-6j
        Zbus_0 = np.linalg.inv(Ybus_0)
        
    Z1_kk = Z2_kk = Zbus_1[k, k]  
    Z0_kk = Zbus_0[k, k]
    
    # Calculate the parallel equivalent of Negative and Zero Sequence networks + fault impedance
    Z0_total = Z0_kk + (3 * Z_fault)
    Z_p = (Z2_kk * Z0_total) / (Z2_kk + Z0_total) if (Z2_kk + Z0_total) != 0 else 1e6
    
    # 2. Prepare Norton Current Injections
    I_sync = np.zeros(num_bus, dtype=complex)
    ibrs = []
    V_pre = 1.0 + 0j 
    
    for i, bus in enumerate(grid.buses):
        for comp in bus.components:
            c_name = type(comp).__name__
            if c_name == 'SynchronousMachine':
                Xd_sub = getattr(comp, 'Xd_sub', 0.2)
                I_sync[i] += V_pre / (1j * Xd_sub)
                
            elif c_name in ['Inverter', 'Battery']:
                limit = getattr(comp, 'I_fault_limit_pu', 1.2)
                S_rated = getattr(comp, 'S_max', getattr(comp, 'P_max', 0.0)) / grid.Sbase
                I_max = limit * S_rated 
                ibrs.append({'bus_idx': i, 'I_max': I_max, 'name': comp.name, 'I_inj': 0j})

    # 3. Iterative Solver Loop
    V1_profile = np.ones(num_bus, dtype=complex) * V_pre
    V2_new = np.zeros(num_bus, dtype=complex)
    V0_new = np.zeros(num_bus, dtype=complex)
    I_f0 = 0j
    converged = False
    
    for iteration in range(max_iter):
        I_1_inj = I_sync.copy()
        
        for ibr in ibrs:
            idx = ibr['bus_idx']
            v1_mag = abs(V1_profile[idx])
            v1_angle = np.angle(V1_profile[idx])
            
            if v1_mag < 0.9:
                ibr['I_inj'] = ibr['I_max'] * np.exp(1j * (v1_angle + np.pi/2))
            else:
                ibr['I_inj'] = 0j
                
            I_1_inj[idx] += ibr['I_inj']

        E_th = Zbus_1 @ I_1_inj
        E_th_k = E_th[k]
        
        # Solve the positive sequence fault current
        I_f1 = E_th_k / (Z1_kk + Z_p)
        
        # Current division to find Negative and Zero sequence fault currents
        # Negative sign indicates they flow OUT of the fault point into the passive networks
        I_f2 = -I_f1 * (Z0_total / (Z2_kk + Z0_total))
        I_f0 = -I_f1 * (Z2_kk / (Z2_kk + Z0_total))
        
        # Calculate sequence voltages at all buses
        V1_new = E_th - (Zbus_1[:, k] * I_f1)
        V2_new = - (Zbus_1[:, k] * I_f2)  
        V0_new = - (Zbus_0[:, k] * I_f0)
        
        v_diff = np.max(np.abs(V1_new - V1_profile))
        V1_profile = V1_new
        
        if v_diff < tol:
            converged = True
            if verbose > 1:
                print(f"   🔄 Iterative DLG Solver converged in {iteration+1} iterations.")
            break

    if not converged and verbose > 0:
        print("   ⚠️ Iterative DLG Solver reached max iterations without strict convergence.")

    # 4. Compile Results
    # Total physical fault current flowing into the ground return path is 3 * I_f0
    I_fault_total = 3 * I_f0
    
    Va_profile = V1_profile + V2_new + V0_new 

    return {
        'I_fault_pu': abs(I_fault_total), 
        'Z_th_mag': abs(Z1_kk + Z_p), 
        'V_profile': [abs(v) for v in Va_profile]
    }


def OpenConductorFault(grid: Grid, from_bus: str, to_bus: str, n: int = 1, max_iter: int = 15, tol: float = 1e-4, verbose: int = 0) -> dict:
    """
    Iterative Solver for Open Conductor (Series) Faults.
    
    Evaluates asymmetric flow configurations caused by broken conductors spanning a transmission line.
    
    Parameters:
    -----------
    n : int
        Number of open conductors. 
        - n=1: One Conductor Open (Yields a parallel sequence network structure like LL faults).
        - n=2: Two Conductors Open (Yields a series sequence network structure like SLG faults).
    """
    if n not in [1, 2]:
        raise ValueError(f"Invalid number of open conductors (n). Must be 1 or 2. Got: {n}")
        
    if verbose > 0: 
        print(f"\n⚡ Initiating Iterative {n}-Open Conductor Fault Analysis (IBR-Enabled) between Bus '{from_bus}' and '{to_bus}'...")
    
    nodes_list = list(grid.nodes)
    if from_bus not in nodes_list or to_bus not in nodes_list:
        raise ValueError(f"❌ Bus '{from_bus}' or '{to_bus}' not found in the grid!")
    
    faulted_branch = None
    for u, v, data in grid.edges(data=True):
        if (u == from_bus and v == to_bus) or (u == to_bus and v == from_bus):
            faulted_branch = data.get('obj')
            break
            
    if faulted_branch is None:
        raise ValueError(f"❌ No branch found connecting Bus '{from_bus}' and '{to_bus}'!")
    
    m = nodes_list.index(from_bus)
    k = nodes_list.index(to_bus)
    num_bus = len(grid.buses)
    
    # 1. Build Y-bus Matrices
    Ybus_1 = grid.build_ybus_pos()
    Ybus_0 = grid.build_ybus_zero()
    
    try:
        Zbus_1 = np.linalg.inv(Ybus_1)
    except np.linalg.LinAlgError:
        Ybus_1 += np.eye(num_bus) * 1e-6j
        Zbus_1 = np.linalg.inv(Ybus_1)
        
    try:
        Zbus_0 = np.linalg.inv(Ybus_0)
    except np.linalg.LinAlgError:
        Ybus_0 += np.eye(num_bus) * 1e-6j
        Zbus_0 = np.linalg.inv(Ybus_0)

    # Calculate Equivalent Thevenin Impedance LOOKING INTO the open port from the grid
    Z1_th = Zbus_1[m, m] + Zbus_1[k, k] - (2 * Zbus_1[m, k])
    Z2_th = Z1_th  # Assumption: Z2 approx equals Z1
    Z0_th = Zbus_0[m, m] + Zbus_0[k, k] - (2 * Zbus_0[m, k])
    
    # 2. Prepare Norton Current Injections
    I_sync = np.zeros(num_bus, dtype=complex)
    ibrs = []
    V_pre = 1.0 + 0j 
    
    for i, bus in enumerate(grid.buses):
        for comp in bus.components:
            c_name = type(comp).__name__
            if c_name == 'SynchronousMachine':
                Xd_sub = getattr(comp, 'Xd_sub', 0.2)
                I_sync[i] += V_pre / (1j * Xd_sub)
                
            elif c_name in ['Inverter', 'Battery']:
                limit = getattr(comp, 'I_fault_limit_pu', 1.2)
                S_rated = getattr(comp, 'S_max', getattr(comp, 'P_max', 0.0)) / grid.Sbase
                I_max = limit * S_rated 
                ibrs.append({'bus_idx': i, 'I_max': I_max, 'name': comp.name, 'I_inj': 0j})

    # 3. Iterative Solver Loop
    V1_profile = np.ones(num_bus, dtype=complex) * V_pre
    V2_new = np.zeros(num_bus, dtype=complex)
    V0_new = np.zeros(num_bus, dtype=complex)
    I_fault_total = 0j
    converged = False
    
    for iteration in range(max_iter):
        I_1_inj = I_sync.copy()
        
        for ibr in ibrs:
            idx = ibr['bus_idx']
            v1_mag = abs(V1_profile[idx])
            v1_angle = np.angle(V1_profile[idx])
            
            if v1_mag < 0.9:
                ibr['I_inj'] = ibr['I_max'] * np.exp(1j * (v1_angle + np.pi/2))
            else:
                ibr['I_inj'] = 0j
                
            I_1_inj[idx] += ibr['I_inj']

        # Calculate pre-fault voltage profile due to ALL grid injections
        V1_unf = Zbus_1 @ I_1_inj
        
        # Driving voltage for the series fault 
        E_th_mk = V_pre 
        
        # Calculate Sequence Currents actively passing through the open point
        if n == 1:
            Z_eq = (Z2_th * Z0_th) / (Z2_th + Z0_th) if (Z2_th + Z0_th) != 0 else 1e6
            I_sequence = E_th_mk / (Z1_th + Z_eq)
            
            I_a1 = I_sequence
            I_a2 = -I_sequence * (Z0_th / (Z2_th + Z0_th))
            I_a0 = -I_sequence * (Z2_th / (Z2_th + Z0_th))
            I_fault_total = I_sequence  
            
        elif n == 2:
            I_sequence = E_th_mk / (Z1_th + Z2_th + Z0_th)
            I_a1 = I_sequence
            I_a2 = I_sequence
            I_a0 = I_sequence
            I_fault_total = I_sequence  
            
        # Distribute sequence voltages globally across all buses due to the series disturbance
        V1_new = V1_unf - (Zbus_1[:, m] - Zbus_1[:, k]) * I_a1
        V2_new = - (Zbus_1[:, m] - Zbus_1[:, k]) * I_a2
        V0_new = - (Zbus_0[:, m] - Zbus_0[:, k]) * I_a0
        
        v_diff = np.max(np.abs(V1_new - V1_profile))
        V1_profile = V1_new
        
        if v_diff < tol:
            converged = True
            if verbose > 1:
                print(f"   🔄 Iterative {n}OP Solver converged in {iteration+1} iterations.")
            break

    if not converged and verbose > 0:
        print(f"   ⚠️ Iterative {n}OP Solver reached max iterations without strict convergence.")

    # 4. Compile Results
    Va_profile = V1_profile + V2_new + V0_new

    return {
        'I_fault_pu': abs(I_fault_total),
        'Z_th_mag': abs(Z1_th + Z2_th + Z0_th) / 3 if n == 2 else abs(Z1_th),
        'V_profile': [abs(v) for v in Va_profile]
    }

def analyze_fault(grid: Grid, path: str | None = None, fault_type: str = '3PH', Z_fault: complex = 0j, verbose: int = 0):
    """
    Unified Batch Analysis Engine. Automatically scans the entire grid and applies the specified 
    fault type sequentially to all valid system components (Buses or Transmission Lines).
    
    Supported Fault Types:
    - Shunt Faults (applied at Bus): '3PH', 'SLG', 'LL', 'DLG'
    - Series Faults (applied at Line): '1OP', '2OP'
    """
    results = []
    bus_names = [bus.name for bus in grid.buses]
    fault_type = fault_type.upper()
    
    print(f"🔍 Scanning grid for {fault_type} Fault...")

    # ==========================================
    # 1. Routing & Solvers Mapping
    # ==========================================
    shunt_solvers = {
        '3PH': ThreePhaseFault,
        'SLG': LineToGroundFault,
        'LL': LineToLineFault,
        'DLG': DoubleLineToGroundFault
    }
    
    series_faults = ['1OP', '2OP'] 

    if fault_type not in shunt_solvers and fault_type not in series_faults:
        raise ValueError(f"Unknown fault type: {fault_type}. Available: {list(shunt_solvers.keys()) + series_faults}")

    # ==========================================
    # 2. Execution Loop
    # ==========================================
    
    # --- A. Shunt Faults (Occur at Buses) ---
    if fault_type in shunt_solvers:
        solver_func = shunt_solvers[fault_type]
        for bus in grid.buses:
            data = solver_func(grid, fault_bus=bus.name, Z_fault=Z_fault, verbose=verbose)
            
            v_base = getattr(bus, 'Vbase', 11.0) 
            s_base = getattr(grid, 'Sbase', 100.0) 
            
            i_base = (s_base * 1000) / (np.sqrt(3) * v_base)
            i_fault_ka = (data['I_fault_pu'] * i_base) / 1000
            s_sc_mva = s_base * data['I_fault_pu']
            
            row = {
                "Fault_Location": f"Bus_{bus.name}",
                "Fault_Type": fault_type,
                "I_fault_pu": round(data['I_fault_pu'], 4),
                "I_fault_kA": round(i_fault_ka, 4),
                "S_sc_MVA": round(s_sc_mva, 2),
                "Eqv_Z_th_mag": round(data['Z_th_mag'], 4)
            }
            
            for i, b_name in enumerate(bus_names):
                row[f"V_{b_name}_pu"] = round(data['V_profile'][i], 4)
            results.append(row)
            
    # --- B. Series Faults (Occur at Branches/Lines) ---
    elif fault_type in series_faults:
        n_open = int(fault_type[0])  # Extract the number 1 or 2
        
        # Loop through all transmission lines (Branches)
        for u, v, data_edge in grid.edges(data=True):
            branch = data_edge.get('obj')
            if branch is None or type(branch).__name__ != 'TransmissionLine':
                continue  # Skip transformers; open conductor faults are typically evaluated strictly on lines
                
            data = OpenConductorFault(grid, from_bus=u, to_bus=v, n=n_open, verbose=verbose)
            
            # Utilize the Base Voltage of the connecting bus to derive physical currents
            u_bus = grid.bus(u)
            v_base = getattr(u_bus, 'Vbase', 11.0) 
            s_base = getattr(grid, 'Sbase', 100.0) 
            
            i_base = (s_base * 1000) / (np.sqrt(3) * v_base)
            i_fault_ka = (data['I_fault_pu'] * i_base) / 1000
            
            row = {
                "Fault_Location": f"Line_{u}-{v}",
                "Fault_Type": fault_type,
                "I_fault_pu": round(data['I_fault_pu'], 4),
                "I_fault_kA": round(i_fault_ka, 4),
                "S_sc_MVA": 0.0,  # Series faults do not directly cause a short circuit MVA drain
                "Eqv_Z_th_mag": round(data['Z_th_mag'], 4)
            }
            
            for i, b_name in enumerate(bus_names):
                row[f"V_{b_name}_pu"] = round(data['V_profile'][i], 4)
            results.append(row)

    # ==========================================
    # 3. Compile Results
    # ==========================================
    df = pd.DataFrame(results)
    
    if path:
        df.to_csv(path, index=False)
        print(f"✅ Engineering Fault Analysis complete. Report saved to {path}")
    else: 
        print(f"✅ Engineering Fault Analysis complete.")
    
    return results


def analyze_faults(grid: Grid, path: str | None = None, fault_types: list[str] = ['3PH', 'SLG', 'LL', 'DLG', '1OP', '2OP'], Z_fault: complex = 0j, verbose: int = 0):
    """
    Comprehensive Batch Analysis Engine.
    Executes multiple fault scenarios simultaneously and generates an aggregated Master Report 
    containing the peak short-circuit currents (kA) mapped out across the entire power grid.
    """
    print("\n" + "="*50)
    print("🌍 INITIATING ULTIMATE FAULT ANALYSIS")
    print("="*50)
    
    master_results = []
    
    for f_type in fault_types:
        type_results = analyze_fault(grid, path=None, fault_type=f_type, Z_fault=Z_fault, verbose=verbose)
        master_results.extend(type_results)
        
    df_master = pd.DataFrame(master_results)
    
    if path:
        df_master.to_csv(path, index=False)
        print(f"\n✅ ULTIMATE REPORT SAVED TO: {path}")
        
        try:
            print("\n📊 FAULT CURRENT SUMMARY (kA):")
            
            # Pivot the data to create a clean comparison matrix
            summary_df = df_master.pivot(index="Fault_Location", columns="Fault_Type", values="I_fault_kA")
            summary_df['MAX_kA'] = summary_df.max(axis=1)
            
            # Enforce strict numerical sorting by extracting the integer ID from the location string
            # This prevents pandas from mis-sorting bus arrays alphabetically (e.g., Bus_1, Bus_10, Bus_2)
            sort_keys = summary_df.index.to_series().str.extract(r'(\d+)').astype(float)
            summary_df['sort_key'] = sort_keys[0]
            summary_df = summary_df.sort_values(by='sort_key').drop(columns=['sort_key'])
            
            print(summary_df.to_string())
        except Exception as e:
            print(f"Failed to generate summary table: {e}")

    return df_master