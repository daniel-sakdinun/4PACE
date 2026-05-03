=========================
Grid Configuration (YAML)
=========================

The ``config.yaml`` file is the core of 4PACE. It defines the network topology, component specifications, economic planning data, and complex dynamic controllers.

This guide serves as an exhaustive reference for all supported structures and parameters.

Global Settings
===============
These parameters define the base values for the entire per-unit (p.u.) system.

.. list-table:: 
   :widths: 20 20 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``Sbase``
     - Float
     - **[Required]** System base apparent power in MVA (e.g., 100.0).

Buses & Components
==================
The ``buses`` list contains all nodes in the network. Each bus can host multiple components.

Bus Attributes
--------------
* ``name`` (String): **[Required]** Unique identifier for the bus.
* ``Vbase`` (Float): **[Required]** Base voltage of the bus in kV.
* ``bus_type`` (String): Reference for Load Flow (e.g., 'Slack', 'PV', 'PQ').
* ``components`` (List): A list of electrical devices attached to this bus.

Component Types
---------------
Inside the ``components`` list, you must specify the ``type`` of the device. Below are the supported components and their parameters:

**SynchronousMachine** (Generators)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Represents generators, motors, or condensers with full dynamic and fault capabilities.

* *Basic & OPF*: ``name``, ``P``, ``Q``, ``Pmin``, ``Pmax``, ``Qmin``, ``Qmax``, ``S_rated``, ``pf``, ``mode`` ('generator', 'motor', 'condenser'), ``a``, ``b``, ``c`` (cost curve parameters).
* *Fault & Impedance*: ``R``, ``X``, ``Xd``, ``Xd_prime``, ``Xd_sub`` (Crucial for Faults), ``X2``, ``X0``.
* *Dynamics*: ``H`` (Inertia), ``Td0_prime``.
* *Controllers*: Nested blocks for ``avr``, ``gov``, and ``pss`` (See Controllers section below).
* *Planning (CEP)*: ``is_candidate``, ``capex_per_mw``, ``max_build_mw``, ``lifetime_years``, ``interest_rate``.

**AsynchronousMachine** (Induction Motor)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* *Basic*: ``name``, ``P_rated``, ``V_rated``, ``poles``, ``freq``.
* *Circuit*: ``Rs``, ``Xs``, ``Rr``, ``Xr``, ``Xm``.
* *State*: ``s`` (slip), ``load_type`` (e.g., 'constant_torque').

**Inverter** (Solar PV / Wind)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* *Basic*: ``name``, ``S_max``, ``P``, ``Q``, ``control_mode``, ``source_type`` ('solar', 'wind', 'bess').
* *Fault Logic*: ``I_fault_limit_pu`` (Max reactive current injection limit during faults, e.g., 1.2).
* *Planning (CEP)*: ``is_candidate``, ``capex_per_mw``, ``max_build_mw``, ``lifetime_years``.

**Battery** (BESS)
^^^^^^^^^^^^^^^^^^
* *Basic*: ``name``, ``P_max``, ``E_max``, ``init_soc``, ``eta`` (efficiency).
* *Fault Logic*: ``I_fault_limit_pu`` (e.g., 1.5).
* *Planning (CEP)*: ``is_candidate``, ``capex_per_mw``, ``capex_per_mwh``, ``max_build_mw``, ``max_build_mwh``.

**Load**
^^^^^^^^
* *Basic*: ``name``, ``model`` ('P', 'I', 'Z'), ``P`` (MW), ``Q`` (MVAr), ``R``, ``X``.

**Shunt** (Capacitors / Reactors)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* *Basic*: ``name``, ``Q_nom`` (MVAr), ``V_nom`` (kV).

Controllers (Nested in SynchronousMachine)
==========================================

**SEXS** (AVR)
--------------
Simplified Excitation System (IEEE Type 1 Equivalent).
* ``type``: "SEXS"
* ``Ka`` (Gain), ``Ta`` (Time constant), ``Efd_min``, ``Efd_max`` (Limits).

**TGOV1** (Governor)
--------------------
IEEE Standard Steam Turbine-Governor Model.
* ``type``: "TGOV1"
* ``R`` (Speed droop), ``T1``, ``T2``, ``T3`` (Time constants), ``Vmax``, ``Vmin`` (Valve limits), ``Dt`` (Damping).

**PSS1A** (Power System Stabilizer)
-----------------------------------
IEEE Standard Power System Stabilizer.
* ``type``: "PSS1A"
* ``K_pss`` (Gain), ``T_w`` (Washout), ``T1``, ``T2``, ``T3``, ``T4`` (Lead-lag time constants), ``V_max``, ``V_min``.

FACTS Devices
=============

**CSVGN1** (SVC - Shunt FACTS)
------------------------------
Standard Static Var Compensator. Resides inside the bus ``components`` list.
* ``type``: "CSVGN1"
* ``name``, ``V_ref`` (Target Voltage), ``K_svc`` (Gain), ``T_1`` (Delay), ``B_max``, ``B_min`` (Susceptance limits).

**STATCOM1** (Shunt FACTS)
--------------------------
Voltage Source Converter based Shunt. Resides inside the bus ``components`` list.
* ``type``: "STATCOM1"
* ``name``, ``V_ref``, ``K_r`` (Gain), ``T_r`` (Delay), ``Iq_max``, ``Iq_min`` (Reactive current limits).

**TCSC1** (Series FACTS)
------------------------
Thyristor Controlled Series Capacitor. Resides in the global ``series_facts`` list (outside buses).
* ``type``: "TCSC1"
* ``name``, ``branch_name`` (Target line to clamp onto), ``P_ref`` (Target active power flow in p.u.).
* ``K_p`` (Gain), ``T_p`` (Delay), ``X_max``, ``X_min`` (Reactance injection limits).

Branches
========
The ``branches`` list defines the connections between buses.

**TransmissionLine**
--------------------
* *Basic*: ``from_bus``, ``to_bus``, ``R``, ``X``, ``B_shunt``, ``S_max``, ``length_km``.
* *Fault Logic*: ``R0``, ``X0``, ``B0_shunt`` (Zero Sequence Parameters).
* *Planning (CEP)*: ``is_candidate``, ``capex_per_mva``, ``max_build_mva``, ``lifetime_years``.

**Transformer**
---------------
* *Basic*: ``from_bus``, ``to_bus``, ``R``, ``X``, ``S_max``.
* *Tap Changer*: ``tap_ratio``, ``phase_shift``, ``auto_tap``, ``target_V``, ``tap_step``, ``tap_min``, ``tap_max``.
* *Fault Logic*: ``R0``, ``X0``, ``connection_type`` ('yg-yg', 'delta-yg', 'yg-delta', 'delta-delta').
* *Planning (CEP)*: ``is_candidate``, ``capex_per_mva``, ``max_build_mva``.

Example: Modified IEEE 14 Bus System
====================================
Here is how multiple components, complex controllers, and FACTS devices come together:

.. code-block:: yaml

    Sbase: 100.0
    buses:
      - name: "1"
        Vbase: 132.0
        bus_type: Slack
        components:
          - type: SynchronousMachine
            name: G1
            a: 0.0
            b: 20.0
            c: 0.043
            Pmax: 332.4
            Pmin: 0.0
            Qmax: 150.0
            Qmin: -20.0
            Xd_sub: 0.15
            X2: 0.15
            X0: 0.05
            # --- Control Systems ---
            avr:
              type: SEXS
              Ka: 200.0
              Ta: 0.02
            gov:
              type: TGOV1
              R: 0.05
              T1: 0.5
            pss:
              type: PSS1A
              K_pss: 10.0

      - name: "2"
        Vbase: 132.0
        bus_type: PV
        components:
          - type: SynchronousMachine
            name: G2
            P: 40.0
            a: 0.0
            b: 20.0
            c: 0.25
            Pmax: 140.0
            Pmin: 0.0
            Qmax: 100.0
            Qmin: -50.0
            Xd: 2.0
            Xd_sub: 0.15
            Td0_prime: 0.5
            X2: 0.20
            X0: 0.07
            # --- Control Systems ---
            avr:
              type: SEXS
              Ka: 200.0
              Ta: 0.02
            gov:
              type: TGOV1
              R: 0.05
              T1: 0.5
            pss:
              type: PSS1A
              K_pss: 10.0

          - type: Load
            name: L2
            P: 21.7
            Q: 12.7

      - name: "3"
        Vbase: 132.0
        bus_type: PV
        components:
          - type: SynchronousMachine
            name: G3
            P: 0.0
            a: 0.0
            b: 40.0
            c: 0.01
            Pmax: 100.0
            Pmin: 0.0
            Qmax: 80.0
            Qmin: -20.0
            Xd_sub: 0.25
            X2: 0.25
            X0: 0.09
            # --- Control Systems ---
            avr:
              type: SEXS
              Ka: 200.0
              Ta: 0.02
            gov:
              type: TGOV1
              R: 0.05
              T1: 0.5
            pss:
              type: PSS1A
              K_pss: 10.0

          - type: Load
            name: L3
            P: 94.2
            Q: 19.0

      - name: "4"
        Vbase: 132.0
        bus_type: PQ
        components:
          - type: Load
            name: L4
            P: 47.8
            Q: -3.9
          - type: CSVGN1
            name: SVC_Bus4
            V_ref: 1.02
            K_svc: 50.0
            T_1: 0.05
            B_max: 2.0
            B_min: -1.0

      - name: "5"
        Vbase: 132.0
        bus_type: PQ
        components:
          - type: Load
            name: L5
            P: 7.6
            Q: 1.6

      - name: "6"
        Vbase: 33.0
        bus_type: PV
        components:
          - type: SynchronousMachine
            name: G6
            P: 0.0
            a: 0.0
            b: 20.0
            c: 0.01
            Pmax: 100.0
            Pmin: 0.0
            Qmax: 60.0
            Qmin: -20.0
            Xd_sub: 0.25
            X2: 0.25
            X0: 0.10
            # --- Control Systems ---
            avr:
              type: SEXS
              Ka: 200.0
              Ta: 0.02
            gov:
              type: TGOV1
              R: 0.05
              T1: 0.5
            pss:
              type: PSS1A
              K_pss: 10.0

          - type: Load
            name: L6
            P: 11.2
            Q: 7.5

      - name: "7"
        Vbase: 33.0
        bus_type: PQ
        components: []

      - name: "8"
        Vbase: 33.0
        bus_type: PV
        components:
          - type: SynchronousMachine
            name: G8
            P: 0.0
            a: 0.0
            b: 20.0
            c: 0.01
            Pmax: 100.0
            Pmin: 0.0
            Qmax: 60.0
            Qmin: -20.0
            Xd_sub: 0.25
            X2: 0.25
            X0: 0.10
            # --- Control Systems ---
            avr:
              type: SEXS
              Ka: 200.0
              Ta: 0.02
            gov:
              type: TGOV1
              R: 0.05
              T1: 0.5
            pss:
              type: PSS1A
              K_pss: 10.0

      - name: "9"
        Vbase: 33.0
        bus_type: PQ
        components:
          - type: Load
            name: L9
            P: 29.5
            Q: 16.6
          - type: STATCOM1
            name: STAT_Bus9
            V_ref: 1.0
            K_r: 40.0
            T_r: 0.02
            Iq_max: 1.5
            Iq_min: -1.5

      - name: "10"
        Vbase: 33.0
        bus_type: PQ
        components:
          - type: Load
            name: L10
            P: 9.0
            Q: 5.8

      - name: "11"
        Vbase: 33.0
        bus_type: PQ
        components:
          - type: Load
            name: L11
            P: 3.5
            Q: 1.8

      - name: "12"
        Vbase: 33.0
        bus_type: PQ
        components:
          - type: Load
            name: L12
            P: 6.1
            Q: 1.6

      - name: "13"
        Vbase: 33.0
        bus_type: PQ
        components:
          - type: Load
            name: L13
            P: 13.5
            Q: 5.8

      - name: "14"
        Vbase: 33.0
        bus_type: PQ
        components:
          - type: Load
            name: L14
            P: 14.9
            Q: 5.0
          - type: Shunt
            name: Cap14
            Q_nom: 25.0
          - type: Battery
            name: BESS_14
            is_candidate: true
            capex_per_mw: 10000.0
            capex_per_mwh: 15000.0
            max_build_mw: 50.0
            max_build_mwh: 200.0
            lifetime_years: 10
            interest_rate: 0.05
            eta: 0.95
            init_soc: 0.5
            P_max: 45.82
            E_max: 200.00
            I_fault_limit_pu: 1.5
          - type: Inverter
            name: PV_14
            is_candidate: true
            source_type: solar
            capex_per_mw: 50000.0
            max_build_mw: 80.0
            lifetime_years: 20
            interest_rate: 0.05
            S_max: 80
            I_fault_limit_pu: 1.2

    branches:
      # Transmission Lines
      - {type: TransmissionLine, name: "L1-2", from_bus: "1", to_bus: "2", R: 0.01938, X: 0.05917, S_max: 250, R0: 0.05814, X0: 0.17751}
      - {type: TransmissionLine, name: "L1-5", from_bus: "1", to_bus: "5", R: 0.05403, X: 0.22304, S_max: 250, R0: 0.16209, X0: 0.66912}
      - {type: TransmissionLine, name: "L2-3", from_bus: "2", to_bus: "3", R: 0.04699, X: 0.19797, S_max: 150, R0: 0.14097, X0: 0.59391}
      - {type: TransmissionLine, name: "L2-4", from_bus: "2", to_bus: "4", R: 0.05811, X: 0.17632, S_max: 100, R0: 0.17433, X0: 0.52896}
      - {type: TransmissionLine, name: "L2-5", from_bus: "2", to_bus: "5", R: 0.05695, X: 0.17388, S_max: 100, R0: 0.17085, X0: 0.52164}
      - {type: TransmissionLine, name: "L3-4", from_bus: "3", to_bus: "4", R: 0.06701, X: 0.17103, S_max: 150, R0: 0.20103, X0: 0.51309}
      - {type: TransmissionLine, name: "L4-5", from_bus: "4", to_bus: "5", R: 0.01335, X: 0.04211, S_max: 180, R0: 0.04005, X0: 0.12633}
      - {type: TransmissionLine, name: "L6-11", from_bus: "6", to_bus: "11", R: 0.09498, X: 0.19890, S_max: 120, R0: 0.28494, X0: 0.59670}
      - {type: TransmissionLine, name: "L6-12", from_bus: "6", to_bus: "12", R: 0.12291, X: 0.25581, S_max: 100, R0: 0.36873, X0: 0.76743}
      - {type: TransmissionLine, name: "L6-13", from_bus: "6", to_bus: "13", R: 0.06615, X: 0.13027, S_max: 150, R0: 0.19845, X0: 0.39081}
      - {type: TransmissionLine, name: "L7-8", from_bus: "7", to_bus: "8", R: 0.0, X: 0.17615, S_max: 100, R0: 0.0, X0: 0.52845}
      - {type: TransmissionLine, name: "L7-9", from_bus: "7", to_bus: "9", R: 0.0, X: 0.11001, S_max: 100, R0: 0.0, X0: 0.33003}
      - {type: TransmissionLine, name: "L9-10", from_bus: "9", to_bus: "10", R: 0.03181, X: 0.08450, S_max: 100, R0: 0.09543, X0: 0.25350}
      - {type: TransmissionLine, name: "L9-14", from_bus: "9", to_bus: "14", R: 0.12711, X: 0.27038, S_max: 100, R0: 0.38133, X0: 0.81114}
      - {type: TransmissionLine, name: "L10-11", from_bus: "10", to_bus: "11", R: 0.08205, X: 0.19207, S_max: 100, R0: 0.24615, X0: 0.57621}
      - {type: TransmissionLine, name: "L12-13", from_bus: "12", to_bus: "13", R: 0.22092, X: 0.19988, S_max: 100, R0: 0.66276, X0: 0.59964}
      - {type: TransmissionLine, name: "L13-14", from_bus: "13", to_bus: "14", R: 0.17093, X: 0.34802, S_max: 100, R0: 0.51279, X0: 1.04406}

      # Transformers
      - {type: Transformer, name: "T4-7", from_bus: "4", to_bus: "7", R: 0.0, X: 0.20912, tap_ratio: 0.978, R0: 0.0, X0: 0.20912, connection_type: "yg-yg"}
      - {type: Transformer, name: "T4-9", from_bus: "4", to_bus: "9", R: 0.0, X: 0.55618, tap_ratio: 0.969, R0: 0.0, X0: 0.55618, connection_type: "yg-yg"}
      - {type: Transformer, name: "T5-6", from_bus: "5", to_bus: "6", R: 0.0, X: 0.25202, tap_ratio: 0.932, R0: 0.0, X0: 0.25202, connection_type: "yg-yg"}

    series_facts:
      # Adding a TCSC to clamp onto Transmission Line L4-5
      - type: TCSC1
        name: TCSC_L4_5
        branch_name: L4-5
        P_ref: 0.5   # Target flow in pu (50 MW)
        K_p: 2.0
        T_p: 0.05
        X_max: 0.05  # Slight inductive capability
        X_min: -0.3