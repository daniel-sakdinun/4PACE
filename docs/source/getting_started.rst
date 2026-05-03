===============
Getting Started
===============

Welcome to 4PACE! This guide will help you set up the environment and run your first modern power system simulation.

Installation
------------
Install 4PACE directly from PyPI using pip:

.. code-block:: bash

    pip install 4pace

Grid Configuration (YAML)
------------

4PACE utilizes a highly intuitive, declarative YAML format to define your power system topology. You can easily configure complex parameters like Sequence Impedances (:math:`X_d''`, :math:`X_0`), IBR fault current limits, and Transformer grounding connections without writing a single line of Python code.

A standard 4PACE configuration file requires three main blocks: ``Sbase``, ``buses``, and ``branches``.

.. code-block:: yaml

    # 🌍 System Base
    Sbase: 100.0  # System Base Power (MVA)

    # 🏢 Buses & Components
    buses:
      - name: "1"
        Vbase: 11.0  # Base Voltage (kV)
        components:
          - type: "SynchronousMachine"
            name: "G1"
            P_max: 300.0
            # ⚡ Sequence Impedances for Fault Analysis
            Xd_sub: 0.15   # Positive/Negative Seq Subtransient Reactance (pu)
            X0: 0.05       # Zero Seq Reactance (pu)

      - name: "2"
        Vbase: 22.0
        components:
          - type: "Inverter"
            name: "Solar_Farm_A"
            S_max: 50.0
            # ☀️ IBR Fault Logic
            I_fault_limit_pu: 1.2  # Limits reactive fault current to 1.2x of rated

          - type: "Battery"
            name: "BESS_1"
            P_max: 20.0
            capacity_mwh: 40.0
            I_fault_limit_pu: 1.5  # BESS can inject up to 1.5x during faults

    # 🌉 Branches (Lines & Transformers)
    branches:
      - type: "TransmissionLine"
        from_from: "1"
        to_bus: "2"
        R: 0.019
        X: 0.059
        # 📏 Zero Sequence Line Data (Crucial for SLG, DLG, Open Conductor)
        X0: 0.177 

      - type: "Transformer"
        from_bus: "2"
        to_bus: "3"
        X: 0.08
        # 🔄 Transformer Grounding Configuration (Controls Zero Sequence Path)
        # Available options: "Yg-Yg", "Delta-Yg", "Yg-Delta", "Delta-Delta"
        connection_type: "Delta-Yg" 
        tap_ratio: 1.0

* **Xd_sub and X0**: Mandatory for ``SynchronousMachine`` to establish Norton Equivalent Sources correctly.
* **I_fault_limit_pu**: Crucial for ``Inverter`` and ``Battery``. If omitted, the system defaults to standard limit (usually 1.2).
* **connection_type**: Essential for ``Transformer``. It acts as a switch determining if Zero Sequence current can pass through the grounding path.

Quick Start: Microgrid Investment Planning & N-1 Security
---------------------------------------------------------
Co-optimize Solar/BESS CapEx against 24-hour OpEx, then run a Multiverse Security-Constrained OPF to ensure the grid survives any single transmission line failure.

.. code-block:: python

    import pandas as pd
    from fourpace.psys import Grid
    from fourpace.pfa import CEP, SCOPF, Validate_N1

    grid = Grid.load('config.yaml')

    # Attach a 24-hour load and solar profile
    load_profile = pd.read_csv('profile.csv')
    grid.attach_profile(load_profile)

    # 1. Capacity Expansion Planning (Find optimal BESS/Solar sizes)
    CEP(grid, solver='CLARABEL')

    # 2. Security-Constrained OPF (Generate rescue plans for N-1 contingencies)
    rescue_plan = SCOPF(grid, solver='CLARABEL')

    # 3. Non-Linear Physics Auditor (Validate AI-generated dispatch with AC Load Flow)
    if rescue_plan:
        Validate_N1(grid, rescue_plan)

Quick Start: The "One-Click" Ultimate Fault Analysis
----------------------------------------------------
Define your grid, sequence impedances (:math:`X_d''`, :math:`X_0`), and IBR current limits directly in your `config.yaml`, then run the Iterative Fault Analyzer.

.. code-block:: python

    from fourpace.psys import Grid
    from fourpace.fault import analyze_faults

    # 1. Load system configurations (includes Sequence Network & IBR data)
    grid = Grid.load('config.yaml')

    # 2. Run Comprehensive Iterative Fault Analysis (Supports Solar & BESS limits!)
    # Analyzes: 3PH, SLG, LL, DLG, 1OP, 2OP
    master_report_df = analyze_faults(grid, path="Master_Fault_Report.csv", verbose=0)

Quick Start: Transient Stability & Critical Clearing Time (CCT)
---------------------------------------------------------------
Simulate multi-machine dynamics equipped with AVRs, Governors, and FACTS devices to find the exact millisecond a network loses synchronism.

.. code-block:: python

    from fourpace.psys import Grid
    from fourpace.dynamics import find_cct

    grid = Grid.load('config.yaml')

    # Launch the Binary Search Engine to find the ultimate CCT for a fault at Bus 2
    find_cct(grid, fault_bus="2", t_min=0.01, t_max=0.50, tol=0.002, path='rotor_angle.csv')