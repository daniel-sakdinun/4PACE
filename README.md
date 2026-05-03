# ⚡ Four-Quadrant Power Analysis & Computational Engine (4PACE)

[![Status](https://img.shields.io/badge/Status-Early%20Access-orange.svg)]()
[![Version](https://img.shields.io/badge/version-v0.1.0a1-blue.svg)]()
[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![Documentation Status](https://readthedocs.org/projects/4pace/badge/?version=latest)](https://4pace.readthedocs.io/en/latest/?badge=latest)

> **⚠️ Early Access Notice (v0.1.0a1)**
> 
> Welcome to the foundation of 4PACE. Version **0.1.0a1** introduces our fully OOP Differential-Algebraic Equation (DAE) core. It comes equipped with the **Ultimate Transient Stability Engine**, an **IBR Iterative Fault Solver**, and seamless integration for **Series/Shunt FACTS**. This is an early alpha release built for rigorous community testing and feedback.

**4PACE** is a high-precision Power System Optimization and Analysis Engine designed for modern electrical grids. Built entirely in Python utilizing robust libraries (`cvxpy`, `networkx`, `numpy`), 4PACE handles everything from microgrid investment planning to dynamic transient stability on complex mesh networks like the IEEE 14 Bus System.

---

## 📚 Comprehensive Documentation

Ready to build your first modern grid? Access the complete guide, including in-depth YAML configurations, parameter references, and API details here:

👉 **[Read the Full 4PACE Documentation on ReadTheDocs](https://4pace.readthedocs.io/)** 👈

---

## 🛠️ Quick Start Guide

### 📝 Human-Centric Grid Configuration (`config.yaml`)
4PACE utilizes a highly intuitive, declarative YAML format to define your power system topology. You can easily configure complex parameters like Sequence Impedances, IBR fault current limits, and Transformer grounding connections without writing a single line of Python code.
```yaml
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
    from_bus: "1"
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
```

### Example 1: Microgrid Investment Planning & N-1 Security
Co-optimize Solar/BESS CapEx against 24-hour OpEx, then run a Multiverse Security-Constrained OPF to ensure the grid survives any single transmission line failure.
```python
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
```

### Example 2: The "One-Click" Ultimate Fault Analysis
Define your grid, sequence impedances ($X_d''$, $X_0$), and IBR current limits directly in your `config.yaml`, then run the Iterative Fault Analyzer.
```python
from fourpace.psys import Grid
from fourpace.fault import analyze_faults

# 1. Load system configurations (includes Sequence Network & IBR data)
grid = Grid.load('config.yaml')

# 2. Run Comprehensive Iterative Fault Analysis (Supports Solar & BESS limits!)
# Analyzes: 3PH, SLG, LL, DLG, 1OP, 2OP
master_report_df = analyze_faults(grid, path="Master_Fault_Report.csv", verbose=0)
```

### Example 3: Transient Stability & Critical Clearing Time (CCT)
Simulate multi-machine dynamics equipped with AVRs, Governors, and FACTS devices to find the exact millisecond a network loses synchronism.
```python
from fourpace.psys import Grid
from fourpace.dynamics import find_cct

grid = Grid.load('config.yaml')

# Launch the Binary Search Engine to find the ultimate CCT for a fault at Bus 2
find_cct(grid, fault_bus="2", t_min=0.01, t_max=0.50, tol=0.002, path='rotor_angle.csv')
```

---

### 🔬 Design Philosophy: From a Passion Project to Open-Source
4PACE didn't start as a grand corporate venture. It began as a personal passion project right after I took a course in Power System Analysis. Over time, it evolved line by line, gradually acquiring the advanced computational features you would typically only find in expensive, commercial-grade software. Seeing its potential, I decided to share it with the global electrical engineering community.

The goal of 4PACE is not to directly compete with industry giants, but to offer a highly capable, completely transparent alternative. Proprietary software often acts as a "black box", but with 4PACE, users can inspect every line of the underlying code, verify the mathematical models, and directly report issues or contribute fixes. 

To keep the engine lightweight and accessible, 4PACE strictly adheres to a **Minimal Dependency** architecture. By relying on just a few essential libraries (`cvxpy`, `numpy`, `scipy`, `networkx`, `pandas`, `pyyaml`), We drastically reduce complexity. You won't spend weeks trying to decipher the source code. This simplicity and transparency are designed to foster community growth, ensuring 4PACE evolves into a robust, powerful tool that anyone can pick up and use for exactly $0.

---

### 📜 License
This project is licensed under the GNU Affero General Public License v3.0 (AGPLv3). See the LICENSE file for details. This ensures that the core mathematical engine remains open, transparent, and beneficial to the entire engineering community.

---