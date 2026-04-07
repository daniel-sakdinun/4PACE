# ⚡ Four-Quadrant Power Analysis & Computational Engine (4PACE)

[![Status](https://img.shields.io/badge/Status-Early%20Access-orange.svg)]()
[![PyPI version](https://img.shields.io/pypi/v/4pace.svg)](https://pypi.org/project/4pace/)
[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **⚠️ Early Access Notice**
> 
> This project is currently in **Early Access (v0.2.0a0)**. While the core engine is functional and capable of solving complex power flow and optimization problems, it is still undergoing significant development. Features may change, and bug reports or feedback are highly encouraged to help reach the stable v1.0.0 release.

**4PACE** is a high-precision Power System Optimization Engine designed for analyzing and determining the most efficient operating points for electrical grids. Developed entirely in Python using only core scientific libraries, it ensures maximum computational efficiency and provides the flexibility needed to scale capabilities for modern grid demands.

---

## 🧐 Why "4PACE"?

The name 4PACE reflects the core pillars of this project:
* **4-Quadrant:** The ability to simulate electrical equipment behavior across all four quadrants of active and reactive power, whether operating as a Generator, Motor, or Condenser.
* **Power Analysis:** Comprehensive steady-state analysis powered by robust algorithms like Fast Decoupled Power Flow (FDPF).
* **Computational Engine:** Driven by the mathematical prowess of NumPy and SciPy to ensure rapid numerical convergence.

---

* **Network Topology Management:** Built on top of `networkx` for robust and intuitive graph-based grid modeling.
* **Advanced Steady-State Analysis:** Highly optimized **Fast Decoupled Power Flow (FDPF)** algorithm for rapid convergence, complete with **Dynamic Q-Limit Enforcement** (automatic PV to PQ bus switching).
* **Optimal Power Flow (OPF):** Full non-linear AC OPF considering $I^2R$ losses, voltage bounds, reactive power (Q), and hard branch flow limits ($S_{max}$).
* **Grid Security & Monitoring:** Automatic transmission line overload checking and load percentage reporting.
* **Advanced Equipment Models:**
  * **Synchronous Machines:** Supports Generator, Motor, and Condenser modes with strict PQ limits.
  * **Asynchronous Machines (NEW):** Induction motor models featuring real-time slip ($s$) and torque-balance dynamics that actively respond to grid voltage variations.
  * **Branch Models:** Pi-model Transmission Lines and **OLTC Transformers** (On-Load Tap Changers with automatic voltage regulation and tap step adjustments).
  * **Shunts (NEW):** Capacitor banks and Reactors with voltage-dependent reactive power injection ($Q \propto V^2$).
  * **Loads:** Supports Constant Power (PQ), Constant Current (I), and Constant Impedance (Z) models.
* **Human-Readable Config:** Effortlessly load entire grid topologies via `YAML` without hardcoding Python scripts.

---

## 📦 Installation

Since 4PACE is now available on PyPI, you can install the latest early access version directly using `pip`:

```bash
pip install 4pace
```

Note: Requires Python 3.10 or higher.

---

## 🚀 Quick Start (QOL Feature)
**4PACE** allows you to separate your grid configuration from your logic.

1. Define your grid (config.yaml)
Create a configuration file to specify your buses, machines, lines, and advanced features like OLTC:

```yaml
Sbase: 100.0
buses:
  - name: A01
    Vbase: 6.9
    bus_type: Slack
    components:
      - type: SynchronousMachine
        name: Gen1
        S_rated: 100
        a: 500
        b: 7
        c: 0.004
  - name: B01
    Vbase: 22
    bus_type: PQ
    components:
      - type: Load
        name: City_Mall
        P: 240
        Q: 80
      - type: Shunt
        name: CapBank1
        Q_nom: 50
branches:
  - type: Transformer
    name: TR1
    from_bus: A01
    to_bus: B01
    R: 0.002
    X: 0.04
    auto_tap: true
    controlled_bus: B01
    target_V: 1.0
```


2. Run the Solver

```python
from fourpace.psys import Grid # Professional Clean Import

# Load your configuration
grid = Grid.load('config.yaml') 

# Run Optimal Power Flow
grid.eco_dispatch()

# Solve Fast Decoupled Power Flow (with OLTC & Q-Limits)
grid.solve()

# Check for Transmission Line Overloads
grid.check_overload()
```

---

## 🗺️ Roadmap
[ ] Integration of Inverter-Based Resources (IBRs).

[ ] Transient and Dynamic Stability Analysis.

[ ] Interactive Web-based Visualization Dashboard.

## 🤝 Contributing
Contributions are welcome! If you find a bug or have a feature request, please open an issue on our GitHub repository.

---

Developed with ❤️ for the Power Engineering Community.
