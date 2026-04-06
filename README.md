# ⚡ Four-Quadrant Power Analysis & Computational Engine (4PACE)

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

4PACE is a high-precision Power System Optimization Engine designed for analyzing and determining the most efficient operating points for electrical grids. Developed entirely in Python using only core scientific libraries, it ensures maximum computational efficiency and provides the flexibility needed to scale capabilities for modern grid demands.

## 🧐 Why "4PACE"?

The name 4PACE reflects the core pillars of this project:
* **4-Quadrant:** The ability to simulate electrical equipment behavior across all four quadrants of active and reactive power, whether operating as a Generator, Motor, or Condenser.
* **Power Analysis:** Comprehensive steady-state analysis powered by a robust Newton-Raphson solver.
* **Computational** Engine: Driven by the mathematical prowess of NumPy and SciPy to ensure rapid numerical convergence.

## ✨ Key Features

* **Network Topology Management:** Built on top of `networkx` for robust and intuitive graph-based grid modeling.
* **Steady-State Analysis:** Accurate **Newton-Raphson Power Flow** algorithm for evaluating grid voltage and power distribution.
* **Optimal Power Flow (OPF):**
  * **AC OPF:** Full non-linear optimization considering $I^2R$ losses, voltage limits, reactive power (Q), and branch flow limits ($S_{max}$).
* **Economic Dispatch:** Lambda iteration method for minimizing generation fuel costs.
* **Advanced Equipment Models:**
  * **Synchronous Machines:** Supports both Generator and Motor modes with PQ limits.
  * **Branch Models:** Pi-model Transmission Lines and Off-nominal Tap / Phase-shifting Transformers.
  * **Loads:** Supports Constant Power (PQ), Constant Current (I), and Constant Impedance (Z) models.
* **Human-Readable Config:** Effortlessly load entire grid topologies via `YAML` without hardcoding Python scripts.

## 📦 Installation

Clone this repository and install the required dependencies. This engine is intentionally kept lightweight and relies only on core scientific libraries.

```bash
git clone [https://github.com/yourusername/python_power_dispatch.git](https://github.com/yourusername/python_power_dispatch.git)
cd python_power_dispatch
pip install -r requirements.txt
