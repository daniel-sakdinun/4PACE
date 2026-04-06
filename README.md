# python_power_dispatch# ⚡ Python Power System Analysis & Dispatch Engine

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A lightweight, high-performance Power System Analysis and Optimization engine built from scratch using pure Python. This project leverages Graph Theory and advanced Numerical Optimization to solve complex power grid problems.

## ✨ Key Features

* **Network Topology Management:** Built on top of `networkx` for robust and intuitive graph-based grid modeling.
* **Steady-State Analysis:** Accurate **Newton-Raphson Power Flow** algorithm for evaluating grid voltage and power distribution.
* **Optimal Power Flow (OPF):** * **DC OPF:** Fast linear optimization for active power dispatch.
  * **AC OPF:** Full non-linear optimization considering $I^2R$ losses, voltage limits, reactive power (Q), and branch flow limits ($S_{max}$).
* **Economic Dispatch:** Lambda iteration method for minimizing generation fuel costs.
* **Advanced Equipment Models:**
  * **Synchronous Machines:** Supports both Generator and Motor modes with PQ limits.
  * **Branch Models:** Pi-model Transmission Lines and Off-nominal Tap / Phase-shifting Transformers.
  * **Loads:** Supports Constant Power (PQ), Constant Current (I), and Constant Impedance (Z) models.
* **Human-Readable Config:** Effortlessly load entire grid topologies via `YAML` or `JSON` without hardcoding Python scripts.

## 📦 Installation

Clone this repository and install the required dependencies. This engine is intentionally kept lightweight and relies only on core scientific libraries.

```bash
git clone [https://github.com/yourusername/python_power_dispatch.git](https://github.com/yourusername/python_power_dispatch.git)
cd python_power_dispatch
pip install numpy scipy networkx pyyaml
