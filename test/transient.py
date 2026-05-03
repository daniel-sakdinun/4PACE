import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from fourpace.psys import Grid
from fourpace.dynamics import find_cct

grid = Grid.load('config.yaml')

find_cct(grid, fault_bus="2", t_min=0.01, t_max=0.50, tol=0.002, path='rotor_angle.csv')