import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from fourpace.psys import Grid
from fourpace.fault import analyze_faults

grid = Grid.load('config.yaml')

analyze_faults(grid, path="IEEE14_Batch_Fault.csv", verbose=-1)