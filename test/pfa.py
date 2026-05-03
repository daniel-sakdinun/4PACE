import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from fourpace.psys import Grid
from fourpace.pfa import CEP, SCOPF, Validate_N1
import pandas as pd

grid = Grid.load('config.yaml')

load_profile = pd.read_csv('profile.csv')
grid.attach_profile(load_profile)

CEP(grid, solver='CLARABEL')

rescue_plan = SCOPF(grid, solver='CLARABEL')

if rescue_plan:
    Validate_N1(grid, rescue_plan)