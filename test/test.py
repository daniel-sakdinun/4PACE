import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from fourpace.psys import Grid

grid = Grid.load('config.yaml')

grid.eco_dispatch()
grid.solve()
grid.check_overload()