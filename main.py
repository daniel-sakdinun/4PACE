from psys import Grid

grid = Grid.load('config.yml')

grid.solve()
grid.eco_dispatch()