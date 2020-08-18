"""Couple landlab tidal-flow with dorado passive particles."""
# This example uses a modified version of Greg Tucker's idealized 2D example.

import matplotlib.pyplot as plt
from landlab.components import TidalFlowCalculator
from landlab import RasterModelGrid
from map_fun import gridded_vars
from plot_fun import group_plot
from particletransport import init_particles
from particletransport import tidal_particles

# 2D raster example
# parameters
nrows = 150
ncols = 250
grid_spacing = 100.0  # m
mean_depth = 2.0  # m
tidal_range = 2.0  # m
roughness = 0.01  # s/m^1/3, i.e., Manning's n
tide_period = 2*60*60  # tidal period in seconds
n_tide_periods = 15  # number of tidal periods to move particles around for

# create and set up the grid
grid = RasterModelGrid((nrows, ncols), xy_spacing=grid_spacing)
z = grid.add_zeros('topographic__elevation', at='node')
z[:] = -mean_depth
grid.set_closed_boundaries_at_grid_edges(False, False, True, True)

# instantiate the TidalFlowCalculator
tfc = TidalFlowCalculator(grid, tidal_range=tidal_range, roughness=0.01)

# run it
tfc.run_one_step()

# get gridded values
gvals = gridded_vars(grid)

# initialize the particle parameters
seed_xloc = list(range(70, 180))
seed_yloc = list(range(50, 110))
Np_tracer = 100
params = init_particles(seed_xloc, seed_yloc, Np_tracer, grid_spacing, gvals)

# move the particles with the tides
walk_data = tidal_particles(params, tide_period, n_tide_periods)

# make figure
group_plot(gvals)
plt.savefig('group_plot.png')
plt.close()
