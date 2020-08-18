"""Couple landlab tidal-flow with dorado passive particles."""
# this example uses Greg Tuckers' 'Straight Channel example'

import matplotlib.pyplot as plt
from landlab.components import TidalFlowCalculator
from landlab import RasterModelGrid
from map_fun import gridded_vars
from plot_fun import group_plot
from particletransport import init_particles
from particletransport import tidal_particles
import numpy as np

from landlab.grid.mappers import map_max_of_link_nodes_to_link

# Straight Channel example
# parameters
nrows = 400
ncols = 200
grid_spacing = 2.0  # m
marsh_height = 1.0  # m
channel_depth = 2.0  # m
tidal_range = 3.1  # m
tidal_period = 12.5 * 3600.0  # s
n_tide_periods = 100  # number of tidal periods
open_nodes = np.arange(94, 105, dtype=np.int)  # IDs of open-boundary nodes (along channel at bottom/south boundary)
roughness_shallow = 0.2  # Manning's n for areas above mean sea level (i.e., the marsh)
roughness_deep = 0.01  # Manning's n for areas below mean sea level (i.e., the channel)

# create and set up the grid
grid = RasterModelGrid((nrows, ncols), xy_spacing=grid_spacing)
z = grid.add_zeros('topographic__elevation', at='node')
z[grid.core_nodes] = marsh_height
channel = np.logical_and(grid.x_of_node >= 188.0, grid.x_of_node <= 208.0)
z[channel] = -channel_depth
grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
grid.status_at_node[open_nodes] = grid.BC_NODE_IS_FIXED_VALUE

# set up roughness field (calculate on nodes, then map to links)
roughness_at_nodes = roughness_shallow + np.zeros(z.size)
roughness_at_nodes[z < 0.0] = roughness_deep
roughness = grid.add_zeros('roughness', at='link')
map_max_of_link_nodes_to_link(grid, roughness_at_nodes, out=roughness)

# instantiate the TidalFlowCalculator
tfc = TidalFlowCalculator(grid, tidal_range=tidal_range,
                          tidal_period=tidal_period, roughness='roughness')

# run it
tfc.run_one_step()

# get gridded values
gvals = gridded_vars(grid)

# initialize the particle parameters
seed_xloc = list(range(96, 103))
seed_yloc = list(range(195, 200))
Np_tracer = 200
params = init_particles(seed_xloc, seed_yloc, Np_tracer, grid_spacing, gvals)

# move the particles with the tides
walk_data = tidal_particles(params, tidal_period/100, n_tide_periods)

# make figure
group_plot(gvals)
plt.savefig('group_plot.png')
plt.close()
