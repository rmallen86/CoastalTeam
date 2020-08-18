"""Landlab tidal-flow-calculator with random vegetated field."""

import matplotlib.pyplot as plt
from landlab.components import TidalFlowCalculator
from landlab import RasterModelGrid
from landlab.grid.mappers import map_max_of_link_nodes_to_link
from map_fun import gridded_vars
from plot_fun import group_plot
from plot_fun import plot_depth
from particletransport import init_particles
from particletransport import tidal_particles
import gstools as gs
import numpy as np

# 2D domain using a square
# parameters
nrows = 250
ncols = nrows
grid_spacing = 1.0  # m
mean_depth = 1.5  # m
tidal_range = 0.5  # m
roughness_low = 0.01  # s/m^1/3, i.e., Manning's n
roughness_high = 0.1  # s/m^1/3, i.e., Manning's n
tide_period = 2*60*60  # tidal period in seconds
n_tide_periods = 50  # number of tidal periods to move particles around for

# generate a random field with gstools
# (https://github.com/GeoStat-Framework/GSTools)
seed = 1  # defines the random seed
len_scale = 10  # length scale - integer = isotropic; list [1, 5] = anisotropic

x = y = range(nrows)
model = gs.Gaussian(dim=2, var=1, len_scale=len_scale)
srf = gs.SRF(model, seed=seed)
srf.structured([x, y])
gs.transform.binary(srf)
# get array info
srf_array = srf.field

# create and set up the grid
grid = RasterModelGrid((nrows, ncols), xy_spacing=grid_spacing)
z = grid.add_zeros('topographic__elevation', at='node')
grid.set_closed_boundaries_at_grid_edges(True, False, True, False)

# set up roughness field (calculate on nodes, then map to links)
roughness_at_nodes = np.zeros_like(z)
roughness_at_nodes[srf_array.flatten() > 0] = roughness_high  # high roughness
roughness_at_nodes[srf_array.flatten() < 0] = roughness_low  # low roughness
roughness = grid.add_zeros('roughness', at='link')
map_max_of_link_nodes_to_link(grid, roughness_at_nodes, out=roughness)

# instantiate the TidalFlowCalculator
tfc = TidalFlowCalculator(grid, tidal_range=tidal_range,
                          tidal_period=tide_period, roughness='roughness')

# run it
tfc.run_one_step()

# get gridded values
gvals = gridded_vars(grid)

# initialize the particle parameters
seed_xloc = list(range(100, 125))
seed_yloc = list(range(100, 125))
Np_tracer = 100
params = init_particles(seed_xloc, seed_yloc, Np_tracer, grid_spacing, gvals)

# move the particles with the tides
walk_data = tidal_particles(params, tide_period/10, n_tide_periods,
                            plot_grid=np.flipud(np.reshape(roughness_at_nodes,
                                                           grid.shape)))

# make figure
group_plot(gvals)
plt.savefig('group_plot.png')
plt.close()

plot_depth(grid)
plt.savefig('depth_plot.png')
plt.close()
