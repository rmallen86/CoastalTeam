"""Mapping functions."""
from landlab.grid.mappers import map_link_vector_components_to_node as mln
import numpy as np


def map_velocity_components_to_nodes(grid):
    """Map the velocity components from the links to the nodes.

    Inputs :
        grid - A landlab grid object

    Returns :
        ebb_vel_x - ebb x component of flow velocity
        ebb_vel_y - ebb y component of flow velocity
        flood_vel_x - flood x component of flow velocity
        flood_vel_y - flood y component of flow velocity

    """
    ebb_vel_x, ebb_vel_y = mln(grid, grid.at_link['ebb_tide_flow__velocity'])
    flood_vel_x = -ebb_vel_x
    flood_vel_y = -ebb_vel_y
    return (ebb_vel_x, ebb_vel_y, flood_vel_x, flood_vel_y)


class gridded_vars:
    """Class to hold the gridded variable data."""

    def __init__(self, grid):
        """Initialize the class.

        Inputs :
            grid - A landlab grid object.

        """
        # get velocity components
        (eb_x, eb_y, fl_x, fl_y) = map_velocity_components_to_nodes(grid)

        # convert components to attributes of the class
        self.ex = np.flipud(np.reshape(eb_x, grid.shape))  # ebb x
        self.ey = np.flipud(np.reshape(eb_y, grid.shape))  # ebb y
        self.fx = np.flipud(np.reshape(fl_x, grid.shape))  # flood x
        self.fy = np.flipud(np.reshape(fl_y, grid.shape))  # flood y
        self.elev = np.reshape(grid.at_node['topographic__elevation'],
                               grid.shape)
        self.depth = np.reshape(grid.at_node['mean_water__depth'],
                                grid.shape)
