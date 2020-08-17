"""Mapping velocity components."""
from landlab.grid.mappers import map_link_vector_components_to_node as mln


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
