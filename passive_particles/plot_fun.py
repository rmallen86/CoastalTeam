"""Plotting functions, adapted from: https://github.com/landlab/landlab/blob/gt/tidal-flow-component/notebooks/tutorials/tidal_flow/tidal_flow_calculator.ipynb"""
import matplotlib.pyplot as plt
from landlab import imshow_grid
from map_fun import map_velocity_components_to_nodes as mvcn
import numpy as np


def plot_depth(grid, resample=1):
    """Plot of water depths.

    Inputs :
        grid - A landlab grid object
        resample - Downsampling value

    Returns :
        Draws a figure that can be rendered with plt.show()
    """
    # depth
    plt.figure()
    imshow_grid(grid, grid.at_node['mean_water__depth'],
                cmap='YlGnBu', color_for_closed='g')
    plt.title('Water depth (m)')
    plt.xlabel('Distance (m)')
    plt.ylabel('Distance (m)')


def plot_ebb_quiver(grid, resample=1):
    """Quiver plot of ebb velocities.

    Inputs :
        grid - A landlab grid object
        resample - Downsampling value

    Returns :
        Draws a figure that can be rendered with plt.show()
    """
    (ebb_x, ebb_y, _, _) = mvcn(grid)

    # down-sample for legible quiver plots if needed
    if resample != 1:
        xr = grid.x_of_node.reshape((grid.number_of_node_rows,
                                     grid.number_of_node_columns))[::resample,
                                                                   ::resample]
        yr = grid.y_of_node.reshape((grid.number_of_node_rows,
                                     grid.number_of_node_columns))[::resample,
                                                                   ::resample]
        ebb_xr = ebb_x.reshape((grid.number_of_node_rows,
                                grid.number_of_node_columns))[::resample,
                                                              ::resample]
        ebb_yr = ebb_y.reshape((grid.number_of_node_rows,
                                grid.number_of_node_columns))[::resample,
                                                              ::resample]
    else:
        xr = grid.x_of_node
        yr = grid.y_of_node
        ebb_xr = ebb_x
        ebb_yr = ebb_y

    # ebb tide
    plt.figure()
    imshow_grid(grid, grid.at_node['topographic__elevation'])
    plt.quiver(xr, yr, ebb_xr, ebb_yr)
    plt.title('Ebb Tide')
    plt.xlabel('Distance (m)')
    plt.ylabel('Distance (m)')


def plot_flood_quiver(grid, resample=1):
    """Quiver plot of flood velocities.

    Inputs :
        grid - A landlab grid object
        resample - Downsampling value

    Returns :
        Draws a figure that can be rendered with plt.show()
    """
    (_, _, flood_x, flood_y) = mvcn(grid)

    # down-sample for legible quiver plots if needed
    if resample != 1:
        xr = grid.x_of_node.reshape((grid.number_of_node_rows,
                                     grid.number_of_node_columns))[::resample,
                                                                   ::resample]
        yr = grid.y_of_node.reshape((grid.number_of_node_rows,
                                     grid.number_of_node_columns))[::resample,
                                                                   ::resample]
        fld_xr = flood_x.reshape((grid.number_of_node_rows,
                                  grid.number_of_node_columns))[::resample,
                                                                ::resample]
        fld_yr = flood_y.reshape((grid.number_of_node_rows,
                                  grid.number_of_node_columns))[::resample,
                                                                ::resample]
    else:
        xr = grid.x_of_node
        yr = grid.y_of_node
        fld_xr = flood_x
        fld_yr = flood_y

    # flood tide
    plt.figure()
    imshow_grid(grid, grid.at_node['topographic__elevation'])
    plt.quiver(xr, yr, fld_xr, fld_yr)
    plt.title('Flood Tide')
    plt.xlabel('Distance (m)')
    plt.ylabel('Distance (m)')


def plot_ebb_magnitudes(grid, resample=1):
    """Plot of ebb velocities.

    Inputs :
        grid - A landlab grid object
        resample - Downsampling value

    Returns :
        Draws a figure that can be rendered with plt.show()
    """
    (ebb_x, ebb_y, _, _) = mvcn(grid)

    ebb_vel_magnitude = np.sqrt(ebb_x * ebb_x + ebb_y * ebb_y)
    plt.figure()
    imshow_grid(grid, ebb_vel_magnitude, cmap='magma', color_for_closed='g')
    plt.title('Ebb Tide Velocity Magnitude (m/s)')
    plt.xlabel('Distance (m)')
    plt.ylabel('Distance (m)')


def plot_flood_magnitudes(grid, resample=1):
    """Plot of flood velocities.

    Inputs :
        grid - A landlab grid object
        resample - Downsampling value

    Returns :
        Draws a figure that can be rendered with plt.show()
    """
    (_, _, flood_x, flood_y) = mvcn(grid)
    plt.figure()
    flood_vel_magnitude = np.sqrt(flood_x * flood_x + flood_y * flood_y)
    imshow_grid(grid, flood_vel_magnitude, cmap='magma', color_for_closed='g')
    plt.title('Flood Tide Velocity Magnitude (m/s)')
    plt.xlabel('Distance (m)')
    plt.ylabel('Distance (m)')


def group_plot(gridded_vars):
    """Set of plots of the gridded variables.

    Inputs :
        gridded_vars - A gridded_vars object from map_fun.py

    Returns :
        Draws a figure that can be rendered with plt.show()
    """
    # make figure
    plt.figure()
    plt.subplot(2, 3, 1)
    plt.imshow(gridded_vars.ex)
    plt.colorbar()
    plt.title('Ebb x')

    plt.subplot(2, 3, 2)
    plt.imshow(gridded_vars.ey)
    plt.colorbar()
    plt.title('Ebb y')

    plt.subplot(2, 3, 3)
    plt.imshow(np.sqrt(gridded_vars.ex**2+gridded_vars.ey**2))
    plt.colorbar()
    plt.title('Ebb Mag')

    plt.subplot(2, 3, 4)
    plt.imshow(gridded_vars.fx)
    plt.colorbar()
    plt.title('Flood x')

    plt.subplot(2, 3, 5)
    plt.imshow(gridded_vars.fy)
    plt.colorbar()
    plt.title('Flood y')

    plt.subplot(2, 3, 6)
    plt.imshow(np.sqrt(gridded_vars.fx**2+gridded_vars.fy**2))
    plt.colorbar()
    plt.title('Flood Mag')

    plt.tight_layout()
