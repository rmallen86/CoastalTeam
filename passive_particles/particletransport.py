"""Functions associated with the passive particle routing."""
from dorado import particle_track as pt
from dorado.routines import plot_state
import matplotlib.pyplot as plt


def init_particles(init_x, init_y, Np_tracer, grid_spacing, gridded_vals):
    """Initialize the particles.

    Inputs :
        init_x - List of starting x locations for the particles
        init_y - List of starting y locations for the particles
        Np_tracer - Number of particles to route
        grid_spacing - Size of the raster grid cells
        gridded_vals - gridded_vals object from map_fun.py

    Outputs :
        params - an initialized params class

    """
    # create class
    params = pt.params()
    # populate with required parameters
    params.seed_xloc = init_y
    params.seed_yloc = init_x
    params.u = gridded_vals.ex
    params.v = gridded_vals.ey
    params.depth = gridded_vals.depth
    params.topography = gridded_vals.elev
    params.Np_tracer = Np_tracer
    params.dx = grid_spacing
    # hold ebb and flood velocity components in the class too
    params.ex = gridded_vals.ex
    params.ey = gridded_vals.ey
    params.fx = gridded_vals.fx
    params.fy = gridded_vals.fy

    return params


def tidal_particles(params, tide_period, n_tide_periods):
    """Route the particles in tides.

    Inputs :
        params - params output from init_particles()
        tide_period - tidal period in seconds
        n_tide_periods - number of tidal periods to iterate particles over

    Returns :
        walk_data - history of particle locations and travel times

        Also saves image of particle locations to disk for each flood/ebb tide.

    """
    # define the particle
    particle = pt.Particle(params)
    # record each 1/2 tidal cycle so each ebb and flood
    for i in range(0, int(2*n_tide_periods)):
        if i == 0:
            # start with ebb tide
            walk_data = particle.run_iteration()
        else:
            if i % 2 != 0:
                # ebb tide
                params.u = params.ex
                params.v = params.ey
                particle = pt.Particle(params)
                walk_data = particle.run_iteration(previous_walk_data=walk_data,
                                                   target_time=tide_period/2*(i+1))
            else:
                # flood tide
                params.u = params.fx
                params.v = params.fy
                particle = pt.Particle(params)
                walk_data = particle.run_iteration(previous_walk_data=walk_data,
                                                   target_time=tide_period/2*(i+1))

        # plot and save particle locations
        plot_state(particle.depth, walk_data, -1, None, 'r')
        if i % 2 != 0:
            plt.title('Ebb Tide, Time = ' + str(tide_period/2*(i+1)) + 's')
        else:
            plt.title('Flood Tide, Time = ' + str(tide_period/2*(i+1)) + 's')
        plt.savefig(str(i) + '.png')
        plt.close()
