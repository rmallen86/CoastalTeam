"""Functions associated with the passive particle routing."""
from dorado import particle_track as pt
from dorado.routines import plot_state
import matplotlib.pyplot as plt


def init_particles(init_x, init_y, Np_tracer, grid_spacing, gridded_vals):
    """Initialize the particles.

    Inputs :
        init_x : `list`
            List of starting x locations for the particles

        init_y : `list`
            List of starting y locations for the particles

        Np_tracer : `int`
            Number of particles to route

        grid_spacing : `int`
            Size of the raster grid cells

        gridded_vals : `obj`
            gridded_vals object from map_fun.py

    Outputs :
        params : `obj`
            an initialized params class

    """
    # create class
    params = pt.modelParams()
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


def tidal_particles(params, tide_period, n_tide_periods, plot_grid=None):
    """Route the particles in tides.

    Inputs :
        params : `obj`
            params output from init_particles()

        tide_period : `int`
            tidal period in seconds

        n_tide_periods : 'int'
            number of tidal periods to iterate particles over

        plot_grid : `numpy.ndarray` (Optional)
            grid on which to plot the particles (e.g. depth)

    Returns :
        walk_data : `list`
            history of particle locations and travel times

        Also saves image of particle locations to disk for each flood/ebb tide.

    """
    # define the particle
    particle = pt.Particles(params)
    # generate a set of particles to route around
    particle.generate_particles(params.Np_tracer,
                                params.seed_xloc,
                                params.seed_yloc)
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
                particle = pt.Particles(params)
                particle.generate_particles(0, [], [],
                                            previous_walk_data=walk_data)
                walk_data = particle.run_iteration(target_time=tide_period/2*(i+1))
            else:
                # flood tide
                params.u = params.fx
                params.v = params.fy
                particle = pt.Particles(params)
                particle.generate_particles(0, [], [],
                                            previous_walk_data=walk_data)
                walk_data = particle.run_iteration(target_time=tide_period/2*(i+1))

        # plot and save particle locations
        if plot_grid is None:
            plot_state(particle.depth, walk_data, -1, None, 'r')
        else:
            plot_state(plot_grid, walk_data, -1, None, 'r')

        # set colorbar
        plt.colorbar()

        if i % 2 != 0:
            plt.title('Ebb Tide, Time = ' + str(tide_period/2*(i+1)) + 's')
        else:
            plt.title('Flood Tide, Time = ' + str(tide_period/2*(i+1)) + 's')

        plt.tight_layout()
        plt.savefig(str(i) + '.png')
        plt.close()

    return walk_data
