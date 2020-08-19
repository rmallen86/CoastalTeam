"""
Total sediment erosion (assumes limitless supply)
Based on totalsedimenterosionmudsine.m by Giulio Mariotti 
Includes some code based on G. Tucker tidal_flow_calculator.py in landlab component
Does not incorporate averaging over multiple tidal cycles or linear increase in critical shear stress with depth
"""

# imports
import numpy as np
import matplotlib.pyplot as plt
from landlab import RasterModelGrid, imshow_grid
from landlab.components import TidalFlowCalculator
from landlab.io import read_esri_ascii
from landlab.grid.mappers import map_mean_of_link_nodes_to_link, map_node_to_cell, map_link_vector_components_to_node, map_min_of_node_links_to_node


def map_velocity_components_to_nodes(grid):
    """Map the velocity components from the links to the nodes, and return the node arrays."""
    ebb_vel_x, ebb_vel_y = map_link_vector_components_to_node(grid, grid.at_link['ebb_tide_flow__velocity'])
    flood_vel_x = -ebb_vel_x
    flood_vel_y = -ebb_vel_y
    ebb_vel = np.sqrt(ebb_vel_x*ebb_vel_x + ebb_vel_y*ebb_vel_y)
    flood_vel = -ebb_vel
    return (ebb_vel_x, ebb_vel_y, flood_vel_x, flood_vel_y, ebb_vel, flood_vel)

def plot_tidal_flow(grid, resample=1):
    (ebb_x, ebb_y, flood_x, flood_y, ebb, flood) = map_velocity_components_to_nodes(grid)

    # depth
    plt.figure()
    imshow_grid(grid, grid.at_node['mean_water__depth'], cmap='YlGnBu', color_for_closed='g')
    plt.title('Water depth (m)')
    plt.xlabel('Distance (m)')
    plt.ylabel('Distance (m)')

    # down-sample for legible quiver plots if needed
    if resample != 1:
        xr = grid.x_of_node.reshape((grid.number_of_node_rows, grid.number_of_node_columns))[::resample, ::resample]
        yr = grid.y_of_node.reshape((grid.number_of_node_rows, grid.number_of_node_columns))[::resample, ::resample]
        ebb_xr = ebb_x.reshape((grid.number_of_node_rows, grid.number_of_node_columns))[::resample, ::resample]
        ebb_yr = ebb_y.reshape((grid.number_of_node_rows, grid.number_of_node_columns))[::resample, ::resample]
        fld_xr = flood_x.reshape((grid.number_of_node_rows, grid.number_of_node_columns))[::resample, ::resample]
        fld_yr = flood_y.reshape((grid.number_of_node_rows, grid.number_of_node_columns))[::resample, ::resample]
    else:
        xr = grid.x_of_node
        yr = grid.y_of_node
        ebb_xr = ebb_x
        ebb_yr = ebb_y
        fld_xr = flood_x
        fld_yr = flood_y
        
    # ebb tide
    plt.figure()
    imshow_grid(grid, grid.at_node['topographic__elevation'])
    plt.quiver(xr, yr, ebb_xr, ebb_yr)
    plt.title('Ebb Tide')
    plt.xlabel('Distance (m)')
    plt.ylabel('Distance (m)')

    ebb_vel_magnitude = ebb
    plt.figure()
    imshow_grid(grid, ebb_vel_magnitude, cmap='magma', color_for_closed='g')
    plt.title('Ebb Tide Velocity Magnitude (m/s)')
    plt.xlabel('Distance (m)')
    plt.ylabel('Distance (m)')

    # flood tide
    plt.figure()
    imshow_grid(grid, grid.at_node['topographic__elevation'])
    plt.quiver(xr, yr, fld_xr, fld_yr)
    plt.title('Flood Tide')
    plt.xlabel('Distance (m)')
    plt.ylabel('Distance (m)')

    plt.figure()
    flood_vel_magnitude = flood
    imshow_grid(grid, flood_vel_magnitude, cmap='magma', color_for_closed='g')
    plt.title('Flood Tide Velocity Magnitude (m/s)')
    plt.xlabel('Distance (m)')
    plt.ylabel('Distance (m)')
    
def map_node2cell_addGrid(grid,var1,var2): #takes a grid, plus the variable you want to map to cell, and the string name
    a = grid.map_node_to_cell(var1)
    return grid.add_field(var2, a, at='cell',clobber=True)
    
def map_link2cell_addGrid(grid,var1,var2): #takes a grid, plus the variable you want to map to cell, and the string name
    a = grid.map_min_of_node_links_to_node(var1)
    b = grid.map_node_to_cell(a)
    return grid.add_field(var2, b, at='cell',clobber=True)
    
def populateGrids(grid, tfc, tau_cr, tau_crv, veg):
    rate = tfc.calc_tidal_inundation_rate()
    grid.add_field('tidal_innundation_rate',rate,at = 'node',units='m/s',clobber=True)
    map_node2cell_addGrid(grid,rate,'tidal_innundation_rate_cell')

    tfc._calc_effective_water_depth()
    grid.add_field('effective_water_depth',tfc._water_depth,at='node',units='m',clobber=True)
    map_node2cell_addGrid(grid,tfc._water_depth,'effective_water_depth_cell')
    
    tfc.run_one_step()

    topo = grid.at_node['topographic__elevation']
    map_node2cell_addGrid(grid,topo,'topographic_elevation_cell')
    
    msl = tfc._mean_sea_level
    dHW = np.maximum(0,topo+msl+tfc._tidal_half_range)
    dHW[topo==999] = 0;

    ftide = np.minimum(1,np.maximum(10^-3, dHW/tfc._tidal_range))
    #ftide[topo==999] = 999
    grid.add_field('hydroperiod',ftide,at='node',units='m',clobber=True)
    grid.add_field('water_depth_at_MHW',dHW,at='node',units='m',clobber=True)
    map_node2cell_addGrid(grid,ftide,'hydroperiod_cell')
    map_node2cell_addGrid(grid,dHW,'water_depth_at_MHW_cell')
    
    lev_an = -topo-msl #water depth with respect to MSL
    grid.add_field('lev_at_node',lev_an,at = 'node',clobber=True)
    lev_atlink = grid.map_mean_of_link_nodes_to_link('lev_at_node')
    map_node2cell_addGrid(grid,lev_an,'lev_at_cell')

    taucr = grid.add_zeros('tau_cr',at='link') + tau_cr
    v = grid.at_link['veg_atlink'] 
    taucr[v==1] = tau_crv
    taucr_node = grid.map_min_of_node_links_to_node(taucr)
    grid.add_field('tau_cr_node',taucr_node,at='node',clobber=True)
    map_link2cell_addGrid(grid,taucr,'tau_cr_cell')
    
    ebb = grid.at_link['ebb_tide_flow__velocity']
    ebb_node = grid.map_min_of_node_links_to_node(ebb)
    grid.add_field('ebb_tide_flow__velocity_node',ebb_node,at='node',clobber=True)
    grid.add_field('flood_tide_flow__velocity_node',-ebb_node,at='node',clobber=True)
    map_node2cell_addGrid(grid,ebb,'ebb_tide_flow__velocity_cell')
    map_node2cell_addGrid(grid,-ebb,'flood_tide_flow__velocity_cell')
    
    rough = grid.at_link['roughness']
    rough_node = grid.map_min_of_node_links_to_node(rough)
    grid.add_field('roughness_node',rough_node,at='node',clobber=True)
    map_link2cell_addGrid(grid,rough,'roughness_cell')
    
    grid.add_field('water_depth_at_link',tfc._water_depth_at_links,at = 'link',units='m',clobber=True)
    wd = tfc._water_depth_at_links
    wd_node = grid.map_min_of_node_links_to_node(wd)
    grid.add_field('water_depth_at_node',wd_node,at='node',clobber=True)
    map_link2cell_addGrid(grid,tfc._water_depth_at_links,'water_depth_at_cell')

def totalsedimenterosion_mudsine(grid, mud_erodability,tr,tcg):

    fupeak = np.pi/2
    #total sed erosion for loop
    ntdcy = 10 #number of tidal cycles
    taucr = grid.at_node['tau_cr_node']
    E = np.zeros(taucr.size)
    
    # lev = grid.at_node['lev_at_node']
    # xi = -lev-tr/2
    # xi[xi<0] = 0
    # taucr += xi*tcg
    
    utide = grid.at_node['flood_tide_flow__velocity_node']*fupeak*np.sin(np.pi/2) #intra-tidal velocity
    rough = grid.at_node['roughness_node']
    h = grid.at_node['water_depth_at_node']

    tauC = 1025*9.81* (rough**2) * (utide**2) * (h**(-1/3))
    E += mud_erodability*(np.sqrt(1+(tauC/taucr)**2)-1)
    
    grid.add_field('erosion',E,at='node',clobber=True)
    grid.add_field('utide',utide,at='node',clobber=True)
    grid.add_field('tauC',tauC,at='node',clobber=True)
    return E
    
def totalsedimenterosion_mudsine_link(grid, mud_erodability):

    fupeak = np.pi/2
    #total sed erosion for loop
    ntdcy = 10 #number of tidal cycles
    E = grid.add_zeros('Erosion',at = 'cell')
    taucr = grid.at_cell['tau_cr_cell']
    
    # get rid of for loop
    #for i in range(ntdcy):
    utide = grid.at_cell['flood_tide_flow__velocity_cell']*fupeak*np.sin(np.pi/2) #intra-tidal velocity
    print(utide.mean())
    tauC = 1025*9.81* grid.at_cell['roughness_cell']**2 * utide**2 * grid.at_cell['water_depth_at_cell']**(-1/3)
    E += mud_erodability*(np.sqrt(1+(tauC/taucr)**2)-1)
    #print(max(E))

def updategrids(grid, tfc): #need to update grids used in totalsedimenterosion_mudsine that were changed by tidal_flow_calculator
    ebb = grid.at_link['ebb_tide_flow__velocity']
    ebb_node = grid.map_min_of_node_links_to_node(ebb)
    grid.add_field('ebb_tide_flow__velocity_node',ebb_node,at='node',clobber=True)
    grid.add_field('flood_tide_flow__velocity_node',-ebb_node,at='node',clobber=True)
    
    rough = grid.at_link['roughness']
    rough_node = grid.map_min_of_node_links_to_node(rough)
    grid.add_field('roughness_node',rough_node,at='node',clobber=True)
    
    grid.add_field('water_depth_at_link',tfc._water_depth_at_links,at = 'link',units='m',clobber=True)
    wd = tfc._water_depth_at_links
    wd_node = grid.map_min_of_node_links_to_node(wd)
    grid.add_field('water_depth_at_node',wd_node,at='node',clobber=True)
    
    msl = tfc._mean_sea_level
    topo = grid.at_node['topographic__elevation']
    lev_an = -topo-msl #water depth with respect to MSL
    grid.add_field('lev_at_node',lev_an,at = 'node',clobber=True)
    
    
