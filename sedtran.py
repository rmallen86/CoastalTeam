# -*- coding: utf-8 -*-
"""
Sedtran -- sediment transport to evolve the bed surface
Based on sedtran.m by ______ (github link ____)

Created on Mon Aug 17 2020

@author: Collin Roland ...
% Sediment transport to evolve the bed surface
% Input:
    % d: 
    % A: cell types (0: out of domain, 1: normal, 2: BC cells)
    % DiffS: coefficient for tidal dispersion [-]
    % h: tidal averaged water depth
    % ho: 
    % E: 
    % WS: 
    % dx: cell size [m]
    % dt: time step size [days]
    % rbulk: dry bulk density [kg/m3]
    % co: sea boundary suspended sediment concentration for mud [g/L]
    % Ux: water velocity in the x-direction
    % Uy: water velocity in the y-direction
    % FLX: sediment flux through open boundary
    % fTide: hydroperiod
    % Ttide: tidal period [day]
    % kro: minimum water depth [m]
% Output:
    % EmD: Values of erosion/deposition with which to update the bed
    % SSM: suspended sediment mass ???
    % FLX: updated sediment flux through open boundary
"""

# imports
import numpy as np
import scipy.io as spio
#import excludeboundarycell

def sedtran(d, A, DiffS, h, ho, E, WS, dx, dt, rbulk, co, Ux, Uy, FLX, fTide, Ttide):
    
    # DELETE BELOW THIS WHEN FINISHED

    mat = spio.loadmat('C:/Users/colli/Documents/Python_Scripts/ESPIn/coastal/MarshMorpho2D/rightbeforeTime.mat')
    A_int = mat['A']
    A = A_int.astype(float)
    d = np.zeros_like(A) # creates an array of zeros the shape of A (in reality this is an array of floats)
    d = d.astype(float)
    DiffS = mat['DiffS']
    h = mat['h']
    ho = mat['ho']
    E = mat['E2'] # there's no E var
    WS = mat['WS']
    dx = mat['dx']
    dt = mat['dt']
    rbulk = mat['rbulk2'] # there's no rbulk var
    co = mat['co2'] # there's no co var
    Ux = mat['Ux']
    Uy = mat['Uy']
    FLX = mat['FLX']
    fTide = mat['fTide']
    Ttide = mat['Ttide']
    # DELETE ABOVE THIS WHEN FINISHED
    
    # Eliminate the cells in which the water depth is too small
    find_A = np.argwhere(ho<=0)
    A[find_A[:,0],find_A[:,1]]= 0 # this indexing doesn't quite work    
    
    # rivermouthfront
    p = np.argwhere(A.flatten(order='F')>0)
    G = 0*d
    G_shape = np.shape(G)
    NN=len(p)
    G2 = G.flatten(order='F')
    G2 = G2.astype(int)
    p2 = np.arange(1,NN+1)
    p3 = p2.reshape((86628,1))
    G2[p] = p3
    G = G2.reshape(G_shape,order='F')

    E2 = E.flatten(order='F')
    
    rhs = E2[p] # in the rhs there is already the addition of the erosion input
    [N, M] = np.shape(G)
    i =np.empty(0); j=np.empty(0);s=np.empty(0);
    S=0*G;
    
    # boundary conditions imposed SSC
    a = np.argwhere(A.flatten(order='F')==2)
    rhs_ind = G2[a]
    rhs_ind = rhs_ind.reshape(15,)
    rhs[rhs_ind] = co*(h.flatten(order='F')[a])*(fTide.flatten(order='F')[a])
    
    Dxx = (DiffS*Ttide/2*(np.abs(Ux*Ux))*(24.0*3600.0)**2.0)*h #.*(ho>kro);%.*(hgross>0.01);%% the h is not the coefficient for diffusion, it the h in the mass balance eq.
    Dyy = (DiffS*Ttide/2*(np.abs(Ux*Ux))*(24.0*3600.0)**2.0)*h
    
    # the factor 24*3600 is used to convert the Ux and Uy from m/s to m/day
    
    [row, col] = np.unravel_index(p, size(A))
    for k in [N -1 1 -N]:
        [a, q] = exludeboundarycell(k, N, M, p)
        
        a =a[A[q]]