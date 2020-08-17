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
import numpy as np
import scipy.io as spio
import excludeboundarycell

def sedtran(d, A, DiffS, h, ho, E, WS, dx, dt, rbulk, co, Ux, Uy, FLX, fTide, Ttide):
    
    # DELETE BELOW THIS WHEN FINISHED

    mat = spio.loadmat('C:/Users/colli/Documents/Python_Scripts/ESPIn/coastal/MarshMorpho2D/rightbeforeTime.mat')
    A = mat['A']
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
    
    find_A = np.argwhere(ho<=0)
    A[find_A] = 0
    
    # rivermouthfront
    
    p = np.argwhere(A>0); # exclude the no land cells
    G = 0*d
    NN=len(p)
    G[p] =[1:NN]
    rhp = E[p] # in the rhs there is already the addition of the erosion input
    [N, M] = size(G)
    i =(); j=();s=();S=0*G;
    
    # boundary conditions imposed SSC
    a = np.argwhere(A==2)
    rhs[G[a]] = co*h[a]*fTide[a]
    
    Dxx = (DiffS*Ttide/2*(np.abs(Ux*Ux))*(24*3600)^2)*h #.*(ho>kro);%.*(hgross>0.01);%% the h is not the coefficient for diffusion, it the h in the mass balance eq.
    Dyy = (DiffS*Ttide/2*(np.abs(Ux*Ux))*(24*3600)^2)*h
    
    # the factor 24*3600 is used to convert the Ux and Uy from m/s to m/day
    
    [row, col] = np.unravel_index(p, size(A))
    for k in [N -1 1 -N]:
        [a, q] = exludeboundarycell(k, N, M, p)
        
        a =a[A[q]]