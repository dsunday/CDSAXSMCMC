# -*- coding: utf-8 -*-
"""
This program will be used to analyse CDSAXS data from qz slices, including MCMC program for uncertainty analysis
"""
import numpy as np
import CDSAXSfunctions as CD
import matplotlib.pyplot as plt
[Qz,Qx,Intensity]=CD.CDSAXSimport('test.txt','212SCANQx.txt',121,26)
Intensity=Intensity*1000000
Trapnumber = 6
Discretization = 6
tpar=np.zeros([Trapnumber+1,3])

tpar[0,0] = 35; tpar[0,1] = 37; tpar[0,2] = 1;
tpar[1,0] = 24; tpar[1,1] = 40; tpar[1,2] = 1;
tpar[2,0] = 20; tpar[2,1] = 50; tpar[2,2] = 1;
tpar[3,0] = 15; tpar[3,1] = 50; tpar[3,2] = 1;
tpar[4,0] = 13; tpar[4,1] = 7; tpar[4,2] = 1;
tpar[5,0] = 12; tpar[5,1] = 15; tpar[5,2] = 1;
tpar[6,0] = 10; tpar[6,1] = 0; tpar[6,2] = 1;
X1 = 32; X2 = 35; X3 = 33;
X=[X1, X2]
ppar=np.zeros([3,3])
ppar[0,0]=0.5;  ppar[0,1]=0.2;  ppar[0,2]=0;
ppar[1,0]=0.5;  ppar[1,1]=0.2;  ppar[1,2]=0;
ppar[2,0]=0.5;  ppar[2,1]=0.2;  ppar[2,2]=0;

Pitch = 135.7;

DW = 1.5
I0 = 0.005
Bk = 33.4

SPAR = [DW,I0, Bk, Pitch]
Coord=CD.SCNCoordAssign(tpar,Trapnumber,X1,X2,X3,Pitch)
SimInt = CD.SCNIntensitySim(Coord,Qx,Qz,Trapnumber,DW,I0,Bk)
Coordp=CD.SCNParabolaCoord(tpar,ppar,Discretization,Trapnumber, X1, X2, Pitch)
SimIntp = CD.SCNIntensitySim(Coordp,Qx,Qz,Trapnumber+Discretization-1,DW,I0,Bk)

CD.PlotQzCut(Qz,SimIntp,Intensity,26)

Chi2=CD.Misfit(Intensity,SimInt)
C=Chi2.sum()