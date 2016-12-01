# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 08:55:11 2016

@author: dfs1
"""
# import libraries
import numpy as np
import CDSAXSfunctions as CD
import matplotlib.pyplot as plt

# Import data
Intensity=np.loadtxt('Sat3N_Intensity.txt')
Qx = np.loadtxt('Sat3N_Qx.txt')
Qz = np.loadtxt('Sat3N_Qz.txt')
CutLengths=np.loadtxt('Sat3N_CutLength.txt')

Trapnumber=5
Discretization=18

Spline=np.zeros([Trapnumber,5])
Pitch = 84 
Offset = np.zeros([7,1])
Offset[0,0]=20
Offset[1,0]=5
Offset[2,0]=10.5
Offset[3,0]=10.5
Offset[4,0]=10.5
Offset[5,0]=10.5
Offset[6,0]=10

Spline[0,0]=0;  Spline[0,1]=0; Spline[0,2]=0; Spline[0,3]=0; Spline[0,4]=0;
Spline[1,0]=5;  Spline[1,1]=3; Spline[1,2]=0; Spline[1,3]=0; Spline[1,4]=0;
Spline[2,0]=15;  Spline[2,1]=3; Spline[2,2]=1; Spline[2,3]=1; Spline[2,4]=1;
Spline[3,0]=21;  Spline[3,1]=3; Spline[3,2]=-1; Spline[3,3]=-0.5; Spline[3,4]=-1;
Spline[4,0]=24;  Spline[4,1]=3; Spline[4,2]=-1; Spline[4,3]=-0.5; Spline[4,4]=-1;

MCoord=np.loadtxt('MCOORDPSPVP1.txt')
(Coord)= CD.PSPVPCoord(20,Spline,MCoord,Trapnumber, Discretization,Pitch, Offset)
# Define spline parameters