# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 15:03:39 2016

@author: dfs1
"""

import numpy as np
import CDSAXSfunctions as CD
import CDplot as CDp
import matplotlib.pyplot as plt
import time
Intensity=np.loadtxt('Sat3N_Intensity.txt')
Qx = np.loadtxt('Sat3N_Qx.txt')
Qz = np.loadtxt('Sat3N_Qz.txt')
CutLengths=np.loadtxt('Sat3N_CutLength.txt')

Trapnumber = 5
Disc=18

DW = 1.69
I0 = 0.2
Bk = 1.2
Pitch = 84
SPAR=np.zeros([4]) 
SPAR[0]=DW; SPAR[1]=I0; SPAR[2]=Bk;SPAR[3]=Pitch;

Spline=np.zeros([Trapnumber,5])

Offset = np.zeros([7,1])
Offset[0,0]=20
Offset[1,0]=5
Offset[2,0]=10.5
Offset[3,0]=10.5
Offset[4,0]=10.5
Offset[5,0]=10.5
Offset[6,0]=10

Spline[0,0]=0;  Spline[0,1]=0; Spline[0,2]=0; Spline[0,3]=0; Spline[0,4]=0;
Spline[1,0]=2.4;  Spline[1,1]=1; Spline[1,2]=0; Spline[1,3]=0; Spline[1,4]=0;
Spline[2,0]=6;  Spline[2,1]=3; Spline[2,2]=0; Spline[2,3]=0; Spline[2,4]=0;
Spline[3,0]=19;  Spline[3,1]=3; Spline[3,2]=0; Spline[3,3]=0; Spline[3,4]=0;
Spline[4,0]=24;  Spline[4,1]=2; Spline[4,2]=0; Spline[4,3]=0; Spline[4,4]=0;

MCoord=np.loadtxt('MCOORDPSPVP1.txt')

(Coord)= CD.PSPVPCoord(Spline,MCoord,Trapnumber, Disc,Pitch, Offset)
(FITPAR,FITPARLB,FITPARUB)=CD.PSPVP_PB(Offset,Spline,SPAR,Trapnumber,Disc)


MCPAR=np.zeros([7])
MCPAR[0] = 2000 # Chainnumber
MCPAR[1] = len(FITPAR)
MCPAR[2] = 100 #stepnumber
MCPAR[3] = 0 #randomchains
MCPAR[4] = 1 # Resampleinterval
MCPAR[5] = 40 # stepbase
MCPAR[6] = 200 # steplength

def MCMCInit_PSPVP(FITPAR,FITPARLB,FITPARUB,MCPAR):
    
    MCMCInit=np.zeros([int(MCPAR[0]),int(MCPAR[1])+1])
    
    for i in range(int(MCPAR[0])):
        if i <MCPAR[3]: #reversed from matlab code assigns all chains below randomnumber as random chains
            for c in range(int(MCPAR[1])-3):
                MCMCInit[i,c]=FITPARLB[c]+(FITPARUB[c]-FITPARLB[c])*np.random.random_sample()
            MCMCInit[i,int(MCPAR[1])-3:int(MCPAR[1])]=FITPAR[int(MCPAR[1])-3:int(MCPAR[1])]
            SimInt=SimInt_PSPVP(MCMCInit[i,:])
            C=np.sum(CD.Misfit(Intensity,SimInt))
            
            MCMCInit[i,int(MCPAR[1])]=C
            
        else:
            MCMCInit[i,0:int(MCPAR[1])]=FITPAR
            SimInt=SimInt_PSPVP(MCMCInit[i,:])
            C=np.sum(CD.Misfit(Intensity,SimInt))
            MCMCInit[i,int(MCPAR[1])]=C
           
    return MCMCInit
    
def MCMCInit_PSPVPUniform(FITPAR,FITPARLB,FITPARUB,MCPAR):
    
    MCMCInit=np.zeros([int(MCPAR[0]),int(MCPAR[1])+1])
    
    for i in range(int(MCPAR[1])-3):
        if FITPARUB[i]==FITPARLB[i]:
            MCMCInit[:,i]=FITPAR[i]
        else:
            A= np.arange(FITPARLB[i],FITPARUB[i]+0.0001,(FITPARUB[i]-FITPARLB[i])/(int(MCPAR[0])-1))
            R=np.random.rand(int(MCPAR[0]))
            ind=R.argsort()
            A=A[ind]
            MCMCInit[:,i]=A
    MCMCInit[:,int(MCPAR[1])-3:int(MCPAR[1])]=FITPAR[int(MCPAR[1])-3:int(MCPAR[1])]        
    
    for i in range(int(MCPAR[0])):
       SimInt=SimInt_PSPVP(MCMCInit[i,:])
       C=np.sum(CD.Misfit(Intensity,SimInt))
       MCMCInit[i,int(MCPAR[1])]=C
        
    return MCMCInit
    
def SimInt_PSPVP(FITPAR):
    T=int(FITPAR[len(FITPAR)-3])
    Disc=int(FITPAR[len(FITPAR)-2])
    Pitch=int(FITPAR[len(FITPAR)-4])
    Spline=np.reshape(FITPAR[0:T*5],(T,5))
    Offset=FITPAR[T*5:T*5+7]
    SPAR=FITPAR[T*5+7:T*5+11]
    Coord=CD.PSPVPCoord(Spline,MCoord,Trapnumber, Disc,Pitch, Offset)
    F1 = CD.FreeFormTrapezoid(Coord[:,:,0],Qx,Qz,Disc+1)
    F2 = CD.FreeFormTrapezoid(Coord[:,:,1],Qx,Qz,Disc+1)
    F3 = CD.FreeFormTrapezoid(Coord[:,:,2],Qx,Qz,Disc+1)
    F4 = CD.FreeFormTrapezoid(Coord[:,:,3],Qx,Qz,Disc+1)
    F5 = CD.FreeFormTrapezoid(Coord[:,:,4],Qx,Qz,Disc+1)
    F6 = CD.FreeFormTrapezoid(Coord[:,:,5],Qx,Qz,Disc+1)
    F7 = CD.FreeFormTrapezoid(Coord[:,:,6],Qx,Qz,Disc+1)
    F8 = CD.FreeFormTrapezoid(Coord[:,:,7],Qx,Qz,Disc+1)    
    Formfactor=(F1+F2+F3+F4+F5+F6+F7+F8)
    M=np.power(np.exp(-1*(np.power(Qx,2)+np.power(Qz,2))*np.power(SPAR[0],2)),0.5)
    Formfactor=Formfactor*M
    Formfactor=abs(Formfactor)
    SimInt = np.power(Formfactor,2)*SPAR[1]+SPAR[2]
    return SimInt
    
start_time = time.perf_counter()
MCMCInitialU=MCMCInit_PSPVPUniform(FITPAR,FITPARLB,FITPARUB,MCPAR)
end_time=time.perf_counter()   
print(end_time-start_time)
MCMCInitialU=MCMCInitialU[MCMCInitialU[:,int(MCPAR[1])].argsort(),]

A=0;
Spline=np.reshape(MCMCInitialU[A,0:Trapnumber*5],(Trapnumber,5))
Offset=MCMCInitialU[A,Trapnumber*5:Trapnumber*5+7]
SPAR=MCMCInitialU[A,Trapnumber*5+7:Trapnumber*5+11]
Coord=CD.PSPVPCoord(Spline,MCoord,Trapnumber, Disc,Pitch, Offset)
plt.figure(1)
CDp.plotPSPVP(Coord,Disc)
Sim=SimInt_PSPVP(MCMCInitialU[A,:])
plt.figure(2)
plt.semilogy(Qz[:,0],Intensity[:,0],'.')
plt.semilogy(Qz[:,0],Sim[:,0])

plt.semilogy(Qz[:,3],Intensity[:,3],'.')
plt.semilogy(Qz[:,3],Sim[:,3])

C=np.sum(CD.Misfit(Intensity,Sim))
#CDp.PlotQzCut(Qz,Sim,Intensity,11)