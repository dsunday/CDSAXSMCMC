# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 16:32:00 2016

@author: dfs1
"""

import numpy as np
import CDSAXSfunctions as CD
from multiprocessing import Pool
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
MCPAR[0] = 24 # Chainnumber
MCPAR[1] = len(FITPAR)
MCPAR[2] = 20 #stepnumber
MCPAR[3] = 0 #randomchains
MCPAR[4] = 1 # Resampleinterval
MCPAR[5] = 100 # stepbase
MCPAR[6] = 200 # steplength

MCMCInitial=np.load('SAt3N_5T_Init.npy')

MCMC_List=[0]*int(MCPAR[0])
for i in range(int(MCPAR[0])):
    MCMC_List[i]=MCMCInitial[i,:]
    
def MCMC_PSPVP(MCMC_List):
    
    MCMCInit=MCMC_List
    
    L = int(MCPAR[1])
    Stepnumber= int(MCPAR[2])
        
    SampleMatrix=np.zeros([Stepnumber,L+1]) 
    SampleMatrix[0,:]=MCMCInit
    Move = np.zeros([L+1])
    
    ChiPrior = MCMCInit[L]
    for step in np.arange(1,Stepnumber,1): 
        Temp = SampleMatrix[step-1,:]
        for p in range(L-3):
            StepControl = MCPAR[5]+MCPAR[6]*np.random.random_sample()
            Move[p] = (FITPARUB[p]-FITPARLB[p])/StepControl*(np.random.random_sample()-0.5) # need out of bounds check
        Temp[0:L-3]=Temp[0:L-3]+Move[0:L-3]
        
        SimPost=SimInt_PSPVP(Temp)
        ChiPost=np.sum(CD.Misfit(Intensity,SimPost))
        if ChiPost < ChiPrior:
            SampleMatrix[step,0:L]=Temp[0:L]
            SampleMatrix[step,L]=ChiPost
            ChiPrior=ChiPost
            
        else:
            MoveProb = np.exp(-0.5*np.power(ChiPost-ChiPrior,2))
            if np.random.random_sample() < MoveProb:
                SampleMatrix[step,0:L]=Temp[0:L]
                SampleMatrix[step,L]=ChiPost
                ChiPrior=ChiPost
            else:
                SampleMatrix[step,:]=SampleMatrix[step-1,:]
                
    return SampleMatrix
    
def SimInt_PSPVP(FITPAR):
    T=int(FITPAR[len(FITPAR)-3])
    Disc=int(FITPAR[len(FITPAR)-2])
    Pitch=int(FITPAR[len(FITPAR)-4])
    Spline=np.reshape(FITPAR[0:T*5],(T,5))
    Spline
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
    M=np.power(np.exp(-1*(np.power(Qx,2)+np.power(Qz,2))*np.power(SPAR[0],2)),0.5);
    Formfactor=Formfactor*M
    SimInt = np.power(abs(Formfactor),2)*SPAR[1]+SPAR[2]
    return SimInt



start_time = time.perf_counter()
if __name__ =='__main__':  
    pool = Pool(processes=4)
    F=pool.map(MCMC_PSPVP,MCMC_List)
    F=tuple(F)
    np.save('Slurmtest',F)
    end_time=time.perf_counter()   
    print(end_time-start_time)    