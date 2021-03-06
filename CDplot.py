# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 13:25:57 2016

@author: dfs1
"""
import matplotlib.pyplot as plt
import numpy as np
import CDSAXSfunctions as CD 
   
def PlotQzCut(Qz,SimInt,ExpI,numbercuts):
    S=SimInt
    I=ExpI    
    for i in range(0,numbercuts):
#        IMin=np.min(SimInt[:,i-1])
#        IMax=np.max(SimInt[:,i])
#        R=2*IMax/IMin
        S[:,i]=S[:,i]/(10.**(i+1))
        I[:,i]=I[:,i]/(10.**(i+1))
    for i in range(numbercuts):
        plt.semilogy(Qz[:,i],I[:,i],'.')
        plt.semilogy(Qz[:,i],S[:,i])
        
def plotPSPVP(Coord,Trapnumber):
    for S in range(8):
        h=0
        Lc= np.zeros([Trapnumber+1,2])
        Rc= np.zeros([Trapnumber+1,2])
        
        for i in range(Trapnumber+1):
            Lc[i,0]=Coord[i,0,S]
            Rc[i,0]=Coord[i,1,S]
            Lc[i,1]=h
            Rc[i,1]=h
            h=h+Coord[i,2,S]
        plt.plot(Lc[:,0],Lc[:,1])
        plt.plot(Rc[:,0],Rc[:,1])
        Cc=np.zeros([2,2])
        for i in range(Trapnumber):
            Cc[0,0]=Lc[i+1,0]
            Cc[0,1]=Lc[i+1,1]
            Cc[1,0]=Rc[i+1,0]
            Cc[1,1]=Rc[i+1,1]
            plt.plot(Cc[:,0],Cc[:,1])
            
    
