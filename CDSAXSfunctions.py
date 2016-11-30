# -*- coding: utf-8 -*-
"""
Functions to be used in analyzing CDSAXS data
"""
import numpy as np
import matplotlib.pyplot as plt
import CDSAXSfunctions as CD

def CDSAXSimport(filenameqz,filenameqx,AngleNumber,cutnumber):
    Qz=np.zeros([AngleNumber,cutnumber])
    Qx=np.zeros([AngleNumber,cutnumber])
    Intensity=np.zeros([AngleNumber,cutnumber])
    dataz=np.loadtxt(filenameqz)
    datax=np.loadtxt(filenameqx)
    for i in range (1,cutnumber+1):
        Qz[0:AngleNumber,i-1]=dataz[0:AngleNumber,i*3-2]
        Intensity[0:AngleNumber,i-1]=dataz[0:AngleNumber,i*3-1]
    for i in range(AngleNumber):
        Qx[i,0:cutnumber]=datax
    Qx=Qx*10        
    Qz=Qz*10
    return (Qz,Qx,Intensity)
    
def SCNCoordAssign(tpar,Trapnumber,X1,X2,X3,Pitch):
    Coord=np.zeros([Trapnumber+1,5,4])
    for T in range (Trapnumber+1):
        if T==0:
            Coord[T,0,0]=0
            Coord[T,1,0]=tpar[0,0]
            Coord[T,2,0]=tpar[0,1]
            Coord[T,3,0]=0
            Coord[T,4,0]=1
        else:
            Coord[T,0,0]=Coord[T-1,0,0]+0.5*(tpar[T-1,0]-tpar[T,0])
            Coord[T,1,0]=Coord[T,0,0]+tpar[T,0]
            Coord[T,2,0]=tpar[T,1]
            Coord[T,3,0]=0
            Coord[T,4,0]=1
    Coord[:,:,1]=Coord[:,:,0]
    Coord[:,0:2,1]=Coord[:,0:2,1]+X1
    Coord[:,:,2]=Coord[:,:,1]
    Coord[:,0:2,2]=Coord[:,0:2,2]+X2
    Coord[:,:,3]=Coord[:,:,2]
    Coord[:,0:2,3]=Coord[:,0:2,3]+X3
    
    for T in range (Trapnumber-2):
        if Coord[T,1,0] > Coord[T,0,1]:
            Center = (Coord[T,1,0]+Coord[T,0,1])/2
            Coord[T,1,0]=Center; Coord[T,0,1]=Center;
        if Coord[T,1,1] > Coord[T,0,2]:
            Center = (Coord[T,1,1]+Coord[T,0,2])/2
            Coord[T,1,1]=Center; Coord[T,0,2]=Center;
        if Coord[T,1,2] > Coord[T,0,3]:
            Center = (Coord[T,1,2]+Coord[T,0,3])/2
            Coord[T,1,2]=Center; Coord[T,0,3]=Center;    
        if Coord[T,1,3] > (Coord[T,0,0]+Pitch):
            Center = (Coord[T,1,3]+Coord[T,0,0]+Pitch)/2
            Coord[T,1,3]= Center; Coord[T,0,0]= Center;
    return Coord
         
         
def SCNParabolaCoord(tpar, ppar, Disc, Trapnumber, X1, X2, Pitch):
    Height = 0   
    Coord=np.zeros([Trapnumber+Disc,5,4])
    for T in range(Disc):
        # stack 1, part 1
        Hstep=tpar[0,1]/Disc
        if T ==0:
            Coord[T,0,0]=0
            Coord[T,1,0]=X1
            Coord[T,2,0]=Hstep
            Coord[T,3,0]=0
            Coord[T,4,0]=tpar[0,2]
        else:
            if Height < ppar[0,2]:
                Coord[T,0,0]=0
            else:
                Coord[T,0,0]=(-1*ppar[0,1]+np.sqrt(np.power(ppar[0,1],2)-4*ppar[0,0]*(-1*(Height-ppar[0,2]))))/(2*ppar[0,0])
            if Height < ppar[1,2]:
                Coord[T,1,0]=X1
            else:
                 Coord[T,1,0]=X1-(-1*ppar[1,1]+np.sqrt(np.power(ppar[1,1],2)-4*ppar[1,0]*(-1*(Height-ppar[1,2]))))/(2*ppar[1,0])
            Coord[T,2,0]=Hstep
            Coord[T,3,0]=0
            Coord[T,4,0]=tpar[0,2]
        # stack 2 part 1
        if T ==0:
            Coord[T,0,1]=X1
            Coord[T,1,1]=X1+X2
            Coord[T,2,1]=Hstep
            Coord[T,3,1]=0
            Coord[T,4,1]=tpar[0,2]
        else:
            if Height < ppar[1,2]:
                Coord[T,0,1]=X1
            else:
                Coord[T,0,1]=X1+(-1*ppar[1,1]+np.sqrt(np.power(ppar[1,1],2)-4*ppar[1,0]*(-1*(Height-ppar[1,2]))))/(2*ppar[1,0])
            if Height < ppar[1,2]:
                Coord[T,1,1]=X1+X2
            else:
                 Coord[T,1,1]=(X1+X2)-(-1*ppar[2,1]+np.sqrt(np.power(ppar[2,1],2)-4*ppar[2,0]*(-1*(Height-ppar[2,2]))))/(2*ppar[2,0])
            Coord[T,2,1]=Hstep
            Coord[T,3,1]=0
            Coord[T,4,1]=tpar[0,2]
            
        # stack 3 part 1
        if T ==0:
            Coord[T,0,2]=X1+X2
            Coord[T,1,2]=2*X1+X2
            Coord[T,2,2]=Hstep
            Coord[T,3,2]=0
            Coord[T,4,2]=tpar[0,2]
        else:
            if Height < ppar[2,2]:
                Coord[T,0,2]=X1+X2
            else:
                Coord[T,0,2]=(X1+X2)+(-1*ppar[2,1]+np.sqrt(np.power(ppar[2,1],2)-4*ppar[2,0]*(-1*(Height-ppar[2,2]))))/(2*ppar[2,0])
            if Height < ppar[1,2]:
                Coord[T,1,2]=2*X1+X2
            else:
                 Coord[T,1,2]=(2*X1+X2)-(-1*ppar[1,1]+np.sqrt(np.power(ppar[1,1],2)-4*ppar[1,0]*(-1*(Height-ppar[1,2]))))/(2*ppar[1,0])
            Coord[T,2,2]=Hstep
            Coord[T,3,2]=0
            Coord[T,4,2]=tpar[0,2]
            
             # stack 4 part 1
        if T ==0:
            Coord[T,0,3]=2*X1+X2
            Coord[T,1,3]=Pitch
            Coord[T,2,3]=Hstep
            Coord[T,3,3]=0
            Coord[T,4,3]=tpar[0,2]
        else:
            if Height < ppar[2,2]:
                Coord[T,0,3]=2*X1+X2
            else:
                Coord[T,0,3]=(2*X1+X2)+(-1*ppar[1,1]+np.sqrt(np.power(ppar[1,1],2)-4*ppar[1,0]*(-1*(Height-ppar[1,2]))))/(2*ppar[1,0])
            if Height < ppar[1,2]:
                Coord[T,1,3]=Pitch
            else:
                 Coord[T,1,3]=Pitch-(-1*ppar[0,1]+np.sqrt(np.power(ppar[0,1],2)-4*ppar[0,0]*(-1*(Height-ppar[0,2]))))/(2*ppar[0,0])
            Coord[T,2,3]=Hstep
            Coord[T,3,3]=0
            Coord[T,4,3]=tpar[0,2]
        Height=Height+Hstep
    for T in range(1,Trapnumber+1):
        if T==1:
            Coord[Disc,0,0]=(-1*ppar[0,1]+np.sqrt(np.power(ppar[0,1],2)-4*ppar[0,0]*(-1*(tpar[0,1]-ppar[0,2]))))/(2*ppar[0,0])
            Coord[Disc,1,0]=X1-(-1*ppar[1,1]+np.sqrt(np.power(ppar[1,1],2)-4*ppar[1,0]*(-1*(tpar[0,1]-ppar[1,2]))))/(2*ppar[1,0])
            Coord[Disc,2,0]=tpar[T,1]
            Coord[Disc,3,0]=0
            Coord[Disc,4,0]=tpar[T,2]
            
            Coord[Disc,0,1]=X1+(-1*ppar[1,1]+np.sqrt(np.power(ppar[1,1],2)-4*ppar[1,0]*(-1*(tpar[0,1]-ppar[1,2]))))/(2*ppar[1,0])
            Coord[Disc,1,1]=(X1+X2)-(-1*ppar[2,1]+np.sqrt(np.power(ppar[2,1],2)-4*ppar[2,0]*(-1*(tpar[0,1]-ppar[2,2]))))/(2*ppar[2,0])
            Coord[Disc,2,1]=tpar[T,1]
            Coord[Disc,3,1]=0
            Coord[Disc,4,1]=tpar[T,2]
            
            Coord[Disc,0,2]=(X1+X2)+(-1*ppar[2,1]+np.sqrt(np.power(ppar[2,1],2)-4*ppar[2,0]*(-1*(tpar[0,1]-ppar[2,2]))))/(2*ppar[2,0])
            Coord[Disc,1,2]=(2*X1+X2)-(-1*ppar[1,1]+np.sqrt(np.power(ppar[1,1],2)-4*ppar[1,0]*(-1*(tpar[0,1]-ppar[1,2]))))/(2*ppar[1,0])
            Coord[Disc,2,2]=tpar[T,1]
            Coord[Disc,3,2]=0
            Coord[Disc,4,2]=tpar[T,2]
            
            Coord[Disc,0,3]=(2*X1+X2)+(-1*ppar[1,1]+np.sqrt(np.power(ppar[1,1],2)-4*ppar[1,0]*(-1*(tpar[0,1]-ppar[1,2]))))/(2*ppar[1,0])
            Coord[Disc,1,3]=(Pitch)-(-1*ppar[0,1]+np.sqrt(np.power(ppar[0,1],2)-4*ppar[0,0]*(-1*(tpar[0,1]-ppar[0,2]))))/(2*ppar[0,0])
            Coord[Disc,2,3]=tpar[T,1]
            Coord[Disc,3,3]=0
            Coord[Disc,4,3]=tpar[T,2]
        else:
            Coord[Disc+T-1,0,0]=Coord[Disc+T-2,0,0]+0.5*((Coord[Disc+T-2,1,0]-Coord[Disc+T-2,0,0])-tpar[T,0])
            Coord[Disc+T-1,1,0]=Coord[Disc+T-1,0,0]+tpar[T,0]
            Coord[Disc+T-1,2,0]=tpar[T,1]
            Coord[Disc+T-1,3,0]=0;
            Coord[Disc+T-1,4,0]=tpar[T,2]
            
            Coord[Disc+T-1,0,1]=Coord[Disc+T-2,0,1]+0.5*((Coord[Disc+T-2,1,1]-Coord[Disc+T-2,0,1])-tpar[T,0])
            Coord[Disc+T-1,1,1]=Coord[Disc+T-1,0,1]+tpar[T,0]
            Coord[Disc+T-1,2,1]=tpar[T,1]
            Coord[Disc+T-1,3,1]=0;
            Coord[Disc+T-1,4,1]=tpar[T,2]
            
            Coord[Disc+T-1,0,2]=Coord[Disc+T-2,0,2]+0.5*((Coord[Disc+T-2,1,2]-Coord[Disc+T-2,0,2])-tpar[T,0])
            Coord[Disc+T-1,1,2]=Coord[Disc+T-1,0,2]+tpar[T,0]
            Coord[Disc+T-1,2,2]=tpar[T,1]
            Coord[Disc+T-1,3,2]=0;
            Coord[Disc+T-1,4,2]=tpar[T,2]
            
            Coord[Disc+T-1,0,3]=Coord[Disc+T-2,0,3]+0.5*((Coord[Disc+T-2,1,3]-Coord[Disc+T-2,0,3])-tpar[T,0])
            Coord[Disc+T-1,1,3]=Coord[Disc+T-1,0,3]+tpar[T,0]
            Coord[Disc+T-1,2,3]=tpar[T,1]
            Coord[Disc+T-1,3,3]=0;
            Coord[Disc+T-1,4,3]=tpar[T,2]
    return Coord
   
            
                        
def FreeFormTrapezoid(Coord,Qx,Qz,Trapnumber):
    H1 = Coord[0,3]
    H2 = Coord[0,3]
    form=np.zeros([len(Qx[:,1]),len(Qx[1,:])])
    for i in range(Trapnumber):
        H2 = H2+Coord[i,2];
        if i > 0:
            H1 = H1+Coord[i-1,2]
        x1 = Coord[i,0]
        x4 = Coord[i,1]
        x2 = Coord[i+1,0]
        x3 = Coord[i+1,1]
        SL = Coord[i,2]/(x2-x1)
        SR = -Coord[i,2]/(x4-x3)
        
        A1 = (np.exp(1j*Qx*((H1-SR*x4)/SR))/(Qx/SR+Qz))*(np.exp(-1j*H2*(Qx/SR+Qz))-np.exp(-1j*H1*(Qx/SR+Qz)))
        A2 = (np.exp(1j*Qx*((H1-SL*x1)/SL))/(Qx/SL+Qz))*(np.exp(-1j*H2*(Qx/SL+Qz))-np.exp(-1j*H1*(Qx/SL+Qz)))
        form=form+(1j/Qx)*(A1-A2)*Coord[i,4]
    return form
    
def SCNIntensitySim(Coord,Qx,Qz,Trapnumber,DW,I0,Bk):
    F1 = FreeFormTrapezoid(Coord[:,:,0],Qx,Qz,Trapnumber)
    F2 = FreeFormTrapezoid(Coord[:,:,1],Qx,Qz,Trapnumber)
    F3 = FreeFormTrapezoid(Coord[:,:,2],Qx,Qz,Trapnumber)
    F4 = FreeFormTrapezoid(Coord[:,:,3],Qx,Qz,Trapnumber)
    Formfactor=(F1+F2+F3+F4)
    M=np.power(np.exp(-1*(np.power(Qx,2)+np.power(Qz,2))*np.power(DW,2)),0.5);
    Formfactor=Formfactor*M
    SimInt = np.power(abs(Formfactor),2)*I0+Bk
    return SimInt
    
def PlotQzCut(Qz,SimInt,ExpI,numbercuts):
    for i in range(1,numbercuts):
        IMin=np.min(ExpI[:,i-1])
        IMax=np.max(ExpI[:,i])
        R=2*IMax/IMin
        SimInt[:,i]=SimInt[:,i]/R
        ExpI[:,i]=ExpI[:,i]/R
    for i in range(numbercuts):
        plt.semilogy(Qz[:,i],ExpI[:,i],'.')
        plt.semilogy(Qz[:,i],SimInt[:,i])
        
def plotIntelShape(tpar,ppar,Disc,Trapnumber,X,Pitch):
    Coord=CD.SCNParabolaCoord(tpar, ppar, Disc, Trapnumber, X[0,0], X[0,1], Pitch)
    Coord[:,:,5]



def Misfit(Exp,Sim):
    D= abs(Exp-Sim)
    ms=np.zeros([len(Exp[:,1]),len(Exp[1,:]),2])
    ms[:,:,0]=Sim
    ms[:,:,1]=Exp
    MS= np.nanmin(ms,2)
    Chi2=np.power((D/MS),2)
    Chi2[np.isnan(Chi2)]=0
    return Chi2
    
    
        
