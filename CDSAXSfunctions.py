# -*- coding: utf-8 -*-
"""
Functions to be used in analyzing CDSAXS data
"""
import numpy as np
import matplotlib.pyplot as plt
import CDSAXSfunctions as CD
from scipy.interpolate import interp1d

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
    for i in range(int(Trapnumber)):
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

def ParBoundSCN(tpar,ppar,SPAR,X):
    XL=X*0.95
    XU=X*0.95
    tparL=tpar*0.9
    tparU=tpar*1.1
    pparL=ppar*0.5
    pparU=ppar*2
    SPARL=SPAR*0.5
    SPARU=SPAR*2
    return(tparL,tparU,pparL,pparU,XL,XU,SPARL,SPARU)
    """
      c=0;
    Parnumber=Trapnumber*2+1+9+3+2
    FITPAR=np.zeros([Parnumber,1])
    FITPARLB=np.zeros([Parnumber,1])
    FITPARUB=np.zeros([Parnumber,1])
    for i in range(Trapnumber-2):
        c=c+1
        FITPAR[c,0]=tpar[i,0]
        FITPARLB[c,0]=tpar[i,0]*0.9
        FITPARUB[c,0]=tpar[i,0]*1.1
        c=c+1
        FITPAR[c,0]=tpar[i,1]
        FITPARLB[c,0]=tpar[i,1]*0.9
        FITPARUB[c,0]=tpar[i,1]*1.1

    c=c+1
    FITPAR[c,0]=tpar[Trapnumber-2,0]
    FITPARLB[c,0]=tpar[Trapnumber-2,0]*0.9
    FITPARUB[c,0]=tpar[Trapnumber-2,0]*1.1
    """

def Misfit(Exp,Sim):
    D= abs(Exp-Sim)
    ms=np.zeros([len(Exp[:,1]),len(Exp[1,:]),2])
    ms[:,:,0]=Sim
    ms[:,:,1]=Exp
    MS= np.nanmin(ms,2)
    Chi2=np.power((D/MS),2)
    Chi2[np.isnan(Chi2)]=0
    return Chi2
    
def PSPVPCoord(Spline,MCoord,Trapnumber, Disc,Pitch, Offset):
    Coord=np.zeros([Disc+2,5,8])
    Coord[0,0,0]=0
    Coord[0,1,0]=Offset[0]
    
    Coord[0,0,1]= Coord[0,1,0]
    Coord[0,1,1]= Coord[0,0,1]+Offset[1]
    
    Coord[0,0,2]= Coord[0,1,1]
    Coord[0,1,2]=Coord[0,0,2]+Offset[2]
    
    Coord[0,0,3]= Coord[0,1,2]
    Coord[0,1,3]=Coord[0,0,3]+Offset[3]
    
    Coord[0,0,4]= Coord[0,1,3]
    Coord[0,1,4]=Coord[0,0,4]+Offset[4]
    
    Coord[0,0,5]= Coord[0,1,4]
    Coord[0,1,5]=Coord[0,0,5]+Offset[5]
    
    Coord[0,0,6]= Coord[0,1,5]
    Coord[0,1,6]=Coord[0,0,6]+Offset[6]
    
    Coord[0,0,7]= Coord[0,1,6]
    Coord[0,1,7]=Coord[0,0,0]+Pitch
    
    # calculates spline points for PS LIne over Template   
    SX1=Coord[0,0,0]+Spline[:,1]
    SX2=Coord[0,1,0]-Spline[:,1]
    SY2=Spline[:,0]
    xx = np.arange(0,Spline[Trapnumber-1,0],Spline[Trapnumber-1,0]/(Disc+1))
    SR1 = interp1d(SY2,SX2,kind='cubic')
    R1 = SR1(xx)
    SL1 = interp1d(SY2,SX1,kind='cubic')
    L1 = SL1(xx)
    
    for i in range(int(Disc)+1):
        Coord[i,0,0]=L1[i]
        Coord[i,1,0]=R1[i]
        Coord[i,2,0]=xx[1]
        Coord[i,4,0]=MCoord[i,0]
    # Calculates Spline poitns for PS Line 2
    SX3=Coord[0,0,2]+Spline[:,2]
    SX4=Coord[0,1,2]-Spline[:,3]

    SR2 = interp1d(SY2,SX3,kind='cubic')
    R2 = SR2(xx)
    SL2 = interp1d(SY2,SX4,kind='cubic')
    L2 = SL2(xx)
    
    for i in range(int(Disc)+1):
        Coord[i,0,2]=R2[i]
        Coord[i,1,2]=L2[i]
        Coord[i,2,2]=xx[1]
        Coord[i,4,2]=MCoord[i,2]
    
    # Calculates Spline poitns for PS Line 3
    SX5=Coord[0,0,4]+Spline[:,4]
    SX6=Coord[0,1,4]-Spline[:,4]

    SR3 = interp1d(SY2,SX5,kind='cubic')
    R3 = SR3(xx)
    SL3 = interp1d(SY2,SX6,kind='cubic')
    L3 = SL3(xx)
    
    for i in range(int(Disc)+1):
        Coord[i,0,4]=R3[i]
        Coord[i,1,4]=L3[i]
        Coord[i,2,4]=xx[1]
        Coord[i,4,4]=MCoord[i,4]
    
    # Calculates Spline poitns for PS Line 4
    SX7=Coord[0,0,6]+Spline[:,3]
    SX8=Coord[0,1,6]-Spline[:,2]

    SR4 = interp1d(SY2,SX7,kind='cubic')
    R4 = SR4(xx)
    SL4 = interp1d(SY2,SX8,kind='cubic')
    L4 = SL4(xx)
    
    for i in range(int(Disc)+1):
        Coord[i,0,6]=R4[i]
        Coord[i,1,6]=L4[i]
        Coord[i,2,6]=xx[1]
        Coord[i,4,6]=MCoord[i,6]
    
    S=1
    for T in range(int(Disc)+1):
        Coord[T,0,S]=Coord[T,1,S-1]
        Coord[T,1,S]=Coord[T,0,S+1]
        Coord[T,2,S]=xx[1]
        Coord[T,4,S]=MCoord[T,S-1]
        
    S=3
    for T in range(int(Disc)+1):
        Coord[T,0,S]=Coord[T,1,S-1]
        Coord[T,1,S]=Coord[T,0,S+1]
        Coord[T,2,S]=xx[1]
        Coord[T,4,S]=MCoord[T,S-1]
        
    S=5
    for T in range(int(Disc)+1):
        Coord[T,0,S]=Coord[T,1,S-1]
        Coord[T,1,S]=Coord[T,0,S+1]
        Coord[T,2,S]=xx[1]
        Coord[T,4,S]=MCoord[T,S-1]
        
    S=7
    for T in range(int(Disc)+1):
        Coord[T,0,S]=Coord[T,1,S-1]
        Coord[T,1,S]=Pitch+Coord[T,0,0]
        Coord[T,2,S]=xx[1]
        Coord[T,4,S]=MCoord[T,S-1]
    
    
    return (Coord)
    
def PSPVP_PB(Offset, Spline,SPAR,Trapnumber,Disc):
     
    SplineLB=Spline-1
    SplineUB=Spline+1
    
    OffsetLB=Offset-1
    OffsetUB=Offset+1
    SPARLB=SPAR[0:3]*0.8
    SPARUB=SPAR[0:3]*1.2

    FITPAR1=Spline.ravel();FITPAR2=Offset.ravel();
    FITPAR= np.append(FITPAR1,FITPAR2)
    FITPAR=np.append(FITPAR,SPAR)
    FITPAR=np.append(FITPAR,Trapnumber)
    FITPAR=np.append(FITPAR,Disc)
    
    
    FITPARLB=np.append(SplineLB.ravel(),OffsetLB.ravel())
    FITPARLB[0]=0
    FITPARLB=np.append(FITPARLB,SPARLB)

    FITPARUB=np.append(SplineUB.ravel(),OffsetUB.ravel())
    FITPARUB=np.append(FITPARUB,SPARUB)    
    
    return (FITPAR,FITPARLB,FITPARUB)
    

    """   
def MCMCInit_PSPVP(FITPAR,FITPARLB,FITPARUB,MCPAR):
    
    MCMCInit=np.zeros([10,int(MCPAR[1])+1])
    
    C=int(MCPAR[0])
    for i in range(C):
        MCMCInit[i,Trapnumber]=5
    return MCMCInit
    """

#def MCMCSeeding(tpar,ppar,SPAR,X,Qx,Qz,Intensity,Disc,Trapnumber,)
        
