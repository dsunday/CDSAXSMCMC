# -*- coding: utf-8 -*-
"""
Functions to be used in analyzing CDSAXS data
"""
import numpy as np


def CDSAXSimport(filenameqz,filenameqx,cutnumber):
    Qz=np.zeros([121,cutnumber])
    Qx=np.zeros([121,cutnumber])
    Intensity=np.zeros([121,cutnumber])
    dataz=np.loadtxt(filenameqz)
    datax=np.loadtxt(filenameqx)
    for i in range (1,cutnumber+1):
        Qz[0:121,i-1]=dataz[0:121,i*3-2]
        Intensity[0:121,i-1]=dataz[0:121,i*3-1]
    for i in range(121):
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
    
    
def Misfit(Exp,Sim):
    D= abs(Exp-Sim)
    ms=np.zeros([len(Exp[:,1]),len(Exp[1,:])])
    ms[:,:,0]=Sim
    ms[:,:,1]=Exp
    MS= np.nanmin(ms,2)
    Chi2=np.power((D/MS),2)
    return Chi2
    
    
        
