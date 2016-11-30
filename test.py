# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 12:41:28 2015

@author: dfs1
"""

import numpy as np
import CDSAXSfunctions as CD
import matplotlib.pyplot as plt
Exp= Intensity
Sim=SimInt
ms=np.zeros([121,26,2])
ms[:,:,0]=Sim
ms[:,:,1]=Exp
MS= np.nanmin(ms,2)