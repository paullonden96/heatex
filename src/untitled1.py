# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 13:21:24 2023

@author: ZimmermannP
"""
import numpy as np
f = [1]

for i in f:
    x = [i*4,i*2]
    L = [i*0.1,i*24]
    
    m1 = np.reshape(x,(len(x),1))
    m2 = np.reshape(L,(len(L),1))
print(m1.shape(m1))
print(m2.shape(m2))