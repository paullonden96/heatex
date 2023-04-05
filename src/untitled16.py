# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 14:31:08 2023

@author: ZimmermannP
"""

from src.air import air

air1 = air(10,0.85,5200)

x_g = air1.x(2)

print(x_g)