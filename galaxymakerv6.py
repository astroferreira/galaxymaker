##############################################################################
## This is a main script that runs routines to create synthetic galaxies	##
## 						Created by Geferson Lucatelli 						##
##						lucatelli_asph_r2@virgophysics						##
##							Physics Laboratory								##
##								03.16.2016									##
##############################################################################
#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
##IMPORT SCRIPTS
##############################################################################
import galaxymakerlibv6
from galaxymakerlibv6 import*
import os
M=5.0
N=5.0
q=0.75
c=0.0
Io=8.0
ro=0.5
n=1.5
nb=0.7
p=1.0
k=2.0
AA=1.0
NN=5.0
BB=1.60
PHI=0.6
m=1.5
gal=spiral_galaxy_model(M,N,p,k,c)
# log_spiral(m,M,N,p,k,c,Io,n,ro,q)