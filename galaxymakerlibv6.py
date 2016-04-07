##############################################################################
## This is a library of many functions to the code galaxy_maker		 		##
## 						Created by Geferson Lucatelli 						##
##						lucatelli_asph_r2@virgophysics						##
##							Physics Laboratory								##
##								03.16.2016									##
##############################################################################

##############################################################################
##  IMPORTING SOME USEFULL MODULES
##############################################################################
from __future__ import division
import numpy as np
from numpy import*
import matplotlib.pyplot as pyplot
import scipy as sp
from scipy import*
import pylab as pl
from pylab import*
import math 
from math import*
import pyfits as pf
from pyfits import*
import os
import requests
import sys

##############################################################################
##  SOME SURFACE BRIGHTNESS PROFILES
##############################################################################
'''
Some surface brightness to model disk galaxies, elliptical galaxies, spiral galaxies.
'''
def exponential_profile(Io,ro,r):
	Iexp=0.1*Io*np.exp(-r/ro)
	return Iexp
def sersic_profile(In,rn,r,n):
	bn=2*n-1./3.
	Iser=In*np.exp(bn*(1-(r/rn)**(1/n)))
	# Io=In*np.exp(bn)
	return Iser
def bar_profile(Inb,rnb,r,nb):
	bn=2*nb-1./3.
	Ibar=Inb*np.exp(bn*(1-(r/rnb)**(1/nb)))
	return Ibar
def grid(q,c,gal_center,N,M):
	# x,y=meshgrid(arange(-M/2,M/2,0.005), arange(-N/2,N/2,0.005))
	x,y=meshgrid(arange(-M/2,M/2,M/255.0), arange(-N/2,N/2,N/255.0))
	r=(abs(x-gal_center[1])**(c+2.0)+((abs(y-gal_center[0]))/(q))**(c+2.0))**(1.0/(c+2.0))+0.01
	return r,x,y
def fourier_modes(a1,a2,a3,a4,a5,a6,ph,q,Ioe,roe):
	r,x,y=grid(q,c)
	O=np.arctan((y)/(x)*q)
	rfm=r*(1+a1*cos(1*(O+ph))+a2*cos(2*(O+ph))+a3*cos(3*(O+ph))+a4*cos(4*(O+ph))+a5*cos(5*(O+ph))+a6*cos(6*(O+ph))) #CREATE FOURIER MODES
	return rfm
'''
The reference of this solution, comes from my question on StackOverflow
http://stackoverflow.com/questions/36095775/creating-a-spiral-structure-in-python-using-hyperbolic-tangent
'''
def spiral_galaxy_model(M,N,p,k,c):
	nro__galaxias=5
	rn=np.random.uniform(80.0,120.0,nro__galaxias)
	q=np.random.uniform(0.6,1.0,nro__galaxias)
	n=np.random.uniform(1.5,3.0,nro__galaxias)
	In=np.random.uniform(100.0,250.0,nro__galaxias)
	Inb=np.random.uniform(20.0,50.0,nro__galaxias)
	rnb=np.random.uniform(10.0,25.0,nro__galaxias)
	nb=np.random.uniform(0.5,1.3,nro__galaxias)
	fig=plt.figure()
	AA=[1.0]
	BB=np.random.uniform(0.01,2.0,nro__galaxias)
	NN=np.random.uniform(1.0,3.0,nro__galaxias)
	# phi_0=2*NN*arctanh(np.exp(-AA[0]))-np.exp(-AA[0])
	phi_0=np.random.uniform(0.5,2.0,nro__galaxias)
	NNN=[]
	BBB=[]
	PHI0=[]
	ARCTANH=[]
	for i in range(len(BB)):
		power=2.0
		gal_center=(0,0)
		r,x,y=grid(q[i],c,gal_center,N,M)
		b=1./(tanh(phi_0[i]/(2.*NN[i])))
		r_break=2.0
		r_soft=1.9
		# *(	(r_break-r)/r_break	)
		phi_r_tan=2.*NN[i]*arctan(np.exp(AA[0]/(r**1.0+0.001) )/b)
		bb=2.65-4.98*(r_break/(r_break-r_soft))
		phi_r_tan_t=0.5*(np.tanh((2-bb)*r/r_break+bb)+1)
		corte=where(r>=r_break)
		I_ser=sersic_profile(In[i],rn[i],r,n[i])
		I_ser_nuc=sersic_profile(25,0.3,r,8.0)
		f_r_tan=1.*(1.0/p)*np.cos(k*(arctan2(x,y)+phi_r_tan))
		I_bar=bar_profile(Inb[i],rnb[i],r,nb[i])
		tan_spiral=((1.0*I_ser+1.0*I_bar+0.0*I_ser_nuc)+In[i]*10*f_r_tan)
		ax1=fig.add_subplot(1,3,1)
		ax1.imshow(((tan_spiral)))
		ax2=fig.add_subplot(1,3,2)
		axis('off')
		ax2.text(0.0, 0.80, r'Galaxy ID:$'+str(i+1)+'$', fontsize=12)
		ax2.text(0.0, 0.75, r'$q$='+'$'+str(format(q[i],'.2f'))+'$', fontsize=12)
		ax2.text(0.0, 0.70, r'$r_0$='+'$'+str(format(rn[i],'.2f'))+'$', fontsize=12)
		ax2.text(0.0, 0.65, r'$I_0$='+'$'+str(format(In[i],'.2f'))+'$', fontsize=12)
		ax2.text(0.0, 0.60, r'$n$='+'$'+str(format(n[i],'.2f'))+'$', fontsize=12)
		ax2.text(0.0, 0.55, r'$r_{0,bar}$='+'$'+str(format(rnb[i],'.2f'))+'$', fontsize=12)
		ax2.text(0.0, 0.50, r'$I_{0,bar}$='+'$'+str(format(Inb[i],'.2f'))+'$', fontsize=12)
		ax2.text(0.0, 0.45, r'$n_{bar}$='+'$'+str(format(nb[i],'.2f'))+'$', fontsize=12)
		ax2.text(0.0, 0.40, r'$N$='+'$'+str(format(NN[i],'.2f'))+'$', fontsize=12)
		ax2.text(0.0, 0.35, r'$\phi_0$='+'$'+str(format(phi_0[i],'.2f'))+'$', fontsize=12)
		ax2.text(0.0, 0.30, r'$B(\phi_0)$='+'$'+str(format(b,'.2f'))+'$', fontsize=12)
		ax2.text(0.0, 0.25, r'$B$='+'$'+str(format(BB[i],'.2f'))+'$', fontsize=12)
		ax3=fig.add_subplot(1,3,3,axisbg='white')
		ax3.plot(tan_spiral[len(tan_spiral)/2,len(tan_spiral)/2:])
		plt.savefig('gall_'+str(i+1)+'.jpg',dpi=200)
		pf.writeto('gall_'+str(i+1)+'.fits', tan_spiral,  clobber=1)
		print 'Galaxy ', i+1,' Completed   -----    BB=', BB[i], ' NN=',NN[i],'  phi_0=',phi_0[i],'  B_phi=', b
		NNN.append(NN[i])
		BBB.append(b)
		PHI0.append(phi_0[i])
		plt.clf()
	ax1=fig.add_subplot(131)
	plt.title('distribuicao $N$')
	ax1.plot(NNN,'o')
	ax2=fig.add_subplot(132)
	plt.title('distribuicao $B(\phi)$')
	ax2.plot(BBB,'o')
	ax3=fig.add_subplot(133)
	plt.title('distribuicao $\phi_0$')
	ax3.plot(PHI0,'o')
	plt.savefig('distribuicao.jpg',dpi=200)
	plt.clf()
	return tan_spiral
def log_spiral(m,M,N,p,k,c,Io,n,ro,q):
	gal_center=(0,0)
	r,x,y=grid(q,c,gal_center,N,M)
	ph_r=((180/np.pi)*(1/tan(m)))*np.log(r/ro)
	I_ser=sersic_profile(Io,ro,r,n)
	I_bar=bar_profile(8.0,0.4,r,0.8)
	spiral=np.cos(2*arctan2(x,y)+2*ph_r)
	gal=1*I_ser+0.5*spiral+I_bar*0
	plt.imshow(gal)
	plt.show()