from __future__ import division
import numpy as np
from scipy.fftpack import fft, ifft
import matplotlib.pyplot as plt
import sys
import linecache
from scipy import fftpack
import pdb 
from  scipy.integrate import simps

#######################

def int_mean(time, signal):
	sig_m = simps(signal, time)/(time[-1] -time[0])
	return sig_m


########################


infile = open(str(sys.argv[1]), 'r')
outfile = open('avg_vel.dat', 'w')
data = infile.read()
data = data.replace("(", "")
data = data.replace(")", "")
outfile.write(data)
outfile.close()
infile.close()


permeab = np.loadtxt('avg_vel.dat', skiprows=3)

t = permeab[:,0]

perm = float(linecache.getline('avg_vel.dat', 4)[15:])


k_x = permeab[:,1]*perm 
k_y = permeab[:,2]*perm 
k_z = permeab[:,3]*perm 

##############################


nu = 1e-5 # cinematic viscosity





##### create the windows for the fft and relatives parameters of the signal


ts = t[-1] - t[-2] # sampling interval ( the delta t of the simulation)

idx = (t>=float(sys.argv[2]))  # cut the part of the signal that is not already periodic

t = t[idx]



k_x = (k_x[idx])
k_y = (k_y[idx])
k_z = (k_z[idx])

U = int_mean(t,u_x)
V = int_mean(t,u_y)
W = int_mean(t,u_z)




#############################



print U, V, W



