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
outfile = open('force.dat', 'w')
data = infile.read()
data = data.replace("(", "")
data = data.replace(")", "")
outfile.write(data)
outfile.close()
infile.close()


infile = open(str(sys.argv[2]), 'r')
outfile = open('avg_vel.dat', 'w')
data = infile.read()
data = data.replace("(", "")
data = data.replace(")", "")
outfile.write(data)
outfile.close()
infile.close()



forces = np.loadtxt('force.dat', skiprows=3)

t = forces[:,0]

Fp_x = forces[:,4]
Fp_y = forces[:,5]
Fp_z = forces[:,6]

Fv_x = forces[:,7]
Fv_y = forces[:,8]
Fv_z = forces[:,9]


F_x = Fp_x + Fv_x
F_y = Fp_y + Fv_y
F_z = Fp_z + Fv_z



velocities = np.loadtxt('avg_vel.dat', skiprows=3)

vol = float(linecache.getline('avg_vel.dat', 4)[15:])


u_x = velocities[:,1]
u_y = velocities[:,2]
u_z = velocities[:,3]

##############################


nu = 1e-5 # cinematic viscosity

d = 2*float(linecache.getline('PyFoamPrepareCaseParameters', 32)[2:-2])*0.001
h = 0.002
area = np.pi*d*h

##### create the windows for the fft and relatives parameters of the signal


ts = t[-1] - t[-2] # sampling interval ( the delta t of the simulation)

idx = (t>=float(sys.argv[3]))  # cut the part of the signal that is not already periodic

t = t[idx]

F_x = F_x[idx]
F_y = F_y[idx]
F_z = F_z[idx]

Fp_x = Fp_x[idx]
Fp_y = Fp_y[idx]
Fp_z = Fp_z[idx]

Fv_x = Fv_x[idx]
Fv_y = Fv_y[idx]
Fv_z = Fv_z[idx]



u_x = (u_x[idx])/vol
u_y = (u_y[idx])/vol
u_z = (u_z[idx])/vol

U = int_mean(t,u_x)
V = int_mean(t,u_y)
W = int_mean(t,u_z)

u = np.sqrt(U**2 +V**2 +W**2)  # reference velocity


Cd = F_x/(0.5*area*u**2)
Cl = F_y/(0.5*area*u**2)
Cz = F_z/(0.5*area*u**2)


#############################


#pdb.set_trace()


n = len(t)  # the numbers of samplings is needed for the normalization
k = np.arange(n)

## the FFT analysisis done with the scipy package


P = np.zeros([n/2 ,3])
F = np.zeros([n/2 ,3])
ST = np.zeros([1,3])

i=0
for sig in  [Cd,Cl,Cz]: 
	sample_frq = fftpack.fftfreq(sig.size, d=ts)
	sig_fft = fftpack.fft(sig)
	pidxs = np.where(sample_frq > 0)
	freqs = sample_frq[pidxs]  # frequencies of the signal
	power = np.abs(sig_fft)[pidxs]/ (n/2)  # normalized amplitude of the signal
	main_sig = ifft(sig_fft)
	## computing the Strouhal number 

	f = freqs[power.argmax()]  # frenquency of the first armonic

	St = (f*d)/u  # Strouhal number
	F[:,i] = freqs
	P[:,i] = power
	ST[0,i] = St
	i = i+1



fig, ([ax1, ax2], [ax3, ax4], [ax5, ax6], [ax7, ax8]) = plt.subplots(4, 2,figsize=(16, 8), dpi=150)

#fig.gca().ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
ax1.plot(t,Cd)
#ax[0].set_xticks((0.92e-5, 1.02e-5))
ax1.set_xlabel('t [s]',fontsize=10)
ax1.set_ylabel(r'$C_d$',fontsize=20)
ax1.grid()
#ax[0].set_ylim((-1,1))
#ax[0].set_xlim((500,700))

ax3.plot(F[:,0], P[:,0], 'r') # plotting the spectrum
ax3.set_xlabel(r'$freq$',fontsize=20)
ax3.set_ylabel(r'$P_{C_d}$',fontsize=20)
ax3.set_xlim((0,2000))
ax3.grid()

ax2.plot(t,Cl)
ax2.set_xlabel('t [s]',fontsize=10)
ax2.set_ylabel(r'$C_l$',fontsize=20)
ax2.grid()

ax4.plot(F[:,1], P[:,1], 'r') # plotting the spectrum
ax4.set_xlabel(r'$freq$',fontsize=20)
ax4.set_ylabel(r'$P_{C_l}$',fontsize=20)
ax4.set_xlim((0,2000))
ax4.grid()

ax5.plot(t,Cz)
ax5.set_xlabel('t [s]',fontsize=10)
ax5.set_ylabel(r'$C_z$',fontsize=20)
ax5.grid()

ax7.plot(F[:,2], P[:,2], 'r') # plotting the spectrum
ax7.set_xlabel(r'$freq$',fontsize=20)
ax7.set_ylabel(r'$P_{C_z}$',fontsize=20)
ax7.set_xlim((0,2000))
ax7.grid()

ax6.plot(t,u_x,'y')
ax6.set_xlabel('t [s]',fontsize=10)
ax6.set_ylabel(r'$<U>$',fontsize=20)
ax6.grid()

ax8.plot(t,u_y,'g', t,u_z,'m')
ax8.set_xlabel('t [s]',fontsize=10)
ax8.set_ylabel(r'$<V> , <W>$',fontsize=20)
ax8.grid()
plt.show()

fig.savefig('force_avg.png', bbox_inches='tight',dpi=500)


###### Compute statistic of the signal #####
Re = (u*d)/nu
theta = float(linecache.getline('PyFoamPrepareCaseParameters', 30)[6:-2])
phi = float(linecache.getline('PyFoamPrepareCaseParameters', 6)[4:-2])
g = float(linecache.getline('PyFoamPrepareCaseParameters', 12)[2:-2])
porosity = float(linecache.getline('PyFoamPrepareCaseParameters', 24)[9:-2])

print theta, phi, g, porosity, int_mean(t,Fp_x), int_mean(t,Fp_y), int_mean(t,Fp_z), int_mean(t,Fv_x), int_mean(t,Fv_y), int_mean(t,Fv_z), U, V, W, Re


"""
print "# points:", len(sig)
print "mean:", np.mean(sig)
print "median:", np.median(sig)
print "std:", np.std(sig)
print "variance:", np.var(sig)
print "3th quartile:",  np.percentile(sig, 75)
print "max value:", np.amax(sig), "@", t[np.argmax(sig)]
print "min value:", np.amin(sig), "@", t[np.argmin(sig)]
"""

