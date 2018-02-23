import numpy as np

fx = g*np.cos(theta*np.pi/180)*np.cos(phi*np.pi/180)
fy = g*np.sin(theta*np.pi/180)*np.cos(phi*np.pi/180)
fz = g*np.sin(phi*np.pi/180)

R = np.sqrt((4*(1-porosity))/(2*np.pi))
