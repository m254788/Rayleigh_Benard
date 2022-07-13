import numpy as np
import sys
sys.path.insert(1,".")
from hdf5Helper import *

params = np.loadtxt("params.txt",dtype=np.float64)
lx = int(params[0])
ly = int(params[1])
total_vis = int(params[2])
delta_x = params[3]


nnodes = lx*ly

dims = (1,int(ly),int(lx))
uz = np.zeros(nnodes)

for Vis_ind in range(total_vis):
	data = np.loadtxt("evolution"+str(Vis_ind)+".txt",dtype=np.float64)
	T = data[0:nnodes]
	ux = data[nnodes:2*nnodes]
	uy = data[2*nnodes:3*nnodes]

	h5_file = "out"+str(Vis_ind)+".h5"
	xmf_file = "data"+str(Vis_ind)+".xmf"
	writeH5(T,ux,uy,uz,h5_file)
	writeXdmf(dims,delta_x,xmf_file,h5_file)
