

"""

runs using 
1) python/3.3.2              2) numpy/1.11.1/python/3.3  		 3) scipy/0.17.1/python/3.3

module load python3/3.3.2 /scratch/mhess/py3p3/3.3

Aim: do MOR steady state solver with efficient matrices

Task1: recover a snapshot, provided as initial solution with a 2-dim ROM

"""



import numpy as np
import time
import Oseen_step_trafo_p
import gen_Mats
import gen_Mats_adv

ngc = 4753
nlc = 5408
nLocBndDofs = 3104
nGlobBndDofs = 1794
num_elem = 32
nvel = 2
nbnd = 48
nint = 121
glodofphys = np.load('glodofphys.npy', 'r')
physglodof = np.load('physglodof.npy', 'r')
tt_ref_p0075 = np.load('tt_ref_p0075.npy', 'r')
tt_ref_p005 = np.load('tt_ref_p005.npy', 'r')
LocGloBndSign = np.load('LocGloBndSign.npy', 'r')
LocGloBndMap = np.load('LocGloBndMap.npy', 'r')
LocGloSignA = np.load('LocGloSignA.npy', 'r')
LocGloMapA = np.load('LocGloMapA.npy', 'r')
LocGloSign = np.load('LocGloSign.npy', 'r')
LocGloMap = np.load('LocGloMap.npy', 'r')
bmap = (np.load('bmap.npy', 'r'))
imap = (np.load('imap.npy', 'r'))
IP = np.load('IP.npy', 'r')
IPp = np.load('IPp.npy', 'r')
bwdtrans = np.load('bwdtrans.npy', 'r')
cartmap0 = np.load('cartmap0.npy', 'r')
cartmap1 = np.load('cartmap1.npy', 'r')
c_f_bnd = np.load("collect_f_bnd.npy", "r")
c_f_p = np.load("collect_f_p.npy", "r")
c_f_int = np.load("collect_f_int.npy", "r")

f_snap = open('c1/A1.txt', 'r')
snap = np.transpose(np.loadtxt(f_snap))
print("snap.shape ", snap.shape)
print("glodofphys.shape ", glodofphys.shape)

snap_x = snap[:, 0:ngc]
snap_y = snap[:, ngc:2*ngc]

nn = np.dot(snap_y, glodofphys)
print(nn[:,378])

param_list = [  0.0025, 0.00255, 0.002625, 0.0027, 0.00275, 0.0028, 0.002875, 0.00295, 0.0030, 0.003125, 0.00325, 0.003375, 0.0035, 0.003625, 0.00375, 0.0040, 0.00425, 0.0045, 0.00475, 0.0050, 0.005125, 0.00525, 0.005375, 0.0055, 0.005625, 0.00575, 0.005875, 0.0060, 0.006125, 0.00625, 0.006375, 0.0065, 0.006625, 0.00675, 0.006875, 0.0070, 0.007125, 0.00725, 0.007375, 0.0075, 0.007625, 0.00775,  0.0080, 0.0085, 0.0090, 0.0095, 0.0100, 0.0110, 0.0120, 0.01375, 0.0175,  0.02125,  0.0250,  0.03125,  0.0375,  0.04375,  0.0500]


print("param_list ", param_list)


for i in range(0, snap.shape[0]):
	curr_snap_x = snap_x[i,:]
	curr_snap_y = snap_y[i,:]
	curr_param = param_list[i]
	curr_snap_x_loc = np.dot(curr_snap_x, glodofphys)
	curr_snap_y_loc = np.dot(curr_snap_y, glodofphys)
	print("init qoi ", curr_snap_y_loc[378])
	print("at param ", curr_param)
	(u_x_new, u_y_new, f_bnd, f_p, f_int) = Oseen_step_trafo_p.Oseen_step_rom(curr_snap_x_loc, curr_snap_y_loc, curr_param)
	if i == 0:
		collect_f_bnd = f_bnd
		collect_f_p = f_p
		collect_f_int = f_int
	else:
		collect_f_bnd = np.c_[collect_f_bnd, f_bnd]
		collect_f_p = np.c_[collect_f_p, f_p]
		collect_f_int = np.c_[collect_f_int, f_int]

np.save('trafo_f_bnd', collect_f_bnd)
np.save('trafo_f_p', collect_f_p)
np.save('trafo_f_int', collect_f_int)
