

"""

runs using 
1) python/3.3.2              2) numpy/1.11.1/python/3.3  		 3) scipy/0.17.1/python/3.3

module load python3/3.3.2 /scratch/mhess/py3p3/3.3

Aim: do MOR steady state solver with efficient matrices

"""



import numpy as np
import time
from ITHACA_SEM_code import ROM_Oseen_Iteration_step
from ITHACA_SEM_code import Generate_Advection_Terms
from ITHACA_SEM_code import ROM_Oseen_Iteration_step


def remove_from_matrix(matrix, columns, rows):
	return [
		[float(matrix[row_num][col_num])
		for col_num in range(len(matrix[row_num]))
		if not col_num in columns]

		for row_num in range(len(matrix))
		if not row_num in rows]

#(Dbnd, Dint, sing_A, sing_B, sing_Btilde, sing_C, M_trafo_no_pp_incl_dbc) = gen_Mats_cav.gen()

def Oseen_Iteration(mu_vec, Dbnd, Dint, sing_A, sing_B, sing_Btilde, sing_C, M_trafo_no_pp_incl_dbc, glodofphys, LocGloMapMatA, LocGloSign, LocGloMap, snap, ngc, nlc, bmap, imap, num_elem, nsize_bndry, nsize_int, nsize_p, nvel, c_f_bnd, c_f_p, c_f_int, NumDirBCs, nGlobHomBndDofs, no_loc_dbc, elem_loc_dbc, no_not_loc_dbc, elem_not_loc_dbc, IP, IPp, bwdtrans, cartmap0, cartmap1, LocGloBndSign, LocGloBndMap, BndCondCoeffsToGlobalCoeffsMap , init_bnd, forcing):

	print("snap.shape ", snap.shape)
	snap_x = snap[:, 0:ngc]
	snap_y = snap[:, ngc:2*ngc]
	print("snap_x.shape ", snap_x.shape)
	print("snap_y.shape ", snap_y.shape)

	snap_x_mean = np.mean(snap_x, axis=0)
	snap_y_mean = np.mean(snap_y, axis=0)

	nbmap = bmap.shape[0]
	nimap = imap.shape[0]
	nbnd = nbmap
	nint = nimap
	nsize_bndry_p1 = nsize_bndry + 1
	nsize_p_m1 = nsize_p - 1

	u_x = snap_x_mean
	u_y = snap_y_mean
	
	supr_f_bnd = np.dot( np.transpose(Dbnd) , c_f_p )
	supr_f_int = np.dot( np.transpose(Dint) , c_f_p )
	
	c_f_all = np.r_[c_f_bnd, c_f_p, c_f_int]
	(c_f_all, R1, R2) = np.linalg.svd(c_f_all, full_matrices=False)	
	c_f_all = c_f_all[:, 0:np.max(np.where(np.cumsum(R1)/np.sum(R1) < .99))+1] 

	c_f_bnd = c_f_all[0:c_f_bnd.shape[0],:]
	c_f_p = c_f_all[c_f_bnd.shape[0]:c_f_bnd.shape[0]+c_f_p.shape[0],:]
	c_f_int = c_f_all[c_f_bnd.shape[0]+c_f_p.shape[0]:c_f_bnd.shape[0]+c_f_p.shape[0]+c_f_int.shape[0],:]


	print("c_f_bnd.shape ", c_f_bnd.shape)
	print("c_f_p.shape ", c_f_p.shape)
	print("c_f_int.shape ", c_f_int.shape)


	nBndDofs = nGlobHomBndDofs + NumDirBCs
	nGlobBndDofs = nBndDofs

	no_of_dbc_in_loc = no_loc_dbc
	loc_dbc_set = elem_loc_dbc

	RBsize = c_f_bnd.shape[1]

	V_f_bnd = np.zeros([num_elem*nsize_bndry - no_of_dbc_in_loc, c_f_bnd.shape[1]])
	cropped_counter = 0
	for i in range(0, num_elem*nsize_bndry):
		if i not in loc_dbc_set:
			V_f_bnd[cropped_counter,:] = c_f_bnd[i,:]
			cropped_counter = cropped_counter + 1

	V_f_p = c_f_p
	V_f_int = c_f_int

	f_bnd_dbc = 0*np.arange((num_elem*nsize_bndry - no_not_loc_dbc)*1.0)
	non_cropped_counter = 0
	for i in range(0,num_elem*nsize_bndry):
		if i in loc_dbc_set:
			f_bnd_dbc[non_cropped_counter] = c_f_bnd[i,0]
			non_cropped_counter += 1

	t = time.time()

	sing_A_M = np.dot(np.transpose(M_trafo_no_pp_incl_dbc),sing_A)
	sing_A_MM = np.dot(M_trafo_no_pp_incl_dbc,sing_A_M)
	Dbnd_M = np.dot(Dbnd, M_trafo_no_pp_incl_dbc)
	Dbnd_MM = np.dot(Dbnd_M, np.transpose(M_trafo_no_pp_incl_dbc))
	sing_B_M = np.dot(np.transpose(M_trafo_no_pp_incl_dbc),sing_B)
	sing_B_MM = np.dot(M_trafo_no_pp_incl_dbc,sing_B_M)
	sing_A_MM_cropped = np.array(remove_from_matrix(sing_A_MM, loc_dbc_set, loc_dbc_set))
	sing_A_MM_Aud = np.array(remove_from_matrix(sing_A_MM, elem_not_loc_dbc, loc_dbc_set))
	Dbnd_crop = np.array(remove_from_matrix(Dbnd, loc_dbc_set, set()))
	sing_Dbnd_Aud = np.array(remove_from_matrix(Dbnd, elem_not_loc_dbc, set()))
	sing_B_MM_cropped = np.array(remove_from_matrix(sing_B_MM, set(), loc_dbc_set))
	Dbnd_MM_crop = np.array(remove_from_matrix(Dbnd_MM, loc_dbc_set, set()))
	sing_Btilde_cropped = np.array(remove_from_matrix(sing_Btilde, set(), loc_dbc_set))
	sing_Btilde_Aud = np.array(remove_from_matrix(sing_Btilde, set(), elem_not_loc_dbc))
	A_11_1 = np.dot(np.dot(np.transpose(V_f_bnd), sing_A_MM_cropped), V_f_bnd)
	A_21_alt = -np.dot(np.dot(np.transpose(V_f_p), Dbnd_crop), V_f_bnd)
	A_31_1 = np.dot(np.dot(np.transpose(V_f_int), np.transpose(sing_Btilde_cropped)), V_f_bnd)
	A_12_alt = -np.dot(np.dot(np.transpose(V_f_bnd), np.transpose(Dbnd_MM_crop)), V_f_p)
	A_32_alt = -np.dot(np.dot(np.transpose(V_f_int), np.transpose(Dint)), V_f_p)
	A_13_1 = np.dot(np.dot(np.transpose(V_f_bnd), sing_B_MM_cropped), V_f_int)
	A_23_alt = np.transpose(A_32_alt)
	A_33_1 = np.dot(np.dot(np.transpose(V_f_int), sing_C), V_f_int)
	r1_1 = np.dot(np.dot(np.transpose(V_f_bnd), sing_A_MM_Aud), f_bnd_dbc)
	r2_1 = -np.dot(np.dot(np.transpose(V_f_p), sing_Dbnd_Aud), f_bnd_dbc)
	r3_1 = np.dot(np.dot(np.transpose(V_f_int), np.transpose(sing_Btilde_Aud)), f_bnd_dbc)

	Dbnd_M = np.dot(Dbnd, M_trafo_no_pp_incl_dbc)
	Dbnd_MM = np.dot(Dbnd_M, np.transpose(M_trafo_no_pp_incl_dbc))
	Dbnd_crop = np.array(remove_from_matrix(Dbnd, loc_dbc_set, set()))
	sing_Dbnd_Aud = np.array(remove_from_matrix(Dbnd, elem_not_loc_dbc, set()))
	Dbnd_MM_crop = np.array(remove_from_matrix(Dbnd_MM, loc_dbc_set, set()))
	row1_2 = -np.dot(np.transpose(c_f_bnd), np.dot(np.transpose(Dbnd), c_f_p))
	row2_1 = -np.dot(np.transpose(c_f_p), np.dot(Dbnd, c_f_bnd))
	row2_3 = -np.dot(np.transpose(c_f_p), np.dot(Dint, c_f_int))
	row3_2 = -np.dot(np.transpose(c_f_int), np.dot(np.transpose(Dint), c_f_p))
	sing_A_pt1_proj = np.dot(np.dot(np.transpose(c_f_bnd), sing_A), c_f_bnd)
	sing_B_pt1_proj = np.dot(np.dot(np.transpose(c_f_bnd), sing_B), c_f_int)
	sing_Btilde_pt1_proj = np.dot(np.dot(np.transpose(c_f_int), np.transpose(sing_Btilde)), c_f_bnd)
	sing_C_pt1_proj = np.dot(np.dot(np.transpose(c_f_int), sing_C), c_f_int)
	elapsed = time.time() - t
	print('time MatMults: ', elapsed)
	t = time.time()
	sing_A_MM_cropped_proj_x = np.zeros([RBsize, RBsize, RBsize])
	sing_B_MM_cropped_proj_x = np.zeros([RBsize, RBsize, RBsize])
	sing_A_MM_Aud_proj_x = np.zeros([RBsize, RBsize])
	sing_A_MM_cropped_proj_y = np.zeros([RBsize, RBsize, RBsize])
	sing_B_MM_cropped_proj_y = np.zeros([RBsize, RBsize, RBsize])
	sing_A_MM_Aud_proj_y = np.zeros([RBsize, RBsize])
	sing_Btilde_cropped_proj_x = np.zeros([RBsize, RBsize, RBsize])
	sing_Btilde_cropped_proj_y = np.zeros([RBsize, RBsize, RBsize])
	sing_Btilde_Aud_proj_x = np.zeros([RBsize, RBsize])
	sing_Btilde_Aud_proj_y = np.zeros([RBsize, RBsize])
	sing_C_proj_x = np.zeros([RBsize, RBsize, RBsize])
	sing_C_proj_y = np.zeros([RBsize, RBsize, RBsize]) 
	sing_A_proj_x = np.zeros([RBsize, RBsize, RBsize])
	sing_A_proj_y = np.zeros([RBsize, RBsize, RBsize]) 
	sing_B_proj_x = np.zeros([RBsize, RBsize, RBsize])
	sing_B_proj_y = np.zeros([RBsize, RBsize, RBsize]) 
	sing_Btilde_proj_x = np.zeros([RBsize, RBsize, RBsize])
	sing_Btilde_proj_y = np.zeros([RBsize, RBsize, RBsize]) 
	nphys = glodofphys.shape[1]
	phys_basis_x = np.zeros([nphys,RBsize])
	phys_basis_y = np.zeros([nphys,RBsize])
	# use unit inverse tt_cd mapping
#	inv_tt_cd = 0*tt_cd[0:ngc]
	inv_tt_cd = np.arange(ngc)
	for i in range(0,ngc):
		inv_tt_cd[i] = int(i)
		#inv_tt_cd[tt_cd[i]] = int(i)
	# compute for all basis functions, the gen_Mats_adv
	nplanecoeffs = nlc
	for curr_basis in range(0, c_f_bnd.shape[1]):
		cnt = 0
		cnt1 = 0
		offset = 0
		fields = np.zeros( [nvel, nplanecoeffs] )
		for curr_elem in range(0, num_elem):
			for j in range(0, nvel):
				for k in range(0, nbnd):
					fields[j, offset + int(bmap[k])] = c_f_bnd[cnt + k, curr_basis]
				for k in range(0, nint):
					fields[j, offset + int(imap[k])] = c_f_int[cnt1+k, curr_basis]
				cnt += nbnd
				cnt1 += nint
			offset += nbnd + nint # is 169
		nvel = 2
		cnt = 0
		offset = 0
		velo0 = np.zeros( [ngc] )
		velo1 = np.zeros( [ngc] )
		for i in range(0, nplanecoeffs):
			velo0 = np.dot(LocGloMapMatA, np.transpose(fields[0,:]))
			velo1 = np.dot(LocGloMapMatA, np.transpose(fields[1,:]))
	#		velo0[int(LocGloMapA[i])] = int(LocGloSignA[i]) * fields[0,i]
	#		velo1[int(LocGloMapA[i])] = int(LocGloSignA[i]) * fields[1,i]  # use here inv_tt_cd ??
		velo0 = velo0[inv_tt_cd[0:ngc]]
		velo1 = velo1[inv_tt_cd[0:ngc]]
		basis_x = np.dot(velo0, glodofphys)
		basis_y = np.dot(velo1, glodofphys)
		phys_basis_x[:,curr_basis] = basis_x
		phys_basis_y[:,curr_basis] = basis_y

	(phys_basis_x, R1, R2) = np.linalg.svd(phys_basis_x, full_matrices=False)	
	(phys_basis_y, R1, R2) = np.linalg.svd(phys_basis_y, full_matrices=False)
	for curr_basis in range(0, c_f_bnd.shape[1]):
		print(curr_basis)
		basis_x = phys_basis_x[:,curr_basis]
		basis_y = phys_basis_y[:,curr_basis]
		ts = time.time()
		(sing_A_MM_cropped_proj_x[curr_basis,:,:], sing_A_MM_Aud_proj_x[curr_basis,:], sing_B_MM_cropped_proj_x[curr_basis,:,:], sing_Btilde_cropped_proj_x[curr_basis,:,:], sing_Btilde_Aud_proj_x[curr_basis,:], sing_C_proj_x[curr_basis,:,:], sing_A_proj_x[curr_basis,:,:], sing_B_proj_x[curr_basis,:,:], sing_Btilde_proj_x[curr_basis,:,:]) = Generate_Advection_Terms.gen(basis_x, 0*basis_y, M_trafo_no_pp_incl_dbc, V_f_bnd, V_f_int, f_bnd_dbc, elem_not_loc_dbc, loc_dbc_set, c_f_bnd, glodofphys, IP, IPp, bwdtrans, cartmap0, cartmap1, LocGloSign, LocGloMap, LocGloMapMatA, bmap, imap, LocGloBndSign, LocGloBndMap, BndCondCoeffsToGlobalCoeffsMap)
		elapsed = time.time() - ts
		print('adv gen time : ', elapsed)
		(sing_A_MM_cropped_proj_y[curr_basis,:,:], sing_A_MM_Aud_proj_y[curr_basis,:], sing_B_MM_cropped_proj_y[curr_basis,:,:], sing_Btilde_cropped_proj_y[curr_basis,:,:], sing_Btilde_Aud_proj_y[curr_basis,:], sing_C_proj_y[curr_basis,:,:], sing_A_proj_y[curr_basis,:,:], sing_B_proj_y[curr_basis,:,:], sing_Btilde_proj_y[curr_basis,:,:]) = Generate_Advection_Terms.gen(0*basis_x, basis_y, M_trafo_no_pp_incl_dbc, V_f_bnd, V_f_int, f_bnd_dbc, elem_not_loc_dbc, loc_dbc_set, c_f_bnd, glodofphys, IP, IPp, bwdtrans, cartmap0, cartmap1, LocGloSign, LocGloMap, LocGloMapMatA, bmap, imap, LocGloBndSign, LocGloBndMap, BndCondCoeffsToGlobalCoeffsMap)


	elapsed = time.time() - t
	print('time advTerms: ', elapsed)

	u_x_glo = u_x
	u_y_glo = u_y
	u_x = np.dot(u_x, glodofphys)
	u_y = np.dot(u_y, glodofphys)

	u_x_proj = np.dot(u_x, phys_basis_x)
	u_y_proj = np.dot(u_y, phys_basis_y)
	u_x_reproj = np.dot(phys_basis_x, u_x_proj)
	u_y_reproj = np.dot(phys_basis_y, u_y_proj)
	print("rel err u_x proj ", np.linalg.norm(u_x - u_x_reproj) / np.linalg.norm(u_x))
	print("rel err u_y proj ", np.linalg.norm(u_y - u_y_reproj) / np.linalg.norm(u_y))

	sing_A_MM_cropped_proj_sum = np.zeros([RBsize, RBsize])
	sing_B_MM_cropped_proj_sum = np.zeros([RBsize, RBsize])
	sing_Btilde_cropped_proj_sum = np.zeros([RBsize, RBsize])
	sing_C_proj_sum = np.zeros([RBsize, RBsize])
	sing_A_proj_sum = np.zeros([RBsize, RBsize])
	sing_B_proj_sum = np.zeros([RBsize, RBsize])
	sing_Btilde_proj_sum = np.zeros([RBsize, RBsize])
	sing_A_MM_Aud_proj_sum = 0.0*np.arange(RBsize*1.0)
	sing_Btilde_Aud_proj_sum = 0.0*np.arange(RBsize*1.0)
	for i in range(0, RBsize):
		sing_A_MM_cropped_proj_sum += u_x_proj[i] * sing_A_MM_cropped_proj_x[i,:,:] + u_y_proj[i] * sing_A_MM_cropped_proj_y[i,:,:]
		sing_B_MM_cropped_proj_sum += u_x_proj[i] * sing_B_MM_cropped_proj_x[i,:,:] + u_y_proj[i] * sing_B_MM_cropped_proj_y[i,:,:]
		sing_Btilde_cropped_proj_sum += u_x_proj[i] * sing_Btilde_cropped_proj_x[i,:,:] + u_y_proj[i] * sing_Btilde_cropped_proj_y[i,:,:]
		sing_C_proj_sum += u_x_proj[i] * sing_C_proj_x[i,:,:] + u_y_proj[i] * sing_C_proj_y[i,:,:]
		sing_A_proj_sum += u_x_proj[i] * sing_A_proj_x[i,:,:] + u_y_proj[i] * sing_A_proj_y[i,:,:]
		sing_B_proj_sum += u_x_proj[i] * sing_B_proj_x[i,:,:] + u_y_proj[i] * sing_B_proj_y[i,:,:]
		sing_Btilde_proj_sum += u_x_proj[i] * sing_Btilde_proj_x[i,:,:] + u_y_proj[i] * sing_Btilde_proj_y[i,:,:]
		sing_A_MM_Aud_proj_sum += u_x_proj[i] * sing_A_MM_Aud_proj_x[i,:] + u_y_proj[i] * sing_A_MM_Aud_proj_y[i,:]
		sing_Btilde_Aud_proj_sum += u_x_proj[i] * sing_Btilde_Aud_proj_x[i,:] + u_y_proj[i] * sing_Btilde_Aud_proj_y[i,:]


	mKinvis = mu_vec[0]
	Gr_factor = mu_vec[1] # relative to 20e3
	#Gr_factor = paracent_end/20 # relative to 20e3
	curr_iter = 1
	t = time.time()
	(u_x_new, u_y_new) = ROM_Oseen_Iteration_step.Oseen_step(u_x, u_y, mKinvis, Gr_factor, curr_iter, f_bnd_dbc, sing_C_proj_sum, c_f_bnd, c_f_p, c_f_int, V_f_bnd, V_f_p, V_f_int, sing_A, sing_B, sing_Btilde, sing_C, sing_A_proj_sum, sing_B_proj_sum, 	sing_Btilde_proj_sum, row1_2, row2_1, row2_3, row3_2, sing_A_pt1_proj, sing_B_pt1_proj, sing_Btilde_pt1_proj, sing_C_pt1_proj, glodofphys, IP, IPp, bwdtrans, cartmap0, cartmap1, LocGloSign, LocGloMap, LocGloMapMatA, bmap, imap, LocGloBndSign, LocGloBndMap, BndCondCoeffsToGlobalCoeffsMap, init_bnd, forcing, nlc, A_11_1, A_21_alt, A_31_1, A_12_alt, A_32_alt, A_13_1, A_23_alt, A_33_1, r1_1, r2_1, r3_1, sing_A_MM_cropped_proj_sum, sing_A_MM_Aud_proj_sum, sing_B_MM_cropped_proj_sum, sing_Btilde_cropped_proj_sum, sing_Btilde_Aud_proj_sum, loc_dbc_set)
	elapsed = time.time() - t
	print('time : ', elapsed)



	#print(inv_tt_cd)

#	print("np.linalg.norm(u_x_glo): ", np.linalg.norm(u_x_glo))
#	print("np.linalg.norm(u_x_new): ", np.linalg.norm(u_x_new))
#	print("np.linalg.norm(u_x_new - u_x_glo): ", np.linalg.norm(u_x_new - u_x_glo[tt_cd[0:ngc]]))
#	print("np.linalg.norm(u_x_new - u_x_glo): ", np.linalg.norm(u_x_new[inv_tt_cd[0:ngc]] - u_x_glo))
#	print("np.linalg.norm(u_y_glo): ", np.linalg.norm(u_y_glo))
#	print("np.linalg.norm(u_y_new): ", np.linalg.norm(u_y_new))
#	print("np.linalg.norm(u_y_new - u_y_glo): ", np.linalg.norm(u_y_new - u_y_glo[tt_cd[0:ngc]]))
#	print("np.linalg.norm(u_y_new - u_y_glo): ", np.linalg.norm(u_y_new[inv_tt_cd[0:ngc]] - u_y_glo))

	#print(u_x[tt_cd[300:400]])
	#print(u_x_new[300:400])

	u_x_new = u_x_new[inv_tt_cd[0:ngc]]
	u_y_new = u_y_new[inv_tt_cd[0:ngc]]


	rel_err = 1	
	u_x_new = snap_x_mean
	u_y_new = snap_y_mean
	while rel_err > 1e-6:
		u_x_glo = u_x_new
		u_y_glo = u_y_new
		u_x = np.dot(u_x_glo, glodofphys)
		u_y = np.dot(u_y_glo, glodofphys)
		u_x_proj = np.dot(u_x, phys_basis_x)
		u_y_proj = np.dot(u_y, phys_basis_y)
		u_x_reproj = np.dot(phys_basis_x, u_x_proj)
		u_y_reproj = np.dot(phys_basis_y, u_y_proj)
		sing_A_MM_cropped_proj_sum = np.zeros([RBsize, RBsize])
		sing_B_MM_cropped_proj_sum = np.zeros([RBsize, RBsize])
		sing_Btilde_cropped_proj_sum = np.zeros([RBsize, RBsize])
		sing_C_proj_sum = np.zeros([RBsize, RBsize])
		sing_A_proj_sum = np.zeros([RBsize, RBsize])
		sing_B_proj_sum = np.zeros([RBsize, RBsize])
		sing_Btilde_proj_sum = np.zeros([RBsize, RBsize])
		sing_A_MM_Aud_proj_sum = 0.0*np.arange(RBsize*1.0)
		sing_Btilde_Aud_proj_sum = 0.0*np.arange(RBsize*1.0)
		for i in range(0, RBsize):
			sing_A_MM_cropped_proj_sum += u_x_proj[i] * sing_A_MM_cropped_proj_x[i,:,:] + u_y_proj[i] * sing_A_MM_cropped_proj_y[i,:,:]
			sing_B_MM_cropped_proj_sum += u_x_proj[i] * sing_B_MM_cropped_proj_x[i,:,:] + u_y_proj[i] * sing_B_MM_cropped_proj_y[i,:,:]
			sing_Btilde_cropped_proj_sum += u_x_proj[i] * sing_Btilde_cropped_proj_x[i,:,:] + u_y_proj[i] * sing_Btilde_cropped_proj_y[i,:,:]
			sing_C_proj_sum += u_x_proj[i] * sing_C_proj_x[i,:,:] + u_y_proj[i] * sing_C_proj_y[i,:,:]
			sing_A_proj_sum += u_x_proj[i] * sing_A_proj_x[i,:,:] + u_y_proj[i] * sing_A_proj_y[i,:,:]
			sing_B_proj_sum += u_x_proj[i] * sing_B_proj_x[i,:,:] + u_y_proj[i] * sing_B_proj_y[i,:,:]
			sing_Btilde_proj_sum += u_x_proj[i] * sing_Btilde_proj_x[i,:,:] + u_y_proj[i] * sing_Btilde_proj_y[i,:,:]
			sing_A_MM_Aud_proj_sum += u_x_proj[i] * sing_A_MM_Aud_proj_x[i,:] + u_y_proj[i] * sing_A_MM_Aud_proj_y[i,:]
			sing_Btilde_Aud_proj_sum += u_x_proj[i] * sing_Btilde_Aud_proj_x[i,:] + u_y_proj[i] * sing_Btilde_Aud_proj_y[i,:]
		(u_x_new, u_y_new) = ROM_Oseen_Iteration_step.Oseen_step(u_x, u_y, mKinvis, Gr_factor, curr_iter, f_bnd_dbc, sing_C_proj_sum, c_f_bnd, c_f_p, c_f_int, V_f_bnd, V_f_p, V_f_int, sing_A, sing_B, sing_Btilde, sing_C, sing_A_proj_sum, sing_B_proj_sum, sing_Btilde_proj_sum, row1_2, row2_1, row2_3, row3_2, sing_A_pt1_proj, sing_B_pt1_proj, sing_Btilde_pt1_proj, sing_C_pt1_proj, glodofphys, IP, IPp, bwdtrans, cartmap0, cartmap1, LocGloSign, LocGloMap, LocGloMapMatA, bmap, imap, LocGloBndSign, LocGloBndMap, BndCondCoeffsToGlobalCoeffsMap, init_bnd, forcing, nlc, A_11_1, A_21_alt, A_31_1, A_12_alt, A_32_alt, A_13_1, A_23_alt, A_33_1, r1_1, r2_1, r3_1, sing_A_MM_cropped_proj_sum, sing_A_MM_Aud_proj_sum, sing_B_MM_cropped_proj_sum, sing_Btilde_cropped_proj_sum, sing_Btilde_Aud_proj_sum, loc_dbc_set)
		u_x_new = u_x_new[inv_tt_cd[0:ngc]]
		u_y_new = u_y_new[inv_tt_cd[0:ngc]]
		rel_err = np.linalg.norm(u_x_new - u_x_glo) / np.linalg.norm(u_x_glo) + np.linalg.norm(u_y_new - u_y_glo) / np.linalg.norm(u_y_glo)
		Sxp_mu_1R = np.dot(u_x_new[0:ngc], glodofphys)
		print(rel_err)
	print(Gr_factor * 20, "  ", Sxp_mu_1R[1296])
	
	return (u_x_new, u_y_new)

