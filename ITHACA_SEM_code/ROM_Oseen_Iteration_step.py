
"""

runs using 
1) python/3.3.2              2) numpy/1.11.1/python/3.3  		 3) scipy/0.17.1/python/3.3

module load python3/3.3.2 /scratch/mhess/py3p3/3.3

Aim: do MOR steady state solver with efficient matrices

"""


import numpy as np
import time

def remove_from_matrix(matrix, columns, rows):
	return [
		[float(matrix[row_num][col_num])
		for col_num in range(len(matrix[row_num]))
		if not col_num in columns]

		for row_num in range(len(matrix))
		if not row_num in rows]

def Oseen_step(u_x, u_y, mKinvis, Gr_factor, curr_iter, f_bnd_dbc, sing_C_proj_sum, c_f_bnd, c_f_p, c_f_int, V_f_bnd, V_f_p, V_f_int, sing_A_pt1, sing_B_pt1, sing_Btilde_pt1, sing_C_pt1, sing_A_proj_sum, sing_B_proj_sum, sing_Btilde_proj_sum, row1_2, row2_1, row2_3, row3_2, sing_A_pt1_proj, sing_B_pt1_proj, sing_Btilde_pt1_proj, sing_C_pt1_proj, glodofphys, IP, IPp, bwdtrans, cartmap0, cartmap1, LocGloSign, LocGloMap, LocGloMapMatA, bmap, imap, LocGloBndSign, LocGloBndMap, BndCondCoeffsToGlobalCoeffsMap, init_bnd, forcing, nlc, A_11_1, A_21_alt, A_31_1, A_12_alt, A_32_alt, A_13_1, A_23_alt, A_33_1, r1_1, r2_1, r3_1, sing_A_MM_cropped_proj_sum, sing_A_MM_Aud_proj_sum, sing_B_MM_cropped_proj_sum, sing_Btilde_cropped_proj_sum, sing_Btilde_Aud_proj_sum, loc_dbc_set):


	np.set_printoptions(threshold=np.inf)	

	#mKinvis = 1

#	glodofphys = np.load('glodofphys.npy', 'r')
#	LocGloSign = np.load('LocGloSign.npy', 'r')
#	LocGloMap = np.load('LocGloMap.npy', 'r')
#	LocGloSignA = np.load('LocGloSignA.npy', 'r')
#	LocGloMapA = np.load('LocGloMapA.npy', 'r')
#	bmap = (np.load('bmap.npy', 'r'))
#	imap = (np.load('imap.npy', 'r'))
#	cartmap0 = np.load('cartmap0.npy', 'r')

	nphys = cartmap0.shape[1]
	npc = IP.shape[0]
	num_elem = npc // nphys	
#	print('num_elem ', num_elem)
	nsize_bndry = (LocGloBndMap.shape[0] - num_elem) // num_elem
	nsize_bndry_p1 = nsize_bndry + 1
	nsize_int = ((LocGloMap.shape[0] - LocGloBndMap.shape[0] + num_elem) // num_elem) * 2
	nsize_p = nsize_int // 2
	nsize_p_m1 = nsize_p - 1
	ngc = glodofphys.shape[0]
	nlc = IP.shape[1] * num_elem # fields[m_velocity[0]]->GetNcoeffs() 
	nvel = 2
	nbmap = bmap.shape[0]
	nimap = imap.shape[0]
	nbnd = nbmap
	nint = nimap
	HelmMatScale = 1
	Ahrows = nsize_bndry_p1

	mK_fac = 1
	HelmMat = np.load('ITHACA_SEM_data/HM_0.npy', 'r')
	HelmMatRows = np.asarray(np.rint(np.sqrt(HelmMat.shape[0])), dtype=int)
	ncoeffs = HelmMatRows
	HelmMats = np.zeros(( num_elem , HelmMatRows*HelmMatRows ))

	for i in range(0,num_elem):
		HM_str = 'ITHACA_SEM_data/HM_' + str(i) + '.npy'
		HelmMats[i, :] = mK_fac*np.load(HM_str, 'r')

#	A_ROM_Mat_oo = sing_A_pt1_proj + sing_A_proj_sum + row1_2 + sing_B_pt1_proj + sing_B_proj_sum + row2_1 + row2_3 + sing_Btilde_pt1_proj + sing_Btilde_proj_sum + row3_2 + sing_C_pt1_proj + sing_C_proj_sum

#	cav_row1 = np.load('ITHACA_SEM_data/cav_row1.npy', 'r')
#	cav_row3 = np.load('ITHACA_SEM_data/cav_row3.npy', 'r')

#	r_add_body_1 = Gr_factor * np.dot(np.transpose(c_f_bnd), cav_row1)
#	r_add_body_2 = Gr_factor * np.dot(np.transpose(c_f_int), cav_row3)

#	A_ROM_rhs_alt = r_add_body_1 + r_add_body_2

#	ROM_solve_alt_oo = np.linalg.solve(A_ROM_Mat_oo, A_ROM_rhs_alt)

	A_11_alt = mKinvis*A_11_1 + sing_A_MM_cropped_proj_sum
	A_31_alt = mKinvis*A_31_1 + sing_Btilde_cropped_proj_sum
	A_13_alt = mKinvis*A_13_1 + sing_B_MM_cropped_proj_sum
	A_33_alt = mKinvis*A_33_1 + sing_C_proj_sum
	A_ROM_alt = A_11_alt + A_21_alt + A_31_alt + A_12_alt + A_32_alt + A_13_alt + A_23_alt + A_33_alt
	A_ROM_rhs_alt = mKinvis*r1_1 + r2_1 + mKinvis*r3_1 + sing_A_MM_Aud_proj_sum + sing_Btilde_Aud_proj_sum
	ROM_solve_alt = np.linalg.solve(A_ROM_alt, A_ROM_rhs_alt)
	reproj_ROM_f_bnd_alt = np.dot(V_f_bnd, ROM_solve_alt)
	reproj_ROM_f_p_alt = np.dot(V_f_p, ROM_solve_alt)
	reproj_ROM_f_int_alt = np.dot(V_f_int, ROM_solve_alt)

#	reproj_ROM_f_bnd_alt_oo = np.dot(c_f_bnd, ROM_solve_alt_oo)
#	reproj_ROM_f_p_alt_oo = np.dot(c_f_p, ROM_solve_alt_oo)
#	reproj_ROM_f_int_alt_oo = np.dot(c_f_int, ROM_solve_alt_oo)

#	f_bnd = reproj_ROM_f_bnd_alt_oo
#	f_int = reproj_ROM_f_int_alt_oo

	cnt = 0
	cnt1 = 0
	offset = 0
	nz_loc = 1
	nvel = 2
	nplanecoeffs = nlc
	nbnd = nbmap
	nint = nimap
	fields = np.zeros( [nvel, nplanecoeffs] )

	f_bnd_ROM = 0.0*np.arange(c_f_bnd.shape[0]*1.0)
	cropped_counter = 0
	non_cropped_counter = 0
	for i in range(0, c_f_bnd.shape[0]):
		if i not in loc_dbc_set:
			f_bnd_ROM[i] = reproj_ROM_f_bnd_alt[cropped_counter]
			cropped_counter = cropped_counter + 1
		else:
			f_bnd_ROM[i] = -f_bnd_dbc[non_cropped_counter]
			non_cropped_counter += 1

	f_bnd = -f_bnd_ROM
	f_int = -reproj_ROM_f_int_alt

	for curr_elem in range(0, num_elem):
		for j in range(0, nvel):
			for k in range(0, nbnd):
				fields[j, offset + int(bmap[k])] = f_bnd[cnt + k]
			for k in range(0, nint):
				fields[j, offset + int(imap[k])] = f_int[cnt1+k]
			cnt += nbnd
			cnt1 += nint
		offset += nbnd + nint # is 169

	cnt = 0
	offset = 0
	velo0 = np.zeros( [ngc] )
	velo1 = np.zeros( [ngc] )

#	for i in range(0, nplanecoeffs):
#		velo0[int(LocGloMapA[i])] = int(LocGloSignA[i]) * fields[0,i]
#		velo1[int(LocGloMapA[i])] = int(LocGloSignA[i]) * fields[1,i]

	velo0 = np.dot(LocGloMapMatA, np.transpose(fields[0,:]))
	velo1 = np.dot(LocGloMapMatA, np.transpose(fields[1,:]))


	return (velo0, velo1)
