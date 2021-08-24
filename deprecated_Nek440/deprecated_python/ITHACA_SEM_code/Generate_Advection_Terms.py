
import numpy as np

def remove_from_matrix(matrix, columns, rows):
	return [
		[float(matrix[row_num][col_num])
		for col_num in range(len(matrix[row_num]))
		if not col_num in columns]

		for row_num in range(len(matrix))
		if not row_num in rows]

def gen(u_x, u_y, M_trafo_no_pp_incl_dbc, V_f_bnd, V_f_int, f_bnd_dbc, elem_not_loc_dbc, loc_dbc_set, c_f_bnd, glodofphys, IP, IPp, bwdtrans, cartmap0, cartmap1, LocGloSign, LocGloMap, LocGloMapMatA, bmap, imap, LocGloBndSign, LocGloBndMap, BndCondCoeffsToGlobalCoeffsMap):
	
#	glodofphys = np.load('glodofphys.npy', 'r')
#	IP = np.load('IP.npy', 'r')
#	IPp = np.load('IPp.npy', 'r')
#	bwdtrans = np.load('bwdtrans.npy', 'r')
#	cartmap0 = np.load('cartmap0.npy', 'r')
#	cartmap1 = np.load('cartmap1.npy', 'r')
#	LocGloSign = np.load('LocGloSign.npy', 'r')
#	LocGloMap = np.load('LocGloMap.npy', 'r')
#	LocGloSignA = np.load('LocGloSignA.npy', 'r')
#	LocGloMapA = np.load('LocGloMapA.npy', 'r')
#	bmap = (np.load('bmap.npy', 'r'))
#	imap = (np.load('imap.npy', 'r'))

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

	Ah_ele = np.zeros( [nsize_bndry_p1, nsize_bndry_p1] )
	Ah_elem = np.zeros( [num_elem, nsize_bndry_p1, nsize_bndry_p1] )

	B_elem = np.zeros( [num_elem, nsize_bndry, nsize_int] )
	C_elem = np.zeros( [num_elem, nsize_bndry, nsize_int] )
	D_elem = np.zeros( [num_elem, nsize_int, nsize_int] )
	Bh_elem = np.zeros( [num_elem, nsize_bndry_p1, nsize_p_m1] )
	Ch_elem = np.zeros( [num_elem, nsize_p_m1, nsize_bndry_p1] )
	Dh_elem = np.zeros( [num_elem, nsize_p_m1, nsize_p_m1] )

	sing_A = np.zeros([num_elem*nsize_bndry,num_elem*nsize_bndry])
	sing_B = np.zeros([num_elem*nsize_bndry,num_elem*nsize_int])
	sing_Btilde = np.zeros([num_elem*nsize_bndry,num_elem*nsize_int])
	sing_C = np.zeros([num_elem*nsize_int,num_elem*nsize_int])

	for curr_elem in range(0, num_elem):
		Ah_ele_vec = 0*np.arange(Ahrows*Ahrows*1.0)
		H1_bnd_ele_vec = 0*np.arange(Ahrows*Ahrows*1.0)
		B_ele_vec = 0*np.arange(nsize_bndry*nsize_int*1.0)
		C_ele_vec = 0*np.arange(nsize_bndry*nsize_int*1.0)
		D_ele_vec = 0*np.arange(nsize_int*nsize_int*1.0)
		loc_bwd_mat = bwdtrans[(curr_elem*ncoeffs) : (curr_elem*ncoeffs + ncoeffs) ,:]
		loc_cm0 = cartmap0[(curr_elem*nphys) : (curr_elem*nphys + nphys) ,:]
		loc_cm1 = cartmap1[(curr_elem*nphys) : (curr_elem*nphys + nphys) ,:]
		loc_IP = IP[(curr_elem*nphys) : (curr_elem*nphys + nphys) ,:]
		loc_IPp = IPp[(curr_elem*nphys) : (curr_elem*nphys + nphys) ,:]
		HelmMat = HelmMats[curr_elem,:]

		for i in range(0, nbmap):
			coeffs = 0*np.arange(ncoeffs*1.0)
			coeffs[int(bmap[i])] = 1.0
			phys = np.dot(coeffs, loc_bwd_mat)
			
			for k in range(0, 2):
				if (k == 0):
					deriv = np.dot(phys, (loc_cm0))
					tmpphys = u_x[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv
					coeffs = np.dot(tmpphys, loc_IP)
					for nv in range(0, 2):
						for j in range(0, nbmap):
							Ah_ele_vec[ j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += coeffs[int(bmap[j])]
						for j in range(0, nimap):
							C_ele_vec[ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += coeffs[int(imap[j])] 
					
				if (k == 1):
					deriv = np.dot(phys, (loc_cm1))
					tmpphys = u_y[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv
					coeffs = np.dot(tmpphys, loc_IP)
					for nv in range(0, 2):
						for j in range(0, nbmap):
							Ah_ele_vec[ j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += coeffs[int(bmap[j])]
						for j in range(0, nimap):
							C_ele_vec[ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += coeffs[int(imap[j])]
		
		for i in range(0, nimap):
			coeffs = 0*np.arange(ncoeffs*1.0)
			coeffs[int(imap[i])] = 1.0
			phys = np.dot(coeffs, loc_bwd_mat)
			for k in range(0, 2):
				if (k == 0):
					deriv = np.dot(phys, (loc_cm0))
					tmpphys = u_x[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv
					coeffs = np.dot(tmpphys, loc_IP)
					for nv in range(0, 2):
						for j in range(0, nbmap):
							B_ele_vec[ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += coeffs[int(bmap[j])]
						for j in range(0, nimap):
							D_ele_vec[ j+nv*nimap + (i+nv*nimap)*nsize_int ] += coeffs[int(imap[j])] 
				if (k == 1):
					deriv = np.dot(phys, (loc_cm1))
					tmpphys = u_y[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv
					coeffs = np.dot(tmpphys, loc_IP)
					for nv in range(0, 2):
						for j in range(0, nbmap):
							B_ele_vec[ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += coeffs[int(bmap[j])]
						for j in range(0, nimap):
							D_ele_vec[ j+nv*nimap + (i+nv*nimap)*nsize_int ] += coeffs[int(imap[j])] 
			
			


		# write Ah_elem
		for i in range(0, nsize_bndry_p1):
			for j in range(0, nsize_bndry_p1):
				Ah_elem[curr_elem, i, j] = Ah_ele_vec[ i + j*Ahrows ]

			
		for i in range(0, nsize_bndry):
			for j in range(0, nsize_int):
				B_elem[curr_elem, i, j] = B_ele_vec[ i + j*nsize_bndry ]			
			
		for i in range(0, nsize_bndry):
			for j in range(0, nsize_int):
				C_elem[curr_elem, i, j] = C_ele_vec[ i + j*nsize_bndry ]				
								 
		for i in range(0, nsize_int):
			for j in range(0, nsize_int):
				D_elem[curr_elem, i, j] = D_ele_vec[ i + j*nsize_int ]


		Ah_elem_cp = Ah_elem.copy()
		sing_A[curr_elem*nsize_bndry:curr_elem*nsize_bndry+nsize_bndry, curr_elem*nsize_bndry:curr_elem*nsize_bndry+nsize_bndry] = Ah_elem_cp[curr_elem, 0:nsize_bndry_p1-1, 0:nsize_bndry_p1-1] 		
		sing_B[curr_elem*nsize_bndry:curr_elem*nsize_bndry+nsize_bndry, curr_elem*nsize_int:curr_elem*nsize_int+nsize_int] = B_elem[curr_elem,:,:]
		sing_Btilde[curr_elem*nsize_bndry:curr_elem*nsize_bndry+nsize_bndry, curr_elem*nsize_int:curr_elem*nsize_int+nsize_int] = C_elem[curr_elem,:,:]
		sing_C[curr_elem*nsize_int:curr_elem*nsize_int+nsize_int, curr_elem*nsize_int:curr_elem*nsize_int+nsize_int] = D_elem[curr_elem,:,:]

	sing_A_M = np.dot(np.transpose(M_trafo_no_pp_incl_dbc), sing_A)
	sing_A_MM = np.dot(M_trafo_no_pp_incl_dbc, sing_A_M)
	sing_A_MM_cropped = np.array(remove_from_matrix(sing_A_MM, loc_dbc_set, loc_dbc_set))
	sing_A_MM_Aud = np.array(remove_from_matrix(sing_A_MM, elem_not_loc_dbc, loc_dbc_set))
	sing_A_MM_cropped_proj = np.dot(np.dot(np.transpose(V_f_bnd), sing_A_MM_cropped), V_f_bnd)
	sing_A_MM_Aud_proj = np.dot(np.transpose(V_f_bnd), np.dot(sing_A_MM_Aud, f_bnd_dbc)) 

	sing_B_M = np.dot(np.transpose(M_trafo_no_pp_incl_dbc),sing_B)
	sing_B_MM = np.dot(M_trafo_no_pp_incl_dbc, sing_B_M)
	sing_B_MM_cropped = np.array(remove_from_matrix(sing_B_MM, set(), loc_dbc_set))
	sing_B_MM_cropped_proj = np.dot(np.dot(np.transpose(V_f_bnd), sing_B_MM_cropped), V_f_int)

	Btilde2_cropped = np.array(remove_from_matrix(sing_Btilde,  set(), loc_dbc_set))
	Btilde2_Aud = np.array(remove_from_matrix(sing_Btilde, set(), elem_not_loc_dbc))
	Btilde2_cropped_proj = np.dot(np.dot(np.transpose(V_f_int), np.transpose(Btilde2_cropped)), V_f_bnd)
	Btilde2_Aud_proj = np.dot(np.transpose(V_f_int), np.dot(np.transpose(Btilde2_Aud), f_bnd_dbc))

	sing_C_proj = np.dot(np.dot(np.transpose(V_f_int), sing_C), V_f_int)
	sing_A_proj = np.dot(np.dot(np.transpose(c_f_bnd), sing_A), c_f_bnd)
	sing_B_proj = np.dot(np.dot(np.transpose(c_f_bnd), sing_B), V_f_int)
	sing_Btilde_proj = np.dot(np.dot(np.transpose(V_f_int), np.transpose(sing_Btilde)), c_f_bnd)


	return (sing_A_MM_cropped_proj, sing_A_MM_Aud_proj, sing_B_MM_cropped_proj, Btilde2_cropped_proj, Btilde2_Aud_proj, sing_C_proj, sing_A_proj, sing_B_proj, sing_Btilde_proj)
