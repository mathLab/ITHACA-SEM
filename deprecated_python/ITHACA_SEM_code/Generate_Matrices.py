
import numpy as np

def gen():

	print("Generating Parameter-independent terms.")

#	Dbnd = np.load('ITHACA_SEM_data/Dbnd_cav.npy', 'r')
#	Dint = np.load('ITHACA_SEM_data/Dint_cav.npy', 'r')
	glodofphys = np.load('ITHACA_SEM_data/glodofphys.npy', 'r')
#	physglodof = np.load('ITHACA_SEM_data/physglodof.npy', 'r')
	LocGloBndSign = np.load('ITHACA_SEM_data/LocGloBndSign.npy', 'r')
	LocGloBndMap = np.load('ITHACA_SEM_data/LocGloBndMap.npy', 'r')
	LocGloSign = np.load('ITHACA_SEM_data/LocGloSign.npy', 'r')
	LocGloMap = np.load('ITHACA_SEM_data/LocGloMap.npy', 'r')
#	LocGloSignA = np.load('ITHACA_SEM_data/LocGloSignA.npy', 'r')
#	LocGloMapA = np.load('ITHACA_SEM_data/LocGloMapA.npy', 'r')
	IP = np.load('ITHACA_SEM_data/IP.npy', 'r')
	IPp = np.load('ITHACA_SEM_data/IPp.npy', 'r')
	bwdtrans = np.load('ITHACA_SEM_data/bwdtrans.npy', 'r')
	cartmap0 = np.load('ITHACA_SEM_data/cartmap0.npy', 'r')
	cartmap1 = np.load('ITHACA_SEM_data/cartmap1.npy', 'r')
	bmap = (np.load('ITHACA_SEM_data/bmap.npy', 'r'))
	imap = (np.load('ITHACA_SEM_data/imap.npy', 'r'))
	BndCondCoeffsToGlobalCoeffsMap = np.load('ITHACA_SEM_data/BndCondCoeffsToGlobalCoeffsMap.npy', 'r')
	numDirBnd = int(np.load('ITHACA_SEM_data/numDirBnd.npy', 'r'))
	LocGloMapMatA = np.load('ITHACA_SEM_data/LocGloMapMatA.npy', 'r')
	force0 = np.load('ITHACA_SEM_data/forcing0.npy', 'r')
	force1 = np.load('ITHACA_SEM_data/forcing1.npy', 'r')
	forcing = np.transpose(np.c_[force0, force1])
	bndcond_k0_i_0 = np.load('ITHACA_SEM_data/bndcond_k0_i_0.npy', 'r')
	bndcond_k0_i_1 = np.load('ITHACA_SEM_data/bndcond_k0_i_1.npy', 'r')
	bndcond_k0_i_2 = np.load('ITHACA_SEM_data/bndcond_k0_i_2.npy', 'r')
	bndcond_k1_i_0 = np.load('ITHACA_SEM_data/bndcond_k1_i_0.npy', 'r')
	bndcond_k1_i_1 = np.load('ITHACA_SEM_data/bndcond_k1_i_1.npy', 'r')
	bndcond_k1_i_2 = np.load('ITHACA_SEM_data/bndcond_k1_i_2.npy', 'r')


	print('BndCondCoeffsToGlobalCoeffsMap.shape ', BndCondCoeffsToGlobalCoeffsMap.shape)
#	print('LocGloBndMap.shape ', LocGloBndMap.shape)
#	print('maxi LocGloBndMap ', np.asarray(np.rint(np.amax(LocGloBndMap)), dtype=int) )
#	print('LocGloMap.shape ', LocGloMap.shape)
#	print('IP.shape ', IP.shape)
#	print('IPp.shape ', IPp.shape)



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
	
#	testsimu_D = (np.load('ITHACA_SEM_data/testsimu_D.npy', 'r'))
#	testsimu_nD = (np.load('ITHACA_SEM_data/testsimu_nD.npy', 'r'))		
#	testsimu_time_D = (np.load('ITHACA_SEM_data/testsimu_time_D.npy', 'r'))		
#	testsimu_time_nD = (np.load('ITHACA_SEM_data/testsimu_time_nD.npy', 'r'))	
#	print('testsimu_D.shape ', testsimu_D.shape)
#	print('testsimu_nD.shape ', testsimu_nD.shape)
#	print('testsimu_time_D.shape ', testsimu_time_D.shape)
#	print('testsimu_time_nD.shape ', testsimu_time_nD.shape)
#	print('glodofphys.shape ', glodofphys.shape)
	ngc = glodofphys.shape[0]	
	
#	f_tt_cd = open('ITHACA_SEM_data/tt_corr_dofs.txt', 'r')
#	tt_cd = np.loadtxt(f_tt_cd)
#	tt_cd = np.asarray(np.rint(tt_cd), dtype=int)
	# find inverse tt_cd mapping
#	inv_tt_cd = 0*tt_cd[0:ngc]
#	for i in range(0,ngc):
#		inv_tt_cd[tt_cd[i]] = int(i)	
	
#	f_tt_check = open('ITHACA_SEM_data/time_traj.txt', 'r')
#	tt_check = np.loadtxt(f_tt_check)	
#	print('tt_check.shape ', tt_check.shape)
#	nn = np.dot(tt_check[:,tt_cd[0:ngc]], glodofphys)
#	print(nn[:,1296])
	
#	curr_snap_y_loc = np.dot(testsimu[inv_tt_cd[ngc:2*ngc]], glodofphys)

	
#	curr_snap_x_loc_D = np.dot(testsimu_D[0:ngc], glodofphys)
#	print("init qoi D ", curr_snap_x_loc_D[1296])	
#	curr_snap_x_loc_nD = np.dot(testsimu_nD[0:ngc], glodofphys)
#	print("init qoi nD ", curr_snap_x_loc_nD[1296])	

#	testsimu_time_D = np.asarray(np.rint(testsimu_time_D), dtype=int)
#	testsimu_time_nD = np.asarray(np.rint(testsimu_time_nD), dtype=int)
	# find inverse tt_cd mapping
#	inv_tt_D = 0*testsimu_time_D[0:ngc]
#	inv_tt_nD = 0*testsimu_time_nD[0:ngc]
#	for i in range(0,ngc):
#		print("i ",i)
#		print(testsimu_time_nD[i])	
#		inv_tt_D[testsimu_time_D[i]-1] = int(i)
#		inv_tt_nD[testsimu_time_nD[i]-1] = int(i)

#	curr_snap_x_loc_time_inv_D = np.dot(testsimu_time_D[inv_tt_D[0:ngc]], glodofphys)
#	print("init qoi inv D ", curr_snap_x_loc_time_inv_D[1296])	
#	curr_snap_x_loc_time_inv_nD = np.dot(testsimu_time_nD[inv_tt_nD[0:ngc]], glodofphys)
#	print("init qoi inv nD ", curr_snap_x_loc_time_inv_nD[1296])	

	
#	f_snap = open('ITHACA_SEM_data/Acav1.txt', 'r')
#	snap = np.transpose(np.loadtxt(f_snap))

	# call to generate here the c_f_bnd, c_f_p and c_f_int
#	c_f_bnd = np.load("ITHACA_SEM_data/trafo_f_bnd_cav.npy", "r")
#	c_f_p = np.load("ITHACA_SEM_data/trafo_f_p_cav.npy", "r")
#	c_f_int = np.load("ITHACA_SEM_data/trafo_f_int_cav.npy", "r")

	NumDirBCs = numDirBnd

#	nGlobHomBndDofs = 1161 # for cavity
#	nBndDofs = nGlobHomBndDofs + NumDirBCs

#	print("nBndDofs ", nBndDofs)

	nBndDofs = np.asarray(np.rint(np.amax(LocGloBndMap)), dtype=int) + 1
	nGlobHomBndDofs = nBndDofs - NumDirBCs

	nGlobBndDofs = nBndDofs
	rows = nGlobHomBndDofs

	init_bnd = 0.0*np.arange(nBndDofs*1.0) 
	bndcnt = 0
	for j in range(0, bndcond_k0_i_0.shape[0]):
		init_bnd[int(BndCondCoeffsToGlobalCoeffsMap[bndcnt])] = bndcond_k0_i_0[j]
		bndcnt += 1
	for j in range(0, bndcond_k0_i_1.shape[0]):
		init_bnd[int(BndCondCoeffsToGlobalCoeffsMap[bndcnt])] = bndcond_k0_i_1[j]
		bndcnt += 1
	for j in range(0, bndcond_k0_i_2.shape[0]):
		init_bnd[int(BndCondCoeffsToGlobalCoeffsMap[bndcnt])] = bndcond_k0_i_2[j]
		bndcnt += 1
	for j in range(0, bndcond_k1_i_0.shape[0]):
		init_bnd[int(BndCondCoeffsToGlobalCoeffsMap[bndcnt])] = bndcond_k1_i_0[j]
		bndcnt += 1
	for j in range(0, bndcond_k1_i_1.shape[0]):
		init_bnd[int(BndCondCoeffsToGlobalCoeffsMap[bndcnt])] = bndcond_k1_i_1[j]
		bndcnt += 1
	for j in range(0, bndcond_k1_i_2.shape[0]):
		init_bnd[int(BndCondCoeffsToGlobalCoeffsMap[bndcnt])] = bndcond_k1_i_2[j]
		bndcnt += 1



	cnt = 0
	offset = 0
	no_loc_dbc = 0
	no_not_loc_dbc = 0
	elem_loc_dbc = set()
	elem_not_loc_dbc = set()
	for curr_elem in range(0, num_elem):
		for j in range(0, nvel):
			for k in range(0, nbnd):
				if (int(LocGloMap[offset+j*nbnd+k]) < NumDirBCs):
					no_loc_dbc = no_loc_dbc + 1
					elem_loc_dbc.update({cnt+k})
				else:
					no_not_loc_dbc = no_not_loc_dbc + 1
					elem_not_loc_dbc.update({cnt+k})
			cnt += nbnd
		offset += nvel*nbnd + nint

	no_of_dbc_in_loc = no_loc_dbc
	loc_dbc_set = elem_loc_dbc


	M_trafo_no_pp_incl_dbc = np.zeros([ num_elem*nsize_bndry , nBndDofs ])
	for curr_elem in range(0, num_elem):
		cnt = curr_elem*nsize_bndry_p1
		cnt_no_pp = curr_elem*nsize_bndry
		for i in range(0, nsize_bndry_p1):
			gid1 = int(LocGloBndMap[cnt + i])
			sign1 = int(LocGloBndSign[cnt + i])
			if (gid1 >= 0):
				if (i < nsize_bndry):
					M_trafo_no_pp_incl_dbc[ cnt_no_pp + i, gid1 ] = sign1




#	f_tt_cd = open('ITHACA_SEM_data/tt_corr_dofs.txt', 'r')
#	tt_cd = np.loadtxt(f_tt_cd)
#	tt_cd = np.asarray(np.rint(tt_cd), dtype=int)
	# find inverse tt_cd mapping
#	inv_tt_cd = 0*tt_cd[0:ngc]
#	for i in range(0,ngc):
#		inv_tt_cd[tt_cd[i]] = int(i)

	Ah_ele = np.zeros( [nsize_bndry_p1, nsize_bndry_p1] )
	Ah_elem = np.zeros( [num_elem, nsize_bndry_p1, nsize_bndry_p1] )

	B_elem = np.zeros( [num_elem, nsize_bndry, nsize_int] )
	C_elem = np.zeros( [num_elem, nsize_bndry, nsize_int] )
	D_elem = np.zeros( [num_elem, nsize_int, nsize_int] )
	Bh_elem = np.zeros( [num_elem, nsize_bndry_p1, nsize_p_m1] )
	Ch_elem = np.zeros( [num_elem, nsize_p_m1, nsize_bndry_p1] )
	Dh_elem = np.zeros( [num_elem, nsize_p_m1, nsize_p_m1] )
	Dbnd_elem = np.zeros( [num_elem, nsize_p, nsize_bndry] )
	Dint_elem = np.zeros( [num_elem, nsize_p, nsize_int] )

	sing_A = np.zeros([num_elem*nsize_bndry,num_elem*nsize_bndry])
	sing_B = np.zeros([num_elem*nsize_bndry,num_elem*nsize_int])
	sing_Btilde = np.zeros([num_elem*nsize_bndry,num_elem*nsize_int])
	sing_C = np.zeros([num_elem*nsize_int,num_elem*nsize_int])
	sing_Dbnd = np.zeros([num_elem*nsize_p,num_elem*nsize_bndry])
	sing_Dint = np.zeros([num_elem*nsize_p,num_elem*nsize_int])

	for curr_elem in range(0, num_elem):
		Ah_ele_vec = 0*np.arange(Ahrows*Ahrows*1.0)
		H1_bnd_ele_vec = 0*np.arange(Ahrows*Ahrows*1.0)
		B_ele_vec = 0*np.arange(nsize_bndry*nsize_int*1.0)
		C_ele_vec = 0*np.arange(nsize_bndry*nsize_int*1.0)
		D_ele_vec = 0*np.arange(nsize_int*nsize_int*1.0)
		Dbnd_ele_vec = 0*np.arange(nsize_bndry*nsize_p*1.0)
		Dint_ele_vec = 0*np.arange(nsize_int*nsize_p*1.0)
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
				for j in range(0, nbmap):
					Ah_ele_vec[ i+k*nbmap + (j+k*nbmap)*Ahrows ] += HelmMat[int(bmap[i]) + HelmMatRows*int(bmap[j])]
				for j in range(0, nimap):
					B_ele_vec[i+k*nbmap + (j+k*nimap)*nsize_bndry] += HelmMat[int(bmap[i]) + HelmMatRows*int(imap[j])]

			for k in range(0, 2):
				if (k == 0):
					deriv = np.dot(phys, (loc_cm0))
					pcoeff = np.dot(deriv, loc_IPp)
					Dbnd_ele_vec[ (k*nbmap + i)*nsize_p : ((k*nbmap + i)*nsize_p + nsize_p) ] = pcoeff
					
				if (k == 1):
					deriv = np.dot(phys, (loc_cm1))
					pcoeff = np.dot(deriv, loc_IPp) 
					Dbnd_ele_vec[ (k*nbmap + i)*nsize_p : ((k*nbmap + i)*nsize_p + nsize_p) ] = pcoeff			

		for i in range(0, nimap):
			coeffs = 0*np.arange(ncoeffs*1.0)
			coeffs[int(imap[i])] = 1.0
			phys = np.dot(coeffs, loc_bwd_mat)
			for k in range(0, 2):
				for j in range(0, nbmap):
					C_ele_vec[j+k*nbmap + (i+k*nimap)*nsize_bndry] += HelmMat[int(imap[i])+HelmMatRows*int(bmap[j])]
				for j in range(0, nimap):
					D_ele_vec[i+k*nimap + (j+k*nimap)*nsize_int] += HelmMat[int(imap[i])+HelmMatRows*int(imap[j])];
					
			for k in range(0, 2):
				if (k == 0):
					deriv = np.dot(phys, (loc_cm0))
					pcoeff = np.dot(deriv, loc_IPp)
					Dint_ele_vec[ (k*nimap + i)*nsize_p : ((k*nimap + i)*nsize_p + nsize_p) ] = pcoeff
				if (k == 1):
					deriv = np.dot(phys, (loc_cm1))
					pcoeff = np.dot(deriv, loc_IPp)
					Dint_ele_vec[ (k*nimap + i)*nsize_p : ((k*nimap + i)*nsize_p + nsize_p) ] = pcoeff

						
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
			
		for i in range(0, nsize_p):
			for j in range(0, nsize_bndry):
				Dbnd_elem[curr_elem, i, j] = Dbnd_ele_vec[ i + j*nsize_p ]	
			
		for i in range(0, nsize_p):
			for j in range(0, nsize_int):
				Dint_elem[curr_elem, i, j] = Dint_ele_vec[ i + j*nsize_p ]

			
		Ah_elem_cp = Ah_elem.copy()
		sing_A[curr_elem*nsize_bndry:curr_elem*nsize_bndry+nsize_bndry, curr_elem*nsize_bndry:curr_elem*nsize_bndry+nsize_bndry] = Ah_elem_cp[curr_elem, 0:nsize_bndry_p1-1, 0:nsize_bndry_p1-1] 		
		sing_B[curr_elem*nsize_bndry:curr_elem*nsize_bndry+nsize_bndry, curr_elem*nsize_int:curr_elem*nsize_int+nsize_int] = B_elem[curr_elem,:,:]
		sing_Btilde[curr_elem*nsize_bndry:curr_elem*nsize_bndry+nsize_bndry, curr_elem*nsize_int:curr_elem*nsize_int+nsize_int] = C_elem[curr_elem,:,:]
		sing_C[curr_elem*nsize_int:curr_elem*nsize_int+nsize_int, curr_elem*nsize_int:curr_elem*nsize_int+nsize_int] = D_elem[curr_elem,:,:]
		sing_Dbnd[curr_elem*nsize_p:curr_elem*nsize_p+nsize_p, curr_elem*nsize_bndry:curr_elem*nsize_bndry+nsize_bndry] = Dbnd_elem[curr_elem,:,:]
		sing_Dint[curr_elem*nsize_p:curr_elem*nsize_p+nsize_p, curr_elem*nsize_int:curr_elem*nsize_int+nsize_int] = Dint_elem[curr_elem,:,:]



	return (sing_Dbnd, sing_Dint, sing_A, sing_B, sing_Btilde, sing_C, M_trafo_no_pp_incl_dbc, glodofphys, LocGloSign, LocGloMap, ngc, nlc, bmap, imap, num_elem, nsize_bndry, nsize_int, nsize_p, nvel, NumDirBCs, nGlobHomBndDofs, no_loc_dbc, elem_loc_dbc, no_not_loc_dbc, elem_not_loc_dbc, IP, IPp, bwdtrans, cartmap0, cartmap1, LocGloBndSign, LocGloBndMap, BndCondCoeffsToGlobalCoeffsMap, LocGloMapMatA, init_bnd, forcing)
