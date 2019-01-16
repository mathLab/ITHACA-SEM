
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



def Oseen_step_rom(u_x, u_y, curr_para, Dbnd, Dint, sing_A, sing_B, sing_Btilde, sing_C, M_trafo_no_pp_incl_dbc, glodofphys, LocGloSign, LocGloMap, ngc, nlc, bmap, imap, num_elem, nsize_bndry, nsize_int, nsize_p, nvel, NumDirBCs, nGlobHomBndDofs, no_loc_dbc, elem_loc_dbc, no_not_loc_dbc, elem_not_loc_dbc, IP, IPp, bwdtrans, cartmap0, cartmap1, LocGloBndSign, LocGloBndMap, BndCondCoeffsToGlobalCoeffsMap, LocGloMapMatA, init_bnd, forcing_in):

	mKinvis = curr_para[0]
	Gr_factor = curr_para[1]
	
	forcing = np.multiply(Gr_factor, forcing_in)

	np.set_printoptions(threshold=np.inf)	

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
#	print("nlc ", nlc)
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
	Ah_elem_pt1 = np.zeros( [num_elem, nsize_bndry_p1, nsize_bndry_p1] )
	Ah_elem_pt2 = np.zeros( [num_elem, nsize_bndry_p1, nsize_bndry_p1] )
	Ah_elem_pt21 = np.zeros( [num_elem, nsize_bndry_p1, nsize_p_m1] )
	Ah_elem_pt22 = np.zeros( [num_elem, nsize_p_m1, nsize_p_m1] )
	Ah_elem_pt22_noK = np.zeros( [num_elem, nsize_p_m1, nsize_p_m1] )
	Ah_elem_pt23 = np.zeros( [num_elem, nsize_p_m1, nsize_bndry_p1] )
	Ah_elem_pt1_bf = np.zeros( [num_elem, nsize_bndry_p1-1, nsize_bndry_p1-1] )
	Ah_elem_pt1_bf_pt1 = np.zeros( [num_elem, nsize_bndry_p1-1, nsize_bndry_p1-1] )
	Ah_elem_pt1_bf_pt21 = np.zeros( [num_elem, nsize_bndry, nsize_int] )
	Ah_elem_pt1_bf_pt22 = np.zeros( [num_elem, nsize_int, nsize_int] )
	Ah_elem_pt1_bf_pt23 = np.zeros( [num_elem, nsize_int, nsize_bndry] )
	Ah_elem_pt21_bf = np.zeros( [num_elem, nsize_bndry, nsize_p] )
	Ah_elem_pt22_bf = np.zeros( [num_elem, nsize_p, nsize_p] )
	Ah_elem_pt23_bf = np.zeros( [num_elem, nsize_p, nsize_bndry] )

	B_elem = np.zeros( [num_elem, nsize_bndry, nsize_int] )
	C_elem = np.zeros( [num_elem, nsize_bndry, nsize_int] )
	D_elem = np.zeros( [num_elem, nsize_int, nsize_int] )
	DnoK_elem = np.zeros( [num_elem, nsize_int, nsize_int] )
	Dadv_elem = np.zeros( [num_elem, nsize_int, nsize_int] )
	Dadv_elem_alt = np.zeros( [num_elem, nsize_int, nsize_int] )
	Dbnd_elem = np.zeros( [num_elem, nsize_p, nsize_bndry] )
	Dint_elem = np.zeros( [num_elem, nsize_p, nsize_int] )

	Bh_elem = np.zeros( [num_elem, nsize_bndry_p1, nsize_p_m1] )
	Ch_elem = np.zeros( [num_elem, nsize_p_m1, nsize_bndry_p1] )
	Dh_elem = np.zeros( [num_elem, nsize_p_m1, nsize_p_m1] )

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
		DnoK_ele_vec = 0*np.arange(nsize_int*nsize_int*1.0)
		Dadv_ele_vec = 0*np.arange(nsize_int*nsize_int*1.0)
		H1_int_ele_vec = 0*np.arange(nsize_int*nsize_int*1.0)
		D1_ele_vec = 0*np.arange(nsize_int*nsize_int*1.0)
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
					Ah_ele_vec[ i+k*nbmap + (j+k*nbmap)*Ahrows ] += mKinvis *  HelmMat[int(bmap[i]) + HelmMatRows*int(bmap[j])]
				for j in range(0, nimap):
					B_ele_vec[i+k*nbmap + (j+k*nimap)*nsize_bndry] += mKinvis * HelmMat[int(bmap[i]) + HelmMatRows*int(imap[j])]
			
			for k in range(0, 2):
				if (k == 0):
					deriv = np.dot(phys, (loc_cm0))
					tmpphys = u_x[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv
					coeffs = np.dot(tmpphys, loc_IP)
					for nv in range(0, 2):
						for j in range(0, nbmap):
							Ah_ele_vec[ j+nv*nbmap + (i+nv*nbmap)*Ahrows ] = Ah_ele_vec[ j+nv*nbmap + (i+nv*nbmap)*Ahrows ]  + coeffs[int(bmap[j])]

						for j in range(0, nimap):
							C_ele_vec[ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += coeffs[int(imap[j])] 
					pcoeff = np.dot(deriv, loc_IPp)
					Dbnd_ele_vec[ (k*nbmap + i)*nsize_p : ((k*nbmap + i)*nsize_p + nsize_p) ] = pcoeff
					
				if (k == 1):
					deriv = np.dot(phys, (loc_cm1))
					tmpphys = u_y[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv
					coeffs = np.dot(tmpphys, loc_IP)
					for nv in range(0, 2):
						for j in range(0, nbmap):
							Ah_ele_vec[ j+nv*nbmap + (i+nv*nbmap)*Ahrows ] = Ah_ele_vec[ j+nv*nbmap + (i+nv*nbmap)*Ahrows ]  + coeffs[int(bmap[j])]
						for j in range(0, nimap):
							C_ele_vec[ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += coeffs[int(imap[j])]
					pcoeff = np.dot(deriv, loc_IPp) 
					Dbnd_ele_vec[ (k*nbmap + i)*nsize_p : ((k*nbmap + i)*nsize_p + nsize_p) ] = pcoeff
				

		
		for i in range(0, nimap):
			coeffs = 0*np.arange(ncoeffs*1.0)
			coeffs[int(imap[i])] = 1.0
			phys = np.dot(coeffs, loc_bwd_mat)
			for k in range(0, 2):
				for j in range(0, nbmap):
					C_ele_vec[j+k*nbmap + (i+k*nimap)*nsize_bndry] += mKinvis*HelmMat[int(imap[i])+HelmMatRows*int(bmap[j])]
				for j in range(0, nimap):
					D_ele_vec[i+k*nimap + (j+k*nimap)*nsize_int] += mKinvis*HelmMat[int(imap[i])+HelmMatRows*int(imap[j])];
					DnoK_ele_vec[i+k*nimap + (j+k*nimap)*nsize_int] += HelmMat[int(imap[i])+HelmMatRows*int(imap[j])];
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
							Dadv_ele_vec[ j+nv*nimap + (i+nv*nimap)*nsize_int ] += coeffs[int(imap[j])]
					pcoeff = np.dot(deriv, loc_IPp)
					Dint_ele_vec[ (k*nimap + i)*nsize_p : ((k*nimap + i)*nsize_p + nsize_p) ] = pcoeff
				if (k == 1):
					deriv = np.dot(phys, (loc_cm1))
					tmpphys = u_y[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv
					coeffs = np.dot(tmpphys, loc_IP)
					for nv in range(0, 2):
						for j in range(0, nbmap):
							B_ele_vec[ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += coeffs[int(bmap[j])]
						for j in range(0, nimap):
							D_ele_vec[ j+nv*nimap + (i+nv*nimap)*nsize_int ] += coeffs[int(imap[j])] 
							Dadv_ele_vec[ j+nv*nimap + (i+nv*nimap)*nsize_int ] += coeffs[int(imap[j])]
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
				DnoK_elem[curr_elem, i, j] = DnoK_ele_vec[ i + j*nsize_int ]
				Dadv_elem[curr_elem, i, j] = Dadv_ele_vec[ i + j*nsize_int ]

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
		sing_Dint[curr_elem*nsize_p:curr_elem*nsize_p+nsize_p, curr_elem*nsize_int:curr_elem*nsize_int+nsize_int] = Dint_elem[curr_elem,:,:]
		sing_Dbnd[curr_elem*nsize_p:curr_elem*nsize_p+nsize_p, curr_elem*nsize_bndry:curr_elem*nsize_bndry+nsize_bndry] = Dbnd_elem[curr_elem,:,:]
		
		
		
		Ah_elem_pt1_bf_pt22[curr_elem, :, :] = D_elem[curr_elem, :, :]
		Ah_elem_pt1_bf_pt21[curr_elem, :, :] = B_elem[curr_elem, :, :]
		Ah_elem_pt1_bf_pt23[curr_elem, :, :] = np.transpose(C_elem[curr_elem, :, :])
		Ah_elem_pt1_bf_pt1[curr_elem, :, :] = Ah_elem[curr_elem, 0:nsize_bndry_p1-1, 0:nsize_bndry_p1-1]


		D_elem[curr_elem, :, :] = np.linalg.inv( D_elem[curr_elem, :, :] )
		
		B_elem[curr_elem, :, :] = np.dot(B_elem[curr_elem, :, :] , D_elem[curr_elem, :, :])
		Ah_elem[curr_elem, 0:nsize_bndry_p1-1, 0:nsize_bndry_p1-1] = Ah_elem[curr_elem, 0:nsize_bndry_p1-1, 0:nsize_bndry_p1-1] - np.dot( B_elem[curr_elem, :, :] , np.transpose(C_elem[curr_elem, :, :]) )


	
		Elem_Cinv = D_elem[curr_elem, :, :]
		Elem_BCinv = B_elem[curr_elem, :, :]
		Elem_Btilde = C_elem[curr_elem, :, :]
		Elem_DintCinvDTint = np.dot(np.dot(Dint_elem[curr_elem, :, :], Elem_Cinv), np.transpose(Dint_elem[curr_elem, :, :]))
		Elem_BCinvDTint_m_DTbnd = np.dot(Elem_BCinv, np.transpose(Dint_elem[curr_elem, :, :])) - np.transpose(Dbnd_elem[curr_elem, :, :])
		Elem_DintCinvBTtilde_m_Dbnd = np.dot(np.dot( Dint_elem[curr_elem, :, :] , Elem_Cinv ), np.transpose(Elem_Btilde)) - Dbnd_elem[curr_elem, :, :]
	
		Ah_elem_pt22_bf[curr_elem, :, :] = - Elem_DintCinvDTint
		Ah_elem_pt21_bf[curr_elem, :, :] = Elem_BCinvDTint_m_DTbnd
		Ah_elem_pt23_bf[curr_elem, :, :] = Elem_DintCinvBTtilde_m_Dbnd
		Ah_elem_pt1_bf[curr_elem, :, :] = Ah_elem[curr_elem, 0:nsize_bndry_p1-1, 0:nsize_bndry_p1-1]

		

		Bh_curr_ele = np.zeros( [nsize_bndry_p1, nsize_p_m1] )
		Ch_curr_ele = np.zeros( [nsize_p_m1, nsize_bndry_p1] )
		Dh_curr_ele = np.zeros( [nsize_p_m1, nsize_p_m1] )
		#DhnoK_curr_ele = np.zeros( [nsize_p_m1, nsize_p_m1] )
	
		for i in range(0, nsize_p_m1):
			for j in range(0, nsize_p_m1):
				Dh_curr_ele[ i , j ] = - Elem_DintCinvDTint[i+1,j+1]
			
		for i in range(0, nsize_bndry_p1-1):
			Ah_elem[curr_elem, i, nsize_bndry_p1-1] = Elem_BCinvDTint_m_DTbnd[i, 0]
			Ah_elem[curr_elem, nsize_bndry_p1-1, i] = Elem_DintCinvBTtilde_m_Dbnd[0, i]
	
		Ah_elem[curr_elem, nsize_bndry_p1-1, nsize_bndry_p1-1] = -Elem_DintCinvDTint[0,0]
	
		for j in range(0, nsize_p_m1):
			for i in range(0, nsize_bndry_p1-1):
				Bh_curr_ele[i, j] = Elem_BCinvDTint_m_DTbnd[i,j+1]
				Ch_curr_ele[j, i] = Elem_DintCinvBTtilde_m_Dbnd[ j+1 , i ]
	
		for j in range(0, nsize_p_m1):
			Bh_curr_ele[nsize_bndry_p1-1, j] = - Elem_DintCinvDTint[0,j+1]
			Ch_curr_ele[j,nsize_bndry_p1-1] = - Elem_DintCinvDTint[j+1,0]
	
		Ah_elem_pt22[curr_elem, :, :] =  Dh_curr_ele
		Dh_curr_ele = np.linalg.inv( Dh_curr_ele )
		Bh_curr_ele = np.dot(Bh_curr_ele , Dh_curr_ele)
		Ah_elem_pt23[curr_elem, :, :] =  Ch_curr_ele 
		Ah_elem[curr_elem, :, :] = Ah_elem[curr_elem, :, :] - np.dot( Bh_curr_ele , Ch_curr_ele )
		Bh_elem[curr_elem, :, :] = Bh_curr_ele
		Ch_elem[curr_elem, :, :] = Ch_curr_ele
		Dh_elem[curr_elem, :, :] = Dh_curr_ele
 
	# now loading additionally Db1, Di1, A1, B1, Btilde1, C1, A2, B2, Btilde2, C2
	#print("np.linalg.norm(sing_A) ", np.linalg.norm(sing_A))
	#print("np.linalg.norm(mKinvis*A1 + A2) ", np.linalg.norm(mKinvis*A1 + A2))
	#print("np.linalg.norm(sing_A - mKinvis*A1 - A2) ", np.linalg.norm(sing_A - mKinvis*A1 - A2))
	#print("np.linalg.norm(sing_B) ", np.linalg.norm(sing_B))
	#print("np.linalg.norm(mKinvis*B1 + B2) ", np.linalg.norm(mKinvis*B1 + B2))
	#print("np.linalg.norm(sing_B - mKinvis*B1 - B2) ", np.linalg.norm(sing_B - mKinvis*B1 - B2))
	#print("np.linalg.norm(sing_Btilde) ", np.linalg.norm(sing_Btilde))
	#print("np.linalg.norm(mKinvis*Btilde1 + Btilde2) ", np.linalg.norm(mKinvis*Btilde1 + Btilde2))
	#print("np.linalg.norm(sing_Btilde - mKinvis*Btilde1 - Btilde2) ", np.linalg.norm(sing_Btilde - mKinvis*Btilde1 - Btilde2))
	#print("np.linalg.norm(sing_C) ", np.linalg.norm(sing_C))
	#print("np.linalg.norm(mKinvis*C1 + C2) ", np.linalg.norm(mKinvis*C1 + C2))
	#print("np.linalg.norm(sing_C - mKinvis*C1 - C2) ", np.linalg.norm(sing_C - mKinvis*C1 - C2))


	nBndDofs = nGlobHomBndDofs + NumDirBCs
	rows = nBndDofs - NumDirBCs
	my_Gmat = np.zeros( [rows, rows] )

	# build the transformation M
	M_trafo = np.zeros([num_elem*nsize_bndry + num_elem, nGlobHomBndDofs])
	M_trafo_no_pp = np.zeros([num_elem*nsize_bndry, nGlobHomBndDofs])

	for curr_elem in range(0, num_elem):
		cnt = curr_elem*nsize_bndry_p1
		cnt_no_pp = curr_elem*nsize_bndry
		loc_Ah = Ah_elem[curr_elem, :, :]
		for i in range(0, nsize_bndry_p1):
			gid1 = int(LocGloBndMap[cnt + i]) - NumDirBCs
			sign1 = int(LocGloBndSign[cnt + i])
			if (gid1 >= 0):
				M_trafo[ cnt + i, gid1 ] = sign1
				if (i < nsize_bndry):
					M_trafo_no_pp[ cnt_no_pp + i, gid1 ] = sign1
				for j in range(0, nsize_bndry_p1):
					gid2 = int(LocGloBndMap[cnt + j]) - NumDirBCs
					sign2 = int(LocGloBndSign[cnt + j])
					if (gid2 >= 0):
						my_Gmat[gid1,gid2] += sign1*sign2*loc_Ah[i,j]
				

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


	print(num_elem*nsize_bndry)
	print(num_elem*nsize_int)
	print(num_elem*nsize_p)
	nGlobBndDofs = nGlobHomBndDofs + NumDirBCs

	f_bnd = 0*np.arange(num_elem*nsize_bndry*1.0) 
	f_int = 0*np.arange(num_elem*nsize_int*1.0) 
	f_bnd_rhs = 0*np.arange(num_elem*nsize_bndry*1.0) 
	f_int_rhs = 0*np.arange(num_elem*nsize_int*1.0) 
	f_p = 0*np.arange(num_elem*nsize_p*1.0) 
	f_p_rhs = 0*np.arange(num_elem*nsize_p*1.0) 
	bnd = 0*np.arange((num_elem*nsize_p*1.0 - num_elem) + nGlobBndDofs*1.0) 
	fh_bnd = 0*np.arange((num_elem*nsize_p*1.0 - num_elem) + nGlobBndDofs*1.0) # when no body forcing
	# Attention : this is hard-coded for the Channel example
#	dbnd = np.array([-2.5524342894999079e-15, 0.25, 0.2500000000000005, -7.6106763859384379e-16, 1.3910837861917781e-15, -1.23261736339018e-15, 1.8773266417199297e-15, -2.0894592931587442e-15, 3.0297952086575296e-15, -2.5524342894999079e-15, 2.3770867748962876e-15, -1.7612724401277993e-15, 1.1116566410042821e-15, 0.25, -2.5524342894999079e-15, 0.24999999999999914, -7.5756095477040149e-16, -1.7677313319122266e-15, -2.3149731017852412e-15, -3.9705237119397291e-15, -4.0243436284763416e-15, -3.9821690566654888e-15, -2.9506472566747979e-15, -3.2681356861139003e-15, -2.1123409087112298e-15, -2.8114906724929497e-15])
	#print("dbnd.shape ", dbnd.shape)
#	k = 0
#	bndcnt = 234
#	for j in range(0, dbnd.shape[0]):
#		bnd[int(BndCondCoeffsToGlobalCoeffsMap[bndcnt])] = dbnd[j]
	#	print("writing to bnd location: ", int(BndCondCoeffsToGlobalCoeffsMap[bndcnt]))
#		bndcnt += 1

	cnt = 0
	cnt1 = 0
	nvel = 2
	offset = 0
	for i in range (0, num_elem):
		for j in range(0, nvel):
			for k in range(0, nbmap):
				f_bnd_rhs[cnt+k] = forcing[j, offset + int(bmap[k])]
			for k in range(0, nimap):
				f_int_rhs[cnt1+k] = forcing[j, offset + int(imap[k])]
			cnt  += nbmap
			cnt1 += nimap
		offset += nbmap + nimap


	for curr_elem in range(0, num_elem):
		cnt = curr_elem*B_elem[curr_elem, :, :].shape[0]
		cnt1 = curr_elem*B_elem[curr_elem, :, :].shape[1]
		cnt2 = curr_elem*Dint_elem[curr_elem, :, :].shape[0]
		loc_BCinv = B_elem[curr_elem, :, :]
		loc_Cinv = D_elem[curr_elem, :, :]
		loc_Dint = Dint_elem[curr_elem, :, :]
		f_bnd_rhs[cnt : cnt + B_elem[curr_elem, :, :].shape[0]] = f_bnd_rhs[cnt : cnt + B_elem[curr_elem, :, :].shape[0]] - np.dot(loc_BCinv,  f_int_rhs[cnt1:cnt1 + B_elem[curr_elem, :, :].shape[1]])
		f_p_rhs[cnt2 : cnt2 + Dint_elem[curr_elem, :, :].shape[0]] = np.dot( loc_Dint , np.dot(loc_Cinv,  f_int_rhs[cnt1:cnt1 + B_elem[curr_elem, :, :].shape[1]]))

	offset = 0
	cnt = 0
	for curr_elem in range(0, num_elem):
		for j in range(0, nvel):
			for k in range(0, nbmap):
				fh_bnd[ int(LocGloMap[offset+j*nbmap+k]) ] += int(LocGloSign[offset+j*nbmap+k]) * f_bnd_rhs[cnt+k]
			cnt += nbmap
		offset += nsize_p + nvel*nbmap

	offset = 0
	cnt1 = 0
	for curr_elem in range(0, num_elem):
		for j in range(0, nsize_p):
			fh_bnd[ int(LocGloMap[offset + nvel*nbmap + j]) ] = f_p_rhs[cnt1 + j]
		cnt1 += nsize_p
		offset += nvel*nbmap + nsize_p



	bnd[0:init_bnd.shape[0]] = init_bnd
	# should double check the init_bnd
	print("np.linalg.norm(bnd) ", np.linalg.norm(bnd))
	print("np.linalg.norm(init_bnd) ", np.linalg.norm(init_bnd))
	print("bnd.shape ", bnd.shape)
	print("init_bnd.shape ", init_bnd.shape)

#	vdbc = ngc - 4512
	nLocBndDofs = num_elem*nsize_bndry + num_elem

#	nGlobHomBndDofs = 1312
	nGlobDofs = (num_elem*nsize_p*1.0 - num_elem) + nGlobBndDofs
	nDirBndDofs = NumDirBCs
	loc_dbc = 0*np.arange(nLocBndDofs*1.0)
	V_GlobBnd = bnd[0:nGlobBndDofs]



	for i in range(0,nLocBndDofs):
		loc_dbc[i] = int(LocGloBndSign[i]) * V_GlobBnd[int(LocGloBndMap[i])]

	my_V_GlobBnd_copy = V_GlobBnd.copy()
	loc_dbc_regenForcing = loc_dbc.copy()
	loc_dbc_regenForcing_no_pp = 0*np.arange(num_elem*nsize_bndry*1.0)
	m_p_regenForcing_no_pp = 0*np.arange(num_elem*nsize_p*1.0)



	for curr_elem in range(0, num_elem):
		cnt = curr_elem*nsize_bndry_p1
		cnt_no_pp = curr_elem*nsize_bndry
		cnt_m_p = curr_elem*nsize_p
		for i in range(0, nsize_bndry_p1):
			if (i < nsize_bndry):
				loc_dbc_regenForcing_no_pp[cnt_no_pp + i] = loc_dbc_regenForcing[cnt + i]
			else:
				m_p_regenForcing_no_pp[cnt_m_p] = loc_dbc_regenForcing[cnt + i]   # is zero


	for curr_elem in range(0, num_elem):
		cnt = curr_elem*nsize_bndry_p1
		loc_Ah = Ah_elem[curr_elem, :, :]
		loc_dbc[cnt:cnt+nsize_bndry_p1] = np.dot(loc_Ah,  loc_dbc[cnt:cnt+nsize_bndry_p1]) 

	for curr_elem in range(0, num_elem):
		cnt = curr_elem*nsize_bndry_p1
		cnt_no_pp = curr_elem*nsize_bndry
		cnt_m_p = curr_elem*nsize_p

	tmp_compare = np.dot(np.transpose(M_trafo), loc_dbc)



	V_GlobHomBndTmp = 0*np.arange(nGlobHomBndDofs*1.0)
	tmp = 0*np.arange(nGlobBndDofs*1.0)
	tmp_origForcing = 0*np.arange(nGlobBndDofs*1.0)
	for i in range(0,nLocBndDofs):
		tmp[int(LocGloBndMap[i])] += int(LocGloBndSign[i]) * loc_dbc[i]
		tmp_origForcing[int(LocGloBndMap[i])] += int(LocGloBndSign[i]) * loc_dbc_regenForcing[i]


	offset = nDirBndDofs

	V_GlobHomBndTmp = tmp[offset:nGlobBndDofs]

	V_GlobHomBndTmp = -V_GlobHomBndTmp + fh_bnd[NumDirBCs : nGlobHomBndDofs + NumDirBCs]
	my_sys_in = V_GlobHomBndTmp

#	my_sys_in = -V_GlobHomBndTmp

		
	# non-projected system
	my_Asolution = np.linalg.solve(my_Gmat, my_sys_in)
	
	
# add back initial conditions onto difference
	my_bnd_after = bnd.copy()
	
	my_bnd_after[offset:nGlobBndDofs] += my_Asolution
	
	V_GlobBnd = my_bnd_after[0:nGlobBndDofs]
	
	for i in range(0,nLocBndDofs):
		loc_dbc[i] = int(LocGloBndSign[i]) * V_GlobBnd[int(LocGloBndMap[i])]


	F_int = 0*np.arange((num_elem*nsize_p - num_elem)*1.0)
	for curr_elem in range(0, num_elem):
		cnt = curr_elem*nsize_bndry_p1
		loc_Ch = Ch_elem[curr_elem, :, :]
		F_int[curr_elem*nsize_p_m1:curr_elem*nsize_p_m1+nsize_p_m1] = fh_bnd[nGlobHomBndDofs + NumDirBCs + curr_elem*nsize_p_m1: nGlobHomBndDofs + NumDirBCs + curr_elem*nsize_p_m1+nsize_p_m1] - np.dot(loc_Ch,  loc_dbc[cnt:cnt+nsize_bndry_p1])


	for curr_elem in range(0, num_elem):
		cnt = curr_elem*nsize_p_m1
		loc_Dh = Dh_elem[curr_elem, :, :]
		F_int[cnt:cnt+nsize_p_m1] = np.dot(loc_Dh,  F_int[cnt:cnt+nsize_p_m1]) 


	my_bnd_after[nGlobBndDofs:my_bnd_after.shape[0]] = F_int

# first write f_bnd
	nz_loc = 1
#	nbnd = 48
#	nint = 121
	nvel = 2
	cnt = 0
	offset = 0
	no_loc_dbc = 0
	no_not_loc_dbc = 0
	elem_loc_dbc = set()
	elem_not_loc_dbc = set()
	#find = np.c_[0]
	for curr_elem in range(0, num_elem):
		for j in range(0, nvel):
			for k in range(0, nbnd):
				f_bnd[cnt+k] = int(LocGloSign[offset+j*nbnd+k]) * my_bnd_after[int(LocGloMap[offset+j*nbnd+k])]
				if (int(LocGloMap[offset+j*nbnd+k]) < NumDirBCs):
					no_loc_dbc = no_loc_dbc + 1
					elem_loc_dbc.update({cnt+k})
				else:
					no_not_loc_dbc = no_not_loc_dbc + 1
					elem_not_loc_dbc.update({cnt+k})
			cnt += nbnd
		offset += nvel*nbnd + nint


	
	offset = 0
	cnt = 0
	cnt1 = 0
	num_pressure_coeffs = num_elem*nsize_p
	p_coeffs = 0*np.arange(num_pressure_coeffs*1.0)

	for curr_elem in range(0, num_elem):
		for j in range(0, nint):
			f_p[cnt+j] = my_bnd_after[int(LocGloMap[offset+nvel*nbnd+j])]
			p_coeffs[cnt1+j] = f_p[cnt+j]
		cnt += nint
		offset += nvel*nbnd + nint 
		cnt1 += nsize_p
		

	print("f_p.shape ", f_p.shape)

	
#	nsize_bndry = 96
#	nsize_bndry_p1 = 97
#	nsize_int = 242
#	nsize_p = 121

	
	for curr_elem in range(0, num_elem):
		cnt_Dint = curr_elem*nsize_p
		cnt_C = curr_elem*nsize_bndry
		cnt_D = curr_elem*nsize_int
		loc_Dint = Dint_elem[curr_elem, :, :]
		loc_C = C_elem[curr_elem, :, :]
		loc_D = D_elem[curr_elem, :, :]
		f_int[curr_elem*nsize_int:curr_elem*nsize_int + nsize_int] = f_int_rhs[curr_elem*nsize_int:curr_elem*nsize_int + nsize_int] + np.dot(np.transpose(loc_Dint),  f_p[cnt_Dint:cnt_Dint+nsize_p]) - np.dot(np.transpose(loc_C),  f_bnd[cnt_C:cnt_C+nsize_bndry])
		f_int[curr_elem*nsize_int:curr_elem*nsize_int + nsize_int] = np.dot(loc_D, f_int[curr_elem*nsize_int:curr_elem*nsize_int + nsize_int]) 
	
	
	cnt = 0
	cnt1 = 0
	offset = 0
	nplanecoeffs = nlc # fields[m_velocity[0]]->GetNcoeffs() 
	fields = np.zeros( [nvel, nplanecoeffs] )

	for curr_elem in range(0, num_elem):
		for j in range(0, nvel):
			for k in range(0, nbnd):
				fields[j, offset + int(bmap[k])] = f_bnd[cnt + k]
			for k in range(0, nint):
				fields[j, offset + int(imap[k])] = f_int[cnt1+k]
			cnt += nbnd
			cnt1 += nint
		offset += nbnd + nint # is 169

	nz_loc = 1
#	nbnd = 48
#	nint = 121
	nvel = 2
	cnt = 0
	offset = 0
	velo0 = np.zeros( [ngc] )
	velo1 = np.zeros( [ngc] )
	velo0_ref = np.zeros( [ngc] )
	velo1_ref = np.zeros( [ngc] )

# double check the LocGloMapA by matrix replacement -- some dofs got reordered but this is probably fine
#	LocGloSignA = np.load('LocGloSignA.npy', 'r')
#	LocGloMapA = np.load('LocGloMapA.npy', 'r')
#	print("LocGloMapMatA.shape ", LocGloMapMatA.shape)
#	print("fields.shape ", fields.shape)

#	for i in range(0, nplanecoeffs):
#		velo0_ref[int(LocGloMapA[i])] = int(LocGloSignA[i]) * fields[0,i]
#		velo1_ref[int(LocGloMapA[i])] = int(LocGloSignA[i]) * fields[1,i]
	
	velo0 = np.dot(LocGloMapMatA, np.transpose(fields[0,:]))
	velo1 = np.dot(LocGloMapMatA, np.transpose(fields[1,:]))

	row1 = np.dot(sing_A, f_bnd ) - np.dot(np.transpose(sing_Dbnd), f_p) + np.dot(sing_B, f_int)
	row2 = -np.dot(sing_Dbnd, f_bnd) - np.dot(sing_Dint, f_int)
	row3 = np.dot(np.transpose(sing_Btilde), f_bnd) - np.dot(np.transpose(sing_Dint), f_p) + np.dot(sing_C, f_int)


	print(M_trafo_no_pp_incl_dbc.shape)
	print(row1.shape)
	print("np.linalg.norm(row1) ", np.linalg.norm(np.dot(row1, M_trafo_no_pp_incl_dbc)))
	print("np.linalg.norm(row1) ", np.linalg.norm(row1))
	print("np.linalg.norm(row2) ", np.linalg.norm(row2))
	print("np.linalg.norm(row3) ", np.linalg.norm(row3))

#	print( "np.linalg.norm(velo0) ", np.linalg.norm(velo0) )
#	print( "np.linalg.norm(velo1) ", np.linalg.norm(velo1) )
#	print( "np.linalg.norm(velo0_ref) ", np.linalg.norm(velo0_ref) )
#	print( "np.linalg.norm(velo1_ref) ", np.linalg.norm(velo1_ref) )
#	print( "np.linalg.norm(velo0_ref - velo0) ", np.linalg.norm(velo0_ref - velo0) )
#	print( "np.linalg.norm(velo1_ref - velo1) ", np.linalg.norm(velo1_ref - velo1) )


	return (velo0, velo1, f_bnd, f_p, f_int)
