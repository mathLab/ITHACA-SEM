
import numpy as np

"""

runs using 
  1) python/3.3.2              2) numpy/1.11.1/python/3.3  		 3) scipy/0.17.1/python/3.3

module load python3/3.3.2 /scratch/mhess/py3p3/3.3

"""

def convert():

	print("Generating Parameter-independent terms.")



	f_IP = open('Nektar_data/IP.txt', 'r')
	f_IPp = open('Nektar_data/IPp.txt', 'r')
	f_cartmap0 = open('Nektar_data/cartmap0.txt', 'r')
	f_cartmap1 = open('Nektar_data/cartmap1.txt', 'r')
	IP = np.loadtxt(f_IP)
	IPp = np.loadtxt(f_IPp)
	cartmap0 = np.loadtxt(f_cartmap0)
	cartmap1 = np.loadtxt(f_cartmap1)
	np.save('ITHACA_SEM_data/IP', IP)
	np.save('ITHACA_SEM_data/IPp', IPp)
	np.save('ITHACA_SEM_data/cartmap0', cartmap0)
	np.save('ITHACA_SEM_data/cartmap1', cartmap1)

	nphys = cartmap0.shape[1]
	npc = IP.shape[0]
	num_elem = npc // nphys	

	for i in range(0,num_elem):
		HM_str = 'Nektar_data/HelmMat_' + str(i) + '.txt'
		f_HM = open(HM_str)
		HM = np.loadtxt(f_HM)
		HM_str = 'ITHACA_SEM_data/HM_' + str(i)
		np.save(HM_str, HM)
		

	f_bmap = open('Nektar_data/bmap.txt', 'r')
	f_imap = open('Nektar_data/imap.txt', 'r')
	bmap = np.loadtxt(f_bmap)
	np.save('ITHACA_SEM_data/bmap', bmap)
	imap = np.loadtxt(f_imap)
	np.save('ITHACA_SEM_data/imap', imap)
	
	f_bwdtrans = open('Nektar_data/bwdtrans.txt', 'r')
	f_LocGloBndMap = open('Nektar_data/LocGloBndMap.txt', 'r')
	f_LocGloBndSign = open('Nektar_data/LocGloBndSign.txt', 'r')
	f_BndCondCoeffsToGlobalCoeffsMap = open('Nektar_data/BndCondCoeffsToGlobalCoeffsMap.txt', 'r')
	f_LocGloMap = open('Nektar_data/LocGloMap.txt', 'r')
	f_LocGloSign = open('Nektar_data/LocGloSign.txt', 'r')
	f_glodofphys = open('Nektar_data/glo_dof_phys.txt', 'r')
	f_numDirBnd = open('Nektar_data/NumGlobalDirBndCoeffs.txt', 'r')
	f_LocGloMapMatA = open('Nektar_data/LocGloMapMatA.txt', 'r')
	f_forcing0 = open('Nektar_data/forcing0.txt', 'r')
	f_forcing1 = open('Nektar_data/forcing1.txt', 'r')
	f_bndcond_k0_i_0 = open('Nektar_data/bndcond_k0_i_0.txt', 'r')
	f_bndcond_k0_i_1 = open('Nektar_data/bndcond_k0_i_1.txt', 'r')
	f_bndcond_k0_i_2 = open('Nektar_data/bndcond_k0_i_2.txt', 'r')
	f_bndcond_k1_i_0 = open('Nektar_data/bndcond_k1_i_0.txt', 'r')
	f_bndcond_k1_i_1 = open('Nektar_data/bndcond_k1_i_1.txt', 'r')
	f_bndcond_k1_i_2 = open('Nektar_data/bndcond_k1_i_2.txt', 'r')

#	f_LocGloMapA = open('LocGloMapA.txt', 'r')
#	f_LocGloSignA = open('LocGloSignA.txt', 'r')
	
	bwdtrans = np.loadtxt(f_bwdtrans)
	LocGloBndMap = np.loadtxt(f_LocGloBndMap)
	LocGloBndSign = np.loadtxt(f_LocGloBndSign)
	BndCondCoeffsToGlobalCoeffsMap = np.loadtxt(f_BndCondCoeffsToGlobalCoeffsMap)
	LocGloMap = np.loadtxt(f_LocGloMap)
	LocGloSign = np.loadtxt(f_LocGloSign)
	glodofphys = np.loadtxt(f_glodofphys)
	numDirBnd = np.loadtxt(f_numDirBnd)
	LocGloMapMatA = np.loadtxt(f_LocGloMapMatA)
	forcing0 = np.loadtxt(f_forcing0)
	forcing1 = np.loadtxt(f_forcing1)
	bndcond_k0_i_0 = np.loadtxt(f_bndcond_k0_i_0)
	bndcond_k0_i_1 = np.loadtxt(f_bndcond_k0_i_1)
	bndcond_k0_i_2 = np.loadtxt(f_bndcond_k0_i_2)
	bndcond_k1_i_0 = np.loadtxt(f_bndcond_k1_i_0)
	bndcond_k1_i_1 = np.loadtxt(f_bndcond_k1_i_1)
	bndcond_k1_i_2 = np.loadtxt(f_bndcond_k1_i_2)



	np.save('ITHACA_SEM_data/bwdtrans', bwdtrans)
	np.save('ITHACA_SEM_data/LocGloBndMap', LocGloBndMap)
	np.save('ITHACA_SEM_data/LocGloBndSign', LocGloBndSign)
	np.save('ITHACA_SEM_data/BndCondCoeffsToGlobalCoeffsMap', BndCondCoeffsToGlobalCoeffsMap)
	np.save('ITHACA_SEM_data/LocGloMap', LocGloMap)
	np.save('ITHACA_SEM_data/LocGloSign', LocGloSign)
	np.save('ITHACA_SEM_data/glodofphys', glodofphys)
	np.save('ITHACA_SEM_data/numDirBnd', numDirBnd)
	np.save('ITHACA_SEM_data/LocGloMapMatA', LocGloMapMatA)
	np.save('ITHACA_SEM_data/forcing0', forcing0)
	np.save('ITHACA_SEM_data/forcing1', forcing1)
	np.save('ITHACA_SEM_data/bndcond_k0_i_0', bndcond_k0_i_0)
	np.save('ITHACA_SEM_data/bndcond_k0_i_1', bndcond_k0_i_1)
	np.save('ITHACA_SEM_data/bndcond_k0_i_2', bndcond_k0_i_2)
	np.save('ITHACA_SEM_data/bndcond_k1_i_0', bndcond_k1_i_0)
	np.save('ITHACA_SEM_data/bndcond_k1_i_1', bndcond_k1_i_1)
	np.save('ITHACA_SEM_data/bndcond_k1_i_2', bndcond_k1_i_2)

#	f_filetxt = open('Nektar_data/cavity_poi_Oseen_ROM.txt', 'r')
#	f_filetxt = open('Nektar_data/testsimu.txt', 'r')
#	np.save('ITHACA_SEM_data/testsimu', filetxt)

#	f_filetxt = open('Nektar_data/cavity_poi_Oseen_nD.txt', 'r')
#	filetxt = np.loadtxt(f_filetxt)
#	np.save('ITHACA_SEM_data/testsimu_nD', filetxt)

#	f_filetxt = open('Nektar_data/cavity_poi_Oseen_D.txt', 'r')
#	filetxt = np.loadtxt(f_filetxt)
#	np.save('ITHACA_SEM_data/testsimu_D', filetxt)

#	f_filetxt = open('Nektar_data/cav_tt_corr_D.txt', 'r')
#	filetxt = np.loadtxt(f_filetxt)
#	np.save('ITHACA_SEM_data/testsimu_time_D', filetxt)

#	f_filetxt = open('Nektar_data/cav_tt_corr_nD.txt', 'r')
#	filetxt = np.loadtxt(f_filetxt)
#s	np.save('ITHACA_SEM_data/testsimu_time_nD', filetxt)

#	f_physglodof = open('phys_glo_dof.txt', 'r')
#	physglodof = np.loadtxt(f_physglodof)
#	np.save('physglodof', physglodof)

#	f_f0 = open('forcing0.txt', 'r')
#	force0 = np.loadtxt(f_f0)
#	np.save('force0', force0)
#	f_f1 = open('forcing1.txt', 'r')
#	force1 = np.loadtxt(f_f1)
#	np.save('force1', force1)
#	LocGloMapA = np.loadtxt(f_LocGloMapA)
#	LocGloSignA = np.loadtxt(f_LocGloSignA)
#	np.save('LocGloSignA', LocGloSignA)
#	np.save('LocGloMapA', LocGloMapA)


