
# module load python3/3.3.2 /scratch/mhess/py3p3/3.3

import numpy as np
import time
#import gen_channel_narrow_func
#import gen_channel_narrowROM_func

#f_param_set = open('Re_CVT_Channel_recompute.txt', 'r')
#param_set = np.loadtxt(f_param_set)

num_samples = 40
fname = 'append_to_xml.txt'
out_file = open(fname, 'w')

for i in range(0, num_samples):
	out_file.write('<FUNCTION NAME="TestSnap' + str(i+1) + '">\n')
	out_file.write('<F FILE="link_to_data/xml_channel_narrowROM_' + str(i) + '.fld" />\n')
	out_file.write('</FUNCTION>\n')  
out_file.close()


#for i in range(0, num_samples):
#	fname = 'xml_channel_narrow_' + str(i) + '.xml'	
#	gen_channel_narrow_func.gen_channel_func(fname, param_set[i])
#	fname = 'xml_channel_narrowROM_' + str(i) + '.xml'	
#	advname = 'xml_channel_narrow_' + str(i)	
#	gen_channel_narrowROM_func.gen_channel_func(fname, param_set[i], advname)
#	fname = 'xml_channel_narrowROM_repeat_' + str(i) + '.xml'	
#	advname = 'xml_channel_narrowROM_' + str(i)	
#	gen_channel_narrowROM_func.gen_channel_func(fname, param_set[i], advname)

