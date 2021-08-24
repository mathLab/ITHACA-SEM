# Andrew Gloster
# Summer 2016
# Python Version of Communication Model Implementation
# Error Check Functions

#------------------------------------
# Import relevant modules
#------------------------------------

import os

#------------------------------------
# Import relevant functions and classes
#------------------------------------


#------------------------------------
# New Function
#------------------------------------

def Primary_Error_Out(Message):

    # Check for file path and create if necessary
    error_path = 'Error'
    if not os.path.exists(error_path):
        os.mkdir(error_path)

    error_txt = open('Error/error_proc_input.txt', 'w')

    for i in range(0, len(Message)):
        error_txt.write(Message[i])
        error_txt.write("\n")

    error_txt.close()
       
#------------------------------------
# New Function
#------------------------------------

# Check the Processor, Element and Mode inputs for Errors
# Run this whenever using find_topologies from function_main.py
def Error_Check_Proc_Input(PROC_TOT, Num_Core_Per_Socket, Num_Sock_Per_Node, N_Modes):

	# Default no errors
	check = True

	# Check for file path and create if necessary
	error_path = 'Error'
	if not os.path.exists(error_path):
		os.mkdir(error_path)

	# Open file to write errors to.
	error_txt = open('Error/error_proc_input.txt', 'w')

	# PROC_TOT must be an int
	if(type(PROC_TOT) is not int):
		error_txt.write("Please input an integer for the number of Processors")
		error_txt.write('\n')
		check = False

	# PROC_TOT must be a positive number
	if(PROC_TOT <= 0):
		error_txt.write("Please give a positive number of Processors")
		error_txt.write('\n')
		check = False

	# PROC_TOT does not allow for an odd number of cores to be used.
	if(PROC_TOT % 2 != 0):
		error_txt.write("Odd Number of Processors Provided, Povide an Even Number")
		error_txt.write('\n')
		check = False

	# Num_Core_Per_Socket must be an int
	if(type(Num_Core_Per_Socket) is not int):
		error_txt.write("Please input an integer for the number of Cores Per Socket")
		error_txt.write('\n')
		check = False

	# Num_Core_Per_Socket must be a positive number
	if(Num_Core_Per_Socket <= 0):
		error_txt.write("Please give a positive number of Cores Per Socket")
		error_txt.write('\n')
		check = False

	# Num_Core_Per_Socket must be an int
	if(type(Num_Sock_Per_Node) is not int):
		error_txt.write("Please input an integer for the number of Sockets Per Node")
		error_txt.write('\n')
		check = False

	# Num_Core_Per_Socket must be a positive number
	if(Num_Sock_Per_Node <= 0):
		error_txt.write("Please give a positive number of Sockets Per Node")
		error_txt.write('\n')
		check = False

	# Num_Core_Per_Socket must be an int
	if(type(N_Modes) is not int):
		error_txt.write("Please input an integer for the number of Modes")
		error_txt.write('\n')
		check = False

	# Num_Core_Per_Socket must be a positive number
	if(N_Modes <= 0):
		error_txt.write("Please give a positive number of Modes")
		error_txt.write('\n')
		check = False

	# Close the file and return result
	error_txt.close()

	if(check is True):
		return(check, 'No errors in Processor, Element and Mode inputs')
	if(check is False):
		return(check, 'Errors found in Processor, Element and Mode inputs - Please see Error/error_proc_input.txt')	 

#------------------------------------
# New Function
#------------------------------------

# Check that the bandwidth and latency input are valid
def Error_Check_Comm_Input(BW_Node_To_Node, LAT_Node_To_Node, BW_Socket_To_Socket, LAT_Socket_To_Socket, BW_Core_To_Core, LAT_Core_To_Core):

	# Default no errors
	check = True

	# Check for file path and create if necessary
	error_path = 'Error'
	if not os.path.exists(error_path):
		os.mkdir(error_path)

	# Open file to write errors to.
	error_txt = open('Error/error_comm_input.txt', 'w')

	# BW_Node_To_Node must be an float
	if(type(BW_Node_To_Node) is not float):
		error_txt.write("Please input an float for Bandwith Node To Node")
		error_txt.write('\n')
		check = False

	# BW_Node_To_Node must be a positive number
	if(BW_Node_To_Node <= 0):
		error_txt.write("Please give a positive number for Bandwith Node To Node")
		error_txt.write('\n')
		check = False

	# LAT_Node_To_Node must be an float
	if(type(LAT_Node_To_Node) is not float):
		error_txt.write("Please input an float for Latency Node To Node")
		error_txt.write('\n')
		check = False

	# LAT_Node_To_Node must be a positive number
	if(LAT_Node_To_Node <= 0):
		error_txt.write("Please give a positive number for Latency Node To Node")
		error_txt.write('\n')
		check = False

	# BW_Socket_To_Socket must be an float
	if(type(BW_Socket_To_Socket) is not float):
		error_txt.write("Please input an float for Bandwith Socket To Socket")
		error_txt.write('\n')
		check = False

	# BW_Socket_To_Socket must be a positive number
	if(BW_Socket_To_Socket <= 0):
		error_txt.write("Please give a positive number for Bandwith Socket To Socket")
		error_txt.write('\n')
		check = False

	# LAT_Socket_To_Socket must be an float
	if(type(LAT_Socket_To_Socket) is not float):
		error_txt.write("Please input an float for Latency Socket To Socket")
		error_txt.write('\n')
		check = False

	# LAT_Socket_To_Socket must be a positive number
	if(LAT_Socket_To_Socket <= 0):
		error_txt.write("Please give a positive number for Latency Socket To Socket")
		error_txt.write('\n')
		check = False

	# BW_Core_To_Core must be an float
	if(type(BW_Core_To_Core) is not float):
		error_txt.write("Please input an float for Bandwith Core To Core")
		error_txt.write('\n')
		check = False

	# BW_Socket_To_Socket must be a positive number
	if(BW_Core_To_Core <= 0):
		error_txt.write("Please give a positive number for Bandwith Core To Core")
		error_txt.write('\n')
		check = False

	# LAT_Core_To_Core must be an float
	if(type(LAT_Core_To_Core) is not float):
		error_txt.write("Please input an float for Latency Core To Core")
		error_txt.write('\n')
		check = False

	# LAT_Core_To_Core must be a positive number
	if(LAT_Core_To_Core <= 0):
		error_txt.write("Please give a positive number for Latency Core To Core")
		error_txt.write('\n')
		check = False

	if(check is True):
		return(check, 'No errors in Processor, Element and Mode inputs')
	if(check is False):
		return(check, 'Errors found in Communication inputs - Please see Error/error_comm_input.txt')	 




