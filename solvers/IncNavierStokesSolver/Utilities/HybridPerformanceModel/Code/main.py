# Andrew Gloster
# Summer 2016
# Python Version of Communication Model Implementation
# Main Program

#------------------------------------
# Import relevant modules
#------------------------------------

import os
from subprocess import Popen, PIPE

#------------------------------------
# Import relevant main functions
#------------------------------------

from functions_main import PBS_Job_Parse
from functions_main import PBS_Benchmark_Parse
from functions_main import Parse_Benchmark
from functions_main import Parse_Nektar_Output
from functions_main import Parse_Nektar_CG_Benchmark_Output
from functions_main import Find_Nektar_Elements
from functions_main import Find_Nektar_Files
from functions_main import Find_Hardware
from functions_main import Find_Conditions
from functions_main import Filename_Generate

#------------------------------------
# Import function for serial computation
#------------------------------------

from serial import Run_Serial_Fit

#------------------------------------
# Import function for parallel model
#------------------------------------

from parallel import Run_Parallel_Model
from parallel import Run_Parallel_Comparison

#------------------------------------
# Import relevant error checking functions
#------------------------------------

from error_main import Primary_Error_Out
from error_main import Error_Check_Proc_Input
from error_main import Error_Check_Comm_Input

#------------------------------------
# Main Program Begin
#------------------------------------

#------------------------------------
# User Model Choices
#------------------------------------

# Input the mesh file name, stored in Input/Mesh/
Mesh = 'cyl-small.xml'

# Maximum value of N_Z
Max_N_Z = 'output_72.txt'

# Input the conditions file name stored in Input/Conditions/
Conditions_File = 'conditions_80.xml'

# Choose the numerical scheme you're using, IterativeFull or IterativeStaticCond
Scheme = 'IterativeStaticCond'

# Select the values of N_Z you wish to calibrate the model with
Consider_Modes = [80]

# Choose number of constant to fit to the serial model model, 1 or 2
Num_Constants = 2

# Do you wish to run the clearcomparison between the calibrated serial model prediction and the collected data.
# Set to false if your serial data was collected to a different core than the one you are calibrating for.
Compare_Serial = False

# Do you wish to run the parallel portion of the model
Parallel = True

# Choose the Parallelisation method for the model to adopt, Modal, Elemental, Hybrid_Socket or Hybrid_Node 
Parallelisation = 'Hybrid_Socket'

# Do you wish to run the comparison between the full model prediction and the collected data
Compare_Parallel = True

# Input number of nodes you are using (Can replace this with PBS_Job_Parse)
Num_Node = 1

#------------------------------------
# Generate Input Filename Locations
#------------------------------------

(Mesh_File, Input_Nektar_Max, Conditions, Loc_Serial_Timing_Files, Loc_Parallel_Timing_Files, Benchmark_PBS, MPI_Benchmark, Node_Map) = Filename_Generate(Mesh, Max_N_Z, Conditions_File)

#------------------------------------
# Generate Output Locations
#------------------------------------

output_path = 'Output/Figures/'
if os.path.exists(output_path):
    cmd_string_clear = 'rm -r Output/Figures/ \n'
    process = Popen([cmd_string_clear],shell=True, stdout=PIPE, stdin=PIPE)
    process.wait()
    os.mkdir(output_path)
    
if not os.path.exists(output_path):
    os.mkdir(output_path)
    
#------------------------------------
# Parse Max Modes Nektar Input
#------------------------------------

(Pressure, Velocity_1, Velocity_2, Velocity_3) = Parse_Nektar_CG_Benchmark_Output(Input_Nektar_Max)

#------------------------------------
# Pharse Serial Nektar Inputs
#------------------------------------

Nektar_Serial_Elements = Find_Nektar_Elements(Mesh_File)

(Nektar_Modes, Timing_Files) = Find_Nektar_Files(Loc_Serial_Timing_Files)

Timings = {}
for i in range(0, len(Timing_Files)):
    Timings[str(Nektar_Modes[i])] = Parse_Nektar_Output(Loc_Serial_Timing_Files + Timing_Files[i])

#------------------------------------
# Pharse Condition Inputs
#------------------------------------

(Num_Modes, P, Error, Message) = Find_Conditions(Conditions)

# Print error to file if needed
if (Error is False):
    Primary_Error_Out(Message)
    quit()

#------------------------------------
# Fit Serial Model
#------------------------------------

# Fit the model to the supplied data and return the FLOPs
Fit = Run_Serial_Fit(Compare_Serial, Consider_Modes, Num_Constants, P, Nektar_Serial_Elements, Nektar_Modes, Timings, Pressure, Velocity_1, Velocity_2, Velocity_3, 'IterativeStaticCond')

# Stop here if we do not wish to run the parallel model.
if Parallel is False:
    quit()

#------------------------------------
# Pharse Hardware Inputs
#------------------------------------

(Num_Core_Per_Node, Num_Core_Per_Socket, Num_Sock_Per_Node, Error, Message) = Find_Hardware(Node_Map)

# Print error to file if needed
if (Error is False):
    Primary_Error_Out(Message)
    quit()

# Calculate PROC_TOT using outputs of the PBS_Pharse and Find_Hardware functions
PROC_TOT = Num_Core_Per_Node * Num_Node

# #------------------------------------
# # Pharse Communication Inputs
# #------------------------------------

(PROC_BENCHMARK, Error, Message) = PBS_Benchmark_Parse(Benchmark_PBS)

# Print error to file if needed
if (Error is False):
    Primary_Error_Out(Message)
    quit()

(BW_Node_To_Node, LAT_Node_To_Node, BW_Socket_To_Socket, LAT_Socket_To_Socket, BW_Core_To_Core, LAT_Core_To_Core) = Parse_Benchmark(MPI_Benchmark, PROC_BENCHMARK, Num_Core_Per_Socket, Num_Sock_Per_Node)

print(BW_Node_To_Node, LAT_Node_To_Node, BW_Socket_To_Socket, LAT_Socket_To_Socket, BW_Core_To_Core, LAT_Core_To_Core)

#------------------------------------
# Error check all inputs
#------------------------------------

# Error check the inputs
(Error, Message) = Error_Check_Proc_Input(PROC_TOT, Num_Core_Per_Socket, Num_Sock_Per_Node, Num_Modes)

# Print error to file if needed
if (Error is False):
    Primary_Error_Out(Message)
    quit()

(Error, Message) = Error_Check_Comm_Input(BW_Node_To_Node, LAT_Node_To_Node, BW_Socket_To_Socket, LAT_Socket_To_Socket, BW_Core_To_Core, LAT_Core_To_Core)

# Print error to file if needed
if (Error is False):
    Primary_Error_Out(Message)
    quit()

(PROC_Z, PROC_XY, Total) = Run_Parallel_Model(Parallelisation, Scheme, Mesh_File, Num_Modes, P, Num_Constants, Fit, BW_Node_To_Node, LAT_Node_To_Node, BW_Socket_To_Socket, 
    LAT_Socket_To_Socket, BW_Core_To_Core, LAT_Core_To_Core, Num_Core_Per_Socket, Num_Sock_Per_Node, PROC_TOT, Pressure, Velocity_1, Velocity_2, Velocity_3)

if Compare_Parallel is False:
    quit()

Run_Parallel_Comparison(Loc_Parallel_Timing_Files, Parallelisation, PROC_Z, PROC_XY, Total)

#------------------------------------
# Main Program End
#------------------------------------



