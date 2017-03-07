# Andrew Gloster
# Summer 2016
# Python Version of Communication Model Implementation
# Functions relating to full communication and computation model

#------------------------------------
# Import relevant modules
#------------------------------------

import matplotlib.pyplot as plt
import numpy as np

#------------------------------------
# Import relevant main functions
#------------------------------------

from functions_main import Find_Topologies
from functions_main import Partition
from functions_main import Find_Nektar_Files
from functions_main import Parse_Nektar_Output

#------------------------------------
# Import the required class
#------------------------------------

from class_topology import Topology

#------------------------------------
# New Function
#------------------------------------

# Function to run the full communication model
def Run_Parallel_Model(Parallelisation, Scheme, Mesh_File, Num_Modes, P, Num_Constants, Fit, BW_Node_To_Node, LAT_Node_To_Node, BW_Socket_To_Socket, 
	LAT_Socket_To_Socket, BW_Core_To_Core, LAT_Core_To_Core, Num_Core_Per_Socket, Num_Sock_Per_Node, PROC_TOT, Pressure, Velocity_1, Velocity_2,
	Velocity_3):

    # Find possible topologies when running in hybrid
    if (Parallelisation is 'Hybrid_Socket' or Parallelisation is 'Hybrid_Node'):
        (PROC_XY, PROC_Z) = Find_Topologies(PROC_TOT, Num_Modes)

    # Modal parallelisation values
    if (Parallelisation is 'Modal'):
        PROC_Z = [2, 4, 5, 8, 10, 20]
        PROC_XY = [1, 1, 1, 1, 1, 1]
    
    # Elemental parallelisation values
    if (Parallelisation is 'Elemental'):
        PROC_Z = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        PROC_XY = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]

#------------------------------------
# Topology Class Input and Results
#------------------------------------

    Pairwise = []
    Allreduce = [] 
    Alltoall = []

    Communication = []
    Serial = []

    Serial_Pairwise = []
    Serial_Allreduce = []
    Serial_Alltoall = []

    Total = []

    for i in range(0, len(PROC_Z)):

        # Set the communication count to 0
        Communication_Count = 0.0

        # Create class instance with basic input
        Simulation = Topology(PROC_Z[i], PROC_XY[i], Num_Core_Per_Socket, Num_Sock_Per_Node, Scheme)

        # Find the partition for the number of PROC_XY input
        (Num_Element_Msg, Num_Elements) = Partition(Mesh_File, PROC_XY[i])

        # Input the benchmarked communication data
        Simulation.Input_Communication(BW_Node_To_Node, LAT_Node_To_Node, BW_Socket_To_Socket, LAT_Socket_To_Socket, BW_Core_To_Core, LAT_Core_To_Core)

        # Distribute the modes and elements
        Simulation.Distribute_Modes(Num_Modes);
        Simulation.Distribute_Elements(Num_Element_Msg, Num_Elements);
    
        # Read in the CG iterations
        Simulation.CG_Iterations(Pressure, Velocity_1, Velocity_2, Velocity_3)

        # Input P to genererate the data sizes
        Simulation.Data_Size(P)

        # Input the fit from the serial model
        Simulation.Hardware_Constant(Num_Constants, Fit)

        # Count pairwise exchange
        Pairwise.append(Simulation.Communication_Pairwise_Exchange())
        Communication_Count += Pairwise[i]
    
        # Count Allreduce
        Allreduce.append(Simulation.Communication_Allreduce())
        Communication_Count += Allreduce[i]

        # Count Alltoall
        Alltoall.append(Simulation.Communication_Alltoall())
        Communication_Count += Alltoall[i]

        # Count Serial
        Serial.append(Simulation.Serial_Compute())

        # Count Serial with each type of communication
        Serial_Pairwise.append(Serial[i] + Pairwise[i])
        Serial_Allreduce.append(Serial[i] + Allreduce[i])
        Serial_Alltoall.append(Serial[i] + Alltoall[i])

        # Record the total communication
        Communication.append(Communication_Count)

        # Calculate the total length of a time step
        Total.append(Communication[i] + Serial[i])
    
    if (Parallelisation == 'Hybrid'):
        fig, ax = plt.subplots()
        ax.plot(PROC_Z, Total, label = 'Model')
        ax.set_xlabel('$ R_Z $')
        ax.set_ylabel('Timestep (s)')
        ax.set_title('Length of Single Timestep: Model')
        plt.legend(loc=1)
        fig.savefig("Output/Figures/Model_Hybrid.png")

    if (Parallelisation == 'Modal'):
        fig, ax = plt.subplots()
        ax.plot(PROC_Z, Total, label = 'Model')
        ax.set_xlabel('$ R_Z $')
        ax.set_ylabel('Timestep (s)')
        ax.set_title('Length of Single Timestep: Model')
        plt.legend(loc=1)
        fig.savefig("Output/Figures/Model_Modal.png")

    if (Parallelisation == 'Elemental'):
        fig, ax = plt.subplots()
        ax.plot(PROC_XY, Total, label = 'Model')
        ax.set_xlabel('$ R_{XY} $')
        ax.set_ylabel('Timestep (s)')
        ax.set_title('Length of Single Timestep: Model')
        plt.legend(loc=1)
        fig.savefig("Output/Figures/Model_Elemental.png")

    return(PROC_Z, PROC_XY, Total)

# Run a series of plots showing outputs of the model vs data
def Run_Parallel_Comparison(Loc_Parallel_Timing_Files, Parallelisation, PROC_Z, PROC_XY, Total):

    if (Parallelisation is 'Modal'):
        Timings_Fin = {}
        Data = []

        (Nektar_Modes_Fin, Timing_Files_Fin) = Find_Nektar_Files(Loc_Parallel_Timing_Files)

        for i in range(0, len(Timing_Files_Fin)):
            Data.append(np.mean(Parse_Nektar_Output(Loc_Parallel_Timing_Files + Timing_Files_Fin[i]))/10) 

        # Calculate the mean, standard deviation and variance of the difference between the data and the fitted model
        difference = []

        for i in range(0, len(Timing_Files_Fin)):
            difference.append(abs(Data[i] - Total[i]))
        
        mean_diff = np.mean(difference) 
        std_dev_diff = np.std(difference)
        var_diff = np.var(difference)
        
        # Print these results for the user to see
        print('The mean of the differences between the Data and the Model is ' + str(mean_diff))
        print('The standard deviation of the differences between the Data and the Model is ' + str(std_dev_diff))
        print('The variance of the differences between the Data and the Model is ' + str(var_diff))
            
        fig, ax = plt.subplots()
        ax.plot(PROC_Z, Data, label = 'Data')
        ax.plot(PROC_Z, Total, label = 'Model')
        ax.set_xlabel('$ R_Z $')
        ax.set_ylabel('Timestep (s)')
        ax.set_title('Length of Single Timestep: Model vs Data')
        plt.legend(loc=1)
        fig.savefig("Output/Figures/Mode_Full.png")

    if (Parallelisation is 'Elemental'):
        Timings_Fin = {}
        Data = []

        (Nektar_Modes_Fin, Timing_Files_Fin) = Find_Nektar_Files(Loc_Parallel_Timing_Files)

        for i in range(0, len(Timing_Files_Fin)):
            Data.append(np.mean(Parse_Nektar_Output(Loc_Parallel_Timing_Files + Timing_Files_Fin[i]))/10)

        # Calculate the mean, standard deviation and variance of the difference between the data and the fitted model
        difference = []

        for i in range(0, len(Timing_Files_Fin)):
            difference.append(abs(Data[i] - Total[i]))
        
        mean_diff = np.mean(difference) 
        std_dev_diff = np.std(difference)
        var_diff = np.var(difference)
        
        # Print these results for the user to see
        print('The mean of the differences between the Data and the Model is ' + str(mean_diff))
        print('The standard deviation of the differences between the Data and the Model is ' + str(std_dev_diff))
        print('The variance of the differences between the Data and the Model is ' + str(var_diff))

        fig, ax = plt.subplots()
        ax.plot(PROC_XY, Data, label = 'Data')
        ax.plot(PROC_XY, Total, label = 'Model')
        ax.set_xlabel('$ R_{XY} $')
        ax.set_ylabel('Timestep (s)')
        ax.set_title('Length of Single Timestep: Model vs Data')
        plt.legend(loc=1)
        fig.savefig("Output/Figures/Element_Full.png")

    if (Parallelisation is 'Hybrid_Socket' or Parallelisation is 'Hybrid_Node'):
        Data = []

        (Nektar_Modes_Fin, Timing_Files_Fin) = Find_Nektar_Files(Loc_Parallel_Timing_Files)

        for i in range(0, len(Timing_Files_Fin)):
            Data.append(np.mean(Parse_Nektar_Output(Loc_Parallel_Timing_Files + Timing_Files_Fin[i]))/10)

        # Calculate the mean, standard deviation and variance of the difference between the data and the fitted model
        difference = []

        for i in range(0, len(Timing_Files_Fin)):
            difference.append(abs(Data[i] - Total[i]))
        
        mean_diff = np.mean(difference) 
        std_dev_diff = np.std(difference)
        var_diff = np.var(difference)
        
        # Print these results for the user to see
        print('The mean of the differences between the Data and the Model is ' + str(mean_diff))
        print('The standard deviation of the differences between the Data and the Model is ' + str(std_dev_diff))
        print('The variance of the differences between the Data and the Model is ' + str(var_diff))

        fig, ax = plt.subplots()
        ax.plot(PROC_Z, Data, label = 'Data')
        ax.plot(PROC_Z, Total, label = 'Model')
        ax.set_xlabel('$ R_Z $')
        ax.set_ylabel('Timestep (s)')
        ax.set_title('Length of Single Timestep: Model vs Data')
        plt.legend(loc=1)
        fig.savefig("Output/Figures/Hybrid_Full.png")
