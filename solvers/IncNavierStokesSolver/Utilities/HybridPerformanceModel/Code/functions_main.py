# Andrew Gloster
# Summer 2016
# Python Version of Communication Model Implementation
# Main Functions

#------------------------------------
# Import relevant modules
#------------------------------------

import os
from subprocess import Popen, PIPE

#------------------------------------
# Import relevant functions and classes
#------------------------------------

#------------------------------------
# New Function
#------------------------------------

def Filename_Generate(Mesh, Max_N_Z, Conditions_File):

    # Mesh File Location
    Mesh_File = 'Input/Mesh/' + Mesh

    # File containing the output of the maximum value of N_Z, for parsing CG iteration data.
    Input_Nektar_Max = 'Input/Serial_Input/' + Max_N_Z
 
    # File containing the conditions.
    Conditions = 'Input/Conditions/' + Conditions_File

    # Serial timing file location
    Loc_Serial_Timing_Files = 'Input/Serial_Input/'

    # Parallel timing file location
    Loc_Parallel_Timing_Files = 'Input/Parallel_Input/'

    # Hardware benchmarking files
    Benchmark_PBS = 'Input/Benchmark/Benchmark.pbs'
    MPI_Benchmark = 'Input/Benchmark/Benchmark.txt'
    Node_Map = 'Input/Benchmark/node.xml'
    
    return(Mesh_File, Input_Nektar_Max, Conditions, Loc_Serial_Timing_Files, Loc_Parallel_Timing_Files, Benchmark_PBS, MPI_Benchmark, Node_Map)

#------------------------------------
# New Function
#------------------------------------

# Parse the PBS script provided to find the number of nodes
def PBS_Job_Parse(Input_Filename):

    # Open the input file in .pbs format
    f = open('Input/' + Input_Filename, 'r')

    # Error check default to True
    Error = True
    Message = []

    # Iterate over the file looking for the number of nodes chosen by the user
    for line in f:
        a = line.split()
        for i in range(0, len(a)):
            b = a[i].split(':', 2)
            for j in range(0, len(b)):
                c = b[j].split('=', 1)
                if (c[0] == 'select'):
                    try:
                        Num_Node = int(c[1])
                    except:
                        Num_Node = 0
                        Error = False
                        Message.append('Unable to find number of nodes from ' + Input_Filename)

    # Close the file
    f.close()

    # Return the desired values and error information
    return(Num_Node, Error, Message)

#------------------------------------
# New Function
#------------------------------------

# Function to find the total number of elements for the serial calibration.
def Find_Nektar_Elements(Input_Filename):
    
    # Create a file to hold files required for the model to run
    output_path = 'Temporary_Files'
    if os.path.exists(output_path):
        cmd_string_clear = 'rm -r Temporary_Files/ \n'
        process = Popen([cmd_string_clear],shell=True, stdout=PIPE, stdin=PIPE)
        process.wait()
        os.mkdir(output_path)
    
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # Uncompress data if needs be using Nekmesh
    cmd_string_uncompress = 'Nekmesh ' + Input_Filename + ' Temporary_Files/Nektar_Serial_Mesh.xml:xml:uncompress' + " \n"

    # Run the uncompress command using Python'subprocess module
    process = Popen([cmd_string_uncompress],shell=True, stdout=PIPE, stdin=PIPE)
    process.wait()
    
    # Open the uncompressed mesh and count the total number of elements
    f = open('Temporary_Files/Nektar_Serial_Mesh.xml', "r")
    element_count = 0
    for line in f:
        a = line.split()

        # Record an element when strings match
        if(a[0] == '<Q'):
            element_count += 1
    f.close()

    return(element_count)

#------------------------------------
# New Function
#------------------------------------

# Find the files that make up the Serial Data being provided to the model
def Find_Nektar_Files(Input_Filename):

    # Now find the directory for the generated partition files
    current_directory = os.getcwd()

    directory = current_directory + '/' + Input_Filename

    # Holding list for file names
    Timing_Files = []

    # Find the file names in the given directory
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".txt"):
                Timing_Files.append(os.path.join(file))

    Nektar_Modes = []

    # Loop over the the timing file names to parse the nektar modes
    for i in range(0, len(Timing_Files)):
        a = Timing_Files[i].split('_', 1)
        b = a[1].split('.txt', 1)
        Nektar_Modes.append(int(b[0])) 

    # Sort the output by plane number
    (Nektar_Modes, Timing_Files) = zip(*sorted(zip(Nektar_Modes, Timing_Files)))

    # Return the file names
    return(Nektar_Modes, Timing_Files)

#------------------------------------
# New Function
#------------------------------------

# Pharse the Nektar file to find timing information.
def Parse_Nektar_Output(Input_Filename):

    # Open file to be pharsed
    f = open(Input_Filename, "r")
    
    # List for outputs
    times = []

    # Pharse desired data
    for line in f:
        a = line.split()
        for i in range(0, len(a)):
            if (a[i] == 'CPU'):
                b = a[i + 2].split('s', 1)
                times.append(float(b[0]))
   
    # Remove first two entries, not representative of the timings
    times.pop(0)
    times.pop(0)
    
    # Return the list of times
    return(times)

#------------------------------------
# New Function
#------------------------------------

# Parse the Nektar file to find CG information for final timestep
def Parse_Nektar_CG_Benchmark_Output(Input_Filename):

    # Open file to be pharsed
    f = open(Input_Filename, "r")

    # Pharse desired data
    for line in f:
        a = line.split()
        for i in range(0, len(a)):
            
            # Reset the data each time as we only want the last entry
            if (a[i] == 'Pressure'):
                Pressure = {}
                Velocity_1 = {}
                Velocity_2 = {}
                Velocity_3 = {}
                var = 1
                continue

            if (a[i] == 'Velocity'):
                var = 2
                continue

            if (a[i] == 'Plane'):
                plane = int(a[i + 1]) + 1
                plane = str(plane)
                continue
            
            # Append each value of CG to the dictionaries
            if (a[i] == 'CG'):
                if (var == 1):
                    if plane in Pressure.keys():
                        Pressure[plane].append(int(a[i + 4]))
                        continue
                    else:
                        Pressure[plane] = [int(a[i + 4])]

                if (var == 2):
                    if plane in Velocity_1.keys():
                        Velocity_1[plane].append(int(a[i + 4]))
                        var = var + 1
                        continue
                    else:
                        Velocity_1[plane] = [int(a[i + 4])]

                if (var == 3):
                    if plane in Velocity_2.keys():
                        Velocity_2[plane].append(int(a[i + 4]))
                        var = var + 1
                        continue
                    else:
                        Velocity_2[plane] = [int(a[i + 4])]

                if (var == 4):
                    if plane in Velocity_3.keys():
                        Velocity_3[plane].append(int(a[i + 4]))
                        continue
                    else:
                        Velocity_3[plane] = [int(a[i + 4])]

    # Return the dictionaries of CG iteration
    return(Pressure, Velocity_1, Velocity_2, Velocity_3)

#------------------------------------
# New Function
#------------------------------------

# Pharse the lstopo generated file to count socket and core related quantities
def Find_Hardware(Input_Filename):

    # Open the input file containing information generated by lstopo
    f = open(Input_Filename, 'r')

    # Declare initial counters
    count_cores = 0
    count_sockets = 0

    # Error check default to True
    Error = True
    Message = []

    # Count cores and socket by pharsing the lstopo file
    for line in f:
        a = line.split()
        for i in range(0, len(a)):
            if (a[i] == 'type="Core"'):
                count_cores += 1
            if (a[i] == 'type="Socket"'):
                count_sockets += 1

    # Close the file
    f.close()
    
    # PROC_TOT is one of the desired outputs
    Num_Core_Per_Node = count_cores

    # Error check the pharsing
    if(Num_Core_Per_Node == 0):
        Error = False
        Message.append('Unable to find any cores in ' + Input_Filename)

    if(count_sockets == 0):
        Error = False
        Message.append('Unable to find any sockets in ' + Input_Filename)
    
    # Find other desired quantities 
    Num_Core_Per_Socket = Num_Core_Per_Node/count_sockets
    Num_Sock_Per_Node = count_sockets

    # Return information along with error check
    return(Num_Core_Per_Node, Num_Core_Per_Socket, Num_Sock_Per_Node, Error, Message)

#------------------------------------
# New Function
#------------------------------------

# Pharse the condtions input file for P and N_Modes.
def Find_Conditions(Input_Filename):

	# Open the file to be pharsed
    f = open(Input_Filename, "r")

    # Error checking, default to True
    Error = True
    Message = []

    # Iterate over the file to find P and N_Modes
    for line in f:
        a = line.split()
        for i in range(0, len(a)):

            if (a[i] == 'HomModesZ'):
                try: 
                    N_Modes = int(a[i + 2])/2
                except:
                    Error = False

            b = a[i].split('=', 1)
            for j in range(0, len(b)):
                if(b[j] == 'NUMMODES'):
                    c = b[j + 1].split('"',2)
                    try:
                        P = int(c[1]) - 1
                    except:
                        Error = False

    # Close the file
    f.close()

    # Return the errors and desired values
    return(N_Modes, P, Error, Message)

#------------------------------------
# New Function
#------------------------------------

# Function to find PROC_BENCHMARK, number of processes used in Benchmarking
def PBS_Benchmark_Parse(Input_Filename):

    # Open the input file in .pbs format
    f = open(Input_Filename, 'r')

    # Error check default to True
    Error = True
    Message = []

    # Iterate over the file looking for the number of cores chosen by the user
    for line in f:
        a = line.split()
        for i in range(0, len(a)):
            b = a[i].split(':', 2)
            for j in range(0, len(b)):
                c = b[j].split('=', )
                if (c[0] == 'select'):
                    try:
                        Num_Node = int(c[1])
                    except:
                        Num_Node = 0
                        Error = False
                        Message.append('Unable to find number of nodes from ' + Input_Filename)
                if (c[0] == 'ncpus'):
                    try:
                        Num_Cores = int(c[1])
                    except:
                        Num_Cores = 0
                        Error = False
                        Message.append('Unable to find number of cores from ' + Input_Filename)
    
    # Calculate desired quantity
    PROC_BENCHMARK = Num_Cores * Num_Node

    # Close the file
    f.close()

    # Return the desired values and error information
    return(PROC_BENCHMARK, Error, Message)

#------------------------------------
# New Function
#------------------------------------

# Pharse the IBM/Intel MPI Benchmarking file to find bandwidths and latencies
def Parse_Benchmark(Input_Filename, PROC_BENCHMARK, Num_Core_Per_Socket, Num_Sock_Per_Node):

    # Calculate number of cores per node
    Num_Core_Per_Node = Num_Core_Per_Socket * Num_Sock_Per_Node

    # Usable grouping sizes to be checked
    Check = [4, 8, 16, 32]

    # Checks to be used in loop to decide group size choice
    Check_Node = False
    Check_Socket = False

    # Loop to find size choice, require the group that crosses socket and node divisions
    for i in range(0, len(Check)):
        if (Num_Core_Per_Socket % Check[i] != 0):
            Check_Socket = True

        if (Num_Core_Per_Node % Check[i] != 0):
            Check_Node = True
        
        # Break when found
        if (Check_Node is True and Check_Socket is True):
            PROC_PER_GROUP = Check[i]
            break
        
        # Reset the checks
        Check_Node = False
        Check_Socket = False
    
    # Number of groups 
    Num_Groups = PROC_BENCHMARK/PROC_PER_GROUP

    # Lists to hold groupings of groups
    concerned_node_groups = []
    concerned_socket_groups = []
    concerned_core_groups = []
    
    # List to hold groupings
    Groups = []

    # Counter used in iteration
    count = 0
    
    # Reset the checks, now used to confirm which group concerns which communication combination
    Check_Node = False
    Check_Socket = False

    # Loop over each group and categorise
    for i in range(0, Num_Groups):
        Groups.append([])

        for j in range(0, PROC_PER_GROUP):
            Groups[i].append(count)
            count += 1
    
        for j in range(1, PROC_PER_GROUP - 1):
            if (Groups[i][j] % Num_Core_Per_Node == 0):
                concerned_node_groups.append(i)
                Check_Node = True
                continue
                
            if (Groups[i][j] % Num_Core_Per_Socket == 0 and Check_Node is False):
                concerned_socket_groups.append(i)
                Check_Socket = True
                continue
                
        if (Check_Node is False and Check_Socket is False):
            concerned_core_groups.append(i)

        Check_Node = False
        Check_Socket = False
    
    # Open the file to be pharsed
    f = open(Input_Filename, "r")
    
    # List to hold pharsed data
    data = []

    # True/False statements and count set to default requirements
    count = -1
    Finder_0 = False
    Finder_1 = False
    Finder_2 = False
    Finder_3 = False

    # Pharse the data as required, the strings checked are those produced by the benchmarking tool.
    # Clunky but effective
    for line in f:
        pharsed = line.split()
        if (pharsed == ['#', 'Benchmarking', 'Multi-Exchange']):
            Finder_0 = True
            continue

        if (pharsed == ['#', '(', str(Num_Groups), 'groups', 'of', str(PROC_PER_GROUP), 'processes', 'each', 'running', 'simultaneous', ')'] and Finder_0 is True):
            Finder_1 = True
            continue

        if (pharsed == ['Group', '#bytes', '#repetitions', 't_min[usec]', 't_max[usec]', 't_avg[usec]', 'Mbytes/sec'] and Finder_1 is True):
            Finder_2 = True
            continue

        if (Finder_1 is True and Finder_2 is True):
            if (pharsed == []):
                count += 1
                data.append([[],[],[],[],[],[],[]])
                continue
            if (pharsed == ['#------------------------------------------------------------------------------------------']):
                break
    
            data[count][5].append(float(pharsed[5]))
            data[count][6].append(float(pharsed[6]))

    # Calculate Latencies using the groups
    count_lat_node = 0.0
    count_lat_socket = 0.0
    count_lat_core = 0.0
    
    for i in range(0, len(concerned_node_groups)):
        index = concerned_node_groups[i]
        count_lat_node += data[0][5][index]

    LAT_Node_To_Node = (count_lat_node/len(concerned_node_groups)) * 1e-06

    for i in range(0, len(concerned_socket_groups)):
        index = concerned_socket_groups[i]
        count_lat_socket += data[0][5][index]
    LAT_Socket_To_Socket = (count_lat_socket/len(concerned_socket_groups)) * 1e-06


    for i in range(0, len(concerned_core_groups)):
        index = concerned_core_groups[i]
        count_lat_core += data[0][5][index]

    LAT_Core_To_Core = (count_lat_core/len(concerned_core_groups)) * 1e-06


    # Calculate Bandwidth using the groups, memory size chosen by hand, adjust the minus to choose different size
    count_band_node = 0.0
    count_band_socket = 0.0
    count_band_core = 0.0
    
    for i in range(0, len(concerned_node_groups)):
        index = concerned_node_groups[i]
        count_band_node += data[count - 3][6][index]

    BW_Node_To_Node = (count_band_node/len(concerned_node_groups)) * 1e06

    for i in range(0, len(concerned_socket_groups)):
        index = concerned_socket_groups[i]
        count_band_socket += data[count - 3][6][index]

    BW_Socket_To_Socket = (count_band_socket/len(concerned_socket_groups)) * 1e06


    for i in range(0, len(concerned_core_groups)):
        index = concerned_core_groups[i]
        count_band_core += data[count - 3][6][index]

    BW_Core_To_Core = (count_band_core/len(concerned_core_groups)) * 1e06

    # Return the desired values
    return(BW_Node_To_Node, LAT_Node_To_Node, BW_Socket_To_Socket, LAT_Socket_To_Socket, BW_Core_To_Core, LAT_Core_To_Core)

#------------------------------------
# New Function
#------------------------------------

# Input the filename and number of processors you wish to partition by
def Partition(Input_Filename, PROC_XY):

    # Create a file to hold files required for the model to run
    output_path = 'Temporary_Files'
    if os.path.exists(output_path):
        cmd_string_clear = 'rm -r Temporary_Files/ \n'
        process = Popen([cmd_string_clear],shell=True, stdout=PIPE, stdin=PIPE)
        process.wait()
        os.mkdir(output_path)

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # Uncompress data if needs be using Nekmesh
    cmd_string_uncompress = 'Nekmesh ' + Input_Filename + ' Temporary_Files/uncompressed_mesh.xml:xml:uncompress' + " \n"

    # Run the uncompress command using Python'subprocess module
    process = Popen([cmd_string_uncompress],shell=True, stdout=PIPE, stdin=PIPE)
    process.wait()

    if(PROC_XY == 1):
        f = open('Temporary_Files/uncompressed_mesh.xml', "r")
        element_count = 0
        for line in f:
            a = line.split()

            # Record an element when strings match
            if(a[0] == '<Q'):
                element_count += 1
        f.close()
        return([0], [element_count])

    # Run partitioning part of IncNavierStokesSolver to find how METIS splits elements across processes
    cmd_string_partition = 'IncNavierStokesSolver Temporary_Files/uncompressed_mesh.xml --part-only ' + str(PROC_XY) + " \n"

    # Run the partition command using Python's subprocess module
    process = Popen([cmd_string_partition],shell=True, stdout=PIPE, stdin=PIPE)
    process.wait()

    # Now find the directory for the generated partition files
    current_directory = os.getcwd()
    mesh_part_folder = '/Temporary_Files/uncompressed_mesh_xml' 

    directory = current_directory + mesh_part_folder

    # Holding list for file names
    mesh_file = []

    # Find the file names in the given directory
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".xml"):
                mesh_file.append(os.path.join(file))

    # Initialise list for recording edges present on each core
    edges = []

    # Initialise list for recording elements present on each core
    Num_Elements = []

    # Stores the number of messages process i will send to process j
    dictionary_edges = []

    # Pharse all the edge data from the xml files and store them by each process in edges[]
    for i in range(0, len(mesh_file)):

        # Open each file
        f = open('Temporary_Files/uncompressed_mesh_xml/' + mesh_file[i], "r")

        # Append and update variables used for storage
        edges.append([])
        dictionary_edges.append({})
        element_count = 0

        # Iterate over the file, splitting each line into strings
        for line in f:
            a = line.split()

            # Record an edge when strings match
            if(a[0] == '<E'):
                b = a[1].split('"', 2)
                try:
                    edges[i].append(int(b[1]))

            # Python complains about indentation if we leave out an except:
            # So lets add 1 and 1 and call it Kwyjibo, the word made up by Bart Simpson
            # to win a game of scrable by using all his letters in one go
                except:
                    Kwyjibo = 1 + 1

            # Record an element when strings match
            if(a[0] == '<Q'):
                element_count += 1

        Num_Elements.append(element_count)

    # Initialise dictionary counters for cores
    for i in range(0, len(mesh_file)):
        for j in range(0, len(mesh_file)):
            if(j == i):
                dictionary_edges[i][str(j)] = 'Self'
                continue
            dictionary_edges[i][str(j)] = 0

    # Now compare edge lists between processes to find matches.
    # These correspond to neighbouring elements that must communicate.
    # We have +1 message recorded for a match between process i and k
    for i in range(0, len(mesh_file)):
        for k in range(0, len(mesh_file)):
            if(i == k):
                continue
            for j in range(0, len(edges[i])):
                a = edges[i][j]
                for n in range(0, len(edges[k])):
                    if(a == edges[k][n]):
                        dictionary_edges[i][str(k)] += 1

    # Put the counted edges into lists for general use later.
    Num_Element_Msg = []

    for i in range(0, len(mesh_file)):
        Num_Element_Msg.append([])
        for k in range(0, len(mesh_file)):
            Num_Element_Msg[i].append(dictionary_edges[i][str(k)])

    # Return the dictionary of values to be used later
    return (Num_Element_Msg, Num_Elements)

#------------------------------------
# New Function
#------------------------------------

# Find the possible combinations of PROC_Z and PROC_XY such that a cartesian grid is formed.
def Find_Topologies(PROC_TOT, Num_Modes):
    PROC_Z = []
    PROC_XY = []

    PROC_TOT_HALF = PROC_TOT/2

    for i in range(1, PROC_TOT_HALF + 1):
        if(PROC_TOT % i == 0 and Num_Modes % i == 0):
            PROC_Z.append(i)
            PROC_XY.append(PROC_TOT/i)

    if((Num_Modes % PROC_TOT) == 0):
        PROC_Z.append(PROC_TOT)
        PROC_XY.append(1)

    return (PROC_XY, PROC_Z)

#--=---------------------------------
# End of Functions
#------------------------------------