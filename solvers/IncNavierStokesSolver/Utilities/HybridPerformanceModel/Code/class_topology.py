# Andrew Gloster
# Summer 2016
# Python Version of Communication Model Implementation
# Topology Class

#------------------------------------
# Import relevant modules
#------------------------------------

import os

#------------------------------------
# Import relevant functions
#------------------------------------

from serial import Serial_Computation

#------------------------------------
# Class Definitions
#------------------------------------

class Topology:

#------------------------------------
# New Function
#------------------------------------

	# Initialisation of the Class
    def __init__(self, PROC_Z, PROC_XY, Num_Core_Per_Socket, Num_Sock_Per_Node, Scheme):

    	# Store read in data
        self.PROC_Z = PROC_Z
        self.PROC_XY = PROC_XY
        self.Num_Core_Per_Socket = Num_Core_Per_Socket
        self.Num_Sock_Per_Node = Num_Sock_Per_Node
        self.Scheme = Scheme

       	# Calculate PROC_TOT
       	self.PROC_TOT = self.PROC_XY * self.PROC_Z

        # Initialise storage lists
        self.Core = []
        self.Socket = []
        self.Node = []

        # Begin counters used in for loops
        count_core = 0
        count_socket_old = 0
        count_socket_new = 0
        count_node = 0

        # Iterate over cores assigning core, socket and node number
        for i in range(0, self.PROC_Z):
            for j in range(0, self.PROC_XY):
                self.Core.append(count_core)
                self.Socket.append(count_socket_new)
                self.Node.append(count_node)

                count_core += 1

                # Use modulo to find when to iterate socker number based on bumber of cores assigned
                if(count_core % Num_Core_Per_Socket == 0):
                    count_socket_new = count_socket_new + 1
            
                # Use modulo to find when to iterate node number based on bumber of sockets assigned
                if(count_socket_old != count_socket_new):
                    if(count_socket_new % Num_Sock_Per_Node == 0):
                        count_node = count_node + 1;
                
            	# Update socket counter
                count_socket_old = count_socket_new;

        # Flags to track operations that have occured, default to False
        self.Elements_Disributed = False
        self.Modes_Disributed = False
        self.Data_Size_Input = False
        self.Data_Communication_Input = False


#------------------------------------
# New Function
#------------------------------------

    # Print out hardware layout to user when desired, also used for error checking
    def Print_Hardware(self):

    	# Check that there is a file to put outputs into, if not make one
        output_path = 'Output'
        if not os.path.exists(output_path):
            os.mkdir(output_path)

    	# Open new file and write out configurations
        f = open('Output/Hardware_Topology.txt', "w")

        f.write("Core Distribution")
        f.write("\n")

        for i in range(0, self.PROC_Z):
            for j in range(0, self.PROC_XY):
                f.write(repr(self.Core[i*self.PROC_XY + j]).rjust(5))
            f.write("\n")
        

        f.write("\n")


        f.write("Socket Distribution")
        f.write("\n")
        for i in range(0, self.PROC_Z):
            for j in range(0, self.PROC_XY):
                f.write(repr(self.Socket[i*self.PROC_XY + j]).rjust(5))
            f.write("\n")


        f.write("\n")


        f.write("Node Distribution")
        f.write("\n")
        for i in range(0, self.PROC_Z):
            for j in range(0, self.PROC_XY):
                f.write(repr(self.Node[i*self.PROC_XY + j]).rjust(5))
            f.write("\n")

        # Close the file
        f.close()

#------------------------------------
# New Function
#------------------------------------
	
	# Distribute the elements as specified over the computational domain
    def Distribute_Elements(self, Num_Element_Msg, Num_Elements):
        self.Num_Elements = []
    	for i in range(0, self.PROC_Z):
            self.Num_Elements.append([])
            for j in range(0, self.PROC_XY):
                self.Num_Elements[i].append(Num_Elements[j])

        self.Num_Element_Msg = Num_Element_Msg
    
    	self.Elements_Disributed = True

#------------------------------------
# New Function
#------------------------------------

	# Distribute the modes as specified over the computational domain
    def Distribute_Modes(self, Num_Modes):
        self.Num_Modes = Num_Modes
        self.N_Z = Num_Modes * 2
        self.Split =  self.Num_Modes / self.PROC_Z
        self.Modes = []
        self.Planes = []
        self.Plane_Num = []
        for i in range(0, self.PROC_TOT):
            self.Modes.append(self.Split)
            self.Planes.append(self.Split * 2)
        
        for i in range(0, self.PROC_Z):
            self.Plane_Num.append([])
            for j in range(0, self.Split * 2):
                self.Plane_Num[i].append(i * (self.Split * 2) + j)

        self.Modes_Disributed = True

#------------------------------------
# New Function
#------------------------------------

    # Read in the fitted constants from the serial model
    def Hardware_Constant(self, Num_Constants, constants):
        self.Num_Constants = Num_Constants
        self.constants = constants


#------------------------------------
# New Function
#------------------------------------
	
	# Print the element distribution to file, also used for error checking
    def Print_Elements(self, Index):

    	if (self.Elements_Disributed is True):

    	    # Check that there is a file to put outputs into, if not make one
            output_path = 'Output'
            if not os.path.exists(output_path):
                os.mkdir(output_path)

    	    # Open new file and write out configurations
            f = open('Output/Element_Distribution_' + str(Index) + '.txt', "w")

            f.write("Core Distribution")
            f.write("\n")

            for i in range(0, self.PROC_Z):
                for j in range(0, self.PROC_XY):
                    f.write(repr(self.Core[i * self.PROC_XY + j]).rjust(5))
                f.write("\n")

            f.write("\n")

            f.write("Element Distribution")
            f.write("\n")

            for i in range(0, self.PROC_Z):
                for j in range(0, self.PROC_XY):
                    f.write(repr(self.Num_Elements[i][j]).rjust(5))
                f.write("\n")

        else:
        	# Check that there is a file to put outputs into, if not make one
            output_path = 'Output'
            if not os.path.exists(output_path):
                os.mkdir(output_path)

    	    # Open new file and write out configurations
            f = open('Output/Element_Distribution_' + str(Index) + '.txt', "w")

            f.write('Elements have no yet been distributed by the program.')

            f.close()


#------------------------------------
# New Function
#------------------------------------

	# Print the mode distribution to file, also used for error checking
    def Print_Modes(self, Index):

    	if (self.Modes_Disributed is True):

    	    # Check that there is a file to put outputs into, if not make one
            output_path = 'Output'
            if not os.path.exists(output_path):
                os.mkdir(output_path)

    	    # Open new file and write out configurations
            f = open('Output/Mode_Distribution.txt', "w")

            f.write("Core Distribution")
            f.write("\n")

            for i in range(0, self.PROC_Z):
                for j in range(0, self.PROC_XY):
                    f.write(repr(self.Core[i*self.PROC_XY + j]).rjust(5))
                f.write("\n")

            f.write("\n")

            f.write("Mode Distribution")
            f.write("\n")
            for i in range(0, self.PROC_Z):
                for j in range(0, self.PROC_XY):
                    f.write(repr(self.Modes[i*self.PROC_XY + j]).rjust(5))
                f.write("\n")

        else:
        	# Check that there is a file to put outputs into, if not make one
            output_path = 'Output'
            if not os.path.exists(output_path):
                os.mkdir(output_path)

    	    # Open new file and write out configurations
            f = open('Output/Mode_Distribution_' + str(Index) + '.txt', "w")

            f.write('Modes have no yet been distributed by the program.')

            f.close()
#------------------------------------
# New Function
#------------------------------------
    
    # Read in the CG iterations for the planes
    def CG_Iterations(self, Pressure, Velocity_1, Velocity_2, Velocity_3):
        self.Pressure = Pressure
        self.Velocity_1 = Velocity_1
        self.Velocity_2 = Velocity_2
        self.Velocity_3 = Velocity_3

#------------------------------------
# New Function
#------------------------------------

    # Calulate data size in bytes for each element and mode, elements and modes assumed to have homogeneous P
    # 8 bytes per double
    def Data_Size(self, P):
        self.P = P
    	self.Data_Element = (P + 1) * 8
    	self.Data_Mode = 8 * ((P+1) ** 2)

    	# Update tracker flag
        self.Data_Size_Input = True

#------------------------------------
# New Function
#------------------------------------

    # Communication inputs as parsed from MPI Benchmarking
    def Input_Communication(self, BW_Node_To_Node, LAT_Node_To_Node, BW_Socket_To_Socket, LAT_Socket_To_Socket, BW_Core_To_Core, LAT_Core_To_Core):
    	self.BW_Node_To_Node = BW_Node_To_Node
    	self.BW_Socket_To_Socket = BW_Socket_To_Socket
    	self.BW_Core_To_Core = BW_Core_To_Core
    	self.LAT_Node_To_Node = LAT_Node_To_Node
    	self.LAT_Socket_To_Socket = LAT_Socket_To_Socket
    	self.LAT_Core_To_Core = LAT_Core_To_Core

		# Update tracker flag
        self.Data_Communication_Input = True

#------------------------------------
# New Function
#------------------------------------

    # Check the relation between any two cores input to the function, whether they are on the same node/socket etc.
    def Check_Neighbour(self, core_1, core_2):
        if(core_1 == core_2):
            return 0
    
        if(self.Socket[core_1] == self.Socket[core_2]):
            return 1
    
        if(self.Node[core_1] == self.Node[core_2]):
            return 2;
    
        else:
            return 3;

#------------------------------------
# New Function
#------------------------------------

    # Model of the exchange of element data pairwise, using Alltoall style
    def Communication_Pairwise_Exchange(self):
        comm_max = 0.0

        for i in range(0, self.PROC_Z):

            current_row = i

            CG_Iter = 0 

            for j in range(0, len(self.Plane_Num[current_row])):
                try:
                    Plane = self.Plane_Num[current_row][j] + 1
                    CG_Iter += self.Pressure[str(Plane)][0]
                    CG_Iter += self.Velocity_1[str(Plane)][0]
                    CG_Iter += self.Velocity_2[str(Plane)][0]
                    CG_Iter += self.Velocity_3[str(Plane)][0]
                except:
                    continue

            for j in range(0, self.PROC_XY):
                comm_count = 0.0
                current = i * self.PROC_XY + j

                for k in range(0, self.PROC_XY):
                    send_loc = i * self.PROC_XY + k

                    relation = self.Check_Neighbour(current, send_loc)

                    if (relation == 1):
                        comm_count += self.LAT_Core_To_Core + (self.Num_Element_Msg[j][k] * self.Data_Element)/(self.BW_Core_To_Core)

                    if (relation == 2):
                        comm_count += self.LAT_Socket_To_Socket + (self.Num_Element_Msg[j][k] * self.Data_Element)/(self.BW_Socket_To_Socket)

                    if (relation == 3):
                        comm_count += self.LAT_Node_To_Node + (self.Num_Element_Msg[j][k] * self.Data_Element)/(self.BW_Node_To_Node)
                
                comm_count = comm_count * CG_Iter

                if (comm_count >= comm_max):
                    comm_max = comm_count

        return(comm_max)

#------------------------------------
# New Function
#------------------------------------
    
    # Model to describe MPI_AlltoALL for mode communcation
    def Communication_Allreduce(self):

        # We count the number of 
        comm_max = 0.0

        # Set max relation size tracking to False
        relation_1 = False
        relation_2 = False
        relation_3 = False

        for i in range(0, self.PROC_Z):
            current = self.Core[i * self.PROC_XY]

            current_row = i

            CG_Iter = 0 
            for j in range(0, len(self.Plane_Num[current_row])):
                try:
                    Plane = self.Plane_Num[current_row][j] + 1
                    CG_Iter += self.Pressure[str(Plane)][0]
                    CG_Iter += self.Velocity_1[str(Plane)][0]
                    CG_Iter += self.Velocity_2[str(Plane)][0]
                    CG_Iter += self.Velocity_3[str(Plane)][0]
                except:
                    continue
            
            comm_count = 0.0

            for j in range(0, self.PROC_XY - 1):
                current_right = current + 1

                relation = self.Check_Neighbour(current, current_right)

                if (relation == 1):
                    comm_count += self.LAT_Core_To_Core + (3 * 8)/(self.BW_Core_To_Core)
                    relation_1 = True

                if (relation == 2):
                    comm_count += self.LAT_Socket_To_Socket + (3 * 8)/(self.BW_Socket_To_Socket)
                    relation_2 = True

                if (relation == 3):
                    comm_count += self.LAT_Node_To_Node + (3 * 8)/(self.BW_Node_To_Node)
                    relation_3 = True

                current += 1
            
            comm_count = comm_count * CG_Iter
            if (comm_count >= comm_max):
                comm_max = comm_count
                CG_MAX = CG_Iter
        
        if (relation_1 is True and relation_2 is False and relation_3 is False):
            comm_max += (self.LAT_Core_To_Core + (3 * 8)/(self.BW_Core_To_Core)) * CG_MAX

        if (relation_1 is True and relation_2 is True and relation_3 is False):
            comm_max += (self.LAT_Socket_To_Socket + (3 * 8)/(self.BW_Socket_To_Socket)) * CG_MAX

        if (relation_1 is True and relation_2 is True and relation_3 is True):
            comm_max += (self.LAT_Node_To_Node + (3 * 8)/(self.BW_Node_To_Node)) * CG_MAX

        return(comm_max)

#------------------------------------
# New Function
#------------------------------------
    
    # Model to describe MPI_AlltoALL for mode communcation
    def Communication_Alltoall(self):
        
        # Find how long each core requires to send its information to each level, core, socker or node.
        MSG_Node_Mode = []
        MSG_Socket_Mode = []
        MSG_Core_Mode = []
        
        # Divide by PROC_Z as Alltoall involves this data size
        for i in range(0, self.PROC_Z):
            for j in range(0, self.PROC_XY):
                MSG_Node_Mode.append(self.LAT_Node_To_Node + (self.Planes[i] * self.Num_Elements[i][j] * self.Data_Mode)/((self.BW_Node_To_Node) * self.PROC_Z))
                MSG_Socket_Mode.append(self.LAT_Socket_To_Socket + (self.Planes[i] * self.Num_Elements[i][j] * self.Data_Mode)/((self.BW_Socket_To_Socket) * self.PROC_Z))
                MSG_Core_Mode.append(self.LAT_Core_To_Core + (self.Planes[i] * self.Num_Elements[i][j] * self.Data_Mode)/((self.BW_Core_To_Core) * self.PROC_Z))
        
        # We are looking for the core that has the most communication to do in this Alltoall set up.
        # Each core performs PROC_Z - 1 sends (A send contains the recieve in this model)
        comm_max = 0.0

        for i in range(0, self.PROC_Z):
            for j in range(0, self.PROC_XY):

                # Counter for the communication for this the current core
                comm_count = 0.0

                # Current column and core
                current = self.Core[i*self.PROC_XY + j]
                current_column = current % self.PROC_XY
                
                # Loop looking for neighbours and counting the communication
                for k in range(0, self.PROC_TOT):
                    if (k % self.PROC_XY == current_column):
                        relation = self.Check_Neighbour(current, k)

                        if (relation == 0):
                            continue

                        if (relation == 1):
                            comm_count += MSG_Core_Mode[current]
                            continue

                        if (relation == 2):
                            comm_count += MSG_Socket_Mode[current]
                            continue

                        if (relation == 3):
                            comm_count += MSG_Node_Mode[current]
                            continue
                
                # Update the largest communication cost
                if (comm_count >= comm_max):
                    comm_max = comm_count

        # Return the result
        return(18 * comm_max)

#------------------------------------
# New Function
#------------------------------------

    def Serial_Compute(self):
        
        serial_max = 0.0

        for i in range(0, self.PROC_Z):
            for j in range(0, self.PROC_XY):

                serial_count = 0.0

                current = self.PROC_XY * i + j

                N_P = 0
                N_V_1 = 0
                N_V_2 = 0
                N_V_3 = 0
                
                try: 
                    length = len(self.Plane_Num[i])
                except:
                    length = 1

                for k in range(0, length):
                    
                    Current_Plane = str(1 + self.Plane_Num[i][k])

                    try:
                        N_P += self.Pressure[Current_Plane][0]
                    except:
                        Turing = 'King of Computers'

                    try:
                        N_V_1 += self.Velocity_1[Current_Plane][0]
                    except:
                        Turing = 'King of Computers'

                    try:
                        N_V_2 += self.Velocity_2[Current_Plane][0]
                    except:
                        Turing = 'King of Computers'

                    try:
                        N_V_3 += self.Velocity_3[Current_Plane][0]
                    except:
                        Turing = 'King of Computers'

                serial_count = Serial_Computation(self.P, self.Num_Elements[i][j], self.Planes[current], N_P, N_V_1, N_V_2, N_V_3, self.Num_Constants, self.constants, self.Scheme)

                if (serial_count >= serial_max):
                    serial_max = serial_count

        return(serial_max)


#------------------------------------
# End of Class File
#------------------------------------
