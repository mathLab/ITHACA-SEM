# Andrew Gloster
# Summer 2016
# Python Version of Communication Model Implementation
# Functions relating to serial computations and fitting of serial element of the model

#------------------------------------
# Import relevant modules
#------------------------------------

from scipy.optimize import fmin
import matplotlib.pyplot as plt

#------------------------------------
# Import relevant functions and classes
#------------------------------------

import math as m
import numpy as np
from scipy.optimize import minimize

#------------------------------------
# New Function
#------------------------------------

# Compute length of one timestep for given inputs (Model must be fitted prior to this function)
def Serial_Computation(P, Num_Elements, Num_Modes, N_P, N_V_1, N_V_2, N_V_3, Num_Constants, constants, Scheme):

# Operation Counts as per the thesis
    O_A_1 = 9 * Num_Elements * (P + 1) ** 2 * Num_Modes * m.log(Num_Modes, 2)
    O_A_2 = Num_Elements * Num_Modes * (P + 1) ** 2
    O_A_3 = 6 * Num_Elements * (P + 1) ** 4 * Num_Modes
    O_A_4 = 15 * Num_Elements * (P + 1) ** 2 * Num_Modes

    T_A = O_A_1 + O_A_2 + O_A_3 + O_A_4

    if (Scheme == 'IterativeFull'):
        O_E_1 = 8 * Num_Elements * (P + 1) ** 2
        O_E_2 = Num_Elements * (P + 1) ** 2
        O_E_3 = Num_Elements * ((4 * P ** 3) + (18 * P ** 2) + (26 * P) + 12)
        O_E_4 = 6 * Num_Elements * (P + 1) ** 2

        t_e = O_E_1 + O_E_2 + O_E_3 + O_E_4

        T_E = (N_P + N_V_1 + N_V_2 + N_V_3) * t_e

    if (Scheme == 'IterativeStaticCond'):
        O_E_1 = 8 * ((P - 1) ** 2) * 4 * P * Num_Elements * Num_Modes
        O_E_2 = (N_P + N_V_1 + N_V_2 + N_V_3) * 2 * 16 * (P ** 2) * Num_Elements
        O_E_3 = (8 * 4 * ((P - 1) ** 2)) + (4 * ((P - 1) ** 2)) * Num_Elements * Num_Modes
        O_E_4 = 8 * ((((P - 1) ** 2)) ** 2) * Num_Elements * Num_Modes

        T_E = O_E_1 + O_E_2 + O_E_3 + O_E_4

# Divide by operations per second.
    
    if (Num_Constants == 1):
        T_A = T_A / constants
        T_E = T_E / constants

    if (Num_Constants == 2):
        T_A = T_A / constants[0]
        T_E = T_E / constants[1]

    Time = T_A + T_E
    
    return(Time)

#------------------------------------
# New Function
#------------------------------------

def Operation_Count(P, Num_Elements, Num_Modes, N_P, N_V_1, N_V_2, N_V_3, Scheme):

# Operation Counts as per the thesis
    O_A_1 = 9 * Num_Elements * (P + 1) ** 2 * Num_Modes * m.log(Num_Modes, 2)
    O_A_2 = Num_Elements * Num_Modes * (P + 1) ** 2
    O_A_3 = 6 * Num_Elements * (P + 1) ** 4 * Num_Modes
    O_A_4 = 15 * Num_Elements * (P + 1) ** 2 * Num_Modes

    T_A = O_A_1 + O_A_2 + O_A_3 + O_A_4

    if (Scheme == 'IterativeFull'):
        O_E_1 = 8 * Num_Elements * (P + 1) ** 2
        O_E_2 = Num_Elements * (P + 1) ** 2
        O_E_3 = Num_Elements * ((4 * P ** 3) + (18 * P ** 2) + (26 * P) + 12)
        O_E_4 = 6 * Num_Elements * (P + 1) ** 2

        t_e = O_E_1 + O_E_2 + O_E_3 + O_E_4

        T_E = (N_P + N_V_1 + N_V_2 + N_V_3) * t_e

    if (Scheme == 'IterativeStaticCond'):
        O_E_1 = 8 * ((P - 1) ** 2) * 4 * P * Num_Elements * Num_Modes
        O_E_2 = (N_P + N_V_1 + N_V_2 + N_V_3) * 2 * 16 * (P ** 2) * Num_Elements
        O_E_3 = (8 * 4 * ((P - 1) ** 2)) + (4 * ((P - 1) ** 2)) * Num_Elements * Num_Modes
        O_E_4 = 8 * ((((P - 1) ** 2)) ** 2) * Num_Elements * Num_Modes

        T_E = O_E_1 + O_E_2 + O_E_3 + O_E_4

    return(T_A, T_E)

#------------------------------------
# New Function
#------------------------------------

# Find the L_2 norm of the difference between the data and model, used for fitting.
def compare_data(constants, Num_Constants, Data, T_A, T_E):
    
    # Lemons make nice Lemonade (Needed a list to hold results)
    lemons = []

    # Loop over data calculating difference between data and model, this will be minimised by the constants.
    for i in range(0, len(Data)):

        if (Num_Constants == 1):
            lemons.append(Data[i] - (T_A[i] + T_E[i])/constants)

        if (Num_Constants == 2):
            lemons.append(Data[i] - (T_A[i] / constants[0] + T_E[i] / constants[1]))

    # Calculate the L_2 norm of the difference, this is the quantity to be minimised by the optimisation step
    L_2_norm = np.linalg.norm(lemons, 2)

    return(L_2_norm) 

#------------------------------------
# New Function
#------------------------------------

# Function to run the fmin optimisation algorithm to fit the model to the data.
def Fit_Model(Num_Constants, Data, T_A, T_E):

    # Numpy arrays to hold the results and data.
    Data = np.array(Data)
    T_A = np.array(T_A)
    T_E = np.array(T_E)
    
    # Set the initial values depending on how many constants are required
    if (Num_Constants == 1):
        inital = 1e6
        
    if (Num_Constants == 2):
        inital = np.array([1e6, 1e06])
    
    # Run the optimisation algorithm
    Fit = fmin(compare_data, inital, args=(Num_Constants, Data, T_A, T_E), xtol=0.0001, ftol=1, maxiter=1e04, maxfun=1e09)
 
    # Return the fitted constants.
    return(Fit)

#------------------------------------
# New Function
#------------------------------------

# Function to run all the steps in order to fit the serial model to serially produced data.
def Run_Serial_Fit(Compare_Serial, Consider_Modes, Num_Constants, P, Num_Elements, Nektar_Modes, Timings, Pressure, Velocity_1, Velocity_2, Velocity_3, Scheme):
    
    # List to hold data from Nektar simulation
    Data = []

    # Loop over modes to be considered in calibrating the model
    for i in range(0, len(Consider_Modes)):
        Data.append(np.mean(Timings[str(Consider_Modes[i])])/10)
    
    # Lists to hold the operation counts to be fed into the optimisation
    T_A = []
    T_E = []
    
    # Loop over modes to be used in the model
    for i in range(0, len(Consider_Modes)):

        # Initialise CG counters to 0
        N_P = 0
        N_V_1 = 0
        N_V_2 = 0
        N_V_3 = 0
        
        # Loop over all the modes up to the mode being considered adding the CG iterations from each mode to the counter
        for j in range(1, Consider_Modes[i] + 1):
            
            # Skip the 2nd mode as we do not get any data about it from Nektar
            if (j == 2):
                continue 
            
            # We don't necessarily have data for all the CG iterations hence we use try
            # The except case has to be included to keep python happy hence we pay homage to Alan Turing
            try:
                N_P += Pressure[str(j)][0]
            except:
                Turing = 'King of Computers'

            try:
                N_V_1 += Velocity_1[str(j)][0]
            except:
                Turing = 'King of Computers'

            try:
                N_V_2 += Velocity_2[str(j)][0]
            except:
                Turing = 'King of Computers'

            try:
                N_V_3 += Velocity_3[str(j)][0]
            except:
                Turing = 'King of Computers'
         
        # Perform the operation counts for each mode and store these to be fed into the optimisation step
        (t_a, t_e) = Operation_Count(P, Num_Elements, Consider_Modes[i], N_P, N_V_1, N_V_2, N_V_3, Scheme)
        T_A.append(t_a)
        T_E.append(t_e)

    # Now fit the model to the data
    Fit = Fit_Model(Num_Constants, Data, T_A, T_E)
    
    # Print the fitted FLOPs
    print(Fit)
    
    # Output the comparison the results if the user desires
    if Compare_Serial is True:
        
        # Pharse the Nektar data as before
        Data = []

        for i in range(1, len(Nektar_Modes)):
            Data.append(np.mean(Timings[str(Nektar_Modes[i])])/10)

        Time = []

        for i in range(1, len(Nektar_Modes)):
            N_P = 0
            N_V_1 = 0
            N_V_2 = 0
            N_V_3 = 0

            for j in range(1, Nektar_Modes[i] + 1):

                if (j == 2):
                    continue 
        
                try:
                    N_P += Pressure[str(j)][0]
                except:
                    Turing = 'King of Computers'

                try:
                    N_V_1 += Velocity_1[str(j)][0]
                except:
                    Turing = 'King of Computers'

                try:
                    N_V_2 += Velocity_2[str(j)][0]
                except:
                    Turing = 'King of Computers'

                try:
                    N_V_3 += Velocity_3[str(j)][0]
                except:
                    Turing = 'King of Computers'

            Time.append(Serial_Computation(P, Num_Elements, Nektar_Modes[i], N_P, N_V_1, N_V_2, N_V_3, Num_Constants, Fit, Scheme))

        Nektar_Modes = list(Nektar_Modes)
        Nektar_Modes.pop(0)
        
        # Calculate the mean, standard deviation and variance of the difference between the data and the fitted model
        difference = []

        for i in range(0, len(Nektar_Modes)):
            difference.append(abs(Data[i] - Time[i]))
        
        mean_diff = np.mean(difference) 
        std_dev_diff = np.std(difference)
        var_diff = np.var(difference)
        
        # Print these results for the user to see
        print('The mean of the differences between the Data and the Model is ' + str(mean_diff))
        print('The standard deviation of the differences between the Data and the Model is ' + str(std_dev_diff))
        print('The variance of the differences between the Data and the Model is ' + str(var_diff))

        # Plot the two data sets together for a visual comparison
        fig, ax = plt.subplots()
        ax.plot(Nektar_Modes, Data, label = 'Data')
        ax.errorbar(Nektar_Modes, Time, label = 'Model')
        ax.set_xlabel('$ N_Z $')
        ax.set_ylabel('Timestep (s)')
        ax.set_title('Length of Single Timestep: Model vs Data')
        plt.legend(loc=4)
        fig.savefig("Output/Figures/Model_vs_Data.png")

    return(Fit)

#------------------------------------
# End of Functions
#------------------------------------