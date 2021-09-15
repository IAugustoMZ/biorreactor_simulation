import os
import warnings
import pandas as pd
import numpy as np

# ignore warnings
warnings.filterwarnings('ignore')

# class definition
class BioReactor:

    # constructor method
    def __init__(self, Cc0: float, Cs0: float, production: float, batch_time: float) -> None:
        """
        builds the BioReactor class and defines some important constants
        
        args: initial concentration of cells, initial concentration of substract
            desired production rate of ethanol, desired batch time [h]

        returns: None
        """

        # important constant and definitions
        self.cp_s = 93                          # ethanol's concentration that inhibits cell growth [g/dm3]
        self.n = 0.52                           # empiric constant to express product inhibition
        self.mu_max = 0.33                      # maximum cell growth velocity [1/h]
        self.Ks = 1.7                           # Michaelis-Menten constant [g/dm3]
        self.m = 0.03                           # substract consumption to maintain cell mass [g_subs/g_cell.h]
        self.y_c_s = 0.08                       # cell to substract yield [g_cell / g_subs]
        self.y_p_s = 0.45                       # product to substract yield [g_product / g_subs]
        self.y_p_c = 5.6                        # product to cell yield [g_product / g_cell]
        self.kd = 0.01                          # cell death rate [1/h]
        self.prod = production                  # desired amount of ethanol production in a year [m³]
        self.batch_time = batch_time            # total desired batch time (operation) [h]
        self.ethanol_d = 0.79                   # ethanol density [g/cm3]
        self.setup_time = 6                     # estimated time of setup [h]
        self.volume = None                      # total necessary volume [m3]
        self.operational_days= 300              # number of operational days  by year [days]

        # creating dictionary of states to store iterations results
        self.states_dict = {
            'Cc': Cc0,
            'Cs': Cs0,
            'Cp': 0,
        }

        pass

    def rate_laws(self) -> list:
        """
        calculates the rate laws that rule the fermentation chemical system

        args: None

        returns: list of rate laws at a specific state of the system
        """

        # rate of product inhibition
        if self.states_dict['Cp'] >= self.cp_s:
            k_obs = 0
        else:
            k_obs = (1-(self.states_dict['Cp']/self.cp_s))**(self.n)

        # rate of new cells generation (rg) - [g_cell / h]
        rg = self.mu_max*k_obs*((self.states_dict['Cc']*self.states_dict['Cs'])/(self.Ks+self.states_dict['Cs']))

        # rate of cell death (rd) - [g_cell / h]
        rd = self.kd*self.states_dict['Cc']

        # rate of substract consumption to cell mass maintenance (rsm) - [g_subst / h]
        rsm = self.m*self.states_dict['Cc']

        # rate of product production (rp) - [g_product / h]
        rp = self.y_p_c*rg

        return [rg, rd, rsm, rp]

    def calculate_k(self) -> list:
        """
        calculates the Runge-Kutta's k parameters by applying the mass balance 
        to cells, substract and product concentrations

        args: None

        returns: list of Runge-Kutta's k parameters
        """  

        # calculate the rate law at present state
        r = self.rate_laws()

        # derivative of the concentration of cells
        dCc = r[0] - r[1]

        # derivative of the concentration of substract
        dCs = (1/self.y_c_s)*(-r[0]) - r[2]

        # derivative of the concentration of product
        dCp = r[3]

        return [dCc, dCs, dCp]

    def runge_kutta(self) -> pd.DataFrame:
        """
        applies the runge-kutta integration method to estimate the concentration
        profiles in a batch

        args: None

        returns: dataframe with concentration profiles
        """

        # definition of integration grid size
        N = 500

        # integration grid
        t = np.linspace(0, self.batch_time, N)

        # integration step
        h = t[1] - t[0]

        # profiles dataframe creation
        cols = ['batch_time', 'cell_conc', 'substract_conc', 'product_conc']
        results = pd.DataFrame(columns=cols)

        # extract values of states
        Cc = self.states_dict['Cc']
        Cs = self.states_dict['Cs']
        Cp = self.states_dict['Cp']

        for i in range(len(t)):

            # allocate present state in dataframe
            results.loc[i, cols] = t[i], Cc, Cs, Cp

            # calculate k1
            k1_list = self.calculate_k()

            # update intermediate state dictionary
            self.states_dict['Cc'] = Cc + ((h/2)*k1_list[0])
            self.states_dict['Cs'] = Cs + ((h/2)*k1_list[1])
            self.states_dict['Cp'] = Cp + ((h/2)*k1_list[2])

            # calculate k2
            k2_list = self.calculate_k()

            # update intermediate state dictionary
            self.states_dict['Cc'] = Cc + ((h/2)*k2_list[0])
            self.states_dict['Cs'] = Cs + ((h/2)*k2_list[1])
            self.states_dict['Cp'] = Cp + ((h/2)*k2_list[2])

            # calculate k3
            k3_list = self.calculate_k()

            # update intermediate state dictionary
            self.states_dict['Cc'] = Cc + (h*k3_list[0])
            self.states_dict['Cs'] = Cs + (h*k3_list[1])
            self.states_dict['Cp'] = Cp + (h*k3_list[2])

            # calculate k4
            k4_list = self.calculate_k()

            # update the true value of next state
            Cc += (h/6)*(k1_list[0] + (2*k2_list[0]) + (2*k3_list[0]) + k4_list[0])
            Cs += (h/6)*(k1_list[1] + (2*k2_list[1]) + (2*k3_list[1]) + k4_list[1])
            Cp += (h/6)*(k1_list[2] + (2*k2_list[2]) + (2*k3_list[2]) + k4_list[2])

        return results

    def calculate_volume(self):
        """
        calculates the necessary volume to attend the demanded ethanol production

        args: None

        returns: None
        """

        # execute simulation to calculate the outlet concentration
        results = self.runge_kutta()

        # extract concentration of ethanol at outlet of reactor
        Cet = results.tail(1)['product_conc'].values[0]
        
        # calculate volume of reaction
        self.volume = (self.prod*self.ethanol_d*(self.batch_time+self.setup_time))/(Cet*self.operational_days*24*1000)

    def change_cells(self, new_Ks, new_u):
        """
        method to simulate the modification of microrganism used in fermentation

        args: new Monod and maximum velocity parameters

        returns: None
        """

        # change the microrganism parameters
        self.Ks = new_Ks
        self.mu_max = new_u