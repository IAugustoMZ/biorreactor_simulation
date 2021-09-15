import os
import sys
import warnings
import pandas as pd
import numpy as np
sys.path.append('./')
from modules.biorreactor import BioReactor

# ignore warnings
warnings.filterwarnings('ignore')

# experimental deviations
error_cel = 0.000025
error_subst = 0.0015

# class instancing
batch_time_exp = 3.0

# number of experiments
N = 10

# define sampling time
sample_freq = 10                # 10 min
sample_time = np.arange(0, batch_time_exp*60.0+sample_freq, sample_freq)

# execute simulations
experimental_data = []
for i in range(N):

    br = BioReactor(Cc0 = 1.0, Cs0 = 250.0, production = 0.0, batch_time = batch_time_exp)
    results = br.runge_kutta()

    # calculate time in minutes and in seconds
    results['time_min'] = (results['batch_time']*60).astype(int)
    results['time_delta'] = ((results['batch_time']*60) - results['time_min'])*60

    # estimate reaction rate of cell formation
    results['rg'] = (br.mu_max*results['cell_conc']*results['substract_conc'])/(br.Ks + results['substract_conc'])

    # drop NaN
    results = results.dropna()

    # select only sample times from generated dataframe
    samples = []
    for index in results.index:
        if(results.loc[index, 'time_min'] in sample_time):
            samples.append(index)

    results = results.loc[samples,:]
    results = results.loc[results['time_delta'] < 30]
    results = results.drop_duplicates(subset=['time_min'])

    # column projection for desired columns
    results = results[['rg','cell_conc', 'substract_conc']]

    # add gaussian distributed noise
    results['cell_conc'] += np.random.normal(loc = 0.0, scale= error_cel*results['cell_conc'], size=results.shape[0])
    results['substract_conc'] += np.random.normal(loc = 0.0, scale= error_subst*results['substract_conc'], 
                                                    size=results.shape[0])
    
    # append to dataframe list
    experimental_data.append(results)

# concatenate all datasets
full_exp_data = pd.concat(experimental_data, axis=0).reset_index().drop(['index'], axis = 1)

# add gaussian distributed noise to rate of reaction estimation
full_exp_data['rg'] += np.random.normal(loc = 0.0, scale= 0*error_cel*full_exp_data['cell_conc'], 
                                            size=full_exp_data.shape[0])

# save experimental dataset
try:
    full_exp_data.to_csv(os.path.join('./data', 'experimental_data.csv'), index=False)
except:
    os.makedirs('./data')
    full_exp_data.to_csv(os.path.join('./data', 'experimental_data.csv'), index = False)