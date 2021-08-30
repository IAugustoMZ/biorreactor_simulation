import os
import sys
import warnings
import numpy as np
import pandas as pd
sys.path.append('./')
from modules.biorreactor import BioReactor

# ignore warnings
warnings.filterwarnings('ignore')

# class instance
br = BioReactor(Cc0=1.0, Cs0=250.0, production=1000.0, batch_time=12)

# testing the construction method
# print(br.states_dict)

# testing the rate laws calculation
# print(br.rate_laws())

# testing Runge_Kutta's k parameter calculation
# print(br.calculate_k())

# testing the Runge-Kutta's implementation
results = br.runge_kutta()
print(results)