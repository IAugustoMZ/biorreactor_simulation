# Simulation of Glicose to Ethanol Fermentation Bioreactor

This project aims to implement a simulation of a bioreactor fermentation of glycose to produce ethanol.

The simulation is implemented by integration of derived rate laws of cell growth, substract consumption and product formation by using the numerical integration algorithm of 4th order Runge Kutta.

The cell growth is described by a modified Monod's rate law expression, which considers the inhibition that ethanol promotes in cell growth.

## Front-End and Model Deploy

The results of simulations are deployed in a local server built using components of Dash libraries and its support for Bootstrap components. The application runs locally in a Flask served environment.

The application consists in a dashboard where the user can select some parameters of the process and automatically visualize its effects on the desired batch profile (cell, substract and product concentrations). It also calculates the necessary installed volume of reaction to attend the inputted volume of ethanol.

To see this dashboard, the user must run the Jupyter notebook called "simulator.ipynb", inside the scripts folder.