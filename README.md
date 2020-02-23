# Multiscale-Solver
Multiscale model solver of biomechanical active response of artereis

This code is solving the model published by:
Murtada, Sae-II, Jay D. Humphrey, and Gerhard A. Holzapfel.
"Multiscale and multiaxial mechanics of vascular smooth muscle." Biophysical journal 113.3 (2017): 714-727.

The paper is published in: https://www.sciencedirect.com/science/article/pii/S0006349517306677.

Run a simulation:
1. Run the file Model_Solver.m.
2. Enter the type of loading you would like simulate (currently biaxial or uniaxial, inflation and extension or ring test, respectively).
3. You can also change the loading conditions on the Model_Solver.m file; inner pressure, axial stretch or circumferential stretch.
4. To change arterial wall material parameters; both active and passive, see the varaibles in the artery class in Artery.m.

Plot Results:
To plot the results use Plot_Results.m. This file was used during my research and it is a bit scary at first glance, but there are comments within the code that explain most of the plotting options.
I am planning to organize it in the near future and make it more user friendly.

Simulation_Program.m offers automatic simulations to conduct sensitivity tests. Plot_Results.m knows to analyze and visualize these results.
