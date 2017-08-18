# Event Identification Simulation Datasets
1, All the simulation datasets are based on the 68-bus power system with the models as follows: 
This is a 16-machine system with 86 transmission lines and
% 68 buses. Data are extracted from the GE final report
% entitled "Singular Perturbations, Coherency and
% Aggregation of Dynamic Systems," pp.6-42, July 1981.
% detailed generator models
% dc exciters on generators 1 to 9
% constant current active loads

2, 6 types of events are simulated including: Line Trip, faults(including asymetric faults of line to line events and symetric faults of three phase short circuits events), load change, generator trip, capacitor banking, inductive motor tripping.

3, Setup of simulations:
1, simulation time step is 0.01 second but sampling rate is 34 per second;
2, Events occur at 0.5 second and the faults are cleared after 0.2 second;
3, Line trip events means tripping the line without faults;
4, Faults events often are accompanied with tripping the line after faults are cleared, which are different from the line trip described in 3;
5, Load change events are simulated by adding a sudden disturbance into the load control input by the function ml_sig;
6, One second of post-event data matrix is selected to identify the type of events, and notice that the data matrix of faults events are selected after faults clearance; 
7, The model of capacitor banking is based on paper[];
8, All the voltage measurements are in the 'bus_v' and its absolute values or the voltage magnitudes are utilized in our algorithm.
