# Traffic Simulation: Cellular Automata

This is a simulation of two lane traffic using cellular automata. It is based on the algorithm described in this paper:

M. Rickert, K. Nagel, M. Schreckenberg, A. Latour,  
Two lane traffic simulations using cellular automata,  
Physica A: Statistical Mechanics and its Applications,  
Volume 231, Issue 4,  
1996,  
Pages 534-550,  
ISSN 0378-4371,  
https://doi.org/10.1016/0378-4371(95)00442-4.  
(http://www.sciencedirect.com/science/article/pii/0378437195004424)

traffic_serial.c implements the algorithm serially.  
traffic_parallel.c parallelises that implementation.

Usage:  
./e_traffic_serial\[parallel\] density p_change p_brake seed  
density: float that determines number of cars.  
p_change: float, probability of lane change.  
p_brake: float, probability of a car randomly braking.  
seed: integer used to seed the PRNG.

Lines 78, 79, 84, 85 from traffic_serial.c  
and 90, 91, 96, 97 from traffic_parallel.c can be uncommented for command line animation for testing purposes.
