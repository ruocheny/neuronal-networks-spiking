The supplementary codes contains two parts, spiking_simulation for Izhikevich simulation on artificial neuronal networks, and goal_directed_network for spiking neural networks training.

## spiking_simulation

The code is developed with MATLAB and can be run with MATLAB R2024b. 

Please run the source code in the following order:

1. Generate simulated neuronal networks: ```main_genNet.m```
2. Generate spikes with Izhikeivich model on the simulated networks: ```main_genSpikes.m``` (takes ~1 min)
3. Run MFDFA on simulated spike trains: ```main_MFDFA.m```
4. Plot Hurst exponent and multifractal spectrum: ```ana_MFDFA.m```


## goal_directed_network

The code is developed with Python and tension package (https://github.com/zhenruiliao/tension). The ipynb notebook trains spiking neural networks with FORCE algorithm to achieve computational task of integration, differentiation and delay.

This repository includes code modified from the [Tension](https://github.com/zhenruiliao/tension) package, which is licensed under the MIT License. A copy of the original license is provided in ```goal_directed_network/tension```.