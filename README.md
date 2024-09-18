# Distance-Estimation-in-Molecular-Communication

**Overview**

This repository contains MATLAB codes developed for two subsequent papers on distance estimation in molecular communication (MC) systems. The first paper focuses on machine learning-based distance estimation methods, while the second paper introduces a novel Fluid Dynamics-based Distance Estimation (FDDE) algorithm. The codes in this repository are primarily from the FDDE paper, with the main components from the machine learning-based distance estimation methods also included for comparison purposes. The code file names correspond to the figures in the FDDE paper (Figures 4-8).

Additionally, the repository includes data files for flow rate measurements (Flow_Rate_Data.xlsx) and distance estimation measurements (T21.mat).

**Citation**

Please cite the following papers when using or building upon the codes in this repository:

Gulec, F., & Atakan, B. (2020). "Distance estimation methods for a practical macroscale molecular communication system." Nano Communication Networks, 24, 100300. https://doi.org/10.1016/j.nancom.2020.100300

Gulec, F., & Atakan, B. (2021). "Fluid dynamics-based distance estimation algorithm for macroscale molecular communication." Nano Communication Networks, 28, 100351. https://doi.org/10.1016/j.nancom.2021.100351

**Background**

Molecular communication (MC) is a cutting-edge communication paradigm where information is transmitted through the release of molecules into the environment. In practical macroscale MC systems, one of the key challenges is estimating the distance between the transmitter (TX) and the receiver (RX).

**Fluid Dynamics-based Distance Estimation (FDDE)**

In the second paper, a novel approach is introduced based on fluid dynamics. The FDDE algorithm models transmitted molecules as evaporating droplets, with their diameter being updated at each time step. The model accounts for the volume fraction of droplets in a mixture of air and liquid, the beamwidth of the transmitter, and the evaporation rate. The FDDE algorithm is validated using experimental data, demonstrating accurate distance estimation in practical MC systems.

**Machine Learning-based Distance Estimation**

In the first paper, machine learning techniques are applied to estimate distances in MC systems. This method serves as a comparison with the fluid dynamics-based approach, providing insights into the advantages and limitations of each method.

**Code Overview**

The following MATLAB codes are included in the repository:

FDDE Algorithm Codes (Figures 4-8):
MATLAB code files corresponding to the figures in the FDDE paper, simulating the fluid dynamics approach for distance estimation.
Experimental data validation is included.
Machine Learning-based Distance Estimation:
Main components from the first paper on machine learning-based distance estimation, included for comparison with the FDDE approach.
Data Files:
Flow_Rate_Data.xlsx: Contains flow rate measurements used to estimate the flow rate.
T21.mat: Contains distance estimation measurements used in the simulations.

**Results**

The results from the simulations and experiments show the following:

**Fluid Dynamics-based Distance Estimation (FDDE):**

The FDDE algorithm provides reliable distance estimation by modeling the propagation of evaporating molecular droplets and accounting for factors such as the volume fraction of droplets and the transmitter (TX) beamwidth.
FDDE performs better than the power-based (PBE) and combined estimation (CE) methods for shorter distances (up to 1.5 meters). Its performance is comparable to machine learning-based methods like linear regression (LR) and neural network regression (NNR) at these distances.
For longer distances (after 1.6 meters), the error in FDDE increases due to a drop in droplet velocity and the influence of diffusion. Machine learning methods show better performance in these scenarios, but they come at the cost of longer data collection and feature extraction times.

**Comparison:**

The FDDE algorithm offers a significant advantage in terms of computational simplicity. Unlike machine learning approaches that require extensive data collection, feature extraction, and learning processes, FDDE relies on physically measurable parameters and simpler calculations.
While machine learning methods are more accurate for long-distance estimation, FDDE can achieve good results with lower computational complexity and without needing large datasets.

**How to Cite**

If you use or build upon these codes in your research, please ensure that you cite both papers:

Gulec, F., & Atakan, B. (2020). "Distance estimation methods for a practical macroscale molecular communication system." Nano Communication Networks, 24, 100300. https://doi.org/10.1016/j.nancom.2020.100300

Gulec, F., & Atakan, B. (2021). "Fluid dynamics-based distance estimation algorithm for macroscale molecular communication." Nano Communication Networks, 28, 100351. https://doi.org/10.1016/j.nancom.2021.100351
