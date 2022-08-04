# Nested hybrid filters (NHFs)

In this repository is included MATLAB code of the NHFs introduced in [(1)](#references). Other algorithms, which have been compared to NHFs in this paper, are also included:

* /NHFs : 4 implementations of NHFs using sequential Monte Carlo (SMC) or sequential quasi-Monte Carlo (SQMC) in the first layer, and extended Kalman filters (EKFs) or ensemble Kalman filters (EnKFs) in the second layer (i.e., SMC-EKF, SMC-EnKF, SQMC-EKF and SQMC-EnKF).
* /NPF : nested particles filter (NPF)
* /2StageF : two-stage filter 


# References
[1] Pérez-Vieites, S., Mariño, I. P., & Míguez, J. (2018). Probabilistic scheme for joint parameter estimation and state prediction in cojmplex dynamical systems. Physical Review E, 98(6), 063305.

