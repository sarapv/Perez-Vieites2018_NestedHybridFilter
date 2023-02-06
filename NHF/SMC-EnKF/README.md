# SMC-EnKF: an implementation of the nested hybrid filter (NHF)

Code for running a nested hybrid filter (NHF) in a two-scale Lorenz 96 model. In particular, we run an SMC-EnKF, i.e,
  * An SMC for parameter estimation (that estimates the parameter F and the 2 coefficients of the ansatz),
  * A bank of ensemble Kalman filters (EnKFs) that estimate the state.

# Run

The main script is NHF_SMCEnKF_Lorenz96.m
