# SQMC-EnKF: an implementation of the nested hybrid filter (NHF)

Code for running a nested hybrid filter (NHF) in a two-scale Lorenz 96 model. In particular, we run an SQMC-EnKF, i.e,

* An SQMC for parameter estimation (that estimates the parameter F and the 2 coefficients of the ansatz),
* A bank of EnKFs that estimate the state.

# Run

The main script is main_SQMCEnKF.m
