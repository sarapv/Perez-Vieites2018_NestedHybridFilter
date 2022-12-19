# Nested hybrid filters (NHFs)

In this repository is included MATLAB code of the NHFs introduced in [[2]](#references). Other algorithms, which have been compared to NHFs in this paper, are also included:

* /NHF : 4 implementations of NHFs using sequential Monte Carlo (SMC) or sequential quasi-Monte Carlo (SQMC) in the first layer, and extended Kalman filters (EKFs) or ensemble Kalman filters (EnKFs) in the second layer (i.e., SMC-EKF, SMC-EnKF, SQMC-EKF and SQMC-EnKF).
* /NPF : implementation of nested particle filter (NPF) [[1]](#references).
* /two-stage filter : implementation of two-stage filter [[3]](#references).


# References
[1] Crisan, D., & Miguez, J. (2018). Nested particle filters for online parameter estimation in discrete-time state-space Markov models. Bernoulli, 24(4A), 3039-3086.

[2] Pérez-Vieites, S., Mariño, I. P., & Míguez, J. (2018). Probabilistic scheme for joint parameter estimation and state prediction in cojmplex dynamical systems. Physical Review E, 98(6), 063305.

[3] Santitissadeekorn, N., & Jones, C. (2015). Two-stage filtering for joint state-parameter estimation. Monthly Weather Review, 143(6), 2028-2042.


# MIT License

Copyright (c) 2022 Sara Pérez Vieites

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
