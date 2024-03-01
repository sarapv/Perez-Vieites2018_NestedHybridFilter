import numpy as np
import matplotlib.pyplot as plt
from lorenz96model_2scale import lorenz96_2scale

# Ground truth of the model : stochastic Lorenz 96 model


dx = 20
fps = 10
h = 2e-4
t_final = 10
t_obs = 200

x, z, y, u, param, noise_var, T = lorenz96_2scale(dx, fps, h, t_final, integration_method="runge-kutta4")

# Figures
plt.figure(1)
idx_plotted = [1, 2]
auxT = T
for i, idx in enumerate(idx_plotted):
    plt.subplot(len(idx_plotted), 1, i + 1)
    plt.plot(np.arange(0, auxT)*h, x[idx, 0:auxT], 'k', label='state' if i == 0 else "")
    plt.plot(np.arange(0, auxT, t_obs)*h, y[idx, 0:auxT:t_obs], '*', label='observations' if i == 0 else "")
    if i == 0:  # Add legend only to the first subplot to avoid repetition
        plt.legend()
plt.suptitle('Slow state and observations')

plt.figure(2)
idx_plotted_z = [1, 2]
for i, idx in enumerate(idx_plotted_z):
    plt.subplot(len(idx_plotted_z), 1, i + 1)
    plt.plot(np.arange(0, auxT)*h, z[idx, 0:auxT], 'k', label='state' if i == 0 else "")
    plt.plot(np.arange(0, auxT, t_obs)*h, u[idx, 0:auxT:t_obs], '*', label='observations' if i == 0 else "")
    if i == 0:  # Add legend only to the first subplot to avoid repetition
        plt.legend()
plt.suptitle('Fast state and observations')
plt.show()

# Call the NPF : probably we could show the whole algorithm in a high level way (with small functions everywhere)
