
import numpy as np
from functools import partial

from numpy import ndarray, dtype

def truncnorm(mean, std_dev, bounds):
    new_mean = np.zeros_like(mean)
    for i in range(mean.shape[0]):
        idx = np.arange(mean.shape[1])
        while True:
            new_mean[i, idx] = np.random.normal(mean[i, idx], std_dev[i])
            condition = (new_mean[i, :] < bounds[i, 0]) | (new_mean[i, :] > bounds[i, 1])
            idx = np.where(condition)[0]
            if len(idx) == 0:
                break
    return new_mean

def nested_particle_filter(x ,y ,N ,M ,K ,h ,T ,range ,jittering_var ,param):

    dx = x.shape[0]
    dy = y.shape[0]
    dparam = len(param)

    # Estimates
    x_est = np.zeros_like(x)
    param_est = np.zeros((dparam ,T))

    weights_param = np.ones(N ) *( 1 /N)

    Param_particles = np.zeros((dparam ,N))
    for i in range(dparam):
        Param_particles[i ,:] = np.random.uniform(range[i ,0], range[i ,1], N)

    X_all_particles = np.zeros((N, dx, M))
    X_estimates_second_layer = np.zeros((dx, N))

    for t in np.arange(1, T):
        # jittering of Param_particles
        Param_particles = truncnorm(Param_particles, np.sqrt(jittering_var), range)

        log_weights_param = np.zeros(N)
        for m in range(M):
            # TODO: call particle filter for the state give parameters
            print(m)

        unnormalized_weights_param = np.exp(log_weights_param - max(log_weights_param))
        weights_param = unnormalized_weights_param / np.sum(unnormalized_weights_param)

        # Estimates
        x_est[: ,t] = np.cross(X_estimates_second_layer ,weights_param.T)
        param_est[: ,t] = np.cross(Param_particles ,weights_param.T)

        # TODO: resampling




    result = 1

    return result


def particle_filter_for_state_estimation(state ,N ,M ,K ,h):
    result = 1

    return result
