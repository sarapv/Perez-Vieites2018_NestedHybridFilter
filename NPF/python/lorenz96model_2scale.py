import time

import numpy as np


def fx_fz_lorenz96(x_old, z_old, dx, dz, fps, param):

    F = param[0]
    H = param[1]
    B = param[2]
    C = param[3]

    fx_new = np.zeros_like(x_old)
    fz_new = np.zeros_like(z_old)

    # slow variables
    fx_new[0] = (-x_old[dx - 1] * (x_old[dx - 2] - x_old[1])) - x_old[0] + F - ((H * C) / B) * np.sum(z_old[0:fps])
    fx_new[1] = (-x_old[0] * (x_old[dx - 1] - x_old[2])) - x_old[1] + F - ((H * C) / B) * np.sum(
        z_old[fps + np.arange(0, fps)])

    for i in np.arange(2, dx - 1):
        fx_new[i] = (-x_old[i - 1] * (x_old[i - 2] - x_old[i + 1])) - x_old[i] + F - ((H * C) / B) * np.sum(
            z_old[(i - 1) * fps + np.arange(0, fps)])

    fx_new[dx - 1] = (-x_old[dx - 2] * (x_old[dx - 3] - x_old[0])) - x_old[dx - 1] + F - ((H * C) / B) * np.sum(
        z_old[(dx - 2) * fps + np.arange(0, fps)])

    # fast variables
    xaux = np.repeat(x_old, fps)
    fz_new[0] = (C * B * z_old[1] * (z_old[dz - 1] - z_old[2]) - C * z_old[0] + (F * C / B) + (C * H / B) * xaux[0])

    fz_new[1:dz - 3] = C * B * z_old[2:dz - 2] * (z_old[0:dz - 4] - z_old[3:dz - 1]) - C * z_old[1:dz - 3] + (
            F * C / B) + (C * H / B) * xaux[1:dz - 3]

    fz_new[dz - 2] = C * B * z_old[dz - 1] * (z_old[dz - 3] - z_old[0]) - C * z_old[dz - 2] + (F * C / B) + (
            C * H / B) * xaux[dz - 2]

    fz_new[dz - 1] = (C * B * z_old[0] * (z_old[dz - 2] - z_old[1]) - C * z_old[dz - 1] + (F * C / B)
                      + (C * H / B) * xaux[dz - 1])

    return fx_new, fz_new


def runge_kutta4_lorenz96_1timestep(x_old, z_old, dx, dz, fps, h, param, noise_var):
    # Placeholder for the Euler method implementation
    #print("4th order Runge-Kutta method called")

    s2x = noise_var[0]
    s2z = noise_var[2]

    fx1, fz1 = fx_fz_lorenz96(x_old, z_old, dx, dz, fps, param)

    x1_rk = x_old + 0.5 * (h * fx1 + np.sqrt(h*s2x) * np.random.randn(dx))
    z1_rk = z_old + 0.5 * (h * fz1 + np.sqrt(h*s2z) * np.random.randn(dz))
    fx2, fz2 = fx_fz_lorenz96(x1_rk, z1_rk, dx, dz, fps, param)

    x2_rk = x_old + 0.5 * (h * fx2 + np.sqrt(h*s2x) * np.random.randn(dx))
    z2_rk = z_old + 0.5 * (h * fz2 + np.sqrt(h*s2z) * np.random.randn(dz))
    fx3, fz3 = fx_fz_lorenz96(x2_rk, z2_rk, dx, dz,fps, param)

    x3_rk = x_old + 1 * (h * fx3 + np.sqrt(h*s2x) * np.random.randn(dx))
    z3_rk = z_old + 1 * (h * fz3 + np.sqrt(h*s2z) * np.random.randn(dz))
    fx4, fz4 = fx_fz_lorenz96(x3_rk, z3_rk, dx, dz, fps, param)

    x_new = x_old + (1 / 6) * (h * (fx1 + fx2 + fx3 + fx4) + np.sqrt(h*s2x) * np.random.randn(dx))
    z_new = z_old + (1 / 6) * (h * (fz1 + fz2 + fz3 + fz4) + np.sqrt(h*s2z) * np.random.randn(dz))

    return x_new, z_new


def euler_lorenz96_1timestep(x_old, z_old, dx, dz, fps, h, param, noise_var):
    # Placeholder for the Runge-Kutta 4 method implementation
    #print("Euler method called")
    s2x = noise_var[0]
    s2z = noise_var[2]

    fx, fz = fx_fz_lorenz96(x_old, z_old, dx, dz, fps, param)

    x_new = x_old + h * fx + np.sqrt(h*s2x) * np.random.randn(dx)
    z_new = z_old + h * fz + np.sqrt(h*s2z) * np.random.randn(dz)

    return x_new, z_new


def lorenz96_2scale(dx, fps, h, t_final, integration_method):
    """
    Simulates a Lorenz96 model :true_param integration_method: "euler" or "runge-kutta4" are available to integrate the
    SDEs of Lorenz96 :true_param dx: dimension of the state x, e.g., dx = 20 :true_param fps: number of fast variables per each
    slow variables x (i.e., dimension of z = fps * dx), e.g., fps = 10 :true_param h: step-size for the integration of the
    SDEs (in natural units), e.g., h = 2e-4 :true_param t_final: length of the simulation in natural time units (i.e.,
    t_final/h discrete time steps), e.g., t_final = 10

    :return: x: slow state variables from t = 0, ... , T
    :return: y: (partial) observations of the slow state variables from t = 0, ... , T
    :return: z: fast state variables from t = 0, ... , T
    :return: u: (partial) observations of the fast state variables from t = 0, ... , T
    :return: true_param: contains the parameters of the Lorenz96 model
    """

    # Dimension of fast variables
    dz = fps * dx  # no. of fast variables (dz)

    # Static parameters (Lorenz 96)
    F = 8  # forcing parameter
    H = 0.75  # coupling between slow and fast variables
    C = 10  # time scale of variables y
    B = 15  # inverse amplitude of the fast variables
    param = [F, H, B, C]

    # Noise parameters (variance)
    s2y = 4  # variance of the observations: slow variables
    s2u = 1 / 10  # variance of the observations: fast variables
    s2x = 1 / 2  # variance of the state noise: slow variables
    s2z = 1 / 8  # variance of the state noise: fast variables
    noise_var = [s2x, s2y, s2z, s2u]

    # Time variables
    T = np.fix(t_final / h).astype(int)  # no. of discrete time steps

    # Initialize variables
    x = np.zeros((dx, T))
    z = np.zeros((dz, T))

    # Initial conditions
    x[:, 0] = np.random.rand(dx)
    z[:, 0] = (1 / (C * B)) * np.random.rand(dz) - (1 / (2 * C * B))

    # Capture start time
    t0 = time.time()
    ok = 0
    while not ok:
        for t in np.arange(1, T):
            if integration_method == "euler":
                x[:, t], z[:, t] = euler_lorenz96_1timestep(x[:, t - 1], z[:, t - 1], dx, dz, fps, h, param, noise_var)

            elif integration_method == "runge-kutta4":
                x[:, t], z[:, t] = runge_kutta4_lorenz96_1timestep(x[:, t - 1], z[:, t - 1], dx, dz, fps, h, param,
                                                                   noise_var)

            else:
                raise ValueError("Invalid method. Choose 'euler' or 'runge-kutta4'.")

        ok = not np.any(np.isnan(x[0, :]) | np.isinf(x[0, :])) | np.any(np.isnan(z[0, :]) | np.isinf(z[0, :]))

    y = x + np.sqrt(s2y) * np.random.randn(dx, T)
    u = z + np.sqrt(s2u) * np.random.randn(dz, T)

    # Compute elapsed time
    elapsed_time = time.time() - t0
    print(f"time: {elapsed_time:5.4f} s")

    return x, z, y, u, param, noise_var, T
