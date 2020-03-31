import numpy as np
from scipy.integrate import ode


def solveode(diffeq, y0, tspan, args):
    sol = ode(diffeq)
    sol.set_integrator(
        'vode', method='bdf', with_jacobian=True,
        atol=1e-9, rtol=1e-9, min_step=1e-8
    )
    sol.set_initial_value(y0, tspan[0])
    sol.set_f_params(args)

    T = [tspan[0]]
    Y = [y0]

    while sol.successful() and sol.t < tspan[-1]:
        sol.integrate(sol.t+1.)
        T.append(sol.t)
        Y.append(sol.y)

    return np.array(T), np.array(Y)


def get_steady_state(diffeq, y0, tspan, args, steady_state_eps=1e-6):
    iter_ = 0
    while iter_ < 100:
        (T, Y) = solveode(diffeq, y0, tspan, args)
        if T[-1] < tspan[-1] or np.all(np.abs(Y[-1, :] - y0) < steady_state_eps):
            break
        else:
            y0 = Y[-1, :].tolist()
            iter_ += 1
    
    return T[-1], y0
