import numpy as np
from scipy.integrate import ode


def solveode(diffeq, y0, tspan, args):
    sol = ode(diffeq)
    sol.set_integrator(
        'vode', method='bdf', with_jacobian=True, min_step=1e-8
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


def get_steady_state(diffeq, y0, tspan, args,
                     steady_state_time=1000000, steady_state_eps=1e-6):
    sol = ode(diffeq)
    sol.set_integrator(
        'vode', method='bdf', with_jacobian=True, min_step=1e-8
    )
    sol.set_initial_value(y0, 0)
    sol.set_f_params(args)

    T = [0]
    Y = [y0]

    while sol.successful() and sol.t < steady_state_time:
        sol.integrate(steady_state_time, step=True)
        if tspan[-1] < sol.t and np.all(np.abs(sol.y - Y[-1]) < steady_state_eps):
            break
        else:
            T.append(sol.t)
            Y.append(sol.y)

    return T[-1], Y[-1]