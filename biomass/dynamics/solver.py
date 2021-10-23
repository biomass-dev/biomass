import math
import sys
import time
from typing import Callable, List, Optional, Tuple, Union

import numpy as np
from scipy.integrate import ode, solve_ivp
from scipy.integrate._ivp.ivp import OdeResult

__all__ = ["solve_ode", "get_steady_state"]


def solve_ode(
    diffeq: Callable,
    y0: Union[list, np.ndarray],
    t: Union[range, List[int]],
    f_params: Tuple[float, ...],
    *,
    method: str = "BDF",
    vectorized: bool = False,
    options: Optional[dict] = None,
) -> Optional[OdeResult]:
    """
    Solve a system of ordinary differential equations using ``scipy.integrate.solve_ivp()``.

    Parameters
    ----------
    diffeq : callable f(t, y, *x)
        Right-hand side of the differential equation.
    y0 : array
        Initial condition on y (can be a vector).
    t : array
        A sequence of time points for which to solve for y.
    f_params : tuple
        Model parameters.
    method : str (default: "BDF")
        Integration method to use.
    vectorized : bool (default: :obj:`False`)
        Whether `diffeq` is implemented in a vectorized fashion.
    options : dict, optional
        Options passed to a chosen solver.

    Returns
    -------
    sol : OdeResult
        Represents the solution of ODE.
    """
    if options is None:
        options = {}
    options.setdefault("rtol", 1e-8)
    options.setdefault("atol", 1e-8)
    try:
        sol = solve_ivp(
            diffeq,
            (t[0], t[-1]),
            y0,
            method=method,
            t_eval=t,
            vectorized=vectorized,
            args=f_params,
            **options,
        )
        return sol if sol.success else None
    except ValueError:
        return None


def get_steady_state(
    diffeq: Callable,
    y0: list,
    f_params: tuple,
    *,
    dt: float = 1,
    atol: float = 1e-8,
    rtol: float = 1e-8,
    allclose_kws: Optional[dict] = None,
    maximum_wait_time: Union[int, float] = 60.0,
) -> List[float]:
    """
    Simulate a model from given initial conditions until it reaches steady state.

    Parameters
    ----------
    diffeq : callable f(t, y, *x)
        Right-hand side of the differential equation.
    y0 : array
        Initial condition on y (can be a vector).
    f_params : tuple
        Model parameters.
    dt : float (default: 1.0)
        The step size used to calculate the steady state.
    atol : float (default: 1e-8)
        Absolute tolerance for solution.
    rtol : float (default: 1e-8)
        Relative tolerance for solution.
    allclose_kws : dict, optional
        Keyword arguments to pass to ``numpy.allclose()``.
    maximum_wait_time : int or float (default: 60.0 = 1 min.)
        The longest time a user can wait for the system to reach the steady state.

    Returns
    -------
    steady_state : List[float]
        Steady state concentrations of all species.
        Return an empty list if simulation failed.
    """
    if allclose_kws is None:
        allclose_kws = {}
    allclose_kws.setdefault("rtol", 1e-3)
    sol = ode(lambda t, y, f_args: diffeq(t, y, *f_args))
    sol.set_integrator(
        "zvode",
        method="bdf",
        with_jacobian=True,
        atol=atol,
        rtol=rtol,
    )
    sol.set_initial_value(y0, 0)
    sol.set_f_params(f_params)
    ys = [y0]
    start = time.time()
    while sol.successful():
        sol.integrate(sol.t + dt)
        if (
            np.iscomplex(np.real_if_close(sol.y)).any()
            or (time.time() - start) > maximum_wait_time
        ):
            return []
        elif np.allclose(sol.y, ys[-1], **allclose_kws):
            break
        else:
            ys.append(sol.y)
    steady_state = np.real(sol.y).tolist() if sol.successful() else []
    for i, val in enumerate(steady_state):
        if math.fabs(val) < sys.float_info.epsilon:
            steady_state[i] = 0.0
    return steady_state
