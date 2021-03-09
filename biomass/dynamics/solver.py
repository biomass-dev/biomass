import time
import warnings
from typing import Callable, List, Optional, Tuple, Union

import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate._ivp.ivp import OdeResult

__all__ = ["solve_ode", "get_steady_state"]


def solve_ode(
    diffeq: Callable,
    y0: Union[list, np.ndarray],
    t: Union[range, List[int]],
    f_params: Tuple[float, ...],
    method: str = "BDF",
    options: Optional[dict] = None,
) -> Optional[OdeResult]:
    """
    Solve a system of ordinary differential equations.

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
    step_size: int = 1,
    maximum_wait_time: Union[int, float] = 60,
    eps: float = 1e-6,
) -> List[float]:
    """
    Find the steady state for the untreated condition.

    Parameters
    ----------
    diffeq : callable f(t, y, *x)
        Right-hand side of the differential equation.

    y0 : array
        Initial condition on y (can be a vector).

    f_params : tuple
        Model parameters.

    step_size : int (default: 1)
        The step size used to calculate the steady state.

    maximum_wait_time : int or float (default: 60)
        The longest time a user can wait for the system to reach steady state.
        Default value is 60 sec. = 1 min.

    eps : float (default: 1e-6)
        Run until a time t for which the maximal absolute value of the
        regularized relative derivative was smaller than eps.

    Returns
    -------
    y0 : array
        Steady state concentrations of all species.
        Return an empty list if simulation failed.
    """
    start = time.time()
    while True:
        sol = solve_ode(
            diffeq,
            y0,
            range(step_size + 1),
            f_params,
        )
        elapsed_time = time.time() - start
        if maximum_wait_time < elapsed_time:
            warnings.warn(
                f"couldn't reach steady state within maximum_wait_time(={maximum_wait_time}).",
                RuntimeWarning,
            )
            return []
        elif (
            sol is None
            or np.max(np.abs((sol.y[:, -1] - y0) / ([(val + eps) for val in y0]))) < eps
        ):
            break
        else:
            y0 = sol.y[:, -1]

    return [] if sol is None else sol.y[:, -1].tolist()
