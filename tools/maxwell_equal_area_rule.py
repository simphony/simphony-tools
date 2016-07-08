"""Computation of liquid-vapour densities.

The densities are computed by utilizing the Maxwell's equal area rule
(aka Maxwell construction).
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, fminbound
from scipy.signal import argrelmin, argrelmax


def _foff(arg1, arg2, f, C):
    return f(arg1, arg2) - C


def _finters(a, b, f, arg2, C):
    """Find an intersection of a 2 argument function f and a constant C.

    Parameters
    ----------
    a : float
        Lower limit for the search interval
    b : float
        Upper limit for the search interval
    f : function of 2 arguments
    arg2 : float
        Second argument of f
    C : float
        Constant

    Returns
    -------
    The intersection : float

    Raises
    ------
    RuntimeError : If an intersection is not found.

    """
    # brentq from scipy
    x0, r = brentq(_foff, a, b, (arg2, f, C), full_output=True)

    if not r.converged:
        raise RuntimeError("An intersection not found: abort!")

    return x0


def _area_calc(Pr, Tr, reosf, head_ab, body_ab, tail_ab):
    """Compute the Maxwell's areas.

    That is, compute the areas demarcated
    by the isotherm and a reduced pressure.

    Parameters
    ----------
    Pr : float
        Reduced pressure, aka constant pressure
    Tr : float
        Reduced temperature
    reosf : function
        Reduced equation of state
    head_ab : iterable of float (at least two elements)
        Head part of the isotherm (xstart and xend), aka head
    body_ab : iterable of float (at least two elements)
        Body part of the isotherm (xstart and xend), aka body
    tail_ab : iterable of float (at least two elements)
        Tail part of the isotherm (xstart and xend), aka tail

    Returns
    -------
    A tuple with two elements:
        1. A tuple including the intersections (floats) between the constant
           pressure and the head (1.1.), body (1.2.), and tail (1.3.) part.

        2. A tuple including the area
           2.1. below the constant pressure (negative float),
           2.2. above the constant pressure (positive float).

    """
    h0 = _finters(head_ab[0], head_ab[1], reosf, Tr, Pr)
    b0 = _finters(body_ab[0], body_ab[1], reosf, Tr, Pr)
    t0 = _finters(tail_ab[0], tail_ab[1], reosf, Tr, Pr)

    # quad from scipy
    a1 = quad(reosf, h0, b0, (Tr,))[0] - (b0 - h0)*Pr
    a2 = quad(reosf, b0, t0, (Tr,))[0] - (t0 - b0)*Pr

    return ((h0, b0, t0), (a1, a2))


def _area_diff(Pr, Tr, reosf, head_ab, body_ab, tail_ab):
    """Compute the difference between Maxwell's areas.

    That is, compute the difference between areas
    demarcated by the isotherm and a reduced pressure.

    Parameters
    ----------
    Pr : float
        Reduced pressure, aka constant pressure
    Tr : float
        Reduced temperature
    reosf : function
        Reduced equation of state
    head_ab : iterable of float (at least two elements)
        Head part of the isotherm (xstart and xend), aka head
    body_ab : iterable of float (at least two elements)
        Body part of the isotherm (xstart and xend), aka body
    tail_ab : iterable of float (at least two elements)
        Tail part of the isotherm (xstart and xend), aka tail

    Returns
    -------
    The difference between Maxwell's areas (absolute value).

    """
    (intersects, areas) = _area_calc(Pr, Tr, reosf, head_ab, body_ab, tail_ab)
    return np.fabs(areas[0] + areas[1])


def liquid_vapour_density(Trs, reos, iSearchMax=40.0, RhorMin=1e-12):
    """Compute liquid-vapour densities for the given reduced temperatures.

    The densities are computed by utilizing the Maxwell's equal area rule.

    Parameters
    ----------
    Trs : iterable of float
        Reduced temperatures
    reos :  ReducedEquationOfState
    iSearchMax : float
        Upper limit for the reduced molar volume, Vrmol, in the
        search of the primary extrema (default value 40.0)
    RhorMin : float
        Lower limit for the reduced (vapour) density; limits the range
        of solutions, i.e. a constraint for the equal area optimization
        problem (default value 1e-12; VrmolMax := 1/RhorMin)

    Returns
    -------
    List of tuples, where each tuple has 6 elements (floats):
        1. Reduced temperature, and associated
        2. Reduced equilibrium pressure,
        3. Reduced molar volume for the vapour phase,
        4. Reduced molar volume for the liquid phase,
        5. Reduced density for the vapour phase, and
        6. Reduced density for the liquid phase.

    Raises
    ------
    ValueError :
        If a reduced temperature parameter is >= 1.
    RuntimeError :
        If a primary minimum or maximum of an isotherm is not found.
    RuntimeError :
        If a reduced equilibrium pressure is not found (i.e. solution
        to the equal area optimization problem is not converged).

    Algorithm
    ---------
    1. Find the primary minimum and maximum of an isotherm.
        - The two extrema are referred to as priMinPr and priMaxPr.
        - Locations of the extrema are referred to as priMinVr and priMaxVr.
        - The search of extrema locations is executed in two steps:
              A. Delimit the locations to finite intervals; initial search
                 from the range ]VrmolMin, iSearchMax (parameter)], where
                 VrmolMin is defined by the reduced equation of state.
              B. Pinpoint the locations in the intervals (refined search).

    2. Decompose the isotherm into three parts.
        - Head: ]VrmolMin, priMinVr].
        - Body: ]priMinVr, priMaxVr].
        - Tail: ]priMaxVr, VrmolMax := 1/RhorMin (parameter)].

    3. Find the reduced equilibrium pressure, PrEq, from the body part.
        - Utilize the Maxwell's equal area rule (aka Maxwell construction).
        - Search PrEq from the range
          ]max(priMinPr, 0, Pr(VrmolMax)), priMaxPr[.

    4. Find the liquid-vapour densities associated with PrEq.
        - The intersection of PrEq and head defines the
          reduced molar volume (and density) for the liquid phase.
        - The intersection of PrEq and tail defines the
          reduced molar volume (and density) for the vapour phase.
    """
    # Exclude the lower limit for the reduced molar volume, as defined by
    # the reduced equation of state, in order to avoid division by zero
    iSearchRes = 1e-3
    VrmolMin, VrmolMax = reos.VrmolMin + iSearchRes, 1.0/RhorMin
    Vrmols = np.arange(VrmolMin, iSearchMax, iSearchRes)

    vcnt = Vrmols.shape[0]
    reosf = reos.Pr
    results = []

    for Tr in Trs:
        if Tr >= 1.0:
            raise ValueError('Reduced temperature paramater must be < 1.')

        Prs = np.fromfunction(lambda i: reosf(Vrmols[i], Tr),
                              (vcnt,), dtype=int)

        # 1.A Initial search (argrelmin and argrelmax from scipy)
        priMinArr = argrelmin(Prs)[0]
        priMaxArr = argrelmax(Prs)[0]

        if len(priMinArr) == 0 or len(priMaxArr) == 0:
            msg = 'Initial search of a primary minimum or maximum failed.\n'
            msg += 'A possible remedy: try to extend the search interval '
            msg += 'for finding primary extrema,\ni.e. increase the '
            msg += 'value of the iSearchMax parameter.'
            raise RuntimeError(msg)

        mina, minb = Vrmols[priMinArr[0]-1], Vrmols[priMinArr[0]+1]
        maxa, maxb = Vrmols[priMaxArr[0]-1], Vrmols[priMaxArr[0]+1]

        # 1.B refined search (fminbound from scipy)
        priMinVr, priMinPr, e1, nf2 = fminbound(reosf, mina, minb, (Tr,),
                                                xtol=1e-10, full_output=True)
        priMaxVr, priMaxPr, e2, nf2 = fminbound(reosf, maxa, maxb, (Tr,),
                                                xtol=1e-10, full_output=True)

        if e1 != 0 or e2 != 0:
            msg = 'Refined search of a primary minimum or maximum failed '
            msg += 'due to unknown reason.'
            raise RuntimeError(msg)

        # 2. Decompose the isotherm into three parts.
        head_ab = [VrmolMin, priMinVr]
        body_ab = [priMinVr, priMaxVr]
        tail_ab = [priMaxVr, VrmolMax]

        # 3. Find the reduced equilibrium pressure from the isotherm body
        PrVrmolMax = reosf(VrmolMax, Tr)
        priMinPr = max(priMinPr, 0, PrVrmolMax)

        if PrVrmolMax < priMinPr:
            tail_ab[1] = _finters(priMaxVr, VrmolMax, reosf, Tr, priMinPr)

        PrEq, fv, e3, nf3 = fminbound(_area_diff, priMinPr, priMaxPr,
                                      (Tr, reosf, head_ab, body_ab, tail_ab),
                                      xtol=1e-10, full_output=True)

        if e3 != 0:
            msg = 'Reduced equilibrium pressure not found, i.e. solution\n'
            msg += 'to the Maxwell\'s equal area optimization problem '
            msg += 'not converged.'
            raise RuntimeError(msg)

        # 4. Find the liquid-vapour densities associated with PrEq.
        VrmolLiq = _finters(head_ab[0], head_ab[1], reosf, Tr, PrEq)
        VrmolVap = _finters(tail_ab[0], tail_ab[1], reosf, Tr, PrEq)
        RhorVap, RhorLiq = 1.0/VrmolVap, 1.0/VrmolLiq

        results.append((Tr, PrEq, VrmolVap, VrmolLiq, RhorVap, RhorLiq))

    return results
