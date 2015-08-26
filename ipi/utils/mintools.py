"""Contains simple helper algorithms for minimization.

Copyright (C) 2013, Joshua More and Michele Ceriotti

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http.//www.gnu.org/licenses/>.

Algorithms implemented by Michele Ceriotti and Benjamin Helfrecht, 2015

Functions: 
        bracket: Determines the 3 points that bracket the function minimum
        min_brent:  Does one-D minimization (line search) based on bisection 
            method with derivatives. Uses 'bracket' function.
        min_approx: Does approximate n-D minimization (line search) based 
            on sufficient function decrease in the search direction
        BFGS: Constructs an approximate inverse Hessian to determine 
            new search directions. Minimizes the function using 'min_approx' function.
        L-BFGS: Uses the limited memory BFGS algorithm (L-BFGS) to 
            compute new search directions. Minimizes using 'min_approx'
        L-BFGS_nls: L-BFGS algorithm without line search
            *** This function is less stable than L-BFGS and not any more efficient ***
        bracket_neb: Modified 'bracket' routine to make 
            compatible with functions with unknown gradient
        min_brent_neb: Modified 'min_brent' routine to make 
            compatible with functions with unknown gradient

        bracket, bracket_neb, min_brent, min_brent_neb,and BFGS subroutines adapted from: 
            Press, W. H., Teukolsky, S. A., Vetterling, W. T., and Flannery, B. P. (1992). 
            Numerical Recipes in C: The Art of Scientific Computing. 
            Cambridge: Cambridge University Press

        LBFGS subroutine adapted from:
            Nocedal, J. (1980). Updating Quasi-Newton Matrices with
            Limited Storage. Mathematics of Computation, 35, 773-782.
            DOI: http://dx.doi.org/10.1090/S0025-5718-1980-0572855-7
"""

#TODO: CLEAN UP BFGS, L-BFGS, L-BFGS_nls TO NOT EXIT WITHIN MINTOOLS.PY BUT USE UNIVERSAL SOFTEXIT

__all__ = [ "min_brent" ]

import numpy as np
import math
from ipi.utils.messages import verbosity, warning, info

# Bracketing function
def bracket(fdf, fdf0=None, x0=0.0, init_step=1.0e-3): 

    """Given an initial point, determines the initial bracket for the minimum
     Arguments:
            fdf: function to minimize, derivative of function to minimize
            x0: initial point
            fdf0: value of function and its derivative at x0
            init_step: intial step size
    """

    # Constants
    gold = 1.618034 # Golden ratio
    glimit = 100.0 # Limit for magnification of parabolic fit step
    tiny = 1.0e-20 # Prevent division by zero

    # x0: initial point of evaluation (e.g. initial atomic position)
    # ax, bx, cx: bracketing points with ax < bx < cx
    # fa, fb, fc: value of function at ax, bx, cx
    if fdf0 is None: fdf0 = fdf(x0)
    ax = x0 
    fa, dfa = fdf0 
    bx = x0 + init_step
    fb, dfb = fdf(bx)
    info(" @BRACKET: Started bracketing", verbosity.debug)
    info(" @BRACKET: Evaluated first step", verbosity.debug)

    # Switch direction to move downhill, if necessary, and rearrange
    if fb > fa:
        tmp = ax
        ax = bx
        bx = tmp
        tmp = fb
        fb = fa
        fa = tmp
        tmp = dfb
        dfb = dfa
        dfa = tmp

    # Initial guess for third bracketing point
    cx = bx + gold * (bx - ax)
    fc, dfc = fdf(cx)
    info(" @BRACKET: Evaluated initial bracket: (%f:%f, %f:%f, %f:%f)" % (ax, fa, bx, fb, cx, fc), verbosity.debug) 

    # Loop until acceptable bracketing condition is achieved
    # u is a point between two of the bracketing points
    # Use parabolic extrapolation to find u. "tiny" prevents possible division by zero
    while fb > fc:
        r = (bx - ax) * (fb - fc)
        q = (bx - cx) * (fb - fa)
        u = bx - ((bx - cx) * q - (bx - ax) * r) / (2.0 * math.copysign(max(abs(q - r), tiny), (q - r))) # Point from parabolic fit
        ulim = bx + glimit * (cx - bx) # Limit for parabolic fit point; *Can test various possibilities*

        # Find minimums between b and c or a and u
        # If parabolic fit unsuccessful, use default step magnification
        # Otherwise:
        # - Parabolic fit between c and its allowed limit
        # - Limit u to maximum allowed value
        # - Use default magnification
        if ((bx - u) * (u - cx)) > 0.0:
            fu, dfu = fdf(u)
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug) 
            if fu < fc:
                ax = bx
                bx = u
                fa = fb
                fb = fu
                dfa = dfb
                dfb = dfu
                info(" @BRACKET: Bracketing completed: (%f:%f, %f:%f, %f:%f)" % (ax, fa, bx, fb, cx, fc), verbosity.debug)
                return (ax, bx, cx, fb, dfb)

            elif fu > fb:
                cx = u
                fc = fu
                dfc = dfu
                info(" @BRACKET: Bracketing completed", verbosity.debug)
                return (ax, bx, cx, fb, dfb)

            u = cx + gold * (cx - bx)
            fu, dfu = fdf(u)
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug) 
        elif ((cx - u) * (u - ulim)) > 0.0:
            fu, dfu = fdf(u)
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug) 
            if fu < fc:
                bx = cx
                cx = u
                u = cx + gold * (cx - bx)
                fb = fc
                fc = fu
                dfb = dfc
                dfc = dfu
                fu, dfu = fdf(u)
                info(" @BRACKET: Evaluated new bracket point", verbosity.debug) 
        elif ((u - ulim) * (ulim - cx)) >= 0.0:
            u = ulim
            fu, dfu = fdf(u)
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug) 
        else:
            u = cx + gold * (cx - bx)
            fu, dfu = fdf(u)
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug) 

        # Shift points
        ax = bx
        bx = cx
        cx = u
        fa = fb
        fb = fc
        fc = fu
        dfa = dfb
        dfb = dfc
        dfc = dfu

    info(" @BRACKET: Bracketing completed: (%f:%f, %f:%f, %f:%f)" % (ax, fa, bx, fb, cx, fc), verbosity.debug)
    return (ax, bx, cx, fb, dfb)

# One dimensional minimization function using function derivatives
# and Brent's method
def min_brent(fdf, fdf0=None, x0=0.0, tol=1.0e-6, itmax=100, init_step=1.0e-3):

    """Given a maximum number of iterations and a convergence tolerance,
     minimizes the specified function 
     Arguments:
            x0: initial x-value
            fdf: function to minimize
            fdf0: initial function value
            tol: convergence tolerance
            itmax: maximum allowed iterations
            init_step: initial step size
    """

    # Initializations and constants
    gold = 0.3819660 # Golden ratio
    zeps = 1.0e-10 # Safeguard against trying to find fractional precision for min that is exactly zero
    e = 0.0 # Size of step before last

    # Call initial bracketing routine
    (ax, bx, cx, fb, dfb) = bracket(fdf, fdf0, x0, init_step)
    
    # Set bracket points
    if ax < cx:
        a = ax
    else:
        a = cx
    if ax > cx:
        b = ax
    else: 
        b = cx

    # Initial points to evaluate
    # f* is evaluation of function
    # df* is the evaluation of the derivative
    # x = point with least function value so far
    # w = point with 2nd least function value
    # v = previous value of w
    # u = point at which function was evaluated most recently
    x = w = v = bx
    fw = fv = fx = fb # Function
    dfw = dfv = dfx = dfb # Function derivative
    
    # Main loop
    j = 1
    info(" @MINIMIZE: Started 1D minimization", verbosity.debug)
    while j <= itmax:

        # Determine tolerance
        xm = 0.5 * (a + b)
        tol1 = tol * abs(x) + zeps
        tol2 = 2.0 * tol1

        # Test for satisfactory completion
        if abs(x - xm) <= (tol2 - 0.5 * (b - a)):
            info(" @MINIMIZE: Finished minimization, energy = %f" % fx, verbosity.debug)
            return (x, fx)

        # Initialize d values (used to determine step size) to outside of bracket
        if abs(e) > tol1:
            d1 = 2.0 * (b - a)
            d2 = d1

            # Secant method with both d points
            if dfw != dfx:
                d1 = (w - x) * dfx / (dfx - dfw)
            if dfv != dfx:
                d2 = (v - x) * dfx / (dfx - dfv)

            # Choose estimate based on derivative at x and move distance on step
            # before last
            u1 = x + d1
            u2 = x + d2
            ok1 = ((a - u1) * (u1 - b) > 0.0) and (dfx * d1 <= 0.0)
            ok2 = ((a - u2) * (u2 - b) > 0.0) and (dfx * d2 <= 0.0)
            olde = e
            e = d

            # Take an acceptable d; if both are acceptable, choose smallest
            if ok1 or ok2:
                if ok1 and ok2:
                    if abs(d1) < abs(d2):
                        d = d1
                    else:
                        d = d2
                elif ok1:
                    d = d1
                else:
                    d = d2
                if abs (d) <= abs (0.5 * olde):
                    u = x + d
                    if ((u - a) < tol2) or ((b - u) < tol2):
                         d = math.copysign(tol1, (xm - x))
                    else:
                        if dfx >= 0.0:
                            e = a - x
                        else:
                            e = b - x
                        d = 0.5 * e
            else:
                if dfx >= 0.0:
                    e = a - x
                else:
                    e = b - x
                d = 0.5 * e
        else:
            if dfx >= 0.0:
                e = a - x
            else:
                e = b - x
            d = 0.5 * e
        if abs(d) >= tol1:
            u = x + d
            fu, dfu = fdf(u)
        else:
            u = x + math.copysign(tol1, d)
            fu, dfu = fdf(u)
         
            # If minimum step goes uphill, minimum has been found
            if fu > fx:
                info(" @MINIMIZE: Finished minimization, energy = %f" % fx, verbosity.debug)
                return (x, fx)
                
        if fu <= fx:
            if u >= x:
                a = x
            else:
                b = x
            v = w
            fv = fw
            dfv = dfw
            w = x
            fw = fx
            dfw = dfx
            x = u
            fx = fu
            dfx = dfu
        else:
            if u < x:
                a = u
            else:
                b = u
            if (fu <= fw) or (w == x):
                v = w
                fv = fw
                dfv = dfw
                w = u
                fw = fu
                dfw = dfu
            elif (fu < fv) or (v == x) or (v == w):
                v = u
                fv = fu
                dfv = dfu
        j += 1
    
    # Exit if maximum number of iterations exceeded
    info(" @MINIMIZE: Error -- maximum iterations for minimization (%d) exceeded, exiting minimization" % itmax, verbosity.low)
    info(" @MINIMIZE: Finished minimization, energy = %f" % fx, verbosity.debug)
    return (x, fx)

# Approximate line search
def min_approx(fdf, x0, fdf0=None, d0=None, max_step=100.0, tol=1.0e-6, itmax=100):
    
    """Given an n-dimensional function and its gradient, and an 
    initial point and a direction, finds a new point where the function
    is thought to be 'sufficiently' minimized, i.e. carries out an 
    approximate minimization. 
        Arguments:
            fdf: function and its gradient
            fdf0: initial function and gradient value
            d0: n-dimensional initial direction
            x0: n-dimensional initial point
            max_step: maximum step size
            tol: tolerance for exiting line search
            itmax: maximum number of iterations for the line search
    """
    
    # Initializations and constants
    info(" @MINIMIZE: Started approx. line search", verbosity.debug) 
    n = len(x0.flatten())
    if fdf0 is None: fdf0 = fdf(x0)
    f0, df0 = fdf0
    if d0 is None: d0 = -df0 / np.sqrt(np.dot(df0.flatten(), df0.flatten()))
    x = np.zeros(n)
    alf = 1.0e-4

    # Step size
    stepsum = np.sqrt(np.dot(d0.flatten(), d0.flatten()))

    # Scale if attempted step is too large
    if stepsum > max_step:
        info(" @MINIMIZE: Scaled step size for line search", verbosity.debug)
        d0 = np.multiply(d0, max_step / stepsum)

    slope = np.dot(df0.flatten(), d0.flatten())

    if slope >= 0.0:
        info(" @MINIMIZE: Warning -- gradient is >= 0 (%f)" % slope, verbosity.low)

    test = np.amax(np.divide(np.absolute(d0.flatten()), np.maximum(np.absolute(x0.flatten()), np.ones(n))))

    # Setup to try Newton step first
    alamin = tol / test
    alam = 1.0

    # Minimization Loop
    i = 1
    while i < itmax:
        x = np.add(x0, (alam * d0))
        fx, dfx = fdf(x)
        info(" @MINIMIZE: Calculated energy", verbosity.debug)
        
        # Check for convergence on change in x
        if alam < alamin:
            x = x0
            info(" @MINIMIZE: Convergence in position, exited line search", verbosity.debug)
            return (x, fx, dfx)

        # Sufficient function decrease
        elif fx <= (f0 + alf * alam * slope):
            info(" @MINIMIZE: Sufficient function decrease, exited line search", verbosity.debug)
            return (x, fx, dfx)

        # No convergence; backtrack
        else:
            info(" @MINIMIZE: No convergence on step; backtrack to find point", verbosity.debug)

            # First backtrack
            if alam == 1.0:
                tmplam = -slope / (2.0 * (fx - f0 - slope))
            
            # Subsequent backtracks
            else:
                rhs1 = fx - f0 - alam * slope
                rhs2 = f2 - f0 - alam2 * slope
                a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2)
                b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2)
                if a == 0.0:
                    tmplam = -slope / (2.0 * b)
                
                else:
                    disc = b * b - 3.0 * a * slope
                    if disc < 0.0:
                        tmplam = 0.5 * alam

                    elif b <= 0.0:
                        tmplam = (-b + np.sqrt(disc)) / (3.0 * a)

                    else:
                        tmplam = -slope / (b + np.sqrt(disc))

                    # Coefficient less than 0.5 * lambda_1
                    if tmplam > (0.5 * alam):
                        tmplam = 0.5 * alam

        alam2 = alam
        f2 = fx

        # Coefficient greater than 0.1 * lambda_1
        alam = max(tmplam, 0.1 * alam)

        i += 1

    info(" @MINIMIZE: Error - maximum iterations for line search (%d) exceeded, exiting search" % itmax, verbosity.low)
    info(" @MINIMIZE: Finished minimization, energy = %f" % fx, verbosity.debug)
    return (x, fx, dfx)
        
# BFGS algorithm with approximate line search
def BFGS(x0, d0, fdf, fdf0=None, invhessian=None, max_step=100, tol=1.0e-6, itmax=100):
    
    """BFGS minimization. Uses approximate line minimizations.
    Does one step.
        Arguments:
            fdf: function and gradient
            fdf0: initial function and gradient value
            d0: initial direction for line minimization
            x0: initial point
            max_step: limit on step length
            tol: convergence tolerance
            itmax: maximum number of allowed iterations
    """
    
    # Original function value, gradient, other initializations
    zeps = 1.0e-10
    if fdf0 is None: fdf0 = fdf(x0)
    f0, df0 = fdf0
    n = len(x0.flatten())
    if invhessian is None: invhessian = np.eye(n)
    dg = np.zeros(n)
    g = df0.flatten()
    hdg = np.zeros(n)
    x = np.zeros(n)
    linesum = np.dot(x0.flatten(), x0.flatten())
    
    # Initial line direction
    xi = d0

    # Maximum step size
    max_step = max_step * max(np.sqrt(linesum), n)

    # Perform approximate line minimization in direction d0
    x, fx, dfx = min_approx(fdf, x0, fdf0, xi, max_step, tol, itmax) 

    info(" @MINIMIZE: Started BFGS", verbosity.debug)

    # Update line direction (xi) and current point (x0)
    xi = np.subtract(x, x0).flatten()
    x0 = x

    # Store old gradient
    dg = g

    # Get new gradient      
    g = dfx
    info(" @MINIMIZE: Updated gradient", verbosity.debug)
    g = g.flatten()

    # Compute difference of gradients
    dg = np.subtract(g, dg)

    # Difference of gradients times current matrix
    hdg = np.dot(invhessian, dg)

    fac = np.dot(dg.flatten(), xi.flatten())
    fae = np.dot(dg.flatten(), hdg.flatten())
    sumdg = np.dot(dg.flatten(), dg.flatten())
    sumxi = np.dot(xi.flatten(), xi.flatten())

    # Skip update if not 'fac' sufficiently positive
    if fac > np.sqrt(zeps * sumdg * sumxi):
        fac = 1.0 / fac
        fad = 1.0 / fae

        # Compute BFGS term
        dg = np.subtract(fac * xi, fad * hdg)

        invhessian = invhessian + np.outer(xi, xi) * fac - np.outer(hdg, hdg) * fad + np.outer(dg, dg) * fae        
        info(" @MINIMIZE: Updated hessian", verbosity.debug)
    else:
        info(" @MINIMIZE: Skipped hessian update; direction x gradient insufficient", verbosity.debug)
    
    # Update direction
    xi = np.dot(invhessian, -g)
    info(" @MINIMIZE: Updated search direction", verbosity.debug)
    return (x, fx, xi, invhessian)

# L-BFGS algorithm with approximate line search
def L_BFGS(x0, d0, fdf, qlist, glist, fdf0=None, max_step=100, tol=1.0e-6, itmax=100, m=0, k=0):
    
    """L-BFGS minimization. Uses approximate line minimizations.
    Does one step.
        Arguments:
            fdf = function and gradient
            fdf0 = initial function and gradient value
            d0 = initial direction for line minimization
            x0 = initial point
            qlist = list of previous positions used for reduced inverse Hessian construction
            glist = list of previous gradients used for reduced inverse Hessian construction
            m = number of corrections to store and use
            k = iteration (MD step) number
            max_step = limit on step length
            tol = convergence tolerance
            itmax = maximum number of allowed iterations
    """
    
    # Original function value, gradient, other initializations
    zeps = 1.0e-10
    if fdf0 is None: fdf0 = fdf(x0)
    f0, df0 = fdf0
    n = len(x0.flatten())
    dg = np.zeros(n)
    g = df0
    x = np.zeros(n)
    linesum = np.dot(x0.flatten(), x0.flatten())
    alpha = np.zeros(m)
    beta = np.zeros(m)
    rho = np.zeros(m)
    q = np.zeros(n)
    
    # Initial line direction
    xi = d0

    # Maximum step size
    max_step = max_step * max(np.sqrt(linesum), n)

    # Perform approximate line minimization in direction d0
    x, fx, dfx = min_approx(fdf, x0, fdf0, xi, max_step, tol, itmax)

    info(" @MINIMIZE: Started L-BFGS", verbosity.debug)

    # Update line direction (xi) and current point (x0)
    xi = np.subtract(x, x0)

    # Build list of previous positions
    if k < m:
        qlist[k] = xi.flatten()
    else:
        qlist = np.roll(qlist, -1, axis=0)
        qlist[m - 1] = xi.flatten()
    
    # Update current point
    x0 = x

    # Store old gradient
    dg = g

    # Get new gradient      
    g = dfx
    info(" @MINIMIZE: Updated gradient", verbosity.debug)

    # Compute difference of gradients
    q = g.flatten()
    dg = np.subtract(g, dg)

    # Build list of previous gradients
    if k < m:
        glist[k] = dg.flatten()
    else:
        glist = np.roll(glist, -1, axis=0)
        glist[m - 1] = dg.flatten()

    fac = np.dot(dg.flatten(), xi.flatten())
    sumdg = np.dot(dg.flatten(), dg.flatten())
    sumxi = np.dot(xi.flatten(), xi.flatten())

    # Determine bounds for L-BFGS 'two loop recursion'
    if k < (m - 1):
        bound1 = k
        bound2 = k + 1
    else:
        bound1 = m - 1
        bound2 = m
        
    # Skip update if not 'fac' sufficiently positive
    if fac > np.sqrt(zeps * sumdg * sumxi):
        
        # Begin two loop recursion:
        # First loop
        for j in range(bound1, -1,  -1):
            rho[j] = 1.0 / np.dot(glist[j], qlist[j])
            alpha[j] = rho[j] * np.dot(qlist[j], q)
            q = q - alpha[j] * glist[j]

        info(" @MINIMIZE: First L-BFGS loop recursion completed", verbosity.debug)

        # Two possiblities for scaling: using first or most recent 
        # members of the gradient and position lists
        hk = np.dot(glist[bound1], qlist[bound1]) / np.dot(glist[bound1], glist[bound1])
        #hk = np.dot(glist[0], qlist[0]) / np.dot(glist[0], glist[0])
        xi = hk * q

        # Second loop
        for j in range(0, bound2, 1):
            beta[j] = rho[j] * np.dot(glist[j], xi)
            xi = xi + qlist[j] * (alpha[j] - beta[j])

        info(" @MINIMIZE: Second L-BFGS loop recursion completed", verbosity.debug)

    else:
        info(" @MINIMIZE: Skipped direction update; direction * gradient insufficient", verbosity.debug)

    # Update direction xi #TODO: MOVE THIS TO OUTSIDE OF IF/ELSE SO RUNS EVERY TIME
    xi = -1.0 * xi.reshape(d0.shape)
    info(" @MINIMIZE: Updated search direction", verbosity.debug)
    return (x, fx, xi, qlist, glist)

# Bracketing for NEB, TODO: DEBUG THIS IF USING SD OR CG OPTIONS FOR NEB
def bracket_neb(fdf, fdf0=None, x0=0.0, init_step=1.0e-3): 

    """Given an initial point, determines the initial bracket for the minimum
     Arguments:
            fdf: function to minimize 
            x0: initial point
            fdf0: value of function at x0
            init_step: initial step size
    """

    # Constants
    gold = 1.618034 # Golden ratio
    glimit = 100.0 # Limit for magnification of parabolic fit step
    tiny = 1.0e-20 # Prevent division by zero

    # x0: initial point of evaluation (e.g. initial atomic position)
    # ax, bx, cx: bracketing points with ax < bx < cx
    # fa, fb, fc: value of function at ax, bx, cx
    if fdf0 is None: fdf0 = fdf(x0)
    ax = x0 
    fa = fdf0 
    bx = x0 + init_step
    fb = fdf(bx)[1]
    info(" @BRACKET: Started bracketing", verbosity.debug)
    info(" @BRACKET: Evaluated first step", verbosity.debug)

    # Switch direction to move downhill, if necessary
    if fb > fa:
        tmp = ax
        ax = bx
        bx = tmp
        tmp = fb
        fb = fa
        fa = tmp

    # Initial guess for third bracketing point
    cx = bx + gold * (bx - ax)
    fc = fdf(cx)[1]
    info(" @BRACKET: Evaluated initial bracket: (%f:%f, %f:%f, %f:%f)" % (ax, fa, bx, fb, cx, fc), verbosity.debug) 

    # Loop until acceptable bracketing condition is achieved
    # u is a point between two of the bracketing points
    # Use parabolic extrapolation to find u. "tiny" prevents possible division by zero
    while fb > fc:
        r = (bx - ax) * (fb - fc)
        q = (bx - cx) * (fb - fa)
        u = bx - ((bx - cx) * q - (bx - ax) * r) / (2.0 * math.copysign(max(abs(q - r), tiny), (q - r))) # Point from parabolic fit
        ulim = bx + glimit * (cx - bx) # Limit for parabolic fit point; *Can test various possibilities*

        # Find minimums between b and c or a and u
        # If parabolic fit unsuccessful, use default step magnification
        # Otherwise:
        # - Parabolic fit between c and its allowed limit
        # - Limit u to maximum allowed value
        # - Use default magnification
        if ((bx - u) * (u - cx)) > 0.0:
            fu = fdf(u)[1]
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug) 
            if fu < fc:
                ax = bx
                bx = u
                fa = fb
                fb = fu
                info(" @BRACKET: Bracketing completed: (%f:%f, %f:%f, %f:%f)" % (ax, fa, bx, fb, cx, fc), verbosity.debug)
                return (ax, bx, cx, fb)

            elif fu > fb:
                cx = u
                fc = fu
                info(" @BRACKET: Bracketing completed", verbosity.debug)
                return (ax, bx, cx, fb)

            u = cx + gold * (cx - bx)
            fu = fdf(u)[1]
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug) 
        elif ((cx - u) * (u - ulim)) > 0.0:
            fu = fdf(u)[1]
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug) 
            if fu < fc:
                bx = cx
                cx = u
                u = cx + gold * (cx - bx)
                fb = fc
                fc = fu
                fu = fdf(u)[1]
                info(" @BRACKET: Evaluated new bracket point", verbosity.debug) 
        elif ((u - ulim) * (ulim - cx)) >= 0.0:
            u = ulim
            fu = fdf(u)[1]
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug) 
        else:
            u = cx + gold * (cx - bx)
            fu = fdf(u)[1]
            info(" @BRACKET: Evaluated new bracket point", verbosity.debug) 

        # Shift points
        ax = bx
        bx = cx
        cx = u
        fa = fb
        fb = fc
        fc = fu

    info(" @BRACKET: Bracketing completed: (%f:%f, %f:%f, %f:%f)" % (ax, fa, bx, fb, cx, fc), verbosity.debug)
    return (ax, bx, cx, fb)


# Minimize using only forces; for NEB
def min_brent_neb(fdf, fdf0=None, x0=0.0, tol=1.0e-6, itmax=100, init_step=1.0e-3):

    """Given a maximum number of iterations and a convergence tolerance,
     minimizes the specified function 
     Arguments:
            x0: initial x-value
            fdf: function to minimize
            fdf0: initial function value
            tol: convergence tolerance
            itmax: maximum allowed iterations
            init_step: initial step size
    """

    # Initializations and constants
    gold = 0.3819660
    zeps = 1e-10    
    e = 0.0 # Step size for step before last

    (ax, bx, cx, fb) = bracket_neb(fdf, fdf0, x0, init_step)
    
    # Set bracket points
    if ax < cx:
        a = ax
    else:
        a = cx
    if ax > cx:
        b = ax
    else: 
        b = cx

    # Initial points to evaluate
    # f* is evaluation of arbitrary function
    # x = point with least function value so far
    # w = point with 2nd least function value
    # v = previous value of w
    # u = point at which function was evaluated most recently
    # d = used to determine step size
    x = w = v = bx
    fw = fv = fx = fb 

    # Main loop
    j = 1
    while j <= itmax:

        # Determine tolerance
        xm = 0.5 * (a + b)
        tol1 = tol * abs(x) + zeps
        tol2 = 2.0 * tol1

        # Test for satisfactory completion
        if abs(x - xm) <= (tol2 - 0.5 * (b - a)):
            xmin = x
            return xmin, fx

        # Complete an iteration if error is greater than tolerance
        # and construct parabolic fit from parameters
        if abs(e) > tol1:
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2.0 * (q - r)
            if q > 0.0:
                p = -p
            q = abs(q)
            etmp = e
            e = d

            # Determine acceptability of parabolic fit
            if (abs(p) >= abs(0.5 * q * etmp)) or (p <= (q * (a-x))) or (p >= (q * (b-x))):
                if x >= xm:
                    e = a - x
                else:
                    e = b - x 
                d = gold * e

            # Take parabolic step
            else:
                d = p / q
                u = x + d
                if ((u - a) < tol2) or ((b - u) < tol2):
                    d = abs(tol1) * (xm - x) / abs(xm - x)
        else:
            if x < xm:
                e = a - x
            else:
                e = b - x
            d = gold * e
        if abs(d) >= tol1:
            u = x + d
        else:
            u = x + abs(tol1) * d / abs(d)

        fu = fdf(u)[1]

        if fu <= fx:
            if (u >= x):
                a = x
            else:
                b = x
            
            # Reassign bracketing points
            v = w
            w = x
            x = u
            fv = fw
            fw = fx
            fx = fu
        else:
            if u < x:
                a = u
            else:
                b = u
            if fu <= fw or w == x:
                v = w
                w = u
                fv = fw
                fw = fu
            elif (fu <= fv) or (v == x) or (v == w):
                v = u
                fv = fu
        j += 1
    
    # Exit if maximum number of iterations exceeded
    xmin = x
    return xmin, fx 

# L-BFGS without line search; WARNING: UNSTABLE
def L_BFGS_nls(x0, d0, fdf, qlist, glist, fdf0=None, max_step=100, tol=1.0e-6, itmax=100, init_step=1.0e-3, m=0, k=0):
    
    """L-BFGS minimization without line search
    Does one step.
        Arguments:
            fdf: function and gradient
            fdf0: initial function and gradient value
            d0: initial direction for line minimization
            x0: initial point
            qlist: list of previous positions used for reduced inverse Hessian construction
            glist: list of previous gradients used for reduced inverse Hessian construction
            m: number of corrections to store and use
            k: iteration (MD step) number
            max_step: limit on step length
            tol: convergence tolerance
            itmax: maximum number of allowed iterations
            init_step: initial step size
    """
    
    # Original function value, gradient, other initializations
    zeps = 1.0e-10
    if fdf0 is None: fdf0 = fdf(x0)
    f0, df0 = fdf0
    n = len(x0.flatten())
    dg = np.zeros(n)
    g = df0 
    x = np.zeros(n)
    linesum = np.dot(x0.flatten(), x0.flatten())
    alpha = np.zeros(m)
    beta = np.zeros(m)
    rho = np.zeros(m)
    q = np.zeros(n)
    
    # Initial line direction
    xi = d0
    dg = df0 

    # Step size
    stepsize = np.sqrt(np.dot(d0.flatten(), d0.flatten()))

    # First iteration; use initial step
    if k == 0:
        scale = 1.0
        while np.sqrt(np.dot(g.flatten(), g.flatten())) >= np.sqrt(np.dot(df0.flatten(), df0.flatten()))\
                or np.isnan(np.sqrt(np.dot(g.flatten(), g.flatten()))) == True\
                or np.isinf(np.sqrt(np.dot(g.flatten(), g.flatten()))) == True:
            x = np.add(x0, (scale * init_step * d0 / np.sqrt(np.dot(d0.flatten(), d0.flatten()))))
            scale *= 0.1
            fx, g = fdf(x)
    else:

        # Scale if attempted step is too large
        if stepsize > max_step:
            d0 = max_step * d0 / np.sqrt(np.dot(d0.flatten(), d0.flatten()))
            info(" @MINIMIZE: Scaled step size", verbosity.debug)

        x = np.add(x0, d0)       
        print "step size:", np.sqrt(np.dot(d0.flatten(), d0.flatten()))
        fx, g = fdf(x)

    info(" @MINIMIZE: Started L-BFGS", verbosity.debug)
    info(" @MINIMIZE: Updated gradient", verbosity.debug)

    # Update line direction (xi) and current point (x0)
    xi = np.subtract(x, x0)

    # Build list of previous positions
    if k < m:
        qlist[k] = xi.flatten()
    else:
        qlist = np.roll(qlist, -1, axis=0)
        qlist[m - 1] = xi.flatten()
    
    # Update current point
    x0 = x

    # Compute difference of gradients
    q = g.flatten()
    dg = np.subtract(g, dg)

    # Build list of previous gradients
    if k < m:
        glist[k] = dg.flatten()
    else:
        glist = np.roll(glist, -1, axis=0)
        glist[m - 1] = dg.flatten()

    # Determine bounds for L-BFGS 'two loop recursion'
    if k < (m - 1):
        bound1 = k
        bound2 = k + 1
    else:
        bound1 = m - 1
        bound2 = m
        
    # Begin two loop recursion:
    # First loop
    for j in range(bound1, -1,  -1):
        rho[j] = 1.0 / np.dot(glist[j], qlist[j])
        alpha[j] = rho[j] * np.dot(qlist[j], q)
        q = q - alpha[j] * glist[j]

    info(" @MINIMIZE: First L-BFGS loop recursion completed", verbosity.debug)

    # Two possiblities for scaling: using first or most recent 
    # members of the gradient and position lists
    hk = np.dot(glist[bound1], qlist[bound1]) / np.dot(glist[bound1], glist[bound1])
    #hk = np.dot(glist[0], qlist[0]) / np.dot(glist[0], glist[0])
    xi = hk * q

    # Second loop
    for j in range(0, bound2, 1):
        beta[j] = rho[j] * np.dot(glist[j], xi)
        xi = xi + qlist[j] * (alpha[j] - beta[j])

    # Update direction xi
    xi = -xi.reshape(d0.shape)

    info(" @MINIMIZE: Second L-BFGS loop recursion completed", verbosity.debug)
    info(" @MINIMIZE: Updated search direction", verbosity.debug)

    return (x, fx, xi, qlist, glist)

