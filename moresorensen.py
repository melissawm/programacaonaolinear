from scipy.linalg import norm, cholesky, solve
import numpy as np

def ms(g,H,delta):
    """
    Método de Moré e Sorensen para a solução do subproblema de região de confiança na norma 2.
    """
    theta = 0.0001
    eps_D = 0.1
    n = g.shape[0]
    value = 0
    gnorm = norm(g)
    g2D = gnorm/delta
    Hnorminf = norm(H,np.inf)
    HnormF = norm(H,'fro')
    lower = np.maximum(0, g2D - np.minimum(Hnorminf,HnormF))
    upper = np.maximum(0, g2D + np.minimum(Hnorminf,HnormF))
    Dlower = (1-eps_D)*delta
    Dupper = (1+eps_D)*delta
    
    if lower < 1.0e-12:
        lbd = 0.0
    else:
        lbd = np.maximum(np.sqrt(lower*upper),lower+theta*(upper-lower))

    identity = np.eye(n)
    i = 0
    while i < 101:
        new_lbd = -1
        Hlambda = H+lbd*identity
        try:
            # Cholesky returns a lower triangular factor L such that Hlambda = L*L.T
            R = cholesky(Hlambda)
            isPosDef = True
        except:
            isPosDef = False

        if isPosDef:
            s = - solve(R, solve(R.T,g))
            norms = norm(s)
            if ((lbd < 1.0e-12 and norms <= Dupper) or (norms >= Dlower and norms<= Dupper)):
                # Then we are finished. s is the solution.
                updates = i
                i = 101
                if (norms >= Dlower and norms<= Dupper):
                    text = "(Solution found in the boundary of the trust region.)"
                else:
                    text = "(Solution found in the interior of the trust region.)"

            w = solve(R.T,s)
            normw2 = norm(w)**2
            new_lbd = lbd + ((norms-delta)/delta)*(norms**2/normw2)
            if norms > Dupper:
                lower = lbd
            else:
                upper = lbd
            
            theta_range = theta * ( upper - lower )
            if new_lbd > lower + theta_range and new_lbd < upper - theta_range:
                lbd = new_lbd
            else:
                lbd = np.maximum(np.sqrt(lower*upper),lower+theta_range)
        else:
            lower = lbd
            lbd = np.maximum(np.sqrt(lower*upper),lower+theta*(upper-lower))

        i = i + 1
            
    if i == 100:
        print("Error: 100 iterations in MS.")

    return s

if __name__ == "__main__":
    s = ms(np.asarray([-1., 1.]), np.asarray([[1., 0.], [0., 1.]]), 1.)