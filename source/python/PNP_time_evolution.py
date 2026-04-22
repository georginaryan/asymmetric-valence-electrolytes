import os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

def solve_pnp_time(
    zp=1.0, #cation valence
    zn=1.0, #anion valence
    eps=0.1, #small parameter
    Ip=0.2, #cation flux
    In=0.1, #anion flux
    V=1.0, #potential difference
    Nx=250, #mesh points
    t_final=50.0, #final time
):

    x = np.linspace(0, 1, Nx)
    dx = x[1] - x[0]

    # Initial conditions
    p0 = np.ones(Nx) # uniform unit cation concentration
    n0 = np.ones(Nx) # uniform unit anion concentration
    # calculate potential later from p0, n0 and BCs. 

    y0 = np.concatenate([p0, n0])
    
    ## Poisson equation as an algebraic constraint
    def solve_poisson(p, n):
        #construct Laplacian (central finite differences)
        main = -2 * np.ones(Nx)
        off = np.ones(Nx - 1)
        A = diags([off, main, off], [-1, 0, 1]) / dx**2
        
        #right hand side
        rhs = (-p + n) / eps**2
    
        A = A.tolil() # make A easily editable
        
        # set phi(0)=V
        A[0, :] = 0
        A[0, 0] = 1
        rhs[0] = V

        #set phi(1)=0
        A[-1, :] = 0
        A[-1, -1] = 1
        rhs[-1] = 0
    
        return spsolve(A.tocsr(), rhs)
        
    ## Nernst-Planck equations, to solve as an IVP
    def rhs(t, y):
        p = y[:Nx] #extract cation concentration
        n = y[Nx:2 * Nx] #extract anion concentration

        phi = solve_poisson(p, n) #solve Poisson as algebraic constraint

        dpdx = np.gradient(p, dx)
        dndx = np.gradient(n, dx)
        dphidx = np.gradient(phi, dx)

        Qp = dpdx + zp * p * dphidx #cation Nernst-Planck flux
        Qn = dndx - zn * n * dphidx #anion Nernst-Planck flux

        # Flux boundary conditions
        Qp[0] = -Ip
        Qp[-1] = -Ip
        Qn[0] = In
        Qn[-1] = In

        #Nernst-Planck continuity equations
        dpdt = np.gradient(Qp, dx)
        dndt = np.gradient(Qn, dx)

        return np.concatenate([dpdt, dndt])
    
    def export_pnp_csv(x, p, n, phi, zp, zn, Ip, In, eps):
        data = np.column_stack((x, p, n, phi))
        # Project root
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

        output_dir = os.path.join(project_root, "data")
        os.makedirs(output_dir, exist_ok=True)

        filename = (
            f"pnp_zp{zp}_zn{zn}_"
            f"Ip{Ip:.2f}_In{In:.2f}_"
            f"eps{eps:.3f}.csv"
        )

        filepath = os.path.join(output_dir, filename)

        np.savetxt(
            filepath,
            data,
            delimiter=",",
            header="x,p,n,phi",
            comments=""
        )
        
    sol = solve_ivp(
        rhs,
        (0, t_final),
        y0,
        method="BDF",
        rtol=1e-6,
        atol=1e-8,
    )
    
    p_ss = sol.y[:Nx, -1]
    n_ss = sol.y[Nx:2 * Nx, -1]
    phi_ss = solve_poisson(p_ss, n_ss) #calculate final potential from final cation and anion concentrations
    
    export_pnp_csv(x, p_ss, n_ss, phi_ss, zp, zn, Ip, In, eps)
    
    return x, p_ss, n_ss, phi_ss
