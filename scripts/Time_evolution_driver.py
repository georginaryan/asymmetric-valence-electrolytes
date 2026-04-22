import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from source.python.PNP_time_evolution import solve_pnp_time

#For Figure 2
zn = 1
Ip = 2
In = 1
eps = 0.05
V = 1.0
x, p11, n11, phi11 = solve_pnp_time(1, zn, eps, Ip, In, V)
x, p21, n21, phi21 = solve_pnp_time(2, zn, eps, Ip, In, V)
x, p31, n31, phi31 = solve_pnp_time(3, zn, eps, Ip, In, V)

#For Figure 3
zp=2
zn = 1
Ip = 1
In = 1
V = 1.0
x, p11, n11, phi11 = solve_pnp_time(zp, zn, 0.2, Ip, In, V)
x, p21, n21, phi21 = solve_pnp_time(zp, zn, 0.01, Ip, In, V)