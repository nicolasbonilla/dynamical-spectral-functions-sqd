L=14 ground-state checkpoint for arXiv:2608.16436
=================================================

Paper:      N. Bonilla Vargas, "Captured weight and boundary leakage bound the error of
            sample-based spectral functions", arXiv:2608.16436 (v3).
Code:       https://github.com/nicolasbonilla/dynamical-spectral-functions-sqd

System: one-dimensional Hubbard ring, periodic boundary conditions, U/t = 8, half filling,
L = 14 sites (N-electron sector nup = nd = 7).

Files
-----
L14_gs.npz     (94,230,224 bytes)
               psi  float64[11778624]  ground-state vector of the N-electron sector
               E0   float64            ground-state energy, -4.613102625985336 t
               FAF  float64            one-body fermionic magic F1 of the ground state, 22.50497728150742
               nup, nd  int32          7, 7
               sha256 0fc2c9eb8e0e8c53f9d396582c9b22e1d0cc931cddd0c56c38b3d37950884903

L14_order.npy  (41,225,312 bytes)
               int32[10306296]  the stored ranking of the (N+1) sector (dimension 10,306,296) used for
               the L = 14 certificate point; it is verified to be a permutation of the sector, and the
               ranking protocol was not re-run at this size (Sec. S4 of the paper)
               sha256 0df36e6eb8086cbae06e52f33de34b2bef433556757dc78db97dc37d33e0737c

Use
---
In the repository both files sit in ckpt/, where the scripts look by default:

    python src/frontier/run_L14.py
    python src/frontier/z_suppK_sym.py 14

Elsewhere, set P2_CKPT to the folder that holds them.

They reproduce the L = 14 row of Table III (vacuity factors 8.61, 5.59, 4.61, 3.26) and the
L = 14 Krylov support count (10,305,696 at an amplitude threshold of 1e-12). The residual of the
stored ground state is 2.7e-9.

License: CC BY 4.0
