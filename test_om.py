from gpaw import GPAW
from gpaw.spinorbit import get_L_vlmm
import numpy as np

Lx_ii = get_L_vlmm()[0][2]
Ly_ii = get_L_vlmm()[1][2]
Lz_ii = get_L_vlmm()[2][2]

calc = GPAW('gs_soc.gpw')
#print(calc.get_occupation_numbers())

mx = 0
my = 0
mz = 0
f_km = []
for wf in calc.wfs.kpt_u:
    P_msi = wf.projections[0]
    P0_mi = P_msi[:, 0, 7:12]
    P1_mi = P_msi[:, 1, 7:12]
    f_m = wf.f_n
    f_km.append(f_m)
    A0_mm = P0_mi.conj().dot(Lx_ii).dot(P0_mi.T)
    A1_mm = P1_mi.conj().dot(Lx_ii).dot(P1_mi.T)
    mx += np.dot(f_m, np.diag(A0_mm)).real
    mx += np.dot(f_m, np.diag(A1_mm)).real

    A0_mm = P0_mi.conj().dot(Ly_ii).dot(P0_mi.T)
    A1_mm = P1_mi.conj().dot(Ly_ii).dot(P1_mi.T)
    my += np.dot(f_m, np.diag(A0_mm)).real
    my += np.dot(f_m, np.diag(A1_mm)).real

    A0_mm = P0_mi.conj().dot(Lz_ii).dot(P0_mi.T)
    A1_mm = P1_mi.conj().dot(Lz_ii).dot(P1_mi.T)
    mz += np.dot(f_m, np.diag(A0_mm)).real
    mz += np.dot(f_m, np.diag(A1_mm)).real

print(2*mx)
print(2*my)
print(2*mz)
print(np.sum(f_km, axis=0))
print(np.shape(f_km))
print(f_km[0])
