from gpaw.response.pair import PairDensity, KPointPair
from gpaw import GPAW
import numpy as np
from gpaw.spinorbit import BZWaveFunctions, soc_eigenstates
from ase.dft.kpoints import bandpath
from ase.units import Bohr, Hartree, kJ, _me, _e, _hplanck
from gpaw.mpi import world, serial_comm, broadcast, rank

def parallelisation_sizes(kpts_k, rank=None):
    if rank is None:
        rank = world.rank
    nK = len(kpts_k) 
    myKsize = -(-nK // world.size)
    myKrange = range(rank * myKsize,
                     min((rank + 1) * myKsize, nK))
    myKsize = len(myKrange)
    return myKrange, myKsize


def get_orbital_magnetization(calc, kpts_k=None, n1=0, n2=None, theta=0.0, phi=0.0, scale=1.0):
    kpts_kc = calc.get_ibz_k_points()
    if kpts_k == None:
        kpts_k = range(len(kpts_kc))

    if n2 == None:
        n2 = calc.get_number_of_bands()
    bands_n = range(n1, n2)

    pdensity = PairDensity(calc, txt=None, communicator=serial_comm)

    myKrange, myKsize = parallelisation_sizes(kpts_k)

    soc = soc_eigenstates(calc, scale=scale, theta=theta, phi=phi, n1=n1, n2=n2) 
    wfs_k = soc.wfs
    v_kmsn = soc.eigenvectors()
    e_km = soc.eigenvalues() / Hartree
    Ns = calc.get_number_of_spins()
    Nn = len(bands_n)
    chempot = calc.get_chemical_potential()
    
    berryc_xy_k = np.zeros(len(kpts_k))
    berryc_yz_k = np.zeros(len(kpts_k))
    berryc_zx_k = np.zeros(len(kpts_k))

    threshold = 1

    for i in myKrange:
        ik = kpts_k[i]
        k_c = kpts_kc[ik]

        f_m = wfs_k[ik].f_m
        f_mm = (f_m[:, np.newaxis] - f_m[:])
        
        aeps_mm = e_km[ik,:, np.newaxis] + e_km[ik, :]
        deps_mm = e_km[ik,:, np.newaxis] - e_km[ik, :]
        deps_mm[deps_mm == 0.0] = np.inf

        smallness_mm = np.abs(-1e-3 / deps_mm)
        inds_mm = (np.logical_and(np.inf > smallness_mm,
                                  smallness_mm > threshold))

        frac_mm = aeps_mm*(f_mm / np.square(deps_mm))
        frac_mm[inds_mm] = 0
        frac_mm = np.nan_to_num(frac_mm)

        rho_snnv = np.zeros((2, Nn, Nn, 3), complex)
        for s in range(Ns):
            e0_n = calc.get_eigenvalues(kpt=ik, spin=s)[bands_n] / Hartree
            deps0_nn = e0_n[:, np.newaxis] - e0_n[:]
            kpt = pdensity.get_k_point(s, k_c, n1, n2)
            for n in bands_n:
                #rho_snnv[s, n] = pdensity.optical_pair_velocity(n, bands_n, kpt, kpt)
                rho_snnv[s, n] = pdensity.optical_pair_velocity(n, kpt, kpt)

        if Ns == 1:
            rho_snnv[1] = rho_snnv[0]
        
        v_msn = v_kmsn[ik]
        rho_mmv = np.dot(v_msn[:, 0].conj(),                                       
                           np.dot(v_msn[:, 0], rho_snnv[0]))                             
        rho_mmv += np.dot(v_msn[:, 1].conj(),                                       
                          np.dot(v_msn[:, 1], rho_snnv[1]))                             

        
        A_mm = rho_mmv[:, :, 0].conj() * rho_mmv[:, :, 1]
        berryc_xy_k[i] = np.einsum('mn, mn', A_mm, frac_mm).imag
        
        A_mm = rho_mmv[:, :, 1].conj() * rho_mmv[:, :, 2]
        berryc_yz_k[i] = np.einsum('mn, mn', A_mm, frac_mm).imag
        
        A_mm = rho_mmv[:, :, 2].conj() * rho_mmv[:, :, 0]
        berryc_zx_k[i] = np.einsum('mn, mn', A_mm, frac_mm).imag

    berryc_xy_k = berryc_xy_k * Bohr**2
    berryc_yz_k = berryc_yz_k * Bohr**2
    berryc_zx_k = berryc_zx_k * Bohr**2

    world.sum(berryc_xy_k)
    world.sum(berryc_yz_k)
    world.sum(berryc_zx_k)

    ## Units of Bohr magnetons
    #prefactor = _me/(_hplanck/(2*np.pi))**2
    ## Si units to eV and Ångstrom
    #m_to_Ångstrom = 1e10
    #J_to_eV = kJ/1000
    #prefactor = prefactor/((m_to_Ångstrom**2)*J_to_eV)
    #N = len(atom)

    ## Convert
    #M_x = Omega_yz*prefactor
    #M_y = Omega_zx*prefactor
    #M_z = Omega_xy*prefactor

    M_x_k = (_e/2*_hplanck)*berryc_yz_k
    M_y_k = (_e/2*_hplanck)*berryc_zx_k
    M_z_k = (_e/2*_hplanck)*berryc_xy_k

    return M_x, M_y, M_z
