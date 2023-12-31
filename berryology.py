import numpy as np
from ase.units import Bohr, Hartree, _e, _hplanck

from gpaw import GPAW
from gpaw.occupations import fermi_dirac
from gpaw.spinorbit import soc_eigenstates
from gpaw.response.pair import PairDensity
from gpaw.mpi import world, serial_comm, rank

#def get_orbital_magnetization(calc_file, kpts_k=None, 
#                              n1=0, n2=None, mi=None, mf=None, 
#                              width=None, theta=0.0, phi=0.0, scale=1.0):

#    calc1 = GPAW(calc_file, txt=None)
#    kpts_kc = calc1.get_bz_k_points()
    
#    if width is not None:
#        calc1.wfs.occupations._width = width

#    if n2 == None:
#        n2 = calc1.get_number_of_bands()
#        eigenvalues = None

#    bands_n = range(n1, n2)
#    Ns = calc1.get_number_of_spins()
#    Nn = len(bands_n)
    #if not n2 == None:
    #    eigenvalues = np.ndarray(shape=(Ns, len(calc1.get_ibz_k_points()), n2 -n1), dtype=float)
    #    for s in range(Ns):
    #        for k in range(len(calc1.get_ibz_k_points())):
    #            eigenvalues[s][k][:] = calc1.get_eigenvalues(kpt=k , spin = s)[:n2] 
                
#    soc = soc_eigenstates(calc1, scale=scale, theta=theta, phi=phi, n1=n1, n2=n2) 
#    e_km = soc.eigenvalues() / Hartree
#    efermi = soc.fermi_level / Hartree
#    f_km = soc.occupation_numbers()

#    calc2 = GPAW(calc_file, txt=None, communicator=serial_comm)
#    pdensity = PairDensity(calc2, txt=None, world=serial_comm)

#    M_x_k = np.zeros(len(kpts_k))
#    M_y_k = np.zeros(len(kpts_k))
#    M_z_k = np.zeros(len(kpts_k))

#    threshold = 1

#    del calc1
#    for ik, wfs in soc.wfs.items():
#        v_msn = wfs.v_msn
#        k_c = kpts_kc[ik]

#        if mf is not None:
#            f_m[:mf+1-2*n1] = 1
#            f_m[mf+1-2*n1:] = 0
#        if mi is not None:
#            f_m[:mi-2*n1] = 0

#        f_m = f_km[ik]
#        f_mm = (f_m[:, np.newaxis] - f_m[:])
        
#        aeps_mm = (e_km[ik,:, np.newaxis] + e_km[ik, :] - 2*efermi)
#        deps_mm = e_km[ik,:, np.newaxis] - e_km[ik, :]
#        deps_mm[deps_mm == 0.0] = np.inf

#        smallness_mm = np.abs(-1e-3 / deps_mm)
#        inds_mm = (np.logical_and(np.inf > smallness_mm,
#                               smallness_mm > threshold))

#        frac_mm = aeps_mm*(f_mm / np.square(deps_mm))
#        frac_mm[inds_mm] = 0
#        frac_mm = np.nan_to_num(frac_mm)

#        rho_snnv = np.zeros((2, Nn, Nn, 3), complex)
#        for s in range(Ns):
#            kpt = pdensity.get_k_point(s, k_c, n1, n2)
#            for n in bands_n:
#                rho_snnv[s, n] = pdensity.optical_pair_velocity(n, kpt, kpt)

#        if Ns == 1:
#            rho_snnv[1] = rho_snnv[0]
        
#        rho_mmv = np.dot(v_msn[:, 0].conj(),                                       
#                           np.dot(v_msn[:, 0], rho_snnv[0]))                             
#        rho_mmv += np.dot(v_msn[:, 1].conj(),                                       
#                          np.dot(v_msn[:, 1], rho_snnv[1]))                             

#        A_mm = rho_mmv[:, :, 1].conj() * rho_mmv[:, :, 2]
#        M_x_k[ik] = np.einsum('mn, mn', A_mm, frac_mm).imag
        
#        A_mm = rho_mmv[:, :, 2].conj() * rho_mmv[:, :, 0]
#        M_y_k[ik] = np.einsum('mn, mn', A_mm, frac_mm).imag

#        A_mm = rho_mmv[:, :, 0].conj() * rho_mmv[:, :, 1]
#        M_z_k[ik] = np.einsum('mn, mn', A_mm, frac_mm).imag


#    world.sum(M_x_k)
#    world.sum(M_y_k)
#    world.sum(M_z_k)

#    M_x_k = 0.5*M_x_k
#    M_y_k = 0.5*M_y_k
#    M_z_k = 0.5*M_z_k

#    return M_x_k, M_y_k, M_z_k

def get_berrycurvature(calc_file,
                       n1=0, n2=None, mi=None, mf=None,
                       theta=0.0, phi=0.0, scale=1.0, orb_mag=False):

    calc1 = GPAW(calc_file, txt=None)
    kpts_kc = calc1.get_bz_k_points()

    if n2 is None:
        n2 = calc1.get_number_of_bands()
    bands_n = range(n1, n2)
    Nn = len(bands_n)
    Ns = calc1.get_number_of_spins()

    #if not n2 == None:
    #    eigenvalues = np.ndarray(shape=(Ns, len(calc1.get_ibz_k_points()), n2 -n1), dtype=float)
    #    for s in range(Ns):
    #        for k in range(len(calc1.get_ibz_k_points())):
    #            eigenvalues[s][k][:] = calc1.get_eigenvalues(kpt=k , spin = s)[:n2] 

    soc = soc_eigenstates(calc1, scale=scale, theta=theta, phi=phi, n1=n1, n2=n2)
    e_km = soc.eigenvalues() / Hartree
    f_km = soc.occupation_numbers()
    if orb_mag:
        efermi = soc.fermi_level / Hartree
    
    calc2 = GPAW(calc_file, txt=None, communicator=serial_comm)
    pdensity = PairDensity(calc2, txt=None, world=serial_comm)

    berryc_xy_k = np.zeros(len(kpts_kc))
    berryc_yz_k = np.zeros(len(kpts_kc))
    berryc_zx_k = np.zeros(len(kpts_kc))

    threshold = 1
    
    del calc1
    for ik, wfs in soc.wfs.items():
        v_msn = wfs.v_msn
        k_c = kpts_kc[ik]

        if mf is not None:
            f_m[:mf+1-2*n1] = 1
            f_m[mf+1-2*n1:] = 0
        if mi is not None:
            f_m[:mi-2*n1] = 0

        f_m = f_km[ik]
        f_mm = (f_m[:, np.newaxis] - f_m[:])
        
        #deps_mm = e_km[ik,:, np.newaxis] - e_km[ik, :]
        #deps_mm[deps_mm == 0.0] = np.inf

        #smallness_mm = np.abs(-1e-3 / deps_mm)
        #inds_mm = (np.logical_and(np.inf > smallness_mm,
        #                          smallness_mm > threshold))

        if orb_mag:
            deps_mm = e_km[ik, :, np.newaxis] - e_km[ik, :]
            aeps_mm = (e_km[ik, :, np.newaxis] + e_km[ik, :] - 2 * efermi)
        else:
            deps_mm = e_km[ik, :, np.newaxis] - e_km[ik, :]
            aeps_mm = 1

        smallness_mm = np.abs(1e-3 / deps_mm)
        inds_mm = (np.logical_and(np.inf > smallness_mm,
                                  smallness_mm > threshold))

        #print(np.square(deps_mm))
        frac_mm = aeps_mm*f_mm / np.square(deps_mm)
        frac_mm[inds_mm] = 0
        frac_mm = np.nan_to_num(frac_mm)

        rho_snnv = np.zeros((2, Nn, Nn, 3), complex)
        for s in range(Ns):
            kpt = pdensity.get_k_point(s, k_c, n1, n2)
            #for n in bands_n:
            for n in range(0, Nn):
                rho_snnv[s, n] = pdensity.optical_pair_velocity(n, kpt, kpt)

        if Ns == 1:
            rho_snnv[1] = rho_snnv[0]
        
        rho_mmv = np.dot(v_msn[:, 0].conj(),                                       
                           np.dot(v_msn[:, 0], rho_snnv[0]))                             
        rho_mmv += np.dot(v_msn[:, 1].conj(),                                       
                          np.dot(v_msn[:, 1], rho_snnv[1]))                             

        A_mm = rho_mmv[:, :, 0].conj() * rho_mmv[:, :, 1]
        berryc_xy_k[ik] = np.einsum('mn, mn', A_mm, frac_mm).imag
                
        A_mm = rho_mmv[:, :, 1].conj() * rho_mmv[:, :, 2]
        berryc_yz_k[ik] = np.einsum('mn, mn', A_mm, frac_mm).imag
        
        A_mm = rho_mmv[:, :, 2].conj() * rho_mmv[:, :, 0]
        berryc_zx_k[ik] = np.einsum('mn, mn', A_mm, frac_mm).imag
        
    soc.kpt_comm.sum(berryc_xy_k)
    soc.kpt_comm.sum(berryc_yz_k)
    soc.kpt_comm.sum(berryc_zx_k)

    if not orb_mag:
        berryc_xy_k = berryc_xy_k * Bohr**2
        berryc_yz_k = berryc_yz_k * Bohr**2
        berryc_zx_k = berryc_zx_k * Bohr**2
    if orb_mag:
        berryc_xy_k = berryc_xy_k * 0.5
        berryc_yz_k = berryc_yz_k * 0.5
        berryc_zx_k = berryc_zx_k * 0.5

    #berryc_xy_k = berryc_xy_k * Bohr**2
    #berryc_yz_k = berryc_yz_k * Bohr**2
    #berryc_zx_k = berryc_zx_k * Bohr**2

    return berryc_xy_k, berryc_yz_k, berryc_zx_k

def get_hall_conductivity(*args, **kwargs):
    berryc_xy_k, berryc_yz_k, berryc_zx_k = get_berrycurvature(*args, **kwargs, orb_mag=False)

    sigma_xy = (_e**2/_hplanck)*berryc_xy_k
    sigma_yz = (_e**2/_hplanck)*berryc_yz_k
    sigma_zx = (_e**2/_hplanck)*berryc_zx_k

    return sigma_xy, sigma_yz, sigma_zx


def get_orbital_magnetization(*args, **kwargs):
    M_x, M_y, M_z = get_berrycurvature(*args,
                                       **kwargs,
                                       orb_mag=True)
    
    return M_x, M_y, M_z
