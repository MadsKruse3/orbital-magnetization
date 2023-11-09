import pytest
from gpaw.response.berryology import get_orbital_magnetization

### Slow test uses gpw file with 64 kpts.. ### 
#@pytest.mark.response
#def test_ahc_fe_bcc(in_tmp_dir, gpw_files):
    #calc_file = gpw_files['fe_pw_wfs']

    ## Check results with default options (with the magnetization pointing along the z-direction). 
    #sigma_xy_w_so, sigma_yz_w_so, sigma_zx_w_so = get_hall_conductivity(calc_file)

    #print(sigma_xy_w_so)
    #print(sigma_yz_w_so)
    #print(sigma_zx_w_so)

    # ------------------ Check results without spin-orbit coupling. -----------------------#
    #sigma_xy_wo_so, sigma_yz_wo_so, sigma_zx_wo_so = get_hall_conductivity(calc_file, scale=0.0)

    
    # ------------------ Check results when M is oriented along the x-axis. -----------------------#
    #sigma_xy_M_x, sigma_yz_M_x, sigma_zx_M_x = get_hall_conductivity(calc_file, theta=90, phi=0)

    # ------------------ Check results when M is oriented along the y-axis. -----------------------#
    #sigma_xy_M_y, sigma_yz_M_y, sigma_zx_M_z = get_hall_conductivity(calc_file, theta=90, phi=90)

    #----------------- Test results when the contribution from 3 bands (occupied and unoccupied) are included in matrix elements etc. ---------------------#
    #n1 = 1
    #n2 = 4
    #sigma_xy_wo_so, sigma_yz_wo_so, sigma_zx_wo_so = get_hall_conductivity(calc_file, n1=n1, n2=n2)

    #----------------- Test results when occupation numbers are changed such that only 2 bands are occupied ---------------------#
    #mi = 2
    #mf = 3
    #sigma_xy_wo_so, sigma_yz_wo_so, sigma_zx_wo_so = get_hall_conductivity(calc_file, mi=mi, mf=mf)


from gpaw import GPAW, PW, FermiDirac
from ase.build import bulk
from ase.units import Bohr
@pytest.mark.response
def test_fe_bcc(in_tmp_dir):

    a = 5.42
    a = a*Bohr
    atoms = bulk('Fe', 'bcc', a=a)
    atoms.set_initial_magnetic_moments([3.75])

    xc = 'LDA'
    kpts = 2
    occw = 0.01
    nbands = 6
    pw = 200
    conv = {'density': 1e-3,
           'forces': 1e-2}
    
    calc = GPAW(xc=xc,
                mode=PW(pw),
                kpts={'size': (kpts, kpts, kpts), 'gamma': True},
                occupations=FermiDirac(occw),
                convergence=conv,
                nbands=nbands)


    atoms.calc = calc
    atoms.get_potential_energy()
    calc.write('gs_Fe_test.gpw', mode='all')

    calc_file = 'gs_Fe_test.gpw'

    # ------------------ Check results with default parameters. -----------------------#
    #test_M_x_def, test_M_y_def, test_M_z_def = get_orbital_magnetization(calc_file)

    #M_x_def = [-2.76987700e-32, 6.77625564e-10, 1.26921221e-03, -5.44642998e-10, -5.44642870e-10, -1.26921221e-03, -6.77625612e-10, -2.22312060e-03]
    #M_y_def = [7.60665920e-31, -6.77625595e-10, 5.44642799e-10, -1.26921221e-03, 1.26921221e-03, 5.44642954e-10, -6.77625576e-10, 2.22312057e-03]
    #M_z_def = [1.37139630e-30, -3.18938286e-03, -1.21292465e-03, -1.21292465e-03, -1.21292465e-03, -1.21292465e-03, -3.18938286e-03, 6.09181302e-04]

    #assert test_M_x_def == pytest.approx(M_x_def, rel=1e-12, abs=1e-6)
    #assert test_M_y_def == pytest.approx(M_y_def, rel=1e-12, abs=1e-6)
    #assert test_M_z_def == pytest.approx(M_z_def, rel=1e-12, abs=1e-6)

    #print(test_M_x_def)
    #print(test_M_y_def)
    #print(test_M_z_def)

    # ------------------ Check results with reduced spin-orbit coupling. -----------------------#
    #test_M_x_wr_so, test_M_y_wr_so, test_M_z_wr_so = get_orbital_magnetization(calc_file, scale=0.5)

    #M_x_wr_so = [3.46095953e-32, -6.14684541e-11, 6.31021165e-04, -2.97312911e-10, -2.97313249e-10, -6.31021165e-04, 6.14683833e-11, -1.13285376e-03]
    #M_y_wr_so = [2.47346006e-32, 6.14684503e-11, 2.97313019e-10, -6.31021165e-04, 6.31021165e-04, 2.97312691e-10, 6.14683171e-11, 1.13285375e-03]
    #M_z_wr_so = [-8.02610117e-32, -1.45348416e-03, -6.16515199e-04, -6.16515199e-04, -6.16515199e-04, -6.16515199e-04, -1.45348416e-03, 2.89134048e-04]

    #assert test_M_x_wr_so == pytest.approx(M_x_wr_so, rel=1e-12, abs=1e-6)
    #assert test_M_y_wr_so == pytest.approx(M_y_wr_so, rel=1e-12, abs=1e-6)
    #assert test_M_z_wr_so == pytest.approx(M_z_wr_so, rel=1e-12, abs=1e-6)

#    print(test_M_x_wr_so)
#    print(test_M_y_wr_so)
#    print(test_M_z_wr_so)

    # ------------------ Check results when spins are oriented along the x-axis. -----------------------#
    #test_M_x_s_x, test_M_y_s_x, test_M_z_s_x = get_orbital_magnetization(calc_file, theta=90, phi=0)

    #M_x_s_x = [ 2.23625528e-30, -1.21292465e-03, -1.21292465e-03, -3.18938286e-03, -3.18938286e-03, -1.21292465e-03, -1.21292465e-03, 6.09181248e-04]
    #M_y_s_x = [ 9.93132815e-31, 1.26921221e-03, -5.44642805e-10, 6.77625600e-10, -6.77625563e-10, 5.44642852e-10, -1.26921221e-03, -2.22312060e-03]
    #M_z_s_x = [-4.67493368e-32, -5.44642993e-10, 1.26921221e-03, 6.77625530e-10, 6.77625527e-10, -1.26921221e-03, 5.44642990e-10, -2.22312060e-03]

    #assert test_M_x_s_x == pytest.approx(M_x_s_x, rel=1e-12, abs=1e-6)
    #assert test_M_y_s_x == pytest.approx(M_y_s_x, rel=1e-12, abs=1e-6)
    #assert test_M_z_s_x == pytest.approx(M_z_s_x, rel=1e-12, abs=1e-6)

#    print(test_M_x_s_x)
#    print(test_M_y_s_x)
#    print(test_M_z_s_x)

    # ------------------ Check results when M is oriented along the y-axis. -----------------------#
    #test_M_x_s_y, test_M_y_s_y, test_M_z_s_y = get_orbital_magnetization(calc_file, theta=90, phi=90)


    #M_x_s_y = [1.81873705e-30, 1.26921221e-03, 6.77625647e-10, -5.44642845e-10, 5.44642964e-10, -6.77625521e-10, -1.26921221e-03, -2.22312060e-03]
    #M_y_s_y = [9.95981597e-31, -1.21292465e-03, -3.18938286e-03, -1.21292465e-03, -1.21292465e-03, -3.18938286e-03, -1.21292465e-03, 6.09181302e-04]
    #M_z_s_y = [1.71447934e-30, 5.44642787e-10, -6.77625637e-10, -1.26921221e-03, 1.26921221e-03, -6.77625553e-10, 5.44642878e-10, 2.22312057e-03]
    
    #assert test_M_x_s_y == pytest.approx(M_x_s_y, rel=1e-12, abs=1e-6)
    #assert test_M_y_s_y == pytest.approx(M_y_s_y, rel=1e-12, abs=1e-6)
    #assert test_M_z_s_y == pytest.approx(M_z_s_y, rel=1e-12, abs=1e-6)

    #print(test_M_x_s_y)
    #print(test_M_y_s_y)
    #print(test_M_z_s_y)

    #----------------- Test results when the contribution from 3 bands (occupied and unoccupied) are included in matrix elements etc. ---------------------

    n1 = 1
    n2 = 6
    test_M_x_n6, test_M_y_n6, test_M_z_n6 = get_orbital_magnetization(calc_file, n1=n1, n2=n2)

    M_x_n6 = [-7.04288559e-45, 1.35633762e-04, -1.10562561e-01, -3.15288233e-03, -3.15288804e-03, 1.11595571e-01, -1.11511764e-04, -7.87927469e-03]
    M_y_n6 = [3.92252677e-46, 1.11464347e-04, 3.15288804e-03, 1.11595571e-01, -1.10562561e-01, 3.15288233e-03, 1.35586345e-04, 7.84404289e-03]
    M_z_n6 = [6.74127924e-46, -1.17532148e-01, -1.15124694e-01, -1.16196726e-01, -1.15124694e-01, -1.16196726e-01, -1.17532143e-01, 2.30705817e-03]

    assert test_M_x_n6 == pytest.approx(M_x_n6, rel=1e-12, abs=1e-3)
    assert test_M_y_n6 == pytest.approx(M_y_n6, rel=1e-12, abs=1e-3)
    assert test_M_z_n6 == pytest.approx(M_z_n6, rel=1e-12, abs=1e-3)

    #print(test_M_x_n6)
    #print(test_M_y_n6)
    #print(test_M_z_n6)
    
    #----------------- Test results when occupation numbers are changed such that only 2 bands are occupied ---------------------#
    mi = 2
    mf = 3
    test_M_x_bandocc2, test_M_y_bandocc2, test_M_z_bandocc2 = get_orbital_magnetization(calc_file, mi=mi, mf=mf)

    M_x_bandocc2 = [6.03860123e-28, -4.41656305e-06, 6.75622496e-02, -3.75404894e-06, -3.75404879e-06, -6.75622496e-02, 4.41656305e-06, 2.67070952e-05]
    M_y_bandocc2 = [-7.29449523e-28, 4.41656304e-06, 3.75404879e-06, -6.75622496e-02, 6.75622496e-02, 3.75404893e-06, 4.41656305e-06, -2.67070951e-05]
    M_z_bandocc2 = [8.79333319e-29, 1.33821069e+00, -2.60610423e-02, -2.60610423e-02, -2.60610423e-02, -2.60610423e-02, 1.33821069e+00, 1.10350807e-03]

    assert test_M_x_bandocc2 == pytest.approx(M_x_bandocc2, rel=1e-12, abs=1e-6)
    assert test_M_y_bandocc2 == pytest.approx(M_y_bandocc2, rel=1e-12, abs=1e-6)
    assert test_M_z_bandocc2 == pytest.approx(M_z_bandocc2, rel=1e-12, abs=1e-6)

#    print(test_M_x_bandocc2)
#    print(test_M_y_bandocc2)
#    print(test_M_z_bandocc2)


#@pytest.mark.response
#def test_ahc_default_params_bcc_Fe(in_tmp_dir):
#    calc_file = 'gs_Fe_test.gpw'
    # ------------------ Check results with default parameters. -----------------------#
#    test_sigma_xy_def, test_sigma_yz_def, test_sigma_zx_def = get_hall_conductivity(calc_file)

#    sigma_xy_def = [-7.03325011e-32, 2.55683148e-07, 3.47873186e-07, 3.47873186e-07, 3.47873186e-07, 3.47873186e-07,  2.55683148e-07, 6.52075901e-08]
#    sigma_yz_def = [6.95671942e-32, -4.70688597e-12, -2.72226471e-07, -6.85077200e-12, -6.85077190e-12,  2.72226471e-07,  4.70688593e-12, -1.04109688e-07]
#    sigma_zx_def = [5.01550705e-32,  4.70688597e-12,  6.85077195e-12,  2.72226471e-07, -2.72226471e-07,  6.85077200e-12,  4.70688588e-12,  1.04109687e-07]

#    assert test_sigma_xy_def == pytest.approx(sigma_xy_w_so, rel=1e-12, abs=1e-6)
#    assert test_sigma_yz_def == pytest.approx(sigma_yz_w_so, rel=1e-12, abs=1e-6)
#    assert test_sigma_zx_def == pytest.approx(sigma_zx_w_so, rel=1e-12, abs=1e-6)


#@pytest.mark.response
#def test_ahc_wr_so_bcc_Fe(in_tmp_dir):
#    calc_file = 'gs_Fe_test.gpw'

    # ------------------ Check results with reduced spin-orbit coupling. -----------------------#
#    test_sigma_xy_wr_so, test_sigma_yz_wr_so, test_sigma_zx_wr_so = get_hall_conductivity(calc_file, scale=0.5)
    
#    sigma_xy_wr_so = [2.65945571e-35, 9.09123982e-21, -1.51837123e-21, -8.20896458e-21, -2.39734970e-21, -3.15988374e-21, -3.30625648e-21, 5.45657831e-18]
#    sigma_yz_wr_so = [-2.34745867e-35, -1.22830874e-21, 1.50917611e-21, 4.12587328e-22, -1.15692559e-20, -3.16713755e-21, 1.21494356e-20, -2.41266576e-18]
#    sigma_zx_wr_so = [1.09714122e-34, 1.23861601e-21, -1.69421491e-20, -8.21831398e-21, 2.38723898e-21, 5.37227383e-21, 1.21423014e-20, -7.90097766e-18]

#    assert test_sigma_xy_wr_so == pytest.approx(sigma_xy_wr_so, rel=1e-12, abs=1e-6)
#    assert test_sigma_yz_wr_so == pytest.approx(sigma_yz_wr_so, rel=1e-12, abs=1e-6)
#    assert test_sigma_zx_wr_so == pytest.approx(sigma_zx_wr_so, rel=1e-12, abs=1e-6)


#@pytest.mark.response
#def test_ahc_M_x_bcc_Fe(in_tmp_dir):
#    calc_file = 'gs_Fe_test.gpw'

    # ------------------ Check results when M is oriented along the x-axis. -----------------------#
#    test_sigma_xy_M_x, test_sigma_yz_M_x, test_sigma_zx_M_x = get_hall_conductivity(calc_file, theta=90, phi=0)

#    sigma_xy_M_x = [2.02056944e-32, -6.85077191e-12, -2.72226471e-07, -4.70688596e-12, -4.70688592e-12, 2.72226471e-07, 6.85077194e-12, -1.04109688e-07]
#    sigma_yz_M_x = [-4.69041620e-33, 3.47873186e-07, 3.47873186e-07, 2.55683148e-07, 2.55683148e-07, 3.47873186e-07, 3.47873186e-07, 6.52075888e-08]
#    sigma_zx_M_x = [2.11868104e-32, -2.72226471e-07, -6.85077201e-12, -4.70688590e-12, 4.70688592e-12, 6.85077202e-12, 2.72226471e-07, -1.04109688e-07]
    

#    assert test_sigma_xy_M_x == pytest.approx(sigma_xy_M_x, rel=1e-12, abs=1e-6)
#    assert test_sigma_yz_M_x == pytest.approx(sigma_yz_M_x, rel=1e-12, abs=1e-6)
#    assert test_sigma_zx_M_x == pytest.approx(sigma_zx_M_x, rel=1e-12, abs=1e-6)


#@pytest.mark.response
#def test_ahc_M_y_bcc_Fe(in_tmp_dir):
#    calc_file = 'gs_Fe_test.gpw'

    # ------------------ Check results when M is oriented along the y-axis. -----------------------#
#    test_sigma_xy_M_y, test_sigma_yz_M_y, test_sigma_zx_M_y = get_hall_conductivity(calc_file, theta=90, phi=90)

#    sigma_xy_M_y = [8.42146901e-32, 6.85077202e-12, 4.70688589e-12, 2.72226471e-07, -2.72226471e-07, 4.70688593e-12, 6.85077193e-12, 1.04109687e-07]
#    sigma_yz_M_y = [-3.20258099e-32, -2.72226471e-07, -4.70688585e-12, -6.85077203e-12, 6.85077200e-12, 4.70688586e-12, 2.72226471e-07, -1.04109688e-07]
#    sigma_zx_M_y = [-7.88077284e-32, 3.47873186e-07, 2.55683148e-07, 3.47873186e-07, 3.47873186e-07, 2.55683148e-07, 3.47873186e-07, 6.52075902e-08]
    
#    assert test_sigma_xy_M_y == pytest.approx(sigma_xy_M_y, rel=1e-12, abs=1e-6)
#    assert test_sigma_yz_M_y == pytest.approx(sigma_yz_M_y, rel=1e-12, abs=1e-6)
#    assert test_sigma_zx_M_y == pytest.approx(sigma_zx_M_y, rel=1e-12, abs=1e-6)


#@pytest.mark.response
#def test_ahc_fewbands_bcc_Fe(in_tmp_dir):
#    calc_file = 'gs_Fe_test.gpw'

    #----------------- Test results when the contribution from 3 bands (occupied and unoccupied) are included in matrix elements etc. ---------------------#
#    n1 = 1
#    n2 = 6
#    test_sigma_xy_n6, test_sigma_yz_n6, test_sigma_zx_n6 = get_hall_conductivity(calc_file, n1=n1, n2=n2)

#    sigma_xy_n6 = [-1.57755951e-49, 3.83301771e-04, 2.98900587e-04, 3.02384790e-04, 2.98900587e-04, 3.02384790e-04, 3.83301756e-04, 9.56402602e-07]
#    sigma_yz_n6 = [ 1.82606166e-48,-5.54517452e-07, 3.01298109e-04, 1.21287525e-06, 1.21289425e-06, -3.04773305e-04, 5.68102614e-07, -4.30542259e-07]
#    sigma_zx_n6 = [-1.03075846e-49, -5.67835032e-07, -1.21289425e-06, -3.04773305e-04, 3.01298109e-04, -1.21287525e-06, -5.54249870e-07, 4.17036218e-07]
    
#    assert test_sigma_xy_n6 == pytest.approx(sigma_xy_n6, rel=1e-12, abs=1e-5)
#    assert test_sigma_yz_n6 == pytest.approx(sigma_yz_n6, rel=1e-12, abs=1e-5)
#    assert test_sigma_zx_n6 == pytest.approx(sigma_zx_n6, rel=1e-12, abs=1e-5)


#@pytest.mark.response
#def test_ahc_few_occ_bands_bcc_Fe(in_tmp_dir):
#    calc_file = 'gs_Fe_test.gpw'

    #----------------- Test results when occupation numbers are changed such that only 2 bands are occupied ---------------------#
#    mi = 2
#    mf = 3
#    test_sigma_xy_bandocc2, test_sigma_yz_bandocc2, test_sigma_zx_bandocc2 = get_hall_conductivity(calc_file, mi=mi, mf=mf)

#    sigma_xy_bandocc2 = [-9.42778557e-33, -3.48313049e-04,  6.83259931e-06,  6.83259931e-06, 6.83259931e-06,  6.83259931e-06, -3.48313049e-04,  8.26705496e-08]
#    sigma_yz_bandocc2 = [-6.45456155e-32,  1.14912373e-09, -1.78199461e-05,  9.88020731e-10, 9.88020692e-10,  1.78199461e-05, -1.14912373e-09,  2.52643923e-09]
#    sigma_zx_bandocc2 = [ 7.79251365e-32, -1.14912372e-09, -9.88020693e-10,  1.78199461e-05, -1.78199461e-05, -9.88020729e-10, -1.14912373e-09, -2.52643925e-09]

#    assert test_sigma_xy_bandocc2 == pytest.approx(sigma_xy_bandocc2, rel=1e-12, abs=1e-6)
#    assert test_sigma_yz_bandocc2 == pytest.approx(sigma_yz_bandocc2, rel=1e-12, abs=1e-6)
#    assert test_sigma_zx_bandocc2 == pytest.approx(sigma_zx_bandocc2, rel=1e-12, abs=1e-6)







