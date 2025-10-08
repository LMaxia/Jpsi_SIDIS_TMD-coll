#####---- Importing packages ----#######################
# System control
import os, sys
#### ---------------------------------------------- ####
# Numerical manipulation
import numpy as np
from scipy.integrate import quad
import warnings
#### ---------------------------------------------- ####
import yaml
basepath = os.path.dirname(__file__)
with open(os.path.join(basepath, 'input.yaml'), 'r') as infile:
    input_dict = yaml.safe_load(infile)
pol = input_dict['pol']
Q = input_dict['Q']; Q2 = Q**2
M = input_dict['Mjp']; M2 = M**2
ildme = input_dict['ildme']
ipdf = input_dict['ipdf']
wave = input_dict['wave']
iq = input_dict['iq']
ig = input_dict['ig']
xB = input_dict['xB']
inum = input_dict['inum']
qT_sp = input_dict['qT_sp']
qT_ep = input_dict['qT_ep']
scale_err = input_dict['scale_err']
ldme_err = input_dict['ldme_err']
#### ---------------------------------------------- ####
# Additional packages
from packages.hard_part import FUU
from packages.ldme import sel as ldme_sel
import lhapdf
########################################################


#####---- Factorization scale selection ----############
# Other options can be included here and referenced⮐
# in the input.yaml file
mu2 = (M2 + Q2)
mu2_info = '(M2 + Q2)'
#### ------------------- ####
mu = np.sqrt(mu2) 
########################################################

#####---- LDME and PDF selection ----###################
# LDME selection
# To add a set, include it in the ldme.py file and⮐
# reference it in the input.yaml file
ldme = ldme_sel(ildme)
print(f'\nLDME set: {ldme.info}')
#### ---------------------------------------------- ####
# PDF selection
# To add a set, download it from LHAPDF, and include it⮐
# here and in the input.yaml file
if ipdf == 1:
    pdf_set = lhapdf.mkPDF("NNPDF40_lo_as_01180")
    pdf_info = 'NNPDF40'
elif ipdf == 2:
    pdf_set = lhapdf.mkPDF("MSHT20lo_as130")
    pdf_info = 'MSHT20 (LO & as = 0.130)'
elif ipdf == 3:
    pdf_set = lhapdf.mkPDF("MSHT20nlo_as118")
    pdf_info = 'MSHT20 (NLO & as = 0.118)'
elif ipdf == 4:
    pdf_set = lhapdf.mkPDF("NNPDF40_nlo_as_01180")
    pdf_info = r'NNPDF40 (NLO & as = 0.118)'
elif ipdf == 5:
    pdf_set = lhapdf.mkPDF("MSTW2008lo68cl")
    pdf_info = r'MSTW2008 (LO & 68\% of c.l.)'
else:
    print('ERROR: other pdf sets not implemented')
    sys.exit()
########################################################

#####---- qT array ----#################################
qT_array = np.linspace(qT_sp, qT_ep, num=inum)
########################################################

#####---- xB function ----##############################
def xmax(qT): 
    return Q2/(M2 + Q2 + 2*M*qT)
########################################################


#####---- LHAPDF class ----#############################
# Determination of different PDF flavors from⮐
# LHAPDF database
class LHApdf():
    def __init__(self, xi, _muh2):
        ####----------------#### 
        fdo = pdf_set.xfxQ2(1, xi, _muh2)/xi
        fup = pdf_set.xfxQ2(2, xi, _muh2)/xi
        fst = pdf_set.xfxQ2(3, xi, _muh2)/xi
        fdob = pdf_set.xfxQ2(-1, xi, _muh2)/xi
        fupb = pdf_set.xfxQ2(-2, xi, _muh2)/xi
        fstb = pdf_set.xfxQ2(-3, xi, _muh2)/xi
        ####----------------#### 
        self.qqb_PDF = (fup + fdo + fst + fupb + fdob + fstb)
        self.qqb_PDF_weight = (fup + fupb + (fdo + fst + fdob + fstb)/4)
        self.gluon_PDF = pdf_set.xfxQ2(21, xi, _muh2)/xi
        ####----------------#### 
        self.als = pdf_set.alphasQ2(_muh2)  
########################################################


#####---- Function definition ----######################
def integrand(xh, qT, qT2, _muh2, err_ldme):  
    sqrt_arg = ((M2 + Q2)*xh - Q2)**2 - 4*M2*qT2*xh**2
    if sqrt_arg < 0:
        if np.abs(sqrt_arg) < 10**(-6):
            sqrt_arg = np.abs(sqrt_arg)
            print('!!Square-root argument is set positive by hand!!')
        else:
            print('!!ERROR!!')
            print(f'Square-root argument is negative at qT = {qT} GeV')
            print(f'xh = {xh}')
            print(f'xmax = {xmax(qT)}')
            print(f'xmin = {xB}')
            print(f'Sqrt argument = {sqrt_arg}')
            sys.exit()
    ####---------------------------------------#### 
    z_delta = np.sqrt(sqrt_arg)
    if z_delta < 0:
        print('Delta is negative after applying the sqrt!')
        print(f'qT: {qT} GeV')
        print(f'z_Delta: {z_delta}')
        sys.exit()
    ####---------------------------------------#### 
    # Fixed value of zh after integration over delta
    # Note that the delta depends quadratically⮐ 
    # on zh, thus two solutions are included here
    zhp = (Q2*(1 - xh) + M2*xh + z_delta)/(2*(Q2*(1 - xh) + qT2*xh)) 
    zhm = (Q2*(1 - xh) + M2*xh - z_delta)/(2*(Q2*(1 - xh) + qT2*xh)) 
    if zhp < 0 or zhp > 1:
        print("ERROR: zh plus solution is out of range!")
        print(f"zhp: {zhp}")
        sys.exit()
    if zhm < 0 or zhm > 1:
        print("ERROR: zh minus solution is out of range!")
        print(f"zhm: {zhm}")
        sys.exit()
    ####---------------------------------------#### 
    # Prefactor coming from delta evaluation over z
    facp = zhp**2*Q2/z_delta
    facm = zhm**2*Q2/z_delta   
    ####---------------------------------------#### 
    pdf = LHApdf(xB/xh, _muh2)
    ####---------------------------------------#### 
    structure_function_p = FUU(xh, zhp, pdf, ldme, err_ldme)
    structure_function_m = FUU(xh, zhm, pdf, ldme, err_ldme)
    ####---------------------------------------#### 
    return facp*structure_function_p + facm*structure_function_m
########################################################


#####--- File string ----###############################
str_fixed_values = f"Q {Q}GeV - xB {int(xB*1000)/1000}"
#### ------------------- ####
if pol == 'T':
    str_pol = 'transv'
    out_pol = 'F_UUT'
    print('\nPolarization: Transverse')
elif pol == 'L':
    str_pol = 'long'
    out_pol = 'F_UUL'
    print('\nPolarization: Longitudindal')
elif pol == 'A':
    str_pol = '2phi'
    out_pol = 'F_UU^cos(2phi)'
    print('\nPolarization: Asymmetry (2phi)')
#### ------------------- ####
if iq == 0:
    part_str = ' g'
    print('Gluon contribution')
elif ig == 0:
    part_str = ' q'
    print('Quark contribution')
else:
    part_str = ''
if wave == 0:
    mod_str = 'NRQCD (1s0wave)'
    print('Decomposition: 1s0 wave')
elif wave == 1:
    if ig == 1 and iq == 1:
        print('ERROR! 3s1 wave requires to select only one partonic channel\nPlease, select either ig=0 or iq=0')
        sys.exit()
    mod_str = 'NRQCD (3s1wave)'
    print('Decomposition: 3s1 wave')
elif wave == 2:
    mod_str = 'NRQCD (3pJwave)'
    print('Decomposition: 3pJ wave')
else:
    mod_str = f'NRQCD ({ldme.s})' 
#### ---------------------------------------------- ####
print(f'Q value: {Q} GeV\nxB value: {xB}\n')
#### ---------------------------------------------- ####
flnm_str =  f"{mod_str}{part_str} - pol {str_pol} - {str_fixed_values} - mu2 {mu2_info} - ipdf {ipdf}"
if ipdf == 3 or ipdf == 4:
    flnm_str += " (nlo)"
#### ---------------------------------------------- ####
filename = f"Jpsi unpol - ep (SIDIS) - plot vs qT - xh integr - {flnm_str}.tsv"
########################################################


#####---- Script ----###################################
if __name__ == "__main__":
    warnings.filterwarnings('ignore')   
    ####--------------------------------------------#### 
    if pol == 'T':
        pol_folder = 'transverse'
    elif pol == 'L':
        pol_folder = 'longitudinal'
    elif pol == 'A':
        pol_folder = 'asymmetries'
    os.chdir(f'{basepath}/output')
    ####--------------------------------------------#### 
    f_out = open(filename, 'w')
    f_out.write('# ep -> J/psi X\n')
    f_out.write(f'# {out_pol} full collinear\n')
    f_out.write(f'\n# ldme -> {ldme.info}\n')
    f_out.write(f'# pdf -> {pdf_info}\n')
    f_out.write('\n# Plot vs qT\n')
    f_out.write(f'# fixed variables:\t{str_fixed_values}\n')
    f_out.write(f'# Regularisation: mu2 = {mu2_info}\n')
    f_out.write('\n# qT\tF_mid\tF[mu_f/2]\tF[2*mu_f]\tF[mu_psi/2]\tF[2mu_psi]\n')
    ####--------------------------------------------#### 
    var_hard = [1/2, 2]
    var_ldme = [1/2, 2]
    ####--------------------------------------------#### 
    index = 0
    num = 0
    ####--------------------------------------------#### 
    for qT in qT_array:    
        qT2 = qT**2  
        index += 1
        num += 1
        print(f'Point {num}/{inum}')
        ####----------------------------------------####
        xh_min = xB
        xh_max = xmax(qT)
        ####----------------------------------------####
        coll_mid = quad(integrand, xh_min, xh_max, epsabs=1.e-20, limit=50000, args = (qT, qT2, mu2, 1))
        ####----------------------------------------####
        coll_var_list = ''
        if scale_err:
            coll_var_hard1 = quad(integrand, xh_min, xh_max, epsabs=1.e-20, limit=50000, args=(qT, qT2, mu2*var_hard[0]**2, 1))[0]
            coll_var_hard2 = quad(integrand, xh_min, xh_max, epsabs=1.e-20, limit=50000, args=(qT, qT2, mu2*var_hard[1]**2, 1))[0]
        else:
            coll_var_hard1 = 0
            coll_var_hard2 = 0
        if ldme_err:
            coll_var_ldme1 = quad(integrand, xh_min, xh_max, epsabs=1.e-20, limit=50000, args=(qT, qT2, mu2, var_ldme[0]))[0]
            coll_var_ldme2 = quad(integrand, xh_min, xh_max, epsabs=1.e-20, limit=50000, args=(qT, qT2, mu2, var_ldme[1]))[0]
        else:
            coll_var_ldme1 = 0
            coll_var_ldme2 = 0
        coll_var_list += f'\t{coll_var_hard1}\t{coll_var_hard2}\t{coll_var_ldme1}\t{coll_var_ldme2}' 
        ####----------------------------------------####
        if qT == qT_array[0]:
            print(f'qT value: {qT}')
            print(f'Full collinear: {coll_mid}')
        elif index == int(inum/6):
            print(f'qT value: {qT}')
            print(f'Full Collinear: {coll_mid}')
            index = 0
        ####----------------------------------------####
        f_out.write(f'{qT}\t{coll_mid[0]}{coll_var_list}\n')
    ####--------------------------------------------#### 
    f_out.close()
########################################################