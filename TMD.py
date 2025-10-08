#####---- Importing packages ----#######################
# System control
import sys
import os
#### ---------------------------------------------- ####
# Numerical manipulation
import numpy as np
import scipy.special as sp
from scipy.integrate import quad
import itertools
import warnings
#### ---------------------------------------------- ####
from scipy.optimize import root
#### ---------------------------------------------- ####
import yaml
basepath = os.path.dirname(__file__)
with open(os.path.join(basepath, 'input.yaml'), 'r') as infile:
    input_dict = yaml.safe_load(infile)
e_c = input_dict['e_c']
alpha_QED = input_dict['alpha_QED']
pol = input_dict['pol']
Q = input_dict['Q']; Q2 = Q**2
M = input_dict['Mjp']; M2 = M**2
Mz = input_dict['Mz']
ildme = input_dict['ildme']
ipdf = input_dict['ipdf']
wave = input_dict['wave']
iq = input_dict['iq']
ig = input_dict['ig']
xB = input_dict['xB']
qT_sp = input_dict['qT_sp']
qT_ep = input_dict['qT_ep']
inum = input_dict['inum']
acc = input_dict['acc']
iBCO = input_dict['iBCO']
gpsi = input_dict['gpsi']
scale_err = input_dict['scale_err']
ldme_err = input_dict['ldme_err']
np_err = input_dict['np_err']
npoints = input_dict['npoints']
#### ---------------------------------------------- ####
# Additional packages
from packages.sigma_factor import small_qT_expansion as sqT_exp
from packages.ldme import sel as ldme_sel
import lhapdf
########################################################


#####---- Factorization scale selection ----############
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


#####---- TMD variables ----############################
b0 = 1.12292 
bT_max = 1.5
nf = 5
CA = 3
beta0 = (11*CA - 2*nf)/3
beta1 = (102*CA - 38*nf)/3
#### ---------------------------------------------- ####
# numerical determination of Lambda_QCD
def nsolve(func, x0, method='hybr', tol=1e-10):
    solution = root(func, x0, method=method, tol=tol)
    if solution.success:
        return solution.x[0]
    else:
        raise ValueError(f"Solver failed: {solution.message}")
#### ---------------------------------------------- ####
def alpha_Mz_1loop(Lambda):
    return pdf_set.alphasQ(Mz) - 4*np.pi/(beta0*np.log(Mz**2/Lambda**2))
LQCD_1loop_guess = 0.16675389136254443
LQCD_1loop = nsolve(alpha_Mz_1loop, LQCD_1loop_guess)
print(f'\nLambda_QCD from 1-loop evolution: {LQCD_1loop}')
#### ---------------------------------------------- ####
if acc == '1-loop' or acc == 'NLL':
    LQCD = LQCD_1loop
else:
    print('Uncorrect accuracy selected.\nPlease choose between [1-loop] or [NLL]\n')
    sys.exit()
print(f'alpha_s(M_z): {pdf_set.alphasQ(Mz)}')
#### ---------------------------------------------- ####
# S_NP by Aybat and Rogers
g1 = 0.184
g2 = 0.201
g3 = -0.129
xc = 0.009
A = [.05, .414, 0.8]            
def B(_x): 
    # 9/4 due to Casimir scaling (CA/CF)
    # 0.5 since this BNP is constructed from the Drell-Yan nonperturbative Sudakov with 2 protons
    return 0.5*(9/4)*(g2*(1 + 2*g3*np.log((10*_x*xc)/(xc + _x))) - g1*np.log(2))
########################################################
    

#####---- TMD perturbative tails ----###################
def x(qT):
    xmax = Q2/(M2 + Q2 + 2*M*qT)
    return xB/xmax
####---------------------------#### 
# Only gluon TMD contributes and determined from LHAPDF
def fgl(_x, mub): 
    return pdf_set.xfxQ(21, _x, mub)/_x
########################################################


#####---- Sudakov Class ----############################
A1 = CA/2; A2 = (CA/4)*(CA*(67/18 - np.pi**2/6) - 5*nf/9)
B1 = -(CA/2)*(beta0/6)
####---------------------------#### 
# Choice of the total B_CO contribution
# These options assume ζf * ζShF = M**2 + Q**2
if iBCO == 0:
    B1_CO = 0
    BCO_info = '0'
elif iBCO == 1:
    B1_CO = -CA/2
    BCO_info = '1'
elif iBCO == 2:
    B1_CO = -(CA/2)*(1 + np.log(M2/(M2 + Q2)))
    BCO_info = 'full'
else:
    print('ERROR: other B_CO contributions not implemented')
    sys.exit() 
####---------------------------#### 
def Delta_A2(C1): 
    return (CA/4)*beta0*np.log(C1)
def Delta_B1(C1, C2): 
    return -CA*np.log(C2/C1)
#### ---------------------------------------------- ####
def L(_eta2):
    return np.log(_eta2/LQCD**2)
#### ---------------------------------------------- ####
class Sudakov():
    def __init__(self, bT, _x, _mub2, _muh2, C1, C2, err_np=1):
        ####----------------------------------------#### 
        # Perturbative Sudakov
        one_loop = -(4/beta0)*(A1*np.log(_muh2/_mub2) - (B1 + B1_CO + Delta_B1(C1, C2) + A1*L(_muh2))*np.log(L(_muh2)/L(_mub2)))
        nll = one_loop - (16/beta0**2)*(A2 + Delta_A2(C1))*(np.log(L(_muh2)/L(_mub2)) - (1/L(_mub2))*np.log(_muh2/_mub2))
        ####-------------------#### 
        if acc == '1-loop':
            self.pert = one_loop
        elif acc == 'NLL':
            self.pert = nll
        ####----------------------------------------#### 
        # Non-Perturbative Sudakov
        Q_NP = 1.6
        self.np = (A[err_np]*np.log(np.sqrt(_muh2)/Q_NP) + B(_x) + gpsi)*bT**2
########################################################


#####---- Function definition ----######################
def convolution(bT, qT, x, _muh2, C1, C2, C3, C4, err_np=1):
    ####--------------------------------#### 
    # Soft scale
    bpCS = (bT/((1 + (bT/bT_max)**2)**0.5))
    mub = b0/((bpCS**2 + (b0/np.sqrt(_muh2))**2)**0.5)
    mub2 = mub**2
    ####--------------------------------#### 
    # TMD tails and scale
    pdf_tmd = sp.j0(bT*qT)*fgl(x, C3*mub)
    ####--------------------------------#### 
    # Sudakov 
    S = Sudakov(bT, x, C1**2*mub2, _muh2, C1, C2, err_np)
    Sud_full = np.exp(-S.pert)*np.exp(-S.np)
    ####--------------------------------#### 
    # Kinematical and dinamical factors
    # Hard alpha_s is always evaluated from the PDF dataset
    alpha_s = pdf_set.alphasQ2(_muh2)
    fac = (2*np.pi**2)*e_c**2*alpha_QED*alpha_s*sqT_exp(ldme, C4).sigma
    ####--------------------------------#### 
    # Fourier Transform
    return np.nan_to_num((bT/(2*np.pi))*fac*Sud_full*pdf_tmd)
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
print(f'Q value: {Q} GeV\nxB value: {xB}')
print(f'accuracy: {acc}')
print(f'B_CO: {BCO_info}\n')
print(f'PDF set: {pdf_info}\n')
print(f'qT array: ({np.round(qT_sp,6)}, {np.round(qT_ep,3)})\n')
#### ---------------------------------------------- ####
flnm_str = f"{mod_str}{part_str} - pol {str_pol} - {str_fixed_values} - {acc} - {npoints} points - iBCO {iBCO} - mu2 {mu2_info} - ipdf {ipdf}"
if ipdf == 3 or ipdf == 4:
    flnm_str += " (nlo)"
if gpsi != 0:
    flnm_str += f" - gpsi {gpsi}"
#### ---------------------------------------------- ####
filename = f"Jpsi unpol - ep (SIDIS) - plot vs qT (TMD evolved) - {flnm_str}.tsv"
########################################################


####################---- Script ----####################
if __name__ == "__main__":
    warnings.filterwarnings('ignore')   
    ####--------------------------------------------#### 
    if pol == 'T':
        pol_folder = 'transverse'
    elif pol == 'L':
        pol_folder = 'longitudinal'
    os.chdir(f'{basepath}/output')
    ####--------------------------------------------#### 
    f_out = open(filename, 'w')
    f_out.write('# ep -> J/psi X\n')
    f_out.write(f'# {out_pol} (TMD evolved)\n')
    f_out.write(f'\n# ldme -> {ldme.info}\n')
    f_out.write(f'# pdf -> {pdf_info}\n')
    f_out.write('\n# Plot vs qT\n')
    f_out.write(f'# fixed variables:\t{str_fixed_values}\n')
    f_out.write(f'# Regularisation: mu2 = {mu2_info}\n')
    ####--------------------------------------------#### 
    var_hard = [1, 1/2, 2]
    var_soft = [1, 1/2, 2]
    var_pdf = [1, 1/2, 2]
    var_ldme = [1, 1/2, 2]
    ####--------------------------------------------####
    list_var_com = np.array(list(itertools.product(var_hard, var_soft, var_pdf)))
    if npoints == 7:
        f_out.write('\n# qT\tF_mid\t'+
                    'F[C1=2,C2=1/2,C3=1/2]\tF[C1=2,C2=2,C3=2]\t'+
                    'F[C1=0.5,C2=1,C3=1]\tF[C1=0.5,C2=0.5,C3=0.5]\t'+
                    'F[C1=2,C2=1,C3=1]\tF[C1=2,C2=2,C3=2]\t'+
                    'F[mu_psi/2]\tF[2*mu_psi]\tF[A_np=0.05]\tF[A_np=0.8]\n')
        # First element of the table is the mid value
        index_rm = [0]
        # Remove extremes from list
        for i in range(len(list_var_com)):
            check_05 = 0
            check_2 = 0
            for j in range(len(list_var_com[0])):
                if list_var_com[i][j] == 1/2:
                    check_05 += 1
                elif list_var_com[i][j] == 2:
                    check_2 += 1
                else:
                    None
            if (check_05 and check_2) or list_var_com[i][2] != list_var_com[i][1]: 
                index_rm.append(i)
            else:
                None
        list_var_com = np.delete(list_var_com, index_rm, 0)
    elif npoints == 9:
        f_out.write('\n# qT\tF_mid\t'+
                    'F[C1=1,C2=1/2,C3=1/2]\tF[C1=1,C2=2,C3=2]\t'+
                    'F[C1=0.5,C2=1,C3=1]\tF[C1=0.5,C2=0.5,C3=0.5]\tF[C1=0.5,C2=2,C3=2]\t'+
                    'F[C1=2,C2=1,C3=1]\tF[C1=2,C2=0.5,C3=0.5]\tF[C1=2,C2=2,C3=2]\t'+
                    'F[mu_psi/2]\tF[2*mu_psi]\tF[A_np=0.05]\tF[A_np=0.8]\n')
        # First element of the table is the mid value
        index_rm = [0]
        for i in range(len(list_var_com)):
            if list_var_com[i][2] != list_var_com[i][1]: 
                index_rm.append(i)
            else:
                None
        list_var_com = np.delete(list_var_com, index_rm, 0)
    elif npoints == 14:
        f_out.write('\n# qT\tF_mid\t'+
                    'F[C1=1,C2=1,C3=0.5]\tF[C1=1,C2=1,C3=2]\tF[C1=1,C2=0.5,C3=1]\tF[C1=1,C2=0.5,C3=0.5]\tF[C1=1,C2=2,C3=1]\tF[C1=1,C2=2,C3=2]\t'+
                    'F[C1=0.5,C2=1,C3=1]\tF[C1=0.5,C2=1,C3=0.5]\tF[C1=0.5,C2=0.5,C3=1]\tF[C1=0.5,C2=0.5,C3=0.5]\t'+
                    'F[C1=2,C2=1,C3=1]\tF[C1=2,C2=1,C3=2]\tF[C1=2,C2=2,C3=1]\tF[C1=2,C2=2,C3=2]\t'+
                    'F[mu_psi/2]\tF[2*mu_psi]\tF[A_np=0.05]\tF[A_np=0.8]\n')
        # First element of the table is the mid value
        index_rm = [0]
        # Remove extremes from list
        for i in range(len(list_var_com)):
            check_05 = 0
            check_2 = 0
            for j in range(len(list_var_com[0])):
                if list_var_com[i][j] == 1/2:
                    check_05 += 1
                elif list_var_com[i][j] == 2:
                    check_2 += 1
                else:
                    None
            if check_05 and check_2: 
                index_rm.append(i)
            else:
                None
        list_var_com = np.delete(list_var_com, index_rm, 0)
    elif npoints == 21:
        f_out.write('\n# qT\tF_mid\t'+
                    'F[C1=1,C2=1,C3=0.5]\tF[C1=1,C2=1,C3=2]\tF[C1=1,C2=0.5,C3=1]\tF[C1=1,C2=0.5,C3=0.5]\tF[C1=1,C2=2,C3=1]\tF[C1=1,C2=2,C3=2]\t'+
                    'F[C1=0.5,C2=1,C3=1]\tF[C1=0.5,C2=1,C3=0.5]\tF[C1=0.5,C2=1,C3=2]\tF[C1=0.5,C2=0.5,C3=1]\tF[C1=0.5,C2=0.5,C3=0.5]\tF[C1=0.5,C2=2,C3=1]\tF[C1=0.5,C2=2,C3=2]\t'+
                    'F[C1=2,C2=1,C3=1]\tF[C1=2,C2=1,C3=0.5]\tF[C1=2,C2=1,C3=2]\tF[C1=2,C2=0.5,C3=1]\tF[C1=2,C2=0.5,C3=0.5]\tF[C1=2,C2=2,C3=1]\tF[C1=2,C2=2,C3=2]\t'+
                    'F[mu_psi/2]\tF[2*mu_psi]\tF[A_np=0.05]\tF[A_np=0.8]\n')
        # First element of the table is the mid value
        index_rm = [0]
        # Remove extremes from list
        for i in range(len(list_var_com)):
            check_05 = 0
            check_2 = 0
            for j in range(len(list_var_com[0])-1):
                if list_var_com[i][j+1] == 1/2:
                    check_05 += 1
                elif list_var_com[i][j+1] == 2:
                    check_2 += 1
                else:
                    None
            if check_05 and check_2: 
                index_rm.append(i)
            else:
                None
        list_var_com = np.delete(list_var_com, index_rm, 0)
    elif npoints == 27:
        f_out.write('\n# qT\tF_mid\t'+
                    'F[C1=1,C2=1,C3=0.5]\tF[C1=1,C2=1,C3=2]\tF[C1=1,C2=0.5,C3=1]\tF[C1=1,C2=0.5,C3=0.5]\tF[C1=1,C2=0.5,C3=2]\tF[C1=1,C2=2,C3=1]\tF[C1=1,C2=2,C3=0.5]\tF[C1=1,C2=2,C3=2]\t'+
                    'F[C1=0.5,C2=1,C3=1]\tF[C1=0.5,C2=1,C3=0.5]\tF[C1=0.5,C2=1,C3=2]\tF[C1=0.5,C2=0.5,C3=1]\tF[C1=0.5,C2=0.5,C3=0.5]\tF[C1=0.5,C2=0.5,C3=2]\tF[C1=0.5,C2=2,C3=1]\tF[C1=0.5,C2=2,C3=0.5]\tF[C1=0.5,C2=2,C3=2]\t'+
                    'F[C1=2,C2=1,C3=1]\tF[C1=2,C2=1,C3=0.5]\tF[C1=2,C2=1,C3=2]\tF[C1=2,C2=0.5,C3=1]\tF[C1=2,C2=0.5,C3=0.5]\tF[C1=2,C2=0.5,C3=2]\tF[C1=2,C2=2,C3=1]\tF[C1=2,C2=2,C3=0.5]\tF[C1=2,C2=2,C3=2]\t'+
                    'F[mu_psi/2]\tF[2*mu_psi]\tF[A_np=0.05]\tF[A_np=0.8]\n')
        # First element of the table is the mid value
        index_rm = [0]
        list_var_com = np.delete(list_var_com, index_rm, 0)
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
        if not x(qT) <= 1:
            print(f'Change xB value (current {xB})\nLight-cone fraction exceed physical bound: x={x(qT)}')
            sys.exit()
        ####----------------------------------------####
        # To evaluate asymmetry, another intergration method is advised
        tmd_ev = quad(convolution, 0, np.infty, epsabs=1.e-20, limit=50000, args=(qT, x(qT), mu2, 1, 1, 1, 1))
        tmd_ev_var_list = ''
        if scale_err:
            for i in range(len(list_var_com)):
                tmd_ev_var = quad(convolution, 0, np.infty, epsabs=1.e-20, limit=50000, args=(qT, x(qT), mu2*list_var_com[i][0]**2, list_var_com[i][1], list_var_com[i][0], list_var_com[i][2], 1))[0]
                tmd_ev_var_list += f'\t{tmd_ev_var}'
        if ldme_err:
            tmd_var_ldme1 = quad(convolution, 0, np.infty, epsabs=1.e-20, limit=50000, args=(qT, x(qT), mu2, 1, 1, 1, 1/2))[0]
            tmd_var_ldme2 = quad(convolution, 0, np.infty, epsabs=1.e-20, limit=50000, args=(qT, x(qT), mu2, 1, 1, 1, 2))[0]
            tmd_ev_var_list += f'\t{tmd_var_ldme1}\t{tmd_var_ldme2}' 
        else:
            tmd_ev_var_list += '\t0\t0' 
        if np_err:
            tmd_var_np1 = quad(convolution, 0, np.infty, epsabs=1.e-20, limit=50000, args=(qT, x(qT), mu2, 1, 1, 1, 1, 0))[0]
            tmd_var_np2 = quad(convolution, 0, np.infty, epsabs=1.e-20, limit=50000, args=(qT, x(qT), mu2, 1, 1, 1, 1, 2))[0]
            tmd_ev_var_list += f'\t{tmd_var_np1}\t{tmd_var_np2}' 
        else:
            tmd_ev_var_list += '\t0\t0' 
        ####----------------------------------------####
        if qT == qT_array[0]:
            print(f'qT value: {qT}')
            print(f'TMD evolved: {tmd_ev}')
        elif index == int(inum/6):
            print(f'qT value: {qT}')
            print(f'TMD evolved: {tmd_ev}')
            index = 0
        ####----------------------------------------####
        f_out.write(f'{qT}\t{tmd_ev[0]}{tmd_ev_var_list}\n')
    ####--------------------------------------------#### 
    f_out.close()
########################################################