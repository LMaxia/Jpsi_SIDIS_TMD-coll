import sys
#####---- Charm mass ----###############################
mc = 1.5
mc2 = mc**2
########################################################


#####---- LDME class ----###############################
class LDME:
   def __init__(self, j3s1CS, j1s0CO, j3s1CO, j3p0CO, j1s0CO_err=0, j3s1CO_err=0, j3p0CO_err=0):
      # scale is given in terms of M (quarkonium mass)
      self.j3s1CS = j3s1CS
      self.j1s0CO = j1s0CO
      self.j3s1CO = j3s1CO
      self.j3p0CO = j3p0CO
      self.j3p1CO = 3.0*j3p0CO
      self.j3p2CO = 5.0*j3p0CO
   ##### -------------------------- #####
      self.j1s0CO_err = j1s0CO_err
      self.j3s1CO_err = j3s1CO_err
      self.j3p0CO_err = j3p0CO_err
#### ---------------------------------------------- ####
   def j1s0CO_var(self, scale):
      if scale==1:
         return self.j1s0CO
      ##### ----------------------- #####
      elif scale==1/2:
         return self.j1s0CO - self.j1s0CO_err
      ##### ----------------------- #####
      elif scale==2:
         return self.j1s0CO + self.j1s0CO_err
      ##### ----------------------- #####
      else:
         print(f'scale variation is uncorrect: {scale}')
         sys.exit()
#### ---------------------------------------------- ####
   def j3s1CO_var(self, scale):
      if scale==1:
         return self.j3s1CO
      ##### ----------------------- #####
      elif scale==1/2:
         return self.j3s1CO - self.j3s1CO_err
      ##### ----------------------- #####
      elif scale==2:
         return self.j3s1CO + self.j3s1CO_err
      ##### ----------------------- #####
      else:
         print(f'scale variation is uncorrect: {scale}')
         sys.exit()
#### ---------------------------------------------- ####
   def j3p0CO_var(self, scale):
      if scale==1:
         return self.j3p0CO
      ##### ----------------------- #####
      elif scale==1/2:
         return self.j3p0CO - self.j3p0CO_err
      ##### ----------------------- #####
      elif scale==2:
         return self.j3p0CO + self.j3p0CO_err
      ##### ----------------------- #####
      else:
         print(f'scale variation is uncorrect: {scale}')
         sys.exit()
########################################################


#####---- LDME selection ----###########################
def sel(ildme):
   if ildme == 0:
      return ldmeCS
   elif ildme == 1:
      return ldmeCM
   elif ildme == 2:
      return ldmeBKfd
   elif ildme == 3:
      return ldmeBK
   elif ildme == 4:
      return ldmeSY
   elif ildme == 5:
      return ldmeGW
   elif ildme == 6:
      return ldmeSV
   elif ildme == 7:
      return ldmeKZ
   elif ildme == 8:
      return ldmeBC
########################################################


#####---- LDME sets info ----###########################
#### ------------------- set 0 -------------------- ####
ldmeCS = LDME(1.4, 0, 0, 0)
ldmeCS.info = 'CS channel'
#### ---------------------------------------------- ####

#### ------------------- set 1 -------------------- ####
ldmeCM = LDME(1.16, 0.089, 0.003, 0.0056*mc2, 0.0098, 0.00012, 0.0021*mc2)
ldmeCM.s = 'CM12'
ldmeCM.info = 'CM12 - 10.1103/PhysRevLett.108.242004'
# Authors: K.-T. Chao, Y.-Q. Ma, H.-S. Shao, K. Wang, Y.-J. Zhang
# Doi: https://link.aps.org/doi/10.1103/PhysRevLett.108.242004
# Details: Fit with polarisation production included
#### ---------------------------------------------- ####

#### ------------------- set 2 -------------------- ####
ldmeBKfd = LDME(1.32, 0.0304, 0.00168, -0.00908, 0.0035, 0.00046, 0.0161)
ldmeBKfd.s = 'BK11fd'
ldmeBKfd.info = 'BK11fd - 0.1103/PhysRevD.84.051501 (with feed down)'
# Authors: M. Butenschoen, B.A. Kniehl
# Doi: https://link.aps.org/doi/10.1103/PhysRevD.84.051501
# Details: NLO WITH feed-down contributions
#### ---------------------------------------------- ####

#### ------------------- set 3 -------------------- ####
ldmeBK = LDME(1.32, 0.0497, 0.00224, -0.0161, 0.0044, 0.00059, 0.0020)
ldmeBK.s = 'BK11'
ldmeBK.info = 'BK11 - 0.1103/PhysRevD.84.051501 (prompt)'
# Authors: M. Butenschoen, B.A. Kniehl
# Doi: https://link.aps.org/doi/10.1103/PhysRevD.84.051501
# Details: NLO prompt Jpsi production (Table I)
#### ---------------------------------------------- ####

#### ------------------- set 4 -------------------- ####
ldmeSY = LDME(0., 0.1423, -0.0093, -0.0175*mc2)
ldmeSY.s = 'SYY13'
ldmeSY.info = 'SYY13 - 10.1103/PhysRevD.88.054008'
# Authors: P. Sun, C.-P. Yuan, F. Yuan
# Doi: https://link.aps.org/doi/10.1103/PhysRevD.88.054008
# Details: Obtained with CSS resummation; j3s1CS is set as 0
#### ---------------------------------------------- ####

#### ------------------- set 5 -------------------- ####
ldmeGW = LDME(1.16, 0.097, -0.0046, -0.0095*mc2)
ldmeGW.s = 'G13'
ldmeGW.info = 'G13 - 10.1103/PhysRevLett.110.042002'
# Authors: B. Gong, L.-P. Wan, J.-X. Wang, H.-F. Zhang
# Doi: https://link.aps.org/doi/10.1103/PhysRevLett.110.042002
# Details: --
#### ---------------------------------------------- ####

#### ------------------- set 6 -------------------- ####
ldmeSV = LDME(1.2, 0.018, 0.0013, 0.018*mc2, 0.0087, 0.0013, 0.0087*mc2)
ldmeSV.s = 'SV13'
ldmeSV.info = 'SV13 - 10.1103/PhysRevC.87.044905'
# Authors: R. Sharma, I. Vitev
# Doi: https://link.aps.org/doi/10.1103/PhysRevC.87.044905
# Details: LO extraction
#### ---------------------------------------------- ####

#### ------------------- set 7 -------------------- ####
#-----From F. Kniehl, L. Zwirner, Nuclear Physics B 621, 337-358 (2002)
#      HERA data for CTEQ5L, for kappa=0.5
ldmeKZ = LDME(1.4, 0.033, 0.0039, 0.033*mc**2/3.4)
ldmeKZ.s = 'KZ02'
ldmeKZ.info = 'KZ02 - 10.1016/S0550-3213(01)00564-8 (k=0.5)'
# Authors: B.A. Kniehl, L. Zwirner
# Doi: https://doi.org/10.1016/S0550-3213(01)00564-8
# Details: valid for HERA and CTEQ5L; values are obtained with kappa=0.5
#### ---------------------------------------------- ####

#### ------------------- set 8 -------------------- ####
ldmeBC = LDME(1.18, -0.0476, 0.017, 0.03*mc2)
ldmeBC.s = 'BCVW22'
ldmeBC.info = 'BCVW22 - 10.1103/PhysRevD.105.L111503'
# Authors: N. Brambilla, H. S. Chung, A. Vairo, X.-P. Wang
# Doi: https://doi.org/10.1103/PhysRevD.105.L111503
# Details: Recent global fit
#### ---------------------------------------------- ####
########################################################