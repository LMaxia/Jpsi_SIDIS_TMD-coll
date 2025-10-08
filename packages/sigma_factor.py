#####---- Importing variables ----######################
import os
import yaml
basepath = os.path.dirname(os.path.dirname(__file__))
with open(os.path.join(basepath, 'input.yaml'), 'r') as infile:
    input_dict = yaml.safe_load(infile)
M = input_dict['Mjp']
Q = input_dict['Q']
pol = input_dict['pol']
wave = input_dict['wave']
########################################################


#####---- Hard structure functions selection ----#######
def sigma_UU2phi(ldme, ldme_scale):
   return (1/(M*(M**2 + Q**2)))*(-ldme.j1s0CO_var(ldme_scale) + (4*ldme.j3p0CO_var(ldme_scale))*(3*M**2 - Q**2)/(M**2*(M**2 + Q**2)))
#### ---------------------------------------------- ####
if wave == 10:
   class small_qT_expansion:
      def __init__(self, ldme, ldme_scale):
         sigma_UUT = (1/(M*(M**2 + Q**2)))*(ldme.j1s0CO_var(ldme_scale) + (4*ldme.j3p0CO_var(ldme_scale))*(7*M**4 + 2*M**2*Q**2 + 3*Q**4)/(M**2*(M**2 + Q**2)**2))
         sigma_UUL = (1/(M*(M**2 + Q**2)))*(16*Q**2/((M**2 + Q**2)**2)*ldme.j3p0CO_var(ldme_scale))
         # sigma_UU2phi = (4/(M*(M**2 + Q**2)))*(-ldme.j1s0CO_var(ldme_scale) + (4*ldme.j3p0CO_var(ldme_scale))*(3*M**2 - Q**2)/(M**2*(M**2 + Q**2)))
         sigma_UU2phi = (1/(M*(M**2 + Q**2)))*(-ldme.j1s0CO_var(ldme_scale) + (4*ldme.j3p0CO_var(ldme_scale))*(3*M**2 - Q**2)/(M**2*(M**2 + Q**2)))

         if pol == 'T':
            self.sigma = sigma_UUT
         elif pol == 'L':
            self.sigma = sigma_UUL
         elif pol == 'A':
            self.sigma = sigma_UU2phi
#### ---------------------------------------------- ####

#### ---------------------------------------------- ####
elif wave == 0:
   class small_qT_expansion:
      def __init__(self, ldme, ldme_scale):
         sigma_UUT = (1/(M*(M**2+Q**2)))
         sigma_UUL = 0.
         sigma_UU2phi = (2/(M*(M**2+Q**2)))*(-1)

         if pol == 'T':
            self.sigma = sigma_UUT
         elif pol == 'L':
            self.sigma = sigma_UUL
         elif pol == 'A':
            self.sigma = sigma_UU2phi
#### ---------------------------------------------- ####

#### ---------------------------------------------- ####
elif wave == 2:
   class small_qT_expansion:
      def __init__(self, ldme, ldme_scale):
         sigma_UUT = (1/(M*(M**2+Q**2)))*((4)*(7*M**4 + 2*M**2*Q**2 + 3*Q**4)/(M**2*(M**2 + Q**2)**2))
         sigma_UUL = (1/(M*(M**2+Q**2)))*(16*Q**2/((M**2 + Q**2)**2))
         sigma_UU2phi = (2/(M*(M**2+Q**2)))*((4)*(3*M**2 - Q**2)/(M**2 * (M**2 + Q**2)))

         if pol == 'T':
            self.sigma = sigma_UUT
         elif pol == 'L':
            self.sigma = sigma_UUL
         elif pol == 'A':
            self.sigma = sigma_UU2phi
#### ---------------------------------------------- ####

#### ---------------------------------------------- ####
else: 
   class small_qT_expansion:
      def __init__(self, ldme, ldme_scale):
         self.sigma = 0.
#### ---------------------------------------------- ####s