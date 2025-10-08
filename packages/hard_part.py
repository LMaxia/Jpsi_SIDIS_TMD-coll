#####---- Importing variables ----######################
import os
import yaml
basepath = os.path.dirname(os.path.dirname(__file__))
with open(os.path.join(basepath, 'input.yaml'), 'r') as infile:
    input_dict = yaml.safe_load(infile)
e_c = input_dict['e_c']
alpha_QED = input_dict['alpha_QED']
M = input_dict['Mjp']
Q = input_dict['Q']
pol = input_dict['pol']
ig = input_dict['ig']
iq = input_dict['iq']
wave = input_dict['wave']
########################################################


#####---- Hard structure functions selection ----#######
if pol == 'T':
   def FUU(xh, zh, pdf, ldme, ldme_scale):
      return FUUT(xh, zh, pdf, ldme, ldme_scale)
elif pol == 'L':
   def FUU(xh, zh, pdf, ldme, ldme_scale):
      return FUUL(xh, zh, pdf, ldme, ldme_scale)
elif pol == 'A':
   def FUU(xh, zh, pdf, ldme, ldme_scale):
      return FUUcos2phi(xh, zh, pdf, ldme, ldme_scale)
########################################################
   

#####---- Hard structure functions variables ----#######
#### ---------------------------------------------- ####
#### ---- full ----####
if wave == 10: 
   ####---- (gluon) contribution ----####
   if iq == 0:    
      ####---- Transverse ----#### 
      def FUUT(xh, zh, pdf, ldme, ldme_scale):
         FUUT1s0g = (6*xh*(M**10*xh**5 + M**8*Q**2*xh**4*(-3*zh + xh*(2 + 3*zh)) + M**6*Q**4*xh**3*(5*zh**2 - 4*xh*zh*(1 + 2*zh) + xh**2*(2 + 4*zh + 4*zh**2)) + \
            M**4*Q**6*xh**2*(-5*zh**3 - 2*xh**2*zh*(2 + 2*zh + 5*zh**2) + xh*zh*(2 + 2*zh + 11*zh**2) + 2*xh**3*(1 + 3*zh**2 + zh**3)) + \
            Q**10*(-1 + xh)*zh*(-2*xh*zh**2*(1 - zh + zh**2) + (1 - zh + zh**2)**2 + xh**4*(1 - 2*zh + 2*zh**2) - 2*xh**3*zh*(1 - 2*zh + 2*zh**2) + xh**2*zh*(2 - zh + 2*zh**3)) + \
            M**2*Q**8*xh*(1 - 2*zh + 3*zh**2 - 2*zh**3 + 3*zh**4 + 2*xh*zh**2*(-2 + zh - 4*zh**2) - 4*xh**3*zh*(1 - zh + 2*zh**2 + zh**3) + xh**4*(1 + 4*zh**3) + xh**2*zh*(2 + 3*zh + 10*zh**3))))/ \
            (M*Q**6*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*(-1 + zh)*zh**2)
         
         FUUT3s1g = (32*xh**2*(Q**6*(-1 + xh)*xh*(-1 + zh)*zh + M**6*xh**2*(1 - zh + zh**2) + \
            M**2*Q**4*(1 + (-2 + 2*xh - 5*xh**2)*zh + (3 - 6*xh + 10*xh**2)*zh**2 + (-2 + 4*xh - 6*xh**2)*zh**3 + (1 - 2*xh + 2*xh**2)*zh**4) + \
            M**4*Q**2*xh*(zh*(-3 + 3*zh - 2*zh**2) + xh*(-3 + 7*zh - 4*zh**2 + 2*zh**3))))/(27.*M*Q**2*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*zh**2)
         
         FUUT3pjg = (24*xh*(7*M**16*xh**8*zh + M**14*Q**2*xh**7*(8 + (-5 + 23*xh)*zh + (-31 + 21*xh)*zh**2) + \
            M**12*Q**4*xh**6*(8 + 3*(-11 + 18*xh + 7*xh**2)*zh + (19 - 126*xh + 71*xh**2)*zh**2 + 4*(15 - 19*xh + 7*xh**2)*zh**3 + 2*zh**4) + \
            M**10*Q**6*xh**5*(20 - 63*zh + 87*zh**2 - 50*zh**3 - 56*zh**4 - 8*zh**5 + xh**3*zh*(-3 + 67*zh + 110*zh**2 + 14*zh**3) - xh**2*(48 - 147*zh + 77*zh**2 + 262*zh**3 + 92*zh**4) + \
            xh*(-8 + 87*zh - 259*zh**2 + 304*zh**3 + 112*zh**4 + 4*zh**5)) + \
            Q**16*(-1 + xh)**2*zh**2*(xh**6*(3 - 6*zh + 6*zh**2) + xh**5*(4 - 20*zh + 34*zh**2 - 24*zh**3) + 3*(-1 + 2*zh - 2*zh**2 + zh**3)**2 + \
            xh**4*(3 - 14*zh + 44*zh**2 - 60*zh**3 + 36*zh**4) - 2*xh**3*(2 - 7*zh + 3*zh**2 + 11*zh**3 - 18*zh**4 + 12*zh**5) - \
            2*xh*(-2 + 7*zh - 10*zh**2 + 4*zh**3 + 4*zh**4 - 6*zh**5 + 3*zh**6) + xh**2*(3 + 4*zh - 28*zh**2 + 36*zh**3 - 20*zh**4 + 2*zh**5 + 6*zh**6)) + \
            M**8*Q**8*xh**4*(20 - 135*zh + 255*zh**2 - 224*zh**3 + 114*zh**4 + 14*zh**5 + 12*zh**6 + xh**4*zh*(-7 - 19*zh + 154*zh**2 + 62*zh**3) - \
            2*xh**3*(32 - 34*zh - 89*zh**2 + 120*zh**3 + 161*zh**4 + 20*zh**5) + 2*xh**2*(-20 + 134*zh - 203*zh**2 + 15*zh**3 + 226*zh**4 + 58*zh**5 + 2*zh**6) - \
            2*xh*(-10 - 28*zh + 164*zh**2 - 262*zh**3 + 207*zh**4 + 38*zh**5 + 6*zh**6)) + \
            M**6*Q**10*xh**3*(4 - 79*zh + 305*zh**2 - 446*zh**3 + 318*zh**4 - 134*zh**5 + 12*zh**6 - 8*zh**7 + xh**5*zh*(9 - 61*zh + 76*zh**2 + 108*zh**3) - \
            xh**4*(24 + 63*zh - 199*zh**2 - 84*zh**3 + 392*zh**4 + 144*zh**5) + xh*zh*(-35 - 169*zh + 560*zh**2 - 614*zh**3 + 356*zh**4 + 18*zh**5 + 12*zh**6) + \
            2*xh**3*(-12 + 68*zh + 42*zh**2 - 297*zh**3 + 195*zh**4 + 188*zh**5 + 24*zh**6) + xh**2*(12 + 88*zh - 322*zh**2 + 224*zh**3 + 252*zh**4 - 482*zh**5 - 64*zh**6 - 8*zh**7)) + \
            M**2*Q**14*xh*zh*(-5 + 23*zh - 60*zh**2 + 104*zh**3 - 98*zh**4 + 48*zh**5 - 11*zh**6 - zh**7 + xh**7*(3 + 5*zh - 26*zh**2 + 38*zh**3) + \
            xh**6*(1 - 43*zh + 66*zh**2 + 20*zh**3 - 112*zh**4) + xh**5*(-7 + 39*zh + 58*zh**2 - 222*zh**3 + 132*zh**4 + 112*zh**5) + \
            xh**4*(3 + 7*zh - 74*zh**2 + 16*zh**3 + 210*zh**4 - 232*zh**5 - 40*zh**6) + xh**3*(-1 + 17*zh - 146*zh**2 + 296*zh**3 - 268*zh**4 + 54*zh**5 + 110*zh**6 + 2*zh**7) + \
            xh*(-3 + 3*zh + 20*zh**2 - 142*zh**3 + 222*zh**4 - 154*zh**5 + 53*zh**6 + 3*zh**7) - xh**2*(-9 + 51*zh - 162*zh**2 + 110*zh**3 + 86*zh**4 - 172*zh**5 + 112*zh**6 + 4*zh**7)) + \
            M**4*Q**12*xh**2*(4 - 23*zh + 107*zh**2 - 272*zh**3 + 322*zh**4 - 192*zh**5 + 64*zh**6 - 5*zh**7 + 2*zh**8 + xh**6*zh*(11 - 23*zh - 16*zh**2 + 92*zh**3) - \
            2*xh**5*zh*(21 + zh - 106*zh**2 + 82*zh**3 + 96*zh**4) + xh**4*zh*(-19 + 253*zh - 342*zh**2 - 156*zh**3 + 408*zh**4 + 120*zh**5) - \
            2*xh**3*(-6 + 14*zh - 53*zh**2 + 223*zh**3 - 343*zh**4 + 139*zh**5 + 124*zh**6 + 12*zh**7) - \
            2*xh*(-2 + 7*zh - 10*zh**2 - 113*zh**3 + 252*zh**4 - 211*zh**5 + 92*zh**6 + zh**7 + 2*zh**8) + \
            xh**2*(12 - 61*zh - 29*zh**2 + 76*zh**3 + 156*zh**4 - 378*zh**5 + 316*zh**6 + 16*zh**7 + 4*zh**8))))/ \
            (M**3*Q**6*(Q**2*(-1 + xh) + M**2*xh)**3*(M**2*xh + Q**2*(1 + xh - zh))**4*(-1 + zh)*zh**3)
         
         FUUTg = ((FUUT1s0g*ldme.j1s0CO_var(ldme_scale)) + (FUUT3s1g*((15/8)*ldme.j3s1CO_var(ldme_scale) + ldme.j3s1CS)) + (FUUT3pjg*ldme.j3p0CO_var(ldme_scale)))*pdf.gluon_PDF

         return e_c**2*alpha_QED*pdf.als**2*FUUTg

      ####---- Longitudinal ----#### 
      def FUUL(xh, zh, pdf, ldme, ldme_scale):
         FUUL1s0g = (6*xh**2*(2*M**8*xh**4 + Q**8*(-1 + xh)**2*zh**2*(1 + 2*xh**2 + xh*(2 - 4*zh) - 2*zh + 2*zh**2) + 2*M**6*Q**2*xh**3*(-3*zh + 2*xh*(1 + zh)) + \
            2*M**2*Q**6*(-1 + xh)*xh*zh*(2*xh**2*(1 + zh) + zh*(-1 + 3*zh) + xh*(1 - 4*zh - 2*zh**2)) + M**4*Q**4*xh**2*(-2*xh*zh*(4 + 5*zh) + zh*(-2 + 9*zh) + 2*xh**2*(1 + 4*zh + zh**2))))/ \
            (M*Q**4*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*zh**2)
         
         FUUL3s1g = (16*xh**2*(Q**4*(-1 + xh)**2*zh**2 + 2*M**2*Q**2*(-1 + xh)*xh*zh*(-1 + 6*zh - 6*zh**2 + 2*zh**3) + M**4*xh**2*(-2 + 10*zh - 11*zh**2 + 4*zh**3)))/ \
            (27.*M*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*zh**2)
         
         FUUL3pjg = (48*xh**2*(M**14*xh**7*zh*(-5 + 7*zh) + M**12*Q**2*xh**6*(-4 - 3*(-5 + 7*xh)*zh + (7 + 19*xh)*zh**2 + 2*(-13 + 7*xh)*zh**3) + \
            Q**14*(-1 + xh)**3*(-1 + zh)*zh**3*(2 + 3*xh**4 + xh**3*(8 - 12*zh) - 6*zh + 11*zh**2 - 8*zh**3 + 3*zh**4 + xh**2*(11 - 24*zh + 18*zh**2) + xh*(6 - 22*zh + 24*zh**2 - 12*zh**3)) + \
            M**10*Q**4*xh**5*(xh**2*zh*(-34 + 2*zh + 55*zh**2 + 7*zh**3) - 2*xh*(6 - 18*zh - 31*zh**2 + 40*zh**3 + 23*zh**4) + 2*(-2 + 15*zh - 34*zh**2 + 9*zh**3 + 19*zh**4 + zh**5)) + \
            M**2*Q**12*(-1 + xh)**2*xh*zh**2*(xh**4*(-4 - 13*zh + 19*zh**2) - 2*xh**3*(7 - 9*zh - 24*zh**2 + 28*zh**3) + xh**2*(-16 + 49*zh - 37*zh**2 - 46*zh**3 + 56*zh**4) + \
            zh*(14 - 50*zh + 69*zh**2 - 47*zh**3 + 15*zh**4 + zh**5) - 2*xh*(5 - 31*zh + 52*zh**2 - 36*zh**3 + 2*zh**4 + 10*zh**5)) + \
            M**8*Q**6*xh**4*(xh**3*zh*(-26 - 42*zh + 77*zh**2 + 31*zh**3) - xh**2*(12 - 18*zh - 138*zh**2 + 43*zh**3 + 161*zh**4 + 20*zh**5) + \
            2*xh*(-4 + 29*zh - 48*zh**2 - 55*zh**3 + 80*zh**4 + 29*zh**5 + zh**6) - 2*(2 - 17*zh + 53*zh**2 - 79*zh**3 + 38*zh**4 + 10*zh**5 + 3*zh**6)) + \
            M**4*Q**10*xh**2*zh*(7 - 61*zh + 172*zh**2 - 236*zh**3 + 180*zh**4 - 72*zh**5 + 4*zh**6 - 2*zh**7 + xh**5*(-1 - 25*zh - 8*zh**2 + 46*zh**3) + \
            xh**4*(-9 + 23*zh + 124*zh**2 - 82*zh**3 - 96*zh**4) + 2*xh**3*(-1 + 28*zh - 79*zh**2 - 48*zh**3 + 102*zh**4 + 30*zh**5) - \
            2*xh**2*(1 - 21*zh + 72*zh**2 - 144*zh**3 + 54*zh**4 + 62*zh**5 + 6*zh**6) + xh*(7 - 35*zh + 14*zh**2 + 80*zh**3 - 180*zh**4 + 136*zh**5 + 8*zh**6 + 2*zh**7)) + \
            M**6*Q**8*xh**3*(-4 + 39*zh - 139*zh**2 + 246*zh**3 - 244*zh**4 + 116*zh**5 - 4*zh**6 + 6*zh**7 + xh**4*zh*(-9 - 53*zh + 38*zh**2 + 54*zh**3) - \
            4*xh**3*(1 + 3*zh - 28*zh**2 - 23*zh**3 + 49*zh**4 + 18*zh**5) + 2*xh**2*(-2 + 13*zh + 10*zh**2 - 141*zh**3 + 62*zh**4 + 94*zh**5 + 12*zh**6) - \
            2*xh*(2 - 16*zh + 36*zh**2 - 19*zh**3 - 87*zh**4 + 96*zh**5 + 16*zh**6 + 2*zh**7))))/(M**3*(Q**2*(-1 + xh) + M**2*xh)**3*(M**2*Q*xh + Q**3*(1 + xh - zh))**4*(-1 + zh)*zh**3)

         FUULg = ((FUUL1s0g*ldme.j1s0CO_var(ldme_scale)) + (FUUL3s1g*((15/8)*ldme.j3s1CO_var(ldme_scale) + ldme.j3s1CS)) + (FUUL3pjg*ldme.j3p0CO_var(ldme_scale)))*pdf.gluon_PDF
      
         return e_c**2*alpha_QED*pdf.als**2*FUULg
      
      ####---- Cos2phi ----####
      def FUUcos2phi(xh, zh, pdf, ldme, ldme_scale):
         FUUcos2phi1s0g = 4*(3*xh*(2*M**8*xh**5*(-1 + zh) + M**6*Q**2*xh**3*((-2 + zh)*zh - 6*xh*(-1 + zh)*zh + 4*xh**2*(-1 + zh**2)) + \
            M**4*Q**4*xh**2*(2*xh**2*zh*(4 + zh - 5*zh**2) + zh*(2 + 2*zh - zh**2) + xh*zh*(-2 - 10*zh + 9*zh**2) + 2*xh**3*(-1 - 3*zh + 3*zh**2 + zh**3)) + \
            Q**8*(-1 + xh)**2*zh**2*(2*xh**3*(-1 + zh) + zh + xh**2*(-2 + 6*zh - 4*zh**2) + xh*(-2 + 3*zh - 4*zh**2 + 2*zh**3)) + \
            M**2*Q**6*(-1 + xh)*xh*zh*(3*zh + 4*xh**3*(-1 + zh**2) - 2*xh**2*(1 - 5*zh + 2*zh**2 + 2*zh**3) + xh*(-2 + zh - 8*zh**2 + 6*zh**3))))/ \
            (M*Q**4*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*(-1 + zh)*zh**2)
         
         FUUcos2phi3s1g = 4*(-16*M*xh**3*(-1 + zh)*(M**2*xh + Q**2*(-1 + xh)*zh)*(M**2 + Q**2*(-1 + 4*zh - 2*zh**2)))/ \
            (27.*Q**2*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*zh**2)
         
         FUUcos2phi3pjg = 4*(12*xh*(2*M**14*xh**7*(-1 + zh)*(-4 + (3 + 7*xh)*zh) + M**12*Q**2*xh**6*(8 - 38*zh + 50*zh**2 - 17*zh**3 + 2*xh**2*zh*(-31 + 17*zh + 14*zh**2) + \
            xh*(16 - 8*zh + 38*zh**2 - 46*zh**3)) + M**10*Q**4*xh**5*(8 - 38*zh + 68*zh**2 - 65*zh**3 + 15*zh**4 + 2*xh**3*zh*(-54 - 8*zh + 55*zh**2 + 7*zh**3) - \
            2*xh**2*zh*(-18 - 91*zh + 63*zh**2 + 46*zh**3) + xh*(8 - 10*zh - 62*zh**2 + 19*zh**3 + 55*zh**4 + 4*zh**5)) + \
            Q**14*(-1 + xh)**3*zh**3*(6*xh**5*(-1 + zh) + (-1 + zh)**2*zh - 8*xh**4*(2 - 5*zh + 3*zh**2) + xh**3*(-22 + 69*zh - 84*zh**2 + 36*zh**3) + \
            xh**2*(-16 + 59*zh - 90*zh**2 + 72*zh**3 - 24*zh**4) + xh*(-4 + 21*zh - 38*zh**2 + 37*zh**3 - 22*zh**4 + 6*zh**5)) + \
            M**2*Q**12*(-1 + xh)**2*xh*zh**2*(2*xh**5*(-6 - 13*zh + 19*zh**2) + xh**4*(-38 + 54*zh + 96*zh**2 - 112*zh**3) + zh*(7 - 15*zh + 11*zh**2 - 3*zh**3) + \
            xh**3*(-44 + 131*zh - 109*zh**2 - 92*zh**3 + 112*zh**4) - xh**2*(38 - 151*zh + 239*zh**2 - 174*zh**3 + 8*zh**4 + 40*zh**5) + \
            xh*(-12 + 67*zh - 133*zh**2 + 151*zh**3 - 103*zh**4 + 30*zh**5 + 2*zh**6)) + \
            M**8*Q**6*xh**4*(8 - 42*zh + 92*zh**2 - 86*zh**3 + 49*zh**4 - 3*zh**5 + 2*xh**4*zh*(-46 - 62*zh + 77*zh**2 + 31*zh**3) - \
            2*xh*zh*(-22 + 105*zh - 167*zh**2 + 92*zh**3 + 8*zh**4 + 6*zh**5) - 2*xh**3*(8 - 8*zh - 182*zh**2 + zh**3 + 161*zh**4 + 20*zh**5) + \
            xh**2*(-8 + 94*zh - 186*zh**2 - 244*zh**3 + 249*zh**4 + 116*zh**5 + 4*zh**6)) + \
            M**4*Q**10*xh**2*zh*(2*xh**6*(-3 - 35*zh - 8*zh**2 + 46*zh**3) + zh*(-2 - 11*zh + 36*zh**2 - 33*zh**3 + 13*zh**4) - \
            2*xh**5*(12 - 31*zh - 159*zh**2 + 82*zh**3 + 96*zh**4) + xh**4*zh*(128 - 363*zh - 288*zh**2 + 408*zh**3 + 120*zh**4) - \
            2*xh**3*zh*(-23 + 145*zh - 327*zh**2 + 81*zh**3 + 124*zh**4 + 12*zh**5) - \
            4*xh*(-4 + 26*zh - 70*zh**2 + 110*zh**3 - 93*zh**4 + 37*zh**5 - 2*zh**6 + zh**7) + xh**2*(14 - 60*zh + 82*zh**2 + 110*zh**3 - 393*zh**4 + 263*zh**5 + 16*zh**6 + 4*zh**7)) + \
            M**6*Q**8*xh**3*(2*xh**5*zh*(-19 - 73*zh + 38*zh**2 + 54*zh**3) - zh*(12 - 58*zh + 108*zh**2 - 84*zh**3 + 33*zh**4 + zh**5) - \
            2*xh**4*(4 + 19*zh - 145*zh**2 - 146*zh**3 + 196*zh**4 + 72*zh**5) + 2*xh**3*(-4 + 33*zh + 17*zh**2 - 316*zh**3 + 68*zh**4 + 188*zh**5 + 24*zh**6) - \
            2*xh**2*(4 - 41*zh + 122*zh**2 - 87*zh**3 - 173*zh**4 + 167*zh**5 + 32*zh**6 + 4*zh**7) + \
            xh*(-8 + 68*zh - 252*zh**2 + 534*zh**3 - 554*zh**4 + 265*zh**5 - 17*zh**6 + 12*zh**7))))/ \
            (M**3*(Q**2*(-1 + xh) + M**2*xh)**3*(M**2*Q*xh + Q**3*(1 + xh - zh))**4*(-1 + zh)*zh**3)

         FUUcos2phig = ((FUUcos2phi1s0g*ldme.j1s0CO_var(ldme_scale)) + (FUUcos2phi3s1g*((15/8)*ldme.j3s1CO_var(ldme_scale) + ldme.j3s1CS)) + (FUUcos2phi3pjg *ldme.j3p0CO_var(ldme_scale)))*pdf.gluon_PDF
      
         return e_c**2*alpha_QED*pdf.als**2*FUUcos2phig

   ####---- (quark) contribution ----####
   elif ig == 0:  
      ####---- Transverse ----#### 
      def FUUT(xh, zh, pdf, ldme, ldme_scale):
         FUUT1s0q = (-4*xh*(M**4*xh**2 + 2*M**2*Q**2*(-1 + xh)*xh*zh + Q**4*(1 + zh**2 - 2*xh*(1 - zh + zh**2) + xh**2*(1 - 2*zh + 2*zh**2))))/(3.*M*Q**4*(M**2*xh + Q**2*(1 + xh - zh))**2*(-1 + zh))

         FUUT3s1q = (-2*(M**6*xh**3 - M**4*Q**2*xh**2*(2 + xh*(-2 + zh) + zh) + M**2*Q**4*xh*(-4*xh*zh + (1 + zh)**2 + xh**2*(-1 + 2*zh)) + Q**6*(-1 + xh)*zh*(1 + zh**2 - 2*xh*zh**2 + xh**2*(1 - 2*zh + 2*zh**2))))/ \
            (9.*M**3*Q**4*(-1 + xh)**2*(M**2*xh - Q**2*zh)**2)  
         
         FUUT3pjq = (-16*xh*(7*M**8*xh**4 + 2*M**6*Q**2*xh**3*(1 + xh - 8*zh + 7*xh*zh) + 2*M**4*Q**4*xh**2*(-3 + 3*zh + 7*zh**2 + xh*(6 - 6*zh - 9*zh**2) + xh**2*(-3 + 3*zh + 7*zh**2)) + \
            2*M**2*Q**6*xh*(1 + 4*zh - zh**2 - 4*zh**3 + xh**3*(1 - 7*zh + 10*zh**2) - xh**2*(9 - 22*zh + 16*zh**2 + 2*zh**3) + xh*(-9 + 5*zh + zh**2 + 5*zh**3)) + \
            Q**8*(3*(-1 + zh)**2*(1 + zh**2) + xh**4*(3 - 6*zh + 6*zh**2) - 2*xh**3*(2 - 5*zh**2 + 6*zh**3) - 2*xh*(2 - 8*zh + 7*zh**2 - 4*zh**3 + 3*zh**4) + 2*xh**2*(1 - 2*zh + zh**3 + 3*zh**4))))/ \
            (3.*M**3*(M**2*Q*xh + Q**3*(1 + xh - zh))**4*(-1 + zh))
         
         FUUTq = ((FUUT1s0q*ldme.j1s0CO_var(ldme_scale)) + (FUUT3pjq*ldme.j3p0CO_var(ldme_scale)))*pdf.qqb_PDF + (FUUT3s1q*ldme.j3s1CO_var(ldme_scale))*pdf.qqb_PDF_weight

         return e_c**2*alpha_QED*pdf.als**2*FUUTq

      ####---- Longitudinal ----#### 
      def FUUL(xh, zh, pdf, ldme, ldme_scale):
         FUUL1s0q = (-8*xh**2*(M**2*xh + Q**2*(-1 + xh)*zh))/(3.*M*(M**2*Q*xh + Q**3*(1 + xh - zh))**2)

         FUUL3s1q = (-4*Q**2*xh*(-1 + zh)*zh**2)/(9.*M**3*(M**2*xh - Q**2*zh)**2)

         FUUL3pjq = (-32*xh**2*(M**6*xh**3*(-5 + 7*zh) + Q**6*(-1 + xh)*(-1 + zh)*zh*(3 + 3*xh**2 + xh*(2 - 6*zh) - 2*zh + 3*zh**2) + M**4*Q**2*xh**2*(2 + 3*zh - 9*zh**2 + xh*(-6 + 3*zh + 7*zh**2)) + \
            M**2*Q**4*xh*(7 - 13*zh + 5*zh**2 + 5*zh**3 + xh**2*(-1 - 7*zh + 10*zh**2) - 2*xh*(-1 - 6*zh + 8*zh**2 + zh**3))))/(3.*M**3*Q**2*(M**2*xh + Q**2*(1 + xh - zh))**4*(-1 + zh))
         
         FUULq = ((FUUL1s0q*ldme.j1s0CO_var(ldme_scale)) + (FUUL3pjq*ldme.j3p0CO_var(ldme_scale)))*pdf.qqb_PDF + (FUUL3s1q*ldme.j3s1CO_var(ldme_scale))*pdf.qqb_PDF_weight

         return e_c**2*alpha_QED*pdf.als**2*FUULq
      
      ####---- Cos2phi ----#### 
      def FUUcos2phi(xh, zh, pdf, ldme, ldme_scale):         
         FUUcos2phi1s0q = 4*(-4*xh*(1 + xh*(-1 + zh))*(M**2*xh + Q**2*(-1 + xh)*zh))/(3.*M*(M**2*Q*xh + Q**3*(1 + xh - zh))**2*(-1 + zh))

         FUUcos2phi3s1q = 4*(2*xh*(-1 + zh)*(M**2 - Q**2*zh)*(M**2*xh + Q**2*(-1 + xh)*zh))/(9.*M**3*(-1 + xh)*(M**2*Q*xh - Q**3*zh)**2)

         FUUcos2phi3pjq = 4*(-16*xh*(M**2*xh + Q**2*(-1 + xh)*zh)*(M**4*xh**2*(-3 + 7*xh*(-1 + zh)) + 2*M**2*Q**2*xh*(-1 + 5*xh**2*(-1 + zh) + zh - xh*zh**2) + \
            Q**4*(3*xh**3*(-1 + zh) + (-1 + zh)**2 + xh**2*(-1 + 8*zh - 6*zh**2) + xh*(3 - zh - 5*zh**2 + 3*zh**3))))/ \
            (3.*M**3*Q**2*(M**2*xh + Q**2*(1 + xh - zh))**4*(-1 + zh))

         FUUcos2phiq = ((FUUcos2phi1s0q * ldme.j1s0CO_var(ldme_scale)) + (FUUcos2phi3pjq * ldme.j3p0CO_var(ldme_scale))) * pdf.qqb_PDF + ( FUUcos2phi3s1q * ldme.j3s1CO_var(ldme_scale)) * pdf.qqb_PDF_weight

         return e_c**2 * alpha_QED * pdf.als**2 * FUUcos2phiq

   ####----(gluon + quark) contribution ----####
   else:      
      ####---- Transverse ----#### 
      def FUUT(xh, zh, pdf, ldme, ldme_scale):
         FUUT1s0g = (6*xh*(M**10*xh**5 + M**8*Q**2*xh**4*(-3*zh + xh*(2 + 3*zh)) + M**6*Q**4*xh**3*(5*zh**2 - 4*xh*zh*(1 + 2*zh) + xh**2*(2 + 4*zh + 4*zh**2)) + \
            M**4*Q**6*xh**2*(-5*zh**3 - 2*xh**2*zh*(2 + 2*zh + 5*zh**2) + xh*zh*(2 + 2*zh + 11*zh**2) + 2*xh**3*(1 + 3*zh**2 + zh**3)) + \
            Q**10*(-1 + xh)*zh*(-2*xh*zh**2*(1 - zh + zh**2) + (1 - zh + zh**2)**2 + xh**4*(1 - 2*zh + 2*zh**2) - 2*xh**3*zh*(1 - 2*zh + 2*zh**2) + xh**2*zh*(2 - zh + 2*zh**3)) + \
            M**2*Q**8*xh*(1 - 2*zh + 3*zh**2 - 2*zh**3 + 3*zh**4 + 2*xh*zh**2*(-2 + zh - 4*zh**2) - 4*xh**3*zh*(1 - zh + 2*zh**2 + zh**3) + xh**4*(1 + 4*zh**3) + xh**2*zh*(2 + 3*zh + 10*zh**3))))/ \
            (M*Q**6*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*(-1 + zh)*zh**2)
         FUUT3s1g = (32*xh**2*(Q**6*(-1 + xh)*xh*(-1 + zh)*zh + M**6*xh**2*(1 - zh + zh**2) + \
            M**2*Q**4*(1 + (-2 + 2*xh - 5*xh**2)*zh + (3 - 6*xh + 10*xh**2)*zh**2 + (-2 + 4*xh - 6*xh**2)*zh**3 + (1 - 2*xh + 2*xh**2)*zh**4) + \
            M**4*Q**2*xh*(zh*(-3 + 3*zh - 2*zh**2) + xh*(-3 + 7*zh - 4*zh**2 + 2*zh**3))))/(27.*M*Q**2*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*zh**2)
         FUUT3pjg = (24*xh*(7*M**16*xh**8*zh + M**14*Q**2*xh**7*(8 + (-5 + 23*xh)*zh + (-31 + 21*xh)*zh**2) + \
            M**12*Q**4*xh**6*(8 + 3*(-11 + 18*xh + 7*xh**2)*zh + (19 - 126*xh + 71*xh**2)*zh**2 + 4*(15 - 19*xh + 7*xh**2)*zh**3 + 2*zh**4) + \
            M**10*Q**6*xh**5*(20 - 63*zh + 87*zh**2 - 50*zh**3 - 56*zh**4 - 8*zh**5 + xh**3*zh*(-3 + 67*zh + 110*zh**2 + 14*zh**3) - xh**2*(48 - 147*zh + 77*zh**2 + 262*zh**3 + 92*zh**4) + \
            xh*(-8 + 87*zh - 259*zh**2 + 304*zh**3 + 112*zh**4 + 4*zh**5)) + \
            Q**16*(-1 + xh)**2*zh**2*(xh**6*(3 - 6*zh + 6*zh**2) + xh**5*(4 - 20*zh + 34*zh**2 - 24*zh**3) + 3*(-1 + 2*zh - 2*zh**2 + zh**3)**2 + \
            xh**4*(3 - 14*zh + 44*zh**2 - 60*zh**3 + 36*zh**4) - 2*xh**3*(2 - 7*zh + 3*zh**2 + 11*zh**3 - 18*zh**4 + 12*zh**5) - \
            2*xh*(-2 + 7*zh - 10*zh**2 + 4*zh**3 + 4*zh**4 - 6*zh**5 + 3*zh**6) + xh**2*(3 + 4*zh - 28*zh**2 + 36*zh**3 - 20*zh**4 + 2*zh**5 + 6*zh**6)) + \
            M**8*Q**8*xh**4*(20 - 135*zh + 255*zh**2 - 224*zh**3 + 114*zh**4 + 14*zh**5 + 12*zh**6 + xh**4*zh*(-7 - 19*zh + 154*zh**2 + 62*zh**3) - \
            2*xh**3*(32 - 34*zh - 89*zh**2 + 120*zh**3 + 161*zh**4 + 20*zh**5) + 2*xh**2*(-20 + 134*zh - 203*zh**2 + 15*zh**3 + 226*zh**4 + 58*zh**5 + 2*zh**6) - \
            2*xh*(-10 - 28*zh + 164*zh**2 - 262*zh**3 + 207*zh**4 + 38*zh**5 + 6*zh**6)) + \
            M**6*Q**10*xh**3*(4 - 79*zh + 305*zh**2 - 446*zh**3 + 318*zh**4 - 134*zh**5 + 12*zh**6 - 8*zh**7 + xh**5*zh*(9 - 61*zh + 76*zh**2 + 108*zh**3) - \
            xh**4*(24 + 63*zh - 199*zh**2 - 84*zh**3 + 392*zh**4 + 144*zh**5) + xh*zh*(-35 - 169*zh + 560*zh**2 - 614*zh**3 + 356*zh**4 + 18*zh**5 + 12*zh**6) + \
            2*xh**3*(-12 + 68*zh + 42*zh**2 - 297*zh**3 + 195*zh**4 + 188*zh**5 + 24*zh**6) + xh**2*(12 + 88*zh - 322*zh**2 + 224*zh**3 + 252*zh**4 - 482*zh**5 - 64*zh**6 - 8*zh**7)) + \
            M**2*Q**14*xh*zh*(-5 + 23*zh - 60*zh**2 + 104*zh**3 - 98*zh**4 + 48*zh**5 - 11*zh**6 - zh**7 + xh**7*(3 + 5*zh - 26*zh**2 + 38*zh**3) + \
            xh**6*(1 - 43*zh + 66*zh**2 + 20*zh**3 - 112*zh**4) + xh**5*(-7 + 39*zh + 58*zh**2 - 222*zh**3 + 132*zh**4 + 112*zh**5) + \
            xh**4*(3 + 7*zh - 74*zh**2 + 16*zh**3 + 210*zh**4 - 232*zh**5 - 40*zh**6) + xh**3*(-1 + 17*zh - 146*zh**2 + 296*zh**3 - 268*zh**4 + 54*zh**5 + 110*zh**6 + 2*zh**7) + \
            xh*(-3 + 3*zh + 20*zh**2 - 142*zh**3 + 222*zh**4 - 154*zh**5 + 53*zh**6 + 3*zh**7) - xh**2*(-9 + 51*zh - 162*zh**2 + 110*zh**3 + 86*zh**4 - 172*zh**5 + 112*zh**6 + 4*zh**7)) + \
            M**4*Q**12*xh**2*(4 - 23*zh + 107*zh**2 - 272*zh**3 + 322*zh**4 - 192*zh**5 + 64*zh**6 - 5*zh**7 + 2*zh**8 + xh**6*zh*(11 - 23*zh - 16*zh**2 + 92*zh**3) - \
            2*xh**5*zh*(21 + zh - 106*zh**2 + 82*zh**3 + 96*zh**4) + xh**4*zh*(-19 + 253*zh - 342*zh**2 - 156*zh**3 + 408*zh**4 + 120*zh**5) - \
            2*xh**3*(-6 + 14*zh - 53*zh**2 + 223*zh**3 - 343*zh**4 + 139*zh**5 + 124*zh**6 + 12*zh**7) - \
            2*xh*(-2 + 7*zh - 10*zh**2 - 113*zh**3 + 252*zh**4 - 211*zh**5 + 92*zh**6 + zh**7 + 2*zh**8) + \
            xh**2*(12 - 61*zh - 29*zh**2 + 76*zh**3 + 156*zh**4 - 378*zh**5 + 316*zh**6 + 16*zh**7 + 4*zh**8))))/ \
            (M**3*Q**6*(Q**2*(-1 + xh) + M**2*xh)**3*(M**2*xh + Q**2*(1 + xh - zh))**4*(-1 + zh)*zh**3)

         FUUT1s0q = (-4*xh*(M**4*xh**2 + 2*M**2*Q**2*(-1 + xh)*xh*zh + Q**4*(1 + zh**2 - 2*xh*(1 - zh + zh**2) + xh**2*(1 - 2*zh + 2*zh**2))))/(3.*M*Q**4*(M**2*xh + Q**2*(1 + xh - zh))**2*(-1 + zh))
         FUUT3s1q = (-2*(M**6*xh**3 - M**4*Q**2*xh**2*(2 + xh*(-2 + zh) + zh) + M**2*Q**4*xh*(-4*xh*zh + (1 + zh)**2 + xh**2*(-1 + 2*zh)) + Q**6*(-1 + xh)*zh*(1 + zh**2 - 2*xh*zh**2 + xh**2*(1 - 2*zh + 2*zh**2))))/ \
            (9.*M**3*Q**4*(-1 + xh)**2*(M**2*xh - Q**2*zh)**2)  
         FUUT3pjq = (-16*xh*(7*M**8*xh**4 + 2*M**6*Q**2*xh**3*(1 + xh - 8*zh + 7*xh*zh) + 2*M**4*Q**4*xh**2*(-3 + 3*zh + 7*zh**2 + xh*(6 - 6*zh - 9*zh**2) + xh**2*(-3 + 3*zh + 7*zh**2)) + \
            2*M**2*Q**6*xh*(1 + 4*zh - zh**2 - 4*zh**3 + xh**3*(1 - 7*zh + 10*zh**2) - xh**2*(9 - 22*zh + 16*zh**2 + 2*zh**3) + xh*(-9 + 5*zh + zh**2 + 5*zh**3)) + \
            Q**8*(3*(-1 + zh)**2*(1 + zh**2) + xh**4*(3 - 6*zh + 6*zh**2) - 2*xh**3*(2 - 5*zh**2 + 6*zh**3) - 2*xh*(2 - 8*zh + 7*zh**2 - 4*zh**3 + 3*zh**4) + 2*xh**2*(1 - 2*zh + zh**3 + 3*zh**4))))/ \
            (3.*M**3*(M**2*Q*xh + Q**3*(1 + xh - zh))**4*(-1 + zh))
         
         FUUTg = ((FUUT1s0g*ldme.j1s0CO_var(ldme_scale)) + (FUUT3s1g*((15/8)*ldme.j3s1CO_var(ldme_scale) + ldme.j3s1CS)) + (FUUT3pjg*ldme.j3p0CO_var(ldme_scale)))*pdf.gluon_PDF
         FUUTq = ((FUUT1s0q*ldme.j1s0CO_var(ldme_scale)) + (FUUT3pjq*ldme.j3p0CO_var(ldme_scale)))*pdf.qqb_PDF + (FUUT3s1q*ldme.j3s1CO_var(ldme_scale))*pdf.qqb_PDF_weight
         
         return e_c**2*alpha_QED*pdf.als**2*(FUUTg + FUUTq)

      ####---- Longitudinal ----#### 
      def FUUL(xh, zh, pdf, ldme, ldme_scale):
         FUUL1s0g = (6*xh**2*(2*M**8*xh**4 + Q**8*(-1 + xh)**2*zh**2*(1 + 2*xh**2 + xh*(2 - 4*zh) - 2*zh + 2*zh**2) + 2*M**6*Q**2*xh**3*(-3*zh + 2*xh*(1 + zh)) + \
            2*M**2*Q**6*(-1 + xh)*xh*zh*(2*xh**2*(1 + zh) + zh*(-1 + 3*zh) + xh*(1 - 4*zh - 2*zh**2)) + M**4*Q**4*xh**2*(-2*xh*zh*(4 + 5*zh) + zh*(-2 + 9*zh) + 2*xh**2*(1 + 4*zh + zh**2))))/ \
            (M*Q**4*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*zh**2)
         FUUL3s1g = (16*xh**2*(Q**4*(-1 + xh)**2*zh**2 + 2*M**2*Q**2*(-1 + xh)*xh*zh*(-1 + 6*zh - 6*zh**2 + 2*zh**3) + M**4*xh**2*(-2 + 10*zh - 11*zh**2 + 4*zh**3)))/ \
            (27.*M*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*zh**2)
         FUUL3pjg = (48*xh**2*(M**14*xh**7*zh*(-5 + 7*zh) + M**12*Q**2*xh**6*(-4 - 3*(-5 + 7*xh)*zh + (7 + 19*xh)*zh**2 + 2*(-13 + 7*xh)*zh**3) + \
            Q**14*(-1 + xh)**3*(-1 + zh)*zh**3*(2 + 3*xh**4 + xh**3*(8 - 12*zh) - 6*zh + 11*zh**2 - 8*zh**3 + 3*zh**4 + xh**2*(11 - 24*zh + 18*zh**2) + xh*(6 - 22*zh + 24*zh**2 - 12*zh**3)) + \
            M**10*Q**4*xh**5*(xh**2*zh*(-34 + 2*zh + 55*zh**2 + 7*zh**3) - 2*xh*(6 - 18*zh - 31*zh**2 + 40*zh**3 + 23*zh**4) + 2*(-2 + 15*zh - 34*zh**2 + 9*zh**3 + 19*zh**4 + zh**5)) + \
            M**2*Q**12*(-1 + xh)**2*xh*zh**2*(xh**4*(-4 - 13*zh + 19*zh**2) - 2*xh**3*(7 - 9*zh - 24*zh**2 + 28*zh**3) + xh**2*(-16 + 49*zh - 37*zh**2 - 46*zh**3 + 56*zh**4) + \
            zh*(14 - 50*zh + 69*zh**2 - 47*zh**3 + 15*zh**4 + zh**5) - 2*xh*(5 - 31*zh + 52*zh**2 - 36*zh**3 + 2*zh**4 + 10*zh**5)) + \
            M**8*Q**6*xh**4*(xh**3*zh*(-26 - 42*zh + 77*zh**2 + 31*zh**3) - xh**2*(12 - 18*zh - 138*zh**2 + 43*zh**3 + 161*zh**4 + 20*zh**5) + \
            2*xh*(-4 + 29*zh - 48*zh**2 - 55*zh**3 + 80*zh**4 + 29*zh**5 + zh**6) - 2*(2 - 17*zh + 53*zh**2 - 79*zh**3 + 38*zh**4 + 10*zh**5 + 3*zh**6)) + \
            M**4*Q**10*xh**2*zh*(7 - 61*zh + 172*zh**2 - 236*zh**3 + 180*zh**4 - 72*zh**5 + 4*zh**6 - 2*zh**7 + xh**5*(-1 - 25*zh - 8*zh**2 + 46*zh**3) + \
            xh**4*(-9 + 23*zh + 124*zh**2 - 82*zh**3 - 96*zh**4) + 2*xh**3*(-1 + 28*zh - 79*zh**2 - 48*zh**3 + 102*zh**4 + 30*zh**5) - \
            2*xh**2*(1 - 21*zh + 72*zh**2 - 144*zh**3 + 54*zh**4 + 62*zh**5 + 6*zh**6) + xh*(7 - 35*zh + 14*zh**2 + 80*zh**3 - 180*zh**4 + 136*zh**5 + 8*zh**6 + 2*zh**7)) + \
            M**6*Q**8*xh**3*(-4 + 39*zh - 139*zh**2 + 246*zh**3 - 244*zh**4 + 116*zh**5 - 4*zh**6 + 6*zh**7 + xh**4*zh*(-9 - 53*zh + 38*zh**2 + 54*zh**3) - \
            4*xh**3*(1 + 3*zh - 28*zh**2 - 23*zh**3 + 49*zh**4 + 18*zh**5) + 2*xh**2*(-2 + 13*zh + 10*zh**2 - 141*zh**3 + 62*zh**4 + 94*zh**5 + 12*zh**6) - \
            2*xh*(2 - 16*zh + 36*zh**2 - 19*zh**3 - 87*zh**4 + 96*zh**5 + 16*zh**6 + 2*zh**7))))/(M**3*(Q**2*(-1 + xh) + M**2*xh)**3*(M**2*Q*xh + Q**3*(1 + xh - zh))**4*(-1 + zh)*zh**3)

         FUUL1s0q = (-8*xh**2*(M**2*xh + Q**2*(-1 + xh)*zh))/(3.*M*(M**2*Q*xh + Q**3*(1 + xh - zh))**2)
         FUUL3s1q = (-4*Q**2*xh*(-1 + zh)*zh**2)/(9.*M**3*(M**2*xh - Q**2*zh)**2)
         FUUL3pjq = (-32*xh**2*(M**6*xh**3*(-5 + 7*zh) + Q**6*(-1 + xh)*(-1 + zh)*zh*(3 + 3*xh**2 + xh*(2 - 6*zh) - 2*zh + 3*zh**2) + M**4*Q**2*xh**2*(2 + 3*zh - 9*zh**2 + xh*(-6 + 3*zh + 7*zh**2)) + \
            M**2*Q**4*xh*(7 - 13*zh + 5*zh**2 + 5*zh**3 + xh**2*(-1 - 7*zh + 10*zh**2) - 2*xh*(-1 - 6*zh + 8*zh**2 + zh**3))))/(3.*M**3*Q**2*(M**2*xh + Q**2*(1 + xh - zh))**4*(-1 + zh))

         FUULg = ((FUUL1s0g*ldme.j1s0CO_var(ldme_scale)) + (FUUL3s1g*((15/8)*ldme.j3s1CO_var(ldme_scale) + ldme.j3s1CS)) + (FUUL3pjg*ldme.j3p0CO_var(ldme_scale)))*pdf.gluon_PDF
         FUULq = ((FUUL1s0q*ldme.j1s0CO_var(ldme_scale)) + (FUUL3pjq*ldme.j3p0CO_var(ldme_scale)))*pdf.qqb_PDF + (FUUL3s1q*ldme.j3s1CO_var(ldme_scale))*pdf.qqb_PDF_weight
         
         return e_c**2*alpha_QED*pdf.als**2*(FUULg + FUULq)

      ####---- Cos2phi ----#### 
      def FUUcos2phi(xh, zh, pdf, ldme, ldme_scale):
         FUUcos2phi1s0g = 4*(3*xh*(2*M**8*xh**5*(-1 + zh) + M**6*Q**2*xh**3*((-2 + zh)*zh - 6*xh*(-1 + zh)*zh + 4*xh**2*(-1 + zh**2)) + \
            M**4*Q**4*xh**2*(2*xh**2*zh*(4 + zh - 5*zh**2) + zh*(2 + 2*zh - zh**2) + xh*zh*(-2 - 10*zh + 9*zh**2) + 2*xh**3*(-1 - 3*zh + 3*zh**2 + zh**3)) + \
            Q**8*(-1 + xh)**2*zh**2*(2*xh**3*(-1 + zh) + zh + xh**2*(-2 + 6*zh - 4*zh**2) + xh*(-2 + 3*zh - 4*zh**2 + 2*zh**3)) + \
            M**2*Q**6*(-1 + xh)*xh*zh*(3*zh + 4*xh**3*(-1 + zh**2) - 2*xh**2*(1 - 5*zh + 2*zh**2 + 2*zh**3) + xh*(-2 + zh - 8*zh**2 + 6*zh**3))))/ \
            (M*Q**4*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*(-1 + zh)*zh**2)
         FUUcos2phi3s1g = 4*(-16*M*xh**3*(-1 + zh)*(M**2*xh + Q**2*(-1 + xh)*zh)*(M**2 + Q**2*(-1 + 4*zh - 2*zh**2)))/ \
            (27.*Q**2*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*zh**2)
         FUUcos2phi3pjg = 4*(12*xh*(2*M**14*xh**7*(-1 + zh)*(-4 + (3 + 7*xh)*zh) + M**12*Q**2*xh**6*(8 - 38*zh + 50*zh**2 - 17*zh**3 + 2*xh**2*zh*(-31 + 17*zh + 14*zh**2) + \
            xh*(16 - 8*zh + 38*zh**2 - 46*zh**3)) + M**10*Q**4*xh**5*(8 - 38*zh + 68*zh**2 - 65*zh**3 + 15*zh**4 + 2*xh**3*zh*(-54 - 8*zh + 55*zh**2 + 7*zh**3) - \
            2*xh**2*zh*(-18 - 91*zh + 63*zh**2 + 46*zh**3) + xh*(8 - 10*zh - 62*zh**2 + 19*zh**3 + 55*zh**4 + 4*zh**5)) + \
            Q**14*(-1 + xh)**3*zh**3*(6*xh**5*(-1 + zh) + (-1 + zh)**2*zh - 8*xh**4*(2 - 5*zh + 3*zh**2) + xh**3*(-22 + 69*zh - 84*zh**2 + 36*zh**3) + \
            xh**2*(-16 + 59*zh - 90*zh**2 + 72*zh**3 - 24*zh**4) + xh*(-4 + 21*zh - 38*zh**2 + 37*zh**3 - 22*zh**4 + 6*zh**5)) + \
            M**2*Q**12*(-1 + xh)**2*xh*zh**2*(2*xh**5*(-6 - 13*zh + 19*zh**2) + xh**4*(-38 + 54*zh + 96*zh**2 - 112*zh**3) + zh*(7 - 15*zh + 11*zh**2 - 3*zh**3) + \
            xh**3*(-44 + 131*zh - 109*zh**2 - 92*zh**3 + 112*zh**4) - xh**2*(38 - 151*zh + 239*zh**2 - 174*zh**3 + 8*zh**4 + 40*zh**5) + \
            xh*(-12 + 67*zh - 133*zh**2 + 151*zh**3 - 103*zh**4 + 30*zh**5 + 2*zh**6)) + \
            M**8*Q**6*xh**4*(8 - 42*zh + 92*zh**2 - 86*zh**3 + 49*zh**4 - 3*zh**5 + 2*xh**4*zh*(-46 - 62*zh + 77*zh**2 + 31*zh**3) - \
            2*xh*zh*(-22 + 105*zh - 167*zh**2 + 92*zh**3 + 8*zh**4 + 6*zh**5) - 2*xh**3*(8 - 8*zh - 182*zh**2 + zh**3 + 161*zh**4 + 20*zh**5) + \
            xh**2*(-8 + 94*zh - 186*zh**2 - 244*zh**3 + 249*zh**4 + 116*zh**5 + 4*zh**6)) + \
            M**4*Q**10*xh**2*zh*(2*xh**6*(-3 - 35*zh - 8*zh**2 + 46*zh**3) + zh*(-2 - 11*zh + 36*zh**2 - 33*zh**3 + 13*zh**4) - \
            2*xh**5*(12 - 31*zh - 159*zh**2 + 82*zh**3 + 96*zh**4) + xh**4*zh*(128 - 363*zh - 288*zh**2 + 408*zh**3 + 120*zh**4) - \
            2*xh**3*zh*(-23 + 145*zh - 327*zh**2 + 81*zh**3 + 124*zh**4 + 12*zh**5) - \
            4*xh*(-4 + 26*zh - 70*zh**2 + 110*zh**3 - 93*zh**4 + 37*zh**5 - 2*zh**6 + zh**7) + xh**2*(14 - 60*zh + 82*zh**2 + 110*zh**3 - 393*zh**4 + 263*zh**5 + 16*zh**6 + 4*zh**7)) + \
            M**6*Q**8*xh**3*(2*xh**5*zh*(-19 - 73*zh + 38*zh**2 + 54*zh**3) - zh*(12 - 58*zh + 108*zh**2 - 84*zh**3 + 33*zh**4 + zh**5) - \
            2*xh**4*(4 + 19*zh - 145*zh**2 - 146*zh**3 + 196*zh**4 + 72*zh**5) + 2*xh**3*(-4 + 33*zh + 17*zh**2 - 316*zh**3 + 68*zh**4 + 188*zh**5 + 24*zh**6) - \
            2*xh**2*(4 - 41*zh + 122*zh**2 - 87*zh**3 - 173*zh**4 + 167*zh**5 + 32*zh**6 + 4*zh**7) + \
            xh*(-8 + 68*zh - 252*zh**2 + 534*zh**3 - 554*zh**4 + 265*zh**5 - 17*zh**6 + 12*zh**7))))/ \
            (M**3*(Q**2*(-1 + xh) + M**2*xh)**3*(M**2*Q*xh + Q**3*(1 + xh - zh))**4*(-1 + zh)*zh**3)

         FUUcos2phi1s0q = 4*(-4*xh*(1 + xh*(-1 + zh))*(M**2*xh + Q**2*(-1 + xh)*zh))/(3.*M*(M**2*Q*xh + Q**3*(1 + xh - zh))**2*(-1 + zh))
         FUUcos2phi3s1q = 4*(2*xh*(-1 + zh)*(M**2 - Q**2*zh)*(M**2*xh + Q**2*(-1 + xh)*zh))/(9.*M**3*(-1 + xh)*(M**2*Q*xh - Q**3*zh)**2)
         FUUcos2phi3pjq = 4*(-16*xh*(M**2*xh + Q**2*(-1 + xh)*zh)*(M**4*xh**2*(-3 + 7*xh*(-1 + zh)) + 2*M**2*Q**2*xh*(-1 + 5*xh**2*(-1 + zh) + zh - xh*zh**2) + \
            Q**4*(3*xh**3*(-1 + zh) + (-1 + zh)**2 + xh**2*(-1 + 8*zh - 6*zh**2) + xh*(3 - zh - 5*zh**2 + 3*zh**3))))/ \
            (3.*M**3*Q**2*(M**2*xh + Q**2*(1 + xh - zh))**4*(-1 + zh))
         
         FUUcos2phig = ((FUUcos2phi1s0g*ldme.j1s0CO_var(ldme_scale)) + (FUUcos2phi3s1g*((15/8)*ldme.j3s1CO_var(ldme_scale) + ldme.j3s1CS)) + (FUUcos2phi3pjg *ldme.j3p0CO_var(ldme_scale)))*pdf.gluon_PDF
         FUUcos2phiq = ((FUUcos2phi1s0q*ldme.j1s0CO_var(ldme_scale)) + (FUUcos2phi3pjq*ldme.j3p0CO_var(ldme_scale)))*pdf.qqb_PDF + (FUUcos2phi3s1q*ldme.j3s1CO_var(ldme_scale))*pdf.qqb_PDF_weight
         
         return e_c**2*alpha_QED*pdf.als**2*(FUUcos2phig + FUUcos2phiq)
#### ---------------------------------------------- ####


#### ---------------------------------------------- ####
#### ---- 1S0 CO ----####
elif wave == 0: 
   ####---- (gluon) contribution ----####
   if iq == 0:
      ####---- Transverse ----#### 
      def FUUT(xh, zh, pdf, ldme, ldme_scale):
         FUUT1s0g = (6*xh*(M**10*xh**5 + M**8*Q**2*xh**4*(-3*zh + xh*(2 + 3*zh)) + M**6*Q**4*xh**3*(5*zh**2 - 4*xh*zh*(1 + 2*zh) + xh**2*(2 + 4*zh + 4*zh**2)) + \
            M**4*Q**6*xh**2*(-5*zh**3 - 2*xh**2*zh*(2 + 2*zh + 5*zh**2) + xh*zh*(2 + 2*zh + 11*zh**2) + 2*xh**3*(1 + 3*zh**2 + zh**3)) + \
            Q**10*(-1 + xh)*zh*(-2*xh*zh**2*(1 - zh + zh**2) + (1 - zh + zh**2)**2 + xh**4*(1 - 2*zh + 2*zh**2) - 2*xh**3*zh*(1 - 2*zh + 2*zh**2) + xh**2*zh*(2 - zh + 2*zh**3)) + \
            M**2*Q**8*xh*(1 - 2*zh + 3*zh**2 - 2*zh**3 + 3*zh**4 + 2*xh*zh**2*(-2 + zh - 4*zh**2) - 4*xh**3*zh*(1 - zh + 2*zh**2 + zh**3) + xh**4*(1 + 4*zh**3) + xh**2*zh*(2 + 3*zh + 10*zh**3))))/ \
            (M*Q**6*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*(-1 + zh)*zh**2)
        
         FUUTg = FUUT1s0g*pdf.gluon_PDF

         return e_c**2*alpha_QED*pdf.als**2*FUUTg

      ####---- Longitudinal ----#### 
      def FUUL(xh, zh, pdf, ldme, ldme_scale):
         FUUL1s0g = (6*xh**2*(2*M**8*xh**4 + Q**8*(-1 + xh)**2*zh**2*(1 + 2*xh**2 + xh*(2 - 4*zh) - 2*zh + 2*zh**2) + 2*M**6*Q**2*xh**3*(-3*zh + 2*xh*(1 + zh)) + \
            2*M**2*Q**6*(-1 + xh)*xh*zh*(2*xh**2*(1 + zh) + zh*(-1 + 3*zh) + xh*(1 - 4*zh - 2*zh**2)) + M**4*Q**4*xh**2*(-2*xh*zh*(4 + 5*zh) + zh*(-2 + 9*zh) + 2*xh**2*(1 + 4*zh + zh**2))))/ \
            (M*Q**4*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*zh**2)
         
         FUULg = FUUL1s0g*pdf.gluon_PDF

         return e_c**2*alpha_QED*pdf.als**2*FUULg
      
      ####---- Cos2phi ----#### 
      def FUUcos2phi(xh, zh, pdf, ldme, ldme_scale):
         FUUcos2phi1s0g = 4*(3*xh*(2*M**8*xh**5*(-1 + zh) + M**6*Q**2*xh**3*((-2 + zh)*zh - 6*xh*(-1 + zh)*zh + 4*xh**2*(-1 + zh**2)) + \
            M**4*Q**4*xh**2*(2*xh**2*zh*(4 + zh - 5*zh**2) + zh*(2 + 2*zh - zh**2) + xh*zh*(-2 - 10*zh + 9*zh**2) + 2*xh**3*(-1 - 3*zh + 3*zh**2 + zh**3)) + \
            Q**8*(-1 + xh)**2*zh**2*(2*xh**3*(-1 + zh) + zh + xh**2*(-2 + 6*zh - 4*zh**2) + xh*(-2 + 3*zh - 4*zh**2 + 2*zh**3)) + \
            M**2*Q**6*(-1 + xh)*xh*zh*(3*zh + 4*xh**3*(-1 + zh**2) - 2*xh**2*(1 - 5*zh + 2*zh**2 + 2*zh**3) + xh*(-2 + zh - 8*zh**2 + 6*zh**3))))/ \
            (M*Q**4*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*(-1 + zh)*zh**2)
         
         FUUcos2phig = FUUcos2phi1s0g*pdf.gluon_PDF

         return e_c**2*alpha_QED*pdf.als**2*FUUcos2phig
      
   ####---- (quark) contribution ----####     
   elif ig == 0:  
      ####---- Transverse ----####    
      def FUUT(xh, zh, pdf, ldme, ldme_scale):
         FUUT1s0q = (-4*xh*(M**4*xh**2 + 2*M**2*Q**2*(-1 + xh)*xh*zh + Q**4*(1 + zh**2 - 2*xh*(1 - zh + zh**2) + xh**2*(1 - 2*zh + 2*zh**2))))/(3.*M*Q**4*(M**2*xh + Q**2*(1 + xh - zh))**2*(-1 + zh))
            
         FUUTq = (FUUT1s0q)*pdf.qqb_PDF

         return e_c**2*alpha_QED*pdf.als**2*FUUTq

      ####---- Longitudinal ----#### 
      def FUUL(xh, zh, pdf, ldme, ldme_scale):
         FUUL1s0q = (-8*xh**2*(M**2*xh + Q**2*(-1 + xh)*zh))/(3.*M*(M**2*Q*xh + Q**3*(1 + xh - zh))**2)
         
         FUULq = (FUUL1s0q)*pdf.qqb_PDF

         return e_c**2*alpha_QED*pdf.als**2*FUULq
      
      ####---- Cos2phi ----#### 
      def FUUcos2phi(xh, zh, pdf, ldme, ldme_scale):         
         FUUcos2phi1s0q = 4*(-4*xh*(1 + xh*(-1 + zh))*(M**2*xh + Q**2*(-1 + xh)*zh))/(3.*M*(M**2*Q*xh + Q**3*(1 + xh - zh))**2*(-1 + zh))
         
         FUUcos2phiq = (FUUcos2phi1s0q)*pdf.qqb_PDF

         return e_c**2*alpha_QED*pdf.als**2*FUUcos2phiq

   ####----(gluon + quark) contribution ----####
   else:          
      ####---- Transverse ----#### 
      def FUUT(xh, zh, pdf, ldme, ldme_scale):
         FUUT1s0g = (6*xh*(M**10*xh**5 + M**8*Q**2*xh**4*(-3*zh + xh*(2 + 3*zh)) + M**6*Q**4*xh**3*(5*zh**2 - 4*xh*zh*(1 + 2*zh) + xh**2*(2 + 4*zh + 4*zh**2)) + \
            M**4*Q**6*xh**2*(-5*zh**3 - 2*xh**2*zh*(2 + 2*zh + 5*zh**2) + xh*zh*(2 + 2*zh + 11*zh**2) + 2*xh**3*(1 + 3*zh**2 + zh**3)) + \
            Q**10*(-1 + xh)*zh*(-2*xh*zh**2*(1 - zh + zh**2) + (1 - zh + zh**2)**2 + xh**4*(1 - 2*zh + 2*zh**2) - 2*xh**3*zh*(1 - 2*zh + 2*zh**2) + xh**2*zh*(2 - zh + 2*zh**3)) + \
            M**2*Q**8*xh*(1 - 2*zh + 3*zh**2 - 2*zh**3 + 3*zh**4 + 2*xh*zh**2*(-2 + zh - 4*zh**2) - 4*xh**3*zh*(1 - zh + 2*zh**2 + zh**3) + xh**4*(1 + 4*zh**3) + xh**2*zh*(2 + 3*zh + 10*zh**3))))/ \
            (M*Q**6*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*(-1 + zh)*zh**2)
         
         FUUT1s0q = (-4*xh*(M**4*xh**2 + 2*M**2*Q**2*(-1 + xh)*xh*zh + Q**4*(1 + zh**2 - 2*xh*(1 - zh + zh**2) + xh**2*(1 - 2*zh + 2*zh**2))))/(3.*M*Q**4*(M**2*xh + Q**2*(1 + xh - zh))**2*(-1 + zh))
         
         FUUTg = (FUUT1s0g)*pdf.gluon_PDF
         FUUTq = (FUUT1s0q)*pdf.qqb_PDF

         return e_c**2*alpha_QED*pdf.als**2*(FUUTg + FUUTq)

      ####---- Longitudinal ----#### 
      def FUUL(xh, zh, pdf, ldme, ldme_scale):
         FUUL1s0g = (6*xh**2*(2*M**8*xh**4 + Q**8*(-1 + xh)**2*zh**2*(1 + 2*xh**2 + xh*(2 - 4*zh) - 2*zh + 2*zh**2) + 2*M**6*Q**2*xh**3*(-3*zh + 2*xh*(1 + zh)) + \
            2*M**2*Q**6*(-1 + xh)*xh*zh*(2*xh**2*(1 + zh) + zh*(-1 + 3*zh) + xh*(1 - 4*zh - 2*zh**2)) + M**4*Q**4*xh**2*(-2*xh*zh*(4 + 5*zh) + zh*(-2 + 9*zh) + 2*xh**2*(1 + 4*zh + zh**2))))/ \
            (M*Q**4*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*zh**2)

         FUUL1s0q = (-8*xh**2*(M**2*xh + Q**2*(-1 + xh)*zh))/(3.*M*(M**2*Q*xh + Q**3*(1 + xh - zh))**2)
         
         FUULg = (FUUL1s0g)*pdf.gluon_PDF
         FUULq = (FUUL1s0q)*pdf.qqb_PDF

         return e_c**2*alpha_QED*pdf.als**2*(FUULg + FUULq)
      
      ####---- Cos2phi ----#### 
      def FUUcos2phi(xh, zh, pdf, ldme, ldme_scale):
         FUUcos2phi1s0g = 4*(3*xh*(2*M**8*xh**5*(-1 + zh) + M**6*Q**2*xh**3*((-2 + zh)*zh - 6*xh*(-1 + zh)*zh + 4*xh**2*(-1 + zh**2)) + \
            M**4*Q**4*xh**2*(2*xh**2*zh*(4 + zh - 5*zh**2) + zh*(2 + 2*zh - zh**2) + xh*zh*(-2 - 10*zh + 9*zh**2) + 2*xh**3*(-1 - 3*zh + 3*zh**2 + zh**3)) + \
            Q**8*(-1 + xh)**2*zh**2*(2*xh**3*(-1 + zh) + zh + xh**2*(-2 + 6*zh - 4*zh**2) + xh*(-2 + 3*zh - 4*zh**2 + 2*zh**3)) + \
            M**2*Q**6*(-1 + xh)*xh*zh*(3*zh + 4*xh**3*(-1 + zh**2) - 2*xh**2*(1 - 5*zh + 2*zh**2 + 2*zh**3) + xh*(-2 + zh - 8*zh**2 + 6*zh**3))))/ \
            (M*Q**4*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*(-1 + zh)*zh**2)
         
         FUUcos2phi1s0q = 4*(-4*xh*(1 + xh*(-1 + zh))*(M**2*xh + Q**2*(-1 + xh)*zh))/(3.*M*(M**2*Q*xh + Q**3*(1 + xh - zh))**2*(-1 + zh))
         
         FUUcos2phig = (FUUcos2phi1s0g)*pdf.gluon_PDF
         FUUcos2phiq = (FUUcos2phi1s0q)*pdf.qqb_PDF

         return e_c**2*alpha_QED*pdf.als**2*(FUUcos2phig + FUUcos2phiq)


#### ---------------------------------------------- ####
#### ---- 3s1  ----####
if wave == 1: 
   ####---- (gluon) contribution ----####
   if iq == 0:    
      ####---- Transverse ----####   
      def FUUT(xh, zh, pdf, ldme, ldme_scale):
         FUUT3s1g = (32*xh**2*(Q**6*(-1 + xh)*xh*(-1 + zh)*zh + M**6*xh**2*(1 - zh + zh**2) + \
            M**2*Q**4*(1 + (-2 + 2*xh - 5*xh**2)*zh + (3 - 6*xh + 10*xh**2)*zh**2 + (-2 + 4*xh - 6*xh**2)*zh**3 + (1 - 2*xh + 2*xh**2)*zh**4) + \
            M**4*Q**2*xh*(zh*(-3 + 3*zh - 2*zh**2) + xh*(-3 + 7*zh - 4*zh**2 + 2*zh**3))))/(27.*M*Q**2*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*zh**2)
         
         FUUTg = (FUUT3s1g)*pdf.gluon_PDF

         return e_c**2*alpha_QED*pdf.als**2*FUUTg

      ####---- Longitudinal ----#### 
      def FUUL(xh, zh, pdf, ldme, ldme_scale):
         FUUL3s1g = (16*xh**2*(Q**4*(-1 + xh)**2*zh**2 + 2*M**2*Q**2*(-1 + xh)*xh*zh*(-1 + 6*zh - 6*zh**2 + 2*zh**3) + M**4*xh**2*(-2 + 10*zh - 11*zh**2 + 4*zh**3)))/ \
            (27.*M*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*zh**2)

         FUULg = (FUUL3s1g)*pdf.gluon_PDF

         return e_c**2*alpha_QED*pdf.als**2*FUULg
      
      ####---- Cos2phi ----####
      def FUUcos2phi(xh, zh, pdf, ldme, ldme_scale):
         FUUcos2phi3s1g = 4*(-16*M*xh**3*(-1 + zh)*(M**2*xh + Q**2*(-1 + xh)*zh)*(M**2 + Q**2*(-1 + 4*zh - 2*zh**2)))/ \
            (27.*Q**2*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*zh**2)

         FUUcos2phig = (FUUcos2phi3s1g)*pdf.gluon_PDF

         return e_c**2*alpha_QED*pdf.als**2*FUUcos2phig

   ####---- (quark) contribution ----####    
   elif ig == 0:  
      ####---- Transverse ----#### 
      def FUUT(xh, zh, pdf, ldme, ldme_scale):
         FUUT3s1q = (-2*(M**6*xh**3 - M**4*Q**2*xh**2*(2 + xh*(-2 + zh) + zh) + M**2*Q**4*xh*(-4*xh*zh + (1 + zh)**2 + xh**2*(-1 + 2*zh)) + Q**6*(-1 + xh)*zh*(1 + zh**2 - 2*xh*zh**2 + xh**2*(1 - 2*zh + 2*zh**2))))/ \
            (9.*M**3*Q**4*(-1 + xh)**2*(M**2*xh - Q**2*zh)**2)
         
         FUUTq = (FUUT3s1q)*pdf.qqb_PDF_weight

         return e_c**2*alpha_QED*pdf.als**2*FUUTq

      ####---- Longitudinal ----#### 
      def FUUL(xh, zh, pdf, ldme, ldme_scale):
         FUUL3s1q = (-4*Q**2*xh*(-1 + zh)*zh**2)/(9.*M**3*(M**2*xh - Q**2*zh)**2)

         FUULq = (FUUL3s1q)*pdf.qqb_PDF_weight

         return e_c**2*alpha_QED*pdf.als**2*FUULq
      
      ####---- Cos2phi ----####
      def FUUcos2phi(xh, zh, pdf, ldme, ldme_scale):
         FUUcos2phi3s1q = 4*(2*xh*(-1 + zh)*(M**2 - Q**2*zh)*(M**2*xh + Q**2*(-1 + xh)*zh))/(9.*M**3*(-1 + xh)*(M**2*Q*xh - Q**3*zh)**2)

         FUUcos2phiq = (FUUcos2phi3s1q)*pdf.qqb_PDF_weight

         return e_c**2*alpha_QED*pdf.als**2*FUUcos2phiq

   ####----(gluon + quark) contribution ----####
   else:        
      ####---- Transverse ----####     
      def FUUT(xh, zh, pdf, ldme, ldme_scale):
         FUUT3s1g = (32*xh**2*(Q**6*(-1 + xh)*xh*(-1 + zh)*zh + M**6*xh**2*(1 - zh + zh**2) + \
            M**2*Q**4*(1 + (-2 + 2*xh - 5*xh**2)*zh + (3 - 6*xh + 10*xh**2)*zh**2 + (-2 + 4*xh - 6*xh**2)*zh**3 + (1 - 2*xh + 2*xh**2)*zh**4) + \
            M**4*Q**2*xh*(zh*(-3 + 3*zh - 2*zh**2) + xh*(-3 + 7*zh - 4*zh**2 + 2*zh**3))))/(27.*M*Q**2*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*zh**2)
         
         FUUT3s1q = (-2*(M**6*xh**3 - M**4*Q**2*xh**2*(2 + xh*(-2 + zh) + zh) + M**2*Q**4*xh*(-4*xh*zh + (1 + zh)**2 + xh**2*(-1 + 2*zh)) + Q**6*(-1 + xh)*zh*(1 + zh**2 - 2*xh*zh**2 + xh**2*(1 - 2*zh + 2*zh**2))))/ \
            (9.*M**3*Q**4*(-1 + xh)**2*(M**2*xh - Q**2*zh)**2)
         
         FUUTg = (FUUT3s1g)*pdf.gluon_PDF
         FUUTq = (FUUT3s1q)*pdf.qqb_PDF_weight

         return e_c**2*alpha_QED*pdf.als**2*(FUUTg + FUUTq)

      ####---- Longitudinal ----#### 
      def FUUL(xh, zh, pdf, ldme, ldme_scale):
         FUUL3s1g = (16*xh**2*(Q**4*(-1 + xh)**2*zh**2 + 2*M**2*Q**2*(-1 + xh)*xh*zh*(-1 + 6*zh - 6*zh**2 + 2*zh**3) + M**4*xh**2*(-2 + 10*zh - 11*zh**2 + 4*zh**3)))/ \
            (27.*M*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*zh**2)

         FUUL3s1q = (-4*Q**2*xh*(-1 + zh)*zh**2)/(9.*M**3*(M**2*xh - Q**2*zh)**2)

         FUULg = (FUUL3s1g)*pdf.gluon_PDF
         FUULq = (FUUL3s1q)*pdf.qqb_PDF_weight

         return e_c**2*alpha_QED*pdf.als**2*(FUULg + FUULq)
      
      ####---- Cos2phi ----#### 
      def FUUcos2phi(xh, zh, pdf, ldme, ldme_scale):
         FUUcos2phi3s1g = 4*(-16*M*xh**3*(-1 + zh)*(M**2*xh + Q**2*(-1 + xh)*zh)*(M**2 + Q**2*(-1 + 4*zh - 2*zh**2)))/ \
            (27.*Q**2*(Q**2*(-1 + xh) + M**2*xh)**2*(M**2*xh + Q**2*(1 + xh - zh))**2*zh**2)
         
         FUUcos2phi3s1q = 4*(2*xh*(-1 + zh)*(M**2 - Q**2*zh)*(M**2*xh + Q**2*(-1 + xh)*zh))/(9.*M**3*(-1 + xh)*(M**2*Q*xh - Q**3*zh)**2)
         
         FUUcos2phig = (FUUcos2phi3s1g)*pdf.gluon_PDF
         FUUcos2phiq = (FUUcos2phi3s1q)*pdf.qqb_PDF_weight

         return e_c**2*alpha_QED*pdf.als**2*(FUUcos2phig + FUUcos2phiq)


#### ---------------------------------------------- ####
#### ---- 3pJ CO ----####
if wave == 2: 
   ####---- (gluon) contribution ----####
   if iq == 0:  
      ####---- Transverse ----####     
      def FUUT(xh, zh, pdf, ldme, ldme_scale):
         FUUT3pjg = (24*xh*(7*M**16*xh**8*zh + M**14*Q**2*xh**7*(8 + (-5 + 23*xh)*zh + (-31 + 21*xh)*zh**2) + \
            M**12*Q**4*xh**6*(8 + 3*(-11 + 18*xh + 7*xh**2)*zh + (19 - 126*xh + 71*xh**2)*zh**2 + 4*(15 - 19*xh + 7*xh**2)*zh**3 + 2*zh**4) + \
            M**10*Q**6*xh**5*(20 - 63*zh + 87*zh**2 - 50*zh**3 - 56*zh**4 - 8*zh**5 + xh**3*zh*(-3 + 67*zh + 110*zh**2 + 14*zh**3) - xh**2*(48 - 147*zh + 77*zh**2 + 262*zh**3 + 92*zh**4) + \
            xh*(-8 + 87*zh - 259*zh**2 + 304*zh**3 + 112*zh**4 + 4*zh**5)) + \
            Q**16*(-1 + xh)**2*zh**2*(xh**6*(3 - 6*zh + 6*zh**2) + xh**5*(4 - 20*zh + 34*zh**2 - 24*zh**3) + 3*(-1 + 2*zh - 2*zh**2 + zh**3)**2 + \
            xh**4*(3 - 14*zh + 44*zh**2 - 60*zh**3 + 36*zh**4) - 2*xh**3*(2 - 7*zh + 3*zh**2 + 11*zh**3 - 18*zh**4 + 12*zh**5) - \
            2*xh*(-2 + 7*zh - 10*zh**2 + 4*zh**3 + 4*zh**4 - 6*zh**5 + 3*zh**6) + xh**2*(3 + 4*zh - 28*zh**2 + 36*zh**3 - 20*zh**4 + 2*zh**5 + 6*zh**6)) + \
            M**8*Q**8*xh**4*(20 - 135*zh + 255*zh**2 - 224*zh**3 + 114*zh**4 + 14*zh**5 + 12*zh**6 + xh**4*zh*(-7 - 19*zh + 154*zh**2 + 62*zh**3) - \
            2*xh**3*(32 - 34*zh - 89*zh**2 + 120*zh**3 + 161*zh**4 + 20*zh**5) + 2*xh**2*(-20 + 134*zh - 203*zh**2 + 15*zh**3 + 226*zh**4 + 58*zh**5 + 2*zh**6) - \
            2*xh*(-10 - 28*zh + 164*zh**2 - 262*zh**3 + 207*zh**4 + 38*zh**5 + 6*zh**6)) + \
            M**6*Q**10*xh**3*(4 - 79*zh + 305*zh**2 - 446*zh**3 + 318*zh**4 - 134*zh**5 + 12*zh**6 - 8*zh**7 + xh**5*zh*(9 - 61*zh + 76*zh**2 + 108*zh**3) - \
            xh**4*(24 + 63*zh - 199*zh**2 - 84*zh**3 + 392*zh**4 + 144*zh**5) + xh*zh*(-35 - 169*zh + 560*zh**2 - 614*zh**3 + 356*zh**4 + 18*zh**5 + 12*zh**6) + \
            2*xh**3*(-12 + 68*zh + 42*zh**2 - 297*zh**3 + 195*zh**4 + 188*zh**5 + 24*zh**6) + xh**2*(12 + 88*zh - 322*zh**2 + 224*zh**3 + 252*zh**4 - 482*zh**5 - 64*zh**6 - 8*zh**7)) + \
            M**2*Q**14*xh*zh*(-5 + 23*zh - 60*zh**2 + 104*zh**3 - 98*zh**4 + 48*zh**5 - 11*zh**6 - zh**7 + xh**7*(3 + 5*zh - 26*zh**2 + 38*zh**3) + \
            xh**6*(1 - 43*zh + 66*zh**2 + 20*zh**3 - 112*zh**4) + xh**5*(-7 + 39*zh + 58*zh**2 - 222*zh**3 + 132*zh**4 + 112*zh**5) + \
            xh**4*(3 + 7*zh - 74*zh**2 + 16*zh**3 + 210*zh**4 - 232*zh**5 - 40*zh**6) + xh**3*(-1 + 17*zh - 146*zh**2 + 296*zh**3 - 268*zh**4 + 54*zh**5 + 110*zh**6 + 2*zh**7) + \
            xh*(-3 + 3*zh + 20*zh**2 - 142*zh**3 + 222*zh**4 - 154*zh**5 + 53*zh**6 + 3*zh**7) - xh**2*(-9 + 51*zh - 162*zh**2 + 110*zh**3 + 86*zh**4 - 172*zh**5 + 112*zh**6 + 4*zh**7)) + \
            M**4*Q**12*xh**2*(4 - 23*zh + 107*zh**2 - 272*zh**3 + 322*zh**4 - 192*zh**5 + 64*zh**6 - 5*zh**7 + 2*zh**8 + xh**6*zh*(11 - 23*zh - 16*zh**2 + 92*zh**3) - \
            2*xh**5*zh*(21 + zh - 106*zh**2 + 82*zh**3 + 96*zh**4) + xh**4*zh*(-19 + 253*zh - 342*zh**2 - 156*zh**3 + 408*zh**4 + 120*zh**5) - \
            2*xh**3*(-6 + 14*zh - 53*zh**2 + 223*zh**3 - 343*zh**4 + 139*zh**5 + 124*zh**6 + 12*zh**7) - \
            2*xh*(-2 + 7*zh - 10*zh**2 - 113*zh**3 + 252*zh**4 - 211*zh**5 + 92*zh**6 + zh**7 + 2*zh**8) + \
            xh**2*(12 - 61*zh - 29*zh**2 + 76*zh**3 + 156*zh**4 - 378*zh**5 + 316*zh**6 + 16*zh**7 + 4*zh**8))))/ \
            (M**3*Q**6*(Q**2*(-1 + xh) + M**2*xh)**3*(M**2*xh + Q**2*(1 + xh - zh))**4*(-1 + zh)*zh**3)

         FUUTg = (FUUT3pjg)*pdf.gluon_PDF

         return e_c**2*alpha_QED*pdf.als**2*FUUTg

      ####---- Longitudinal ----#### 
      def FUUL(xh, zh, pdf, ldme, ldme_scale):
         FUUL3pjg = (48*xh**2*(M**14*xh**7*zh*(-5 + 7*zh) + M**12*Q**2*xh**6*(-4 - 3*(-5 + 7*xh)*zh + (7 + 19*xh)*zh**2 + 2*(-13 + 7*xh)*zh**3) + \
            Q**14*(-1 + xh)**3*(-1 + zh)*zh**3*(2 + 3*xh**4 + xh**3*(8 - 12*zh) - 6*zh + 11*zh**2 - 8*zh**3 + 3*zh**4 + xh**2*(11 - 24*zh + 18*zh**2) + xh*(6 - 22*zh + 24*zh**2 - 12*zh**3)) + \
            M**10*Q**4*xh**5*(xh**2*zh*(-34 + 2*zh + 55*zh**2 + 7*zh**3) - 2*xh*(6 - 18*zh - 31*zh**2 + 40*zh**3 + 23*zh**4) + 2*(-2 + 15*zh - 34*zh**2 + 9*zh**3 + 19*zh**4 + zh**5)) + \
            M**2*Q**12*(-1 + xh)**2*xh*zh**2*(xh**4*(-4 - 13*zh + 19*zh**2) - 2*xh**3*(7 - 9*zh - 24*zh**2 + 28*zh**3) + xh**2*(-16 + 49*zh - 37*zh**2 - 46*zh**3 + 56*zh**4) + \
            zh*(14 - 50*zh + 69*zh**2 - 47*zh**3 + 15*zh**4 + zh**5) - 2*xh*(5 - 31*zh + 52*zh**2 - 36*zh**3 + 2*zh**4 + 10*zh**5)) + \
            M**8*Q**6*xh**4*(xh**3*zh*(-26 - 42*zh + 77*zh**2 + 31*zh**3) - xh**2*(12 - 18*zh - 138*zh**2 + 43*zh**3 + 161*zh**4 + 20*zh**5) + \
            2*xh*(-4 + 29*zh - 48*zh**2 - 55*zh**3 + 80*zh**4 + 29*zh**5 + zh**6) - 2*(2 - 17*zh + 53*zh**2 - 79*zh**3 + 38*zh**4 + 10*zh**5 + 3*zh**6)) + \
            M**4*Q**10*xh**2*zh*(7 - 61*zh + 172*zh**2 - 236*zh**3 + 180*zh**4 - 72*zh**5 + 4*zh**6 - 2*zh**7 + xh**5*(-1 - 25*zh - 8*zh**2 + 46*zh**3) + \
            xh**4*(-9 + 23*zh + 124*zh**2 - 82*zh**3 - 96*zh**4) + 2*xh**3*(-1 + 28*zh - 79*zh**2 - 48*zh**3 + 102*zh**4 + 30*zh**5) - \
            2*xh**2*(1 - 21*zh + 72*zh**2 - 144*zh**3 + 54*zh**4 + 62*zh**5 + 6*zh**6) + xh*(7 - 35*zh + 14*zh**2 + 80*zh**3 - 180*zh**4 + 136*zh**5 + 8*zh**6 + 2*zh**7)) + \
            M**6*Q**8*xh**3*(-4 + 39*zh - 139*zh**2 + 246*zh**3 - 244*zh**4 + 116*zh**5 - 4*zh**6 + 6*zh**7 + xh**4*zh*(-9 - 53*zh + 38*zh**2 + 54*zh**3) - \
            4*xh**3*(1 + 3*zh - 28*zh**2 - 23*zh**3 + 49*zh**4 + 18*zh**5) + 2*xh**2*(-2 + 13*zh + 10*zh**2 - 141*zh**3 + 62*zh**4 + 94*zh**5 + 12*zh**6) - \
            2*xh*(2 - 16*zh + 36*zh**2 - 19*zh**3 - 87*zh**4 + 96*zh**5 + 16*zh**6 + 2*zh**7))))/(M**3*(Q**2*(-1 + xh) + M**2*xh)**3*(M**2*Q*xh + Q**3*(1 + xh - zh))**4*(-1 + zh)*zh**3)

         FUULg = (FUUL3pjg)*pdf.gluon_PDF

         return e_c**2*alpha_QED*pdf.als**2*FUULg
      
      ####---- Cos2phi ----####
      def FUUcos2phi(xh, zh, pdf, ldme, ldme_scale):
         FUUcos2phi3pjg = 4*(12*xh*(2*M**14*xh**7*(-1 + zh)*(-4 + (3 + 7*xh)*zh) + M**12*Q**2*xh**6*(8 - 38*zh + 50*zh**2 - 17*zh**3 + 2*xh**2*zh*(-31 + 17*zh + 14*zh**2) + \
            xh*(16 - 8*zh + 38*zh**2 - 46*zh**3)) + M**10*Q**4*xh**5*(8 - 38*zh + 68*zh**2 - 65*zh**3 + 15*zh**4 + 2*xh**3*zh*(-54 - 8*zh + 55*zh**2 + 7*zh**3) - \
            2*xh**2*zh*(-18 - 91*zh + 63*zh**2 + 46*zh**3) + xh*(8 - 10*zh - 62*zh**2 + 19*zh**3 + 55*zh**4 + 4*zh**5)) + \
            Q**14*(-1 + xh)**3*zh**3*(6*xh**5*(-1 + zh) + (-1 + zh)**2*zh - 8*xh**4*(2 - 5*zh + 3*zh**2) + xh**3*(-22 + 69*zh - 84*zh**2 + 36*zh**3) + \
            xh**2*(-16 + 59*zh - 90*zh**2 + 72*zh**3 - 24*zh**4) + xh*(-4 + 21*zh - 38*zh**2 + 37*zh**3 - 22*zh**4 + 6*zh**5)) + \
            M**2*Q**12*(-1 + xh)**2*xh*zh**2*(2*xh**5*(-6 - 13*zh + 19*zh**2) + xh**4*(-38 + 54*zh + 96*zh**2 - 112*zh**3) + zh*(7 - 15*zh + 11*zh**2 - 3*zh**3) + \
            xh**3*(-44 + 131*zh - 109*zh**2 - 92*zh**3 + 112*zh**4) - xh**2*(38 - 151*zh + 239*zh**2 - 174*zh**3 + 8*zh**4 + 40*zh**5) + \
            xh*(-12 + 67*zh - 133*zh**2 + 151*zh**3 - 103*zh**4 + 30*zh**5 + 2*zh**6)) + \
            M**8*Q**6*xh**4*(8 - 42*zh + 92*zh**2 - 86*zh**3 + 49*zh**4 - 3*zh**5 + 2*xh**4*zh*(-46 - 62*zh + 77*zh**2 + 31*zh**3) - \
            2*xh*zh*(-22 + 105*zh - 167*zh**2 + 92*zh**3 + 8*zh**4 + 6*zh**5) - 2*xh**3*(8 - 8*zh - 182*zh**2 + zh**3 + 161*zh**4 + 20*zh**5) + \
            xh**2*(-8 + 94*zh - 186*zh**2 - 244*zh**3 + 249*zh**4 + 116*zh**5 + 4*zh**6)) + \
            M**4*Q**10*xh**2*zh*(2*xh**6*(-3 - 35*zh - 8*zh**2 + 46*zh**3) + zh*(-2 - 11*zh + 36*zh**2 - 33*zh**3 + 13*zh**4) - \
            2*xh**5*(12 - 31*zh - 159*zh**2 + 82*zh**3 + 96*zh**4) + xh**4*zh*(128 - 363*zh - 288*zh**2 + 408*zh**3 + 120*zh**4) - \
            2*xh**3*zh*(-23 + 145*zh - 327*zh**2 + 81*zh**3 + 124*zh**4 + 12*zh**5) - \
            4*xh*(-4 + 26*zh - 70*zh**2 + 110*zh**3 - 93*zh**4 + 37*zh**5 - 2*zh**6 + zh**7) + xh**2*(14 - 60*zh + 82*zh**2 + 110*zh**3 - 393*zh**4 + 263*zh**5 + 16*zh**6 + 4*zh**7)) + \
            M**6*Q**8*xh**3*(2*xh**5*zh*(-19 - 73*zh + 38*zh**2 + 54*zh**3) - zh*(12 - 58*zh + 108*zh**2 - 84*zh**3 + 33*zh**4 + zh**5) - \
            2*xh**4*(4 + 19*zh - 145*zh**2 - 146*zh**3 + 196*zh**4 + 72*zh**5) + 2*xh**3*(-4 + 33*zh + 17*zh**2 - 316*zh**3 + 68*zh**4 + 188*zh**5 + 24*zh**6) - \
            2*xh**2*(4 - 41*zh + 122*zh**2 - 87*zh**3 - 173*zh**4 + 167*zh**5 + 32*zh**6 + 4*zh**7) + \
            xh*(-8 + 68*zh - 252*zh**2 + 534*zh**3 - 554*zh**4 + 265*zh**5 - 17*zh**6 + 12*zh**7))))/ \
            (M**3*(Q**2*(-1 + xh) + M**2*xh)**3*(M**2*Q*xh + Q**3*(1 + xh - zh))**4*(-1 + zh)*zh**3)

         FUUcos2phig = (FUUcos2phi3pjg)*pdf.gluon_PDF

         return e_c**2*alpha_QED*pdf.als**2*FUUcos2phig

   ####---- (quark) contribution ----####   
   elif ig == 0: 
      ####---- Transverse ----####  
      def FUUT(xh, zh, pdf, ldme, ldme_scale):
         FUUT3pjq = (-16*xh*(7*M**8*xh**4 + 2*M**6*Q**2*xh**3*(1 + xh - 8*zh + 7*xh*zh) + 2*M**4*Q**4*xh**2*(-3 + 3*zh + 7*zh**2 + xh*(6 - 6*zh - 9*zh**2) + xh**2*(-3 + 3*zh + 7*zh**2)) + \
            2*M**2*Q**6*xh*(1 + 4*zh - zh**2 - 4*zh**3 + xh**3*(1 - 7*zh + 10*zh**2) - xh**2*(9 - 22*zh + 16*zh**2 + 2*zh**3) + xh*(-9 + 5*zh + zh**2 + 5*zh**3)) + \
            Q**8*(3*(-1 + zh)**2*(1 + zh**2) + xh**4*(3 - 6*zh + 6*zh**2) - 2*xh**3*(2 - 5*zh**2 + 6*zh**3) - 2*xh*(2 - 8*zh + 7*zh**2 - 4*zh**3 + 3*zh**4) + 2*xh**2*(1 - 2*zh + zh**3 + 3*zh**4))))/ \
            (3.*M**3*(M**2*Q*xh + Q**3*(1 + xh - zh))**4*(-1 + zh))
         
         FUUTq = (FUUT3pjq)*pdf.qqb_PDF

         return e_c**2*alpha_QED*pdf.als**2*FUUTq

      ####---- Longitudinal ----#### 
      def FUUL(xh, zh, pdf, ldme, ldme_scale):
         FUUL3pjq = (-32*xh**2*(M**6*xh**3*(-5 + 7*zh) + Q**6*(-1 + xh)*(-1 + zh)*zh*(3 + 3*xh**2 + xh*(2 - 6*zh) - 2*zh + 3*zh**2) + M**4*Q**2*xh**2*(2 + 3*zh - 9*zh**2 + xh*(-6 + 3*zh + 7*zh**2)) + \
            M**2*Q**4*xh*(7 - 13*zh + 5*zh**2 + 5*zh**3 + xh**2*(-1 - 7*zh + 10*zh**2) - 2*xh*(-1 - 6*zh + 8*zh**2 + zh**3))))/(3.*M**3*Q**2*(M**2*xh + Q**2*(1 + xh - zh))**4*(-1 + zh))
         
         FUULq = (FUUL3pjq)*pdf.qqb_PDF

         return e_c**2*alpha_QED*pdf.als**2*FUULq
      
      ####---- Cos2phi ----####
      def FUUcos2phi(xh, zh, pdf, ldme, ldme_scale):
         FUUcos2phi3pjq = 4*(-16*xh*(M**2*xh + Q**2*(-1 + xh)*zh)*(M**4*xh**2*(-3 + 7*xh*(-1 + zh)) + 2*M**2*Q**2*xh*(-1 + 5*xh**2*(-1 + zh) + zh - xh*zh**2) + \
            Q**4*(3*xh**3*(-1 + zh) + (-1 + zh)**2 + xh**2*(-1 + 8*zh - 6*zh**2) + xh*(3 - zh - 5*zh**2 + 3*zh**3))))/ \
            (3.*M**3*Q**2*(M**2*xh + Q**2*(1 + xh - zh))**4*(-1 + zh))

         FUUcos2phiq = (FUUcos2phi3pjq)*pdf.qqb_PDF

         return e_c**2*alpha_QED*pdf.als**2*FUUcos2phiq

   ####----(gluon + quark) contribution ----####
   else:         
      ####---- Transverse ----####   
      def FUUT(xh, zh, pdf, ldme, ldme_scale):
         FUUT3pjg = (24*xh*(7*M**16*xh**8*zh + M**14*Q**2*xh**7*(8 + (-5 + 23*xh)*zh + (-31 + 21*xh)*zh**2) + \
            M**12*Q**4*xh**6*(8 + 3*(-11 + 18*xh + 7*xh**2)*zh + (19 - 126*xh + 71*xh**2)*zh**2 + 4*(15 - 19*xh + 7*xh**2)*zh**3 + 2*zh**4) + \
            M**10*Q**6*xh**5*(20 - 63*zh + 87*zh**2 - 50*zh**3 - 56*zh**4 - 8*zh**5 + xh**3*zh*(-3 + 67*zh + 110*zh**2 + 14*zh**3) - xh**2*(48 - 147*zh + 77*zh**2 + 262*zh**3 + 92*zh**4) + \
            xh*(-8 + 87*zh - 259*zh**2 + 304*zh**3 + 112*zh**4 + 4*zh**5)) + \
            Q**16*(-1 + xh)**2*zh**2*(xh**6*(3 - 6*zh + 6*zh**2) + xh**5*(4 - 20*zh + 34*zh**2 - 24*zh**3) + 3*(-1 + 2*zh - 2*zh**2 + zh**3)**2 + \
            xh**4*(3 - 14*zh + 44*zh**2 - 60*zh**3 + 36*zh**4) - 2*xh**3*(2 - 7*zh + 3*zh**2 + 11*zh**3 - 18*zh**4 + 12*zh**5) - \
            2*xh*(-2 + 7*zh - 10*zh**2 + 4*zh**3 + 4*zh**4 - 6*zh**5 + 3*zh**6) + xh**2*(3 + 4*zh - 28*zh**2 + 36*zh**3 - 20*zh**4 + 2*zh**5 + 6*zh**6)) + \
            M**8*Q**8*xh**4*(20 - 135*zh + 255*zh**2 - 224*zh**3 + 114*zh**4 + 14*zh**5 + 12*zh**6 + xh**4*zh*(-7 - 19*zh + 154*zh**2 + 62*zh**3) - \
            2*xh**3*(32 - 34*zh - 89*zh**2 + 120*zh**3 + 161*zh**4 + 20*zh**5) + 2*xh**2*(-20 + 134*zh - 203*zh**2 + 15*zh**3 + 226*zh**4 + 58*zh**5 + 2*zh**6) - \
            2*xh*(-10 - 28*zh + 164*zh**2 - 262*zh**3 + 207*zh**4 + 38*zh**5 + 6*zh**6)) + \
            M**6*Q**10*xh**3*(4 - 79*zh + 305*zh**2 - 446*zh**3 + 318*zh**4 - 134*zh**5 + 12*zh**6 - 8*zh**7 + xh**5*zh*(9 - 61*zh + 76*zh**2 + 108*zh**3) - \
            xh**4*(24 + 63*zh - 199*zh**2 - 84*zh**3 + 392*zh**4 + 144*zh**5) + xh*zh*(-35 - 169*zh + 560*zh**2 - 614*zh**3 + 356*zh**4 + 18*zh**5 + 12*zh**6) + \
            2*xh**3*(-12 + 68*zh + 42*zh**2 - 297*zh**3 + 195*zh**4 + 188*zh**5 + 24*zh**6) + xh**2*(12 + 88*zh - 322*zh**2 + 224*zh**3 + 252*zh**4 - 482*zh**5 - 64*zh**6 - 8*zh**7)) + \
            M**2*Q**14*xh*zh*(-5 + 23*zh - 60*zh**2 + 104*zh**3 - 98*zh**4 + 48*zh**5 - 11*zh**6 - zh**7 + xh**7*(3 + 5*zh - 26*zh**2 + 38*zh**3) + \
            xh**6*(1 - 43*zh + 66*zh**2 + 20*zh**3 - 112*zh**4) + xh**5*(-7 + 39*zh + 58*zh**2 - 222*zh**3 + 132*zh**4 + 112*zh**5) + \
            xh**4*(3 + 7*zh - 74*zh**2 + 16*zh**3 + 210*zh**4 - 232*zh**5 - 40*zh**6) + xh**3*(-1 + 17*zh - 146*zh**2 + 296*zh**3 - 268*zh**4 + 54*zh**5 + 110*zh**6 + 2*zh**7) + \
            xh*(-3 + 3*zh + 20*zh**2 - 142*zh**3 + 222*zh**4 - 154*zh**5 + 53*zh**6 + 3*zh**7) - xh**2*(-9 + 51*zh - 162*zh**2 + 110*zh**3 + 86*zh**4 - 172*zh**5 + 112*zh**6 + 4*zh**7)) + \
            M**4*Q**12*xh**2*(4 - 23*zh + 107*zh**2 - 272*zh**3 + 322*zh**4 - 192*zh**5 + 64*zh**6 - 5*zh**7 + 2*zh**8 + xh**6*zh*(11 - 23*zh - 16*zh**2 + 92*zh**3) - \
            2*xh**5*zh*(21 + zh - 106*zh**2 + 82*zh**3 + 96*zh**4) + xh**4*zh*(-19 + 253*zh - 342*zh**2 - 156*zh**3 + 408*zh**4 + 120*zh**5) - \
            2*xh**3*(-6 + 14*zh - 53*zh**2 + 223*zh**3 - 343*zh**4 + 139*zh**5 + 124*zh**6 + 12*zh**7) - \
            2*xh*(-2 + 7*zh - 10*zh**2 - 113*zh**3 + 252*zh**4 - 211*zh**5 + 92*zh**6 + zh**7 + 2*zh**8) + \
            xh**2*(12 - 61*zh - 29*zh**2 + 76*zh**3 + 156*zh**4 - 378*zh**5 + 316*zh**6 + 16*zh**7 + 4*zh**8))))/ \
            (M**3*Q**6*(Q**2*(-1 + xh) + M**2*xh)**3*(M**2*xh + Q**2*(1 + xh - zh))**4*(-1 + zh)*zh**3)

         FUUT3pjq = (-16*xh*(7*M**8*xh**4 + 2*M**6*Q**2*xh**3*(1 + xh - 8*zh + 7*xh*zh) + 2*M**4*Q**4*xh**2*(-3 + 3*zh + 7*zh**2 + xh*(6 - 6*zh - 9*zh**2) + xh**2*(-3 + 3*zh + 7*zh**2)) + \
            2*M**2*Q**6*xh*(1 + 4*zh - zh**2 - 4*zh**3 + xh**3*(1 - 7*zh + 10*zh**2) - xh**2*(9 - 22*zh + 16*zh**2 + 2*zh**3) + xh*(-9 + 5*zh + zh**2 + 5*zh**3)) + \
            Q**8*(3*(-1 + zh)**2*(1 + zh**2) + xh**4*(3 - 6*zh + 6*zh**2) - 2*xh**3*(2 - 5*zh**2 + 6*zh**3) - 2*xh*(2 - 8*zh + 7*zh**2 - 4*zh**3 + 3*zh**4) + 2*xh**2*(1 - 2*zh + zh**3 + 3*zh**4))))/ \
            (3.*M**3*(M**2*Q*xh + Q**3*(1 + xh - zh))**4*(-1 + zh))
         
         FUUTg = (FUUT3pjg)*pdf.gluon_PDF
         FUUTq = (FUUT3pjq)*pdf.qqb_PDF

         return e_c**2*alpha_QED*pdf.als**2*(FUUTg + FUUTq)

      ####---- Longitudinal ----#### 
      def FUUL(xh, zh, pdf, ldme, ldme_scale):
         FUUL3pjg = (48*xh**2*(M**14*xh**7*zh*(-5 + 7*zh) + M**12*Q**2*xh**6*(-4 - 3*(-5 + 7*xh)*zh + (7 + 19*xh)*zh**2 + 2*(-13 + 7*xh)*zh**3) + \
            Q**14*(-1 + xh)**3*(-1 + zh)*zh**3*(2 + 3*xh**4 + xh**3*(8 - 12*zh) - 6*zh + 11*zh**2 - 8*zh**3 + 3*zh**4 + xh**2*(11 - 24*zh + 18*zh**2) + xh*(6 - 22*zh + 24*zh**2 - 12*zh**3)) + \
            M**10*Q**4*xh**5*(xh**2*zh*(-34 + 2*zh + 55*zh**2 + 7*zh**3) - 2*xh*(6 - 18*zh - 31*zh**2 + 40*zh**3 + 23*zh**4) + 2*(-2 + 15*zh - 34*zh**2 + 9*zh**3 + 19*zh**4 + zh**5)) + \
            M**2*Q**12*(-1 + xh)**2*xh*zh**2*(xh**4*(-4 - 13*zh + 19*zh**2) - 2*xh**3*(7 - 9*zh - 24*zh**2 + 28*zh**3) + xh**2*(-16 + 49*zh - 37*zh**2 - 46*zh**3 + 56*zh**4) + \
            zh*(14 - 50*zh + 69*zh**2 - 47*zh**3 + 15*zh**4 + zh**5) - 2*xh*(5 - 31*zh + 52*zh**2 - 36*zh**3 + 2*zh**4 + 10*zh**5)) + \
            M**8*Q**6*xh**4*(xh**3*zh*(-26 - 42*zh + 77*zh**2 + 31*zh**3) - xh**2*(12 - 18*zh - 138*zh**2 + 43*zh**3 + 161*zh**4 + 20*zh**5) + \
            2*xh*(-4 + 29*zh - 48*zh**2 - 55*zh**3 + 80*zh**4 + 29*zh**5 + zh**6) - 2*(2 - 17*zh + 53*zh**2 - 79*zh**3 + 38*zh**4 + 10*zh**5 + 3*zh**6)) + \
            M**4*Q**10*xh**2*zh*(7 - 61*zh + 172*zh**2 - 236*zh**3 + 180*zh**4 - 72*zh**5 + 4*zh**6 - 2*zh**7 + xh**5*(-1 - 25*zh - 8*zh**2 + 46*zh**3) + \
            xh**4*(-9 + 23*zh + 124*zh**2 - 82*zh**3 - 96*zh**4) + 2*xh**3*(-1 + 28*zh - 79*zh**2 - 48*zh**3 + 102*zh**4 + 30*zh**5) - \
            2*xh**2*(1 - 21*zh + 72*zh**2 - 144*zh**3 + 54*zh**4 + 62*zh**5 + 6*zh**6) + xh*(7 - 35*zh + 14*zh**2 + 80*zh**3 - 180*zh**4 + 136*zh**5 + 8*zh**6 + 2*zh**7)) + \
            M**6*Q**8*xh**3*(-4 + 39*zh - 139*zh**2 + 246*zh**3 - 244*zh**4 + 116*zh**5 - 4*zh**6 + 6*zh**7 + xh**4*zh*(-9 - 53*zh + 38*zh**2 + 54*zh**3) - \
            4*xh**3*(1 + 3*zh - 28*zh**2 - 23*zh**3 + 49*zh**4 + 18*zh**5) + 2*xh**2*(-2 + 13*zh + 10*zh**2 - 141*zh**3 + 62*zh**4 + 94*zh**5 + 12*zh**6) - \
            2*xh*(2 - 16*zh + 36*zh**2 - 19*zh**3 - 87*zh**4 + 96*zh**5 + 16*zh**6 + 2*zh**7))))/(M**3*(Q**2*(-1 + xh) + M**2*xh)**3*(M**2*Q*xh + Q**3*(1 + xh - zh))**4*(-1 + zh)*zh**3)

         FUUL3pjq = (-32*xh**2*(M**6*xh**3*(-5 + 7*zh) + Q**6*(-1 + xh)*(-1 + zh)*zh*(3 + 3*xh**2 + xh*(2 - 6*zh) - 2*zh + 3*zh**2) + M**4*Q**2*xh**2*(2 + 3*zh - 9*zh**2 + xh*(-6 + 3*zh + 7*zh**2)) + \
            M**2*Q**4*xh*(7 - 13*zh + 5*zh**2 + 5*zh**3 + xh**2*(-1 - 7*zh + 10*zh**2) - 2*xh*(-1 - 6*zh + 8*zh**2 + zh**3))))/(3.*M**3*Q**2*(M**2*xh + Q**2*(1 + xh - zh))**4*(-1 + zh))
         
         FUULg = (FUUL3pjg)*pdf.gluon_PDF
         FUULq = (FUUL3pjq)*pdf.qqb_PDF

         return e_c**2*alpha_QED*pdf.als**2*(FUULg + FUULq)
      
      ####---- Cos2phi ----####
      def FUUcos2phi(xh, zh, pdf, ldme, ldme_scale):
         FUUcos2phi3pjg = 4*(12*xh*(2*M**14*xh**7*(-1 + zh)*(-4 + (3 + 7*xh)*zh) + M**12*Q**2*xh**6*(8 - 38*zh + 50*zh**2 - 17*zh**3 + 2*xh**2*zh*(-31 + 17*zh + 14*zh**2) + \
            xh*(16 - 8*zh + 38*zh**2 - 46*zh**3)) + M**10*Q**4*xh**5*(8 - 38*zh + 68*zh**2 - 65*zh**3 + 15*zh**4 + 2*xh**3*zh*(-54 - 8*zh + 55*zh**2 + 7*zh**3) - \
            2*xh**2*zh*(-18 - 91*zh + 63*zh**2 + 46*zh**3) + xh*(8 - 10*zh - 62*zh**2 + 19*zh**3 + 55*zh**4 + 4*zh**5)) + \
            Q**14*(-1 + xh)**3*zh**3*(6*xh**5*(-1 + zh) + (-1 + zh)**2*zh - 8*xh**4*(2 - 5*zh + 3*zh**2) + xh**3*(-22 + 69*zh - 84*zh**2 + 36*zh**3) + \
            xh**2*(-16 + 59*zh - 90*zh**2 + 72*zh**3 - 24*zh**4) + xh*(-4 + 21*zh - 38*zh**2 + 37*zh**3 - 22*zh**4 + 6*zh**5)) + \
            M**2*Q**12*(-1 + xh)**2*xh*zh**2*(2*xh**5*(-6 - 13*zh + 19*zh**2) + xh**4*(-38 + 54*zh + 96*zh**2 - 112*zh**3) + zh*(7 - 15*zh + 11*zh**2 - 3*zh**3) + \
            xh**3*(-44 + 131*zh - 109*zh**2 - 92*zh**3 + 112*zh**4) - xh**2*(38 - 151*zh + 239*zh**2 - 174*zh**3 + 8*zh**4 + 40*zh**5) + \
            xh*(-12 + 67*zh - 133*zh**2 + 151*zh**3 - 103*zh**4 + 30*zh**5 + 2*zh**6)) + \
            M**8*Q**6*xh**4*(8 - 42*zh + 92*zh**2 - 86*zh**3 + 49*zh**4 - 3*zh**5 + 2*xh**4*zh*(-46 - 62*zh + 77*zh**2 + 31*zh**3) - \
            2*xh*zh*(-22 + 105*zh - 167*zh**2 + 92*zh**3 + 8*zh**4 + 6*zh**5) - 2*xh**3*(8 - 8*zh - 182*zh**2 + zh**3 + 161*zh**4 + 20*zh**5) + \
            xh**2*(-8 + 94*zh - 186*zh**2 - 244*zh**3 + 249*zh**4 + 116*zh**5 + 4*zh**6)) + \
            M**4*Q**10*xh**2*zh*(2*xh**6*(-3 - 35*zh - 8*zh**2 + 46*zh**3) + zh*(-2 - 11*zh + 36*zh**2 - 33*zh**3 + 13*zh**4) - \
            2*xh**5*(12 - 31*zh - 159*zh**2 + 82*zh**3 + 96*zh**4) + xh**4*zh*(128 - 363*zh - 288*zh**2 + 408*zh**3 + 120*zh**4) - \
            2*xh**3*zh*(-23 + 145*zh - 327*zh**2 + 81*zh**3 + 124*zh**4 + 12*zh**5) - \
            4*xh*(-4 + 26*zh - 70*zh**2 + 110*zh**3 - 93*zh**4 + 37*zh**5 - 2*zh**6 + zh**7) + xh**2*(14 - 60*zh + 82*zh**2 + 110*zh**3 - 393*zh**4 + 263*zh**5 + 16*zh**6 + 4*zh**7)) + \
            M**6*Q**8*xh**3*(2*xh**5*zh*(-19 - 73*zh + 38*zh**2 + 54*zh**3) - zh*(12 - 58*zh + 108*zh**2 - 84*zh**3 + 33*zh**4 + zh**5) - \
            2*xh**4*(4 + 19*zh - 145*zh**2 - 146*zh**3 + 196*zh**4 + 72*zh**5) + 2*xh**3*(-4 + 33*zh + 17*zh**2 - 316*zh**3 + 68*zh**4 + 188*zh**5 + 24*zh**6) - \
            2*xh**2*(4 - 41*zh + 122*zh**2 - 87*zh**3 - 173*zh**4 + 167*zh**5 + 32*zh**6 + 4*zh**7) + \
            xh*(-8 + 68*zh - 252*zh**2 + 534*zh**3 - 554*zh**4 + 265*zh**5 - 17*zh**6 + 12*zh**7))))/ \
            (M**3*(Q**2*(-1 + xh) + M**2*xh)**3*(M**2*Q*xh + Q**3*(1 + xh - zh))**4*(-1 + zh)*zh**3)
         
         FUUcos2phi3pjq = 4*(-16*xh*(M**2*xh + Q**2*(-1 + xh)*zh)*(M**4*xh**2*(-3 + 7*xh*(-1 + zh)) + 2*M**2*Q**2*xh*(-1 + 5*xh**2*(-1 + zh) + zh - xh*zh**2) + \
            Q**4*(3*xh**3*(-1 + zh) + (-1 + zh)**2 + xh**2*(-1 + 8*zh - 6*zh**2) + xh*(3 - zh - 5*zh**2 + 3*zh**3))))/ \
            (3.*M**3*Q**2*(M**2*xh + Q**2*(1 + xh - zh))**4*(-1 + zh))
         
         FUUcos2phig = (FUUcos2phi3pjg)*pdf.gluon_PDF
         FUUcos2phiq = (FUUcos2phi3pjq)*pdf.qqb_PDF

         return e_c**2*alpha_QED*pdf.als**2*(FUUcos2phig + FUUcos2phiq)
########################################################  
  
