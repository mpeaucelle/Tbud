############################################################
# Ebalance function
# Author: Marc Peaucelle, 2020/06
# Description:
# This function computes the bud energy balance as described in Hamer et al. 1985, and Landsberg et al. 1974
# and is adapted from the tealeaves R package (Muir 2019)
# References: 
#     - Landsberg, J. J., Butler, D. R. & Thorpe, M. R. Apple bud and blossom temperatures. Journal of Horticultural Science 49, 227–239; 10.1080/00221589.1974.11514574 (1974).
#     - HAMER, P. The heat balance of apple buds and blossoms. Part I. Heat transfer in the outdoor environment. Agricultural and Forest Meteorology 35, 339–352; 10.1016/0168-1923(85)90094-2 (1985).
#     - Muir, C. D. tealeaves: an R package for modelling leaf temperature using energy budgets. AoB PLANTS 11, plz054; 10.1093/aobpla/plz054 (2019).
#     - Monteith, J. L. & Unsworth, M. H. Principles of environmental physics. Plants, animals, and the atmosphere. 4th ed. (Elsevier/Academic Press, Amsterdam, Boston, 2013).
#     - Jones, H. G. Plants and microclimate. A quantitative approach to environmental plant physiology (Cambridge university press, Cambridge, 2013). --> provides absroptivity values for leaves and needles
############################################################
# Inputs:
# Tbud: leaf temperature (K)
# T_air: air temperature (K) 
# P: atmospheric pressure (kPa)
# RH: relative humidity (0-1)
# Sw: incoming shortwave radiation (W m-2)
# Lw: incoming longwave radiation (W m-2)
# wind: wind speed (m s-1)
# bud_d: bud diameter
# Ebal_only: Flag to return the energy balance only or all the components
############################################################
# outputs: 
# Ebal: Energy balance (W m-2)
# Rabs; Absorbed radiation (W m-2)
# Sr: Longwave radiation emitted by the bud(W m-2)
# H: Sensible heat flux (W m-2)

############################################################

Ebalance<-function(Tbud,T_air,P,RH, Sw,Lw,wind,bud_d,abs_sw=0.5,Ebal_only=FALSE){
  
  # constants
  c_p =  1.01 # Heat capacity of air J / (g * K)
  D_h0 = 1.9e-5 # Diffusion coefficient for heat in air at 0°C m ^ 2 / s
  D_m0 = 13.3e-6 # Diffusion coefficient for momentum in air at 0°C m ^ 2 / s
  D_w0 = 21.2e-6 # Diffusion coefficient for water vapour in air at 0°C m ^ 2 / s
  eT = 1.75 # Exponent for temperature dependence of diffusion 
  G = 9.8 # Gravitational acceleration m / s ^ 2
  R = 8.3144598 # Ideal gas constant J / (mol * K)
  R_air = 287.058 # Specific gas constant for dry air J / (kg * K)
  s = 5.67e-08 # Stefan-Boltzmann constant W / (m ^ 2 * K ^ 4)
  #abs_sw = 0.5 # leaf absorptivity for shortwave radiation defined by the user
  abs_lw = 0.97 # leaf absorptivity for longwave radiation
  
  ## 1. Estimates energy balance for given environmental conditions
  ## 1.1 Total absorbed radiation in W m^-2 
  # Rabs = Rsw + Rlw 
  # with Rlw =  (La+Lg)/2 for a spherical object as defined in Hammer 1985
  # La and Lg are the atmosphere and ground emitted radiation respectively
  # Here Lg is estimated from air temperature, while soil temperature is generally higher than air

  # r = soil albedo defined to 0.2 as in Muir 2019 
  r = 0.2
  R_abs = abs_sw * (1 + r) * Sw + ((abs_lw * Lw ) + (abs_lw * s *  T_air ^ 4))/2 # Can be replaced by real soil T instead of T_air

  
  ## 1.2 re-emitted longwave radiation in W m^-2 : Foster and Smith 1986
  # Tbud in K, I kept the notation from Muir 2019 for an easier comparison with leaf energy balance. 
  # S_r corresponds to LW_bud in Peaucelle, Penuelas & Verbeeck 2020.
  S_r =  s * abs_lw * Tbud ^ 4
  
  ## 1.3 Sensible heat flux in W m^-2
  ## 1.3.1  Density of dry air in g m-3
  # P in kPa
  # T in K
  # multiply by 1e6 to get g m-3
  P_a = 1e6 * P / (R_air * (T_air + Tbud) / 2)
  
  #### Now estimate boundary layer conductance to heat
  ## 1.3.2 Diffusion coefficient to heat
  # T in K
  # P in kPa
  D_h = D_h0 * (((T_air+Tbud)/2) / 273.15) ^ eT * (101.3246 / P)     
  
  ## 1.3.3 Calculate Grashof number
  # P in kPa
  # T in K
  Tv_bud = Virtual_Temp(Tbud,P)
  Tv_air = Virtual_Temp(T_air,P,RH)
  D_m = D_m0 * (((T_air+Tbud)/2) / 273.15) ^ eT * (101.3246 / P)
  Gr = (1 / T_air) * G * bud_d ^ 3 * abs(Tv_bud - Tv_air) / D_m ^ 2
  
  ## 1.3.4  Calculate Reynolds number
  Re =  wind * bud_d/D_m
  
  ## 1.3.5  Calculate Nusselt number
  
  D_w = D_w0 * (((T_air+Tbud)/2) / 273.15) ^ eT * (101.3246 / P)
  
  # Nu for a smooth sphere is described in Monteith 2013
  # see nu_constant function at the end of the script. 
  # a,b,c,d are constants that depend on whether flow is laminar or turbulent and the direction of flow in the case of free convection (see below).
  # Nu_free, Gr^0:25<220 Nu = 2 + 0.54*Gr^0.25
  # Nu_forced, Re<300 Nu = 2 + 0.54*Re^0.5
  # Nu_forced, Re>50 Nu = 0.34Re^0.6
  
  # In general, when the Archimedes number (also called the Richardson number) Ar=Gr/Re2≪0.1`, free convection dominates; when Ar=Gr/Re2≫10`, forced convection dominates (Nobel 2009). 
  
  nu_cst = nu_constant(Re,1,T_air,Tbud,1)
  Nu_forced = nu_cst$a * Re ^ nu_cst$b 
  nu_cst = nu_constant(Re,2,T_air,Tbud,1)
  Nu_free =  nu_cst$a * Gr ^ nu_cst$b

  # N Checked equation here for sphere. 
  Nu = (Nu_forced ^ 3.5 + Nu_free ^ 3.5) ^ (1 / 3.5)

 
  ## 1.3.6 Boundary layer conductance to heat
  g_h = D_h*Nu/bud_d


  ## 1.3.7 Final sensible heat flux
  H = P_a * c_p * g_h * (Tbud - T_air) 

  ## 1.4 Energy balance
  Ebal = R_abs - ( S_r + H )
  
  if(Ebal_only){  
    out<-Ebal
  } else {
    out<-list(R_abs=R_abs,
              S_r=S_r,
              H=H,
              Ebal=Ebal)
  }
  return(out)
}

nu_constant<-function(Re,type,T_air,Tbud){
  # Nu for a smooth sphere is described in Monteith 2013, p.
  # see nu_constant function at the end of the script. 
  # a,b,c,d are constants that depend on whether flow is laminar or turbulent and the direction of flow in the case of free convection (see below).
  # Nu_free, Gr^0:25<220 Nu = 2 + 0.54*Gr^0.25
  # Nu_forced, Re<300 Nu = 2 + 0.54*Re^0.5
  # Nu_forced, Re>50 Nu = 0.34Re^0.6
  if (type==1) { # type = forced
    # laminar flow
    if (Re <= 300){
      a = 0.54
      b = 0.5
    } else {
      # turbulent flow
      a = 0.34
      b = 0.6
    }
  } else{ # type = free
      a = 0.54
      b = 0.25
  }
  ############### Values for leaves avilable in the tealeaves package, Muir 2019 ##################
  #  if (type==1) { # type = forced
  #   # laminar flow
  #   if (Re <= 4000){
  #     a = 0.6
  #     b = 0.5
  #   } else {
  #     # turbulent flow
  #     a = 0.032
  #     b = 0.8
  #   }
  # } else{ # type = free
  #   if (((surface==1)&(Tbud>T_air))|
  #       ((surface==2)&(Tbud<T_air))) {
  #     
  #     a = 0.5
  #     b = 0.25
  #   } else {
  #     a = 0.23
  #     b = 0.25
  #   }
  # }
  
  return (list(a=a,
               b=b))
}

Virtual_Temp<-function (temp,P,RH=1){
  # P = 101.3246 # kPa
  # temp in K
  epsilon = 0.622 # ratio of water to air molar masses
  p_s = 10 ^ (-7.90298 * (373.16 / temp - 1) + 
                5.02808 * log10(373.16 / temp) - 
                1.3816e-7 * (10 ^ (11.344 * (1 - temp / 373.16) - 1)) + 
                8.1328e-3 * (10 ^ (-3.49149 * (373.16 / temp - 1)) - 1) + 
                log10(P*10)) # P*10 --> hPa
  p_s = p_s*RH
  # virtual temperature for leaves Eq. 2.35 in Monteith & Unsworth (2013) applied here for buds
  Tv = temp / (1 - (1 - epsilon) * (p_s*0.1 / P)) #p_s*0.1--> kPa
  return(Tv)
}

sh_constant = function(type, unitless) {
    a = 0.33 #forced
    if (type==2) a= 0.25
    return(a)
}
