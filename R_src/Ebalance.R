############################################################
# Ebalance function
# Author: Marc Peaucelle, 2020/06
# Description:
# This function computes the bud energy balance as described in Hamer et al. 1985, and Landsberg et al. 1974
# and is adapted from the tealeaves R package (Muir 2019)
#
# Important Note: This model is simplified to represent idealized conditions.
# To have robust simulations of bud temperature and to use this approach to model phenology, 
# a proper calibration of model parameters, including bud traits (e.g. solar absorptivity, size and shape), but also the surrounding environment (e.g. ground albedo and soil temperature)
# Running the model with micrometeorological data should also produce better results for bud temperature
#
# References: 
#     - Landsberg, J. J., Butler, D. R. & Thorpe, M. R. Apple bud and blossom temperatures. Journal of Horticultural Science 49, 227â€“239; 10.1080/00221589.1974.11514574 (1974).
#     - HAMER, P. The heat balance of apple buds and blossoms. Part I. Heat transfer in the outdoor environment. Agricultural and Forest Meteorology 35, 339â€“352; 10.1016/0168-1923(85)90094-2 (1985).
#     - Muir, C. D. tealeaves: an R package for modelling leaf temperature using energy budgets. AoB PLANTS 11, plz054; 10.1093/aobpla/plz054 (2019).
#     - Monteith, J. L. & Unsworth, M. H. Principles of environmental physics. Plants, animals, and the atmosphere. 4th ed. (Elsevier/Academic Press, Amsterdam, Boston, 2013).
#     - Jones, H. G. Plants and microclimate. A quantitative approach to environmental plant physiology (Cambridge university press, Cambridge, 2013). --> provides absroptivity values for leaves and needles
############################################################
# Inputs:
# Tbud: Bud temperature (K)
# T_air: air temperature (K) 
# P: atmospheric pressure (kPa)
# RH: relative humidity (0-1)
# Sw: incoming shortwave radiation (W m-2)
# Lw: incoming longwave radiation (W m-2)
# wind: wind speed (m s-1)
# bud_d: bud diameter
# abs_sw; solar ab
# Ebal_only: Flag to return the energy balance only or all the components
############################################################
# outputs: 
# Ebal: Energy balance (W m-2)
# Rabs; Absorbed radiation (W m-2)
# Sr: Longwave radiation emitted by the bud(W m-2)
# H: Sensible heat flux (W m-2)

############################################################

Ebalance<-function(Tbud, T_air, T_ground, P, RH, qair, Sw, Lw, wind, rain, bud_d=0.005, abs_sw=0.5, r=0.2, 
                   evap = TRUE, md = NULL, max_dw = 0.15, tstep = 1800, kE=1, Ebal_only=FALSE){
  
  # constants
  c_p =  1.01 # Heat capacity of air J / (g * K)
  c_w =  1.82 # Heat capacity of water J / (g * K)
  c_b = 2820.12 # Heat capacity of bud in J kg-1 C-1
  rho_b = 782.93 # density of bud in kg m-3
  D_h0 = 1.9e-5 # Diffusion coefficient for heat in air at 0Â°C m ^ 2 /² s
  D_m0 = 13.3e-6 # Diffusion coefficient for momentum in air at 0Â°C m ^ 2 / s
  #  D_w0 = 21.2e-6 # Diffusion coefficient for water vapour in air at 0Â°C m ^ 2 / s
  eT = 1.75 # Exponent for temperature dependence of diffusion 
  G = 9.8 # Gravitational acceleration m / s ^ 2
  R = 8.3144598 # Ideal gas constant J / (mol * K)
  R_air = 287.058 # Specific gas constant for dry air J / (kg * K)
  s = 5.67e-08 # Stefan-Boltzmann constant W / (m ^ 2 * K ^ 4)
  #abs_sw = 0.5 # bud absorptivity for shortwave radiation defined by the user
  abs_lw = 0.97 # bud absorptivity for longwave radiation
  
  # Initialize latent heat and storage
  E = 0
  
  if (rain > 0 & evap){
    # Sl = surface area of the bud in m2
    Sl = 4 * pi * (bud_d/2)**2
    
    # max_md = maximum dew mass on bud (in kg) = max water density (in mm = kg/m²) * area
    max_md = max_dw * Sl
    md = min(max_md,(md+rain*Sl))
  }
  
  # if buds are wet, change the albedo 
  if(md>0 & evap){
    abs_sw = 0.9
  }
  ## 1. Estimates energy balance for given environmental conditions
  ## 1.1 Total absorbed radiation in W m^-2 
  # Rabs = Rsw + Rlw 
  # with Rlw =  (La+Lg)/2 for a spherical object as defined in Hammer 1985
  # La and Lg are the atmosphere and ground emitted radiation respectively
  # Here Lg is estimated from air temperature, while soil temperature is generally higher than air
  
  # r = soil albedo. Default defined to 0.2 as in Muir 2019 
  # average radiation sky and ground -> /2
  R_abs = abs_sw * (1 + r) * Sw /2 + ((abs_lw * Lw ) + (abs_lw * s *  T_ground ^ 4)) /2 # Should be replaced by real soil T instead of T_air
  
  
  ## 1.2 re-emitted longwave radiation in W m^-2 : Foster and Smith 1986
  # Tbud in K, I kept the notation from Muir 2019 for an easier comparison with leaf energy balance. 
  # S_r corresponds to LW_bud in Peaucelle, Penuelas & Verbeeck 2020.
  S_r =  s * abs_lw * Tbud ^ 4
  
  ## 1.3 Sensible heat flux in W m^-2
  ## 1.3.1  Density of air in g m-3
  ## 1.3.1.1  Density of dry air in g m-3# P in kPa
  # T in K
  # multiply by 1e6 to get g m-3
  P_a = 1e6 * P / (R_air * (T_air + Tbud) / 2)
  
  ## 1.3.1.2  Density of air rho_air
  # R_water = 461.55 J kg-1 K-1
  # R_air  = 287.058 J kg-1 K-1
  # R_water/R_air = 1.608
  
  rho_air = P_a * (1 + qair)/(1 + 1.608 * qair)
  
  
  #### Now estimate boundary layer conductance to heat
  ## 1.3.2 Diffusion coefficient to heat
  # T in K
  # P in kPa
  # if(md>0 & evap){
  #   D_h0 = 0.144 * 1e-6 
  #   D_h0 = 10 * 1e-6
  # }
  D_h = D_h0 * (((T_air+Tbud)/2) / 273.15) ^ eT * (101.3246 / P)     
  
  ## 1.3.3 Calculate Grashof number
  # P in kPa
  # T in K
  Tv_bud = Virtual_Temp(Tbud,P)
  Tv_air = Virtual_Temp(T_air,P,RH)
  # if(md>0 & evap){
  #   D_m0 = 1.01 * 1e-6 
  # }
  D_m = D_m0 * (((T_air+Tbud)/2) / 273.15) ^ eT * (101.3246 / P)
  
  # D_m from Schymanski et al. 2013 and Gerlein-safdi
  # D_m = (T_air+Tbud)/2 * 9 * 1e-8 - 1.13 * 1e-5
  
  Gr = (1 / T_air) * G * bud_d ^ 3 * abs(Tv_bud - Tv_air) / D_m ^ 2
  
  ## 1.3.4  Calculate Reynolds number
  Re =  wind * bud_d/D_m
  
  ## 1.3.5  Calculate Nusselt number
  # D_w = D_w0 * (((T_air+Tbud)/2) / 273.15) ^ eT * (101.3246 / P)
  
  # Nu for a smooth sphere is described in Monteith 2013
  # see nu_constant function at the end of the script. 
  # a,b,c,d are constants that depend on whether flow is laminar or turbulent and the direction of flow in the case of free convection (see below).
  # Nu_free, Gr^0:25<220 Nu = 2 + 0.54*Gr^0.25
  # Nu_forced, Re<300 Nu = 2 + 0.54*Re^0.5
  # Nu_forced, Re>50 Nu = 0.34Re^0.6
  
  # In general, when the Archimedes number (also called the Richardson number) Ar=Gr/Re2â‰ª0.1`, free convection dominates; when Ar=Gr/Re2â‰«10`, forced convection dominates (Nobel 2009). 
  # Chose here the shape of the object for Nusselt constants
  #shape = "sphere"
  shape = "cylinder"
  nu_cst = nu_constant(Re,1,T_air,Tbud,Gr,shape)
  Nu_forced = nu_cst$a1 + nu_cst$a * Re ^ nu_cst$b 
  nu_cst = nu_constant(Re,2,T_air,Tbud,Gr,shape)
  Nu_free =  nu_cst$a1 + nu_cst$a * Gr ^ nu_cst$b
  
  # Nu mixed equation
  
  Nu = (Nu_forced ^ 3.5 + Nu_free ^ 3.5) ^ (1 / 3.5)
  
  # Nu from gerlrein-safdi
  Npr = 0.71 # Prandtl number
  # Nu = 0.664 * (Re**(1/2)) * (Npr**(1/3))
  
  ## 1.3.6 Boundary layer conductance to heat
  g_h = D_h*Nu/bud_d
  
  ## 1.3.7 Final sensible heat flux
  # heat capacity of moist air cs
  c_s = c_p + c_w * qair  
  
  H = rho_air * c_s * g_h * (Tbud - T_air) 
  
  gw = 0
  ## 1.4 Evaporative cooling
  if (evap) {
    
    # Sl = surface area of the bud in m2
    Sl = 4 * pi * (bud_d/2)**2
    
    # max_md = maximum dew mass on bud (in kg) = max water density (in mm = kg/m²) * area
    max_md = max_dw * Sl
    
    
    # alapha_a =  thermal diffusivity of air (m² s-1)
    # Dva = diffusivity of water vapor in air (m² s-1)
    Tb = (T_air + Tbud)/2
    alpha_a = Tb * 1.32e-7 - 1.73e-5
    Dva = Tb * 1.49e-7 - 1.96e-5
    # Nle = Lewis number
    Nle = alpha_a / Dva
    
    # ka = thermal conductivity of air in the boundary layer in W K-1 m-1
    # L = characteristic length of leave in m 
    # Nnu = Nusselt number (different function than in Muir)
    
    ka = Tb * 6.84e-5 + 5.62e-3
    hc = ka * Nu/bud_d
    
    gw = hc / (rho_air * c_s * Nle^(2/3)) #XXX test *2
    
    # Dmd = dmd/dt in kg s-1
    # md is the dew mass in kg 
    # ec = vapor pressure in Pa (Tair)
    # esat = saturation vapor pressure in Pa (Tbud)
    esat <-  6.112 * exp((17.67 * (Tbud-273.15))/((Tbud -273.15)+ 243.5))
    ec <- qair * P * 10 / (0.378 * qair + 0.622)
    
    Dmd = 0.622 * Sl * rho_air * gw * ((ec-esat)/P) 
    
    # gamma_v = gamma_c = the latent heat of vaporization/condensation in J kg-1
    gamma_v = 1e6 * (2.5 - 2.36 * 1e-3 * (Tbud-273.15))
    
    # formule Muir in J mol-1
    # 56847.68250 - 43.12514 * T[K] (Nobel, 2019)
    # 1 mole  = 18.02 grams H2O --> *1000/18.02 in J kg-1
    #gamma_v = 3154699 - 2393.182 * Tbud
    
    
    # integrate over time
    if (Dmd >= 0){
      # if dew is accumulating, no latent heat
      E = 0
    } else {
      if((md + Dmd * tstep) > 0){
        E =  gamma_v * Dmd * tstep
      } else {
        E = gamma_v * md
      }
    }
    
    # adjust md 
    md = max(c(0,min(c(md + Dmd,max_md * tstep))))
  }
  
  ## 1.5 Energy balance
  
  Ebal = R_abs - ( S_r + H + kE*E ) #XXX
  
  if(Ebal_only){  
    out<-Ebal
  } else {
    out<-list(R_abs=R_abs,
              S_r=S_r,
              H=H,
              E=E,
              Tbud=Tbud,
              md = md,
              gw=gw,
              g_h=g_h,
              Nu=Nu,
              Ebal=Ebal)
  }
  return(out)
}

nu_constant<-function(Re,type,T_air,Tbud,Gr,shape){
  # Nu for a smooth sphere is described in Monteith 2013, p.
  # see nu_constant function at the end of the script. 
  # a,b,c,d are constants that depend on whether flow is laminar or turbulent and the direction of flow in the case of free convection (see below).
  # Nu_free, Gr^0:25<220 Nu = 2 + 0.54*Gr^0.25
  # Nu_forced, Re<300 Nu = 2 + 0.54*Re^0.5
  # Nu_forced, Re>50 Nu = 0.34Re^0.6
  if (shape == "sphere"){
    if (type==1) { # type = forced
      # laminar flow
      if (Re <= 300){
        a1 = 2
        a = 0.54
        b = 0.5
      } else {
        # turbulent flow
        a1 = 0
        a = 0.34
        b = 0.6
      }
    } else{ # type = free
      a1 = 2
      a = 0.54
      b = 0.25
    }
  } else if (shape == "cylinder") {
    # cylinder
    
    if (type==1) { # type = forced
      # Wide range
      # laminar flow
      # if (Re <= 1000){
      #   a1 = 0.32
      #   a = 0.51
      #   b = 0.55
      # } else {
      #   # turbulent flow
      #   a1 = 0
      #   a = 0.24
      #   b = 0.6
      # }
      
      if (Re <= 4){
        a1 = 0
        a = 0.89
        b = 0.33
      } else if (4 < Re & Re <= 40) {
        a1 = 0
        a = 0.82
        b = 0.39
      } else if (40 < Re & Re <= 4000) {
        a1 = 0
        a = 0.62
        b = 0.47
      } else if (4000 < Re & Re <= 400000) {
        a1 = 0
        a = 0.17
        b = 0.62
      } else {
        a1 = 0
        a = 0.024
        b = 0.81
      }
      
    } else{ # type = free for horizontal cylinder
      if(Gr < 10^9){
        a1 = 0
        a = 0.48
        b = 0.25
      } else {
        a1 = 0
        a = 0.09
        b = 0.33
      }
    }
  } else {
    print("Indicate shape = sphere or shape = cylinder")
  }
  ############### Values for leaves available in the tealeaves package, Muir 2019 ##################
  # Note that in Muir 2019, an additional arguement "surface" is needed, corresponding to the upper or lower surface of the leaf
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
  
  return (list(a1=a1,
               a=a,
               b=b))
}

Virtual_Temp<-function (temp,P,RH=1){
  # P = 101.3246 # kPa
  # temp in K
  # RH (0-1), assumed to be 100% in the bud
  epsilon = 0.622 # ratio of water to air molar masses
  p_s = 10 ^ (-7.90298 * (373.16 / temp - 1) +
                5.02808 * log10(373.16 / temp) -
                1.3816e-7 * (10 ^ (11.344 * (1 - temp / 373.16) - 1)) +
                8.1328e-3 * (10 ^ (-3.49149 * (373.16 / temp - 1)) - 1) +
                log10(P*10)) # P*10 --> hPa
  p_s = p_s*RH
  # virtual temperature of leaf Eq. 2.35 in Monteith & Unsworth (2013)
  Tv = temp / (1 - (1 - epsilon) * (p_s*0.1 / P)) #p_s*0.1--> kPa
  return(Tv)
}


