#####################################################################
# GHG.R
#
# R function that calculates dissolved pCO2, pCH4, and pN2O from raw GC data from headspace
# equilibration.  This function uses only freshwater constants, updating the Henry’s law 
# constants from to Sander 2015.

GHG <-  function(...){
  arguments <- list(...)
  
  # test arguments and initialize variables
  if (is.data.frame(arguments[[1]])) {
    input.table=arguments[[1]]
    if (dim(input.table)[2]!=13){
      stop("You should input a data frame with 13 columns", call.=FALSE)
    }else{
      Site = as.character(input.table$site)
      Timestamp = input.table$datetime.EST
      mCH4_headspace = input.table$HS.mCH4.before #the mCH4 (ppmv) of the headspace "before" equilibration
      mCO2_headspace = input.table$HS.mCO2.before #the mCO2 (ppmv) of the headspace "before" equilibration
      mN2O_headspace = input.table$HS.mN2O.before #the mN2O (ppmv) of the headspace "before" equilibration
      mCH4_eq = input.table$HS.mCH4.after #the measured mCH4 (ppmv) of the headspace "after" equilibration
      mCO2_eq = input.table$HS.mCO2.after #the measured mCO2 (ppmv) of the headspace "after" equilibration
      mN2O_eq = input.table$HS.mN2O.after #the measured mN2O (ppmv) of the headspace "after" equilibration
      temp_insitu = input.table$Temp.insitu #in situ water temperature in degrees celsius
      temp_eq = input.table$Temp.equil #the water temperature after equilibration in degree celsius
      #alk = input.table$Alkalinity.measured #Total alkalinity (micro eq/L) of the water sample
      vol_gas = input.table$Volume.gas #Volume of gas in the headspace vessel (mL)
      vol_water = input.table$Volume.water #Volume of water in the headspace vessel (mL)   
      Bar.pressure = input.table$Bar.pressure #Barometric pressure at field conditions in kPa. 101.325 kPa = 1 atm   
    } 
  } else if (length(arguments)==13) {
    Site = as.character(arguments[[1]])
    Timestamp = arguments[[2]]
    mCH4_headspace = arguments[[3]] #the mCH4 (ppmv) of the headspace "before" equilibration
    mCO2_headspace = arguments[[4]] #the mCO2 (ppmv) of the headspace "before" equilibration
    mN2O_headspace = arguments[[5]] #the mN2O (ppmv) of the headspace "before" equilibration
    mCH4_eq = arguments[[6]] #the measured mCH4 (ppmv) of the headspace "after" equilibration
    mCO2_eq = arguments[[7]] #the measured mCO2 (ppmv) of the headspace "after" equilibration
    mN2O_eq = arguments[[8]] #the measured mN2O (ppmv) of the headspace "after" equilibration
    temp_insitu = arguments[[9]] #in situ water temperature in degrees celsius
    temp_eq = arguments[[10]] #the water temperature after equilibration in degree celsius
   # alk = arguments[[11]] #Total alkalinity (micro eq/L) of the water sample
    vol_gas = arguments[[11]] #Volume of gas in the headspace vessel (mL)
    vol_water = arguments[[12]] #Volume of water in the headspace vessel (mL)   
    Bar.pressure = arguments[[13]] #Barometric pressure at field conditions in kPa. 101.325 kPa = 1 atm   
  } else {
    stop("You should input either a data frame or a vector of 13 values", call.=FALSE)
  }
  
  
  
  #initialization of variables
  pGHG <- data.frame(matrix(NA,length(mCO2_headspace),11))
  names(pGHG) <- c("Site",
                        "Timestamp",
                        "mCO2 simple headspace (ppmv)", 
                        "pCO2 simple headspace (micro-atm)", 
                        "CO2 concentration simple headspace (micro-mol/L)", 
                        "mCH4 simple headspace (ppmv)", 
                        "pCH4 simple headspace (micro-atm)", 
                        "CH4 concentration simple headspace (micro-mol/L)",
                        "mN2O simple headspace (ppmv)", 
                        "pN2O simple headspace (micro-atm)", 
                        "N2O concentration simple headspace (micro-mol/L)")
  
  R <- 0.082057338 #L atm K-1 mol-1
  
  #the function uniroot cannot handle vectors, so we need a loop
  for (i in 1:length(mCO2_headspace)){ 
    
    #AT = alk[i]*(1e-6) #conversion to mol/L
    
    #Constants of the carbonate equilibrium
    # Kw = the dissociation constant of H2O into H+ and OH-
    # Kh.co2 = the solubility of CO2 in water - equilibration conditions
    # Kh2.co2 = the solubility of CO2 in water - in situ field conditions
    # K1 = the equilibrium constant between CO2 and HCO3-
    # K2 = the equilibrium constant between HCO3- and CO3 2-
    # Kh.ch4 = the solubility of CH4 in water - equilibration conditions
    # Kh2.ch4 = the solubility of CH4 in water - in situ field conditions
    # Kh.n2o = the solubility of N2O in water - equilibration conditions
    # Kh2.n2o = the solubility of N2O in water - in situ field conditions
    
    # Solubility coefficients from Sander 2015 (1.4 x 10^-5 mol m^-3 Pa ^-1 and 1900 K for CH4, 
    # 3.3 x 10^-4 mol m^-3 Pa ^-1 and 2400 K for CO2, 
    # and 2.4 x 10^-4 mol m^-3 Pa ^-1 and 2700 K for N2O)
    # Dissociation of water from Dickson and Riley (1979)
    # K1 and K2 from Millero, F. (1979).  
    K1=10^-(-126.34048+6320.813/(temp_eq[i]+273.15)+19.568224*log(temp_eq[i]+273.15))
    K2=10^-(-90.18333+5143.692/(temp_eq[i]+273.15)+14.613358*log(temp_eq[i]+273.15))
    Kw = exp(148.9652-13847.26/(temp_eq[i]+273.15)-23.6521*log(273.15+temp_eq[i]))
    
    Kh.co2 = 0.00033*exp(2400*(1/(273.15+temp_eq[i])-1/298.15))*101325/1000 # mol/L/atm equilibration conditions
    Kh2.co2 = 0.00033*exp(2400*(1/(273.15+temp_insitu[i])-1/298.15))*101325/1000 # mol/L/atm original conditions
    Kh.ch4 = 0.000014*exp(1900*(1/(273.15+temp_eq[i])-1/298.15))*101325/1000# mol/L/atm equilibration conditions
    Kh2.ch4 = 0.000014*exp(1900*(1/(273.15+temp_insitu[i])-1/298.15))*101325/1000 # mol/L/atm original conditions
    Kh.n2o = 0.00024*exp(2700*(1/(273.15+temp_eq[i])-1/298.15))*101325/1000 # mol/L/atm equilibration conditions
    Kh2.n2o = 0.00024*exp(2700*(1/(273.15+temp_insitu[i])-1/298.15))*101325/1000 # mol/L/atm original conditions
    
    HS.ratio <- vol_gas[i]/vol_water[i] #Headspace ratio (=vol of gas/vol of water)
    
    #The following calculations assume 1 atm, this is corrected later for measured pressure in the field.
    
    pGHG[i,1] <- as.character(Site[i])
    pGHG[i,2] <- Timestamp[i]

    
    
    #Calculation not accounting for alkalinity effects and associated error
    #concentration and total mass in the water sample assuming ideal gas from the pCO2 measured at the headspace
    CO2_solution <- mCO2_eq[i]/1000000*Kh.co2 #mol/L
    CO2_solution_mass <- CO2_solution * vol_water[i]/1000 #mol
    
    #mass of CO2 in the measured headspace
    final_C_headspace_mass <- mCO2_eq[i]/1000000*(vol_gas[i]/1000) / (R * (temp_eq[i]+273.15)) #mol
    
    mols_headspace <- mCO2_headspace[i]/1000000*(vol_gas[i]/1000)/(R * (temp_eq[i]+273.15)) #mol PV / RT = n
    
    #implication: mass, concentration, and partial pressure of CO2 in the original sample (aount in sample and headspace after equilibration minus original mass in the headspace)
    Sample_CO2_mass <- CO2_solution_mass + final_C_headspace_mass - mols_headspace #mol
    Sample_CO2_conc <- Sample_CO2_mass/(vol_water[i]/1000) #mol/L
    pGHG[i,3] <- Sample_CO2_conc/Kh2.co2*1000000 #ppmv
    pGHG[i,4] <- Sample_CO2_conc/Kh2.co2*1000000*Bar.pressure[i]/101.325 # micro-atm
    pGHG[i,5] <- Sample_CO2_conc*1000000 # micro-mol/L
    
    #concentration and total mass in the water sample assuming ideal gas from the pCH4 measured at the headspace
    CH4_solution <- mCH4_eq[i]/1000000*Kh.ch4 #mol/L
    CH4_solution_mass <- CH4_solution * vol_water[i]/1000 #mol
    
    #mass of CH4 in the measured headspace
    final_C_headspace_mass.ch4 <- mCH4_eq[i]/1000000*(vol_gas[i]/1000) / (R * (temp_eq[i]+273.15)) #mol
    
    mols_headspace.ch4 <- mCH4_headspace[i]/1000000*(vol_gas[i]/1000)/(R * (temp_eq[i]+273.15)) #mol PV / RT = n
    
    #implication: mass, concentration, and partial pressure of CH4 in the original sample (aount in sample and headspace after equilibration minus original mass in the headspace)
    Sample_CH4_mass <- CH4_solution_mass + final_C_headspace_mass.ch4 - mols_headspace.ch4 #mol
    Sample_CH4_conc <- Sample_CH4_mass/(vol_water[i]/1000) #mol/L
    pGHG[i,6] <- Sample_CH4_conc/Kh2.ch4*1000000 #ppmv
    pGHG[i,7] <- Sample_CH4_conc/Kh2.ch4*1000000*Bar.pressure[i]/101.325 # micro-atm
    pGHG[i,8] <- Sample_CH4_conc*1000000 # micro-mol/L
    
    #concentration and total mass in the water sample assuming ideal gas from the pN2O measured at the headspace
    N2O_solution <- mN2O_eq[i]/1000000*Kh.n2o #mol/L
    N2O_solution_mass <- N2O_solution * vol_water[i]/1000 #mol
    
    #mass of N2O in the measured headspace
    final_N_headspace_mass <- mN2O_eq[i]/1000000*(vol_gas[i]/1000) / (R * (temp_eq[i]+273.15)) #mol
    
    mols_headspace.n2o <- mN2O_headspace[i]/1000000*(vol_gas[i]/1000)/(R * (temp_eq[i]+273.15)) #mol PV / RT = n
    
    #implication: mass, concentration, and partial pressure of N2O in the original sample (aount in sample and headspace after equilibration minus original mass in the headspace)
    Sample_N2O_mass <- N2O_solution_mass + final_N_headspace_mass - mols_headspace.n2o #mol
    Sample_N2O_conc <- Sample_N2O_mass/(vol_water[i]/1000) #mol/L
    pGHG[i,9] <- Sample_N2O_conc/Kh2.n2o*1000000 #ppmv
    pGHG[i,10] <- Sample_N2O_conc/Kh2.n2o*1000000*Bar.pressure[i]/101.325 # micro-atm
    pGHG[i,11] <- Sample_N2O_conc*1000000 # micro-mol/L
  }
  
  
  return(pGHG) #Output data frame
  
}
