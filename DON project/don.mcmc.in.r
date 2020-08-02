#-------------------
# methyleugenol.mcmc.in
#-------------------
Integrate (Lsodes, 1e-4, 1e-6, 1);

MCMC("sim.out","", # name of output and restart file
     "",           # name of data file
     10000,0,      # iterations, print predictions flag,
     1,10000,     # printing frequency, iters to print
     10101010);    # random seed (default )

Level { # top
  Distrib(Bodywt, TruncNormal, 83.2711,20.8715, 33.4669, 182.9474);
  Distrib(Pct_M_fat, TruncNormal, 0.4194,0.1228, 0.0939, 0.6909); # Percentage of fat by total mass (mean, std, min, max)
  Distrib(Pct_M_liv, TruncNormal, 0.02, 0.0065, 0.0041, 0.0483);# liver
  Distrib(Pct_M_rp, TruncNormal, 0.0473, 0.0108, 0.0227, 0.1112);# rapidly perfused tissues
  Distrib(Pct_M_sp, TruncNormal, 0.4486, 0.1088, 0.2305, 0.7624);# slowly perfused tissues
  Distrib(Pct_M_blood, TruncNormal, 0.0647, 0.0102, 0.0397, 0.1003);#blood
  Distrib(Flow_tot, TruncNormal, 349.3047,36.4727, 253.9084, 435.7289);# total cardiac output
  Distrib(Pct_Flow_fat, TruncNormal, 0.0678, 0.0179, 0.0418, 0.0960); # fraction of blood flow to fat 
  Distrib(Pct_Flow_liv, TruncNormal, 0.2153, 0.0118, 0.1886, 0.2440); # fraction of blood flow to liver
  Distrib(Pct_Flow_rp, TruncNormal, 0.4722, 0.0116, 0.4397, 0.5071); # fraction of blood flow to rapidly perfused tissues 
  Distrib(Pct_Flow_sp, TruncNormal, 0.2447, 0.0265, 0.1986, 0.2991); # fraction of blood flow to slowly perfused tissues
  Distrib(Vmax_1HME_1A2, TruncNormal_cv, 0.860, 0.67, 0.1169, 4.5606); 
  Distrib(Vmax_1HME_2C9, TruncNormal_cv, 1.210, 0.54, 0.2321, 4.6118); 
  Distrib(Vmax_1HME_2C19, TruncNormal_cv, 0.2092, 1.05, 0.0115, 1.8638); 
  Distrib(Vmax_1HME_2D6, TruncNormal_cv, 0.0555, 0.61, 0.0090, 0.2411); 
  Distrib(Vmax_1HME_3A4, TruncNormal_cv, 2.0892, 1.14, 0.0892, 20.2598);
  Distrib(Vmax_MEO_2B6, TruncNormal_cv, 0.6721, 0.39, 0.2068, 1.9192);
  Distrib(Vmax_MEO_2E1, TruncNormal_cv, 0.3243, 0.61, 0.0530, 1.4213);
  Distrib(Vmax_1EU_1A2, TruncNormal_cv, 0.1543, 0.67, 0.0217, 0.8276);
  Distrib(Vmax_1HMEG_1A9, TruncNormal_cv, 0.0983, 0.49, 0.0211, 0.3948);
  Distrib(Vmax_1HMEG_2B7, TruncNormal_cv, 0.2923, 0.45, 0.0398, 1.3526);
  Distrib(Vmax_1OME, TruncNormal_cv, 44.4406, 0.27, 19.5618, 94.6727);
  Distrib(Vmax_1HMES, TruncNormal_cv, 0.9982, 0.51, 0.2188, 3.6805);
  Distrib(Km_1OME, TruncNormal_cv, 1567, 0.46, 384, 5251);
  Distrib(Km_1HMES, TruncNormal_cv, 282, 0.76, 27, 1662);
  
  
  
   Likelihood (C_ven, LogNormal, Prediction(C_ven), 1.1);
  
  Simulation {
  PO_dose = 0.002632;
  Print (C_ven, 0.25, 0.5, 1, 2);
  Data (C_ven, 0.000302356, 0.000240638,0.000207592,0.000141524);
  }
} End. 
