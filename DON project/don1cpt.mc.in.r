#-------------------
# don1cpt.mc.in.r
#-------------------
Integrate (Lsodes, 1e-4, 1e-6, 1);

MonteCarlo("sim.out", # name of output and restart file
     50000,           # number of runs
     10101010);    # random seed (default )


  
  
  Simulation {
  PO_dose = 1;
  Distrib(km_d15g, LogUniform, 0.01,10);
  Distrib(km_d3g, LogUniform, 0.01,10);
  Distrib(ku_d15g, LogUniform, 0.01,10);
  Distrib(ku_d3g, LogUniform, 0.01, 10);
  Distrib(ke_don, LogUniform, 0.01, 10);
  Print (Pct_d15g_ex, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,12,13,14,15,16,17,18,19,20,21,22,23,24);
  Data (Pct_d15g_ex, 0.9448, 8.0286,6.9637,5.5285,8.4916,3.7691,4.2321,3.6765,1.2689,0.9911,1.0837,0.6207,0.4355,0.8985,1.6393,1.6393,3.2135,0.7596,1.6856,0.7133,0.7133,0.2503,0.4587,0.2503);
  }
} End. 
