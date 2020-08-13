#-------------------
# don_PBPK.mc.in.r
#-------------------
Integrate (Lsodes, 1e-4, 1e-6, 1);

MonteCarlo("sim.out", # name of output and restart file
     1000,           # number of runs
     102020202);    # random seed (default )

    Distrib(Vmax_vitro, LogUniform, 0.0001, 0.01);
    Distrib(Km_vitro, LogUniform, 10000, 100000);


  Simulation {#1
  PO_dose = 1;
  PrintStep (Q_elim_kid, 1,24,1);
  }
End. 

