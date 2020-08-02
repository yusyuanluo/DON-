# donpbtk1cptumb1.model ------DON and only two metabolites_d15GA & D3GA with bile and feces compartments 

# Constants
# =========
# Molecular weight (in grams/mol)
#MW_don = 296.32;
#Dose   = 0.07/(1000*MW_don); # 1 microgram/kg ----(moles)

# Default parameter values
# ========================
# Units:
# Volumes: liter
# Time:    hour
# Quantity:   nmol

# Conversions between mg and nMol:
#Dose = Dose/(1000*MW_don); in number of moles
#IngDose = Dose * 1e+9; # In nmol

States  = { 
  Q_DON_cpt,    # Quantity od don in central compartment (nmol)
  Q_d3g_cpt,    # Quantity of metabolite d3g in urine (nmol)
  Q_d15g_cpt,   # Quantity of metabolite d15g in urine  (nmol)
  Q_DON_ex,     # Quantity of d3g in excrete (nmol)
  Q_d3g_ex,    # Quantity of d15g in excrete (nmol)
  Q_d15g_ex,   # Quantity of d3g in excrete (nmol)
  };  

Inputs = {PO_dose
}; 

Outputs = {
  C_don_cpt, C_d3g_cpt, C_d15g_cpt, Q_total, Pct_don_ex, Pct_d3g_ex, Pct_d15g_ex
};
# Oral input modeling
PO_dose = 1; ##microgram per kw body weight
Oral_dose;   ### converted to nmol
Period = 1e3; ##study period 
Ka = 0.09; ## absorption rate constant (1/h)
R_ing = PerExp (Oral_dose, Period, 0.0, Ka);


#IngDose    = 236.33; # ingested input (nmol)
#Fgutabs    = 0.56; #
#kgutabs    = 0.09; # Intestinal absortion rate (/h)
#kgblad_mag = 8.00;
#kgblad_time = 15; 
#kgblad =    PerDose(kgblad_mag, 48, kgblad_time, 1);

# Distribution volumes (L)
Vdist = 86.8;

# Body weight (kg)
BW = 70;

# Elimination rate constants (/h)
km_d15g    = 3.29;     #metabolic rate constant for D15GA (nmol/h)
km_d3g     = 0.73;     #metabolic rate constant for D3GA (nmol/h)
ku_d15g    = 2.31;     #Excretion rate constant for D15GA (nmol/h)
ku_d3g     = 1.28;     #Excretion rate constant for D3GA (nmol/h)
ke_don     = 1.28;     #Excretion rate constant for DON (nmol/h)
#kbile_1    = 1.00;     #Biliary excretion Rate constant DON-15G(nmol/h)
#kbile_2    = 1.00;     #Biliary excretion Rate constant DON-3G(nmol/h)
#km_tot      = 4.02;    #Excretion rate constant for DON (nmol/h)
#kelim       = 12.90;    #Excretion rate constant for DON (nmol/h)
#kgutelim   = 0.24; # kgutelim =(kgutabs/Fgutabs) - kgutabs


Initialize {
Oral_dose= PO_dose * BW * 1000 /296.32;
#A_Gi = IngDose;
#km_tot = km_d15g + km_d3g;
#kelim = ke_don + km_tot + kgutelim;
#kgutelim =(kgutabs/Fgutabs) - kgutabs;
}

Dynamics { 
  dt (Q_DON_cpt)  = R_ing * Ka - Q_DON_cpt * (km_d3g + km_d15g + ke_don);
  dt (Q_d3g_cpt) = Q_DON_cpt * km_d3g - Q_d3g_cpt * ku_d3g;
  dt (Q_d15g_cpt) = Q_DON_cpt * km_d15g - Q_d15g_cpt *ku_d15g;
  dt (Q_DON_ex) = Q_DON_cpt * ke_don;
  dt (Q_d3g_ex) = Q_DON_cpt * ku_d3g;
  dt (Q_d15g_ex) = Q_DON_cpt * ku_d15g;
}
CalcOutputs { 
  C_don_cpt = Q_DON_cpt  / Vdist;
  C_d3g_cpt = Q_d3g_cpt / Vdist;
  C_d15g_cpt = Q_d15g_cpt /Vdist;
  Q_total = Q_DON_ex + Q_d3g_ex + Q_d15g_ex;
  Pct_don_ex = (Q_DON_cpt * ke_don / Oral_dose) *100;
  Pct_d3g_ex = (Q_d3g_ex / Oral_dose) *100;
  Pct_d15g_ex = (Q_d15g_ex / Oral_dose) *100;
  
  #ExRate_don = Q_DON_cpt * ke_don; 
  #ExRate_d15g = Q_d15g_cpt * ku_d15g; 
  #ExRate_d3g  = Q_d3g_cpt * ku_d3g;
  #BileExc = A_gblad * kgblad;
  #ExRateFec_don = A_Gi * kgutelim;
  #PctExRate_don  = (Acpt * ke_don/IngDose) * 100;
  #PctExRate_d15g = (Am_1 * ku_d15g/IngDose) * 100; 
  #PctExRate_d3g  = (Am_2 * ku_d3g/IngDose) * 100;
  #PctBileExc = (A_gblad * kgblad/IngDose) * 100;
  #PctFecExc_don = (A_Gi * kgutelim/IngDose) * 100;
}

End.
