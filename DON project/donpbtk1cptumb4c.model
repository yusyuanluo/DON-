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
# Time:    minute
# Quantity:   nmol

# Conversions between mg and nMol:
#Dose = Dose/(1000*MW_don); in number of moles
#IngDose = Dose * 1e+9; # In nmol

States  = { 
  Acpt,      # Quantity in central compartment (nmol)
  Au_don,    # Quantity of don in urine (nmol)
  Au_d15g,   # Quantity of d15g in urine  (nmol)
  Au_d3g,     # Quantity of d3g in kidney (nmol)
  Am_1,       # Quantity of d15g metabolite (nmol)
  Am_2,        # Quantity of d3g metabolite (nmol)
  A_Gi,      # Quantity of in the GI compartment (nmol)
  A_fec,      # Quantity of free DON in feces (nmol)
  A_gblad     # Quantity of DON in the Gall bladder (nmol)
};  

# Gall bladder input modeling

Inputs = { kgblad
}; 

Outputs = {
  Ccpt, ExRate_don, ExRate_d15g, ExRate_d3g, PctExRate_don, PctExRate_d15g, PctExRate_d3g
};
# Oral input modeling
IngDose    = 236.33; # ingested input (nmol)
Fgutabs    = 0.56; #
kgutabs    = 0.09; # Intestinal absortion rate (/h)
kgblad_mag = 8.00;
kgblad_time = 15; 
kgblad =    PerDose(kgblad_mag, 48, kgblad_time, 1);

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
kbile_1    = 1.00;     #Biliary excretion Rate constant DON-15G(nmol/h)
kbile_2    = 1.00;     #Biliary excretion Rate constant DON-3G(nmol/h)
km_tot      = 4.02;    #Excretion rate constant for DON (nmol/h)
kelim       = 12.90;    #Excretion rate constant for DON (nmol/h)
kgutelim   = 0.24; # kgutelim =(kgutabs/Fgutabs) - kgutabs


Initialize {
A_Gi = IngDose;
km_tot = km_d15g + km_d3g;
kelim = ke_don + km_tot + kgutelim;
kgutelim =(kgutabs/Fgutabs) - kgutabs;
}

Dynamics { 
  dt (A_Gi)  = - A_Gi * (kgutabs + kgutelim) + A_gblad * kgblad;
  dt (Acpt) = A_Gi * kgutabs - kelim * Acpt;
  dt (A_fec) = A_Gi * kgutelim;
  dt (A_gblad) = (Am_1 * kbile_1 + Am_2 * kbile_2) - A_gblad * kgblad;
  dt (Au_don) = Acpt * ke_don;
  dt (Am_1) = (Acpt * km_d15g) - Am_1 * (ku_d15g + kbile_1);
  dt (Am_2) = (Acpt * km_d3g) - Am_2 * (ku_d3g + kbile_2);
  dt (Au_d15g) = Am_1 * ku_d15g;
  dt (Au_d3g) = Am_2 * ku_d3g;
}

CalcOutputs { 
  Ccpt = Acpt  / Vdist;
  ExRate_don = Acpt * ke_don; 
  ExRate_d15g = Am_1 * ku_d15g; 
  ExRate_d3g  = Am_2 * ku_d3g;
  #BileExc = A_gblad * kgblad;
  #ExRateFec_don = A_Gi * kgutelim;
  PctExRate_don  = (Acpt * ke_don/IngDose) * 100;
  PctExRate_d15g = (Am_1 * ku_d15g/IngDose) * 100; 
  PctExRate_d3g  = (Am_2 * ku_d3g/IngDose) * 100;
  #PctBileExc = (A_gblad * kgblad/IngDose) * 100;
  #PctFecExc_don = (A_Gi * kgutelim/IngDose) * 100;
}

End.
