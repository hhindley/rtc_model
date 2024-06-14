# rtc model params
L_val = 10.; 
c_val = 0.001; 
kr_val = 0.125; 
Vmax_init_val = 39.51; 
Km_init_val = 250.;
θtscr_val = 160.01;
θtlr_val = 255.73; 
nA_val = 338.; 
nB_val = 408.; 
nR_val = (532.)*6; 
d_val = 0.2; 
krep_val = 137.; 
ktag_val = 9780.;
atp_val = 3000.;
km_a_val = 20.; 
km_b_val = 16.; 
g_max_val = 1260.; 
kdeg_val = 0.001;
kin_val = 0.00022;
ω_ab_val = 1.e-5;  
ω_r_val = 1.e-6; 
kdam_val =  0.0
lam_val = 0.014; 
sf = 1e6/(6.022e23*1e-15);
kb_gm_val = 1/sf 
ku_gm_val = 1. ;
kc_val = 0.6  ;
k_diss_val = 0.006;

# trna params
rh_val = 30;
thr_t_val = 4;
kin_trna_val = 0.00175;
kdeg_trna_val = 0.00001;

k_inhib2_val_trna = 0.0025;
inhib_val_trna = 0.1;

k_inhib1_val_trna = 0.5;
k_inhib1a_val_trna = 0.15;
k_inhib1b_val_trna = 0.25;

k_trna_inhib_vals = [k_inhib1a_val_trna, k_inhib1b_val_trna, k_inhib1_val_trna];


# inhib params for rtc model 
k_inhib1_val = 1;
k_inhib2_val = 0.0025;
inhib_val = 0.1;
k_inhib1a_val = 0.3;
k_inhib1b_val = 0.5;

k_inhib_vals = [k_inhib1a_val, k_inhib1b_val, k_inhib1_val];


tspan = (0, 1e9);