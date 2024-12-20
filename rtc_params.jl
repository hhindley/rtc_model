
const L_val::Float64 = 50.;  
const c_val::Float64 = 0.01;  
const kr_val::Float64 = 0.125; 
const Vmax_init_val::Float64 = 39.51;  
const Km_init_val::Float64 = 250.; 
const θtscr_val::Float64 = 160.01; 
const θtlr_val::Float64 = 255.73;  
const nA_val::Float64 = 338.; 
const nB_val::Float64 = 408.;  
const nR_val::Float64 = (532.)*6;  
const d_val::Float64 = 0.2;  
const krep_val::Float64 = 15.67 
const ktag_val::Float64 = 9780.; 
const atp_val::Float64 = 3000.; 
const km_a_val::Float64 = 20.;  
const km_b_val::Float64 = 16.; 
const g_max_val::Float64 = 1260.;  
const kdeg_val::Float64 = 1.; 
const kinc_val::Float64 = 0.022/100;
const kin_val = kinc_val * g_max_val*atp_val/(θtlr_val+atp_val);
const ω_ab_val::Float64 = 5e-5  
const ω_r_val::Float64 = 1.e-6;  
const kdam_val::Float64 =  0.0;
const lam_val::Float64 = 0.014; 
const sf::Float64 = 1e6/(6.022e23*1e-15);
const kb_gm_val::Float64 = 1;
const ku_gm_val::Float64 = 1.;
const kc_val::Float64 = 0.6;
const k_diss_val::Float64 = 0.006;

const k_inhib2_val::Float64 = 0.0025;
const inhib_val::Float64 = 0.1;
const k_inhib1a_val::Float64 = 0.3;
const k_inhib1b_val::Float64 = 0.5;
const k_inhib1_val::Float64 = 1.;

k_inhib_vals = [k_inhib1a_val, k_inhib1b_val, k_inhib1_val];


tspan = (0, 1e9);