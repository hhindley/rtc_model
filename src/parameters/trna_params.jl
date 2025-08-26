const rh_val = 35;
const thr_t_val = 100;
const kinc_trna_val = 0.0022;
const kdeg_trna_val = 1;
const kin_trna_val = kinc_trna_val * g_max_val*atp_val/(θtlr_val+atp_val);


const k_inhib2_val_trna = 0.0025;
const inhib_val_trna = 0.1;

const k_inhib1a_val_trna = 0.05;
const k_inhib1b_val_trna = 0.1;
const k_inhib1_val_trna = 0.2;

k_trna_inhib_vals = [k_inhib1a_val_trna, k_inhib1b_val_trna, k_inhib1_val_trna];
