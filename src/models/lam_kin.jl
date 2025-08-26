indexof(sym, syms) = findfirst(isequal(sym),syms)

@independent_variables t 
@parameters L c kr Vmax_init Km_init ω_ab ω_r θtscr g_max θtlr km_a km_b d krep kdam ktag kdeg atp na nb nr lam_c kc k_diss rh_max kin_c
species_rtc1 = @syms rm_a(t) rtca(t) rm_b(t) rtcb(t) rm_r(t) rtcr(t) rh(t) rd(t) rt(t) 
species_rtc = [Symbol(i) for i in species_rtc1]
  
D = Differential(t)

@mtkmodel LAM_KIN begin
    @parameters begin
        L 
        c 
        kr
        Vmax_init 
        Km_init 
        ω_ab  
        ω_r 
        θtscr
        g_max
        θtlr 
        km_a 
        km_b 
        d 
        krep 
        kdam 
        ktag 
        kdeg 
        atp 
        na 
        nb 
        nr 
        lam_c
        kc 
        k_diss 
        rh_max
        kin_c
    end
    @variables begin
        rm_a(t) 
        rtca(t) 
        rm_b(t) 
        rtcb(t) 
        rm_r(t) 
        rtcr(t) 
        rh(t) 
        rd(t) 
        rt(t) 

        rhs_rm_a(t) 
        rhs_rtca(t) 
        rhs_rm_b(t) 
        rhs_rtcb(t) 
        rhs_rm_r(t) 
        rhs_rtcr(t) 
        rhs_rh(t) 
        rhs_rd(t) 
        rhs_rt(t) 

        alpha(t)
        fa(t)
        ra(t)
        Voc(t)
        sig_o(t)
        tscr_ab(t)
        tscr_r(t)
        tlr_el(t)
        Vrep(t)
        Vdam(t)
        Vinflux(t)
        Vtag(t)
        lam(t)
        kin(t)
        

    end

    @equations begin
        alpha ~ rt/kr 
        fa ~ (1+alpha)^6/(L*((1+c*alpha)^6)+(1+alpha)^6)
        ra ~ fa*rtcr 
        
        # transcription
        Voc ~ Vmax_init*atp/(Km_init+atp) 
        sig_o ~ ra*Voc/k_diss 
    
        tscr_ab ~ sig_o*ω_ab*atp/(θtscr+atp)
        tscr_r ~ ω_r*atp/(θtscr+atp) 

        tlr_el ~ g_max*atp/(θtlr+atp)

        lam ~ lam_c*(rh_max - rh)
        kin ~ (rh*kin_c)/(1+rh)

        # # ribosomes
        Vrep ~ rtcb*rt*krep/(rt+km_b) 
        Vdam ~ kdam*rh 
        Vinflux ~ kin 
        Vtag ~ rtca*rd*ktag/(rd+km_a) 

        rhs_rm_a ~ tscr_ab - dil(rm_a,lam) - deg(rm_a)
        rhs_rtca ~ tlr(rm_a, na, rh, tlr_el) - dil(rtca,lam)     
        rhs_rm_b ~ tscr_ab - dil(rm_b,lam) - deg(rm_b)
        rhs_rtcb ~ tlr(rm_b, nb, rh, tlr_el) - dil(rtcb,lam)
        rhs_rm_r ~ tscr_r - dil(rm_r,lam) - deg(rm_r)
        rhs_rtcr ~ tlr(rm_r, nr, rh, tlr_el) - dil(rtcr,lam)
        rhs_rh ~ Vrep - Vdam + Vinflux - dil(rh,lam)
        rhs_rd ~ Vdam - Vtag - kdeg*rd - dil(rd,lam)
        rhs_rt ~ Vtag - Vrep - dil(rt,lam)

        D(rm_a) ~ rhs_rm_a
        D(rtca) ~ rhs_rtca
        D(rm_b) ~ rhs_rm_b
        D(rtcb) ~ rhs_rtcb 
        D(rm_r) ~ rhs_rm_r
        D(rtcr) ~ rhs_rtcr
        D(rh) ~ rhs_rh
        D(rd) ~ rhs_rd
        D(rt) ~ rhs_rt
    end
end

@mtkbuild lam_kin = LAM_KIN()

lam_c_val = 0.00068;
rh_max_val = 38.8;
kin_c_val = 0.269;

params_rtc_lamkin = OrderedDict(L=>L_val, c=>c_val, kr=>kr_val, Vmax_init=>Vmax_init_val, Km_init=>Km_init_val, θtscr=>θtscr_val, θtlr=>θtlr_val, na=>nA_val, nb=>nB_val, nr=>nR_val, d=>d_val, krep=>krep_val, ktag=>ktag_val,
atp=>atp_val, km_a=>km_a_val, km_b=>km_b_val, g_max=>g_max_val, kdeg=>kdeg_val, ω_ab=>ω_ab_val, ω_r=>ω_r_val, kdam=>kdam_val, lam_c=>lam_c_val, kc=>kc_val, k_diss=>k_diss_val, rh_max=>rh_max_val, kin_c=>kin_c_val)

init_rtc_lamkin = [lam_kin.rm_a=>0.0,lam_kin.rtca=>0.0,lam_kin.rm_b=>0.0,lam_kin.rtcb=>0.0,lam_kin.rm_r=>0.0,lam_kin.rtcr=>0.0,lam_kin.rh=>11.29,lam_kin.rd=>0.0,lam_kin.rt=>0.0]

ssvals_rtc_lamkin = steady_states(lam_kin, init_rtc_lamkin, params_rtc_lamkin)
ssvals_rtc_lamkin = collect(ssvals_rtc_lamkin)

