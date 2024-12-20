indexof(sym, syms) = findfirst(isequal(sym),syms)

@independent_variables t 
@parameters L c kr Vmax_init Km_init ω_ab ω_r θtscr g_max θtlr km_a km_b d krep kdam ktag kdeg kin atp na nb nr lam kc k_diss 
# species_rtc1 = @syms rm_a(t) rtca(t) rm_b(t) rtcb(t) rm_r(t) rtcr(t) rh(t) rd(t) rt(t) 
species_rtc1_ΔA = @syms rm_b(t) rtcb(t) rm_r(t) rtcr(t) rh(t) rt(t) 
species_rtc_ΔA = [Symbol(i) for i in species_rtc1_ΔA]
  
D = Differential(t)

@mtkmodel RTC_ΔA begin
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
        kin 
        atp 
        na 
        nb 
        nr 
        lam 
        kc 
        k_diss 
    end
    @variables begin
        rm_b(t) 
        rtcb(t) 
        rm_r(t) 
        rtcr(t) 
        rh(t) 
        rt(t) 

        rhs_rm_b(t) 
        rhs_rtcb(t) 
        rhs_rm_r(t) 
        rhs_rtcr(t) 
        rhs_rh(t) 
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

        # # ribosomes
        Vrep ~ rtcb*rt*krep/(rt+km_b)  
        Vdam ~ kdam*rh 

        rhs_rm_b ~ tscr_ab - dil(rm_b,lam) - deg(rm_b)
        rhs_rtcb ~ tlr(rm_b, nb, rh, tlr_el) - dil(rtcb,lam)
        rhs_rm_r ~ tscr_r - dil(rm_r,lam) - deg(rm_r)
        rhs_rtcr ~ tlr(rm_r, nr, rh, tlr_el) - dil(rtcr,lam)
        rhs_rh ~ Vrep - Vdam + kin - dil(rh,lam)
        rhs_rt ~ Vdam - Vrep - dil(rt,lam)

        D(rm_b) ~ rhs_rm_b
        D(rtcb) ~ rhs_rtcb 
        D(rm_r) ~ rhs_rm_r
        D(rtcr) ~ rhs_rtcr
        D(rh) ~ rhs_rh
        D(rt) ~ rhs_rt
    end
end

@mtkbuild rtc_model_ΔA = RTC_ΔA()

init_rtc_ΔA = [rtc_model.rm_b=>0.0,rtc_model.rtcb=>0.0,rtc_model.rm_r=>0.0,rtc_model.rtcr=>0.0,rtc_model.rh=>11.29,rtc_model.rt=>0.0]

params_rtc = OrderedDict(L=>L_val, c=>c_val, kr=>kr_val, Vmax_init=>Vmax_init_val, Km_init=>Km_init_val, θtscr=>θtscr_val, θtlr=>θtlr_val, na=>nA_val, nb=>nB_val, nr=>nR_val, d=>d_val, krep=>krep_val, ktag=>ktag_val,
atp=>atp_val, km_a=>km_a_val, km_b=>km_b_val, g_max=>g_max_val, kdeg=>kdeg_val, kin=>kin_val, ω_ab=>ω_ab_val, ω_r=>ω_r_val, kdam=>kdam_val, lam=>lam_val, kc=>kc_val, k_diss=>k_diss_val)

ssvals_rtc_ΔA = steady_states(rtc_model_ΔA, init_rtc_ΔA, params_rtc)

ssvals_rtc_ΔA = collect(ssvals_rtc_ΔA)

