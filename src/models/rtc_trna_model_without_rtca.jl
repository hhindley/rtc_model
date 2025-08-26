indexof(sym,syms) = findfirst(isequal(sym),syms)
@independent_variables t 
@parameters L c kr Vmax_init Km_init ω_ab ω_r θtscr g_max θtlr km_a km_b d krep kdam ktag kdeg kin atp na nb nr lam kc k_diss rh thr_t
species_trna1_ΔA = @syms rm_b(t) rtcb(t) rm_r(t) rtcr(t) trna(t) rt(t) 
species_rtc_ΔA = [Symbol(i) for i in species_trna1_ΔA]

D = Differential(t)

@mtkmodel RTC_TRNA_ΔA begin
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
        rh
        thr_t
    end
    @variables begin
        rm_b(t) 
        rtcb(t) 
        rm_r(t) 
        rtcr(t) 
        trna(t) 
        rt(t)

        rhs_rm_b(t) 
        rhs_rtcb(t) 
        rhs_rm_r(t) 
        rhs_rtcr(t) 
        rhs_trna(t) 
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

        tlr_el ~ (g_max*atp/(θtlr+atp)) * trna/(thr_t+trna) 

        # # ribosomes
        Vrep ~ rtcb*rt*krep/(rt+km_b) 
        Vdam ~ kdam*trna 

        rhs_rm_b ~ tscr_ab - lam*(rm_b) - deg(rm_b)
        rhs_rm_r ~ tscr_r - lam*(rm_r) - deg(rm_r)

        rhs_rtcb ~ tlr(rm_b, nb, rh, tlr_el) - dil(rtcb,lam)
        rhs_rtcr ~ tlr(rm_r, nr, rh, tlr_el) - dil(rtcr,lam)

        rhs_trna ~ Vrep - Vdam + kin - dil(trna,lam)
        rhs_rt ~ Vdam - Vrep - dil(rt,lam)

        D(rm_b) ~ rhs_rm_b
        D(rtcb) ~ rhs_rtcb 
        D(rm_r) ~ rhs_rm_r
        D(rtcr) ~ rhs_rtcr
        D(trna) ~ rhs_trna
        D(rt) ~ rhs_rt
    end
end


@mtkbuild rtc_trna_model_ΔA = RTC_TRNA_ΔA()

init_trna_ΔA = [rtc_trna_model_ΔA.rm_b=>0.0,rtc_trna_model_ΔA.rtcb=>0.0,rtc_trna_model_ΔA.rm_r=>0.0,rtc_trna_model_ΔA.rtcr=>0.0,rtc_trna_model_ΔA.trna=>135.5,rtc_trna_model_ΔA.rt=>0.0] # tRNA initial conc = 135.5

tspan = (0, 1e9);

params_trna_ΔA = Dict(L=>L_val, c=>c_val, kr=>kr_val*12, Vmax_init=>Vmax_init_val, Km_init=>Km_init_val, θtscr=>θtscr_val, θtlr=>θtlr_val, na=>nA_val, nb=>nB_val, nr=>nR_val, d=>d_val, krep=>krep_val, ktag=>ktag_val,
atp=>atp_val, km_a=>km_a_val, km_b=>km_b_val, g_max=>g_max_val, kdeg=>kdeg_trna_val, kin=>kin_trna_val, ω_ab=>ω_ab_val, ω_r=>ω_r_val, kdam=>kdam_val, lam=>lam_val, kc=>kc_val, k_diss=>k_diss_val, rh=>rh_val, thr_t=>thr_t_val)

ssvals_trna_ΔA = steady_states(rtc_trna_model_ΔA, init_trna, params_trna)
ssvals_trna_ΔA = collect(ssvals_trna_ΔA)
