PATH = "YOURPATH"

include("$PATH/funcs.jl")
include("$PATH/params.jl")

indexof(sym,syms) = findfirst(isequal(sym),syms)
@variables t 
@parameters L c kr Vmax_init Km_init ω_ab ω_r θtscr g_max θtlr km_a km_b d krep kdam ktag kdeg kin atp na nb nr lam kc k_diss rh thr_t k_inhib1 k_inhib2 inhib
species_trna_inhib1 = @syms rm_a(t) rtca(t) rm_b(t) rtcb(t) rm_r(t) rtcr(t) trna(t) rd(t) rt(t) rtc_i(t)
species_trna_inhib = [Symbol(i) for i in species_trna_inhib1]

D = Differential(t)

function build_trna_inhib_model(inhib_protein1)
    inhib_protein = inhib_protein1
    @mtkmodel RTC_TRNA_INHIB begin
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
            k_inhib1
            k_inhib2
            inhib
        end
        @variables begin
            rm_a(t) 
            rtca(t) 
            rm_b(t) 
            rtcb(t) 
            rm_r(t) 
            rtcr(t) 
            trna(t) 
            rd(t) 
            rt(t) 
            rtc_i(t)

            rhs_rm_a(t) 
            rhs_rtca(t) 
            rhs_rm_b(t) 
            rhs_rtcb(t) 
            rhs_rm_r(t) 
            rhs_rtcr(t) 
            rhs_trna(t) 
            rhs_rd(t) 
            rhs_rt(t)
            rhs_rtc_i(t) 

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
            Vinflux ~ kin* g_max*atp/(θtlr+atp)
            Vtag ~ rtca*rd*ktag/(rd+km_a)

            rhs_rm_a ~ tscr_ab - dil(rm_a,lam) - deg(rm_a)
            rhs_rm_b ~ tscr_ab - dil(rm_b,lam) - deg(rm_b)
            rhs_rm_r ~ tscr_r - dil(rm_r,lam) - deg(rm_r)

            rhs_trna ~ Vrep - Vdam + Vinflux - dil(trna,lam)
            rhs_rd ~ Vdam - Vtag - kdeg*rd - dil(rd,lam)
            rhs_rt ~ Vtag - Vrep - dil(rt,lam)

            if inhib_protein == :rtca_inhib
                rhs_rtc_i ~ k_inhib1*rtca*inhib - k_inhib2*rtc_i - dil(rtc_i,lam)            
                rhs_rtca ~ tlr(rm_a, na, rh, tlr_el) - dil(rtca,lam) - k_inhib1*rtca*inhib + k_inhib2*rtc_i   
                rhs_rtcb ~ tlr(rm_b, nb, rh, tlr_el) - dil(rtcb,lam)
                rhs_rtcr ~ tlr(rm_r, nr, rh, tlr_el) - dil(rtcr,lam)
            elseif inhib_protein == :rtcb_inhib
                rhs_rtc_i ~ k_inhib1*rtcb*inhib - k_inhib2*rtc_i - dil(rtc_i,lam)            
                rhs_rtca ~ tlr(rm_a, na, rh, tlr_el) - dil(rtca,lam)
                rhs_rtcb ~ tlr(rm_b, nb, rh, tlr_el) - dil(rtcb,lam) - k_inhib1*rtcb*inhib + k_inhib2*rtc_i   
                rhs_rtcr ~ tlr(rm_r, nr, rh, tlr_el) - dil(rtcr,lam)
            else 
                rhs_rtc_i ~ k_inhib1*rtcr*inhib - k_inhib2*rtc_i - dil(rtc_i,lam)            
                rhs_rtca ~ tlr(rm_a, na, rh, tlr_el) - dil(rtca,lam)
                rhs_rtcb ~ tlr(rm_b, nb, rh, tlr_el) - dil(rtcb,lam) 
                rhs_rtcr ~ tlr(rm_r, nr, rh, tlr_el) - dil(rtcr,lam) - k_inhib1*rtcr*inhib + k_inhib2*rtc_i   
            end

            D(rm_a) ~ rhs_rm_a
            D(rtca) ~ rhs_rtca
            D(rm_b) ~ rhs_rm_b
            D(rtcb) ~ rhs_rtcb 
            D(rm_r) ~ rhs_rm_r
            D(rtcr) ~ rhs_rtcr
            D(trna) ~ rhs_trna
            D(rd) ~ rhs_rd
            D(rt) ~ rhs_rt
            D(rtc_i) ~ rhs_rtc_i
        end
    end

    return @mtkbuild rtc_trna_inhib_model = RTC_TRNA_INHIB()
end 

rtca_trna_inhib_model = build_trna_inhib_model(:rtca_inhib)
rtcb_trna_inhib_model = build_trna_inhib_model(:rtcb_inhib)
rtcr_trna_inhib_model = build_trna_inhib_model(:rtcr_inhib)


params_trna_inhib = Dict(L=>L_val, c=>c_val, kr=>kr_val*12, Vmax_init=>Vmax_init_val, Km_init=>Km_init_val, θtscr=>θtscr_val, θtlr=>θtlr_val, na=>nA_val, nb=>nB_val, nr=>nR_val, d=>d_val, krep=>krep_val, ktag=>ktag_val,
atp=>atp_val, km_a=>km_a_val, km_b=>km_b_val, g_max=>g_max_val, kdeg=>kdeg_trna_val, kin=>kin_trna_val, ω_ab=>ω_ab_val, ω_r=>ω_r_val, kdam=>kdam_val, lam=>lam_val, kc=>kc_val, k_diss=>k_diss_val, rh=>rh_val, thr_t=>thr_t_val, k_inhib1=>k_inhib1_val_trna, k_inhib2=>k_inhib2_val_trna, inhib=>inhib_val_trna)

init_trna_inhib_rtca = [rtca_trna_inhib_model.rm_a=>0.0,rtca_trna_inhib_model.rtca=>0.0,rtca_trna_inhib_model.rm_b=>0.0,rtca_trna_inhib_model.rtcb=>0.0,rtca_trna_inhib_model.rm_r=>0.0,rtca_trna_inhib_model.rtcr=>0.0,rtca_trna_inhib_model.trna=>135.5,rtca_trna_inhib_model.rd=>0.0,rtca_trna_inhib_model.rt=>0.0,rtca_trna_inhib_model.rtc_i=>0.0]
init_trna_inhib_rtcb = [rtcb_trna_inhib_model.rm_a=>0.0,rtcb_trna_inhib_model.rtca=>0.0,rtcb_trna_inhib_model.rm_b=>0.0,rtcb_trna_inhib_model.rtcb=>0.0,rtcb_trna_inhib_model.rm_r=>0.0,rtcb_trna_inhib_model.rtcr=>0.0,rtcb_trna_inhib_model.trna=>135.5,rtcb_trna_inhib_model.rd=>0.0,rtcb_trna_inhib_model.rt=>0.0,rtcb_trna_inhib_model.rtc_i=>0.0]
init_trna_inhib_rtcr = [rtcr_trna_inhib_model.rm_a=>0.0,rtcr_trna_inhib_model.rtca=>0.0,rtcr_trna_inhib_model.rm_b=>0.0,rtcr_trna_inhib_model.rtcb=>0.0,rtcr_trna_inhib_model.rm_r=>0.0,rtcr_trna_inhib_model.rtcr=>0.0,rtcr_trna_inhib_model.trna=>135.5,rtcr_trna_inhib_model.rd=>0.0,rtcr_trna_inhib_model.rt=>0.0,rtcr_trna_inhib_model.rtc_i=>0.0]


ssvals_trna_rtca_inhib = steady_states(rtca_trna_inhib_model, init_trna_inhib_rtca, params_trna_inhib)
ssvals_trna_rtcb_inhib = steady_states(rtcb_trna_inhib_model, init_trna_inhib_rtcb, params_trna_inhib)
ssvals_trna_rtcr_inhib = steady_states(rtcr_trna_inhib_model, init_trna_inhib_rtcr, params_trna_inhib)

