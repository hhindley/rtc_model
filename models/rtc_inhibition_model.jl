PATH = "YOURPATH"

include("$PATH/funcs.jl")
include("$PATH/params.jl")

indexof(sym, syms) = findfirst(isequal(sym),syms)

@variables t 
@parameters L c kr Vmax_init Km_init ω_ab ω_r θtscr g_max θtlr km_a km_b d krep kdam ktag kdeg kin atp na nb nr lam kc k_diss k_inhib1 k_inhib2 inhib
species_inhib1 = @syms rm_a(t) rtca(t) rm_b(t) rtcb(t) rm_r(t) rtcr(t) rh(t) rd(t) rt(t) rtc_i(t)
species_inhib = [Symbol(i) for i in species_inhib1]
    
D = Differential(t)

function build_inhib_model(inhib_protein1)
    inhib_protein = inhib_protein1
    @mtkmodel RTC_INHIB begin
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
            rh(t) 
            rd(t) 
            rt(t) 
            rtc_i(t)

            rhs_rm_a(t) 
            rhs_rtca(t) 
            rhs_rm_b(t) 
            rhs_rtcb(t) 
            rhs_rm_r(t) 
            rhs_rtcr(t) 
            rhs_rh(t) 
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

            tlr_el ~ g_max*atp/(θtlr+atp)

            # # ribosomes
            Vrep ~ rtcb*rt*krep/(rt+km_b) 
            Vdam ~ kdam*rh 
            Vinflux ~ kin* g_max*atp/(θtlr+atp) 
            Vtag ~ rtca*rd*ktag/(rd+km_a)

            rhs_rm_a ~ tscr_ab - dil(rm_a,lam) - deg(rm_a)
            rhs_rm_b ~ tscr_ab - dil(rm_b,lam) - deg(rm_b)
            rhs_rm_r ~ tscr_r - dil(rm_r,lam) - deg(rm_r)

            rhs_rh ~ Vrep - Vdam + Vinflux - dil(rh,lam)
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
            D(rh) ~ rhs_rh
            D(rd) ~ rhs_rd
            D(rt) ~ rhs_rt
            D(rtc_i) ~ rhs_rtc_i
        end
    end

    return @mtkbuild rtc_inhib_model = RTC_INHIB()
end 

rtca_inhib_model = build_inhib_model(:rtca_inhib)
rtcb_inhib_model = build_inhib_model(:rtcb_inhib)
rtcr_inhib_model = build_inhib_model(:rtcr_inhib)


params_inhib = Dict(L=>L_val, c=>c_val, kr=>kr_val, Vmax_init=>Vmax_init_val, Km_init=>Km_init_val, θtscr=>θtscr_val, θtlr=>θtlr_val, na=>nA_val, nb=>nB_val, nr=>nR_val, d=>d_val, krep=>krep_val, ktag=>ktag_val,
atp=>atp_val, km_a=>km_a_val, km_b=>km_b_val, g_max=>g_max_val, kdeg=>kdeg_val, kin=>kin_val, ω_ab=>ω_ab_val, ω_r=>ω_r_val, kdam=>kdam_val, lam=>lam_val, kc=>kc_val, k_diss=>k_diss_val, k_inhib1=>k_inhib1_val, k_inhib2=>k_inhib2_val, inhib=>inhib_val)

init_inhib_rtca = [rtca_inhib_model.rm_a=>0.0,rtca_inhib_model.rtca=>0.0,rtca_inhib_model.rm_b=>0.0,rtca_inhib_model.rtcb=>0.0,rtca_inhib_model.rm_r=>0.0,rtca_inhib_model.rtcr=>0.0,rtca_inhib_model.rh=>11.29,rtca_inhib_model.rd=>0.0,rtca_inhib_model.rt=>0.0,rtca_inhib_model.rtc_i=>0.0]
init_inhib_rtcb = [rtcb_inhib_model.rm_a=>0.0,rtcb_inhib_model.rtca=>0.0,rtcb_inhib_model.rm_b=>0.0,rtcb_inhib_model.rtcb=>0.0,rtcb_inhib_model.rm_r=>0.0,rtcb_inhib_model.rtcr=>0.0,rtcb_inhib_model.rh=>11.29,rtcb_inhib_model.rd=>0.0,rtcb_inhib_model.rt=>0.0,rtcb_inhib_model.rtc_i=>0.0]
init_inhib_rtcr = [rtcr_inhib_model.rm_a=>0.0,rtcr_inhib_model.rtca=>0.0,rtcr_inhib_model.rm_b=>0.0,rtcr_inhib_model.rtcb=>0.0,rtcr_inhib_model.rm_r=>0.0,rtcr_inhib_model.rtcr=>0.0,rtcr_inhib_model.rh=>11.29,rtcr_inhib_model.rd=>0.0,rtcr_inhib_model.rt=>0.0,rtcr_inhib_model.rtc_i=>0.0]

ssvals_rtca_inhib = steady_states(rtca_inhib_model, init_inhib_rtca, params_inhib)
ssvals_rtcb_inhib = steady_states(rtcb_inhib_model, init_inhib_rtcb, params_inhib)
ssvals_rtcr_inhib = steady_states(rtcr_inhib_model, init_inhib_rtcr, params_inhib)
