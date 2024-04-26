using BifurcationKit

# model functions 
ω_rtcBA(θ_x, a, sig_o) = sig_o * w_BA*a/(θ_x + a)

dil(x,lam) = lam*x

deg(x) = d*x 

tlr(rm_x, nx, rh, tlr_el) = (1/nx)*kc*rh*rm_x*tlr_el 

# solving functions 
function steady_states(model, init, params)
    prob = SteadyStateProblem(model, init, params; jac=true)
    return solve(prob, DynamicSS())
end 

function sol(model, init, tspan, params)
    prob = ODEProblem(model, init, tspan, params; jac=true)
    solu = solve(prob, Rodas4())
    return solu
end

# bifurcation functions 
norminf(x) = norm(x, Inf)

function get_br(model, init, params, kdam_max)
    prob = ODEProblem(model, init, tspan, params; jac=true)
    odefun = prob.f
    F = (u,p) -> odefun(u,p,0)
    J = (u,p) -> odefun.jac(u,p,0)
    id_kdam = indexof(kdam, parameters(model))
    par_tm = prob.p
    # Bifurcation Problem
    if nameof(model) == :rtc_model || nameof(model) == :rtc_inhib_model
        prob = BifurcationProblem(F, prob.u0, setproperties(par_tm), (@lens _[id_kdam]); J=J,
        record_from_solution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], rh = x[7], rd = x[8], rt = x[9]),)
        opts_br = ContinuationPar(p_min = 0., p_max = kdam_max, ds = 0.001, a=0.1,
        dsmax = 0.05, 
        # options to detect bifurcations
        detect_bifurcation = 3, n_inversion = 4, max_bisection_steps = 20,
        # number of eigenvalues
        nev = 2, 
        # maximum number of continuation steps
        max_steps = 50000,)
        # continuation of equilibria
        br = continuation(prob, PALC(θ=0.5), opts_br; plot = false, bothside=true, normC = norminf)
    else
        # print("tRNA rtc model")
        prob = BifurcationProblem(F, prob.u0, setproperties(par_tm), (@lens _[id_kdam]); J=J,
        record_from_solution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], trna = x[7], rd = x[8], rt = x[9]),)
        opts_br = ContinuationPar(p_min = 0., p_max = kdam_max, ds = 0.001,
        dsmax = 0.15, dsmin = 0.0001,
        # options to detect bifurcations
        detect_bifurcation = 3, n_inversion = 2, max_bisection_steps = 20, 
        # maximum number of continuation steps
        max_steps = 50000,)
        # continuation of equilibria
        br = continuation(prob, PALC(θ=0.5), opts_br; plot = false, bothside=true, normC = norminf)
    end
    return br
end
