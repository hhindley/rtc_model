# model functions 
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
function get_br(model, init, params; kdam_max=1.5, ds=0.0001, a=0.01, dsmax=0.1, dsmin=0.00001, detect_bifurcation=3, n_inversion=2, max_bisection_steps=20, nev=2, max_steps=50000, θ=0.05, tol=1e-10, tol_stability=1e-8, tol_bisection_eigenvalue=1e-16)
    prob = ODEProblem(model, init, tspan, params; jac=true)
    odefun = prob.f
    F = (u,p) -> odefun(u,p,0)
    J = (u,p) -> odefun.jac(u,p,0)
    par_tm = prob.p[1]

    id_kdam = indexin(0.0, par_tm)[1]

    # Bifurcation Problem
    prob = BifurcationProblem(F, prob.u0, par_tm, id_kdam; J=J)

    opt_newt = NewtonPar(tol = tol)

    opts_br = ContinuationPar(p_min = 0., p_max = kdam_max,
    ds = ds, 
    a = a, 
    dsmax = dsmax, 
    dsmin = dsmin,
    detect_bifurcation = detect_bifurcation, 
    n_inversion = n_inversion, 
    max_bisection_steps = max_bisection_steps,
    nev = nev, 
    max_steps = max_steps,
    newton_options = opt_newt,
    tol_stability = tol_stability,
    tol_bisection_eigenvalue = tol_bisection_eigenvalue,
    )
    br = continuation(prob, PALC(θ=θ), opts_br)

    return br
end