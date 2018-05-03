function rkqs(y::Vector{Float64},dy::Vector{Float64},t::Float64,htry::Float64,
    hdid::Float64,hnext::Float64,ys::Vector{Float64},derivs::Function,sparams::WK3.SolverParams,
    mparams::WK3.ModelParams,k::Float64)
    # quality control parameters
    safety = 0.9;
    pgrow = -0.2;
    pshr = -0.2;
    errcon = 1.89e-4;
    errmax = 0.0;

    # allocators for step size, new final time, error estimate, and solution
    htemp = 0.0;
    tnew = 0.0;
    yerr = zeros(sparams.nvar);
    ytemp = zeros(sparams.nvar);

    # first guess for step size
    h = htry;

    while true
        # attempt an integration step
        rkck(y,dy,sparams.nvar,t,h,ytemp,yerr,derivs,k,mparams);
        # check for error < tolerance (and possibly exit)
        if abs.(elastancefn(t*ts,mparams,k)[1]*(y[2]*Vs - mparams.V0)/Ps - y[1]) > sparams.eps
            errmax = 0.0;
            for i = 1:sparams.nvar
                errmax = max(errmax,abs.(yerr[i]/ys[i]));
            end
            errmax /= sparams.eps;
            if errmax <= 1.0
                break
            end
            # resize step if error > tolerance
            htemp = safety*h*errmax^pshr;
            h = h > 0 ? max(htemp,0.1*h) : min(htemp,0.1*h);
            tnew = t+h;
            if tnew == t
                error("Step size underflow in rkqs.")
            end
        else
            warn("Transvalvular pressure difference under tolerance.
                Plv - Pa = $(abs.(elastancefn(t*ts,mparams,k)[1]*(y[2]*Vs - mparams.V0)/Ps - y[1]))
                at t = $t. Taking reduced-order time step.")
            break
        end
    end

    # guess next step size
    if errmax > errcon
        hnext = safety*h*errmax^pgrow;
    else
        hnext = 5*h;
    end

    # output solution, new time, and step size taken
    y[:] = ytemp;
    t += h;
    hdid = h;
    return t,hdid,hnext
end
