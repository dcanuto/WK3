function rkqs(y::Vector{Float64},dy::Vector{Float64},n::Int64,t::Float64,htry::Float64,
    hdid::Float64,hnext::Float64,eps::Float64,ys::Vector{Float64},derivs::Function,k=0.)
    # quality control parameters
    safety = 0.9;
    pgrow = -0.2;
    pshr = -0.2;
    errcon = 1.89e-4;
    errmax = 0.0;

    # allocators for step size, new final time, error estimate, and solution
    htemp = 0.0;
    tnew = 0.0;
    yerr = zeros(n);
    ytemp = zeros(n);

    # first guess for step size
    h = htry;

    # conversions
    mmHgToPa = 133.32;
    cm3Tom3 = 1e-6;

    # non-dimensional scalings
    ts = 0.8;
    Vs = 125*cm3Tom3;
    Ps = 120*mmHgToPa;

    # model parameters
    V0 = 10*cm3Tom3;

    while true
        # attempt an integration step
        if k != 0
            rkck(y,dy,n,t,h,ytemp,yerr,derivs,k);
        else
            rkck(y,dy,n,t,h,ytemp,yerr,derivs);
        end
        # check for error < tolerance (and possibly exit)
        if k == 0 || abs.(elastancefn(t*ts,[ts],k)[1]*(y[2]*Vs - V0)/Ps - y[1]) > eps
            errmax = 0.0;
            for i = 1:n
                errmax = max(errmax,abs.(yerr[i]/ys[i]));
            end
            errmax /= eps;
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
                Plv - Pa = $(abs.(elastancefn(t*ts,[ts],k)[1]*(y[2]*Vs - V0)/Ps - y[1]))
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
