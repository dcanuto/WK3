function rkqs(y::Vector{Float64},dy::Vector{Float64},n::Int64,t::Float64,htry::Float64,
    hdid::Float64,hnext::Float64,eps::Float64,ys::Vector{Float64},derivs::Function)
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

    while true
        # attempt an integration step
        rkck(y,dy,n,t,h,ytemp,yerr,derivs);
        # check for error < tolerance (and possibly exit)
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
