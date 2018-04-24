using WK3

function main()
    # conversions
    mmHgToPa = 133.32;
    cm3Tom3 = 1e-6;

    # non-dimensional scalings
    Ps = 120*mmHgToPa;
    Vs = 125*cm3Tom3;
    Qs = 500*cm3Tom3;
    ts = 0.8;

    # model parameters
    V0 = 10*cm3Tom3;

    # initialization
    y0 = [80*mmHgToPa/Ps;125*cm3Tom3/Vs;0*cm3Tom3/Qs];
    t0 = 0/ts;
    tf = 0.8/ts;

    # solver parameters
    nvar = length(y0);
    dtsav = 1e-3;
    h0 = 1e-3;
    hmin = 1e-16;
    eps = 1e-8;

    # solver loop
    tic();
    to,yo = odeint(y0,nvar,t0,tf,dtsav,eps,h0,hmin,wk3odes,rkqs);
    toc();

    # reshape output to vectors of individual state variables' time series
    yvec = Vector{Float64}[];
    for i = 1:length(yo[1])
        yi = Float64[];
        for j = 1:length(yo)
            push!(yi,yo[j][i])
        end
        push!(yvec,yi)
    end

    # diagnostic variables
    E = Float64[];
    τ1 = 0.215;
    τ2 = 0.362;
    m1 = 1.32;
    m2 = 27.4;
    Emax = 3.5e8;
    Emin = 3.77e6;
    tm = linspace(0,ts,1e4);
    g1 = (tm/τ1).^m1;
    g2 = (tm/τ2).^m2;
    h1 = g1./(1+g1);
    h2 = 1./(1+g2);
    k = (Emax-Emin)/maximum(h1.*h2);
    for i = 1:length(to)
        g1t = (to[i]*ts/τ1).^m1;
        g2t = (to[i]*ts/τ2).^m2;
        h1t = g1t/(1+g1t);
        h2t = 1/(1+g2t);
        push!(E,k*h1t*h2t+Emin);
    end
    Pv = E.*(yvec[2]*Vs-V0)/Ps;

    # output
    return to,yvec,Pv
end
