function elastancefn(t::Float64,mparams::WK3.ModelParams,k::Float64)
    # two-Hill elastance
    tp = mod(t,sum(mparams.th));
    g1t = (tp/mparams.t1)^mparams.m1;
    g2t = (tp/mparams.t2)^mparams.m2;
    h1t = g1t/(1+g1t);
    h2t = 1/(1+g2t);
    E = k*h1t*h2t + mparams.Emin;

    # dE/dt
    g1p = mparams.m1/mparams.t1*(tp/mparams.t1)^(mparams.m1-1);
    g2p = mparams.m2/mparams.t2*(tp/mparams.t2)^(mparams.m2-1);
    h1p = g1p/((1+g1t)^2);
    h2p = -g2p/((1+g2t)^2);
    dE = k*(h1p*h2t + h2p*h1t);

    return E,dE
end
