function elastancefn(t::Float64,th::Vector{Float64},k::Float64)
    # two-Hill function parameters
    τ1 = 0.215;
    τ2 = 0.362;
    m1 = 1.32;
    m2 = 27.4;
    Emin = 3.77e6;

    # two-Hill elastance
    tp = mod(t,sum(th));
    g1t = (tp/τ1)^m1;
    g2t = (tp/τ2)^m2;
    h1t = g1t/(1+g1t);
    h2t = 1/(1+g2t);
    E = k*h1t*h2t + Emin;

    # dE/dt
    g1p = m1/τ1*(tp/τ1)^(m1-1);
    g2p = m2/τ2*(tp/τ2)^(m2-1);
    h1p = g1p/((1+g1t)^2);
    h2p = -g2p/((1+g2t)^2);
    dE = k*(h1p*h2t + h2p*h1t);

    return E,dE
end
