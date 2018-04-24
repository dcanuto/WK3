function wk3odes(t::Float64,y::Vector{Float64},dy::Vector{Float64},k::Float64)
    # conversions
    mmHgToPa = 133.32;
    cm3Tom3 = 1e-6;

    # non-dimensional scalings
    Ps = 120*mmHgToPa;
    Vs = 125*cm3Tom3;
    Qs = 500*cm3Tom3;
    ts = 0.8;

    # model parameters
    Zc = 2e6; # Pa⋅s/m^3
    R = 1.3e8;
    Rav = 0;
    C = 1.35e-8 # m^3/Pa;
    V0 = 10*cm3Tom3;

    # ventricular elastance and pressure
    E,dE = elastancefn(t*ts,[ts],k);
    Pv = E*(y[2]*Vs-V0)/Ps;

    # derivatives
    if y[3] != 0
        dy[1] = ts/Ps*((1+Zc/R)*C)^-1*((1+Zc/R)*Qs*y[3] + C*Zc/R*(-E*Qs*y[3]+
            dE*(Vs*y[2]-V0)) - Ps*y[1]/R);
    else
        dy[1] = ts/Ps*((1+Zc/R)*C)^-1*(-Ps*y[1]/R);
    end
    if Pv > y[1] || (Pv < y[1] && y[3] > 0)
        dy[2] = -ts/Vs*Qs*y[3];
        dy[3] = ts/(Qs*(Rav+Zc))*(-E*Qs*y[3] + dE*(Vs*y[2]-V0) - (1+Zc/R)*Qs/C*y[3] +
            Ps/(R*C)*y[1]);
    elseif Pv < y[1] && y[3] <= 0
        dy[2] = 0.;
        dy[3] = 0.;
    end

end
