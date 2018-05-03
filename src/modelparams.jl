type ModelParams
    Zc::Float64
    R::Float64
    Rav::Float64
    C::Float64
    V0::Float64
    τ1::Float64
    τ2::Float64
    m1::Float64
    m2::Float64
    Emax::Float64
    Emin::Float64
    th::Vector{Float64}

    function ModelParams()
        this = new()
        this.Zc = 2e6; # Pa⋅s/m^3
        this.R = 1.3e8;
        this.Rav = 0;
        this.C = 1.35e-8 # m^3/Pa;
        this.V0 = 10*cm3Tom3;
        this.τ1 = 0.215;
        this.τ2 = 0.362;
        this.m1 = 1.32;
        this.m2 = 27.4;
        this.Emax = 3.5e8;
        this.Emin = 3.77e6;
        this.th = [0.8];
        return this
    end
end
