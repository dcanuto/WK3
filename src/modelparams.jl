type ModelParams
    Zc::Float64
    R::Float64
    Rav::Float64
    C::Float64
    V0::Float64
    t1::Float64
    t2::Float64
    m1::Float64
    m2::Float64
    Emax::Float64
    Emin::Float64
    th::Vector{Float64}

    function ModelParams(old=Dict("a"=>0),restart="no")
        this = new()
        if restart == "no"
            this.Zc = 2e6; # Pa⋅s/m^3
            this.R = 1.3e8;
            this.Rav = 0;
            this.C = 1.35e-8 # m^3/Pa;
            this.V0 = 10*cm3Tom3;
            this.t1 = 0.215;
            this.t2 = 0.362;
            this.m1 = 1.32;
            this.m2 = 27.4;
            this.Emax = 3.5e8;
            this.Emin = 3.77e6;
            this.th = [0.8];
        elseif restart == "yes"
            this.Zc = old["Zc"];
            this.R = old["R"];
            this.Rav = old["Rav"];
            this.C = old["C"]; # m^3/Pa
            this.V0 = old["V0"];
            this.t1 = old["t1"];
            this.t2 = old["t2"];
            this.m1 = old["m1"];
            this.m2 = old["m2"];
            this.Emax = old["Emax"];
            this.Emin = old["Emin"];
            this.th = [old["th"]];
        end
        return this
    end
end
