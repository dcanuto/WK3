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
            # healthy
            # # this.Zc = 2e5; # Pa⋅s/m^3
            # # this.R = 2e7;
            # this.Zc = 1.186*2e5; # Pa⋅s/m^3
            # this.R = 0.940*3e7;
            # this.Rav = 2e5;
            # # this.Rav = 0;
            # # this.C = 1.3e-7 # m^3/Pa;
            # this.C = 0.901*1.5e-8 # m^3/Pa;
            # this.V0 = 10*cm3Tom3;
            # # this.t1 = 0.215;
            # this.t1 = 0.913*1;
            # # this.t2 = 0.362;
            # this.t2 = 1.308*0.3;
            # this.m1 = 0.821*2;
            # this.m2 = 0.389*27.4;
            # # this.m2 = 7;
            # # this.Emax = 7e7;
            # this.Emax = 0.968*1e8;
            # patient
            this.Zc = 1.007*2e5; # Pa⋅s/m^3
            this.R = 0.963*3e7;
            this.Rav = 2e5;
            this.C = 0.864*1e-8 # m^3/Pa;
            this.V0 = 10*cm3Tom3;
            this.t1 = 0.896*0.25;
            this.t2 = 0.4;
            this.m1 = 1.066*2;
            this.m2 = 1.037*20;
            this.Emax = 0.764*1.4e8;
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
