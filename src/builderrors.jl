type Errors
    a::Float64 # parameter distribution smoothing
    h::Float64

    odev::Vector{Float64} # σ for observation, param distributions
    pdev::Vector{Float64}

    lb::Vector{Float64} # lower/upper bounds for truncated normal param distrs.
    ub::Vector{Float64}

    function Errors(nparams=9)
        this = new();
        δ = 0.985;
        this.a = (3*δ-1)/(2*δ);
        this.h = sqrt.(1-this.a^2);
        this.odev = [7.];
        this.pdev = [2e4,1e6,1e-8,1*cm3Tom3,1e-1,2e-1,2e-1,4.5,2.8e7];
        this.lb = -Inf*ones(nparams);
        this.ub = Inf*ones(nparams);
        this.lb[5] = 0.05;
        this.lb[6] = 0.2;
        this.lb[7] = 0.1;
        this.lb[8] = 4e6;
        this.lb[8] = 2.7;
        this.ub[5] = 0.4;
        this.ub[6] = 0.8;
        this.ub[8] = 2e8;
        return this
    end
end
