type Errors
    a::Float64 # parameter distribution smoothing
    h::Float64

    odev::Vector{Float64} # σ for observation, param distributions
    pdev::Vector{Float64}

    lb::Vector{Float64} # lower/upper bounds for truncated normal param distrs.
    ub::Vector{Float64}

    function Errors(nparams=9)
        this = new();
        δ = 0.99; # ∃ (0, 1, typically within 0.95-0.99), lower values = less parameter dispersion
        this.a = (3*δ-1)/(2*δ);
        this.h = sqrt.(1-this.a^2);
        this.odev = [2.5*mmHgToPa,2.5*cm3Tom3,2.5*cm3Tom3];
        this.pdev = [1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12]; # no params
        # this.pdev = [2e4,3e6,1e-9,1e-4,0.025,0.04,0.2,2,1.4e7]; # all params
        # this.pdev = [2e4,2e6,1e-9,1e-12,0.1,0.02,0.2,2,7e6]; # all params - V0
        # this.pdev = [2e4,3e6,1e-9,1e-4,0.025,1e-12,0.2,2,1.4e7]; # all params - t2
        # this.pdev = [2e4,2e6,1e-12,1e-12,0.1,0.1,0.1,2.7,1e-12]; # all params - C - V0 - Emax
        # this.pdev = [2e4,2e6,1e-12,1e-12,0.1,0.1,0.1,2.7,7e6]; # all params - C - V0
        # this.pdev = [2e4,2e6,1e-9,1e-12,1e-12,1e-12,1e-12,1e-12,1e-5]; # Windkessel only
        # this.pdev = [2e4,2e6,1e-9,1e-12,1e-12,1e-12,1e-12,1e-12,1.4e7]; # Windkessel + Emax
        # this.pdev = [2e4,1e-12,1e-9,1e-12,1e-12,1e-12,1e-12,1e-12,7e6]; # Windkessel - R + Emax
        # this.pdev = [2e4,2e6,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,7e6]; # Windkessel - C + Emax
        # this.pdev = [2e4,2e6,1e-9,1e-12,1e-12,0.1,1e-12,1e-12,1.4e7]; # Windkessel + Emax + τ2
        # this.pdev = [2e4,2e6,1e-8,1e-12,1e-12,1e-12,1e-12,2.7,1e-12]; # Windkessel + m2
        # this.pdev = [2e4,2e6,1e-8,1e-12,1e-12,1e-12,1e-12,2.7,7e6]; # Windkessel + Emax + m2
        # this.pdev = [2e4,2e6,1e-12,1e-12,1e-12,1e-12,1e-12,2.7,7e6]; # Windkessel - C + Emax + m2
        # this.pdev = [2e4,2e6,1e-9,1e-12,1e-12,1e-12,0.2,2,1.4e7]; # Windkessel + Emax + m1 + m2
        # this.pdev = [2e4,2e6,1e-9,1e-12,0.1,1e-12,0.2,2,7e6]; # Windkessel + Emax + t1 + m1 + m2
        # this.pdev = [2e4,1e-12,1e-9,1e-12,1e-12,1e-12,0.1,2.7,7e6]; # Windkessel - R + Emax + m1 + m2
        # this.pdev = [2e4,1e-12,1e-9,1e-12,1e-2,1e-12,0.1,2.7,7e6]; # Windkessel - R + Emax + t1 + m1 + m2
        # this.pdev = [2e4,2e6,1e-12,1e-12,1e-12,1e-12,0.1,2.7,7e6]; # Windkessel - C + Emax + m1 + m2
        # this.pdev = [2e4,2e6,1e-12,1e-12,0.2,1e-12,0.1,2.7,7e6]; # Windkessel - C + Emax + m1 + m2 + t1
        # this.pdev = [2e4,2e6,1e-12,1e-12,1e-12,0.2,0.1,2.7,7e6]; # Windkessel - C + Emax + m1 + m2 + Δt
        # this.pdev = [2e4,2e6,1e-12,1e-12,1e-12,1e-12,0.1,2.7,1e-12]; # Windkessel - C + m1 + m2
        # this.pdev = [1e-12,2e6,1e-9,1e-12,1e-12,1e-12,0.1,1e-12,7e6]; # Windkessel - Zc + Emax + m1
        # this.pdev = [1e-12,2e6,1e-12,1e-12,1e-12,1e-12,0.1,2.7,7e6]; # R + Emax + m1 + m2
        # this.pdev = [1e-12,2e6,1e-12,1e-12,1e-12,3e-1,0.1,2.7,7e6]; # R + Emax + m1 + m2 + Δτ
        # this.pdev = [1e-12,1e-12,1e-12,1e-12,0.1,0.1,0.2,2.7,1e-12]; # t1 + t2 + m1 + m2
        # this.pdev = [2e4,1e-12,1e-12,1e-12,1e-12,1e-12,0.1,2.7,7e6]; # Zc + Emax + m1 + m2
        # this.pdev = [1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,0.1,2.7,7e6]; # Emax + m1 + m2
        # this.pdev = [2e4,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,7e6]; # Zc + Emax
        # this.pdev = [1e-12,2e6,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,7e6]; # R + Emax
        # this.pdev = [1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,0.1,1e-12,7e6]; # Emax + m1
        # this.pdev = [1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,2.7,7e6]; # Emax + m2
        # this.pdev = [1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,7e6]; # Emax
        # this.pdev = [1e-12,1e-12,1e-12,1e-4,1e-12,1e-12,1e-12,1e-12,1e-12]; # V0
        this.lb = -Inf*ones(nparams);
        this.ub = Inf*ones(nparams);
        this.lb[1] = 1e3; # Zc
        this.lb[2] = 1e5; # R
        this.lb[3] = 5e-9; # C
        this.lb[5] = 0.01; # τ1
        this.lb[6] = 0.01; # τ2
        this.lb[7] = 1.; # m1
        this.lb[8] = 1; # m2
        this.lb[9] = 1e6; # Emax
        # this.ub[5] = 0.4; # τ1
        # this.ub[6] = 0.55; # Δτ
        this.ub[9] = 2e8; # Emax
        return this
    end
end
