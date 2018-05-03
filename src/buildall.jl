type CVSystem
    Pa::Vector{Float64}
    Pv::Vector{Float64}
    V::Vector{Float64}
    E::Vector{Float64}
    Q::Vector{Float64}
    t::Vector{Float64}
    sparams::SolverParams
    mparams::ModelParams

    function CVSystem(nvar=0)
        this = new()
        this.Pa = [];
        this.Pv = [];
        this.V = [];
        this.E = [];
        this.Q = [];
        this.t = [];
        this.sparams = SolverParams(nvar);
        this.mparams = ModelParams();
        return this
    end
end
