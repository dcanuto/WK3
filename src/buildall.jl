type CVSystem
    Pa::Vector{Float64}
    Pv::Vector{Float64}
    V::Vector{Float64}
    E::Vector{Float64}
    Q::Vector{Float64}
    t::Vector{Float64}
    sparams::SolverParams
    mparams::ModelParams

    function CVSystem(nvar=0,old=Dict("a"=>0),restart="no")
        this = new()
        if restart == "no"
            this.mparams = ModelParams();
        elseif restart == "yes"
            mparams = old["system"]["mparams"];
            this.mparams = ModelParams(mparams,restart);
        end
        this.Pa = [];
        this.V = [];
        this.Q = [];
        this.Pv = [];
        this.E = [];
        this.t = [];
        this.sparams = SolverParams(nvar);
        return this
    end
end
