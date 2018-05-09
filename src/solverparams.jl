type SolverParams
nvar::Int64
dtsav::Float64
h0::Float64
h0nom::Float64
hmin::Float64
eps::Float64

    function SolverParams(nvar=0)
        this = new()
        this.nvar = nvar;
        this.dtsav = 1e-3;
        this.h0 = 1e-3;
        this.h0nom = 1e-3;
        this.hmin = 1e-16;
        this.eps = 1e-8;
        return this
    end
end
