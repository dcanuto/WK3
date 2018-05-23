module WK3

using MAT
using Interpolations
using Distributions

# conversions
include("constants.jl")

# 5th order Runge-Kutta solver w/ adaptive stepsize (from Numerical Recipes)
include("solverparams.jl")
include("modelparams.jl")
include("rkck.jl")
include("rkqs.jl")
include("odeint.jl")

# type definitions
include("buildall.jl")

# odes for 3-element Windkessel w/ ventricular pressure source & diode valve
include("wk3odes.jl")

# ventricular elastance function
include("elastancefn.jl")

# data assimilation
include("patdatainterp.jl")
include("builderrors.jl")
include("paramwalk.jl")

# Fisher information matrix
include("fdjac.jl")

export odeint
export rkqs
export rkck
export wk3odes
export elastancefn
export SolverParams
export ModelParams
export buildall
export patdatainterp
export builderrors
export paramwalk!
export fdjac

end
