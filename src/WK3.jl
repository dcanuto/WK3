module WK3

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

export odeint
export rkqs
export rkck
export wk3odes
export elastancefn
export SolverParams
export ModelParams
export buildall

end
