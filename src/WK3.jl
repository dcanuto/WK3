module WK3

# odes for 3-element Windkessel w/ ventricular pressure source & diode valve
include("wk3odes.jl")

# ventricular elastance function
include("elastancefn.jl")

# 5th order Runge-Kutta solver w/ adaptive stepsize (from Numerical Recipes)
include("rkck.jl")
include("rkqs.jl")
include("odeint.jl")

export odeint
export rkqs
export rkck
export wk3odes
export elastancefn

end
