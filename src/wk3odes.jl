function wk3odes(t::Float64,y::Vector{Float64},dy::Vector{Float64},mparams::WK3.ModelParams,k::Float64)
    # ventricular elastance and pressure

    E,dE = elastancefn(t*ts,mparams,k);
    Pv = E.*(y[2].*Vs.-mparams.V0)./Ps;

    # derivatives
    # println("ODEs, first loop:")
    if y[3] != 0
        dy[1] = ts/Ps*((1+mparams.Zc/mparams.R)*mparams.C)^-1*((1+mparams.Zc/mparams.R)*Qs*y[3]
            .+ mparams.C*mparams.Zc/mparams.R*(-E*Qs*y[3].+dE*(Vs*y[2]-mparams.V0)) .- Ps*y[1]/mparams.R);
    else
        dy[1] = ts/Ps*((1+mparams.Zc/mparams.R)*mparams.C)^-1.*(-Ps.*y[1]./mparams.R);
    end
    # println("ODEs, second loop:")
    if Pv > y[1] || (Pv < y[1] && y[3] > 0)
        dy[2] = -ts/Vs*Qs*y[3];
        dy[3] = ts/(Qs*(mparams.Rav+mparams.Zc))*(-E*Qs*y[3] .+ dE*(Vs*y[2].-mparams.V0) .-
            (1+mparams.Zc/mparams.R)*Qs/mparams.C*y[3] .+ Ps/(mparams.R*mparams.C)*y[1]);
    elseif Pv < y[1] && y[3] <= 0
        dy[2] = 0.;
        dy[3] = 0.;
    end

    # for i = 1:length(dy)
    #     if dy[i] > 1e6
    #         warn("dy/dt above tolerance at t = $t. dE/dt = $dE.")
    #     end
    # end
end
