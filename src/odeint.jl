function odeint(y0::Vector{Float64},t0::Float64,tf::Float64,sparams::WK3.SolverParams,
    mparams::WK3.ModelParams,derivs::Function,stepper::Function)
    # parameters for dependent variable scaling
    tiny = 1e-30;

    # integration options
    saveflag = true; # save at intervals ≈ dtsav

    # output variables
    tout = Float64[];
    yout = Vector{Float64}[];

    # initialization
    y = zeros(sparams.nvar);
    ys = zeros(sparams.nvar);
    dy = zeros(sparams.nvar);
    hdid = 0.;
    hnext = 0.;
    t = t0;
    y[:] = y0;
    nok = 0;
    nbad = 0;

    # make sure we're stepping in the right direction (forwards vs. backwards)
    h = (tf-t0) >= 0 ? (sparams.h0 >= 0 ? sparams.h0: -sparams.h0) : (sparams.h0 >= 0 ? -sparams.h0 : sparams.h0);

    # assure storage of first step
    if saveflag
        tsav = t-sparams.dtsav*2;
    end

    # calculate elastance scaling
    tm = linspace(0.,mparams.th[end],1e4);
    g1 = (tm/mparams.τ1).^mparams.m1;
    g2 = (tm/mparams.τ2).^mparams.m2;
    h1 = g1./(1+g1);
    h2 = 1./(1+g2);
    k = (mparams.Emax-mparams.Emin)/maximum(h1.*h2);

    while true
        # evaluate derivatives
        derivs(t,y,dy,mparams,k);
        # scaling to monitor accuracy
        ys = abs.(y) + abs.(dy*h)+tiny;
        # store intermediate results
        if saveflag && abs.(t-tsav) > abs(sparams.dtsav)
            push!(tout,t)
            push!(yout,y[:])
            tsav = t;
        end
        # decrease stepsize if overshooting end of integration window
        if (t+h-tf)*(t+h-t0) > 0
            h = tf-t;
        end
        # attempt an integration step
        t,hdid,hnext=stepper(y,dy,t,h,hdid,hnext,ys,derivs,sparams,mparams,k);
        # check valve state, set ventricular flow accordingly
        E,~ = elastancefn(t*ts,mparams,k);
        Pv = E*(y[2]*Vs-mparams.V0)/Ps;
        if Pv < y[1] && y[3] <= 0.
            y[3] = 0.;
        end
        # track times step size was changed
        if hdid == h
            nok+=1;
        else
            nbad+=1;
        end
        # normal exit if finished (saving final step)
        if (t-tf)*(tf-t0) >= 0
            push!(yout,y)
            push!(tout,t)
            return tout,yout
        end
        if abs.(hnext) <= sparams.hmin
            push!(yout,y)
            push!(tout,t)
            println("Step size below minimum allowable in odeint. Aborting. hnext = $hnext. t = $t.")
            return tout,yout
        end
    end
end
