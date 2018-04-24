function odeint(y0::Vector{Float64},nvar::Int64,t0::Float64,tf::Float64,dtsav::Float64,
    eps::Float64,h0::Float64,hmin::Float64,derivs::Function,stepper::Function)
    # parameters for max integration steps, dependent variable scaling
    maxstp = 10000;
    tiny = 1e-30;

    # integration options
    saveflag = true; # save at intervals ≈ dtsav
    changeflag = true; # force changes to state vector outside of integrator

    # output variables
    tout = Float64[];
    yout = Vector{Float64}[];

    # initialization
    y = zeros(nvar);
    ys = zeros(nvar);
    dy = zeros(nvar);
    hdid = 0.;
    hnext = 0.;
    t = t0;
    y[:] = y0;
    nok = 0;
    nbad = 0;

    # make sure we're stepping in the right direction (forwards vs. backwards)
    h = (tf-t0) >= 0 ? (h0 >= 0 ? h0: -h0) : (h0 >= 0 ? -h0 : h0);

    # assure storage of first step
    if saveflag
        tsav = t-dtsav*2;
    end

    # calculate elastance scaling
    τ1 = 0.215;
    τ2 = 0.362;
    m1 = 1.32;
    m2 = 27.4;
    Emax = 3.5e8;
    Emin = 3.77e6;
    th = [0.8];
    tm = linspace(0.,th[end],1e4);
    g1 = (tm/τ1).^m1;
    g2 = (tm/τ2).^m2;
    h1 = g1./(1+g1);
    h2 = 1./(1+g2);
    k = (Emax-Emin)/maximum(h1.*h2);

    for i = 1:maxstp
        # evaluate derivatives
        if changeflag
            derivs(t,y,dy,k);
        else
            derivs(t,y,dy);
        end
        # scaling to monitor accuracy
        ys = abs.(y) + abs.(dy*h)+tiny;
        # store intermediate results
        if saveflag && abs.(t-tsav) > abs(dtsav)
            push!(tout,t)
            push!(yout,y[:])
            tsav = t;
        end
        # decrease stepsize if overshooting end of integration window
        if (t+h-tf)*(t+h-t0) > 0
            h = tf-t;
        end
        # attempt an integration step
        if changeflag
            t,hdid,hnext=stepper(y,dy,nvar,t,h,hdid,hnext,eps,ys,derivs,k);
        else
            t,hdid,hnext=stepper(y,dy,nvar,t,h,hdid,hnext,eps,ys,derivs);
        end
        # force changes to state vector to enforce physical constraints
        if changeflag
            # conversions
            mmHgToPa = 133.32;
            cm3Tom3 = 1e-6;

            # non-dimensional scalings
            Ps = 120*mmHgToPa;
            Vs = 125*cm3Tom3;
            Qs = 500*cm3Tom3;
            ts = 0.8;

            # model parameters
            V0 = 10*cm3Tom3;

            # check valve state, set ventricular flow accordingly
            E,~ = elastancefn(t*ts,[ts],k);
            Pv = E*(y[2]*Vs-V0)/Ps;
            if Pv < y[1] && y[3] <= 0.
                y[3] = 0.;
            end
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
        if abs.(hnext) <= hmin
            push!(yout,y)
            push!(tout,t)
            println("Step size below minimum allowable in odeint. Aborting. hnext = $hnext. t = $t.")
            return tout,yout
        end
    end
    error("Too many integration steps in odeint.")
end
