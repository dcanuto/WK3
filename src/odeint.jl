function odeint(y0::Vector{Float64},nvar::Int64,t0::Float64,tf::Float64,dtsav::Float64,
    eps::Float64,h0::Float64,hmin::Float64,derivs::Function,stepper::Function)
    # parameters for max integration steps, dependent variable scaling, and saving
    maxstp = 10000;
    tiny = 1e-30;
    saveflag = "yes";

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
    if saveflag == "yes"
        tsav = t-dtsav*2;
    end

    for i = 1:maxstp
        # evaluate derivatives
        derivs(t,y,dy);
        # scaling to monitor accuracy
        ys = abs.(y) + abs.(dy*h)+tiny;
        # store intermediate results
        if saveflag == "yes" && abs.(t-tsav) > abs(dtsav)
            push!(tout,t)
            push!(yout,y[:])
            tsav = t;
        end
        # decrease stepsize if overshooting end of integration window
        if (t+h-tf)*(t+h-t0) > 0
            h = tf-t;
        end
        # attempt an integration step
        t,hdid,hnext=stepper(y,dy,nvar,t,h,hdid,hnext,eps,ys,derivs);
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
            error("Step size below minimum allowable in odeint. hnext = $hnext.")
        end
    end
    error("Too many integration steps in odeint.")
end
