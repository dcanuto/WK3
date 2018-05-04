using WK3

function main()
    # option flags
    rstflag = "yes"
    
    # initialization
    numbeats = 10;
    if rstflag == "no"
        y0 = [80*WK3.mmHgToPa/WK3.Ps;125*WK3.cm3Tom3/WK3.Vs;0*WK3.cm3Tom3/WK3.Qs];
        system = WK3.CVSystem(length(y0));
    elseif rstflag == "yes"
        filename = "system.mat";
        vars = MAT.matread(filename);
        system = WK3.CVSystem(3,vars,rstflag);
        y0 = [vars["system"]["Pa"][end];1;vars["system"]["Q"][end]];
    end
    t0 = 0/WK3.ts;
    tf = [0.8/WK3.ts for i in 1:numbeats];

    for k in 1:length(tf)
        # solver loop
        to,yo = odeint(y0,t0,k*tf[k],system.sparams,system.mparams,wk3odes,rkqs);
        t0 += tf[k];
        append!(system.t,to)

        # reshape output to vectors of individual state variables' time series
        for i = 1:length(yo[1])
            for j = 1:length(yo)
                if i == 1
                    push!(system.Pa,yo[j][i])
                elseif i == 2
                    push!(system.V,yo[j][i])
                elseif i == 3
                    push!(system.Q,yo[j][i])
                end
                if j == length(yo)
                    y0[1] = yo[j][1];
                    y0[2] = 1;
                    y0[3] = yo[j][3];
                end
            end
        end
    end

    # diagnostic variables
    tm = linspace(0,WK3.ts,1e4);
    g1 = (tm/system.mparams.t1).^system.mparams.m1;
    g2 = (tm/system.mparams.t2).^system.mparams.m2;
    h1 = g1./(1+g1);
    h2 = 1./(1+g2);
    k = (system.mparams.Emax-system.mparams.Emin)/maximum(h1.*h2);
    for i = 1:length(system.t)
        g1t = (mod(system.t[i]*WK3.ts,system.mparams.th[end])/system.mparams.t1).^system.mparams.m1;
        g2t = (mod(system.t[i]*WK3.ts,system.mparams.th[end])/system.mparams.t2).^system.mparams.m2;
        h1t = g1t/(1+g1t);
        h2t = 1/(1+g2t);
        push!(system.E,k*h1t*h2t+system.mparams.Emin);
    end
    system.Pv = system.E.*(system.V*WK3.Vs-system.mparams.V0)/WK3.Ps;

    # save for post-processing
    filename = "system2.mat";
    file = MAT.matopen(filename,"w");
    write(file,"system",system)
    close(file)

    # output
    return system
end
