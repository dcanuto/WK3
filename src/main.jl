using WK3

function main()
    # option flags
    rstflag = "yes"

    # initialization
    numbeats = 3; # total number of cardiac cycles
    ensemblesize = 3;
    if rstflag == "no"
        y0 = [[80*WK3.mmHgToPa/WK3.Ps;125*WK3.cm3Tom3/WK3.Vs;0*WK3.cm3Tom3/WK3.Qs] for i=1:ensemblesize];
        nvar = [3 for i=1:ensemblesize];
        systems = pmap((a1)->WK3.CVSystem(a1),nvar);
    elseif rstflag == "yes"
        fnames = ["system_$i.mat" for i=1:ensemblesize];
        vars = [MAT.matread(fnames[i]) for i=1:ensemblesize];
        nvar = [3 for i=1:ensemblesize];
        rst = [rstflag for i=1:ensemblesize];
        systems = pmap((a1,a2,a3)->WK3.CVSystem(a1,a2,a3),nvar,vars,rst);
        y0 = [[vars[1]["system"]["Pa"][end];1;vars[1]["system"]["Q"][end]] for i=1:ensemblesize];
    end
    t0 = [0/WK3.ts for i=1:ensemblesize];
    sparams = [systems[i].sparams for i=1:ensemblesize];
    mparams = [systems[i].mparams for i=1:ensemblesize];
    to = Vector{Float64}[];
    yo = Vector{Vector{Float64}}[];
    append!(yo,[[zeros(1)]]);
    append!(to,[zeros(1)]);

    # solver loop
    for n in 1:numbeats
        tf = [n*0.8/WK3.ts for i=1:ensemblesize];
        soln = pmap((a1,a2,a3,a4,a5)->odeint(a1,a2,a3,a4,a5),y0,t0,tf,sparams,mparams);
        to = [soln[i][1] for i in 1:length(soln)];
        yo = [soln[i][2] for i in 1:length(soln)];
        t0 += 0.8/WK3.ts;
        for i=1:ensemblesize
            append!(systems[i].t,to[i])
        end

        # reshape output to vectors of individual state variables' time series
        for k = 1:ensemblesize # k indexes ensemble members
            for i = 1:length(yo[k][1]) # i indexes state variables
                for j = 1:length(yo[k]) # j indexes time steps
                    if i == 1
                        push!(systems[k].Pa,yo[k][j][i])
                    elseif i == 2
                        push!(systems[k].V,yo[k][j][i])
                    elseif i == 3
                        push!(systems[k].Q,yo[k][j][i])
                    end
                    if j == length(yo[k])
                        y0[k][1] = yo[k][j][1];
                        y0[k][2] = 1;
                        y0[k][3] = yo[k][j][3];
                    end
                end
            end
        end
    end

    # diagnostic variables
    tm = linspace(0,WK3.ts,1e4);
    for i=1:ensemblesize
        g1 = (tm/systems[i].mparams.t1).^systems[i].mparams.m1;
        g2 = (tm/systems[i].mparams.t2).^systems[i].mparams.m2;
        h1 = g1./(1+g1);
        h2 = 1./(1+g2);
        k = (systems[i].mparams.Emax-systems[i].mparams.Emin)/maximum(h1.*h2);
        for j = 1:length(systems[i].t)
            g1t = (mod(systems[i].t[j]*WK3.ts,systems[i].mparams.th[end])/systems[i].mparams.t1).^systems[i].mparams.m1;
            g2t = (mod(systems[i].t[j]*WK3.ts,systems[i].mparams.th[end])/systems[i].mparams.t2).^systems[i].mparams.m2;
            h1t = g1t/(1+g1t);
            h2t = 1/(1+g2t);
            push!(systems[i].E,k*h1t*h2t+systems[i].mparams.Emin);
        end
        systems[i].Pv = systems[i].E.*(systems[i].V*WK3.Vs-systems[i].mparams.V0)/WK3.Ps;
    end

    # save for post-processing
    for i=1:ensemblesize
        filename = "system2_$i.mat"; # name of solution save file
        file = MAT.matopen(filename,"w");
        write(file,"system",systems[i])
        close(file)
    end

    # output
    return systems
end
