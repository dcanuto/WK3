using WK3

function main()
    # option flags
    rstflag = "no"

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

    # data assimilation setup
    # load patient data
    filename = "pdata.mat";
    tf,qdata = patdatainterp(systems[1].mparams.th[1],filename);

    # allocators for ensemble augmented state, measurements
    X = [zeros(nvar[1]) for i in (1:length(systems))];
    θ = [zeros(9) for i in (1:length(systems))];
    Y = [zeros(1) for i in (1:length(systems))];

    # allocators for state, forecast measurement mean
    xhat = zeros(nvar[1]);
    θhat = zeros(9);
    yhat = zeros(1);

    # parameter std. deviations
    pdev = Float64[];

    # measurement std. deviations
    odev = Float64[];
    push!(odev,7);

    # allocators for measurement replicates
    yi = zeros(length(systems));

    # output variables
    Zc = Float64[];
    R = Float64[];
    C = Float64[];
    V0 = Float64[];
    t1 = Float64[];
    t2 = Float64[];
    m1 = Float64[];
    m2 = Float64[];
    Emax = Float64[];
    numassims = 0;
    println("Number of assimilations: $(length(tf))")

    # allocators for calculating normalized RMSE ratio
    r1dot = 0;
    r2dot = zeros(ensemblesize);

    # RTPS allocators
    p = 0.95; # relaxation amount, ∃ [0.5, 0.95]
    c = zeros(length(systems));
    σb = zeros(length(xhat));
    σa = zeros(length(xhat));
    σtb = zeros(length(θhat));
    σta = zeros(length(θhat));

    # solver loop
    for n in 1:numbeats
        for nn in 1:length(tf)+1
            # time integration
            if nn <= length(tf)
                tfcur = [(n-1)*0.8/WK3.ts + tf[nn]/WK3.ts for i=1:ensemblesize];
            elseif nn > length(tf)
                tfcur = [n*0.8/WK3.ts for i=1:ensemblesize];
            end
            soln = pmap((a1,a2,a3,a4,a5)->odeint(a1,a2,a3,a4,a5),y0,t0,tfcur,sparams,mparams);
            to = [soln[i][1] for i in 1:length(soln)];
            yo = [soln[i][2] for i in 1:length(soln)];
            t0 = [tfcur[1] for i in 1:length(soln)];
            for i=1:ensemblesize
                append!(systems[i].t,to[i])
            end

            # data assimilation
            if nn <= length(tf)
                println("Assimilating data at t = $(tfcur[1])")

                # dimensional measurement
                y = qdata[nn];

                # non-dimensionalize measurement
                eps = 1e-8; # scaling tolerance for small flowrates
                if abs.(yo[1][end][3]) < eps
                    Qs = eps;
                else
                    Qs = WK3.Qs;
                end
                y /= WK3.Qs;

                # generate measurement replicates
                for i = 1:length(systems)
                    yi[i] = y + rand(Distributions.Normal(0,odev[1]/Qs));
                end

                # forecast parameters and means
                θs = zeros(9);
                for i = 1:length(systems)
                    θ[i][1] = mparams[i].Zc;
                    θ[i][2] = mparams[i].R;
                    θ[i][3] = mparams[i].C;
                    θ[i][4] = mparams[i].V0;
                    θ[i][5] = mparams[i].t1;
                    θ[i][6] = mparams[i].t2;
                    θ[i][7] = mparams[i].m1;
                    θ[i][8] = mparams[i].m2;
                    θ[i][9] = mparams[i].Emax;
                    θs[1] += mparams[i].Zc;
                    θs[2] += mparams[i].R;
                    θs[3] += mparams[i].C;
                    θs[4] += mparams[i].V0;
                    θs[5] += mparams[i].t1;
                    θs[6] += mparams[i].t2;
                    θs[7] += mparams[i].m1;
                    θs[8] += mparams[i].m2;
                    θs[9] += mparams[i].Emax;
                end
                θs /= ensemblesize;
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
                            if tfcur[1]*WK3.ts == systems[1].mparams.th[end];
                                y0[k][2] = 1;
                            else
                                y0[k][2] = yo[k][j][2];
                            end
                            y0[k][3] = yo[k][j][3];
                        end
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

    # # save for post-processing
    # for i=1:ensemblesize
    #     filename = "system2_$i.mat"; # name of solution save file
    #     file = MAT.matopen(filename,"w");
    #     write(file,"system",systems[i])
    #     close(file)
    # end

    # output
    return systems
end
