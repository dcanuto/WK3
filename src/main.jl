using WK3

function main()
    # option flags
    rstflag = "no"

    # initialization
    numbeats = 30; # total number of cardiac cycles
    ensemblesize = 50;
    if rstflag == "no"
        y0 = [[rand(Distributions.Normal(10,1))*WK3.mmHgToPa/WK3.Ps;
            rand(Distributions.Normal(125,12.5))*WK3.cm3Tom3/WK3.Vs;0*WK3.cm3Tom3/WK3.Qs] for i=1:ensemblesize];
        # y0 = [[80*WK3.mmHgToPa/WK3.Ps;125*WK3.cm3Tom3/WK3.Vs;0*WK3.cm3Tom3/WK3.Qs] for i=1:ensemblesize];
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
    nparams = 9;
    X = [zeros(nvar[1]-1) for i in (1:length(systems))];
    θ = [zeros(nparams) for i in (1:length(systems))];
    Y = [zeros(1) for i in (1:length(systems))];

    # parameter scalings
    θs = zeros(nparams);
    for i = 1:length(systems)
        θs[1] += mparams[i].Zc;
        θs[2] += mparams[i].R;
        θs[3] += mparams[i].C;
        θs[4] += mparams[i].V0;
        θs[5] += mparams[i].t1;
        θs[6] += mparams[i].t2 - mparams[i].t1;
        θs[7] += mparams[i].m1;
        θs[8] += mparams[i].m2;
        θs[9] += mparams[i].Emax;
    end
    θs /= ensemblesize;

    # allocators for state, forecast measurement mean
    xhat = zeros(nvar[1]-1);
    θhat = zeros(nparams);
    yhat = zeros(1);

    # allocators for measurement replicates
    yi = zeros(length(systems));

    # output variables for ensemble avg. state, meas., params. w/ times
    tout = Float64[];
    xout = Vector{Float64}[];
    xoutv = Vector{Float64}[];
    yout = Float64[];
    innov = Float64[];
    θout = Vector{Float64}[];
    θoutv = Vector{Float64}[];
    numassims = 0;
    println("Number of assimilations: $(length(tf))")

    # distributions for state, measurement, parameters
    error = WK3.Errors(nparams);

    # allocators for calculating normalized RMSE ratio
    r1dot = 0;
    r2dot = zeros(ensemblesize);

    # RTPS allocators
    p = 0.; # relaxation amount, ∃ [0.5, 0.95]
    c = zeros(length(systems));
    σb = zeros(length(xhat));
    σa = zeros(length(xhat));
    σtb = zeros(length(θhat));
    σta = zeros(length(θhat));

    # solver loop
    tic()
    for n in 1:numbeats
        for nn in 1:length(tf)+1
            # randomize initial ventricular volume if starting new cycle
            if nn == 1 && rstflag == "no"
                for i=1:ensemblesize
                    y0[i][2] = rand(Distributions.Normal(1,0.1));
                end
            end
            # normal time integration
            if nn <= length(tf)
                tfcur = [(n-1)*0.8/WK3.ts + tf[nn]/WK3.ts for i=1:ensemblesize];
            elseif nn > length(tf)
                tfcur = [n*0.8/WK3.ts for i=1:ensemblesize];
            end
            # collect outputs
            soln = pmap((a1,a2,a3,a4,a5)->WK3.odeint(a1,a2,a3,a4,a5),y0,t0,tfcur,sparams,mparams);
            to = [soln[i][1] for i in 1:length(soln)];
            yo = [soln[i][2] for i in 1:length(soln)];
            for i=1:ensemblesize
                append!(systems[i].t,to[i])
            end

            # data assimilation
            if nn <= length(tf)
                println("Assimilating data at t = $(tfcur[1])")

                # get last time step taken, use as time step for assimilation
                if tfcur[1] != t0[1] # ignore if measurement is at t0
                    warn("Assimilating data. Fixing time step.")
                    for i in 1:length(soln)
                        sparams[i].h0 = to[i][end] - to[i][end-1];
                    end
                end

                # non-dimensional measurement
                y = qdata[nn]/WK3.Qs;

                # generate measurement replicates
                for i = 1:length(systems)
                    yi[i] = y + rand(Distributions.Normal(0,error.odev[1]/WK3.Qs));
                end

                # forecast parameters and means
                for i = 1:length(systems)
                    θ[i][1] = mparams[i].Zc;
                    θ[i][2] = mparams[i].R;
                    θ[i][3] = mparams[i].C;
                    θ[i][4] = mparams[i].V0;
                    θ[i][5] = mparams[i].t1;
                    θ[i][6] = mparams[i].t2 - mparams[i].t1;
                    θ[i][7] = mparams[i].m1;
                    θ[i][8] = mparams[i].m2;
                    θ[i][9] = mparams[i].Emax;
                end

                # forecast parameters (dimensional)
                θ = [WK3.paramwalk!(error,θs,θ[i]) for i in 1:length(systems)];

                # non-dimensionalize parameters
                for i = 1:length(systems)
                    θ[i] = θ[i]./θs;
                end

                # RTPS prior standard deviation
                for i in 1:length(θ[1])
                    for j in 1:length(systems)
                        c[j] = θ[j][i];
                    end
                    σtb[i] = std(c;corrected=true);
                end

                # parameters back into model
                for i = 1:length(soln)
                    mparams[i].Zc = θ[i][1]*θs[1];
                    mparams[i].R = θ[i][2]*θs[2];
                    mparams[i].C = θ[i][3]*θs[3];
                    # mparams[i].V0 = θ[i][4]*θs[4];
                    # mparams[i].t1 = θ[i][5]*θs[5];
                    # mparams[i].t2 = mparams[i].t1 + θ[i][6]*θs[6];
                    mparams[i].m1 = θ[i][7]*θs[7];
                    mparams[i].m2 = θ[i][8]*θs[8];
                    mparams[i].Emax = θ[i][9]*θs[9];
                end

                # forecast state w/ forecast parameters (single time step)
                soln = pmap((a1,a2,a3,a4,a5)->WK3.odeint(a1,a2,a3,a4,a5),y0,t0,tfcur,sparams,mparams);
                to = [soln[i][1] for i in 1:length(soln)];
                yoa = [soln[i][2] for i in 1:length(soln)];

                # vector of forecast state vectors
                X = [yoa[i][end][1:2] for i in 1:length(soln)]; # Pa, V only

                # vector of forecast measurements
                Y = [yoa[i][end][3] for i in 1:length(soln)]; # flowrate

                # forecast mean state, parameters, measurement
                xhat = mean(X);
                yhat = mean(Y);
                θhat = mean(θ);
                println("Normalized ̂x, first forecast: $xhat")
                println("Normalized ̂y, first forecast: $yhat")
                println("Normalized ̂θ, first forecast: $θhat")

                # forecast params./meas. cross covariance, meas. covariance
                Pty = zeros(length(θhat),length(yhat))
                Pyy = zeros(length(yhat),length(yhat))
                for i = 1:length(soln)
                    Pty += *((θ[i] - θhat),(Y[i] - yhat)');
                    Pyy += *((Y[i] - yhat),(Y[i] - yhat)');
                end
                Pty /= (ensemblesize);
                Pyy /= (ensemblesize);

                # add meas. noise to meas. covariance (allows invertibility)
                Pyy += diagm((error.odev[1]/WK3.Qs)^2*ones(length(yhat)),0);

                # parameter Kalman gain
                K = Pty*inv(Pyy);

                # parameter analysis step
                for i = 1:length(soln)
                    θ[i][:] += K*(yi[i] - Y[i]);
                end

                # RTPS parameter covariance inflation
                for i in 1:length(θ[1])
                    for j in 1:length(systems)
                        c[j] = θ[j][i];
                    end
                    σta[i] = std(c;corrected=true);
                end
                θhat = mean(θ);
                for i = 1:length(soln)
                    θ[i][:] = θ[i][:] + p*((σtb-σta)./σta).*(θ[i][:]-θhat);
                end

                # # compute FD Jacobian of observation vector w.r.t. parameters (per realization)
                # yvec = [[Y[i]] for i in 1:length(soln)];
                # J = pmap((a1,a2,a3,a4,a5,a6)->WK3.fdjac(a1,a2,a3,a4,a5,a6),y0,t0,tfcur,sparams,mparams,yvec);
                #
                # # non-dimensionalize Jacobian
                # for i in 1:length(soln) # i indexes ensemble members
                #     for j in 1:size(J[1],1) # j indexes observations
                #         J[i][j,:] .*= θs;
                #     end
                # end
                #
                # # compute dimensionless FIM (per realization)
                # I = [J[i]'*((WK3.Qs/error.odev[1])^2)*J[i] for i in 1:length(soln)]
                # println("Minimum magnitude of I[1]: $(minimum(abs.(I[1])))")
                # println("Maximum magnitude of I[1]: $(maximum(abs.(I[1])))")
                # println("Minimum rank of FIM: $(minimum(rank.(I)))")

                # analysis parameters back into ensemble members
                for i = 1:length(soln)
                    mparams[i].Zc = θ[i][1]*θs[1];
                    mparams[i].R = θ[i][2]*θs[2];
                    mparams[i].C = θ[i][3]*θs[3];
                    # mparams[i].V0 = θ[i][4]*θs[4];
                    # mparams[i].t1 = θ[i][5]*θs[5];
                    # mparams[i].t2 = mparams[i].t1 + θ[i][6]*θs[6];
                    mparams[i].m1 = θ[i][7]*θs[7];
                    mparams[i].m2 = θ[i][8]*θs[8];
                    mparams[i].Emax = θ[i][9]*θs[9];
                    systems[i].mparams = mparams[i];
                    if mparams[i].Zc <= error.lb[1]
                        println("Warning: analysis Zc below lower bound for member $i.
                            Setting to lower bound of $(error.lb[1]).")
                        mparams[i].Zc = error.lb[1];
                        θ[i][1] = error.lb[1]/θs[1];
                    end
                    if mparams[i].R <= error.lb[2]
                        println("Warning: analysis R below lower bound for member $i.
                            Setting to lower bound of $(error.lb[2]).")
                        mparams[i].R = error.lb[2];
                        θ[i][2] = error.lb[2]/θs[2];
                    end
                    if mparams[i].C <= error.lb[3]
                        println("Warning: analysis C below lower bound for member $i.
                            Setting to lower bound of $(error.lb[3]).")
                        mparams[i].C = error.lb[3];
                        θ[i][3] = error.lb[3]/θs[3];
                    end
                    if mparams[i].V0 <= error.lb[4]
                        println("Warning: analysis V0 below lower bound for member $i.
                            Setting to lower bound of $(error.lb[4]).")
                        mparams[i].V0 = error.lb[4];
                        θ[i][4] = error.lb[4]/θs[4];
                    end
                    if mparams[i].t1 <= error.lb[5]
                        println("Warning: analysis τ1 below lower bound for member $i.
                            Setting to lower bound of $(error.lb[5]).")
                        mparams[i].t1 = error.lb[5];
                        θ[i][5] = error.lb[5]/θs[5];
                    end
                    if mparams[i].t2 <= error.lb[6]
                        println("Warning: analysis τ2 below lower bound for member $i.
                            Setting to lower bound of $(error.lb[6]).")
                        mparams[i].t2 = error.lb[6];
                        θ[i][6] = error.lb[6]/θs[6];
                    end
                    if mparams[i].m1 <= error.lb[7]
                        println("Warning: analysis m1 below lower bound for member $i.
                            Setting to lower bound of $(error.lb[7]).")
                        mparams[i].m1 = error.lb[7];
                        θ[i][7] = error.lb[7]/θs[7];
                    end
                    if mparams[i].m2 <= error.lb[8]
                        println("Warning: analysis m2 below lower bound for member $i.
                            Setting to lower bound of $(error.lb[8]).")
                        mparams[i].m2 = error.lb[8];
                        θ[i][8] = error.lb[8]/θs[8];
                    end
                    if mparams[i].Emax <= error.lb[9]
                        println("Warning: analysis Emax below lower bound for member $i.
                            Setting to lower bound of $(error.lb[9]).")
                        mparams[i].Emax = error.lb[9];
                        θ[i][9] = error.lb[9]/θs[9];
                    end
                end

                # recalculate and output parameter averages (non-dimensional)
                θhat = mean(θ);
                append!(θout,[θhat])
                println("Normalized analysis parameter averages: $θhat")

                # corrected forecast w/ analysis parameters
                soln = pmap((a1,a2,a3,a4,a5)->WK3.odeint(a1,a2,a3,a4,a5),y0,t0,tfcur,sparams,mparams);
                to = [soln[i][1] for i in 1:length(soln)];
                yoa = [soln[i][2] for i in 1:length(soln)];

                # vector of forecast state vectors
                X = [yoa[i][end][1:2] for i in 1:length(soln)]; # Pa, V only

                # vector of forecast measurements
                Y = [yoa[i][end][3] for i in 1:length(soln)]; # flowrate

                # prior standard deviation
                for i in 1:length(X[1])
                    for j in 1:length(systems)
                        c[j] = X[j][i];
                    end
                    σb[i] = std(c;corrected=true);
                end

                # second forecast mean state, measurement
                xhat = mean(X);
                yhat = mean(Y);
                println("Normalized ̂x, second forecast: $xhat")
                println("Normalized ̂y, second forecast: $yhat")

                # output second forecast measurement
                push!(yout,yhat)
                push!(innov,(yhat-y))

                # forecast state/meas. cross covariance, meas. covariance
                Pxy = zeros(length(xhat),length(yhat))
                Pyy = zeros(length(yhat),length(yhat))
                for i = 1:length(soln)
                    Pxy += *((X[i] - xhat),(Y[i] - yhat)');
                    Pyy += *((Y[i] - yhat),(Y[i] - yhat)');
                end
                Pxy /= (ensemblesize);
                Pyy /= (ensemblesize);

                # add meas. noise to meas. covariance (allows invertibility)
                Pyy += diagm((error.odev[1]/WK3.Qs)^2*ones(length(yhat)),0);

                # Kalman gain
                K = Pxy*inv(Pyy);

                # analysis step, NRR tracking
                for i = 1:length(soln)
                    X[i][:] += K*(yi[i] - Y[i]);
                    r2dot[i] += dot((Y[i]-yi),(Y[i]-yi));
                end
                r1dot += sqrt.(dot((yhat-y),(yhat-y)));

                # RTPS multiplicative covariance inflation
                for i in 1:length(X[1])
                    for j in 1:length(systems)
                        c[j] = X[j][i];
                    end
                    σa[i] = std(c;corrected=true);
                end
                xhat = mean(X);
                # println("Prior standard deviation: $(σb)")
                # println("Posterior standard deviation: $(σa)")
                println("Normalized difference: $((σb-σa)./σa)")
                for i = 1:length(soln)
                    X[i][:] = X[i][:] + p*((σb-σa)./σa).*(X[i][:]-xhat);
                end

                # output ensemble average state
                xhat = mean(X);
                append!(xout,[xhat])

                # analysis back into ensemble members
                # reshape output to vectors of individual state variables' time series
                for k = 1:ensemblesize # k indexes ensemble members
                    # yo[k][end][1:2] = X[k][:]; # replace last step Pa, V w/ analysis
                    # yo[k][end][end] = yoa[k][end][end];
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
                                # reset initial volume if starting ending cardiac cycle
                                if tfcur[1]*WK3.ts == systems[1].mparams.th[end] && nn == length(tf);
                                    y0[k][2] = 1;
                                else
                                    y0[k][2] = yo[k][j][2];
                                end
                                y0[k][3] = yo[k][j][3];
                            end
                        end
                    end
                end

                # reset initial step size
                for i in 1:length(soln)
                    sparams[i].h0 = sparams[i].h0nom;
                end
            elseif nn == length(tf)+1 # make sure we get last set of outputs w/out measurement data
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

            # setup initial time for normal forward time integration
            t0 = [tfcur[1] for i in 1:length(soln)];

            # output measurement time
            if nn <= length(tf)
                push!(tout,t0[1])
            end
        end
    end
    toc()

    # reshape ensemble averages into vectors of individual time series
    for i in 1:length(xout[1]) # i indexes state variables
        xo = Float64[];
        for j in 1:length(xout) # j indexes time steps
            push!(xo,xout[j][i])
        end
        append!(xoutv,[xo])
    end

    for i in 1:length(θout[1]) # i indexes state variables
        θo = Float64[];
        for j in 1:length(θout) # j indexes time steps
            push!(θo,θout[j][i])
        end
        append!(θoutv,[θo])
    end

    # normalized RMSE ratio (optimal ensemble yields NRR ~ 1)
    r2dot = sqrt.(1/systems[1].t[end]*r2dot);
    R1 = 1/systems[1].t[end]*r1dot;
    R2 = sum(r2dot);
    Ra = R1/R2;
    ERa = sqrt.((ensemblesize)/(2*(ensemblesize)));
    println("Normalized RMSE ratio: $(Ra/ERa)")

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
    #     filename = "system_$i.mat"; # name of solution save file
    #     file = MAT.matopen(filename,"w");
    #     write(file,"system",systems[i])
    #     close(file)
    # end

    # output
    return systems,tout,xoutv,yout,innov,θoutv
end
