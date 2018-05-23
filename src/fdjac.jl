function fdjac(y0::Vector{Float64},t0::Float64,tf::Float64,sparams::WK3.SolverParams,
    mparams::WK3.ModelParams,obs::Vector{Float64})
    sqrteps = sqrt(eps()); # ≈ square root of machine precision
    nparams = 9;
    nobs = 1;
    J = zeros(nobs,nparams);
    for j = 1:nparams
        if j == 1
            temp = mparams.Zc;
        elseif j == 2
            temp = mparams.R;
        elseif j == 3
            temp = mparams.C;
        elseif j == 4
            temp = mparams.V0;
        elseif j == 5
            temp = mparams.t1;
        elseif j == 6
            temp = mparams.t2;
        elseif j == 7
            temp = mparams.m1;
        elseif j == 8
            temp = mparams.m2;
        elseif j == 9
            temp = mparams.Emax;
        end
        h = sqrteps*abs.(temp); # scale step size by current value
        if h == 0
            h = sqrteps;
        end
        if j == 1
            mparams.Zc = temp+h; # trick to reduce finite precision error
            h = mparams.Zc - temp;
        elseif j == 2
            mparams.R = temp+h;
            h = mparams.R - temp;
        elseif j == 3
            mparams.C = temp+h;
            h = mparams.C - temp;
        elseif j == 4
            mparams.V0 = temp+h;
            h = mparams.V0 - temp;
        elseif j == 5
            mparams.t1 = temp+h;
            h = mparams.t1 - temp;
        elseif j == 6
            mparams.t2 = temp+h;
            h = mparams.t2 - temp;
        elseif j == 7
            mparams.m1 = temp+h;
            h = mparams.m1 - temp;
        elseif j == 8
            mparams.m2 = temp+h;
            h = mparams.m2 - temp;
        elseif j == 9
            mparams.Emax = temp+h;
            h = mparams.Emax - temp;
        end
        ~,yout = odeint(y0,t0,tf,sparams,mparams);
        if j == 1
            mparams.Zc = temp;
        elseif j == 2
            mparams.R = temp;
        elseif j == 3
            mparams.C = temp;
        elseif j == 4
            mparams.V0 = temp;
        elseif j == 5
            mparams.t1 = temp;
        elseif j == 6
            mparams.t2 = temp;
        elseif j == 7
            mparams.m1 = temp;
        elseif j == 8
            mparams.m2 = temp;
        elseif j == 9
            mparams.Emax = temp;
        end
        obsp = [yout[end][3]]; # dimensionless flowrate from perturbed parameter
        for i = 1:nobs # forward difference formula
            J[i,j] = (obsp[i]-obs[i])/h;
        end
    end
    return J
end
