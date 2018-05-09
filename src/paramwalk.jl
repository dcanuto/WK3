function paramwalk!(err::WK3.Errors,pbar::Vector{Float64},θ::Vector{Float64})
    for i = 1:length(θ)
        θ[i] = rand(Distributions.TruncatedNormal(err.a*θ[i]+(1-err.a)*pbar[i],err.h*
            err.pdev[i],err.lb[i],err.ub[i]));
    end
    return θ
end
