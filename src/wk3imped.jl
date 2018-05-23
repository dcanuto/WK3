function wk3imped(ω::Float64,Rc::Vector{Float64},Rp::Vector{Float64},C::Vector{Float64})
    Zt = Complex{Float64}[];
    magZt = Complex{Float64}[];
    phaseZt = Complex{Float64}[];
    append!(Zt,Rc + Rp./((ω*C.*Rp).^2+1) + im*(1./(ω*C))./(1+1./((ω*Rp.*C).^2)));
    append!(magZt,abs.(Zt))
    append!(phaseZt,angle.(Zt))
    return Zt,magZt,phaseZt
end
