function rkck(y::Vector{Float64},dy::Vector{Float64},n::Int64,t::Float64,
    h::Float64,yout::Vector{Float64},yerr::Vector{Float64},derivs::Function)
    # function evaluation weights
    a2=0.2;a3=0.3;a4=0.6;a5=1.0;a6=0.875;b21=0.2;b31=3.0/40.0;b32=9.0/40.0;
    b41=0.3;b42=-0.9;b43=1.2;b51=-11.0/54.0;b52=2.5;b53=-70.0/27.0;b54=35.0/27.0;
    b61=1631/55296;b62=175/512;b63=575/13824;b64=44275/110592;b65=253/4096;c1=37/378;
    c3=250/621;c4=125/594;c6=512/1771;dc5=-277/14336;dc1=c1-2825/27648;dc3=c3-18575/48384;
    dc4=c4-13525/55296;dc6=c6-0.25;

    # allocators for intermediate derivatives and solution
    ak2=zeros(n);
    ak3=zeros(n);
    ak4=zeros(n);
    ak5=zeros(n);
    ak6=zeros(n);
    ytemp=zeros(n);

    # first step
    ytemp = y + h*b21*dy;
    # second step
    derivs(t+a2*h,ytemp,ak2);
    ytemp = y + h*(b31*dy+b32*ak2);
    # third step
    derivs(t+a3*h,ytemp,ak3);
    ytemp = y + h*(b41*dy+b42*ak2+b43*ak3);
    # fourth step
    derivs(t+a4*h,ytemp,ak4);
    ytemp = y + h*(b51*dy+b52*ak2+b53*ak3+b54*ak4);
    # fifth step
    derivs(t+a5*h,ytemp,ak5);
    ytemp = y + h*(b61*dy+b62*ak2+b63*ak3+b64*ak4+b65*ak5);
    # sixth step
    derivs(t+a6*h,ytemp,ak6);
    # accumulate increments w/ proper weights
    yout[:] = y + h*(c1*dy+c3*ak3+c4*ak4+c6*ak6);
    # estimate error
    yerr[:] = h*(dc1*dy+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6);

end
