function patdatainterp(th::Float64,filename)
    # load patient data
    file = MAT.matread(filename);
    Qdata = file["Q"]; # Q should be in m3/s
    tdata = file["tgood"]; # times should be in s
    Pdata = file["P"]; #  pressure should be in Pa
    Vdata = file["V"]; # volume should be in m3
    thdata = file["th"];
    # vectorize data
    Qdata = vec(Qdata);
    tdata = vec(tdata);
    Pdata = vec(Pdata);
    Vdata = vec(Vdata);
    # create interpolation objects
    Qditp = Interpolations.interpolate(Qdata, Interpolations.BSpline(Interpolations.Linear()), Interpolations.OnGrid());
    Pditp = Interpolations.interpolate(Pdata, Interpolations.BSpline(Interpolations.Linear()), Interpolations.OnGrid());
    Vditp = Interpolations.interpolate(Vdata, Interpolations.BSpline(Interpolations.Linear()), Interpolations.OnGrid());
    # perform interpolations
    t = linspace(0,th,length(Qdata));
    tq = tdata*th/thdata;
    sitp = Interpolations.scale(Qditp,t);
    qq = [sitp[j] for j in tq];
    sitp = Interpolations.scale(Pditp,t);
    pq = [sitp[j] for j in tq];
    sitp = Interpolations.scale(Vditp,t);
    vq = [sitp[j] for j in tq];
    return tq,qq,pq,vq
end
