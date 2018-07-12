function patdatainterp(th::Float64,filename)
    file = MAT.matread(filename);
    Qdata = file["Q"]; # Q should be in m3/s
    tdata = file["tgood"]; # times should be in s
    thdata = file["th"];
    Qdata = vec(Qdata);
    tdata = vec(tdata);
    Qditp = Interpolations.interpolate(Qdata, Interpolations.BSpline(Interpolations.Linear()), Interpolations.OnGrid());
    t = linspace(0,th,length(Qdata));
    sitp = Interpolations.scale(Qditp,t);
    tq = tdata*th/thdata;
    qq = [sitp[j] for j in tq];
    return tq,qq
end
