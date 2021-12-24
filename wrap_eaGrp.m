function Aout = wrap_eaGrp(Aall, fhandle, params)
% wrap per group operation (eg cmk.exp, not group at eg pair level)
% should work on E struc format as well as wrapped nonscalars
Aout = struct();
Snames = fieldnames(Aall);
for si = 1:length(Snames) % eg cmk, dlx
    C = Aall.(Snames{si});
    eifields = fieldnames(C);
    for ei = 1:length(eifields) % eg exp, exp2
        Exp = Aall.(Snames{si}).(eifields{ei});
        Aout.(Snames{si}).(eifields{ei}) =...
            fhandle(Exp, params);            
    end
end
end