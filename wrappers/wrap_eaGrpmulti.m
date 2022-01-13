function Aout = wrap_eaGrpmulti(Aallmulti, fhandle, params)
% same except another layer, eg for multiple Aslp-level inputs eg Dropcts
% and AslpatB
Aout = struct();

Snames = fieldnames(Aallmulti);
for si = 1:length(Snames) % animal type
    C = Aallmulti.(Snames{si});
    eifields = fieldnames(C);
    for ei = 1:length(eifields) % eg exp, exp2
        E = C.(eifields{ei}); % has multi structs
        Aout.(Snames{si}).(eifields{ei}) = fhandle(E, params);            
    end
end
end