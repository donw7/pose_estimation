function Aout = wrap_eaAni(Aall, fhandle, params)
% wrap per animal operation
Aout = struct();
Snames = fieldnames(Aall);
for si = 1:length(Snames) 
    C = Aall.(Snames{si});
    eifields = fieldnames(C);
    for ei = 1:length(eifields) % eg exp, exp2
        Exp = Aall.(Snames{si}).(eifields{ei});
        aifields = fieldnames(Exp);
        for gi = 1:length(Exp) 
            for ai = 1:length(aifields)
                if ~isempty(Exp(gi).(aifields{ai}))
                    % eg ..exp(1).track1
                    Aout.(Snames{si}).(eifields{ei})(gi).(aifields{ai}) =...
                        fhandle(Exp(gi).(aifields{ai}), params);            
                end
            end
        end
    end
end
end