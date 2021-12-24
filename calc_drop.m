function [Dropslp, dropctsAll, totalfrAll, totalfrcum] = calc_drop(Aslp, params)
Dropslp = struct();
dropctsAll = {};
totalfrAll = [];
totalfrcum = 0;
fnames = fieldnames(Aslp);
%%
for grpi = 1:length(Aslp) % group
    for ai = 1:length(fnames) % animal
        if ~isempty(Aslp(grpi).(fnames{ai})) % assumes if first empty, rest all empty
            
            %% index out nonproximal frames
            trackmat = Aslp(grpi).(fnames{ai});
            
            if isfield(params, 'Aslpfridx')
                trackIdx = params.Aslpfridx(grpi).('idx');
            else
                trackIdx = 1:size(trackmat, 1);
            end

            Dropslp(grpi).(fnames{ai}) = count_consecutiveNans(trackmat); 
            dropcts = Dropslp(grpi).(fnames{ai}).dropcts;
            dropctsAll = vertcat(dropctsAll, dropcts);
            totalfrcum = totalfrcum + size(trackmat, 1);
            totalfrAll = vertcat(totalfrAll, size(trackmat, 1));
            
        end
    end
end
end