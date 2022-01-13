% calc & plot temporal distances (velocity)
function [temporaldistcumMat, Aslpvel] = calc_tempdist(Aslpxypairs, numNodes, Aslpfridx)
temporaldistcumMat = [];
Aslpvel = struct();
fnames = fieldnames(Aslpxypairs);

for grpi = 1:length(Aslpxypairs) % group
    for ai = 1:length(fnames) % animal
        if ~isempty(Aslpxypairs(grpi).(fnames{ai})) % assumes if first empty, rest all empty
            
            tracklogvec = Aslpfridx(grpi).('idx');
            trackIdx = find(tracklogvec == 1);
                
            Distances = [];
            for npii = 1:numNodes
                %%
                for fri = 2:length(trackIdx) % frame by frame
                    %%
                    nodeXYmat = [Aslpxypairs(grpi).(fnames{ai}){npii}(trackIdx(fri) - 1,:);...
                        Aslpxypairs(grpi).(fnames{ai}){npii}(trackIdx(fri),:)];
                    
                    Distance = pdist(nodeXYmat);
                    Distances(fri, npii) = Distance;
                    
                end
            end
            temporaldistcumMat = vertcat(temporaldistcumMat, Distances);
            Aslpvel(grpi).(fnames{ai}) = Distances;
        end
    end
end

end