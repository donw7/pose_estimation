% calc & plot histograms of dropped fr counts
function [Aslpspat, spatialdistcumMat] = calc_pwdist(Aslpxypairs, params)
% pass eg nodeIdxbody eg [1 2; 7 7] to nodeIdxcomp
spatialdistcumMat = [];
Aslpspat = struct();
fields = fieldnames(Aslpxypairs);

parfor grpi = 1:length(Aslpxypairs) % group
    for ai = 1:length(fields) % animal
        if ~isempty(Aslpxypairs(grpi).(fields{ai})) % assumes if first empty, rest all empty
            
            Distances = [];
            for npii = 1:size(params.nodeIdxcomp, 2) % eg 1-4, 2-4, etc.
                
                if isfield(params, 'Aslpfridx')
                    tracklogvec = params.Aslpfridx(grpi).('idx');
                    trackIdx = find(tracklogvec == 1);
                else
                    trackIdx = 1:length(Aslpxypairs(grpi).(fields{ai}){1});
                end
            
                for fri = 1:length(trackIdx) % frame by frame
                    % [x-y nodeA; x-y nodeB]
                    nodeAidx = params.nodeIdxcomp(1,npii); % single values to find node column
                    nodeBidx = params.nodeIdxcomp(2,npii); % within same animal
                    
                    nodeXYmat = [Aslpxypairs(grpi).(fields{ai}){nodeAidx}(trackIdx(fri),:);...
                        Aslpxypairs(grpi).(fields{ai}){nodeBidx}(trackIdx(fri),:)];

                    % distance and store
                    Distance = pdist(nodeXYmat);
                    Distances(fri, npii) = Distance;
                end
            end
            spatialdistcumMat = vertcat(spatialdistcumMat, Distances);
            Aslpspat(grpi).(fields{ai}) = Distances;
        end
    end
end

end