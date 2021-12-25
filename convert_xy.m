function Aslpxypairs = convert_xy(Aslp, p)
% converts 3-D mat to x-y pairs
% option for interpolation or not
Aslpxypairs = struct();

for grpi = 1:length(Aslp) % group
    fnames = fieldnames(Aslp);
    for ai = 1:length(fnames) % animal
        if ~isempty(Aslp(grpi).(fnames{ai})) % assumes if first empty, rest all empty
            %% b/w all pairwise relationships of all nodes
            trackmatX = Aslp(grpi).(fnames{ai})(:,:,1);
            trackmatY = Aslp(grpi).(fnames{ai})(:,:,2);

            %% make x-y pairs of each point element w/ vs w/o interpolate
            nodeXY = {}; % cell of length of # nodes
            for ni = 1:length(p.nodeNames)
                nodeXYmat = [];
                switch p.interpTF
                    case 'True'
                        trackmatX(:,ni) = fillmissing(trackmatX(:,ni), 'movmean', 10);
                        trackmatY(:,ni) = fillmissing(trackmatY(:,ni), 'movmean', 10);
                        
                        nodeXYmat = [fillmissing(trackmatX(:,ni), 'linear'),...
                            fillmissing(trackmatY(:,ni), 'linear')]; % x, y coordinates
                        
                    case 'False'
                        nodeXYmat = [trackmatX(:,ni), trackmatY(:,ni)];
                end
                nodeXY{ni} = nodeXYmat;
            end
            Aslpxypairs(grpi).(fnames{ai}) = nodeXY;

        end
    end
end


end