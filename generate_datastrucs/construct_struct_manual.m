function E = construct_struct_manual(numTrials, numAnimals)
%% nest folders into cell array
E = cell(1, numTrials); % for output

% for output
for i = 1:numTrials
    for ai = 1:numAnimals
        E{i}{ai} = {};
    end
end

% D3 = K(mod(K,N)==0)
% D2 = K(mod(K,N)==0) - 1
% D1 = K(mod(K,N)==0) - 2
end