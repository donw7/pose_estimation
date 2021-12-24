function [outStruct] = count_consecutiveNans(inMat)
% finds Nans of input matrix and outputs struct w/ properties

outStruct.Nans = isnan(inMat);
outStruct.startIdx = cell(1, size(inMat, 2));
outStruct.dropcts = cell(1, size(inMat, 2));

%%
for coli = 1:size(inMat, 2) % eg each node
    %%
    x = outStruct.Nans(:,coli)';
    f = find(diff([0,x,0]==1));
    outStruct.startIdx{coli} = f(1:2:end-1);  % Start indices
    outStruct.dropcts{coli} = f(2:2:end)-outStruct.startIdx{coli};  % Consecutive ones’ counts

    
end

end