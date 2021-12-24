function percentdropAll = calc_perc(dropctsAll, totalfrAll)
% plot percentage, by each animal
percentdropAll = [];
for ii = 1:size(dropctsAll, 1) % animal
    counts = cellfun(@sum, dropctsAll(ii,:));
    percent = counts ./ totalfrAll(ii);
    percentdropAll = vertcat(percentdropAll, percent);
end
percentdropAll = percentdropAll .* 100;
end