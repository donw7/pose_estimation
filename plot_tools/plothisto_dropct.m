% plot histograms of dropped fr count instances
function plothisto_dropct(dropctsAll, numNodes, nodeNames, pos, YLim, figname)
figure('Position', pos) % x y w h
dropctsMat = cell(1, numNodes);
for ci = 1:numNodes
    all_concat = [dropctsAll{:,ci}];
    dropctsMat{ci} = all_concat;
    subplot(1, numNodes, ci)
    histogram(all_concat, [1:10])
    ax = gca; % current axes
    ax.FontSize = 13;
    ax.YLim = YLim;
    ylabel('instances of dropped frames events')
    xlabel('duration (# fr)')
    title(nodeNames{ci})
end
saveas(gcf,figname, 'png'); close
end