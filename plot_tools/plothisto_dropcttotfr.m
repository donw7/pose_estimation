% plot histograms of dropped fr counts, total # of fr
function plothisto_dropcttotfr(dropctsAll, numNodes, nodeNames, pos, XLim, YLim, figname)
figure('Position', pos) % x y w h
%%
dropctsMat = cell(1, numNodes);
for ci = 1:numNodes
    %%
    all_concat = [dropctsAll{:,ci}];
    dropctsMat{ci} = all_concat;
    unis = unique(all_concat);
    x = 1:maxk(unis,1); y = zeros(size(x));
    for ui = 1:length(unis)
        y(unis(ui)) = sum(all_concat == unis(ui)) * unis(ui);
    end
    subplot(1, numNodes, ci)
    plot(x, y, 'LineWidth', 1)
    ax = gca; % current axes
    ax.FontSize = 13;
    ax.XLim = XLim; ax.YLim = YLim;
    ylabel('total number of dropped frames')
    xlabel('duration (# fr)')
    title(nodeNames{ci})
end
saveas(gcf,figname, 'png'); close
end