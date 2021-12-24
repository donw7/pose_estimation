function plothisto_spat(spatialdistcumMat, nodeNames, nodeIdxFront, pos, XLim, YLim, figname)
figure('Position', pos) % x y w h
for ni = 1:length(nodeIdxFront)
    %%
    subplot(1, length(nodeIdxFront), ni)
    histogram(spatialdistcumMat(:,ni), 'BinWidth', 25)%, [1:10])
    ax = gca; % current axes
    ax.FontSize = 13;
    ax.XLim = XLim;
    ax.YLim = YLim;
    ylabel('instances')
    xlabel('distance (px)')
    Titlenodepair = strcat(cellstr(nodeNames{nodeIdxFront(1,ni)}), {'-'}, cellstr(nodeNames{nodeIdxFront(2,ni)}));
    title(Titlenodepair)

end
saveas(gcf,figname, 'png'); close
end