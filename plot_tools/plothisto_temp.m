function plothisto_temp(temporaldistcumMat, nodeNames, numNodes, pos, XLim, YLim, figname)
figure('Position', pos) % x y w h
for ni = 1:numNodes
    %%
    subplot(1, numNodes, ni)
    histogram(temporaldistcumMat(:,ni))%, [1:10])
    ax = gca; % current axes
    ax.FontSize = 13;
    ax.YLim = YLim;
    ax.XLim = XLim;
    ylabel('instances')
    xlabel('velocity (px/fr)')
    title(nodeNames{ni})

end
saveas(gcf,figname, 'png'); close
end