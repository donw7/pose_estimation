function plotperc_dropct(percentdropAll, X, leg, pos, YLim, figname)
%%
figure('Position', pos) % x y w h
hold on
cmap = colormap(jet(size(percentdropAll,1)));
for ai = 1:size(percentdropAll, 1)
    p = plot(X, percentdropAll(ai,:), '-o', 'Color', cmap(ai,:))%, 'Color', 'k'); %, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
end
%%
legend(leg)
title('percent dropped frames per animal')
ylabel('percent')
ax = gca; ax.FontSize = 13; ax.YLim = YLim;
saveas(gcf,figname, 'png'); close
end

