function draw_lines_xy(a, b, frvec, axisdims)
% draws animal a, b at same time

figure
aline = animatedline('Color','b');
bline = animatedline('Color','r');
axis(axisdims)

%%
for k = 1:length(frvec)
    addpoints(aline,a(frvec(k),1), a(frvec(k),2)); % x coordinate, y coordinate
    addpoints(bline,b(frvec(k),1), b(frvec(k),2));
    drawnow
end

end