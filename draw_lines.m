function draw_lines(a, b, frvec, axisdims)

figure
aline = animatedline('Color','b');
bline = animatedline('Color','r');
axis(axisdims)

%%
for k = 1:length(frvec)
    addpoints(aline,a(frvec(k),1,1), -a(frvec(k),1,2)); % x coordinate, negative y coordinate
    addpoints(bline,b(frvec(k),1,1), -b(frvec(k),1,2));
    drawnow
end

end