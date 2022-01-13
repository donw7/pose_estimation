function Aslpfridx = getfridx(AslpsocDist, distpx)
Aslpfridx = struct();
for grpi = 1:length(AslpsocDist) % group
    if ~isempty(AslpsocDist(grpi).('socDistance'))
        trackIdx = AslpsocDist(grpi).('socDistance') < distpx;
        Aslpfridx(grpi).('idx') = trackIdx;
    end
end
end