function Aslpfridxspat = getfridxspat(AslpspatB, p)
Aslpfridxspat = struct();
for grpi = 1:length(AslpspatB) % group
    fnames = fieldnames(AslpspatB);
    for ai = 1:length(fnames)
        if ~isempty(AslpspatB(grpi).(fnames{ai}))
            switch p.compTF
                case 'less'
                    trackIdx = AslpspatB(grpi).(fnames{ai}) < p.threshdistpx;
                case 'greater'
                    trackIdx = AslpspatB(grpi).(fnames{ai}) > p.threshdistpx;
            end
            Aslpfridxspat(grpi).(fnames{ai}) = trackIdx;
        end
    end
end
end