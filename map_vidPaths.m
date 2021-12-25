function AvidPaths = map_vidPaths(AgrpFilenames, filenameVidPaths)
filenamesVids = match_filename(filenameVidPaths, 'Pair[A-Z]\d+');
filenameVidIDs = match_filename(filenamesVids, 'Pair[A-Z]\d+');

for fi = 1:length(AgrpFilenames)
    if ~isempty(AgrpFilenames(fi).filename1)
        afilename = AgrpFilenames(fi).filename1; 
        extracted = extract(afilename, regexpPattern('Pair[A-Z]\d+'));
        lv = strcmp(filenameVidIDs, extracted);
        if sum(lv) > 1 % if multiple "Pair" matches, match for "exp"
            extractedA = extract(afilename, regexpPattern('exp\d+'));
            biidx = contains(filenameVidPaths(lv), extractedA);
            fileidxfil = find(lv == 1);
            AvidPaths(fi).path = filenameVidPaths(fileidxfil(biidx));
        else
            AvidPaths(fi).path = filenameVidPaths(find(lv == 1));
        end
    end
end
end