
%%
for id = 1%:length(dvidfolders)
%     folderpath = [dvidfolders(id).folder filesep dvidfolders(id).name];
    folderpath = pwd;
    cd(folderpath)
    dfolder = dir(folderpath);
    mp4idx = find(endsWith({dfolder.name},'.mp4')==1);
    for ii = 1:length(mp4idx)
        vidfullpath = [folderpath filesep dfolder(mp4idx(ii)).name];
        
        %% note approximately 10-15 min per 10 min video 30 fps
        reader = VideoReader(vidfullpath);
        writer = VideoWriter([folderpath filesep dfolder(mp4idx(ii)).name '.avi'], 'Motion JPEG AVI');%'Uncompressed AVI');

        writer.FrameRate = reader.FrameRate;
        open(writer);

        for j = 1:reader.NumFrames
           img = readFrame(reader);
           writeVideo(writer,img);
        end
        
        close(writer);

    end
    

end
