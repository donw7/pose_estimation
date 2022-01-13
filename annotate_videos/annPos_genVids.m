function annPos_genVids(A, params)
%% 
%% eg
% params.Draw.anncircradius = 3;
% params.Draw.inscolors = {'blue', [204, 102, 0], 'magenta', 'green', 'cyan', 'black'}; % up to 5 animal colors
Axy = A.AslpxyInti;
Avec = A.Annveci;

Avpaths = A.AvidPathsi;
ininscolors = params.Draw.inscolors;

if isfield(A, 'AnnOutFr')
    AFr = A.AnnOutFr; % cell array of frames if not writing all
else
    AFr = 'nothing';
end
        
%% save separate videos for each 
for gi = 1:length(Axy)
    fields = fieldnames(Axy(gi));

    if ~isempty(Axy(gi).(fields{1}))
        %%
        S = {}; % eg xy
        V = {}; % eg vec to label
        L = {}; % labels
        
        for ai = 1:length(fields)
            S{ai} = Axy(gi).(fields{ai});
            V{ai} = Avec(gi).(fields{ai}).vec;
            L{ai} = Avec(gi).(fields{ai}).labels;
        end
        vidpath = Avpaths(gi).path{1};
        outVid = VideoWriter([vidpath '_ann']); % prep video container
        tic
            % w/ imshow: 1.5 h for 35k frames
            % w/ im2frame: approx 46 sec for 1000 fr | approx 0.5 h for 35k frames? 
        
        vidread = VideoReader(vidpath);
        
        if isfield(A, 'AnnOutFr')
            fridx = AFr{gi};
        else
            fridx = 1:vidread.NumFrames;
        end
        
        for fri = 1:length(fridx) % each frame
            %%
            I = read(vidread, fridx(fri));
            RGBann = I;
            for ai = 1:length(fields) % draw for each animal
                %%
                pos = zeros(1,2); % for insertMarker [x y] only;
                for ni = 1:length(S{ai})
                    x = S{ai}{ni}(fridx(fri),1); y = S{ai}{ni}(fridx(fri),2);
                    pos(ni,:) = [x y];
                end
                RGBann = insertMarker(RGBann, pos, 'o', 'Color', ininscolors{ai});%,'color',color,'size',10);

                conn = params.Draw.connections;
                for ci = 1:size(conn, 2)
                    RGBann = insertShape(RGBann, 'Line', [pos(conn(1,ci),1) pos(conn(1,ci),2)... % for insertShape line [x1 y1 x2 y2] only
                        pos(conn(2,ci),1) pos(conn(2,ci),2)], 'Color', ininscolors{ai}); 
                end
                % for insertObjectAnnotation [x y radius] for circle; [x y w h] for rectangle
                posxyr = horzcat(pos, repmat(params.Draw.anncircradius, [size(pos, 1), 1]));
                for ni = 1:size(V{ai}, 2)
                    if V{ai}(fridx(fri), ni) == 1 % fr x nodes
    %                 x = S{ai}{ni}(fridx(fri),1); y = S{ai}{ni}(fridx(fri),2);
    %                 pos(ni,:) = [x y];
                        posxyri = posxyr(ni,:);
                        label = strcat(L{ai}{fridx(fri), ni}, '-', params.nodeNames{ni}); 
                        RGBann = insertObjectAnnotation(RGBann, 'circle',...
                            posxyri, label, 'LineWidth', 1, 'FontSize', 8, 'Color', ininscolors{ai});
                    end
                end

            end
            open(outVid)
            I_print = im2frame(RGBann); % faster than getframe(gca)
            I_print = I_print.cdata;
            writeVideo(outVid,I_print);

        end
        fprintf('Printed video for Pair %d\n', gi)
        toc

        close(outVid);
        close all
    end

end


    
end