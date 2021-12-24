%% pipeline for validation and sourcing of key training examples, w/ spatial statistics and video annotations for key mistakes

clearvars %-except

%% assume in working directory
% if need to build paths, see multidirectory_loader
sleapPath = [pwd filesep 'drosophila'];
dh5files = dir('**/*.analysis.h5'); % dir struc obj after search of current and all subfolders
filenamesAnalysis = {dh5files.name}';

% source video paths
vidrootpath = [pwd filesep 'drosophila'];

%% set up options for analysis
clear params
params.numAnimals = 2;
params.numTrials = 1;
params.numNodes = 7;
params.nodeidxHead = 1:4; % specify 4 pts
params.nodeidxBody = 7;
params.nodeidxBodyTargets = 1:6;
params.nodeidxSB = 4;
params.nodeidxSBTargets = [1:3, 5, 7];
disp(params)

% initialize dataStruct
Slp = struct();
Slp = construct_struct_manual(params.numTrials, params.numAnimals);


%% load data by available filenames
for fi = 1:length(filenamesAnalysis) % iterate by available filenames
    aniS = struct();
    aniS.filename = filenamesAnalysis{fi};       
    
    fipath = [dh5files(fi).folder filesep dh5files(fi).name];
    aniS.filepath = fipath;
    aniS.occupancy_matrix = h5read(fipath,'/track_occupancy');
    aniS.tracks_matrix = h5read(fipath,'/tracks');
    aniS.track_names = h5read(fipath,'/track_names');
    aniS.node_names = h5read(fipath,'/node_names');
    for ni = 1:length(aniS.node_names) deblank(aniS.node_names{ni}); end

%     aniS.point_scores = h5read(fipath,'/point_scores');
%     aniS.instance_scores = h5read(fipath,'/instance_scores');
%     aniS.tracking_scores = h5read(fipath,'/tracking_scores');

    % tailmarked is animal 1 ('track_0') usually
    aniS.track = aniS.tracks_matrix(:,:,:,1); % specify tracks (e.g. if A, B is in 1:2 or 2:3 
    Slp{fi}{1} = aniS;

    % animal 2 'track_1'
    aniS.track = aniS.tracks_matrix(:,:,:,2); % specify tracks (e.g. if A, B is in 1:2 or 2:3 
    Slp{fi}{2} = aniS;
% 
%     % consider if need to account for missed tracks - hardcoded here 

end

%%
params.nodeNames = Slp{fi}{1}.node_names;

%% check visually
% aniStruct.track is frames x nodes x 2 (x and y)?
draw_lines(Slp{1}{1}.track, Slp{1}{2}.track, 1:500, [0 1000 -1000 0])

% validation, save figs



%% extract-initialize
Aslp = extract_classification(Slp, {'track', 'track'}); %* manual set analyze by cell type
%
fnames = fieldnames(Aslp);

%% convert to x-y pairs
Aslpxypairs = struct();
for grpi = 1:length(Aslp) % group
    for ai = 1:length(fnames) % animal
        if ~isempty(Aslp(grpi).(fnames{ai})) % assumes if first empty, rest all empty
            %% b/w all pairwise relationships of all nodes
            trackmatX = Aslp(grpi).(fnames{ai})(:,:,1);
            trackmatY = Aslp(grpi).(fnames{ai})(:,:,2);

            %% make x-y pairs of each point element & *interpolate
            nodeXY = {}; % cell of length of # nodes
            sqforms = []; % node x node x fr
            for ni = 1:length(params.nodeNames)
                nodeXYmat = [];
                nodeXYmat = [repnan(trackmatX(:,ni)), repnan(trackmatY(:,ni))]; % x, y coordinates
                nodeXY{ni} = nodeXYmat;
            end
            Aslpxypairs(grpi).(fnames{ai}) = nodeXY;
            
        end
    end
end

% proximity based sampling
AslpsocDist = struct();
nodeidx = 1:4;
socDistCumAll = [];

for grpi = 1:length(Aslp) % group
    if ~isempty(Aslp(grpi).(fnames{ai})) % assumes if first empty, rest all empty
        numFr = size(Aslpxypairs(grpi).(fnames{ai}){1}, 1); % take fr from first node

        % [x-y median of front nodes in animal A; x-y of animal B]
        aniApos = horzcat(Aslpxypairs(grpi).(fnames{1}){nodeidx});
        aniBpos = horzcat(Aslpxypairs(grpi).(fnames{2}){nodeidx});

        %%
        DistanceAnis = [];
        for fri = 1:numFr % frame by frame
            aniAavgX = nanmean(aniApos(fri, [1 3 5 7])); % if 4 points
            aniBavgY = nanmean(aniBpos(fri, [2 4 6 8]));

            nodeXYmat = [aniAavgX; aniBavgY];

            % find distance and store
            Distance = pdist(nodeXYmat);
            DistanceAnis(fri, 1) = Distance;

        end
        AslpsocDist(grpi).('socDistance') = DistanceAnis;
    end
end
socDistCumAll = vertcat(AslpsocDist.('socDistance'));

%% make index that encapsulates all
params.distpx = 3000; % set high to include all fr 
params.Aslpfridx = getfridx(AslpsocDist, params.distpx); %(used to limit frames later)

%% calc/plot histo, percent dropped
% specify position and YLim
[Dropslp, dropctsAll, totalfrAll, totalfrcum] = calc_drop(Aslp, params);

%%
plothisto_dropct(dropctsAll, params.numNodes, params.nodeNames, [100 400 1500 600], [0 100], 'fig1')
X = reordercats(categorical(params.nodeNames), params.nodeNames); % fixed reordering of categories for display
leg = repelem(params.numTrials, params.numTrials); % repeated pair ID legend
percentdropAll = calc_perc(dropctsAll, totalfrAll);
%%
plotperc_dropct(percentdropAll, X, leg, [100 400 800 600], [0 10], 'fig1-perc-test')

%%
plothisto_dropcttotfr(dropctsAll, params.numNodes, params.nodeNames, [100 400 1500 600], [0 750], [0 3000], 'fig1-tot')
%% spatial & temporal diffs (velocity)
%% calc and plot specified pairwise spatial distances
% find node indices for spatial offset comparisons
% Front: pairwise comparisons b/w front nodes (specify)
% body: same except wrt body
idxAnchor = params.nodeidxSB; % wrt scopebase, specify
idxTargets = params.nodeidxSBTargets; % wrt body, specify
idxAll = 1:length(params.nodeNames); % determined auto -
nodepairsIdxAll = combvec(idxAll, idxAnchor);
nodepairsIdx = nodepairsIdxAll(:,idxTargets);
params.nodeIdxcomp = nodepairsIdxAll(:,idxTargets);

[AslpspatSB, spatialdistcumMatSB] = calc_pwdist(Aslpxypairs, params); % SB = wrt scopebase
%% specify position/w/h, XLim, YLim
plothisto_spat(spatialdistcumMatSB, params.nodeNames, nodepairsIdx, [100 400 1800 600], [0 750], [0 300], 'fig3-2')

%% do same for body anchor
idxAnchor = params.nodeidxBody; % wrt body, specify
idxTargets = params.nodeidxBodyTargets; % wrt body, specify
idxAll = 1:length(params.nodeNames); % determined auto -
nodepairsIdxAll = combvec(idxAll, idxAnchor);
nodepairsIdx = nodepairsIdxAll(:,idxTargets);
params.nodeIdxcomp = nodepairsIdxAll(:,idxTargets);

[AslpspatB, spatialdistcumMatB] = calc_pwdist(Aslpxypairs, params); % B = wrt body
%% specify position/w/h, XLim, YLim
plothisto_spat(spatialdistcumMatB, params.nodeNames, nodepairsIdx, [100 400 1800 600], [0 750], [0 300], 'fig4-2')

%% calc and plot specified pairwise temporal distances (velocity)
[temporaldistcumMat, Aslpvel] = calc_tempdist(Aslpxypairs, params.numNodes, params.Aslpfridx);
plothisto_temp(temporaldistcumMat, params.nodeNames, params.numNodes, [100 400 1500 600], [0 100], [0 1000], 'fig5')

%% ------------------------- DO SAME BUT IN SOCIAL PROXIMITY -------------------------------------------------------
%% plot all social distances
histogram(socDistCumAll); ylabel('instances'); xlabel('distance (px)');
title('cumulative histogram of all social distances')
ax = gca; ax.FontSize = 13; ax.YLim = [0 10000];
saveas(gcf,'fig6-soccumdist', 'png'); close

%% dropped nodes social proximity frames
distpx = 50; % set high to include all fr
params.Aslpfridx = getfridx(AslpsocDist, distpx);
[Dropslp, dropctsAll, totalfrAll, totalfrcum] = calc_drop(Aslp, params);

%%
plothisto_dropct(dropctsAll, params.numNodes, params.nodeNames, [100 400 1500 600], [0 100], 'fig1-soc');

X = reordercats(categorical(params.nodeNames), params.nodeNames); % fixed reordering of categories for display
percentdropAll = calc_perc(dropctsAll, totalfrAll);
plotperc_dropct(percentdropAll, X, leg, [100 400 800 600], [0 10], 'fig1-soc-perc')

%%
plothisto_dropcttotfr(dropctsAll, params.numNodes, params.nodeNames, [100 400 1500 600], [0 300], [0 2000], 'fig1-soc-tot')
%% calc and plot specified pairwise spatial distances
% find node indices for spatial offset comparisons
% Front: pairwise comparisons b/w front nodes (specify)
% body: same except wrt body
idxAnchor = params.nodeidxSB; % wrt scopebase, specify
idxTargets = params.nodeidxSBTargets; % wrt body, specify
idxAll = 1:length(params.nodeNames); % determined auto -
nodepairsIdxAll = combvec(idxAll, idxAnchor);
params.nodeIdxcomp = nodepairsIdxAll(:,idxTargets);

[AslpspatSBsoc, spatialdistcumMatSBsoc] = calc_pwdist(Aslpxypairs, params); % SB = wrt scopebase
% specify position/w/h, XLim, YLim
plothisto_spat(spatialdistcumMatSBsoc, params.nodeNames, nodepairsIdx, [100 400 1800 600], [0 500], [0 300], 'fig3-2-soc')

% do same for body anchor
idxAnchor = params.nodeidxBody; % wrt body, specify
idxTargets = params.nodeidxBodyTargets; % wrt body, specify
idxAll = 1:length(params.nodeNames); % determined auto -
nodepairsIdxAll = combvec(idxAll, idxAnchor);
params.nodeIdxcomp = nodepairsIdxAll(:,idxTargets);

[AslpspatBsoc, spatialdistcumMatBsoc] = calc_pwdist(Aslpxypairs, params); % B = wrt body
% specify position/w/h, XLim, YLim
plothisto_spat(spatialdistcumMatBsoc, params.nodeNames, nodepairsIdx, [100 400 1800 600], [0 500], [0 300], 'fig4-2-soc')

% calc and plot specified pairwise temporal distances (velocity)
[temporaldistcumMatsoc, Aslpvelsoc] = calc_tempdist(Aslpxypairs, params.numNodes, Aslpfridx);
plothisto_temp(temporaldistcumMatsoc, params.nodeNames, params.numNodes, [100 400 1500 600], [0 100], [0 1000], 'fig5-soc')

%% ----- CLEAN, INTERPOLATE (rest wrapped) -------------------------------------------------------------------
%% match filenames
AslpFilenames = wrap_eaGrp(Slp, @extract_classification, {'filename'});
AslpPaths = wrap_eaGrp(Slp, @extract_classification, {'filepath'});

cd(vidrootpath)
dvidfiles = dir('*_crop.avi');
filenameVidPaths = join(horzcat({dvidfiles.folder}', {dvidfiles.name}'), filesep);
filenameVidPaths = filenameVidPaths(startsWith({dvidfiles.name}', 'Pair'));
AvidPaths = wrap_eaGrp(AslpFilenames, @map_vidPaths, filenameVidPaths);

%% clean, interpolate
idxAnchor = params.nodeidxBody; % wrt body, specify
idxTargets = params.nodeidxBodyTargets; % wrt body, specify
idxAll = 1:length(params.nodeNames); % determined auto -
nodepairsIdxAll = combvec(idxAll, idxAnchor);
params.nodeIdxcomp = nodepairsIdxAll(:,idxTargets); % comparisons used
params.threshdistpx = 140; % distance to threshold out outliers
params.compTF = 'greater'; % threshold by distance comparator
params.interpTF = 'True'; % whether to interpolate or not, convert_xy()
params.threshdrp = 5; % threshold for looking at frames dropped events

assert(~isfield(params, 'Aslpfridx')) % not a field => consider all fr
AslpAll = wrap_eaGrp(Slp, @extract_classification, {'track', 'track'}); % preserve Aslp name for use above
AslpDropAll = wrap_eaGrp(AslpAll, @calc_drop, params); % conditionally takes Aslpfridx
AslpxyInt = wrap_eaGrp(AslpAll, @convert_xy, params); % interp first w/o delete to find distance
AslpspatBAll = wrap_eaGrp(AslpxyInt, @calc_pwdist, params); % conditionally takes Aslpfridx
AslpfridxspatAll = wrap_eaGrp(AslpspatBAll, @getfridxspat, params); % conditionally takes Aslpfridx

%% merge fr vecs for various criteria for deletion
% set up multi analysis struct
% when considering variable number of data inputs
Amulti = struct();
for si = 1:length(Snames) % for each celltype eg cmk, dlx
    C = AslpspatBAll.(Snames{si});
    expifields = fieldnames(C);
    for expi = 1:length(expifields)
        ExpAnalysis = struct(); % clear
        AslpatB = C.(expifields{expi});
        AslpDropi = AslpDropAll.(Snames{si}).(expifields{expi});
        ExpAnalysis.Dropslp = AslpDropi;
        ExpAnalysis.AslpspatB = AslpatB;
        Amulti.(Snames{si}).(expifields{expi}) = ExpAnalysis;
    end
end
Annvec = wrap_eaGrpmulti(Amulti, @get_Annvec, params);

%% use merged vec to delete
ds = struct();
ds.AslpAll = AslpAll;
ds.Annvec = Annvec;
params.analyzedsnames = {'AslpAll', 'Annvec'}; 
script = 'delete_slpAnnvec';
dso = wrap_eaAniMulti(ds, script, params);

%% interpolate
AslpAlldel = dso;
AslpxyInt = wrap_eaGrp(AslpAlldel, @convert_xy, params);

%% annotate and write videos
params.Draw.anncircradius = 4;
params.Draw.inscolors = {'blue', [204, 102, 0], 'magenta', 'green', 'cyan', 'black'}; % up to 5 animal colors
params.Draw.connections = horzcat(combvec(4, [1:3 5 7]), [6;7]); % [2 by # pw connections to draw]

%%
A = struct();
for si = 1:length(Snames) % for each celltype eg cmk, dlx
    C = Annvec.(Snames{si});
    expifields = fieldnames(C);
    for expi = 1:length(expifields)
        A.Annvec = C.(expifields{expi});
        A.AslpxyInti = AslpxyInt.(Snames{si}).(expifields{expi});
        A.AvidPaths = AvidPaths.(Snames{si}).(expifields{expi});
        annPos_genVids(A, params)
    end
end

%% check videos
% videoReader = VideoReader('PairC1_20210402_XZ70_XZ71_behav_exp.avi_crop.avi_ann.avi');
% videoPlayer = vision.VideoPlayer;
% % while hasFrame(videoReader)
frvec = unionAcrossCols([Annvec(7).track1, Annvec(7).track2]);

% ctrl implay:
H = implay('PairD15_20210603_XZ96_XZ93_behav_exp.avi_crop.avi_ann');
Ctrls = H.DataSource.Controls;
for fi = 1:length(frvec) - 1
    c = frvec(fi);
    while c < frvec(fi+1)
        stepFwd(Ctrls);
        c = c + 1;
    end
    disp('paused until input..')
    pause
end
% play(Ctrls)
% vision.VideoPlayer

%% reformat, and reload interp and dropcts back into S struct
assert(strcmp(params.interpTF, 'True'))
for si = 1:length(Snames) % for each celltype eg cmk, dlx
    expifields = fieldnames(Slp.(Snames{si}));
    for expi = 1:length(expifields)
        AslpxyInti = AslpxyInt.(Snames{si}).(expifields{expi}); % load from above
        AslpDropi = AslpDropAll.(Snames{si}).(expifields{expi});
        AslpspatBAlli = AslpspatBAll.(Snames{si}).(expifields{expi});
        AslpfridxspatAlli = AslpfridxspatAll.(Snames{si}).(expifields{expi}); % load from above
        Annveci = Annvec.(Snames{si}).(expifields{expi});
        for pi = 1:length(AslpxyInti)
            fields = fieldnames(AslpxyInti);
            if ~isempty(AslpxyInti(pi).(fields{ai}))
                for ai = 1:length(fields)
                    xy = AslpxyInti(pi).(fields{ai}); % interp
                    xy4d = zeros(size(xy{:,1}, 1), size(xy, 2), 2); % fr x node x axis (x, y)
                    for ni = 1:size(xy, 2)
                        xy4d(:,ni,1) = xy{:,ni}(:,1); % reformat back to orig
                        xy4d(:,ni,2) = xy{:,ni}(:,2);
                    end
                    % store
                    Slp.(Snames{si}).(expifields{expi}){pi}{ai}.trackInterp = xy4d;
                    Slp.(Snames{si}).(expifields{expi}){pi}{ai}.trackInterpXY = xy;
                    
                    % also add other interp, deletion criteria eg drop, spatial outliers
                    Slp.(Snames{si}).(expifields{expi}){pi}{ai}.DropEvents = AslpDropi(pi).(fields{ai});
                    Slp.(Snames{si}).(expifields{expi}){pi}{ai}.DistancebwNodes = AslpspatBAlli(pi).(fields{ai});
                    Slp.(Snames{si}).(expifields{expi}){pi}{ai}.DistanceOutliers = AslpfridxspatAlli(pi).(fields{ai});
                    Slp.(Snames{si}).(expifields{expi}){pi}{ai}.MergedInterpVec = Annveci(pi).(fields{ai});
                    Slp.(Snames{si}).(expifields{expi}){pi}{ai}.validationParams = params;
                    
                    
                end
            end
        end
    end
end

%% threshold out specific events and find frame indices, write to video only such frames
params.Draw.anncircradius = 3;
params.Draw.inscolors = {'blue', [204, 102, 0], 'magenta', 'green', 'cyan', 'black'}; % up to 5 animal colors
A = struct();
% 
for si = 1:length(Snames) % for each celltype eg cmk, dlx
    expifields = fieldnames(Annvec.(Snames{si}));
    for expi = 1:length(expifields)
        
        Dropslp = AslpDropAll.(Snames{si}).(expifields{expi});
        Aslpfridxspat = AslpfridxspatAll.(Snames{si}).(expifields{expi});
        A.Annveci = Annvec.(Snames{si}).(expifields{expi});
        A.AslpxyInti = AslpxyInt.(Snames{si}).(expifields{expi});
        A.AvidPathsi = AvidPaths.(Snames{si}).(expifields{expi});
        
        for gi = 1:length(Dropslp)
            fields = fieldnames(Dropslp);
            % initialize
            if ~isempty(Dropslp(gi).(fields{1}))
                outFrames = zeros(size(Slp.(Snames{si}).(expifields{expi}){gi}{1}.trackInterp, 1), 1); % n x 1
                Trackmat = zeros(size(Slp.(Snames{si}).(expifields{expi}){gi}{1}.trackInterp)); % orig format trackmat
            end
            for ai = 1:length(fields)
                if ~isempty(Dropslp(gi).(fields{ai}))
                    % load from standard struct
                    Sint = Slp.(Snames{si}).(expifields{expi}){gi}{ai}.trackInterp;
                    % get and threshold dropped fr vectors
                    drp = Dropslp(gi).(fields{ai});
                    threshnans = zeros(size(Dropslp(gi).(fields{ai}).Nans));
                    for ni = 1:length(drp.dropcts)
                        eventidx = drp.dropcts{ni} > params.threshdrp; % limit to events above certain length
                        friidx = drp.startIdx{ni}(eventidx); % fr start indices for such events
                        drplengths = drp.dropcts{ni}(eventidx); % lengths of events
                        for frii = 1:length(friidx) % set only such events to nan
                            threshnans(friidx(frii) : friidx(frii) + drplengths(frii), ni, :) = nan;
                        end
                    end
                    %%
                    [~, DrpIdx] = unionAcrossCols(isnan(threshnans));
                    % also spatially thresholded outliers
                    [~, SpatIdx] = unionAcrossCols(Aslpfridxspat(gi).(fields{ai})); % outlier idx
                    
                    [DrpIdx, SpatIdx] = rmvFrames(DrpIdx, SpatIdx); 
                    [~, mergeIdx] = unionAcrossCols([DrpIdx, SpatIdx]);
                    [~, outFrames] = unionAcrossCols([outFrames, mergeIdx]);% union b/w ani, drpidx, and spatidx
                    
                    % store indices and orig format trackmat
                    % trackmat for making slp file is inverse to fridx for making video
                    NotDrpIdx = ~DrpIdx; % remove nondropped (good) frames
                    NotSpatIdx = ~SpatIdx; % remove non-outlier frames
                    Sint(NotSpatIdx,:,:) = nan;
                    Sint(NotDrpIdx,:,:) = nan;
                    Trackmat(:,:,:,ai) = Sint; % store along 4th dim
                end % if not empty
                
            end % animal
            A.AnnOutFr{gi} = find(outFrames == 1); % merge w/ above params

        end % gi
        % write to video, limit frames
        annPos_genVids(A, params)
            
        %% write to h5, incomplete
        % modify existing so can put back into slp gui and relabel (create - str process error)
        % need to make new paths
%             h5write('PairD11_test.h5', '/tracks', Trackmat)
%             h5write('PairD11_test.h5', '/track_occupancy', Slp.(Snames{si}).exp{gi}{ai}.occupancy_matrix)
    end
end

%% write to h5, incomplete
% modify existing so can put back into slp gui and relabel (create - str process error)
for si = 1:length(Snames) % for each celltype eg cmk, dlx
    for gi = 1:length(Dropslp)
        h5write('PairD11_test.h5', '/tracks', Trackmat)
        h5write('PairD11_test.h5', '/track_occupancy', Slp.(Snames{si}).exp{gi}{ai}.occupancy_matrix)
    end
end























