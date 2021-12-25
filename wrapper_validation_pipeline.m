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
% draw_lines(Slp{1}{1}.track, Slp{1}{2}.track, 1:500, [0 1000 -1000 0])
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
[temporaldistcumMatsoc, Aslpvelsoc] = calc_tempdist(Aslpxypairs, params.numNodes, params.Aslpfridx);
plothisto_temp(temporaldistcumMatsoc, params.nodeNames, params.numNodes, [100 400 1500 600], [0 100], [0 1000], 'fig5-soc')

%% ----- CLEAN, INTERPOLATE (rest wrapped) -------------------------------------------------------------------
%% set params and find vid filenames
% assume single level - for multilevel groups (e.g. various phenotypes or
% genotypes), use wrap_ series e.g. wrap_eaGrp():
% e.g. AslpAll = wrap_eaGrp(Slp, @extract_classification, {'track', 'track'})
% AslpDropAll = wrap_eaGrp(AslpAll, @calc_drop, params);

AslpFilenames = extract_classification(Slp, {'filename'});
idxAnchor = params.nodeidxBody; % wrt body, specify
idxTargets = params.nodeidxBodyTargets; % wrt body, specify
idxAll = 1:length(params.nodeNames); % determined auto -
nodepairsIdxAll = combvec(idxAll, idxAnchor);
params.nodeIdxcomp = nodepairsIdxAll(:,idxTargets); % comparisons used
params.threshdistpx = 90; % distance to threshold out outliers
params.compTF = 'greater'; % threshold by distance comparator
params.interpTF = 'True'; % whether to interpolate or not, convert_xy()
assert(~isfield(params, 'Aslpfridx')) % not a field => consider all fr
AslpxyInt = convert_xy(Aslp, params);
Aslpfridxspat = getfridxspat(AslpspatB, params);

%% match filenames, find vidpaths
cd(vidrootpath)
dvidfiles = dir('*.mp4');
filenameVidPaths = join(horzcat({dvidfiles.folder}', {dvidfiles.name}'), filesep);
filenameVidPaths = filenameVidPaths(startsWith({dvidfiles.name}', 'Pair'));
AvidPaths = map_vidPaths(AslpFilenames, filenameVidPaths);

%% annotate and write whole videos
params.Draw.anncircradius = 3;
params.Draw.inscolors = {'blue', [204, 102, 0], 'magenta', 'green', 'cyan', 'black'}; % up to 5 animal colors
params.Draw.connections = horzcat(combvec(1, [2 3 5]), [5; 4]); % [2 by # pw connections to draw]
params.threshdrp = 10; % threshold for looking at frames dropped events

% get annotation vectors
ExpAnalysis = struct();
ExpAnalysis.Dropslp = Dropslp; % order matters
ExpAnalysis.AslpspatB = AslpspatB;
Annvec = get_Annvec(ExpAnalysis, params);

% use merged ann vec to delete
ds = struct();
ds.Aslp = Aslp;
ds.Annvec = Annvec;
params.analyzedsnames = {'Aslp', 'Annvec'}; 
script = 'delete_slpAnnvec';
sh = str2func(script);

dso = struct();
analyzefields = fieldnames(ds);
firstfields = fieldnames(ds.(analyzefields{1}));

for gi = 1:length(Aslp)
    aifields = fieldnames(Aslp);
    for ai = 1:length(aifields)
        if ~isempty(Aslp(gi).(aifields{ai}))
            anids = struct(); % clear
            % analyzenames
            for anai = 1:length(analyzefields)
                anids.(analyzefields{anai}) =... % workspace for per animal analysis
                    ds.(analyzefields{anai})(gi).(aifields{ai});
            end

            % run
            sh()

            % out
            dso.(analyzefields{1})(gi).(aifields{ai}) = out;
        end
    end

end

%% interpolate
Aslpdel = dso.Aslp;
AslpxyInt = convert_xy(Aslpdel, params);

%% store annvecs
A = struct();
A.Annveci = Annvec;
A.AslpxyInti = AslpxyInt;
A.AvidPathsi = AvidPaths;

%% write complete videos
annPos_genVids(A, params)

%% reformat, and reload interp and dropcts back into S struct
assert(strcmp(params.interpTF, 'True'))
for gi = 1:length(AslpxyInt)
    fields = fieldnames(AslpxyInt);
    if ~isempty(AslpxyInt(gi).(fields{ai}))
        for ai = 1:length(fields)
            xy = AslpxyInt(gi).(fields{ai}); % interp
            xy4d = zeros(size(xy{:,1}, 1), size(xy, 2), 2); % fr x node x axis (x, y)
            for ni = 1:size(xy, 2)
                xy4d(:,ni,1) = xy{:,ni}(:,1); % reformat back to orig
                xy4d(:,ni,2) = xy{:,ni}(:,2);
            end
            % store
            Slp{gi}{ai}.trackInterp = xy4d;
            Slp{gi}{ai}.trackInterpXY = xy;
            
            % also save other categories
            Slp{gi}{ai}.trackDrops = Dropslp(gi).(fields{ai});           
            Slp{gi}{ai}.DistancebwNodes = AslpspatB(gi).(fields{ai});
            Slp{gi}{ai}.DistanceOutliers = Aslpfridxspat(gi).(fields{ai});
            Slp{gi}{ai}.MergedInterpVec = Annvec(gi).(fields{ai});
            Slp{gi}{ai}.validationParams = params;
        end
    end
end

%% threshold out specific events and find frame indices, write to video only such frames
for gi = 1:length(Dropslp)
    fields = fieldnames(Dropslp);

    outFrames = zeros(size(Aslp(gi).(fields{ai}), 1), 1); % n x 1
    Trackmat = zeros(size(Aslp(gi).(fields{ai}))); % orig format trackmat

    for ai = 1:length(fields)
        if ~isempty(Dropslp(gi).(fields{ai}))
            %% 
            Sint = Slp{gi}{ai}.trackInterp;
            %%
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

            [DrpIdx, SpatIdx] = rmvFrames(DrpIdx, SpatIdx); % could be diff lengths?
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
    
    %% write to h5, incomplete
    % modify existing so can put back into slp gui and relabel (create - str process error)
    % need to make new paths
%     h5write('PairD11_test.h5', '/tracks', Trackmat)
%     h5write('PairD11_test.h5', '/track_occupancy', Slp.(Snames{si}).exp{gi}{ai}.occupancy_matrix)

end % gi
% write to video, limit frames
annPos_genVids(A, params)

%% write back to h5
% modify existing so can put back into slp gui and relabel (create - str process error)
for si = 1:length(Snames) % for each celltype eg cmk, dlx
    for gi = 1:length(Dropslp)
        h5write('PairD11_test.h5', '/tracks', Trackmat)
        h5write('PairD11_test.h5', '/track_occupancy', Slp.(Snames{si}).exp{gi}{ai}.occupancy_matrix)
    end
end























