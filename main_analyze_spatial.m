%% main script for analysis of spatial relationships
% inherits from wrapper_validation params and AslpxyInt

%% set parameters
params.reshapeBinlen = 10;
% params.reshapeBinlen = 5; % if diff fps

%% reshape w/ median smooth, unnest, store
% Aslpxyresh = struct();
for gi = 1%:length(AslpxyInt)% iterate by available filenames
    aifields = fieldnames(AslpxyInt(gi));
    for ai = 1:length(aifields)
        if ~isempty(AslpxyInt(gi).(aifields{ai}))
            xy = AslpxyInt(gi).(aifields{ai});
            xy = [xy{:}];
            xymedreshaped = resh_binmeanNan(xy, params.reshapeBinlen);
            Aslpxyresh(gi).(aifields{ai}) = xymedreshaped; % store unnested
            
            % reshape social frames idx to match
            socfridx = params.Aslpfridx(gi).("idx");
            socfridxreshaped = resh_binmeanNan(socfridx, params.reshapeBinlen);
            Aslpfridxresh(gi).("idx") = socfridxreshaped;
        end
    end
end

%% generate heat maps to visualize
% params.res = [2 2];
params.res = [4 4];
for gi = 1:length(Aslpxyresh) %[1 2 8:21]%
    aifields = fieldnames(Aslpxyresh(gi));
    if ~isempty(Aslpxyresh(gi).(aifields{1}))
        %%
        tic
        xymeanAll = cell(1,length(aifields));
        gridvecsocfilAll = cell(1,length(aifields));
        
        for ai = 1:length(aifields)
            xy = Aslpxyresh(gi).(aifields{ai});
            yidx = params.nodeidxHead .* 2; % eg 2 4 6
            xidx = yidx - 1; % eg 1 3 5
            xymean = [nanmean(xy(:,xidx), 2), nanmean(xy(:,yidx), 2)];
            [gridlv, gridvec, gridpos] = convert_xy2gridvec(xymean, params.res);
            % filtering for social events
            socfridx = logical(floor(Aslpfridxresh(gi).("idx")));
            gridvecsocfil = zeros(size(gridvec,1),1);
            gridvecsocfil(socfridx) = gridvec(socfridx);
            xymeanAll{1,ai} = xymean;
            gridvecsocfilAll{1,ai} = gridvecsocfil;
        end
        
        outVid = VideoWriter(['hm_' num2str(gi)]);%AslpID(gi).pairID1{:} '_' IDs.genotypePair{gi}]); % prep video container
        HMcum = zeros(params.res(1), params.res(2));
        open(outVid)
        
        %
        fh = figure;
        set(fh, 'position', [0 500 1500 300]);
        sp1 = subplot(1,3,1);
        alineAll = animatedline('Color','b', 'LineStyle', '--', 'Marker', 'o', 'MarkerSize', 3);
        blineAll = animatedline('Color','r', 'LineStyle', '--', 'Marker', 'o', 'MarkerSize', 3);

        for ti = 1:length(gridvecsocfilAll{1})
            %% all frames
            subplot(1,3,1)
            
            axis([200 900 100 800])
            addpoints(alineAll, xymeanAll{1}(ti,1), xymeanAll{1}(ti,2)); % x coordinate, y coordinate
            addpoints(blineAll, xymeanAll{2}(ti,1), xymeanAll{2}(ti,2));
            drawnow
            title('position, cumulative');

            
            %% only social proximal frames
            sp2 = subplot(1,3,2);
            aline = animatedline('Color','b', 'LineStyle', ':', 'Marker', 'o', 'MarkerSize', 3);
            bline = animatedline('Color','r', 'LineStyle', ':', 'Marker', 'o', 'MarkerSize', 3);
            
            axis([200 900 100 800])
            if gridvecsocfilAll{1}(ti) ~= 0
                addpoints(aline, xymeanAll{1}(ti,1), xymeanAll{1}(ti,2)); % x coordinate, y coordinate
                addpoints(bline, xymeanAll{2}(ti,1), xymeanAll{2}(ti,2));
                drawnow
            end
            title('position, head proximity <40 px only, cumulative');
            
            
            %% heat map, cumulative, only social frames
            sp3 = subplot(1,3,3);
            if gridvecsocfilAll{1}(ti) == 0
                h = heatmap(HMcum, 'colormap', jet);
            else
                [subrow, subcol] = ind2sub(params.res, gridvecsocfilAll{1}(ti));
                HMcum(subrow, subcol) = HMcum(subrow, subcol) + 1;
%                 [subrow, subcol] = ind2sub(params.res, gridvecsocfilAll{2}(ti));
%                 HMcum(subrow, subcol) = HMcum(subrow, subcol) + 1;
                
                h = heatmap(HMcum, 'colormap', jet);
            end
            title('head proximity <40 px only, cumulative');
            I = getframe(gcf);
            writeVideo(outVid, I)
                    
            
        end
        close(outVid);
        fprintf('Printed video for Pair %d\n', gi)
        toc
        close all
        
        % plot real time position
        % save name, geno
        % save HMcum
        % save plot traces
                
    end
end
%HM over time