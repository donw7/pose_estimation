function pairIDs = match_filename(cellStrs, REpattern)
% maps eg filename index onto pairIDs
% takes cell array of strings, searches letters defined by strIdx for
% regexp pattern defined by repattern
% eg cellStrs = filenamesAnalysis
% eg strIdx = 1:7
% eg REpattern = 'Pair[A-Z]'
% eg REpattern = '\d+'
% outputs cell array for strings or mat for numbers

% pairIDs = cell(1, length(cellStrs)); don't preallocate ~ conditional
% format as below

for fi = 1:length(cellStrs)
    filename = cellStrs{fi}; % take first few letters of filename as experiment ID
    extracted = extract(filename, regexpPattern(REpattern));
    isNumber = ~isnan(str2double(extracted{1}));
    if isNumber
        pairIDs(fi) = str2double(extracted{1});
    else
        pairIDs{fi,:} = extracted{1};
end

%% incomplete - consider another case to match by animal IDs
% %% match sleap filenames to expIDs (matches for pair ID but not individual animal ID)
% sleapfilenameIDs = zeros(1, length(filenamesAnalysis));
% counteridx = 1;
% while counteridx <= length(filenamesAnalysis)
%     expID = filenamesAnalysis{counteridx}(1:11); % take first few letters of filename as experiment ID
%     for ei = 1:length(expID_Animal1) 
%         if strcmp(expID_Animal1(1, ei), expID) | strcmp(expID_Animal2(ei), expID)
%             sleapfilenameIDs(counteridx) = ei;
%         end
%     end
%     counteridx = counteridx + 1;
% end

end



