%% multi-directory_loader
% compile paths with various folder structures

%% paths
mainfolderName = 'beh processed data';
% expfolderName = 
% Snames = 
% pairIDnames = {'C', 'D'};
% numExp = [13 16];

if ispc
    rootPathgd = 'G:\My Drive'; % for pc operating systems
    rootPathHD = 'F:\Don';
else
    rootPathgd = '/Volumes/GoogleDrive/My Drive'; % For unix (mac, linux) operating systems
end

sleapPath = [rootPathHD filesep mainfolderName filesep expfolderName];
cd(sleapPath)

dh5filesraw = dir('**/*predictions.analysis.h5'); % searches current and all subfolders
filenamesAnalysis = {dh5files.name}';
filenamesAnalysisraw = {dh5filesraw.name}';

% source video paths
vidrootpath = 'F:\Don\beh processed data';

% match sleap filenames to pairIDs already in filename
celltypeIDs = match_filename(filenamesAnalysis, 'Pair[A-Z]'); % search first 7 letters, map IDs onto sequence found in filenames list
pairIDs = match_filename(filenamesAnalysis, '\d+');
pairIDsnumber = match_filename(filenamesAnalysis, '\d+')';