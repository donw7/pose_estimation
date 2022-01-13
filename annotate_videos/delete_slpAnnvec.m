% script delete_slpAnnvec
% expects anids, per animal workspace
% eg
%

%%
out = []; % keep same name to assign to dso

%%
fields = {};
fields = fieldnames(anids);
vecs = cell(1, size(fields, 1));

trackmat = anids.(analyzefields{1}); % fri by nodes by x or y
annvec = [];
annvec = anids.Annvec.vec;

trackmatdel = [];

%% delete
trackmatX = trackmat(:,:,1);
trackmatY = trackmat(:,:,2);
trackmatX(annvec) = nan;
trackmatY(annvec) = nan;

trackmatdel = cat(3, trackmatX, trackmatY);


%%
out = trackmatdel;