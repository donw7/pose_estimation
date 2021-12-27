function Xresh = resh_binmeanNan(X, q)
% reshape matrix along first dim to binned mean, ignoring nans
% q = bin length

%%
n = size(X, 1);
multiple = n - mod(n,q);
targetDim = multiple / q;
Xcut = X(1:multiple,:);
Xresh = [];
slices = [1:q:multiple];

for fri = 1:length(slices) - 1
    slice = slices(fri): slices(fri+1);
    Xbin = Xcut(slice, :);
    Xresh = vertcat(Xresh, nanmean(Xbin));
    
end

% X3 = reshape(X(1:multiple,:), targetDim, [], q); % [multiple x targetDim x size(X,2)]
% Y1 = nanmean(X3, 3);
% Y = reshape(Y1, targetDim, []);

end
