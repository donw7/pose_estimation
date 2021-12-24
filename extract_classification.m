function extractedS = extract_classification(CMK, Xnames, Ynames)
% extracts from nested structure and outputs into nonscalar struct
% Xnames as in predictor variables, can be variable length eg as cell string 
    % eg {'DM', 'DM'}
% Ynames is response variable, accepts 1 or multiple inputs

%% output nonscalar struct
extractedS = struct();
nameXi = strcat(Xnames{1}, '1');
extractedS(1).(nameXi) = {};
% nameYi = strcat(Ynames{1}, '1');
% extractedS(2).(nameYi) = [];

%% extract from O struct and consolidate into extracted nonscalar structure
if nargin == 2 % CMK and Xnames only
    %%
    for e = 1:size(CMK, 2)
        for x = 1:length(Xnames)
            if ~isempty(CMK{e}) && ~isempty(CMK{e}{1})
                name = strcat(Xnames{x}, num2str(x));
                extractedS(e).(name) = CMK{e}{x}.(Xnames{x});
            end
        end
    end
    
elseif nargin > 2
    %%
    if length(Ynames) == 1

        for e = 1:size(CMK, 2)
            for x = 1:length(Xnames)
                if ~isempty(CMK{e}) && ~isempty(CMK{e}{1})
                    name = strcat(Xnames{x}, num2str(x));
                    extractedS(e).(name) = CMK{e}{x}.(Xnames{x});
                end
            end
            extractedS(e).(Ynames{1}) = CMK{e}{1}.(Ynames{1}); % extracts from first animal only
        end
    elseif length(Ynames) > 1
        %%
        for e = 1:size(CMK, 2)
            for x = 1:length(Xnames)
                if ~isempty(CMK{e}) && ~isempty(CMK{e}{1})
                    name = strcat(Xnames{x}, num2str(x));
                    extractedS(e).(name) = CMK{e}{x}.(Xnames{x});
                end
            end
            for y = 1:length(Ynames)
                if ~isempty(CMK{e}) && ~isempty(CMK{e}{1})
                    name = strcat(Ynames{y}, num2str(y));
                    extractedS(e).(name) = CMK{e}{y}.(Ynames{y}); % extracts from all animals
                end
            end
        end
    end
end
end