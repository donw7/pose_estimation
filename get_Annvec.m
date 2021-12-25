function Annvec = get_Annvec(ExpAnalysis, p)
Annvec = struct();

EAfields = fieldnames(ExpAnalysis);

% based on existence of relevant params fields
if isfield(p, 'threshdrp')
    EAidx = strcmp(EAfields, 'Dropslp');
    Dropslp = ExpAnalysis.(EAfields{EAidx});
    for gi = 1:length(Dropslp)
        fields = fieldnames(Dropslp);
        for ai = 1:length(fields)
            if ~isempty(Dropslp(gi).(fields{ai}))
                matlv = Dropslp(gi).(fields{ai}).Nans(:,:,1); % take x only, isequal to y
                Annvec(gi).(fields{ai}).vec = matlv;
                labels = cell(size(matlv));
                labels(matlv) = deal({'intrp-drp'});
                Annvec(gi).(fields{ai}).labels = labels;
            end
        end
    end
end

if isfield(p, 'threshdistpx')
    EAidx = strcmp(EAfields, 'AslpspatB');
    AslpspatB = ExpAnalysis.(EAfields{EAidx});
    Aslpfridxspat = getfridxspat(AslpspatB, p); % greater or less than
    
    for gi = 1:length(Aslpfridxspat)
        fields = fieldnames(Aslpfridxspat);
        for ai = 1:length(fields)
            if ~isempty(Aslpfridxspat(gi).(fields{ai}))
                matlv = Aslpfridxspat(gi).(fields{ai}); % take x only, isequal to y
                prev = Annvec(gi).(fields{ai}).vec;
                new = prev;
                new(matlv) = 1; % add vecs based on threshold distance
                Annvec(gi).(fields{ai}).vec = new;
                labels = Annvec(gi).(fields{ai}).labels; % load previous
                labels(matlv) = deal({'intrp-dist'});
                Annvec(gi).(fields{ai}).labels = labels;
            end
        end
    end
end