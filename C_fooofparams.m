% initial  analyses of FOOOF data
clear
subfiles = dir("/Volumes/Hera/Abby/FOOOF/SubjectFiles/*.mat");
all_aps = [];
all_peaks = [];
all_fits = [];
all_gauss = [];

for currsub = 1:length(subfiles)
    % get id and date
    filename = subfiles(currsub).name;
    filepath = subfiles(currsub).folder;
    split_name = split(filename,"_");
    id = str2double(split_name{1});
    date = split_name{2};
    split_date = split(date,".");
    date = str2double(split_date{1});
    
    % load subject file
    load(fullfile(filepath,filename))
    
    % add 1 to index columns b/c python does 0 indexing
    ap_params(:,3:4) = ap_params(:,3:4)+1;
    fit_params(:,3:4) = fit_params(:,3:4) + 1;
    peak_params(:,4:6) = peak_params(:,4:6) + 1;
    gauss_params(:,4:6) = gauss_params(:,4:6) + 1;
    trial_score(:,1) = cellfun(@(x) x + 1, trial_score(:,1),'UniformOutput',false);
    
    % numeric code for scores
    score_map = containers.Map({'2_cor','2_incor','2_errcor','2_drop'},{1,0,2,-1});
    trial_score(:,3) = cellfun(@(x) score_map(x), trial_score(:,2), 'UniformOutput',false); 
    
    % add id and date
    temp_aps = [repmat(id,length(ap_params),1), repmat(date,length(ap_params),1), fit_params(:,1:2), ap_params];
    temp_peaks = [repmat(id,length(peak_params),1), repmat(date,length(peak_params),1), gauss_params(:,1:3), peak_params];    
    
    % add trial score to temp matrices
    trial_epochs = cell2mat(trial_score(:,1));
    trial_scores = cell2mat(trial_score(:,3));
    aps_epochs = temp_aps(:,7);
    [~,idx] = ismember(aps_epochs,trial_epochs);
    temp_aps(:,9) = trial_scores(idx(idx~=-0));
    peaks_epochs = temp_peaks(:,10);
    [~,idx] = ismember(peaks_epochs,trial_epochs);
    temp_peaks(:,12) = trial_scores(idx(idx~=0));
    
    % add to big matrix for all subjects
    all_aps = [all_aps; temp_aps];
    all_peaks = [all_peaks; temp_peaks];
    
    % clear temp matricies
    clear temp_aps temp_peaks 
end

%% save as tables and name columns
fit_cols_str = deblank(string(fit_cols))';
gauss_cols_str = deblank(string(gauss_cols))';
all_aps_table = array2table(all_aps, 'VariableNames',["lunaid" "eeg_date" fit_cols_str(1:2) deblank(string(ap_cols))' "score"]);
all_peaks_table = array2table(all_peaks,'VariableNames',["lunaid" "eeg_date" gauss_cols_str(1:3) deblank(string(peak_cols))' "score"]);

writetable(all_aps_table,"/Volumes/Hera/Abby/FOOOF/aps.csv")
writetable(all_peaks_table,"/Volumes/Hera/Abby/FOOOF/peaks.csv")


