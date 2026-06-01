% initial  analyses of FOOOF data
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
    
    % add id and date
    temp_aps = [repmat(id,length(ap_params),1), repmat(date,length(ap_params),1), ap_params];
    temp_peaks = [repmat(id,length(peak_params),1), repmat(date,length(peak_params),1), peak_params];
    temp_gauss = [repmat(id,length(gauss_params),1), repmat(date,length(gauss_params),1), gauss_params];
    temp_fit = [repmat(id,length(fit_params),1), repmat(date,length(fit_params),1), fit_params];
    
    % add to big matrix for all subjects
    all_aps = [all_aps; temp_aps];
    all_peaks = [all_peaks; temp_peaks];
    all_fits = [all_fits; temp_fit];
    all_gauss = [all_gauss; temp_gauss];
    
    % clear temp matricies
    clear temp_aps temp_peaks temp_gauss temp_fit    
end

%% save as tables and name columns
all_aps_table = array2table(all_aps, 'VariableNames',["lunaid" "eeg_date" deblank(string(ap_cols))']);
all_peaks_table = array2table(all_peaks,'VariableNames',["lunaid" "eeg_date" deblank(string(peak_cols))']);
all_gauss_table = array2table(all_gauss,'VariableNames',["lunaid" "eeg_date" deblank(string(gauss_cols))']);
all_fits_table = array2table(all_fits, 'VariableNames',["lunaid" "eeg_date" deblank(string(fit_cols))']);

writetable(all_aps_table,"/Volumes/Hera/Abby/FOOOF/aps.csv")
writetable(all_peaks_table,"/Volumes/Hera/Abby/FOOOF/peaks.csv")
writetable(all_gauss_table,"/Volumes/Hera/Abby/FOOOF/gauss.csv")
writetable(all_fits_table,"/Volumes/Hera/Abby/FOOOF/fits.csv")


