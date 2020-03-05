function Microstim_data_processing(varargin)
% this program is for preparing all data for behavioural microstimulation analysis

%Microstim_data_processing('Curius','Y');
%Microstim_data_processing('Curius','L');
%Microstim_data_processing('Linus','Y');
%Microstim_data_processing('Linus','L');

global analysis_parameters
if nargin==0
    analysis_parameters.monkey='Linus';
    analysis_parameters.DAG_Drive='L';                      % important to get the correct dropbox folder
else
    analysis_parameters.monkey      = varargin{1}; %
    analysis_parameters.DAG_Drive   = varargin{2};                      % important to get the correct dropbox folder
end

analysis_parameters.user=getUserName;                      % important to get the correct dropbox folder
analysis_parameters.dates=[20130000 20160000];  % basically update all!
analysis_parameters.folders.ignore_not_matching_dates =1;
analysis_parameters.overwrite_all_raw_data = true;

%analysis_parameters.mastertable_types={'trialinfo'}; % {'trialinfo','bias_RT','evoked'}
analysis_parameters.mastertable_types={'trialinfo'};

disp(['Drive= ', analysis_parameters.DAG_Drive, 'User= ', analysis_parameters.user, 'Monkey= ', analysis_parameters.monkey]);

switch analysis_parameters.monkey
    case 'Linus'
        analysis_parameters.folders.original_data      =strcat(analysis_parameters.DAG_Drive,':', filesep,'Data', filesep,  analysis_parameters.monkey, filesep);
        analysis_parameters.folders.initial_folder      =strcat(analysis_parameters.DAG_Drive,':', filesep,'Data', filesep,  [analysis_parameters.monkey '_microstim'], filesep);
    case 'Curius'
        analysis_parameters.folders.original_data      =strcat(analysis_parameters.DAG_Drive,':', filesep,'Data', filesep,  analysis_parameters.monkey, filesep, 'setup2_microstim', filesep);
        analysis_parameters.folders.initial_folder      =strcat(analysis_parameters.DAG_Drive,':', filesep,'Data', filesep,  [analysis_parameters.monkey '_microstim'], filesep);
end

session_folders_to_move=dir(analysis_parameters.folders.initial_folder);
session_folders_to_move=session_folders_to_move([session_folders_to_move.isdir]);
for session_idx=3:numel(session_folders_to_move)    
        copyfile (strcat(analysis_parameters.folders.initial_folder, session_folders_to_move(session_idx).name, filesep), strcat(analysis_parameters.folders.original_data, session_folders_to_move(session_idx).name, filesep));
end

analysis_parameters.folders.dropbox              = strcat('C:', filesep,'Users', filesep,  analysis_parameters.user, filesep, 'Dropbox', filesep, 'DAG', filesep, 'microstim_behavior', filesep, analysis_parameters.monkey, '_microstim_dpz', filesep);
analysis_parameters.folders.extended_data        = strcat(analysis_parameters.DAG_Drive,':', filesep,'Data', filesep,  analysis_parameters.monkey, '_microstim_with_parameters', filesep);
analysis_parameters.folders.output               = strcat(analysis_parameters.DAG_Drive,':', filesep,'microstim_behavior', filesep,  analysis_parameters.monkey, '_summaries', filesep);
analysis_parameters.files.original_parameters    = strcat(analysis_parameters.monkey, '_stimulation_parameters.xls');
analysis_parameters.files.updated_parameters     = strcat(analysis_parameters.monkey, '_updated_parameters.xls');
analysis_parameters.folders.code                 = strcat(analysis_parameters.DAG_Drive,':', filesep,'microstim_behavior', filesep, 'code', filesep);
analysis_parameters.files.monkeypsych_analyze    = 'monkeypsych_analyze_working';

%analysis_parameters.folders.original_data      =strcat(current_dir(1),':', filesep,'Data', filesep,  analysis_parameters.monkey, filesep);
%analysis_parameters.folders.code               =strcat(analysis_parameters.DAG_Drive,':', filesep,'microstim_behavior', filesep, 'code', filesep, 'After 20140725', filesep);

copyfile (strcat(analysis_parameters.folders.dropbox, analysis_parameters.files.original_parameters), strcat(analysis_parameters.folders.output, analysis_parameters.files.original_parameters));
[~,parametersheets] = xlsfinfo(strcat(analysis_parameters.folders.output, analysis_parameters.files.original_parameters));

working_folders=dir(analysis_parameters.folders.extended_data);
original_folders=dir(analysis_parameters.folders.original_data);
original_folder_names={original_folders.name};
working_folder_names={working_folders.name};
for session_idx=1:numel(parametersheets)
    if ismember(parametersheets(session_idx), original_folder_names) && (~ismember(parametersheets(session_idx), working_folder_names) ||  analysis_parameters.overwrite_all_raw_data)
        copyfile (strcat(analysis_parameters.folders.original_data, parametersheets{session_idx}, filesep), strcat(analysis_parameters.folders.extended_data, parametersheets{session_idx}, filesep));
        monkeypsych_clean_data(analysis_parameters.folders.extended_data,str2double(parametersheets{session_idx}))
        import_parameters_from_xls_to_trial_task_per_date(strcat(analysis_parameters.folders.extended_data, parametersheets{session_idx}),parametersheets{session_idx});
    end
end
if any(strcmp(analysis_parameters.mastertable_types,'trialinfo'))
    [analysis_parameters.filelist_complete, analysis_parameters.filelist_formatted] = DAG_get_filelist_from_folder;
    %DAG_import_parameters_from_xls_to_trial_task;
    DAG_get_trialinfo_mastertable;
end
for k=1:numel(analysis_parameters.mastertable_types)
    update_stimulation_parameters(analysis_parameters.mastertable_types{k});
end
%feval(analysis_parameters.files.plot_evoked_saccades);
%evoked_saccade_preparation
%evoked_saccade_calculation
%evoked_saccade_plotting


% trial history plots.....
%(994 677)
%(4084 2131)
% [out_comp variable_to_test counter]= monkeypsych_analyze_20140910(analysis_parameters.filelist_formatted,...
%     {'nsacc_max',10,'summary',1,'type',2,'inferential_on',1,'success',1,'display',1,'choice',1,'microstim',NaN,'n_trials_past',2,'runs_as_batches',1,...
%     'trial_history_mode',1,'current_structure','binary','current_parameter','eyetar_l','past_structure','binary','past_parameter','eyetar_l','past_value',NaN,'past_keys',{'microstim',1,'choice',NaN}});
%

% % trial history plots.....
% %(776 362)
% %(2970 1330)
% [out_comp variable_to_test counter]= monkeypsych_analyze_20140814(analysis_parameters.filelist_formatted,...
%     {'nsacc_max',10,'summary',1,'type',2,'inferential_on',0,'success',1,'display',1,'choice',1,'microstim',NaN,...
%     'trial_history_mode',1,'observed_structure','binary','observed_parameter','eyetar_l','past_structure','binary','past_parameter','eyetar_l','past_value',NaN,'past_keys',{'microstim',0,'choice',NaN}});
%
% % trial history plots.....
% %(1138 1671)
% %(4300 6215)
% [out_comp variable_to_test counter]= monkeypsych_analyze_20140814(analysis_parameters.filelist_formatted,...
%     {'nsacc_max',10,'summary',1,'type',2,'inferential_on',0,'success',1,'display',1,'choice',1,'microstim',NaN,...
%     'trial_history_mode',1,'observed_structure','binary','observed_parameter','eyetar_l','past_structure','binary','past_parameter','microstim','past_value',NaN,'past_keys',{'choice',NaN}});
%
% % trial history plots.....
% [out_comp variable_to_test counter]= monkeypsych_analyze_20140814(analysis_parameters.filelist_formatted,...
%     {'nsacc_max',10,'summary',1,'type',2,'inferential_on',0,'success',1,'display',1,'choice',1,'microstim',NaN,'boxplot_on',1,'multicomparison',1,'scatter_on',1,...
%     'trial_history_mode',1,'observed_structure','saccades','observed_parameter','tar_pos','past_structure','saccades','past_parameter','tar_pos','past_value',NaN,'past_keys',{'microstim',0,'choice',NaN}});
%
%
% % trial history plots.....
% %(994 677)
% %(4084 2131)
% DAG_delete_last_empty_trial('Y:\Data\Cornelius\20140910')
% [out_comp variable_to_test counter]= monkeypsych_analyze_20140910({'Y:\Data\Cornelius\20140910',1},{'summary',0});
%


%% missing RT plots and outputs...


% update_stimulation_parameters_20140809
% analysis_parameters.folders.dropbox            =strcat('C:', filesep,'Users', filesep,  analysis_parameters.user, filesep, 'Dropbox', filesep, 'DAG', filesep, 'microstim_behavior', filesep, analysis_parameters.monkey, '_microstim_dpz', filesep);
% analysis_parameters.folders.extended_data      =strcat(current_dir(1),':', filesep,'Data', filesep,  analysis_parameters.monkey, '_microstim_with_parameters', filesep);
% analysis_parameters.folders.original_data      =strcat(current_dir(1),':', filesep,'Data', filesep,  analysis_parameters.monkey, filesep);
% analysis_parameters.folders.output             =strcat(current_dir(1),':', filesep,'Data', filesep,  analysis_parameters.monkey, '_microstim_with_parameters', filesep);
% analysis_parameters.folders.code               =strcat(current_dir(1),':', filesep,'microstim_behavior', filesep, 'code', filesep);
%
%


end

function [filelist_complete filelist_formatted] = DAG_get_filelist_from_folder(varargin)
global analysis_parameters

if numel(varargin)==0
    folder_with_session_days=analysis_parameters.folders.extended_data;
    dates=analysis_parameters.dates;
    disp('get_filelist_from_folder is taking folder and dates from analysis_parameters')
else
    folder_with_session_days=varargin{1};
    dates=varargin{2};
    disp('get_filelist_from_folder is taking folder and dates from input')
end

dir_folder_with_session_days=dir(folder_with_session_days); % dir
if all(isnan(dates)) || all(isempty(dates)) || all(dates==0)
    all_files_in_base_path=DAG_keep_only_numeric_cell_entries({dir_folder_with_session_days.name}); % all files from the main monkey folder
else
    all_files_in_base_path=[];
    ctr=1;
    for k=1: length(dir_folder_with_session_days)
        X=str2double(dir_folder_with_session_days(k).name);
        if X==dates(1) ||  ( X<=  dates(2) && X >  dates(1))
            all_files_in_base_path{ctr}= dir_folder_with_session_days(k).name;
            ctr=ctr+1;
        end
    end
    
end

i_run=1;
for in_folders = 1:length(all_files_in_base_path)
    individual_day_folder = [folder_with_session_days all_files_in_base_path{in_folders}]; % session of interest
    d_individual_day_folder=dir(individual_day_folder); % dir
    files_inside_session_folder={d_individual_day_folder.name}'; % files inside session folders
    for number_of_files = 1:length(files_inside_session_folder) % start looping within the session folder
        if length(files_inside_session_folder{number_of_files}) > 4 && strcmp(files_inside_session_folder{number_of_files}(end-3:end),'.mat')
            filelist_complete(i_run,:)=[individual_day_folder filesep files_inside_session_folder{number_of_files}];
            filelist_formatted(i_run,:)= {filelist_complete(i_run,1:end-21) str2double(filelist_complete(i_run,end-5:end-4))};
            i_run=i_run+1;
        end
    end
    clc;
end
end

function output_cell=DAG_keep_only_numeric_cell_entries(original_cell)
logidx=false(size(original_cell,1),size(original_cell,2));
for idx=1:numel(original_cell)
    if isnan(str2double(original_cell{idx}))
        logidx(idx)=false;
    else
        logidx(idx)=true;
    end
end
output_cell=original_cell(logidx);
end

function [mastertable,out_comp]=DAG_get_trialinfo_mastertable(varargin)
global analysis_parameters
current_date            =datestr(date,'yyyymmdd');

if numel(varargin)==0
    monkey                   =analysis_parameters.monkey;
    monkeypsych_analyze      =analysis_parameters.files.monkeypsych_analyze;
    filelist_formatted       =analysis_parameters.filelist_formatted;
    working_directory        =analysis_parameters.folders.output;
    dates                    =analysis_parameters.dates;
else
    %to be added....
    if numel(varargin)>0
        monkey=varargin{1};
    end
    if numel(varargin)>1
        filelist_formatted=varargin{2};
    end
    if numel(varargin)>2
        working_directory=varargin{3};
    end
    if numel(varargin)>3
        monkeypsych_analyze=varargin{4};
    end
    if numel(varargin)>4
        dates=varargin{5};
    end
end

Sel_all={'display',0,'runs_as_batches',1};
[out_comp,~,~]= feval('monkeypsych_analyze_working',filelist_formatted,Sel_all);
% [out_comp,~,~]= monkeypsych_analyze_working(filelist_formatted,Sel_all);
trialinfo = DAG_most_recent_version(working_directory,strcat(monkey, '_trialinfo_mastertable'));

table_for_updating(1,:)={'Session','Run','total_trials','hits','hits_t_2_e_0','hits_t_2_e_4','task_type','task_effector','demanded_hand','reach_hand','x_positions','x_fixation','x_distance_to_center','y_positions','choice_percentage','microstim_percentage','microstim_start', 'microstim_state','hue', 'percentage_left_chosen_stimulated', 'percentage_left_chosen_baseline', 'RT_R_mean_stimulated', 'RT_R_std_stimulated',...
    'RT_L_mean_stimulated', 'RT_L_std_stimulated', 'RT_R_mean_baseline', 'RT_R_std_baseline', 'RT_L_mean_baseline', 'RT_L_std_baseline', 'min_time_fix_hold','max_time_fix_hold'};

if isempty(trialinfo)
    shift_rows=0;
    original_table=table_for_updating;
else
    load(trialinfo,'mastertable');
    original_table=mastertable;
    title_index=DAG_find_column_index(original_table,'Session');
    if all([original_table{2:end,title_index}] > dates(1)) % ignore_not_matching_dates
        shift_rows=0;
    else
        pot_run_files_start=find([original_table{2:end,title_index}]==dates(1));
        if ~isempty(pot_run_files_start)
            shift_rows=pot_run_files_start(1)-1;
        else
            shift_rows=size(original_table,1)-1;
        end
    end
end

for k=1:numel(out_comp)
    n_trials(k)                 = numel(out_comp{k}.binary);
    hits(k)                     = sum([out_comp{k}.binary.success]);
    hits_t_2_e_0(k)             = sum([out_comp{k}.binary.success] & [out_comp{k}.task.effector]==0 & [out_comp{k}.task.type]==2);
    hits_t_2_e_4(k)             = sum([out_comp{k}.binary.success] & [out_comp{k}.task.effector]==4 & [out_comp{k}.task.type]==2);
    type{k,:}                   = unique([out_comp{k}.task.type]);
    effector{k,:}               = unique([out_comp{k}.task.effector]);
    demanded_hand{k,:}          = unique([out_comp{k}.reaches(~isnan([out_comp{k}.reaches().demanded_hand])).demanded_hand]);
    reach_hand{k,:}             = unique([out_comp{k}.reaches(~isnan([out_comp{k}.reaches().reach_hand])).reach_hand]);
    x_positions{k,:}            = round(real(unique([out_comp{k}.saccades(~isnan(real([out_comp{k}.saccades.tar_pos]))).tar_pos, out_comp{k}.reaches(~isnan(real([out_comp{k}.reaches.tar_pos]))).tar_pos])));
    x_fixation{k,:}             = unique(round(real(unique([out_comp{k}.saccades(~isnan(real([out_comp{k}.saccades.fix_pos]))).fix_pos, out_comp{k}.reaches(~isnan(real([out_comp{k}.reaches.fix_pos]))).fix_pos])))); %%% IT GAVE A VECTOR OF TWO ZEROS, NEED TO CHECK WHY
    x_subtraction_to_center_complex{k,:}   = round(unique([[out_comp{k}.saccades.tar_pos]-[out_comp{k}.saccades.fix_pos], [out_comp{k}.reaches.tar_pos]-[out_comp{k}.reaches.fix_pos]])); 
    x_distance_to_center{k,:}   = round(x_subtraction_to_center_complex{k,:}(~isnan(x_subtraction_to_center_complex{k,:})));
    y_positions{k,:}            = round(imag(unique([out_comp{k}.saccades(~isnan(real([out_comp{k}.saccades.tar_pos]))).tar_pos, out_comp{k}.reaches(~isnan(real([out_comp{k}.reaches.tar_pos]))).tar_pos])));
    
    choice_percentage(k)        = round(nanmean([out_comp{k}.binary.choice])*100);
    microstim_percentage(k)     = round(nanmean([out_comp{k}.binary.microstim])*100);
    
    %% not necessarily correct....
    microstim_start_after_go_withnan  = unique([out_comp{k}.task([out_comp{k}.task.stim_state]==2 | ([out_comp{k}.task.stim_state]==3 & [out_comp{k}.task.type]==1) |[out_comp{k}.task.stim_state]==4 | [out_comp{k}.task.stim_state]==6 | [out_comp{k}.task.stim_state]==9).stim_start]);
    microstim_start_before_go_withnan = unique([out_comp{k}.task([out_comp{k}.binary.success] & (([out_comp{k}.task.stim_state]==3 & [out_comp{k}.task.type]~=1)| [out_comp{k}.task.stim_state]==5 | [out_comp{k}.task.stim_state]==7 | [out_comp{k}.task.stim_state]==10)).stim_to_state_end]);
    microstim_start_before_go_withnan = microstim_start_before_go_withnan(microstim_start_before_go_withnan<=0.2);
    microstim_start{k,:}              = round([sort(microstim_start_before_go_withnan(~isnan(microstim_start_before_go_withnan)).*-1) microstim_start_after_go_withnan(~isnan(microstim_start_after_go_withnan))].*1000);
    microstim_states                  = unique([out_comp{k}.task.stim_state]);
    microstim_state{k,:}              = microstim_states(~isnan(microstim_states));
    hues                                = vertcat(out_comp{k}.saccades(~isnan([out_comp{k}.saccades.tar_pos])).col_dim);
    if ~isempty(hues); hue{k,:}                            = unique(hues(1:2:end,1)); else hue{k,:}=NaN; end;
    percentage_left_chosen_stimulated(k) = round([out_comp{k}.counts.left_choice_percentage_successful_microstim]);
    percentage_left_chosen_baseline(k)   = round([out_comp{k}.counts.left_choice_percentage_successful_baseline]);
    
    RT_R_mean_stimulated(k) =  round(nanmean([out_comp{k}.saccades([out_comp{k}.binary.eyetar_r] & [out_comp{k}.binary.success] & [out_comp{k}.binary.microstim]).lat])*1000);
    RT_R_std_stimulated(k)  =  round(nanstd([out_comp{k}.saccades([out_comp{k}.binary.eyetar_r] & [out_comp{k}.binary.success] & [out_comp{k}.binary.microstim]).lat])*1000);
    RT_L_mean_stimulated(k) =  round(nanmean([out_comp{k}.saccades([out_comp{k}.binary.eyetar_l] & [out_comp{k}.binary.success] & [out_comp{k}.binary.microstim]).lat])*1000);
    RT_L_std_stimulated(k)  =  round(nanstd([out_comp{k}.saccades([out_comp{k}.binary.eyetar_l] & [out_comp{k}.binary.success] & [out_comp{k}.binary.microstim]).lat])*1000);
    
    RT_R_mean_baseline(k) =  round(nanmean([out_comp{k}.saccades([out_comp{k}.binary.eyetar_r] & [out_comp{k}.binary.success] & ~[out_comp{k}.binary.microstim]).lat])*1000);
    RT_R_std_baseline(k)  =  round(nanstd([out_comp{k}.saccades([out_comp{k}.binary.eyetar_r] & [out_comp{k}.binary.success] & ~[out_comp{k}.binary.microstim]).lat])*1000);
    RT_L_mean_baseline(k) =  round(nanmean([out_comp{k}.saccades([out_comp{k}.binary.eyetar_l] & [out_comp{k}.binary.success] & ~[out_comp{k}.binary.microstim]).lat])*1000);
    RT_L_std_baseline(k)  =  round(nanstd([out_comp{k}.saccades([out_comp{k}.binary.eyetar_l] & [out_comp{k}.binary.success] & ~[out_comp{k}.binary.microstim]).lat])*1000);
    min_time_fix_hold(k)     =  round(min([out_comp{k}.timing.fix_time_hold])*1000);
    max_time_fix_hold(k)     =  round(max([out_comp{k}.timing.fix_time_hold])*1000) + round(max([out_comp{k}.timing.fix_time_hold_var])*1000);
    
    table_for_updating(k+1+shift_rows,:)= {str2double(filelist_formatted{k,1}(end-7:end)),filelist_formatted{k,2},n_trials(k),hits(k),hits_t_2_e_0(k),hits_t_2_e_4(k),...
        type{k,:},effector{k,:},demanded_hand{k,:},reach_hand{k,:},x_positions{k,:},x_fixation{k,:},x_distance_to_center{k,:},y_positions{k,:},...
        choice_percentage(k),microstim_percentage(k),microstim_start{k,:},microstim_state{k,:},hue{k,:}...
        percentage_left_chosen_stimulated(k),percentage_left_chosen_baseline(k),...
        RT_R_mean_stimulated(k),RT_R_std_stimulated(k),RT_L_mean_stimulated(k),RT_L_std_stimulated(k),...
        RT_R_mean_baseline(k),RT_R_std_baseline(k),RT_L_mean_baseline(k),RT_L_std_baseline(k),min_time_fix_hold(k),max_time_fix_hold(k)};
    
end

raw_data_rows=[shift_rows+2:shift_rows+2+size(table_for_updating,1)];
mastertable= DAG_update_mastertable_cell(original_table,table_for_updating,raw_data_rows);
save(strcat(working_directory, monkey, '_trialinfo_mastertable_',current_date),'mastertable')
end

% function [full_path name directory date_modified]= DAG_most_recent_version(main_directory,name_fragment,varargin)
% if numel(varargin)>0
%     look_in_subfolders=1;
% else
%     look_in_subfolders=0;    
% end
% directory=main_directory;
% dir_content=dir(main_directory);
% filenames={dir_content.name};
% last_modified_dates=cellstr(datestr({dir_content.date},'yyyymmddHHMMSS'));
% %last_modified_dates=cellfun(@(x) [x(1:8) x(10:15)],last_modified_dates,'UniformOutput',false);
% isdirectory=[dir_content.isdir];
% date_modified=0;
% for k=3:numel(filenames)
%     folderflag=false;
%     if isdirectory(k) && look_in_subfolders
%         [filenames{k} subfolder last_modified_dates{k}]= DAG_most_recent_version([main_directory filesep filenames{k}],name_fragment);
%         folderflag=true;
%     end
%     if ~isempty(strfind(filenames{k},name_fragment))
%         potential_date=str2double(last_modified_dates{k});
%         if potential_date  > date_modified
%             date_modified= potential_date;
%             if folderflag                
%                 directory=subfolder;
%                 name=filenames{k};
%             else
%                 n_file_extension_chars=numel(filenames{k})-strfind(filenames{k},'.');
%                 up_to_char=numel(filenames{k})-n_file_extension_chars-1;
%                 name=filenames{k}(1:up_to_char);
%             end
%         end
%     end
% end
% if date_modified==0    
%     name=[];
%     full_path=[];
% else
%     full_path=[directory filesep name];
% end
% date_modified=num2str(date_modified);
% end

function update_stimulation_parameters(mastertable_type)
global analysis_parameters
monkey=analysis_parameters.monkey;
working_folder=analysis_parameters.folders.output;
dates=analysis_parameters.dates;
original_parameters_excel_file=strcat(working_folder, analysis_parameters.files.original_parameters);
updated_parameters_excel_file=strcat(working_folder,analysis_parameters.files.updated_parameters);

if exist(updated_parameters_excel_file, 'file') && strcmp(mastertable_type,'trialinfo')
    movefile(updated_parameters_excel_file, strcat(analysis_parameters.folders.output, 'Savecopy' ,filesep, analysis_parameters.files.updated_parameters));
elseif exist(updated_parameters_excel_file, 'file') 
    copyfile(updated_parameters_excel_file, strcat(analysis_parameters.folders.output, 'Savecopy' ,filesep, analysis_parameters.files.updated_parameters));
end

mastertable_to_update=DAG_most_recent_version(working_folder,strcat(monkey, '_', mastertable_type, '_mastertable'));
[~,parametersheets] = xlsfinfo(original_parameters_excel_file);
load(mastertable_to_update);
for k=1:numel(parametersheets)
    if isnan(str2double(parametersheets{k}))
        continue
    elseif str2double(parametersheets{k}) >= dates(1) && str2double(parametersheets{k}) <= dates(2)
        stimulationparametersfile=original_parameters_excel_file;
    else
        stimulationparametersfile=updated_parameters_excel_file;
    end
    [~,~,parameter_table] = xlsread(stimulationparametersfile,parametersheets{k});
    titles_update=[mastertable(1,:)];
    Session_column_idx=ismember(titles_update,'Session');
    Session_row_idx=[false; cell2mat([mastertable(2:end,Session_column_idx)])==str2double(parametersheets{k})];
    completed_table = DAG_update_mastertable_cell(parameter_table,[mastertable(1,:); mastertable(Session_row_idx,:)],2:sum(Session_row_idx)+1);
    
    completed_table_xls=completed_table;
    for col=1:numel(completed_table)
        if isempty(completed_table{col}) ||  all(isnan(completed_table{col})) || strcmp(completed_table{col}, 'NaN')  
            completed_table_xls{col}=[]; 
        elseif ~ischar(completed_table{col}) && numel(completed_table{col})>1
            kommas=repmat({','},1,numel(completed_table{col}));
            restructured_cell=vertcat(cellstr(num2str(completed_table{col}(:)))',kommas);
            restructured_string=[restructured_cell{:}];
            restructured_string(end)=[];
            completed_table_xls{col}=restructured_string;
        end
    end
    if size(completed_table_xls,2)>255
        completed_table_xls=completed_table_xls(:,1:255);
    end
    xlswrite(updated_parameters_excel_file,completed_table_xls,parametersheets{k});
end

clear complete_mastertable
sheet_counter=0;
complete_mastertable={};
for k=1:numel(parametersheets)
    if isempty (str2num(parametersheets{k}))
            continue
    else
        sheet_counter=sheet_counter+1;   
    end    
    [~,~,parameter_table] = xlsread(updated_parameters_excel_file,parametersheets{k});
    row_counter=size(complete_mastertable,1);
    n_rows=size(parameter_table,1)-1;
    if sheet_counter==1
        complete_mastertable=parameter_table;        
    else        
        table_for_updating=[parameter_table(1,:); cell(row_counter-1,size(parameter_table,2)); parameter_table(2:end,:)];
        complete_mastertable = DAG_update_mastertable_cell(complete_mastertable,table_for_updating,[row_counter+1:row_counter+n_rows]);
    end
    
end

if size(complete_mastertable,2)>255
   complete_mastertable=complete_mastertable(:,1:255);
   disp('mastertable has too many columns, size has been reduced to 255 columns')
end
xlswrite(updated_parameters_excel_file,complete_mastertable,'mastertable');
end

function [task] = import_parameters_from_xls_to_trial_task_per_date(individual_day_folder,date_as_string)
global analysis_parameters
xls_file=strcat(analysis_parameters.folders.output, analysis_parameters.files.original_parameters);

ii_run=2;
[~,~,raw] = xlsread(xls_file,date_as_string);
d_individual_day_folder=dir(individual_day_folder); % dir
files_inside_session_folder={d_individual_day_folder.name}'; % files inside session folders
for number_of_files = 1:length(files_inside_session_folder) % start looping within the session folder
    if length(files_inside_session_folder{number_of_files}) > 10 && strcmp(files_inside_session_folder{number_of_files}(end-3:end),'.mat')
        load([individual_day_folder filesep files_inside_session_folder{number_of_files}]);
        parameter_index=1;
        while true
            if size(raw,2) < parameter_index
                break;
            end
            current_parameter=raw{1,parameter_index};
            if isnan(current_parameter)
                break;
            else
                task.(current_parameter)=raw{ii_run,parameter_index};
                parameter_index=parameter_index+1;
            end
        end
        ii_run=ii_run+1;
        tmp_field_cell=repmat({task},numel(trial),1);
        [trial.task]=tmp_field_cell{:};
        save([individual_day_folder filesep files_inside_session_folder{number_of_files}],'task','trial','SETTINGS')
    end
end
end

function column_index=DAG_find_column_index(inputcell,title)
for m=1:size(inputcell,2)
    if strcmp(inputcell{1,m},title)
        column_index=m;
    end
end
end

function [task] = DAG_import_parameters_from_xls_to_trial_task(varargin)

%DAG_import_parameters_from_xls_to_trial_task(20140719,'Y:\Data\Curius_microstim_with_parameters','Curius_stimulation_parameters.xls')
%DAG_import_parameters_from_xls_to_trial_task([20140710 20140722],'Y:\Data\Linus_microstim_with_parameters','Linus_stimulation_parameters.xls')

%xls_file='curius_stimulation_parameters_2.ods';
%[~,~,raw] = xlsread(xls_file,num2str(days));

%base_path = 'C:\Users\Adan_Ulises\Desktop\Behavioral_analysis_audv\Pulvinar_project\setup1_microstim';
%base_path = 'C:\Users\lschneider\Desktop\monkey_data\Curius';

global analysis_parameters
if numel(varargin)==0
    dates=analysis_parameters.dates;
    folder_with_session_days=analysis_parameters.folders.extended_data;
    xls_file=strcat(analysis_parameters.folders.output, analysis_parameters.files.original_parameters);
else
    if numel(varargin)>0
        dates=varargin{1};
    end
    if numel(varargin)>1
        folder_with_session_days=varargin{2};
    end
    if numel(varargin)>2
        xls_file=varargin{3};
    end
end

if numel(dates)==1
    dates=[dates dates];
end

dir_folder_with_session_days=dir(folder_with_session_days); % dir
all_files_in_base_path=[];
ctr=1;
for k=1: length(dir_folder_with_session_days)
    X=str2num(dir_folder_with_session_days(k).name);
    if X==dates(1) |  ( X<=  dates(2) & X >  dates(1))
        all_files_in_base_path{ctr}= dir_folder_with_session_days(k).name;
        valid_dates{ctr}=X;
        ctr=ctr+1;
    end
end


for in_folders = 1:length(all_files_in_base_path)
    ii_run=2;
    if length(all_files_in_base_path{in_folders})==8 % take all the folders named as yyyymmdd from the main monkey folder
        [~,~,raw] = xlsread(xls_file,num2str(valid_dates{in_folders}));
        folder_to_look_at{1} = all_files_in_base_path{in_folders}; % start looping through session folders
        individual_day_folder = [folder_with_session_days folder_to_look_at{1}]; % session of interest
        d_individual_day_folder=dir(individual_day_folder); % dir
        files_inside_session_folder={d_individual_day_folder.name}'; % files inside session folders
        for number_of_files = 1:length(files_inside_session_folder) % start looping within the session folder
            SETTINGS=[];
            if length(files_inside_session_folder{number_of_files}) > 10 && strcmp(files_inside_session_folder{number_of_files}(end-3:end),'.mat')
                load([individual_day_folder filesep files_inside_session_folder{number_of_files}]);
                parameter_index=1;
                while true
                    if size(raw,2) < parameter_index
                        break;
                    end
                    current_parameter=raw{1,parameter_index};
                    if isnan(current_parameter)
                        break;
                    else
                        task.(current_parameter)=raw{ii_run,parameter_index};
                        parameter_index=parameter_index+1;
                    end
                end
                ii_run=ii_run+1;
                tmp_field_cell=repmat({task},numel(trial),1);
                [trial.task]=tmp_field_cell{:};
                save([individual_day_folder filesep files_inside_session_folder{number_of_files}],'task','trial','SETTINGS')
            end
        end
    end
end
end

