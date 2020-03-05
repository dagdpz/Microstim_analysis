function evoked_saccade_preparation


run set_settings %?
global analysis_parameters

monkey                      = analysis_parameters.monkey;
batch_processing            = analysis_parameters.batch_processing;
Selection_evoked            = analysis_parameters.Selection_evoked;

dates                       = analysis_parameters.dates;
Evo_amp_TH                  = analysis_parameters.Evo_amp_TH;
Evo_RT_TH                   = analysis_parameters.Evo_RT_TH;

current_date                = analysis_parameters.current_date;

%% Specific analysis_parameters for Batch or single run processing
if batch_processing
%     Aiminput=analysis_parameters.batches.Aiminput;
%     Inputsequal=analysis_parameters.batches.Inputsequal;
%     Inputsrange=analysis_parameters.batches.Inputsrange;
%     drive=analysis_parameters.dag_drive;
%     analysis_parameters.filelist_formatted =
%     DAG_get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,drive,monkey);
%     %% change of inputs!! not drive, but path
    Batch_selection=analysis_parameters.batches.additional_selection;  %Additional keys
    batch_type=analysis_parameters.batches.type; 
    if analysis_parameters.batch_list_rearranged
        Currentselection=[Selection_evoked Batch_selection {'type',batch_type}];
    else
        Currentselection=[Selection_evoked Batch_selection {'type',batch_type} {'runs_as_batches',1}];
    end
    run_file_start=1;
else
    [~, analysis_parameters.filelist_formatted] = DAG_get_filelist_from_folder;
    load(DAG_most_recent_version(pwd,strcat(monkey, '_trialinfo_mastertable')));
    Currentselection        =   [Selection_evoked {'runs_as_batches',1}];
    title_index             =   DAG_find_column_index(mastertable,'Session');
    pot_run_files_start     =   find([mastertable{2:end,title_index}]==dates(1));
    run_file_start          =   pot_run_files_start(1);
end

filelist_formatted          = analysis_parameters.filelist_formatted;
analyzed_output=DAG_most_recent_version(pwd,'evoked_output');

if isempty(analyzed_output)
    out_fix={};
    out_dir={};
    out_mem={};
else
    load(analyzed_output)
end
Sel_all_fix=[Currentselection, {'type',1,'evoked_stimlocked_latency_min',Evo_RT_TH(1),'evoked_stimlocked_latency_max',Evo_RT_TH(2),'evoked_amplitude_min',Evo_amp_TH(1),'consecutive_saccade',0}];
Sel_all_dir=[Currentselection, {'type',2,'evoked_stimlocked_latency_min',Evo_RT_TH(1),'evoked_stimlocked_latency_max',Evo_RT_TH(2),'evoked_amplitude_min',Evo_amp_TH(1),'consecutive_saccade',0}];
Sel_all_mem=[Currentselection, {'type',3,'evoked_stimlocked_latency_min',Evo_RT_TH(1),'evoked_stimlocked_latency_max',Evo_RT_TH(2),'evoked_amplitude_min',Evo_amp_TH(1),'consecutive_saccade',0}];


    
if analysis_parameters.batch_list_rearranged    
    sel_array_fix=repmat({Sel_all_fix},1,size(filelist_formatted,2));
    sel_array_dir=repmat({Sel_all_dir},1,size(filelist_formatted,2));
    sel_array_mem=repmat({Sel_all_mem},1,size(filelist_formatted,2));
    Analyze_input_fix=[filelist_formatted;sel_array_fix];
    Analyze_input_dir=[filelist_formatted;sel_array_dir];
    Analyze_input_mem=[filelist_formatted;sel_array_mem];
    
    [out_new_fix,~,~]= monkeypsych_analyze_working(Analyze_input_fix{:});
    [out_new_dir,~,~]= monkeypsych_analyze_working(Analyze_input_dir{:});
    [out_new_mem,~,~]= monkeypsych_analyze_working(Analyze_input_mem{:});
else
    [out_new_fix,~,~]= monkeypsych_analyze_working(filelist_formatted,Sel_all_fix);
    [out_new_dir,~,~]= monkeypsych_analyze_working(filelist_formatted,Sel_all_dir);
    [out_new_mem,~,~]= monkeypsych_analyze_working(filelist_formatted,Sel_all_mem);
end
out_fix(run_file_start:numel(out_new_fix))=out_new_fix;
out_dir(run_file_start:numel(out_new_dir))=out_new_dir;
out_mem(run_file_start:numel(out_new_mem))=out_new_mem;


save([monkey, '_evoked_output_', current_date],'out_dir','out_fix','out_mem');
end