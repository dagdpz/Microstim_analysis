%clear all

do_from_scratch=1;
recalculate=1;
scatter_on=1;
histograms_on=1;

if do_from_scratch
    f=dir(pwd);
    f={f.name};
    n=find(strcmp(f,'set_settings.m'));
    for k=1:numel(f);
        if k~=n
            delete([pwd filesep f{k}])
        end
    end
    mkdir([pwd filesep 'scatterplots'])
    mkdir([pwd filesep 'histograms'])
    mkdir([pwd filesep 'Savecopy'])
end
global analysis_parameters
%analysis_parameters.batch_list_rearranged=0;
%analysis_parameters.evoked_definition='evoked_and_return'; %'evoked_l';
analysis_parameters.dag_drive='Y';
run set_settings;
%analysis_parameters.folders.mastertable_location;
monkey=analysis_parameters.monkey;
if do_from_scratch
    copyfile (strcat(analysis_parameters.folders.mastertable_location, analysis_parameters.files.original_parameters), strcat(analysis_parameters.folders.output, analysis_parameters.files.original_parameters));
    %current_trialinfo='Linus_trialinfo_mastertable_20140930';
    %copyfile (strcat(analysis_parameters.folders.mastertable_location, current_trialinfo, '.mat'), strcat(analysis_parameters.folders.output, current_trialinfo, '.mat'));
end

%% Specific filelist for batch processing
batch_processing            = analysis_parameters.batch_processing;
switch batch_processing
    case 0
        [~, analysis_parameters.filelist_formatted] = DAG_get_filelist_from_folder;
    case 1
        drive=analysis_parameters.dag_drive;
        analysis_parameters.filelist_formatted={};
        for k=1:numel(analysis_parameters.batches.Aiminput)
            Aiminput=analysis_parameters.batches.Aiminput{k};
            Inputsequal=analysis_parameters.batches.Inputsequal{k};
            Inputsrange=analysis_parameters.batches.Inputsrange{k};
            Inputslist=analysis_parameters.batches.Inputslist{k};
            filelist_formatted_k=get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,drive,monkey);
            analysis_parameters.filelist_formatted =  [analysis_parameters.filelist_formatted;filelist_formatted_k];
            
        end
%             Inputsequal_for_batch={'Session','Electrode_depth'};
%             analysis_parameters.filelist_formatted = DAG_get_batch_input_from_xls(analysis_parameters.filelist_formatted,Inputsequal_for_batch,0,drive,monkey);
    case 2 % if batches are related to batches of other criterions
        filelist_formatted_in_batches = DAG_arrange_baselines(analysis_parameters.batches.pool,analysis_parameters.batches.referred_to);
        analysis_parameters.filelist_formatted =vertcat(filelist_formatted_in_batches{:});
    case 3 % if batches are related to batches of other criterions, only used Sessions/Current_strengths
        filelist_formatted_in_batches = DAG_arrange_batches(analysis_parameters.batches.pool,analysis_parameters.batches.referred_to,...
            [analysis_parameters.folders.mastertable_location analysis_parameters.files.updated_parameters],'Session','Current_strength');
        analysis_parameters.filelist_formatted =vertcat(filelist_formatted_in_batches{:});
end



if do_from_scratch
    if batch_processing
        %if batch_processing==1
        %reduce_stimulation_parameters(vertcat(analysis_parameters.filelist_formatted{:}),strcat(analysis_parameters.folders.output, analysis_parameters.files.original_parameters))
        %else
        reduce_stimulation_parameters(analysis_parameters.filelist_formatted,strcat(analysis_parameters.folders.output, analysis_parameters.files.original_parameters))
            
        %end
        movefile (strcat(analysis_parameters.folders.output, analysis_parameters.files.original_parameters), strcat(analysis_parameters.folders.output, 'Savecopy', filesep, analysis_parameters.files.original_parameters));
        movefile (strcat(analysis_parameters.folders.output, 'reduced_parameters.xls'), strcat(analysis_parameters.folders.output, analysis_parameters.files.original_parameters));
    end
    if analysis_parameters.batch_list_rearranged
        %Inputsequal={'Session','Electrode_depth','Current_strength'};
        Inputsequal={'Session','Electrode_depth'};
        
        
%         filelist_formatted_in_batches = DAG_get_batch_input_from_xls(analysis_parameters.filelist_formatted,Inputsequal,0,analysis_parameters.dag_drive,monkey);
%         analysis_parameters.filelist_formatted =vertcat(filelist_formatted_in_batches{:});
        
        analysis_parameters.filelist_formatted = DAG_get_batch_input_from_xls(analysis_parameters.filelist_formatted,Inputsequal,0,analysis_parameters.dag_drive,monkey);
        
    end
        
    evoked_saccade_preparation
    evoked_saccade_calculation
    
    if batch_processing
        DAG_get_trialinfo_mastertable
    end
    
    get_evoked_mastertable
    update_stimulation_parameters('trialinfo')
    update_stimulation_parameters('evoked')
%     updated_parameters_excel_file=strcat(analysis_parameters.folders.output,analysis_parameters.files.updated_parameters);
%     [~,~,parameter_table] = xlsread(updated_parameters_excel_file,'mastertable');
%     evoked_fix_idx=DAG_find_column_index(parameter_table(1,:),'dP_evoked_fix');
%     evoked_dir_idx=DAG_find_column_index(parameter_table(1,:),'dP_evoked_fix');
%     valid_indexes=cellfun(@(x) ~all(isnan(x)),parameter_table(:,evoked_idx));
%     valid_parameter_table=parameter_table(valid_indexes,:);
%     xlswrite(updated_parameters_excel_file,valid_parameter_table,'mastertable');
end



if scatter_on
    if recalculate && ~ do_from_scratch
        evoked_saccade_calculation
        get_evoked_mastertable
        update_stimulation_parameters('trialinfo')
        update_stimulation_parameters('evoked')
    end
    summary_scatterplot_frenzy;
end

if histograms_on
    if recalculate && batch_processing
        load(DAG_most_recent_version(pwd,'evoked_output'));
        out_fix= concatinate_structure_cell_fields(out_fix,'task','states','binary','saccades','selected');
        out_dir= concatinate_structure_cell_fields(out_dir,'task','states','binary','saccades');
        out_mem= concatinate_structure_cell_fields(out_mem,'task','states','binary','saccades');
        evoked_saccade_calculation(out_fix,out_dir,out_mem);
    end
    % evoked_saccade_plotting;
    evoked_saccade_plotting_RL;
end
