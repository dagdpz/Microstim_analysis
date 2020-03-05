function saccade_histograms=evoked_saccade_calculation(varargin)

run set_settings
global analysis_parameters

monkey                      = analysis_parameters.monkey;
working_data_folder         = analysis_parameters.folders.extended_data;
dates                       = analysis_parameters.dates;
all_locations               = analysis_parameters.all_locations;
location_extensions         = analysis_parameters.location_extensions;
one_saccade_per_trial       = analysis_parameters.one_saccade_per_trial;
amplitudes                  = analysis_parameters.amplitudes;
RT_inc                      = analysis_parameters.RT_inc;
far_definition              = analysis_parameters.far_definition;
current_date                = analysis_parameters.current_date;
batch_processing            = analysis_parameters.batch_processing;
%excentricity_calculated     = analysis_parameters.batches.excentricity;

Evo_amp_TH                  = analysis_parameters.Evo_amp_TH;
Evo_RT_TH                   = analysis_parameters.Evo_RT_TH;
Del_RT_TH                   = analysis_parameters.Del_RT_TH;
Del_amp_TH                  = analysis_parameters.Del_amp_TH;
      
STATE=load_state;
valid_batches=[];         % This variable will keep track of the files that have been analyzed
if nargin>0    
    out_fix= varargin{1};
    out_dir= varargin{2};
    out_mem= varargin{3};
else
    load(DAG_most_recent_version(pwd,'evoked_output'));
end

%% Specific analysis_parameters for Batch or single run processing
run_file_start=1; 
run_file_end=numel(out_dir); 

if batch_processing
    type_modes={analysis_parameters.batches.type_mode};
    batch_starts=repmat({[analysis_parameters.batches.stim_starts]},run_file_end,1);
    batch_states=repmat({analysis_parameters.batches.stim_states},run_file_end,1);
    batch_states_for_starts=repmat({analysis_parameters.batches.stim_states_for_starts},run_file_end,1);
    
else    
    % Load existing mastertables, get relevant run information from trialinfo mastertable 
    type_modes={'fix','dir','mem'};
    load(DAG_most_recent_version(pwd,'trialinfo_mastertable'));
    [~,~,stimstart,stimstate]=get_filelist_from_matfile(working_data_folder,mastertable);    
    % Date selection
    if ~isnan(dates)
        title_index=DAG_find_column_index(mastertable,'Session');
        pot_run_files_start=find([mastertable{2:end,title_index}]==dates(1));
        pot_run_files_end=find([mastertable{2:end,title_index}]==dates(2));
        run_file_start=pot_run_files_start(1);
        if isempty(pot_run_files_end)
            run_file_end=numel(stimstart);
        else         
            run_file_end=pot_run_files_end(end);   
        end
    end
end

for m=1:numel(type_modes)
    type_mode=type_modes{m};
for p = run_file_start:run_file_end;
    if batch_processing
        current_trial_stim_starts   = batch_starts{p};
        current_stim_states         = batch_states{p};
        current_states_for_starts   = batch_states_for_starts{p};
        stim_windows                = [current_states_for_starts',current_trial_stim_starts'/1000,];
    else
        % Assigning Input and Stimulation window parameters for each run
        current_trial_stim_starts   = stimstart{p};
        current_stim_states         = stimstate{p};
        clear stim_windows
        switch type_mode
            case 'mem'               
                stim_windows(:,1)=[6,7,7,7,7,9];
                stim_windows(:,2)=[0,0,0.5,1,1.5,0];
            case 'del'
                
            case 'dir'
                stim_windows(:,1)=[3,3,3,3,3,4,4,4,4,4,4,4,4,4,4];
                stim_windows(:,2)=[-0.2,-0.16,-0.12,-0.08,-0.04,0,0.04,0.08,0.1,0.11,0.12,0.13,0.14,0.15,0.16];
            case 'fix'                
                stim_windows(:,1)=[3,3,3];
                stim_windows(:,2)=[0,0.5,1];
        end
    end  
    
    switch type_mode
        % Switching analysis_parameters for different run types:
        % Sel_all defines keys for monkey_psyche_analyse, taking into account different states to look for saccades and different thresholds for different types
        % Potential stimulation windows are defined for each loop: stim_windows(:,1)= microstim_states, stim_windows(:,2)= microstim_start
        case 'mem'
            out=out_mem{p};
            condition_tags={'L_CH_S','L_IN_S','R_CH_S','R_IN_S','L_CH_B','L_IN_B','R_CH_B','R_IN_B'};
            RT_bins=0:RT_inc:0.4;
            trial_type=3;
            n_locations=numel(all_locations);            
        case 'del'
            
        case 'dir'
            out=out_dir{p};
            condition_tags={'L_CH_S','L_IN_S','R_CH_S','R_IN_S','L_CH_B','L_IN_B','R_CH_B','R_IN_B'};
            RT_bins=0:RT_inc:0.4;
            trial_type=2;
            n_locations=numel(all_locations);            
        case 'fix'
            out=out_fix{p};
            condition_tags={'stim','base'};
            RT_bins=0:RT_inc:1;
            trial_type=1;
            n_locations=1;
    end
    if isempty(out.task)
        continue
    end
    
    for target_locations=1:n_locations % looping for all, close and far targets
        far_targets         =all_locations(target_locations);
        excentricity        =location_extensions{target_locations};
        
        %Can be moved outside this loop? 
        current_stim_starts =current_trial_stim_starts./1000;
        
        %This part is for differenciating stimulation runs and run baselines -> Can be simplified?
        if numel(current_stim_states)>0
            max_states=numel(current_stim_states);
        else
            max_states=1;
            current_stim_states=NaN;
        end
        if numel(current_stim_starts)>0
            max_starts=numel(current_stim_starts);
        else
            max_starts=1;
            current_stim_starts=NaN;
        end
        
        num_stim_window=0; % Counter for current stimulation window (state/start)
        summ_max_starts=0; % Counter for valid stimulation windows in this run (used for summeries later)
        for num_stimstate=1:max_states % looping for microstim_states
            summ_max_starts=summ_max_starts+max_starts;
            for num_stimstart=1:max_starts % looping for microstim_starts                
                if current_stim_starts(num_stimstart)>=0 % positive microstimulation start
                    relevant_start='stim_start';
                    stimstart_abs=current_stim_starts(num_stimstart);
                    % break, if current_stim_start is not found in this run
                    if ~any([out.binary.microstim] & [out.task.stim_state]==current_stim_states(num_stimstate) & ([out.task.(relevant_start)]==stimstart_abs)) 
                        summ_max_starts=summ_max_starts-1; 
                        continue
                    end
                elseif  current_stim_starts(num_stimstart)<0 % negative microstimulation start
                    relevant_start='stim_to_state_end';
                    stimstart_abs=current_stim_starts(num_stimstart)*-1;
                    % break, if current_stim_start is not found in this run
                    if ismember(current_stim_states(num_stimstate), STATE.aquisition) || ~any([out.binary.microstim] & [out.task.stim_state]==current_stim_states(num_stimstate) & ([out.task.(relevant_start)]==stimstart_abs))
                        summ_max_starts=summ_max_starts-1;
                        continue
                    end
                else    % run baseline
                    relevant_start='stim_start';
                    stimstart_abs=0;
                end
                
                %% MAIN CALCULATION FUNCTION !!!
                [n_trials,RTs_histograms_all, Dir_all, Amp_all, RTs_all, Indexes_all, RTs_histograms_L, RTs_histograms_R,RTs_per_amp_L, RTs_per_amp_R,amplitudes_L,amplitudes_R,P, RTs_per_amp_L_with_return, amplitudes_L_with_return]=...
                    calc_saccade_hist(out,amplitudes,RT_bins,condition_tags,current_stim_states(num_stimstate),stimstart_abs,relevant_start,trial_type,far_targets,far_definition,one_saccade_per_trial,Evo_amp_TH,Del_amp_TH);
                
                % Assign number of evoked saccades for each condition for text in plots
                f_names=fieldnames(P.evoked_l);
                for t= 1:numel(f_names)
                    n_evoked_l.(f_names{t})               =P.evoked_l.(f_names{t})*n_trials.(f_names{t});
                    n_evoked_r.(f_names{t})               =P.evoked_r.(f_names{t})*n_trials.(f_names{t});
                    n_return_l.(f_names{t})               =P.return_l.(f_names{t})*n_trials.(f_names{t});
                    n_return_r.(f_names{t})               =P.return_r.(f_names{t})*n_trials.(f_names{t});
                    n_evoked_and_return.(f_names{t})    =P.evoked_and_return.(f_names{t})*n_trials.(f_names{t});
                end                
                f_names=fieldnames(P.delayed);
                for t= 1:numel(f_names)
                    n_delayed.(f_names{t})               =P.delayed.(f_names{t})*n_trials.(f_names{t});
                    %n_delayed_r.(f_names{t})               =P.delayed_r.(f_names{t})*n_trials.(f_names{t});
                end    
                
                % Temporary output storage in saccade_histograms structure
                % NOTE: num_stim_window only increases, if the loop has not been broken before, which means that there is information for the current stimulation start and stimulation state
                num_stim_window=num_stim_window+1;
                type_excentricities=[type_mode excentricity]; %[type_mode state_prefix];
                
                saccade_histograms(p).(type_excentricities).stim_windows                             =stim_windows;
                saccade_histograms(p).(type_excentricities).excentricity                             =excentricity;
                saccade_histograms(p).(type_excentricities).amplitudes                               =amplitudes;
                saccade_histograms(p).(type_excentricities).RT_bins                                  =RT_bins;
                saccade_histograms(p).(type_excentricities).n_trials(num_stim_window)                =n_trials;
                saccade_histograms(p).(type_excentricities).RTs_histograms_all(num_stim_window)      =RTs_histograms_all;
                saccade_histograms(p).(type_excentricities).RTs_histograms_L(num_stim_window)        =RTs_histograms_L;
                saccade_histograms(p).(type_excentricities).RTs_histograms_R(num_stim_window)        =RTs_histograms_R;
                saccade_histograms(p).(type_excentricities).stimulation_onset(num_stim_window)       =current_stim_starts(num_stimstart);
                saccade_histograms(p).(type_excentricities).stimulation_state(num_stim_window)       =current_stim_states(num_stimstate);
                saccade_histograms(p).(type_excentricities).RTs_per_amp_L(num_stim_window)           =RTs_per_amp_L;
                saccade_histograms(p).(type_excentricities).RTs_per_amp_L_with_return(num_stim_window)           =RTs_per_amp_L_with_return;
                saccade_histograms(p).(type_excentricities).RTs_per_amp_R(num_stim_window)           =RTs_per_amp_R;
                saccade_histograms(p).(type_excentricities).amplitudes_L(num_stim_window)            =amplitudes_L;
                saccade_histograms(p).(type_excentricities).amplitudes_L_with_return(num_stim_window)            =amplitudes_L_with_return;
                saccade_histograms(p).(type_excentricities).amplitudes_R(num_stim_window)            =amplitudes_R;
                saccade_histograms(p).(type_excentricities).P.evoked_l(num_stim_window)              =P.evoked_l;
                saccade_histograms(p).(type_excentricities).P.evoked_r(num_stim_window)              =P.evoked_r;
                saccade_histograms(p).(type_excentricities).P.return_l(num_stim_window)              =P.return_l;
                saccade_histograms(p).(type_excentricities).P.return_r(num_stim_window)              =P.return_r;
                saccade_histograms(p).(type_excentricities).P.evoked_and_return(num_stim_window)     =P.evoked_and_return;                 
                saccade_histograms(p).(type_excentricities).n.evoked_l(num_stim_window)              =n_evoked_l;                     
                saccade_histograms(p).(type_excentricities).n.evoked_r(num_stim_window)              =n_evoked_r;                   
                saccade_histograms(p).(type_excentricities).n.return_l(num_stim_window)              =n_return_l;                   
                saccade_histograms(p).(type_excentricities).n.return_r(num_stim_window)              =n_return_r;  
                
                saccade_histograms(p).(type_excentricities).n.delayed(num_stim_window)              =n_delayed;   
                %saccade_histograms(p).(type_excentricities).n.delayed_r(num_stim_window)              =n_delayed_r; 
                
                saccade_histograms(p).(type_excentricities).n.evoked_and_return(num_stim_window)     =n_evoked_and_return;
                saccade_histograms(p).(type_excentricities).Dir_all{num_stim_window}                 =Dir_all;
                saccade_histograms(p).(type_excentricities).Amp_all{num_stim_window}                 =Amp_all;                
                saccade_histograms(p).(type_excentricities).RTs_all{num_stim_window}                 =RTs_all;
                saccade_histograms(p).(type_excentricities).Indexes_all{num_stim_window}             =Indexes_all;
                
                
                %%% to be changed... according to keys used in
                %%% monkeypsych_analyze
                saccade_histograms(p).(type_excentricities).Evo_RT_TH                                =Evo_RT_TH;
                saccade_histograms(p).(type_excentricities).Evo_amp_TH                               =Evo_amp_TH;
                saccade_histograms(p).(type_excentricities).Del_RT_TH                                =Del_RT_TH;
                saccade_histograms(p).(type_excentricities).Del_amp_TH                               =Del_amp_TH;
            end
        end
        logidx_evoked_RTs  = RT_bins>Evo_RT_TH(1)               & RT_bins<=Evo_RT_TH(2);
        logidx_delayed_RTs = RT_bins>=Del_RT_TH(1)               & RT_bins<=Del_RT_TH(2);
        logidx_evoked_amps = [amplitudes>=Evo_amp_TH(1)         & amplitudes<Evo_amp_TH(2), false];
        logidx_delayed_amps = [amplitudes>=Del_amp_TH(1)         & amplitudes<Del_amp_TH(2), false];
        used_del_amps=find(logidx_delayed_amps);
        %% summation part
        if summ_max_starts~=0
            condition_fieldnames=fieldnames(saccade_histograms(p).(type_excentricities).n_trials);            
            % Loop for all conditions
            for n_fname=1:numel(condition_fieldnames)
                saccade_histograms(p).(type_excentricities).summed_detection_probability_histogram.(condition_fieldnames{n_fname})   =0;
                saccade_histograms(p).(type_excentricities).summed_direction_probability_histogram.(condition_fieldnames{n_fname})   =0;
                saccade_histograms(p).(type_excentricities).summed_RTs_histograms_L.(condition_fieldnames{n_fname})                  =0;
                saccade_histograms(p).(type_excentricities).summed_RTs_histograms_R.(condition_fieldnames{n_fname})                  =0;
                saccade_histograms(p).(type_excentricities).summed_dir_all.(condition_fieldnames{n_fname})                           =[];
                saccade_histograms(p).(type_excentricities).summed_amp_all.(condition_fieldnames{n_fname})                           =[];
                saccade_histograms(p).(type_excentricities).summed_RTs_all.(condition_fieldnames{n_fname})                           =[];
                saccade_histograms(p).(type_excentricities).summed_n.trials.(condition_fieldnames{n_fname})                          =0;
                saccade_histograms(p).(type_excentricities).summed_n.evoked_l.(condition_fieldnames{n_fname})                          =0;
                saccade_histograms(p).(type_excentricities).summed_n.evoked_r.(condition_fieldnames{n_fname})                          =0;
                saccade_histograms(p).(type_excentricities).summed_n.return_l.(condition_fieldnames{n_fname})                          =0;
                saccade_histograms(p).(type_excentricities).summed_n.return_r.(condition_fieldnames{n_fname})                          =0;
                saccade_histograms(p).(type_excentricities).summed_n.evoked_and_return.(condition_fieldnames{n_fname})                          =0;
                
                % Loop for all conditions
                for n_stim_start=1:summ_max_starts
                    saccade_histograms(p).(type_excentricities).summed_n.trials.(condition_fieldnames{n_fname})       =  saccade_histograms(p).(type_excentricities).summed_n.trials.(condition_fieldnames{n_fname}) +...
                        saccade_histograms(p).(type_excentricities).n_trials(n_stim_start).(condition_fieldnames{n_fname});
                    if n_fname>numel(condition_fieldnames)/2
                        saccade_histograms(p).(type_excentricities).summed_n.trials.(condition_fieldnames{n_fname})   =  saccade_histograms(p).(type_excentricities).n_trials(n_stim_start).(condition_fieldnames{n_fname});
                    end
                end
                
                for n_stim_start=1:summ_max_starts%
                    % Normalization factor for summary (=number of trials with current stim window divided by total trials, seperately for each condition)
                    if saccade_histograms(p).(type_excentricities).summed_n.trials.(condition_fieldnames{n_fname})==0 || isnan(saccade_histograms(p).(type_excentricities).summed_n.trials.(condition_fieldnames{n_fname}))
                        current_total_trials=1;
                    else
                        current_total_trials=saccade_histograms(p).(type_excentricities).summed_n.trials.(condition_fieldnames{n_fname});
                    end
                    current_normalization=saccade_histograms(p).(type_excentricities).n_trials(n_stim_start).(condition_fieldnames{n_fname})./current_total_trials; %!!!!
                    % Baseline normalization: dived by number of stimulation windows
                    if n_fname>numel(condition_fieldnames)/2
                        current_normalization=current_normalization./summ_max_starts;
                    end
                    
                    %summing up with normalization
                    saccade_histograms(p).(type_excentricities).summed_detection_probability_histogram.(condition_fieldnames{n_fname})   =  saccade_histograms(p).(type_excentricities).summed_detection_probability_histogram.(condition_fieldnames{n_fname}) ...
                        + saccade_histograms(p).(type_excentricities).RTs_histograms_all(n_stim_start).(condition_fieldnames{n_fname}).*current_normalization;
                    saccade_histograms(p).(type_excentricities).summed_RTs_histograms_L.(condition_fieldnames{n_fname})                  =  saccade_histograms(p).(type_excentricities).summed_RTs_histograms_L.(condition_fieldnames{n_fname}) ...
                        + saccade_histograms(p).(type_excentricities).RTs_histograms_L(n_stim_start).(condition_fieldnames{n_fname}).*current_normalization;
                    saccade_histograms(p).(type_excentricities).summed_RTs_histograms_R.(condition_fieldnames{n_fname})                  =  saccade_histograms(p).(type_excentricities).summed_RTs_histograms_R.(condition_fieldnames{n_fname}) ...
                        + saccade_histograms(p).(type_excentricities).RTs_histograms_R(n_stim_start).(condition_fieldnames{n_fname}).*current_normalization;
                    saccade_histograms(p).(type_excentricities).summed_n.evoked_l.(condition_fieldnames{n_fname})                  =  saccade_histograms(p).(type_excentricities).summed_n.evoked_l.(condition_fieldnames{n_fname}) ...
                        + saccade_histograms(p).(type_excentricities).P.evoked_l(n_stim_start).(condition_fieldnames{n_fname}).*saccade_histograms(p).(type_excentricities).n_trials(n_stim_start).(condition_fieldnames{n_fname});
                    saccade_histograms(p).(type_excentricities).summed_n.evoked_r.(condition_fieldnames{n_fname})                  =  saccade_histograms(p).(type_excentricities).summed_n.evoked_r.(condition_fieldnames{n_fname}) ...
                        + saccade_histograms(p).(type_excentricities).P.evoked_r(n_stim_start).(condition_fieldnames{n_fname}).*saccade_histograms(p).(type_excentricities).n_trials(n_stim_start).(condition_fieldnames{n_fname});
                    saccade_histograms(p).(type_excentricities).summed_n.return_l.(condition_fieldnames{n_fname})                  =  saccade_histograms(p).(type_excentricities).summed_n.return_l.(condition_fieldnames{n_fname}) ...
                        + saccade_histograms(p).(type_excentricities).P.return_l(n_stim_start).(condition_fieldnames{n_fname}).*saccade_histograms(p).(type_excentricities).n_trials(n_stim_start).(condition_fieldnames{n_fname});
                    saccade_histograms(p).(type_excentricities).summed_n.return_r.(condition_fieldnames{n_fname})                  =  saccade_histograms(p).(type_excentricities).summed_n.return_r.(condition_fieldnames{n_fname}) ...
                        + saccade_histograms(p).(type_excentricities).P.return_r(n_stim_start).(condition_fieldnames{n_fname}).*saccade_histograms(p).(type_excentricities).n_trials(n_stim_start).(condition_fieldnames{n_fname});
                    saccade_histograms(p).(type_excentricities).summed_n.evoked_and_return.(condition_fieldnames{n_fname})              =  saccade_histograms(p).(type_excentricities).summed_n.evoked_and_return.(condition_fieldnames{n_fname}) ...
                        + saccade_histograms(p).(type_excentricities).P.evoked_and_return(n_stim_start).(condition_fieldnames{n_fname}).*saccade_histograms(p).(type_excentricities).n_trials(n_stim_start).(condition_fieldnames{n_fname});
                   saccade_histograms(p).(type_excentricities).summed_dir_all.(condition_fieldnames{n_fname})                          =[saccade_histograms(p).(type_excentricities).summed_dir_all.(condition_fieldnames{n_fname})...
                        saccade_histograms(p).(type_excentricities).Dir_all{n_stim_start}.(condition_fieldnames{n_fname})];
                    saccade_histograms(p).(type_excentricities).summed_amp_all.(condition_fieldnames{n_fname})                          =[saccade_histograms(p).(type_excentricities).summed_amp_all.(condition_fieldnames{n_fname})...
                        saccade_histograms(p).(type_excentricities).Amp_all{n_stim_start}.(condition_fieldnames{n_fname})];
                     saccade_histograms(p).(type_excentricities).summed_RTs_all.(condition_fieldnames{n_fname})                          =[saccade_histograms(p).(type_excentricities).summed_RTs_all.(condition_fieldnames{n_fname})...
                        saccade_histograms(p).(type_excentricities).RTs_all{n_stim_start}.(condition_fieldnames{n_fname})];
                    
                    saccade_histograms(p).(type_excentricities).p_evoked_per_amp_L.(condition_fieldnames{n_fname})                       = sum(saccade_histograms(p).(type_excentricities).summed_RTs_histograms_L.(condition_fieldnames{n_fname})(logidx_evoked_amps,logidx_evoked_RTs),2);
                    temp_max                                        =amplitudes([saccade_histograms(p).(type_excentricities).p_evoked_per_amp_L.(condition_fieldnames{n_fname})==max(saccade_histograms(p).(type_excentricities).p_evoked_per_amp_L.(condition_fieldnames{n_fname})); false]);
                    saccade_histograms(p).(type_excentricities).amp_evoked_L .(condition_fieldnames{n_fname})    =temp_max(1);
                    
                    saccade_histograms(p).(type_excentricities).p_delayed_per_amp_L.(condition_fieldnames{n_fname}) = sum(saccade_histograms(p).(type_excentricities).summed_RTs_histograms_L.(condition_fieldnames{n_fname})(used_del_amps(1):end,logidx_delayed_RTs),2);
                    saccade_histograms(p).(type_excentricities).p_delayed_per_amp_R.(condition_fieldnames{n_fname}) = sum(saccade_histograms(p).(type_excentricities).summed_RTs_histograms_R.(condition_fieldnames{n_fname})(used_del_amps(1):end,logidx_delayed_RTs),2);
                    saccade_histograms(p).(type_excentricities).p_delayed_per_amp.(condition_fieldnames{n_fname})   = saccade_histograms(p).(type_excentricities).p_delayed_per_amp_R.(condition_fieldnames{n_fname})+saccade_histograms(p).(type_excentricities).p_delayed_per_amp_L.(condition_fieldnames{n_fname});
                    end
            end
            valid_batches=[valid_batches p];            
        end
    end
end
end
valid_batches=unique(valid_batches);
saccade_histograms(1).valid_batches=valid_batches;

%% Save output
if batch_processing
    save([monkey, '_histograms_batched_', current_date],'saccade_histograms');
else
    save([monkey, '_saccade_histograms_', current_date],'saccade_histograms');
end
end

%% calculation function
function [n_trials,RTs_histograms_all, Dir_all, Amp_all, RTs_all, Indexes_all, RTs_histograms_L, RTs_histograms_R, RTs_per_amp_L, RTs_per_amp_R,amplitudes_L,amplitudes_R,P,RTs_per_amp_L_with_return,amplitudes_L_with_return]=...
    calc_saccade_hist(out,amplitudes,Bins,condition_tags,stimstate,stimstart_abs,relevant_start,type,far_targets,far_definition,one_saccade_per_trial,Evo_amp_TH,Del_amp_TH)

n_trials_total=numel(out.saccades);

% Defining logical indexes for the various conditions (excentricities, Left/right, Choice/instructed, stim/baseline)
if far_targets==1
    idx.target_distance=abs(real([out.saccades.tar_pos]))>=far_definition;
elseif far_targets==0
    idx.target_distance=abs(real([out.saccades.tar_pos]))<=far_definition;
else
    idx.target_distance=ones(1,numel([out.saccades.tar_pos]));
end

idx.stim=[out.binary.microstim] & [out.task.stim_state]==stimstate & ([out.task.(relevant_start)]==stimstart_abs);
idx.base=~[out.binary.microstim];

if type~=1
    idx.L=real([out.saccades.tar_pos]-[out.saccades.fix_pos])<0;
    idx.R=real([out.saccades.tar_pos]-[out.saccades.fix_pos])>0;
    idx.CH=[out.binary.choice];
    idx.IN=~[out.binary.choice];
    
    idx.L_CH_S=idx.L & idx.CH & idx.stim;
    idx.L_IN_S=idx.L & idx.IN & idx.stim;
    idx.R_CH_S=idx.R & idx.CH & idx.stim;
    idx.R_IN_S=idx.R & idx.IN & idx.stim;
    
    idx.L_CH_B=idx.L & idx.CH & idx.base;
    idx.L_IN_B=idx.L & idx.IN & idx.base;
    idx.R_CH_B=idx.R & idx.CH & idx.base;
    idx.R_IN_B=idx.R & idx.IN & idx.base;    
end

% sac=[out.saccades];
% idx_evoked             =[sac.evoked]==1;

%         n_evoked_and_return=0;
%         n_base_return=0;
%         n_evoked=0;
%         n_return=0;
for c=1:numel(condition_tags)
    % Loop across all conditions (c)
    Dir_all.(condition_tags{c})         =[];
    Amp_all.(condition_tags{c})         =[];
    RTs_all.(condition_tags{c})         =[];
    Indexes_all.L_with_ret.(condition_tags{c})         =[];
    Indexes_all.R_with_ret.(condition_tags{c})         =[];
    
%     % Calculate evoked probability for current condition
%     idx_e=idx.(condition_tags{c}) & idx_evoked;
%     idx_ne=idx.(condition_tags{c}) & ~idx_evoked;
%     if ~(sum(idx_e)==0 && sum(idx_ne)==0)
%         P_evoked.(condition_tags{c})= sum(idx_e)/(sum(idx_e)+sum(idx_ne));
%     else
%         P_evoked.(condition_tags{c})=0 ;
%     end    
    n_trials.(condition_tags{c})=sum([idx.(condition_tags{c})]);
    for m=1:numel(amplitudes)-1
        % Loop across all amplitudes (m)
        RTs_per_amp_L.(condition_tags{c}){m}=[];
        RTs_per_amp_R.(condition_tags{c}){m}=[];        
        amplitudes_L.(condition_tags{c}){m}=[];
        amplitudes_R.(condition_tags{c}){m}=[];     
        RTs_per_amp_L_with_return.(condition_tags{c}){m}=[];
        amplitudes_L_with_return.(condition_tags{c}){m}=[];
        lat_obs=NaN(n_trials_total,40);
        directions_obs=NaN(n_trials_total,40);
        %Dir_per_amplitude.(condition_tags{c}){m}=[];
        %Dir_evoked_per_amplitude.(condition_tags{c}){m}=[];
        for k=1:n_trials_total
            if idx.(condition_tags{c})(k)
                if one_saccade_per_trial % for looking at the "voluntary" saccades, as default this is set to zero
                    lat_obs(k,:)= [out.saccades(k).lat]; 
                    directions_obs(k,:)=[out.saccades(k).rndpos-out.saccades(k).startpos];
                    amplitudes_obs=abs(directions_obs(k,:)).*sign(real(directions_obs(k,:)));
                else
                    
                    idx_1ao=out.saccades(k).ini_all>out.states(k).start_1ao;
                    idx_obs=out.saccades(k).ini_all>out.states(k).start_obs;
                    idx_1bo=out.saccades(k).ini_all>out.states(k).start_1bo;
                    idx_2bo=out.saccades(k).ini_all>out.states(k).start_2bo;
                                        
                    
                    if c<=numel(condition_tags)/2 
                        % latencies in stimulation conditions
%                         current_ini_1ao=[out.saccades(k).ini_1ao]-out.task(k).ini_mis;
%                         current_ini_obs=[out.saccades(k).ini_obs]-out.task(k).ini_mis;
%                         current_ini_1bo=[out.saccades(k).ini_1bo]-out.task(k).ini_mis;
%                         current_ini_2bo=[out.saccades(k).ini_2bo]-out.task(k).ini_mis;
                        microstim_shift=out.task(k).ini_mis;
                    elseif c>numel(condition_tags)/2 
                        % baseline shifting dependend on relevant_start, stimstate, stimstart_abs
                        if strcmp(relevant_start,'stim_start')
                            if out.states(k).state_obs==stimstate;
                                microstim_shift=out.states(k).start_obs;
                            elseif out.states(k).state_1bo==stimstate;
                                microstim_shift=out.states(k).start_1bo;
                            elseif out.states(k).state_2bo==stimstate;
                                microstim_shift=out.states(k).start_2bo;
                            else % run baseline !
                                microstim_shift=out.states(k).start_obs;
                            end
                            microstim_shift=microstim_shift+stimstart_abs;
                        elseif strcmp(relevant_start,'stim_to_state_end')
                            if out.states(k).state_obs==stimstate;
                                microstim_shift=out.states(k).start_1ao;
                            elseif out.states(k).state_1bo==stimstate;
                                microstim_shift=out.states(k).start_obs;
                            elseif out.states(k).state_2bo==stimstate;
                                microstim_shift=out.states(k).start_1bo;
                            else % run baseline !
                                microstim_shift=out.states(k).start_obs;
                            end
                            microstim_shift=microstim_shift-stimstart_abs;
                        end
                        % latencies in baseline condition
                    end
                               
                        current_ini_1ao=out.saccades(k).ini_all(idx_1ao) - microstim_shift;
                        current_ini_obs=out.saccades(k).ini_all(idx_obs) - microstim_shift;
                        current_ini_1bo=out.saccades(k).ini_all(idx_1bo) - microstim_shift;
                        current_ini_2bo=out.saccades(k).ini_all(idx_2bo) - microstim_shift;         
                        all_shifted_ini=[out.saccades(k).ini_all] - microstim_shift;      
                        all_directions=[out.saccades(k).endpos_all-out.saccades(k).startpos_all];
%                     idx_ini_1ao=current_ini_1ao>=0 & current_ini_1ao<=max(Bins);
%                     idx_ini_obs=current_ini_obs>=0 & current_ini_obs<=max(Bins);
%                     idx_ini_1bo=current_ini_1bo>=0 & current_ini_1bo<=max(Bins);
%                     idx_ini_2bo=current_ini_2bo>=0 & current_ini_2bo<=max(Bins);
%                     
                    timing_index    = out.saccades(k).ini_all - microstim_shift>= 0 & out.saccades(k).ini_all - microstim_shift <= max(Bins);
                    idx_ini_1ao     = idx_1ao & timing_index;
                    idx_ini_obs     = idx_obs & timing_index;
                    idx_ini_1bo     = idx_1bo & timing_index;
                    idx_ini_2bo     = idx_2bo & timing_index;
                    
%                     nanfill=NaN(1,40-sum(idx_ini_1ao)-sum(idx_ini_obs)-sum(idx_ini_1bo)-sum(idx_ini_2bo));
%                     lat_obs(k,:)= [current_ini_1ao(idx_ini_1ao) current_ini_obs(idx_ini_obs) current_ini_1bo(idx_ini_1bo) current_ini_2bo(idx_ini_2bo) nanfill];
%                     directions_obs(k,:)     =   [out.saccades(k).endpos_1ao(idx_ini_1ao)-out.saccades(k).startpos_1ao(idx_ini_1ao)...
%                         out.saccades(k).endpos_obs(idx_ini_obs)-out.saccades(k).startpos_obs(idx_ini_obs)...
%                         out.saccades(k).endpos_1bo(idx_ini_1bo)-out.saccades(k).startpos_1bo(idx_ini_1bo)...
%                         out.saccades(k).endpos_2bo(idx_ini_2bo)-out.saccades(k).startpos_2bo(idx_ini_2bo) nanfill];
%                     amplitudes_obs=[out.saccades(k).amplitudes_1ao(idx_ini_1ao) out.saccades(k).amplitudes_obs(idx_ini_obs)...
%                         out.saccades(k).amplitudes_1bo(idx_ini_1bo) out.saccades(k).amplitudes_2bo(idx_ini_2bo) nanfill];
%                     
                    nanfill                 = NaN(1,40-sum(idx_ini_1ao)-sum(idx_ini_obs)-sum(idx_ini_1bo)-sum(idx_ini_2bo));
                    lat_obs(k,:)            = [all_shifted_ini(idx_ini_1ao) all_shifted_ini(idx_ini_obs) all_shifted_ini(idx_ini_1bo) all_shifted_ini(idx_ini_2bo) nanfill];
                    directions_obs(k,:)     = [all_directions(idx_ini_1ao)  all_directions(idx_ini_obs)  all_directions(idx_ini_1bo)  all_directions(idx_ini_2bo)  nanfill];
                     amplitudes_obs = abs(directions_obs(k,:));
                    
                    
                    
                end
                
                matching_amplitudes_R   = amplitudes(m)<= amplitudes_obs  & amplitudes_obs <=amplitudes(m+1) & real(directions_obs(k,:)) > 0;
                matching_amplitudes_L   = amplitudes(m)<= amplitudes_obs  & amplitudes_obs <=amplitudes(m+1) & real(directions_obs(k,:)) < 0;
                matching_amplitudes_L_return   = (amplitudes(m)<= amplitudes_obs  & amplitudes_obs <=amplitudes(m+1) & real(directions_obs(k,:)) < 0 & lat_obs(k,:) > 0.0 & lat_obs(k,:) < 0.2 & any(real(directions_obs(k,:)) > 0 & lat_obs(k,:) > 0.2 & lat_obs(k,:) < 0.4 )) |...
                                                 (amplitudes(m)<= amplitudes_obs  & amplitudes_obs <=amplitudes(m+1) & real(directions_obs(k,:)) > 0 & lat_obs(k,:) > 0.2 & lat_obs(k,:) < 0.4 & any(real(directions_obs(k,:)) < 0 & lat_obs(k,:) > 0.0 & lat_obs(k,:) < 0.2 ));
                matching_amplitudes_R_return   = (amplitudes(m)<= amplitudes_obs  & amplitudes_obs <=amplitudes(m+1) & real(directions_obs(k,:)) > 0 & lat_obs(k,:) > 0.0 & lat_obs(k,:) < 0.2 & any(real(directions_obs(k,:)) < 0 & lat_obs(k,:) > 0.2 & lat_obs(k,:) < 0.4 )) |...
                                                 (amplitudes(m)<= amplitudes_obs  & amplitudes_obs <=amplitudes(m+1) & real(directions_obs(k,:)) < 0 & lat_obs(k,:) > 0.2 & lat_obs(k,:) < 0.4 & any(real(directions_obs(k,:)) > 0 & lat_obs(k,:) > 0.0 & lat_obs(k,:) < 0.2 ));
                matching_amplitude_all  = matching_amplitudes_R | matching_amplitudes_L;                
                RTs_per_amp_R.(condition_tags{c}){m}      =[RTs_per_amp_R.(condition_tags{c}){m} lat_obs(k,matching_amplitudes_R)];
                RTs_per_amp_L.(condition_tags{c}){m}      =[RTs_per_amp_L.(condition_tags{c}){m} lat_obs(k,matching_amplitudes_L)];
                RTs_per_amp_L_with_return.(condition_tags{c}){m}      =[RTs_per_amp_L_with_return.(condition_tags{c}){m} lat_obs(k,matching_amplitudes_L_return)];
                amplitudes_L_with_return.(condition_tags{c}){m}       =[amplitudes_L_with_return.(condition_tags{c}){m} amplitudes_obs(matching_amplitudes_L_return)];
                amplitudes_L.(condition_tags{c}){m}             =[amplitudes_L.(condition_tags{c}){m} amplitudes_obs(matching_amplitudes_L)];
                amplitudes_R.(condition_tags{c}){m}             =[amplitudes_R.(condition_tags{c}){m} amplitudes_obs(matching_amplitudes_R)];
                Dir_all.(condition_tags{c})        =[Dir_all.(condition_tags{c}) angle(directions_obs(k,matching_amplitude_all))];
                Amp_all.(condition_tags{c})        =[Amp_all.(condition_tags{c}) abs(amplitudes_obs(matching_amplitude_all))];
                RTs_all.(condition_tags{c})        =[RTs_all.(condition_tags{c}) lat_obs(k,matching_amplitude_all)];
                Indexes_all.L_with_ret.(condition_tags{c})        =[Indexes_all.L_with_ret.(condition_tags{c}) matching_amplitudes_L_return(matching_amplitude_all)];
                Indexes_all.R_with_ret.(condition_tags{c})        =[Indexes_all.R_with_ret.(condition_tags{c}) matching_amplitudes_R_return(matching_amplitude_all)];
                
%                 if strcmp(condition_tags{c},'stim') && any(lat_obs(k,real(directions_obs(k,:)) < -0.5)<0.2) && any(lat_obs(k,real(directions_obs(k,:)) > 0.5)>0.2 & lat_obs(k,real(directions_obs(k,:)) > 0.5)<0.4)
%                    n_evoked_and_return = n_evoked_and_return+1;
%                 end
%                 if strcmp(condition_tags{c},'stim') && any(lat_obs(k,real(directions_obs(k,:)) > 0.5)>0.2 & lat_obs(k,real(directions_obs(k,:)) > 0.5)<0.4)
%                    n_return = n_return+1;
%                 end
%                 if strcmp(condition_tags{c},'stim') && any(lat_obs(k,real(directions_obs(k,:)) < -0.5)<0.2)
%                    n_evoked = n_evoked+1;
%                 end
%                 if strcmp(condition_tags{c},'base') && any(lat_obs(k,real(directions_obs(k,:)) > 0.5)>0.2 & lat_obs(k,real(directions_obs(k,:)) > 0.5)<0.4)
%                    
%                 n_base_return= n_base_return+1;
%                 end
            end
        end
        

%         idx_i=any(lat_obs>0.2&lat_obs<0.4,2);
%         ev_and_ret=sum(idx_e(1:size(lat_obs,1))' & any(lat_obs>0.2 & lat_obs<0.4 & real(directions_obs) > 0,2));
%         ev        =sum(idx_e);
%         ret       =sum(any(lat_obs(idx.stim,:)>0.2 & lat_obs(idx.stim,:)<0.4 & real(directions_obs(idx.stim,:)) > 0,2));
%         base_ret  =sum(any(lat_obs(idx.base(1:size(lat_obs,1)),:)>0.2 & lat_obs(idx.base(1:size(lat_obs,1)),:)<0.4 & real(directions_obs(idx.base(1:size(lat_obs,1)),:)) > 0.5,2));
        % sum(idx.stim(1:size(lat_obs,1))' & any(lat_obs>=0 & lat_obs<0.2 & real(directions_obs) <= -0.3,2))
        if n_trials.(condition_tags{c})==0
            current_N=1;
        else
            current_N=n_trials.(condition_tags{c});
        end
        
        % Histograms of Reaction times
        RTs_histograms_R.(condition_tags{c})(m,:)                   =hist(RTs_per_amp_R.(condition_tags{c}){m},Bins)./current_N;
        RTs_histograms_L.(condition_tags{c})(m,:)                   =hist(RTs_per_amp_L.(condition_tags{c}){m},Bins)./current_N;
        RTs_histograms_all.(condition_tags{c})(m,:)                 =RTs_histograms_R.(condition_tags{c})(m,:)+RTs_histograms_L.(condition_tags{c})(m,:);        
    end    
    
            % Calculate evoked probability for current condition
        total       =sum(idx.(condition_tags{c}));
        amp_log_idx_evo=abs(abs(directions_obs(1:end,1:end)))>=  Evo_amp_TH(1) & abs(abs(directions_obs(1:end,1:end)))<=  Evo_amp_TH(2);
        ev_l        =sum(idx.(condition_tags{c})' & any(amp_log_idx_evo & lat_obs>=0   & lat_obs<=0.2 & real(directions_obs) < 0 ,2));
        ev_r        =sum(idx.(condition_tags{c})' & any(amp_log_idx_evo & lat_obs>=0   & lat_obs<=0.2 & real(directions_obs) > 0 ,2));
        ret_l       =sum(idx.(condition_tags{c})' & any(amp_log_idx_evo & lat_obs>=0.2 & lat_obs<=0.4 & real(directions_obs) < 0 ,2));
        ret_r       =sum(idx.(condition_tags{c})' & any(amp_log_idx_evo & lat_obs>=0.2 & lat_obs<=0.4 & real(directions_obs) > 0 ,2));
        ev_and_ret=sum(idx.(condition_tags{c})'   & any(amp_log_idx_evo & lat_obs>=0   & lat_obs<=0.2 & real(directions_obs) < 0 ,2) &...
                                                    any(amp_log_idx_evo & lat_obs>=0.2 & lat_obs<=0.4 & real(directions_obs) > 0 ,2));
%         ev_l_2      =sum(idx.(condition_tags{c})' & sum(amp_log_idx & lat_obs>=0   & lat_obs<=0.2 & real(directions_obs) < 0 ,2)>=2); 
%         ev_r_2      =sum(idx.(condition_tags{c})' & sum(amp_log_idx & lat_obs>=0   & lat_obs<=0.2 & real(directions_obs) > 0 ,2)>=2); 
%         
        amp_log_idx_del=abs(abs(directions_obs(1:end,1:end)))>=  Del_amp_TH(1) & abs(abs(directions_obs(1:end,1:end)))<=  Del_amp_TH(2);
        del        =sum(idx.(condition_tags{c})' & any(amp_log_idx_del & lat_obs>=0.2,2));
        %de_r        =sum(idx.(condition_tags{c})' & any(amp_log_idx_del & lat_obs>=0.2  & real(directions_obs) > 0 ,2));
        %idx_e=idx.(condition_tags{c}) & idx_evoked;
        %idx_ne=idx.(condition_tags{c}) & ~idx_evoked;
        %idx_r
        %idx_e_and_r
        
        if total ~= 0
            P.evoked_l.(condition_tags{c})            = ev_l/total;
            P.evoked_r.(condition_tags{c})            = ev_r/total;
            P.return_l.(condition_tags{c})            = ret_l/total;
            P.return_r.(condition_tags{c})            = ret_r/total;
            P.delayed.(condition_tags{c})             = del/total;
            P.evoked_and_return.(condition_tags{c})   = ev_and_ret/total;
        else
            P.evoked_l.(condition_tags{c})            = 0;
            P.evoked_r.(condition_tags{c})            = 0;
            P.return_l.(condition_tags{c})            = 0;
            P.return_r.(condition_tags{c})            = 0;
            P.delayed.(condition_tags{c})             = 0;
            P.evoked_and_return.(condition_tags{c})   = 0;
        end
        
end
end

% Implemented  tools

function [filelist_formatted,task_type,stimstart,stimstate]=get_filelist_from_matfile(working_data_folder,mastertable)
column_index_session=DAG_find_column_index(mastertable,'Session');
column_index_run=DAG_find_column_index(mastertable,'Run');
column_index_Task_type=DAG_find_column_index(mastertable,'task_type');
column_index_stimstart=DAG_find_column_index(mastertable,'microstim_start');
column_index_stimstate=DAG_find_column_index(mastertable,'microstim_state');

base_path_cell=repmat({[working_data_folder, filesep]},size(mastertable,1)-1,1);
session_cell=cellstr(num2str([mastertable{2:end,column_index_session}]'));
filelist_formatted=[strcat(base_path_cell, session_cell),mastertable(2:end,column_index_run)];

stimstart=[mastertable(2:end,column_index_stimstart)];
stimstate=[mastertable(2:end,column_index_stimstate)];
task_type=[mastertable(2:end,column_index_Task_type)];
end

function column_index=DAG_find_column_index(inputcell,title)
for m=1:size(inputcell,2)
    if strcmp(inputcell{1,m},title)
        column_index=m;
    end
end
end

function STATE=load_state


STATE.INI_TRI           = 1;  % initialize trial
STATE.FIX_ACQ           = 2;  % fixation acquisition
STATE.FIX_HOL           = 3;  % fixation hold
STATE.TAR_ACQ           = 4;  % target acquisition
STATE.TAR_HOL           = 5;  % target hold
STATE.CUE_ON            = 6;  % cue on
STATE.MEM_PER           = 7;  % memory period
STATE.DEL_PER           = 8;  % delay period
STATE.TAR_ACQ_INV       = 9;  % target acquisition invisible
STATE.TAR_HOL_INV       = 10; % target hold invisible
STATE.ABORT             = 19;
STATE.SUCCESS           = 20;
STATE.REWARD            = 21;
STATE.ITI               = 50;
STATE.CLOSE             = 99;
STATE.SUCCESS_ABORT     = -1;

STATE.ALL               =[STATE.INI_TRI STATE.FIX_ACQ  STATE.FIX_HOL  STATE.CUE_ON  STATE.MEM_PER  STATE.DEL_PER STATE.TAR_ACQ_INV  STATE.TAR_HOL_INV...
    STATE.TAR_ACQ  STATE.TAR_HOL STATE.SUCCESS_ABORT  STATE.REWARD  STATE.ITI];
STATE.ALL_NAMES         ={'Trial Initiation', 'Fixation Acquisition', 'Fixation Hold', 'Cue', 'Memory', 'Delay', 'Target Acquisition (inv)',...
    'Target Hold (inv)','Target Acquisition', 'Target Hold', 'Success', 'Reward', 'ITI'};
STATE.ALL_NAMES_Short         ={'INI', 'FIX_ACQ', 'FIX_HOL', 'CUE', 'MEM', 'DEL', 'TAR_ACQ_INV',...
    'TAR_HOL_INV','TAR_ACQ', 'TAR_HOLd', 'SUC', 'REW', 'ITI'};

STATE.aquisition=[STATE.FIX_ACQ STATE.TAR_ACQ STATE.TAR_ACQ_INV];
end


