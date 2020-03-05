function get_evoked_mastertable
run set_settings
global analysis_parameters
excentricity_calculated     = analysis_parameters.batches.excentricity;
monkey                      = analysis_parameters.monkey;
current_date                = analysis_parameters.current_date;
batch_processing            = analysis_parameters.batch_processing;
evoked_definition           = analysis_parameters.evoked_definition;

if batch_processing
    load(DAG_most_recent_version(pwd,'histograms_batched'));
else
    load(DAG_most_recent_version(pwd,'saccade_histograms'));
end



valid_batches=saccade_histograms(1).valid_batches;
load(DAG_most_recent_version(pwd,'trialinfo_mastertable'));
evoked_summary_info=evoked_probability_per_run(saccade_histograms,excentricity_calculated,evoked_definition);
save([monkey, '_evoked_summary_info_', current_date],'evoked_summary_info');
mastertable_extension  = write_summary_to_masterfile_simple(evoked_summary_info);
mastertable = DAG_update_mastertable_cell(mastertable,mastertable_extension,valid_batches+1);
save([monkey, '_evoked_mastertable_', current_date],'mastertable');

end

function evoked_summary_info = evoked_probability_per_run(saccade_histograms,excentricity_calculated,evoked_definition)
evoked_summary_info=[];
Trialtypes=fieldnames(saccade_histograms);

for p=1:size(saccade_histograms,2)
    for t=1:numel(Trialtypes)
        % breaking conditions
        if ~isstruct(saccade_histograms(p).(Trialtypes{t}))
           continue 
        end        
        current_type = saccade_histograms(p).(Trialtypes{t});
        excentricity=current_type.excentricity;
        if isempty(strfind(excentricity,excentricity_calculated)) && isempty(strfind(Trialtypes{t},'fix'))
           continue 
        end
        
        stim_windows        =current_type.stim_windows;        
        Evo_RT_TH           =current_type.Evo_RT_TH;
        Evo_amp_TH          =current_type.Evo_amp_TH;        
        Del_RT_TH           =current_type.Del_RT_TH;
        Del_amp_TH          =current_type.Del_amp_TH;        
        
        logidx_delayed_RTs=[current_type.RT_bins]>Del_RT_TH(1) & [current_type.RT_bins]<=Del_RT_TH(2);
        logidx_delayed_amps=[current_type.amplitudes(2:end)]>=Del_amp_TH(1) & [current_type.amplitudes(1:end-1)]<Del_amp_TH(2);   
        n_trial_fieldnames=fieldnames([current_type.summed_n.trials]);
        for k=1:numel(n_trial_fieldnames)
            
            Evoked_P_per_condition_window_L.(n_trial_fieldnames{k}) = NaN(size(stim_windows,1),1);
            %Evoked_amplitudes_per_window_R.(n_trial_fieldnames{k})  = NaN(size(stim_windows,1),1);
            Evoked_P_per_condition_window_R.(n_trial_fieldnames{k}) = NaN(size(stim_windows,1),1);
            Delayed_P_per_condition_window.(n_trial_fieldnames{k}) = NaN(size(stim_windows,1),1);
            All_evoked_RTs.(n_trial_fieldnames{k})          =[];
            All_evoked_amplitudes.(n_trial_fieldnames{k})   =[];
            
            N_evoked_early.(n_trial_fieldnames{k})      =     0;            
            N_delayed_early.(n_trial_fieldnames{k})     =     0;
            N_early.(n_trial_fieldnames{k})             =     0;            
            N_evoked_late.(n_trial_fieldnames{k})       =     0;            
            N_delayed_late.(n_trial_fieldnames{k})      =     0;
            N_late.(n_trial_fieldnames{k})              =     0;
            
            Evoked_RTs_early.(n_trial_fieldnames{k})    =     [];
            Evoked_amps_early.(n_trial_fieldnames{k})   =     [];
            Evoked_RTs_late.(n_trial_fieldnames{k})     =     [];
            Evoked_amps_late.(n_trial_fieldnames{k})    =     [];
            
            for stim_start=1:size(stim_windows,1)
                current_stim_onset=[current_type.stimulation_onset];
                stim_idx= [current_type.stimulation_state]==stim_windows(stim_start,1) & current_stim_onset==stim_windows(stim_start,2);
                if any(stim_idx)
                    current_L_RTs                    =[current_type.RTs_per_amp_L(stim_idx).(n_trial_fieldnames{k}){:}];      
                    current_L_RTs_with_return        =[current_type.RTs_per_amp_L_with_return(stim_idx).(n_trial_fieldnames{k}){:}];  
                    current_L_amplitudes             =[current_type.amplitudes_L(stim_idx).(n_trial_fieldnames{k}){:}];  
                    current_L_amplitudes_with_return  =[current_type.amplitudes_L_with_return(stim_idx).(n_trial_fieldnames{k}){:}];                  
                    if strcmp(evoked_definition,'evoked_l')
                        idx_evoked     = current_L_RTs>Evo_RT_TH(1) & current_L_RTs<=Evo_RT_TH(2) & current_L_amplitudes > Evo_amp_TH(1) &  current_L_amplitudes < Evo_amp_TH(2);
                        current_L_RTs                    =[current_type.RTs_per_amp_L(stim_idx).(n_trial_fieldnames{k}){:}];
                        current_L_amplitudes             =[current_type.amplitudes_L(stim_idx).(n_trial_fieldnames{k}){:}];
                    elseif strcmp(evoked_definition,'evoked_and_return')
                        current_L_RTs        =[current_type.RTs_per_amp_L_with_return(stim_idx).(n_trial_fieldnames{k}){:}];
                        current_L_amplitudes  =[current_type.amplitudes_L_with_return(stim_idx).(n_trial_fieldnames{k}){:}];
                        idx_evoked     = current_L_RTs_with_return>Evo_RT_TH(1) & current_L_RTs_with_return<=Evo_RT_TH(2) & current_L_amplitudes_with_return > Evo_amp_TH(1) &  current_L_amplitudes_with_return < Evo_amp_TH(2);
                    else
                        current_L_RTs                    =[current_type.RTs_per_amp_L(stim_idx).(n_trial_fieldnames{k}){:}];
                        current_L_amplitudes             =[current_type.amplitudes_L(stim_idx).(n_trial_fieldnames{k}){:}];
                        disp([evoked_definition ' is not a known defintion for ecoked saccades'])
                        idx_evoked     = current_L_RTs>Evo_RT_TH(1) & current_L_RTs<=Evo_RT_TH(2) & current_L_amplitudes > Evo_amp_TH(1) &  current_L_amplitudes < Evo_amp_TH(2);
                    end
                    current_evoked_RT_histograms_L  =[current_type.RTs_histograms_L(stim_idx).(n_trial_fieldnames{k})];
                    current_evoked_RT_histograms_R   =[current_type.RTs_histograms_R(stim_idx).(n_trial_fieldnames{k})];current_all_amplitudes           =[current_type.amplitudes_L(stim_idx).(n_trial_fieldnames{k}){:} current_type.amplitudes_R(stim_idx).(n_trial_fieldnames{k}){:}]; 
                    current_all_RTs                  =[current_type.RTs_per_amp_L(stim_idx).(n_trial_fieldnames{k}){:} current_type.RTs_per_amp_R(stim_idx).(n_trial_fieldnames{k}){:}];
                   
                   
                        
                    idx_delayed    = current_all_RTs>Del_RT_TH(1) & current_all_RTs<=Del_RT_TH(2) & current_all_amplitudes > Del_amp_TH(1) &  current_all_amplitudes < Del_amp_TH(2);
                    All_evoked_RTs.(n_trial_fieldnames{k})           =[All_evoked_RTs.(n_trial_fieldnames{k}) current_L_RTs(idx_evoked)];
                    All_evoked_amplitudes.(n_trial_fieldnames{k})    =[All_evoked_amplitudes.(n_trial_fieldnames{k}) current_L_amplitudes(idx_evoked)];
                    
                    Evoked_mean_RTs_per_condition_window_L.(n_trial_fieldnames{k})(stim_start)   =  nanmean(current_L_RTs(idx_evoked));
                    Evoked_std_RTs_per_condition_window_L.(n_trial_fieldnames{k})(stim_start)    =  nanstd(current_L_RTs(idx_evoked));
                    Evoked_mean_amplitudes_per_window_L.(n_trial_fieldnames{k})(stim_start)      =  nanmean(current_L_amplitudes(idx_evoked));
                    Evoked_std_amplitudes_per_window_L.(n_trial_fieldnames{k})(stim_start)       =  nanstd(current_L_amplitudes(idx_evoked));
%                     Evoked_P_per_condition_window_L.(n_trial_fieldnames{k})(stim_start)          =  sum(sum(current_evoked_RT_histograms_L(logidx_evoked_amps,logidx_evoked_RTs)));
%                     Evoked_P_per_condition_window_R.(n_trial_fieldnames{k})(stim_start)          =  sum(sum(current_evoked_RT_histograms_R(logidx_evoked_amps,logidx_evoked_RTs)));
                    %Delayed_P_per_condition_window.(n_trial_fieldnames{k})(stim_start)           =  sum(sum(current_evoked_RT_histograms_R(logidx_delayed_amps,logidx_delayed_RTs)))+sum(sum(current_evoked_RT_histograms_L(logidx_delayed_amps,logidx_delayed_RTs)));
                    %Delayed_P_per_condition_window.(n_trial_fieldnames{k})(stim_start)           =  sum(current_all_RTs(idx_delayed)>Del_RT_TH(1) & current_all_RTs(idx_delayed)<=Del_RT_TH(2));
                    
                    
                    %Delayed_P_per_condition_window.(n_trial_fieldnames{k})(stim_start)           =  sum(current_all_RTs(idx_delayed)>Del_RT_TH(1) & current_all_RTs(idx_delayed)<=Del_RT_TH(2));
                    %Delayed_P_per_condition_window.(n_trial_fieldnames{k})(stim_start)           =(1-sum(sum(current_evoked_RT_histograms_L(logidx_delayed_amps,logidx_delayed_RTs))) - sum(sum(current_evoked_RT_histograms_R(logidx_delayed_amps,logidx_delayed_RTs)))).*current_type.n_trials(stim_idx).(n_trial_fieldnames{k});
                    % dividing into early and late windows
                    if current_stim_onset(stim_idx) < 0
                        N_evoked_early.(n_trial_fieldnames{k})      =     N_evoked_early.(n_trial_fieldnames{k}) + current_type.n.(evoked_definition)(stim_idx).(n_trial_fieldnames{k});
                        N_early.(n_trial_fieldnames{k})             =     N_early.(n_trial_fieldnames{k})        + current_type.n_trials(stim_idx).(n_trial_fieldnames{k});
                        %N_delayed_early.(n_trial_fieldnames{k})     =     N_delayed_early.(n_trial_fieldnames{k}) + Delayed_P_per_condition_window.(n_trial_fieldnames{k})(stim_start)*current_type.n_trials(stim_idx).(n_trial_fieldnames{k});
                        %N_delayed_early.(n_trial_fieldnames{k})     =     N_delayed_early.(n_trial_fieldnames{k}) + Delayed_P_per_condition_window.(n_trial_fieldnames{k})(stim_start);
                        N_delayed_early.(n_trial_fieldnames{k})     =     N_delayed_early.(n_trial_fieldnames{k}) + current_type.n.delayed(stim_idx).(n_trial_fieldnames{k});
                        
                        Evoked_RTs_early.(n_trial_fieldnames{k})    =     [Evoked_RTs_early.(n_trial_fieldnames{k}) current_L_RTs(idx_evoked)];
                        Evoked_amps_early.(n_trial_fieldnames{k})   =     [Evoked_amps_early.(n_trial_fieldnames{k}) current_L_amplitudes(idx_evoked)];
                    elseif current_stim_onset(stim_idx) >= 0
                        N_evoked_late.(n_trial_fieldnames{k})       =     N_evoked_late.(n_trial_fieldnames{k})  + current_type.n.(evoked_definition)(stim_idx).(n_trial_fieldnames{k});
                        N_late.(n_trial_fieldnames{k})              =     N_late.(n_trial_fieldnames{k})         + current_type.n_trials(stim_idx).(n_trial_fieldnames{k});
                        %N_delayed_late.(n_trial_fieldnames{k})      =     N_delayed_late.(n_trial_fieldnames{k}) + Delayed_P_per_condition_window.(n_trial_fieldnames{k})(stim_start)*current_type.n_trials(stim_idx).(n_trial_fieldnames{k});
                        %N_delayed_late.(n_trial_fieldnames{k})      =     N_delayed_late.(n_trial_fieldnames{k}) + Delayed_P_per_condition_window.(n_trial_fieldnames{k})(stim_start);
                        N_delayed_late.(n_trial_fieldnames{k})      =     N_delayed_late.(n_trial_fieldnames{k}) + current_type.n.delayed(stim_idx).(n_trial_fieldnames{k});
                        
                        Evoked_RTs_late.(n_trial_fieldnames{k})     =     [Evoked_RTs_late.(n_trial_fieldnames{k}) current_L_RTs(idx_evoked)];
                        Evoked_amps_late.(n_trial_fieldnames{k})    =     [Evoked_amps_late.(n_trial_fieldnames{k}) current_L_amplitudes(idx_evoked)];                        
                    end 
                    
                else
                    Evoked_mean_RTs_per_condition_window_L.(n_trial_fieldnames{k})(stim_start)   =  NaN;
                    Evoked_std_RTs_per_condition_window_L.(n_trial_fieldnames{k})(stim_start)    =  NaN;
                    Evoked_mean_amplitudes_per_window_L.(n_trial_fieldnames{k})(stim_start)      =  NaN;
                    Evoked_std_amplitudes_per_window_L.(n_trial_fieldnames{k})(stim_start)       =  NaN;
%                     Evoked_P_per_condition_window_L.(n_trial_fieldnames{k})(stim_start)          =  NaN;
%                     Evoked_P_per_condition_window_R.(n_trial_fieldnames{k})(stim_start)          =  NaN;
                end
            end
        end
        
        P_evoked=current_type.P.(evoked_definition);
        clear tmp_struct
        if strfind(Trialtypes{t},'fix')
            tmp_struct.dP_evoked_fix        = round((nanmean([P_evoked.stim]-[P_evoked.base]))*100)./100;            
            tmp_struct.N_stim_fix           = sum([current_type.n_trials.stim]);
            tmp_struct.N_evoked_fix         = sum([current_type.n.(evoked_definition).stim]);   
            
            tmp_struct.Evoked_P_fix         = round((nanmean([P_evoked.stim]))*100)./100;
            tmp_struct.Evoked_mean_RT_fix   = max(nanmean([All_evoked_RTs.stim]));%,NaN);
            tmp_struct.Evoked_mean_amp_fix  = max(nanmean([All_evoked_amplitudes.stim]));%,NaN);
            
        elseif any(strfind(Trialtypes{t},'dir')) || any(strfind(Trialtypes{t},'mem'))
            tmp_struct.dP_evoked_R_CH           = round((nanmean([P_evoked.R_CH_S]-[P_evoked.R_CH_B]))*100)./100;
            tmp_struct.dP_evoked_R_IN           = round((nanmean([P_evoked.R_IN_S]-[P_evoked.R_IN_B]))*100)./100;            
            tmp_struct.dP_evoked_L_CH           = round((nanmean([P_evoked.L_CH_S]-[P_evoked.L_CH_B]))*100)./100;
            tmp_struct.dP_evoked_L_IN           = round((nanmean([P_evoked.L_IN_S]-[P_evoked.L_IN_B]))*100)./100;
            
            tmp_struct.bias_early     = N_early.L_CH_S / (N_early.L_CH_S + N_early.R_CH_S );
            tmp_struct.bias_late      = N_late.L_CH_S / (N_late.L_CH_S + N_late.R_CH_S );
            tmp_struct.bias_base      = current_type.n_trials(1).L_CH_B / (current_type.n_trials(1).L_CH_B  + current_type.n_trials(1).R_CH_B );
            
            
            tmp_struct.Evoked_P_early_L_IN      = N_evoked_early.L_IN_S / N_early.L_IN_S;
            tmp_struct.Evoked_P_early_R_IN      = N_evoked_early.R_IN_S / N_early.R_IN_S;
            tmp_struct.Evoked_P_early_L_CH      = N_evoked_early.L_CH_S / N_early.L_CH_S;
            tmp_struct.Evoked_P_early_R_CH      = N_evoked_early.R_CH_S / N_early.R_CH_S;
          
            tmp_struct.Evoked_P_late_L_IN       = N_evoked_late.L_IN_S / N_late.L_IN_S;
            tmp_struct.Evoked_P_late_R_IN       = N_evoked_late.R_IN_S / N_late.R_IN_S;
            tmp_struct.Evoked_P_late_L_CH       = N_evoked_late.L_CH_S / N_late.L_CH_S;
            tmp_struct.Evoked_P_late_R_CH       = N_evoked_late.R_CH_S / N_late.R_CH_S;
            
            tmp_struct.Delayed_P_early_L_IN     = N_delayed_early.L_IN_S / N_early.L_IN_S;
            tmp_struct.Delayed_P_early_R_IN     = N_delayed_early.R_IN_S / N_early.R_IN_S;
            tmp_struct.Delayed_P_early_L_CH     = N_delayed_early.L_CH_S / N_early.L_CH_S;
            tmp_struct.Delayed_P_early_R_CH     = N_delayed_early.R_CH_S / N_early.R_CH_S;
          
            tmp_struct.Delayed_P_late_L_IN      = N_delayed_late.L_IN_S / N_late.L_IN_S;
            tmp_struct.Delayed_P_late_R_IN      = N_delayed_late.R_IN_S / N_late.R_IN_S;
            tmp_struct.Delayed_P_late_L_CH      = N_delayed_late.L_CH_S / N_late.L_CH_S;
            tmp_struct.Delayed_P_late_R_CH      = N_delayed_late.R_CH_S / N_late.R_CH_S;
            
            tmp_struct.N_evoked_early_L_IN      = N_evoked_early.L_IN_S;
            tmp_struct.N_evoked_early_R_IN      = N_evoked_early.R_IN_S;
            tmp_struct.N_evoked_early_L_CH      = N_evoked_early.L_CH_S;
            tmp_struct.N_evoked_early_R_CH      = N_evoked_early.R_CH_S;
          
            tmp_struct.N_evoked_late_L_IN       = N_evoked_late.L_IN_S;
            tmp_struct.N_evoked_late_R_IN       = N_evoked_late.R_IN_S;
            tmp_struct.N_evoked_late_L_CH       = N_evoked_late.L_CH_S;
            tmp_struct.N_evoked_late_R_CH       = N_evoked_late.R_CH_S;
            
            tmp_struct.N_stim_early_L_IN         = N_early.L_IN_S;
            tmp_struct.N_stim_early_R_IN         = N_early.R_IN_S;
            tmp_struct.N_stim_early_L_CH         = N_early.L_CH_S;
            tmp_struct.N_stim_early_R_CH         = N_early.R_CH_S;
          
            tmp_struct.N_stim_late_L_IN          = N_late.L_IN_S;
            tmp_struct.N_stim_late_R_IN          = N_late.R_IN_S;
            tmp_struct.N_stim_late_L_CH          = N_late.L_CH_S;
            tmp_struct.N_stim_late_R_CH          = N_late.R_CH_S;
            
            tmp_struct.Evoked_mean_RT_early_L_IN               = nanmean([Evoked_RTs_early.L_IN_S]);
            tmp_struct.Evoked_mean_RT_early_L_CH               = nanmean([Evoked_RTs_early.L_CH_S]);
            tmp_struct.Evoked_mean_RT_early_R_IN               = nanmean([Evoked_RTs_early.R_IN_S]);
            tmp_struct.Evoked_mean_RT_early_R_CH               = nanmean([Evoked_RTs_early.R_CH_S]);            
            
            tmp_struct.Evoked_mean_amp_early_L_IN               = nanmean([Evoked_amps_early.L_IN_S]);
            tmp_struct.Evoked_mean_amp_early_L_CH               = nanmean([Evoked_amps_early.L_CH_S]);
            tmp_struct.Evoked_mean_amp_early_R_IN               = nanmean([Evoked_amps_early.R_IN_S]);
            tmp_struct.Evoked_mean_amp_early_R_CH               = nanmean([Evoked_amps_early.R_CH_S]);
                                               
            tmp_struct.Evoked_mean_RT_late_L_IN               = nanmean([Evoked_RTs_late.L_IN_S]);
            tmp_struct.Evoked_mean_RT_late_L_CH               = nanmean([Evoked_RTs_late.L_CH_S]);
            tmp_struct.Evoked_mean_RT_late_R_IN               = nanmean([Evoked_RTs_late.R_IN_S]);
            tmp_struct.Evoked_mean_RT_late_R_CH               = nanmean([Evoked_RTs_late.R_CH_S]);            
            
            tmp_struct.Evoked_mean_amp_late_L_IN               = nanmean([Evoked_amps_late.L_IN_S]);
            tmp_struct.Evoked_mean_amp_late_L_CH               = nanmean([Evoked_amps_late.L_CH_S]);
            tmp_struct.Evoked_mean_amp_late_R_IN               = nanmean([Evoked_amps_late.R_IN_S]);
            tmp_struct.Evoked_mean_amp_late_R_CH               = nanmean([Evoked_amps_late.R_CH_S]);
            
        end
        % Adding the extension for excentricities
        tmp_fieldnames=fieldnames(tmp_struct);
        for FN=1:numel(tmp_fieldnames)
            evoked_summary_info(p).(Trialtypes{t}).(tmp_fieldnames{FN})=tmp_struct.(tmp_fieldnames{FN});
        end
    end
end
end


function mastertable_extension = write_summary_to_masterfile_simple(info_structure)
% This function converts structure information to cell information in a very particular way: 
% Fieldnames will be stored in the first row, data in the following rows, one row for each element of the input structure info_structure. 
% Note that in this particular case, this is done for subfields of the input structure (Subfieldnames)

titles={};
Subfieldnames=fieldnames(info_structure);
for p=1:size(info_structure,2)
    for t=1:numel(Subfieldnames)
        A=info_structure(p).(Subfieldnames{t});
        if isstruct(A)
            titles_old=titles;
            titles_new=sort(fieldnames(A));
            titles=[titles_old; titles_new(~ismember(titles_new,titles_old))];
            for column=1:numel(titles)
                if isfield(A,titles{column})
                    mastertable_extension{1,column}=titles{column};
                    mastertable_extension{p+1,column}=A.(titles{column});
                end
            end
        end
    end
end
end
