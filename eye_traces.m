function [eye_x eye_y time] = eye_traces(out_comp,target_positions,batch_title,errors)
global GLO

STATE.INI_TRI       = 1; % initialize trial
STATE.FIX_ACQ       = 2; % fixation acquisition
STATE.FIX_HOL       = 3; % fixation hold
STATE.TAR_ACQ       = 4; % target acquisition
STATE.TAR_HOL       = 5; % target hold
STATE.CUE_ON        = 6; % cue on
STATE.MEM_PER       = 7; % memory period
STATE.DEL_PER       = 8; % delay period
STATE.TAR_ACQ_INV   = 9; % target acquisition invisible
STATE.TAR_HOL_INV   = 10; % target hold invisible
STATE.MAT_ACQ       = 11; % target acquisition in sample to match
STATE.MAT_HOL       = 12; % target acquisition in sample to match
STATE.MAT_ACQ_MSK   = 13; % target acquisition in sample to match
STATE.MAT_HOL_MSK   = 14; % target acquisition in sample to match
STATE.SEN_RET       = 15; % return to sensors for poffenberger
STATE.ABORT         = 19;
STATE.SUCCESS       = 20;
STATE.REWARD        = 21;
STATE.ITI           = 50;
STATE.CLOSE         = 99;


GLO.line_stim                       = GLO.All_stim_windows([GLO.windows_to_use])*1000;
GLO.train_duration                  = 200;
% GLO.Labels([1 2 GLO.windows_to_use+2])
Side_chosen_BL_binary                       =   [];
Side_chosen_stim_binary                     =   [];
num_hits                                    =   [];
num_hits_temp                               =   [];

idx=0;
for k = 1:numel(out_comp)
    if out_comp{k}.emptyflag==0
        idx=idx+1;
        notemptybatches(idx)=k;
    end
end

switch GLO.effector_to_use
    case 0
        sac_or_rea='saccades';
    case 4
        sac_or_rea='reaches';
end

lat_80_to_cue=[];
for k = notemptybatches
    switch GLO.effector_to_use
        case 0
            horizontal_distance = real([out_comp{k,1}.saccades.tar_pos] - [out_comp{k,1}.saccades.fix_pos]);
            vertical_distance   = imag([out_comp{k,1}.saccades.tar_pos] - [out_comp{k,1}.saccades.fix_pos]);
        case 4
            horizontal_distance = real([out_comp{k,1}.reaches.tar_pos] - [out_comp{k,1}.reaches.fix_pos]);
            vertical_distance   = imag([out_comp{k,1}.reaches.tar_pos] - [out_comp{k,1}.reaches.fix_pos]);
    end
    switch GLO.excentricity
        case 'all'
            idx_L{k}                            =   horizontal_distance<0;
            idx_R{k}                            =   horizontal_distance>0;
            
        case 'farther'
            
            idx_L{k}                            =   horizontal_distance<-GLO.low_excentricity_threshold;
            idx_R{k}                            =   horizontal_distance>GLO.low_excentricity_threshold;
        case 'closer'
            
            idx_L{k}                            =   horizontal_distance<0 & horizontal_distance>=-15;
            idx_R{k}                            =   horizontal_distance>0 & horizontal_distance<=15;
    end
    
    idx_direct{k}                               =   [out_comp{k,1}.task.type]==2;
    idx_memory{k}                               =   [out_comp{k,1}.task.type]==3;
    
    
    idx_choice{k}                               =   [out_comp{k,1}.binary.choice]==1;
    idx_instructed{k}                           =   [out_comp{k,1}.binary.choice]==0;
    
    idx_success{k}                              =   [out_comp{k,1}.binary.success]==1;
    idx_error{k}                                =   [out_comp{k,1}.binary.success]==0 & ~ismember({out_comp{k,1}.task.abort_code},'ABORT_JAW') & [out_comp{k,1}.states.state_abo]>2; %!!
    
    idx_stim{k}                                 =   [out_comp{k,1}.binary.microstim]==1;
    idx_baseline{k}                             =   [out_comp{k,1}.binary.microstim]==0;
    idx_baseline_cho{k}                         =   [out_comp{k,1}.binary.microstim]==0 & idx_choice{k};
    idx_baseline_ins{k}                         =   [out_comp{k,1}.binary.microstim]==0 & idx_instructed{k};
    
    idx_baseline_succ{k}                        =   idx_baseline{k}  & idx_success{k};
    idx_baseline_cho_succ{k}                    =   [out_comp{k,1}.binary.microstim]==0 & idx_choice{k}  & idx_success{k};
    idx_baseline_ins_succ{k}                    =   [out_comp{k,1}.binary.microstim]==0 & idx_instructed{k}  & idx_success{k};
    
    idx_baseline_unsucc{k}                      =   idx_baseline{k}  & idx_error{k};
    idx_baseline_cho_unsucc{k}                  =   [out_comp{k,1}.binary.microstim]==0 & idx_choice{k}  & idx_error{k};
    idx_baseline_ins_unsucc{k}                  =   [out_comp{k,1}.binary.microstim]==0 & idx_instructed{k}  & idx_error{k};
    
    idx_stim_at_0{k}                            =   [out_comp{k,1}.task.stim_start]==0 & idx_stim{k};
    idx_stim_in_fix{k}                          =   [out_comp{k,1}.task.stim_state]==STATE.FIX_HOL & idx_stim{k};
    idx_stim_in_acq{k}                          =   [out_comp{k,1}.task.stim_state]==STATE.TAR_ACQ & idx_stim{k};
    idx_stim_in_cue{k}                          =   [out_comp{k,1}.task.stim_state]==STATE.CUE_ON & idx_stim{k};
    idx_stim_in_mem{k}                          =   [out_comp{k,1}.task.stim_state]==STATE.MEM_PER & idx_stim{k};
    idx_stim_in_acq_inv{k}                      =   [out_comp{k,1}.task.stim_state]==STATE.TAR_ACQ_INV & idx_stim{k};
    
%     %% NASTY part for getting latencies of saccades towards the cue!!
%     idx_error_80_to_cue=[out_comp{k,1}.task.stim_to_state_end]==0.08 & idx_stim_in_fix{k} & idx_error{k} & idx_instructed{k};
%     %all_saccade_lats_during_cue=vertcat(out_comp{k,1}.saccades(idx_error_80_to_cue).ini_2bo) - repmat(vertcat(out_comp{k,1}.states(idx_error_80_to_cue).start_2bo),1,5);
%     all_saccade_lats_during_cue=vertcat(out_comp{k,1}.saccades(idx_error_80_to_cue).ini_2bo) - repmat(vertcat(out_comp{k,1}.task(idx_error_80_to_cue).ini_mis),1,5);
%     all_saccade_endpos_during_cue=vertcat(out_comp{k,1}.saccades(idx_error_80_to_cue).endpos_2bo) - repmat(vertcat(out_comp{k,1}.saccades(idx_error_80_to_cue).fix_pos),1,5);
%     ischoice_80_to_cue=vertcat(out_comp{k,1}.binary(idx_error_80_to_cue).choice);
%     
%     
%     all_saccade_endpos_during_cue_first=all_saccade_endpos_during_cue(:,1);
%     all_saccade_lats_during_cue_first=all_saccade_lats_during_cue(:,1);
%     trial_ends=cellfun(@(x) x(:,end),{out_comp{k,1}.raw(idx_error_80_to_cue).time_axis},'UniformOutput',false);
%     aborted_after_stim=([out_comp{k,1}.task(idx_error_80_to_cue).stim_start] + [out_comp{k,1}.states(idx_error_80_to_cue).start_mis])< [trial_ends{:}];
%     abort_states_80_to_cue=vertcat(out_comp{k,1}.states(idx_error_80_to_cue).state_abo);   
%     fraction_aborted_during_cue=sum(ismember(abort_states_80_to_cue,[6]))/sum(ismember(abort_states_80_to_cue,[4,5,6,7,9,10]));    
%     lat_80_to_cue=[lat_80_to_cue; all_saccade_lats_during_cue_first(~isnan(all_saccade_lats_during_cue_first))];
%     
%     N_saccades_during_cue=sum(~isnan(all_saccade_lats_during_cue_first))    
%     N_aborts_during_cue=sum(ismember(abort_states_80_to_cue,[6]))    
%     N_aborts_80_to_cue_after_fix=sum(ismember(abort_states_80_to_cue,[4,5,6,7,9,10]))  
%     N_aborts_80_to_cue_after_fix_acq=sum(ismember(abort_states_80_to_cue,[3,4,5,6,7,9,10]) & aborted_after_stim')      
%     N_fix_breaks_right_targets=sum(~isnan(all_saccade_lats_during_cue_first)&real(all_saccade_endpos_during_cue_first>0))    
%     N_saccades_and_aborts_during_cue=sum(~isnan(all_saccade_lats_during_cue_first) & ismember(abort_states_80_to_cue,[6]))     
%     N_saccades_and_aborts_during_cue_choices=sum(~isnan(all_saccade_lats_during_cue_first) & ismember(abort_states_80_to_cue,[6] & ischoice_80_to_cue))
%     
%     mean_lat=nanmean(lat_80_to_cue)
%     std_lat=nanstd(lat_80_to_cue)
    
    
%     %% NASTY part for getting percentage of contraversice undershooting errors!!
%     idx_error_80_after_go=[out_comp{k,1}.task.stim_start]==0.08 & idx_stim_in_acq_inv{k} & idx_error{k} & idx_instructed{k};
%    trial_ends=cellfun(@(x) x(:,end),{out_comp{k,1}.raw(idx_error_80_after_go).time_axis},'UniformOutput',false);
%    aborted_after_stim=([out_comp{k,1}.task(idx_error_80_after_go).stim_start]+ [out_comp{k,1}.states(idx_error_80_after_go).start_mis])< [trial_ends{:}];
%     abort_states_80_after_go=vertcat(out_comp{k,1}.states(idx_error_80_after_go).state_abo); 
%     all_saccade_endpos_after_go=vertcat(out_comp{k,1}.saccades(idx_error_80_after_go).endpos_obs) - repmat(vertcat(out_comp{k,1}.saccades(idx_error_80_after_go).fix_pos),1,5);
%     target_positions_after_go=vertcat(out_comp{k,1}.saccades(idx_error_80_after_go).tar_pos) - repmat(vertcat(out_comp{k,1}.saccades(idx_error_80_after_go).fix_pos),1,1);    
%     all_saccade_endpos_after_go_first=all_saccade_endpos_after_go(:,1);
%     
%     N_aborts_80_after_go_after_fix=sum(ismember(abort_states_80_after_go,[4,5,6,7,9,10])& aborted_after_stim')  
%     N_aborts_80_after_go_after_fix_acq=sum(ismember(abort_states_80_after_go,[3,4,5,6,7,9,10]) & aborted_after_stim')   
%     N_aborts_after_go=sum(ismember(abort_states_80_after_go,[4,5,9,10]))       
%     N_undershooting_left_saccades=sum(real(all_saccade_endpos_after_go_first<0))   
%     N_undershooting_left_saccades_and_abort_after_go=sum(real(all_saccade_endpos_after_go_first<0) & ismember(abort_states_80_after_go,[4,5,9,10]))   
%     N_undershooting_right_saccades=sum(real(all_saccade_endpos_after_go_first>0))   
%     
%     %N_undershooting_left_saccades_and_targets=sum(real(all_saccade_endpos_after_go_first<0) & real(target_positions_after_go<0))   
    
    
    for m=1:numel(GLO.All_stim_windows)
        switch GLO.All_stim_states{m}
            case 'fix'
                current_stim_state_index=idx_stim_in_fix;
            case 'cue'
                current_stim_state_index=idx_stim_in_cue;
            case 'mem'
                current_stim_state_index=idx_stim_in_mem;
            case 'tar_acq_inv'
                current_stim_state_index=idx_stim_in_acq_inv;
            case 'tar_acq'
                current_stim_state_index=idx_stim_in_acq;
        end
        
        if GLO.All_stim_windows(m)<0
            idx_shifting_window_all{m}              =   [out_comp{k,1}.task.stim_to_state_end]==abs(GLO.All_stim_windows(m))    & current_stim_state_index{k};
            idx_shifting_window_cho{m}              =   [out_comp{k,1}.task.stim_to_state_end]==abs(GLO.All_stim_windows(m))    & current_stim_state_index{k} & idx_choice{k};
            idx_shifting_window_ins{m}              =   [out_comp{k,1}.task.stim_to_state_end]==abs(GLO.All_stim_windows(m))    & current_stim_state_index{k} & idx_instructed{k};
            idx_shifting_window{m}                  =   [out_comp{k,1}.task.stim_to_state_end]==abs(GLO.All_stim_windows(m))    & current_stim_state_index{k};
        else
            idx_shifting_window_all{m}              =   [out_comp{k,1}.task.stim_start]==abs(GLO.All_stim_windows(m))           & current_stim_state_index{k};
            idx_shifting_window_cho{m}              =   [out_comp{k,1}.task.stim_start]==abs(GLO.All_stim_windows(m))           & current_stim_state_index{k} & idx_choice{k};
            idx_shifting_window_ins{m}              =   [out_comp{k,1}.task.stim_start]==abs(GLO.All_stim_windows(m))           & current_stim_state_index{k} & idx_instructed{k};
            idx_shifting_window{m}                  =   [out_comp{k,1}.task.stim_start]==abs(GLO.All_stim_windows(m))           & current_stim_state_index{k};
        end
        

            idx_shifting_window_unsucc_all{m}            =  idx_shifting_window_all{m} & idx_error{k};
            idx_shifting_window_unsucc_cho{m}            =  idx_shifting_window_cho{m} & idx_error{k};
            idx_shifting_window_unsucc_ins{m}            =  idx_shifting_window_ins{m} & idx_error{k};

            idx_shifting_window_succ_all{m}              =  idx_shifting_window_all{m} & idx_success{k};
            idx_shifting_window_succ_cho{m}              =  idx_shifting_window_cho{m} & idx_success{k};
            idx_shifting_window_succ_ins{m}              =  idx_shifting_window_ins{m} & idx_success{k};

        
    end

    if errors
        idx_shifting_window_current.all=idx_shifting_window_unsucc_all(GLO.windows_to_use);
        idx_shifting_window_current.cho=idx_shifting_window_unsucc_cho(GLO.windows_to_use);
        idx_shifting_window_current.ins=idx_shifting_window_unsucc_ins(GLO.windows_to_use);
    else
        idx_shifting_window_current.all=idx_shifting_window_succ_all(GLO.windows_to_use);
        idx_shifting_window_current.cho=idx_shifting_window_succ_cho(GLO.windows_to_use);
        idx_shifting_window_current.ins=idx_shifting_window_succ_ins(GLO.windows_to_use);
    end
    
    if errors
        idx_base.all=idx_baseline_unsucc;
        idx_base.cho=idx_baseline_cho_unsucc;
        idx_base.ins=idx_baseline_ins_unsucc;
    else
        idx_base.all=idx_baseline_succ;
        idx_base.cho=idx_baseline_cho_succ;
        idx_base.ins=idx_baseline_ins_succ;
    end
    
    fn= fieldnames(idx_shifting_window_current);
    
    idx_shifting_window     =idx_shifting_window     (GLO.windows_to_use);
    
    number_of_shifting_windows=numel(idx_shifting_window_current.all);
    
    
    % Per target position
    
    switch GLO.effector_to_use
        case 0
            tar_positions=[out_comp{k,1}.saccades.tar_pos]-[out_comp{k,1}.saccades.fix_pos];
            nct_positions=[out_comp{k,1}.saccades.nct_pos]-[out_comp{k,1}.saccades.fix_pos];
        case 4
            tar_positions=[out_comp{k,1}.reaches.tar_pos]-[out_comp{k,1}.reaches.fix_pos];
            nct_positions=[out_comp{k,1}.reaches.nct_pos]-[out_comp{k,1}.reaches.fix_pos];
    end
    
    
    for pos=1:numel(target_positions)
        idx_tar_pos{pos}=abs(tar_positions-target_positions(pos))<=GLO.target_pos_precision;
        for q=1:numel(fn)
            idx_input_baseline=find(idx_base.(fn{q}){k} & idx_tar_pos{pos});
            for idx_trial=1: numel(idx_input_baseline)
                sac_pos{k,1}.(fn{q}).position(pos).window(1).trial{idx_trial}=out_comp{k,1}.saccades(idx_input_baseline(idx_trial)).endpos - real(out_comp{k,1}.saccades(idx_input_baseline(idx_trial)).fix_pos);
                idx_acq_onset=find((idx_direct{k}(idx_input_baseline(idx_trial)) & [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==4]) |(idx_memory{k}(idx_input_baseline(idx_trial)) &   [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==6]));
                idx_go=find((idx_direct{k}(idx_input_baseline(idx_trial)) & [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==4]) |(idx_memory{k}(idx_input_baseline(idx_trial)) &   [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==9]));
                idx_state=(idx_direct{k}(idx_trial) & ([out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==3] | [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==4]  | [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==5])) |...
                    (idx_memory{k}(idx_trial) & ( [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==4] | [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==5] | ...
                    [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==6] | [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==7] | [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==9] | [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==10]));
                eye_x{k,1}.(fn{q}).position(pos).window(1).trial{idx_trial}=out_comp{k,1}.raw(idx_input_baseline(idx_trial)).x_eye(idx_state) - real([out_comp{k,1}.saccades(idx_input_baseline(idx_trial)).fix_pos]);
                eye_y{k,1}.(fn{q}).position(pos).window(1).trial{idx_trial}=out_comp{k,1}.raw(idx_input_baseline(idx_trial)).y_eye(idx_state);
                if isempty(idx_acq_onset)
                    time{k,1}.(fn{q}).position(pos).window(1).trial{idx_trial}=[];
                else
                idx_acq_onset=idx_acq_onset(1);                
                time{k,1}.(fn{q}).position(pos).window(1).trial{idx_trial}=out_comp{k,1}.raw(idx_input_baseline(idx_trial)).time_axis(idx_state)-out_comp{k,1}.raw(idx_input_baseline(idx_trial)).time_axis(idx_acq_onset);
                end
                ini_ms{k,1}.(fn{q}).position(pos).window(1).trial{idx_trial}=out_comp{k,1}.task(idx_input_baseline(idx_trial)).ini_mis;
                if  isempty(idx_go)
                    RT{k,1}.(fn{q}).position(pos).window(1).trial{idx_trial}=NaN;
                else
                    RT{k,1}.(fn{q}).position(pos).window(1).trial{idx_trial}=out_comp{k,1}.saccades(idx_input_baseline(idx_trial)).lat - out_comp{k,1}.raw(idx_input_baseline(idx_trial)).time_axis(idx_acq_onset) + out_comp{k,1}.raw(idx_input_baseline(idx_trial)).time_axis(idx_go(1));
                end
           end
            if isempty(idx_input_baseline)
                eye_x{k,1}.(fn{q}).position(pos).window(1).trial{1}=NaN;
                eye_y{k,1}.(fn{q}).position(pos).window(1).trial{1}=NaN;
                time{k,1}.(fn{q}).position(pos).window(1).trial{1}=NaN;
                RT{k,1}.(fn{q}).position(pos).window(1).trial{1}=NaN;
                ini_ms{k,1}.(fn{q}).position(pos).window(1).trial{1}=NaN;
                sac_pos{k,1}.(fn{q}).position(pos).window(1).trial{1}=NaN;
            end
            for idx_shifting_windows = 1:number_of_shifting_windows
                idx_input=find(idx_shifting_window_current.(fn{q}){idx_shifting_windows} & idx_tar_pos{pos});
                for idx_trial=1: numel(idx_input)
                    sac_pos{k,1}.(fn{q}).position(pos).window(idx_shifting_windows+1).trial{idx_trial}=out_comp{k,1}.saccades(idx_input(idx_trial)).endpos - real(out_comp{k,1}.saccades(idx_input(idx_trial)).fix_pos);
                    idx_acq_onset=find((idx_direct{k}(idx_input(idx_trial)) & [out_comp{k,1}.raw(idx_input(idx_trial)).states==4]) |(idx_memory{k}(idx_input(idx_trial)) &   [out_comp{k,1}.raw(idx_input(idx_trial)).states==6]));
                    idx_go=find((idx_direct{k}(idx_input(idx_trial)) & [out_comp{k,1}.raw(idx_input(idx_trial)).states==4]) |(idx_memory{k}(idx_input(idx_trial)) &   [out_comp{k,1}.raw(idx_input(idx_trial)).states==9]));
                    
                    idx_state=(idx_direct{k}(idx_trial) & ([out_comp{k,1}.raw(idx_input(idx_trial)).states==3] | [out_comp{k,1}.raw(idx_input(idx_trial)).states==4] | [out_comp{k,1}.raw(idx_input(idx_trial)).states==5])) |...
                        (idx_memory{k}(idx_trial) & ([out_comp{k,1}.raw(idx_input(idx_trial)).states==4] | [out_comp{k,1}.raw(idx_input(idx_trial)).states==5] | ...
                        [out_comp{k,1}.raw(idx_input(idx_trial)).states==6] | [out_comp{k,1}.raw(idx_input(idx_trial)).states==7] | [out_comp{k,1}.raw(idx_input(idx_trial)).states==9] | [out_comp{k,1}.raw(idx_input(idx_trial)).states==10]));
                    eye_x{k,1}.(fn{q}).position(pos).window(idx_shifting_windows+1).trial{idx_trial}=out_comp{k,1}.raw(idx_input(idx_trial)).x_eye(idx_state) - real([out_comp{k,1}.saccades(idx_input(idx_trial)).fix_pos]);
                    eye_y{k,1}.(fn{q}).position(pos).window(idx_shifting_windows+1).trial{idx_trial}=out_comp{k,1}.raw(idx_input(idx_trial)).y_eye(idx_state);
                    if isempty(idx_acq_onset)
                        time{k,1}.(fn{q}).position(pos).window(idx_shifting_windows+1).trial{idx_trial}=[];
                    else
                        idx_acq_onset=idx_acq_onset(1);
                        time{k,1}.(fn{q}).position(pos).window(idx_shifting_windows+1).trial{idx_trial}=out_comp{k,1}.raw(idx_input(idx_trial)).time_axis(idx_state)-out_comp{k,1}.raw(idx_input(idx_trial)).time_axis(idx_acq_onset);
                    end
                    if  isempty(idx_go)
                        RT{k,1}.(fn{q}).position(pos).window(idx_shifting_windows+1).trial{idx_trial}=NaN;
                    else
                        RT{k,1}.(fn{q}).position(pos).window(idx_shifting_windows+1).trial{idx_trial}=out_comp{k,1}.saccades(idx_input(idx_trial)).lat - out_comp{k,1}.raw(idx_input(idx_trial)).time_axis(idx_acq_onset) + out_comp{k,1}.raw(idx_input(idx_trial)).time_axis(idx_go(1));
                    end
                    ini_ms{k,1}.(fn{q}).position(pos).window(idx_shifting_windows+1).trial{idx_trial}=out_comp{k,1}.task(idx_input(idx_trial)).ini_mis;
                end
                if isempty(idx_input)
                    eye_x{k,1}.(fn{q}).position(pos).window(idx_shifting_windows+1).trial{1}=NaN;
                    eye_y{k,1}.(fn{q}).position(pos).window(idx_shifting_windows+1).trial{1}=NaN;
                    time{k,1}.(fn{q}).position(pos).window(idx_shifting_windows+1).trial{1}=NaN;
                    RT{k,1}.(fn{q}).position(pos).window(idx_shifting_windows+1).trial{1}=NaN;
                    ini_ms{k,1}.(fn{q}).position(pos).window(idx_shifting_windows+1).trial{1}=NaN;
                    sac_pos{k,1}.(fn{q}).position(pos).window(idx_shifting_windows+1).trial{1}=NaN;
                end
            end
            
        end
    end
    
    GLO.fontsize=20;
    
    col                                                         = jet(numel(GLO.Labels));
    eye_y_offset                                                = 20;
    col(1,:)                                                    = [0.2 0.2 0.2];
    col(2,:)                                                    = [0.5 0.5 0.5];
    if GLO.type_to_use                 == 3
        % elegant memory saccade colors
        col                                                    = [0.2 0.2 0.2;0.5 0.5 0.5;0.859 0.275 0.6; 0.65 0.247 0.6; 0.35 0.3412 0.6471; 0.122 0.255 0.604];
    end
    
    %subplot_assignment=[4 3 5 6 2 1];
    subplot_assignment=[5 6 2 1 4 3];
    %%%                [3 5 1 6 2 4]
    
    
    for q=1:numel(fn)
        plot_1_title='Raw eye traces x versus time';
        summary_1                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_1_title);
        for pos=1:numel(target_positions)
            
            subplot(3,2,subplot_assignment(pos)) % not correctly assigned yet
            hold on
            for idx_shifting_windows = 1:number_of_shifting_windows+1
                for idx_trial=1: numel(eye_x{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial)
                    if idx_shifting_windows>1 && isnan(ini_ms{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}(end)*1000)
                         plot(time{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}*1000,eye_x{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial},'color','k');
                    else
                        plot(time{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}*1000,eye_x{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial},'color',col(idx_shifting_windows+1,:));
                    end
                end
                %             line([GLO.line_stim{idx_shifting_windows}, GLO.line_stim{idx_shifting_windows}+GLO.train_duration], [25-idx_shifting_windows 25-idx_shifting_windows ],'color',col(idx_shifting_windows,:),'linewidth',2)
            end
            if all([idx_memory{k}])
                set(gca,'xlim', [-200,2500],'ylim', [-32,32],'FontSize',GLO.fontsize-6);
            else
                set(gca,'xlim', [-200,500],'ylim', [-32,32],'FontSize',GLO.fontsize-6);
            end
            
            plot(nanmean([RT{k,1}.(fn{q}).position(pos).window(1).trial{:}].*1000),-25,'^','color',col(2,:));
            for idx_shifting_windows = 1:number_of_shifting_windows
                line([GLO.line_stim(idx_shifting_windows), GLO.line_stim(idx_shifting_windows)+GLO.train_duration], [-2-idx_shifting_windows -2-idx_shifting_windows ],'color',col(idx_shifting_windows+2,:),'linewidth',2,'LineStyle', ':')
                line([GLO.line_stim(idx_shifting_windows), GLO.line_stim(idx_shifting_windows)], [0 -2-idx_shifting_windows ],'color',col(idx_shifting_windows+2,:),'linewidth',2,'LineStyle', ':')
                line([GLO.line_stim(idx_shifting_windows)+GLO.train_duration, GLO.line_stim(idx_shifting_windows)+GLO.train_duration], [0 -2-idx_shifting_windows ],'color',col(idx_shifting_windows+2,:),'linewidth',2,'LineStyle', ':')
                plot(nanmean([RT{k,1}.(fn{q}).position(pos).window(idx_shifting_windows + 1).trial{:}].*1000),-25,'^','color',col(idx_shifting_windows+2,:));
            end
            
            %         end
            xlabel('Time relative to GO [ms]', 'fontsize', GLO.fontsize-4);
            ylabel('Horizontal eye position [deg]', 'fontsize', GLO.fontsize-4);
            title(['Position ' num2str(round(target_positions(pos)))], 'fontsize', GLO.fontsize-2);
        end
        mtit(summary_1,  [plot_1_title, ' ', batch_title, ' ', fn{q}, ' batch ', num2str(k), ' plotting errors? ' num2str(errors) ], 'xoff', -0.0, 'yoff', 0.04, 'color', [0 0 0], 'fontsize', GLO.fontsize,'Interpreter', 'none');
        
        if GLO.create_pdf
            export_fig([GLO.folder_to_save batch_title, ' ', fn{q}, ' ' plot_1_title ' batch ' num2str(k) , ' plotting errors ' num2str(errors) ], '-pdf','-transparent') % pdf by run
            close(gcf)
        end
        
        plot_2_title='Raw eye traces y versus x';
        summary_2                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_2_title);
        hold on
        for pos=1:numel(target_positions)
            
            angles=[0:pi/100:2*pi];
            circle_x=cos(angles);
            circle_y=sin(angles);
            stepsize_quadrants=pi/2;
            all_phis=[-pi:stepsize_quadrants:pi];
            stepsize_plot=pi/60;
            
            for idx_shifting_windows = 1:number_of_shifting_windows+1
                for idx_trial=1: numel(eye_x{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial)
                     if idx_shifting_windows>1 && isnan(ini_ms{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}(end)*1000)
                        plot(eye_x{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}(time{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}>=0),...
                            eye_y{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}(time{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}>=0),'color','k');
                        
                     else
                         plot(eye_x{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}(time{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}>=0),...
                             eye_y{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}(time{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}>=0),'color',col(idx_shifting_windows+1,:));
                     end
                     scatter(real(sac_pos{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}),imag(sac_pos{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}),'ro');
%                     center=[real(out_comp{1}.saccades(idx_trial).tar_pos - out_comp{1}.saccades(idx_trial).fix_pos) imag(out_comp{1}.saccades(idx_trial).tar_pos - out_comp{1}.saccades(idx_trial).fix_pos)]; 
%                     plot(5*circle_x+center(1),(5*circle_y+center(2))+eye_y_offset,'-');
                end
            end
            
          
            center=[real(target_positions(pos)) imag(target_positions(pos))];
            hold on
            cl=plot(5*circle_x+center(1),(5*circle_y+center(2))+eye_y_offset,'-');
            set(cl,'Color','r','LineWidth',.5)
            cross=plot(center(1),center(2)+eye_y_offset,'MarkerSize',10,'Marker','+','MarkerEdgeColor','r','MarkerFaceColor','r');
            
        end
        
        
        axis equal
        set(gca,'xlim', [-35,35],'ylim', [5,35],'FontSize',GLO.fontsize-6)
        xlabel('Eye x position', 'fontsize', GLO.fontsize-4);
        ylabel('Eye y position', 'fontsize', GLO.fontsize-4);
        %   title(['position ' num2str(round(target_positions(pos)))], 'fontsize', GLO.fontsize + 2);
        mtit(summary_2,  [plot_2_title, ' ', batch_title, ' ', fn{q}, ' batch ', num2str(k), ' plotting errors? ' num2str(errors) ], 'xoff', -0.0, 'yoff', 0.04, 'color', [0 0 0], 'fontsize', GLO.fontsize,'Interpreter', 'none');
        
        if GLO.create_pdf
            export_fig([GLO.folder_to_save batch_title, ' ', fn{q}, ' ' plot_2_title ' batch ' num2str(k), ' plotting errors ' num2str(errors) ], '-pdf','-transparent') % pdf by run
            close(gcf)
        end
    end
    
end
a=1;

