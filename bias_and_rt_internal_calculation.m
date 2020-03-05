function [choice_bias_trial_baseline,choice_bias_stim,mean_rt,raw_rt,num_hits_temp,mean_x_acc,early_late_index,idx_stim_early,idx_stim_late,per_tar_pos,moving_average,vel,accuracy_xy,accuracy_rad,N] = bias_and_rt_internal_calculation(out_comp,target_positions)
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

for k = notemptybatches
    switch GLO.effector_to_use
        case 0
            horizontal_distance = real([out_comp{k,1}.saccades.tar_pos] - [out_comp{k,1}.saccades.fix_pos]);
            vertical_distance   = imag([out_comp{k,1}.saccades.tar_pos] - [out_comp{k,1}.saccades.fix_pos]);
            %endpositions        = [out_comp{k,1}.saccades.endpos_obs]   - [out_comp{k,1}.saccades.fix_pos];
        case 4
            horizontal_distance = real([out_comp{k,1}.reaches.tar_pos] - [out_comp{k,1}.reaches.fix_pos]);
            vertical_distance   = imag([out_comp{k,1}.reaches.tar_pos] - [out_comp{k,1}.reaches.fix_pos]);
            %endpositions        = [out_comp{k,1}.reaches.endpos_obs]   - [out_comp{k,1}.reaches.fix_pos];
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
    
    idx_choice{k}                               =   [out_comp{k,1}.binary.choice]==1 & ~ismember({out_comp{k,1}.task.abort_code},'ABORT_JAW') & [out_comp{k,1}.states.state_abo]~=2;
    idx_instructed{k}                           =   [out_comp{k,1}.binary.choice]==0 & ~ismember({out_comp{k,1}.task.abort_code},'ABORT_JAW') & [out_comp{k,1}.states.state_abo]~=2;
    
    idx_success{k}                              =   [out_comp{k,1}.binary.success]==1;
    idx_error{k}                                =   [out_comp{k,1}.binary.success]==0 & ~ismember({out_comp{k,1}.task.abort_code},'ABORT_JAW') & [out_comp{k,1}.states.state_abo]~=2;
    
    idx_stim{k}                                 =   [out_comp{k,1}.binary.microstim]==1;
    idx_baseline{k}                             =   [out_comp{k,1}.binary.microstim]==0;
    
    idx_stim_succ{k}                            =   idx_stim{k}      & idx_success{k};
    idx_baseline_succ{k}                        =   idx_baseline{k}  & idx_success{k};
    
    idx_stim_at_0{k}                            =   [out_comp{k,1}.task.stim_start]==0 & idx_stim{k};
    idx_stim_in_fix{k}                          =   [out_comp{k,1}.task.stim_state]==STATE.FIX_HOL & idx_stim{k};
    idx_stim_in_acq{k}                          =   [out_comp{k,1}.task.stim_state]==STATE.TAR_ACQ & idx_stim{k};
    idx_stim_in_cue{k}                          =   [out_comp{k,1}.task.stim_state]==STATE.CUE_ON & idx_stim{k};
    idx_stim_in_mem{k}                          =   [out_comp{k,1}.task.stim_state]==STATE.MEM_PER & idx_stim{k};
    idx_stim_in_acq_inv{k}                      =   [out_comp{k,1}.task.stim_state]==STATE.TAR_ACQ_INV & idx_stim{k};
    
    idx_L_CH_S{k}              =   idx_L{k}                                         & idx_choice{k}             & idx_stim_succ{k};
    idx_L_CH_B{k}              =   idx_L{k}                                         & idx_choice{k}             & idx_baseline_succ{k};
    idx_L_IN_S{k}              =   idx_L{k}                                         & idx_instructed{k}         & idx_stim_succ{k};
    idx_L_IN_B{k}              =   idx_L{k}                                         & idx_instructed{k}         & idx_baseline_succ{k};
    idx_R_CH_S{k}              =   idx_R{k}                                         & idx_choice{k}             & idx_stim_succ{k};
    idx_R_CH_B{k}              =   idx_R{k}                                         & idx_choice{k}             & idx_baseline_succ{k};
    idx_R_IN_S{k}              =   idx_R{k}                                         & idx_instructed{k}         & idx_stim_succ{k};
    idx_R_IN_B{k}              =   idx_R{k}                                         & idx_instructed{k}         & idx_baseline_succ{k};
    
end

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
        idx_shifting_window_succ{m}                  =   [out_comp{k,1}.task.stim_to_state_end]==abs(GLO.All_stim_windows(m))    & current_stim_state_index{k} & idx_success{k};
        idx_shifting_window{m}                       =   [out_comp{k,1}.task.stim_to_state_end]==abs(GLO.All_stim_windows(m))    & current_stim_state_index{k};
    else
        idx_shifting_window_succ{m}                  =   [out_comp{k,1}.task.stim_start]==abs(GLO.All_stim_windows(m))           & current_stim_state_index{k} & idx_success{k};
        idx_shifting_window{m}                       =   [out_comp{k,1}.task.stim_start]==abs(GLO.All_stim_windows(m))           & current_stim_state_index{k};
    end
end

idx_shifting_window_succ=idx_shifting_window_succ(GLO.windows_to_use);
idx_shifting_window     =idx_shifting_window     (GLO.windows_to_use);
k=1;
number_of_shifting_windows=numel(idx_shifting_window_succ);

% Shifting windows
for idx_shifting_windows = 1:number_of_shifting_windows
    left_chosen_successful_stimulated(idx_shifting_windows)=sum(idx_L_CH_S{k} & idx_shifting_window_succ{idx_shifting_windows});
    all_chosen_successful_stimulated(idx_shifting_windows)=sum((idx_L_CH_S{k} | idx_R_CH_S{k})  & idx_shifting_window_succ{idx_shifting_windows});
    moving_average.L_choice{idx_shifting_windows+2}=idx_L_CH_S{k} & idx_shifting_window_succ{idx_shifting_windows};
    moving_average.all_choice{idx_shifting_windows+2} = idx_choice{k} & idx_shifting_window_succ{idx_shifting_windows};
end
moving_average.L_choice{2}=idx_L_CH_B{k};
moving_average.all_choice{2}=idx_choice{k} & idx_baseline_succ{k};
moving_average.session = [out_comp{k,1}.selected.session];
moving_average.run = [out_comp{k,1}.selected.run];

% Per target position

switch GLO.effector_to_use
    case 0
        tar_positions=[out_comp{k,1}.saccades.tar_pos]-[out_comp{k,1}.saccades.fix_pos];
        nct_positions=[out_comp{k,1}.saccades.nct_pos]-[out_comp{k,1}.saccades.fix_pos];
    case 4
        tar_positions=[out_comp{k,1}.reaches.tar_pos]-[out_comp{k,1}.reaches.fix_pos];
        nct_positions=[out_comp{k,1}.reaches.nct_pos]-[out_comp{k,1}.reaches.fix_pos];
end

per_tar_pos.accuracy_xy_IN_raw    = cell(numel(target_positions),number_of_shifting_windows+2);    
per_tar_pos.accuracy_xy_CH_raw    = cell(numel(target_positions),number_of_shifting_windows+2);    
    
for pos=1:numel(target_positions)
     idx_tar_pos_base{pos}=abs(tar_positions-target_positions(pos))<=GLO.target_pos_precision & idx_baseline{k};
     idx_nct_pos_base{pos}=abs(nct_positions-target_positions(pos))<=GLO.target_pos_precision & idx_baseline{k};       
   
    %digits(old_digits); 
    idx_any_pos_base{pos}           = idx_tar_pos_base{pos} | idx_nct_pos_base{pos};
    per_tar_pos.choice_ratio(pos,2) = sum(idx_tar_pos_base{pos} & idx_choice{k}  & idx_success{k})   ./sum(idx_any_pos_base{pos} & idx_choice{k} & idx_success{k});
    
    per_tar_pos.hits_IN(pos,2)      = sum(idx_tar_pos_base{pos} & idx_success{k} & idx_instructed{k});
    per_tar_pos.hits_CH(pos,2)      = sum(idx_tar_pos_base{pos} & idx_success{k} & idx_choice{k});
    per_tar_pos.total_IN(pos,2)     = sum(idx_tar_pos_base{pos} & idx_instructed{k});
    per_tar_pos.total_CH(pos,2)     = sum(idx_tar_pos_base{pos} & idx_choice{k});
    per_tar_pos.hitrate_IN(pos,2)   = per_tar_pos.hits_IN(pos,2)./per_tar_pos.total_IN(pos,2);
    per_tar_pos.hitrate_CH(pos,2)   = per_tar_pos.hits_CH(pos,2)./per_tar_pos.total_CH(pos,2);
    
    % accuracy per position
    idx_success_pos_IN=idx_tar_pos_base{pos} & idx_success{k} & idx_instructed{k};
    idx_success_pos_CH=idx_tar_pos_base{pos} & idx_success{k} & idx_choice{k};     
    accuracy_xy_pos_IN=[out_comp{1}.saccades(idx_success_pos_IN).accuracy_xy];    
    accuracy_xy_pos_CH=[out_comp{1}.saccades(idx_success_pos_CH).accuracy_xy];
    accuracy_rad_pos_IN=[out_comp{1}.saccades(idx_success_pos_IN).accuracy_rad];    
    accuracy_rad_pos_CH=[out_comp{1}.saccades(idx_success_pos_CH).accuracy_rad];
    
    per_tar_pos.accuracy_xy_IN_mean(pos,2)   = nanmean(accuracy_xy_pos_IN);
    per_tar_pos.accuracy_xy_IN_std(pos,2)    = nanstd(real(accuracy_xy_pos_IN)) +  1i*nanstd(imag(accuracy_xy_pos_IN));
    per_tar_pos.accuracy_xy_IN_raw(pos,2)    = {accuracy_xy_pos_IN};    
    per_tar_pos.accuracy_xy_CH_mean(pos,2)   = nanmean(accuracy_xy_pos_CH);
    per_tar_pos.accuracy_xy_CH_std(pos,2)    = nanstd(real(accuracy_xy_pos_CH)) +  1i*nanstd(imag(accuracy_xy_pos_CH));
    per_tar_pos.accuracy_xy_CH_raw(pos,2)    = {accuracy_xy_pos_CH};
    
    
    per_tar_pos.accuracy_rad_IN_mean(pos,2)   = nanmean(accuracy_rad_pos_IN);
    per_tar_pos.accuracy_rad_IN_std(pos,2)    = nanstd(real(accuracy_rad_pos_IN)) +  1i*nanstd(imag(accuracy_rad_pos_IN));
    per_tar_pos.accuracy_rad_IN_raw(pos,2)    = {accuracy_rad_pos_IN};    
    per_tar_pos.accuracy_rad_CH_mean(pos,2)   = nanmean(accuracy_rad_pos_CH);
    per_tar_pos.accuracy_rad_CH_std(pos,2)    = nanstd(real(accuracy_rad_pos_CH)) +  1i*nanstd(imag(accuracy_rad_pos_CH));
    per_tar_pos.accuracy_rad_CH_raw(pos,2)    = {accuracy_rad_pos_CH};
        
    
    for idx_shifting_windows = 1:number_of_shifting_windows        
     idx_tar_pos{pos,idx_shifting_windows}=abs(tar_positions-target_positions(pos))<=GLO.target_pos_precision & idx_shifting_window{idx_shifting_windows};
     idx_nct_pos{pos,idx_shifting_windows}=abs(nct_positions-target_positions(pos))<=GLO.target_pos_precision & idx_shifting_window{idx_shifting_windows};   
     
        idx_any_pos{pos,idx_shifting_windows}=idx_tar_pos{pos,idx_shifting_windows} | idx_nct_pos{pos,idx_shifting_windows};
        per_tar_pos.choice_ratio(pos,idx_shifting_windows+2) = sum(idx_tar_pos{pos,idx_shifting_windows} & idx_choice{k} & idx_success{k})./...
            sum(idx_any_pos{pos,idx_shifting_windows} & idx_choice{k} & idx_success{k});
        per_tar_pos.hits_IN(pos,idx_shifting_windows+2)     = sum(idx_tar_pos{pos,idx_shifting_windows} & idx_success{k} & idx_instructed{k});
        per_tar_pos.hits_CH(pos,idx_shifting_windows+2)     = sum(idx_tar_pos{pos,idx_shifting_windows} & idx_success{k} & idx_choice{k});
        per_tar_pos.total_IN(pos,idx_shifting_windows+2)    = sum(idx_tar_pos{pos,idx_shifting_windows} & idx_instructed{k});
        per_tar_pos.total_CH(pos,idx_shifting_windows+2)    = sum(idx_tar_pos{pos,idx_shifting_windows} & idx_choice{k});
        
        % accuracy per position
        idx_success_pos_IN=idx_tar_pos{pos,idx_shifting_windows} & idx_success{k} & idx_instructed{k};
        idx_success_pos_CH=idx_tar_pos{pos,idx_shifting_windows} & idx_success{k} & idx_choice{k};
        accuracy_xy_pos_IN=[out_comp{1}.(sac_or_rea)(idx_success_pos_IN).accuracy_xy];
        accuracy_xy_pos_CH=[out_comp{1}.(sac_or_rea)(idx_success_pos_CH).accuracy_xy];
        accuracy_rad_pos_IN=[out_comp{1}.(sac_or_rea)(idx_success_pos_IN).accuracy_rad];
        accuracy_rad_pos_CH=[out_comp{1}.(sac_or_rea)(idx_success_pos_CH).accuracy_rad];
        
        per_tar_pos.accuracy_xy_IN_mean(pos,idx_shifting_windows+2)   = nanmean(accuracy_xy_pos_IN);
        per_tar_pos.accuracy_xy_IN_std(pos,idx_shifting_windows+2)    = nanstd(real(accuracy_xy_pos_IN)) +  1i*nanstd(imag(accuracy_xy_pos_IN));
        per_tar_pos.accuracy_xy_IN_raw(pos,idx_shifting_windows+2)    = {accuracy_xy_pos_IN};
        per_tar_pos.accuracy_xy_CH_mean(pos,idx_shifting_windows+2)   = nanmean(accuracy_xy_pos_CH);
        per_tar_pos.accuracy_xy_CH_std(pos,idx_shifting_windows+2)    = nanstd(real(accuracy_xy_pos_CH)) +  1i*nanstd(imag(accuracy_xy_pos_CH));
        per_tar_pos.accuracy_xy_CH_raw(pos,idx_shifting_windows+2)    = {accuracy_xy_pos_CH};        
        
        per_tar_pos.accuracy_rad_IN_mean(pos,idx_shifting_windows+2)   = nanmean(accuracy_rad_pos_IN);
        per_tar_pos.accuracy_rad_IN_std(pos,idx_shifting_windows+2)    = nanstd(real(accuracy_rad_pos_IN)) +  1i*nanstd(imag(accuracy_rad_pos_IN));
        per_tar_pos.accuracy_rad_IN_raw(pos,idx_shifting_windows+2)    = {accuracy_rad_pos_IN};
        per_tar_pos.accuracy_rad_CH_mean(pos,idx_shifting_windows+2)   = nanmean(accuracy_rad_pos_CH);
        per_tar_pos.accuracy_rad_CH_std(pos,idx_shifting_windows+2)    = nanstd(real(accuracy_rad_pos_CH)) +  1i*nanstd(imag(accuracy_rad_pos_CH));
        per_tar_pos.accuracy_rad_CH_raw(pos,idx_shifting_windows+2)    = {accuracy_rad_pos_CH};
        
        % maybe not necessary any more
        per_tar_pos.hitrate_IN(pos,idx_shifting_windows+2)  = per_tar_pos.hits_IN(pos,idx_shifting_windows+2)./per_tar_pos.total_IN(pos,idx_shifting_windows+2);
        per_tar_pos.hitrate_CH(pos,idx_shifting_windows+2)  = per_tar_pos.hits_CH(pos,idx_shifting_windows+2)./per_tar_pos.total_CH(pos,idx_shifting_windows+2);
    end
end

% Trial-wise baseline
choices_left_trial_baseline                     = sum(idx_L_CH_B{k});
choices_left_stim                               = left_chosen_successful_stimulated;
choices_left_stim_early                         = sum(choices_left_stim(GLO.idx_early));
choices_left_stim_late                          = sum(choices_left_stim(GLO.idx_late));
choices_all_trial_baseline                      = sum((idx_L_CH_B{k} | idx_R_CH_B{k}));
choices_all_stim                                = all_chosen_successful_stimulated;
choices_all_stim_early                          = sum(choices_all_stim(GLO.idx_early));
choices_all_stim_late                           = sum(choices_all_stim(GLO.idx_late));

Stimulation_times=1:number_of_shifting_windows;
choice_bias_trial_baseline={choices_left_trial_baseline./choices_all_trial_baseline.*100};
choice_bias_stim={choices_left_stim(Stimulation_times)./choices_all_stim(Stimulation_times).*100};
idx_stim_early={choices_left_stim_early/choices_all_stim_early*100};
idx_stim_late={choices_left_stim_late/choices_all_stim_late*100};
early_late_index.n_early_L_CH         =choices_left_stim_early;
early_late_index.n_early_R_CH         =choices_all_stim_early-choices_left_stim_early;
early_late_index.n_early_all_CH       =choices_all_stim_early;
early_late_index.n_late_L_CH          =choices_left_stim_late;
early_late_index.n_late_R_CH          =choices_all_stim_late-choices_left_stim_late;
early_late_index.n_late_all_CH        =choices_all_stim_late;


% ERROR BARS
plot_matrix={idx_L_CH_S, idx_R_CH_S,  idx_L_IN_S, idx_R_IN_S; idx_L_CH_B,  idx_R_CH_B, idx_L_IN_B, idx_R_IN_B};


for c = 1:size(plot_matrix,2)
    mean_rt(1,c)=NaN;
    raw_rt{1,c}={NaN};
    std_rt(1,c)=NaN;    
    
    vel.mean(1,c)        =NaN;
    vel.raw{1,c}         ={NaN};
    vel.std(1,c)         =NaN;
    
    accuracy_xy.mean(1,c)        =NaN + 1i*NaN;
    accuracy_xy.std(1,c)         =NaN + 1i*NaN;
    accuracy_xy.raw{1,c}         ={NaN + 1i*NaN};  
    
    accuracy_rad.mean(1,c)        =NaN + 1i*NaN;
    accuracy_rad.std(1,c)         =NaN + 1i*NaN;
    accuracy_rad.raw{1,c}         ={NaN + 1i*NaN};   
    accuracy_rad.raw_eu{1,c}      =NaN;   
    accuracy_rad.mean_eu(1,c)     =NaN;  
    accuracy_rad.std_eu(1,c)      =NaN;  
    %N.total{1,c}                  =NaN;  
    N.multistep{1,c}              =NaN;
    
    %N.total{2,c}        =[out_comp{1}.saccades(plot_matrix{2,c}{k}).n_obs];  
    N.multistep{2,c}    =[out_comp{1}.saccades(plot_matrix{2,c}{k}).n_obs]~=1;  
    RT_to_test_stim=[out_comp{1}.(sac_or_rea)(plot_matrix{2,c}{k}).lat];   
    %Velocity_to_test=[out_comp{1}.saccades(plot_matrix{2,c}{k} & [out_comp{1}.saccades.n_obs]==1).velocity];%% ???
    Velocity_to_test=[out_comp{1}.saccades(plot_matrix{2,c}{k}).velocity];%% ???
    accuracy_rad_to_test=[out_comp{1}.(sac_or_rea)(plot_matrix{2,c}{k}).accuracy_rad];
    accuracy_xy_to_test=[out_comp{1}.(sac_or_rea)(plot_matrix{2,c}{k}).accuracy_xy];
        
    accuracy_xy.mean(2,c)        =nanmean(accuracy_xy_to_test);
    accuracy_xy.std(2,c)         =nanstd(accuracy_xy_to_test);
    accuracy_xy.raw{2,c}         ={double(accuracy_xy_to_test)};
    
    accuracy_rad.mean(2,c)        =nanmean(accuracy_rad_to_test);
    accuracy_rad.std(2,c)         =nanstd(accuracy_rad_to_test);
    accuracy_rad.raw{2,c}         ={double(accuracy_rad_to_test)};    
    accuracy_rad.raw_eu{2,c}      ={abs(double(accuracy_rad_to_test))};
    accuracy_rad.mean_eu(2,c)     =nanmean(abs(accuracy_rad_to_test));
    accuracy_rad.std_eu(2,c)      =nanstd(abs(accuracy_rad_to_test));
    
    vel.mean(2,c)        =nanmean(Velocity_to_test);
    vel.std(2,c)         =nanstd(Velocity_to_test);
    vel.raw{2,c}         ={double(Velocity_to_test)};
    
    mean_rt(2,c)        =nanmean(RT_to_test_stim);
    std_rt(2,c)         =nanstd(RT_to_test_stim);
    raw_rt{2,c}         ={double(RT_to_test_stim)};
    %num_hits(c,2)       =sum(~isnan(RT_to_test_stim));
    num_hits_temp(c,1)  =0;
    num_hits_temp(c,2)  =sum([out_comp{1}.binary(plot_matrix{2,c}{k}).success]);
    
    for l=Stimulation_times
        
        %N.total(l+2,c)=[out_comp{1}.saccades(idx_shifting_window_succ{l} & plot_matrix{1,c}{k}).n_obs];  
        N.multistep{l+2,c}=[out_comp{1}.saccades(idx_shifting_window_succ{l} & plot_matrix{1,c}{k}).n_obs]~=1; 
        
        RT_to_test_stim=[out_comp{1}.(sac_or_rea)(idx_shifting_window_succ{l} & plot_matrix{1,c}{k}).lat];
        Velocity_to_test=[out_comp{1}.saccades(idx_shifting_window_succ{l} & plot_matrix{1,c}{k}).velocity];
        accuracy_rad_to_test=[out_comp{1}.(sac_or_rea)(idx_shifting_window_succ{l} & plot_matrix{1,c}{k}).accuracy_rad];
        accuracy_xy_to_test=[out_comp{1}.(sac_or_rea)(idx_shifting_window_succ{l} & plot_matrix{1,c}{k}).accuracy_xy];
         
        accuracy_xy.mean(l+2,c)        =nanmean(accuracy_xy_to_test);
        accuracy_xy.std(l+2,c)         =nanstd(accuracy_xy_to_test);
        accuracy_xy.raw{l+2,c}         ={double(accuracy_xy_to_test)};
        
        accuracy_rad.mean(l+2,c)        =nanmean(accuracy_rad_to_test);
        accuracy_rad.std(l+2,c)         =nanstd(accuracy_rad_to_test);
        accuracy_rad.raw{l+2,c}         ={double(accuracy_rad_to_test)};
        accuracy_rad.raw_eu{l+2,c}      ={abs(double(accuracy_rad_to_test))};
        accuracy_rad.mean_eu(l+2,c)     =nanmean(abs(accuracy_rad_to_test));
        accuracy_rad.std_eu(l+2,c)      =nanstd(abs(accuracy_rad_to_test));
        
        vel.mean(l+2,c)        =nanmean(Velocity_to_test);
        vel.std(l+2,c)         =nanstd(Velocity_to_test);
        vel.raw{l+2,c}         ={double(Velocity_to_test)};
        
        mean_rt(l+2,c)          =nanmean(RT_to_test_stim);
        std_rt(l+2,c)           =nanstd(RT_to_test_stim);
        raw_rt{l+2,c}           ={double(RT_to_test_stim)};
        %num_hits(c,l+2)         =sum(~isnan(RT_to_test_stim));
        num_hits_temp(c,l+2)    =sum([out_comp{1}.binary(idx_shifting_window{l} & plot_matrix{1,c}{k}).success]);
    end
end
clear variable_to_test_stim




for c = 1:size(plot_matrix,2)
    x_acc{1,c}=                    NaN;
    mean_x_acc(1,c)=               NaN;
    std_x_acc(1,c)=                NaN;
    y_acc{1,c}=                    NaN;
    mean_y_acc(1,c)=               NaN;
    std_y_acc(1,c)=                NaN;
    
    switch GLO.effector_to_use
        case 0
%             if rem(c,2) & sum(plot_matrix{2,c}{k}) ~= 0
%                 variable_to_test_stim{2,c}=         [out_comp{1}.saccades(plot_matrix{2,c}{k}).precision_pos_l];
%             else
%                 variable_to_test_stim{2,c}=         [out_comp{1}.saccades(plot_matrix{2,c}{k}).precision_pos_r];
%             end
             variable_to_test_stim{2,c}=         [out_comp{1}.saccades(plot_matrix{2,c}{k}).accuracy_xy];
        case 4
%             if rem(c,2) & sum(plot_matrix{2,c}{k}) ~= 0
%                 variable_to_test_stim{2,c}=         [out_comp{1}.reaches(plot_matrix{2,c}{k}).precision_pos_l];
%             else
%                 variable_to_test_stim{2,c}=         [out_comp{1}.reaches(plot_matrix{2,c}{k}).precision_pos_r];
%             end
             variable_to_test_stim{2,c}=         [out_comp{1}.reaches(plot_matrix{2,c}{k}).accuracy_xy];
    end
    
    if ~isempty(variable_to_test_stim{2,c})
        x_acc{2,c}=                    real(variable_to_test_stim{2,c}).';
        mean_x_acc(2,c)=               nanmean(x_acc{2,c});
        std_x_acc(2,c)=                nanstd(x_acc{2,c});
        y_acc{2,c}=                    imag(variable_to_test_stim{2,c}).';
        mean_y_acc(2,c)=               nanmean(y_acc{2,c});
        std_y_acc(2,c)=                nanstd(y_acc{2,c});
    else
        x_acc{2,c}=                    NaN;
        mean_x_acc(2,c)=               NaN;
        std_x_acc(2,c)=                NaN;
        y_acc{2,c}=                    NaN;
        mean_y_acc(2,c)=               NaN;
        std_y_acc(2,c)=                NaN;
    end
    
    for l=Stimulation_times
        switch GLO.effector_to_use
            case 0
%                 if rem(c,2) & sum(plot_matrix{1,c}{1}) ~= 0 & sum(idx_shifting_window_succ{l}) ~= 0
%                     variable_to_test_stim{l+2,c}=         [out_comp{1}.saccades(idx_shifting_window_succ{l} & plot_matrix{1,c}{k}).precision_pos_l];
%                 else
%                     variable_to_test_stim{l+2,c}=         [out_comp{1}.saccades(idx_shifting_window_succ{l} & plot_matrix{1,c}{k}).precision_pos_r];
%                 end
                  variable_to_test_stim{l+2,c}=         [out_comp{1}.saccades(idx_shifting_window_succ{l} & plot_matrix{1,c}{k}).accuracy_xy];
            case 4
%                 if rem(c,2) & sum(plot_matrix{1,c}{1}) ~= 0 & sum(idx_shifting_window_succ{l}) ~= 0
%                     variable_to_test_stim{l+2,c}=         [out_comp{1}.reaches(idx_shifting_window_succ{l} & plot_matrix{1,c}{k}).precision_pos_l];
%                 else
%                     variable_to_test_stim{l+2,c}=         [out_comp{1}.reaches(idx_shifting_window_succ{l} & plot_matrix{1,c}{k}).precision_pos_r];
%                 end
                variable_to_test_stim{l+2,c}=         [out_comp{1}.reaches(idx_shifting_window_succ{l} & plot_matrix{1,c}{k}).accuracy_xy];
        end
        
        if ~isempty(variable_to_test_stim{l+2,c})
            x_acc{l+2,c}=                    real(variable_to_test_stim{l+2,c}).';
            mean_x_acc(l+2,c)=               nanmean(x_acc{l+2,c});
            std_x_acc(l+2,c)=                nanstd(x_acc{l+2,c});
            y_acc{l+2,c}=                    imag(variable_to_test_stim{l+2,c}).';
            mean_y_acc(l+2,c)=               nanmean(y_acc{l+2,c});
            std_y_acc(l+2,c)=                nanstd(y_acc{l+2,c});
        else
            x_acc{l+2,c}=                    NaN;
            mean_x_acc(l+2,c)=               NaN;
            std_x_acc(l+2,c)=                NaN;
            y_acc{l+2,c}=                    NaN;
            mean_y_acc(l+2,c)=               NaN;
            std_y_acc(l+2,c)=                NaN;
        end
        
    end
end

a=1;

