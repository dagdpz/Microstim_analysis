function batch_calculations(batch_out,monkey,batch_title)
global  GLO

windows                                                         = [0 0 GLO.All_stim_windows(GLO.windows_to_use)];
Conditions                                                      = {'Choice_Left','Choice_Right','Instructed_Left','Instructed_Right'};
positions                                                       = batch_out.target_positions;

n_sessions                                                      = numel(batch_out.raw_rt);
n_windows                                                       = numel(windows);
n_conditions                                                    = numel(Conditions);
n_positions                                                     = numel(positions);

RT_raw                                                          = cell(size(batch_out.raw_rt{1},1),size(batch_out.raw_rt{1},2));
Vel_raw                                                         = cell(size(batch_out.velocities{1}.raw,1),size(batch_out.velocities{1}.raw,2));
Accuracy_raw                                                    = cell(size(batch_out.accuracy_rad{1}.raw,1),size(batch_out.accuracy_rad{1}.raw,2));
RT_raw_first_mode                                               = num2cell(NaN(n_windows, n_conditions));
RT_raw_second_mode                                              = num2cell(NaN(n_windows, n_conditions));
% Moving_average_L_CH                                             = num2cell(NaN(n_windows, n_sessions));
% Moving_average_all_CH                                           = num2cell(NaN(n_windows, n_sessions));
Hit_rate_IN_per_target                                          = cell(n_positions, n_windows);
Hit_rate_CH_per_target                                          = cell(n_positions, n_windows);
Choice_rate_per_target                                          = cell(n_positions, n_windows);
Accuracy_IN_per_target_rad                                      = cell(n_positions, n_windows);
Accuracy_CH_per_target_rad                                      = cell(n_positions, n_windows);
Accuracy_IN_per_target_xy                                       = cell(n_positions, n_windows);
Accuracy_CH_per_target_xy                                       = cell(n_positions, n_windows);

Accuracy_IN_per_target_xy_raw                                   = cell(n_positions, n_windows);
Accuracy_CH_per_target_xy_raw                                   = cell(n_positions, n_windows);
Accuracy_IN_per_target_rad_raw                                  = cell(n_positions, n_windows);
Accuracy_CH_per_target_rad_raw                                  = cell(n_positions, n_windows);

Precision_IN_per_target_rad                                         = cell(n_positions, n_windows);
Precision_CH_per_target_rad                                         = cell(n_positions, n_windows);

run_counter=0;

%% Reassignments (including adding run baseline as first window)
for s=1:n_sessions
    
    % bias
    bias_per_session_run_baseline(s,:)                          = batch_out.bias_B_run{s}{:};
    bias_per_session_trial_baseline(s,:)                        = batch_out.bias_B_trial{s}{:};
    bias_per_session_stimulated(s,:)                            = [batch_out.bias_S{s}{:}];
    bias_per_session(s,:)                                       = [bias_per_session_run_baseline(s) bias_per_session_trial_baseline(s) bias_per_session_stimulated(s,:)];
    bias_per_session_stimulated_early(s,:)                      = batch_out.idx_stim_early{s}{:};
    bias_per_session_stimulated_late(s,:)                       = batch_out.idx_stim_late{s}{:};
    bias_diff_stimulated_early(s,:)                             = bias_per_session_stimulated_early(s) - bias_per_session_trial_baseline(s);
    bias_diff_stimulated_late(s,:)                              = bias_per_session_stimulated_late(s) - bias_per_session_trial_baseline(s);
    
    % adding run baseline as the first window
%     batch_out.per_tar_pos{s}.hits_IN(:,1)                       = batch_out.per_tar_pos_base{s}.hits_IN(:,2);
%     batch_out.per_tar_pos{s}.hits_CH(:,1)                       = batch_out.per_tar_pos_base{s}.hits_CH(:,2);
%     batch_out.per_tar_pos{s}.total_IN(:,1)                      = batch_out.per_tar_pos_base{s}.total_IN(:,2);
%     batch_out.per_tar_pos{s}.total_CH(:,1)                      = batch_out.per_tar_pos_base{s}.total_CH(:,2);
%     batch_out.per_tar_pos{s}.hitrate_IN(:,1)                    = batch_out.per_tar_pos_base{s}.hitrate_IN(:,2);
%     batch_out.per_tar_pos{s}.hitrate_CH(:,1)                    = batch_out.per_tar_pos_base{s}.hitrate_CH(:,2);
%     batch_out.per_tar_pos{s}.choice_ratio(:,1)                  = batch_out.per_tar_pos_base{s}.choice_ratio(:,2);
%     batch_out.per_tar_pos{s}.accuracy_xy_IN_mean(:,1)           = batch_out.per_tar_pos_base{s}.accuracy_xy_IN_mean(:,2);
%     batch_out.per_tar_pos{s}.accuracy_xy_IN_std(:,1)            = batch_out.per_tar_pos_base{s}.accuracy_xy_IN_std(:,2);
%     batch_out.per_tar_pos{s}.accuracy_xy_IN_raw(:,1)            = batch_out.per_tar_pos_base{s}.accuracy_xy_IN_raw(:,2);
%     batch_out.per_tar_pos{s}.accuracy_xy_CH_mean(:,1)           = batch_out.per_tar_pos_base{s}.accuracy_xy_CH_mean(:,2);
%     batch_out.per_tar_pos{s}.accuracy_xy_CH_std(:,1)            = batch_out.per_tar_pos_base{s}.accuracy_xy_CH_std(:,2);
%     batch_out.per_tar_pos{s}.accuracy_xy_CH_raw(:,1)            = batch_out.per_tar_pos_base{s}.accuracy_xy_CH_raw(:,2);
        
    batch_out.moving_average_t{s}.L_choice(1)                   = batch_out.moving_average_r{s}.L_choice(2);
    batch_out.moving_average_t{s}.all_choice(1)                 = batch_out.moving_average_r{s}.all_choice(2);
    
    batch_out.velocities{s}.mean(1,:)                           =batch_out.velocities_B_run{s}.mean(2,:);    
    batch_out.accuracy_rad{s}.mean(1,:)                         =batch_out.accuracy_rad_B_run{s}.mean(2,:); 
    batch_out.accuracy_rad{s}.mean_eu(1,:)                      =batch_out.accuracy_rad_B_run{s}.mean_eu(2,:);
    batch_out.accuracy_xy{s}.mean(1,:)                          =batch_out.accuracy_xy_B_run{s}.mean(2,:); 
    batch_out.accuracy_xy{s}.std(1,:)                           =batch_out.accuracy_xy_B_run{s}.std(2,:); 
    batch_out.N{s}.multistep(1,:)                               =batch_out.N_B_run{s}.multistep(2,:); 
    %batch_out.N{s}.total(1,:)                                   =batch_out.N_B_run{s}.total(2,:); 
    
    %batch_out.accuracy_xy{s}.mean_eu(1,:)                       =batch_out.accuracy_xy_B_run{s}.mean_eu(2,:);
    
    
    batch_out.mean_rt{s}(1,:)                                   =batch_out.mean_rt_B_run{s}(2,:);
    batch_out.mean_x_acc{s}(1,:)                                =batch_out.mean_x_acc_B_run{s}(2,:);    
    
    session_runs=[batch_out.moving_average_t{s}.session; batch_out.moving_average_t{s}.run]';
    unique_session_runs=unique(session_runs,'rows');
    n_runs_in_current_batch=size(unique_session_runs,1);
        
    for w=1:n_windows
        
        % cutting trials in batches according to run
        if w>1
            run_in_batch_start=1;
            for r=1:n_runs_in_current_batch
                run_in_batch_end=sum(ismember(session_runs,unique_session_runs(r,:),'rows'))+run_in_batch_start-1;
                Moving_average_L_CH_per_run{w, run_counter+r}                               = batch_out.moving_average_t{s}.L_choice {w}(run_in_batch_start:run_in_batch_end);
                Moving_average_all_CH_per_run{w, run_counter+r}                             = batch_out.moving_average_t{s}.all_choice {w}(run_in_batch_start:run_in_batch_end);
                run_in_batch_start=run_in_batch_end+1;
            end
        end
        
        Moving_average_L_CH{w,s}                               = batch_out.moving_average_t{s}.L_choice {w};
        Moving_average_all_CH{w,s}                             = batch_out.moving_average_t{s}.all_choice {w};
        
        
        Hits_IN_L(s,w)            = 0;
        Hits_CH_L(s,w)            = 0;
        Total_trials_IN_L(s,w)    = 0;
        Total_trials_CH_L(s,w)    = 0;
        Hits_IN_R(s,w)            = 0;
        Hits_CH_R(s,w)            = 0;
        Total_trials_IN_R(s,w)    = 0;
        Total_trials_CH_R(s,w)    = 0;
        
        for c=1:n_conditions
            batch_out.raw_rt{s}{1,c}                = batch_out.raw_rt_B_run{s}{2,c};
            batch_out.velocities{s}.raw{1,c}        = batch_out.velocities_B_run{s}.raw{2,c};
            batch_out.accuracy_rad{s}.raw{1,c}      = batch_out.accuracy_rad_B_run{s}.raw{2,c};
            batch_out.accuracy_xy{s}.raw{1,c}       = batch_out.accuracy_xy_B_run{s}.raw{2,c};
            RT_raw{w,c}                             = [RT_raw{w,c}, batch_out.raw_rt{s}{w,c}{:}];
            Vel_raw{w,c}                            = [Vel_raw{w,c}, batch_out.velocities{s}.raw{w,c}{:}];
            Accuracy_raw{w,c}                       = [Accuracy_raw{w,c}, batch_out.accuracy_rad{s}.raw{w,c}{:}];
            
            RT_per_session_all{w,c}(s,:)            = batch_out.mean_rt{s}(w,c);
            
            Accuracy_eu_per_session_all{w,c}(s,:)   = batch_out.accuracy_rad{s}.mean_eu(w,c);
            accuracy_S_TB{w,c}(s,:)                 = batch_out.mean_x_acc{s}(w,c);
            
            
            % first and second mode seperation
            temp_vel_raw     =[batch_out.velocities{s}.raw{w,c}{:}];
            temp_RT_raw      =batch_out.raw_rt{s}{w,c}{:};            
            idx_second_mode{w,c}         =temp_RT_raw>windows(w)+0.2;
            if w < 3
                idx_second_mode{w,c}     = true(size(temp_RT_raw));
                idx_first_mode{w,c}      = idx_second_mode{w,c};
            else
                idx_first_mode{w,c}  =~idx_second_mode{w,c};
            end
            if ~isempty(idx_second_mode{w,c})
                RT_raw_first_mode{w,c}     =[RT_raw_first_mode{w,c}  temp_RT_raw(idx_first_mode{w,c})];
                RT_raw_second_mode{w,c}    =[RT_raw_second_mode{w,c} temp_RT_raw(idx_second_mode{w,c})];
            end
            
            Fraction_multistep{w,c}(s,:)            =sum(batch_out.N{s}.multistep{w,c})/numel(batch_out.N{s}.multistep{w,c});
            Fraction_second_mode{w,c}(s,:)          =sum(idx_second_mode{w,c})/(sum(idx_first_mode{w,c}) + sum(idx_second_mode{w,c}));
            
            if ~isempty(nanmean(temp_RT_raw(idx_first_mode{w,c})))
                RT_mean_first_mode{w,c}(s,:)        =   nanmean(temp_RT_raw(idx_first_mode{w,c}));
                Vel_mean_first_mode{w,c}(s,:)       =   nanmean(temp_vel_raw(idx_first_mode{w,c} & ~[batch_out.N{s}.multistep{w,c}]));
                N_first_mode{w,c}(s,:)              =   sum(idx_first_mode{w,c});
            else
                RT_mean_first_mode{w,c}(s,:)        =   NaN;
                Vel_mean_first_mode{w,c}(s,:)       =   NaN;
                N_first_mode{w,c}(s,:)              =   0;
            end
            if ~isempty(nanmean(temp_RT_raw(idx_second_mode{w,c})))
                RT_mean_second_mode{w,c}(s,:)       =   nanmean(temp_RT_raw(idx_second_mode{w,c}));
                Vel_mean_second_mode{w,c}(s,:)      =   nanmean(temp_vel_raw(idx_second_mode{w,c} & ~[batch_out.N{s}.multistep{w,c}]));
                N_second_mode{w,c}(s,:)             =   sum(idx_second_mode{w,c});
            else
                RT_mean_second_mode{w,c}(s,:)       =   NaN;
                Vel_mean_second_mode{w,c}(s,:)      =   NaN;
                N_second_mode{w,c}(s,:)             =   0;
            end
%             if ~isempty(nanmean(temp_vel_raw(~[batch_out.N{s}.multistep{w,c}])))
                Vel_per_session_all{w,c}(s,:)       =  nanmean(temp_vel_raw(~[batch_out.N{s}.multistep{w,c}]));
                
            RT_per_session_diff{w,c}(s,:)           = RT_per_session_all{w,c}(s,:) - batch_out.mean_rt{s}(2,c);            
            RT_mean_first_mode_diff{w,c}(s,:)       = RT_mean_first_mode{w,c}(s,:) - batch_out.mean_rt{s}(2,c);            
            RT_mean_second_mode_diff{w,c}(s,:)      = RT_mean_second_mode{w,c}(s,:) - batch_out.mean_rt{s}(2,c);
%             else
%                 Vel_per_session_all{w,c}(s,:)           =    NaN;
%             end
        end
        bias_per_session_first_mode(s,w)= sum(idx_first_mode{w,1})/(sum(idx_first_mode{w,2}) + sum(idx_first_mode{w,1}))*100;
        bias_per_session_second_mode(s,w)= sum(idx_second_mode{w,1})/(sum(idx_second_mode{w,2}) + sum(idx_second_mode{w,1}))*100;
        
        
        for p=1:n_positions
            Hit_rate_IN_per_target{p,w}     =[Hit_rate_IN_per_target{p,w}         batch_out.per_tar_pos{s}.hitrate_IN(p,w)];
            Hit_rate_CH_per_target{p,w}     =[Hit_rate_CH_per_target{p,w}         batch_out.per_tar_pos{s}.hitrate_CH(p,w)];
            Choice_rate_per_target{p,w}     =[Choice_rate_per_target{p,w}         batch_out.per_tar_pos{s}.choice_ratio(p,w)];
            
%             Accuracy_IN_per_target_rad{p,w}     =[Accuracy_IN_per_target_rad{p,w}         batch_out.per_tar_pos{s}.accuracy_xy_IN_mean(p,w)];
%             Accuracy_CH_per_target_rad{p,w}     =[Accuracy_CH_per_target_rad{p,w}         batch_out.per_tar_pos{s}.accuracy_xy_CH_mean(p,w)];
%             Accuracy_IN_per_target_raw{p,w} =[Accuracy_IN_per_target_raw{p,w}     batch_out.per_tar_pos{s}.accuracy_xy_IN_raw{p,w}];
%             Accuracy_CH_per_target_raw{p,w} =[Accuracy_CH_per_target_raw{p,w}     batch_out.per_tar_pos{s}.accuracy_xy_CH_raw{p,w}];
%             Precision_IN_per_target{p,w}    =[Precision_IN_per_target{p,w}        batch_out.per_tar_pos{s}.accuracy_xy_IN_std(p,w)];
%             Precision_CH_per_target{p,w}    =[Precision_CH_per_target{p,w}        batch_out.per_tar_pos{s}.accuracy_xy_CH_std(p,w)];
                        
            Accuracy_IN_per_target_rad{p,w} =[Accuracy_IN_per_target_rad{p,w}     batch_out.per_tar_pos{s}.accuracy_rad_IN_mean(p,w)];
            Accuracy_CH_per_target_rad{p,w} =[Accuracy_CH_per_target_rad{p,w}     batch_out.per_tar_pos{s}.accuracy_rad_CH_mean(p,w)];
            
            
            Accuracy_IN_per_target_rad_raw{p,w} =[Accuracy_IN_per_target_rad_raw{p,w}     batch_out.per_tar_pos{s}.accuracy_rad_IN_raw{p,w}];
            Accuracy_CH_per_target_rad_raw{p,w} =[Accuracy_CH_per_target_rad_raw{p,w}     batch_out.per_tar_pos{s}.accuracy_rad_CH_raw{p,w}];
                        
            Accuracy_IN_per_target_xy{p,w} =[Accuracy_IN_per_target_xy{p,w}     batch_out.per_tar_pos{s}.accuracy_xy_IN_mean(p,w)];
            Accuracy_CH_per_target_xy{p,w} =[Accuracy_CH_per_target_xy{p,w}     batch_out.per_tar_pos{s}.accuracy_xy_CH_mean(p,w)];
            Accuracy_IN_per_target_xy_raw{p,w} =[Accuracy_IN_per_target_xy_raw{p,w}     batch_out.per_tar_pos{s}.accuracy_xy_IN_raw{p,w}];
            Accuracy_CH_per_target_xy_raw{p,w} =[Accuracy_CH_per_target_xy_raw{p,w}     batch_out.per_tar_pos{s}.accuracy_xy_CH_raw{p,w}];
            
            Precision_IN_per_target_rad{p,w}    =[Precision_IN_per_target_rad{p,w}        batch_out.per_tar_pos{s}.accuracy_rad_IN_std(p,w)];
            Precision_CH_per_target_rad{p,w}    =[Precision_CH_per_target_rad{p,w}        batch_out.per_tar_pos{s}.accuracy_rad_CH_std(p,w)];
            
            if real(positions(p))<0
                Hits_IN_L(s,w)            = Hits_IN_L(s,w) + batch_out.per_tar_pos{s}.hits_IN(p,w);
                Hits_CH_L(s,w)            = Hits_CH_L(s,w) + batch_out.per_tar_pos{s}.hits_CH(p,w);
                Total_trials_IN_L(s,w)    = Total_trials_IN_L(s,w) + batch_out.per_tar_pos{s}.total_IN(p,w);
                Total_trials_CH_L(s,w)    = Total_trials_CH_L(s,w) + batch_out.per_tar_pos{s}.total_CH(p,w);
            elseif real(positions(p))>0
                Hits_IN_R(s,w)            = Hits_IN_R(s,w) + batch_out.per_tar_pos{s}.hits_IN(p,w);
                Hits_CH_R(s,w)            = Hits_CH_R(s,w) + batch_out.per_tar_pos{s}.hits_CH(p,w);
                Total_trials_IN_R(s,w)    = Total_trials_IN_R(s,w) + batch_out.per_tar_pos{s}.total_IN(p,w);
                Total_trials_CH_R(s,w)    = Total_trials_CH_R(s,w) + batch_out.per_tar_pos{s}.total_CH(p,w);
            end
        end
    end
    run_counter=run_counter+n_runs_in_current_batch;
    
        bias_per_session_diff(s,:)                                  = bias_per_session(s,:) - bias_per_session(s,2) ;
        bias_per_session_first_mode_diff(s,:)                       = bias_per_session_first_mode(s,:) - bias_per_session_first_mode(s,2) ;
        bias_per_session_second_mode_diff(s,:)                       = bias_per_session_second_mode(s,:) - bias_per_session_second_mode(s,2) ;
end

%% means and standard deviations (and reformatting seconds to ms, probability to %,...)
for w=1:n_windows    
    %% hitrate for choices is now combined for L and R
    Hit_rate.per_session.L_IN(:,w)    = (Hits_IN_L(:,w)./Total_trials_IN_L(:,w)).*100;
    Hit_rate.per_session.L_CH(:,w)    = (Hits_CH_L(:,w)+Hits_CH_R(:,w))./(Total_trials_CH_L(:,w)+Total_trials_CH_R(:,w)).*100;
    Hit_rate.per_session.R_IN(:,w)    = (Hits_IN_R(:,w)./Total_trials_IN_R(:,w)).*100;
    Hit_rate.per_session.R_CH(:,w)    = (Hits_CH_L(:,w)+Hits_CH_R(:,w))./(Total_trials_CH_L(:,w)+Total_trials_CH_R(:,w)).*100;
    
    Hit_rate.mean.L_IN(w)    = nanmean(Hits_IN_L(:,w)./Total_trials_IN_L(:,w)).*100;
    Hit_rate.mean.L_CH(w)    = nanmean((Hits_CH_L(:,w)+Hits_CH_R(:,w))./(Total_trials_CH_L(:,w)+Total_trials_CH_R(:,w))).*100;
    Hit_rate.mean.R_IN(w)    = nanmean(Hits_IN_R(:,w)./Total_trials_IN_R(:,w)).*100;
    Hit_rate.mean.R_CH(w)    = nanmean((Hits_CH_L(:,w)+Hits_CH_R(:,w))./(Total_trials_CH_L(:,w)+Total_trials_CH_R(:,w))).*100;
    
    Hit_rate.sem.L_IN(w)     = sterr(Hits_IN_L(:,w)./Total_trials_IN_L(:,w)).*100;
    Hit_rate.sem.L_CH(w)     = sterr((Hits_CH_L(:,w)+Hits_CH_R(:,w))./(Total_trials_CH_L(:,w)+Total_trials_CH_R(:,w))).*100;
    Hit_rate.sem.R_IN(w)     = sterr(Hits_IN_R(:,w)./Total_trials_IN_R(:,w)).*100;
    Hit_rate.sem.R_CH(w)     = sterr((Hits_CH_L(:,w)+Hits_CH_R(:,w))./(Total_trials_CH_L(:,w)+Total_trials_CH_R(:,w))).*100;
    
    bias.mean.first_mode(w)     =nanmean(bias_per_session_first_mode(:,w));
    bias.mean.second_mode(w)    =nanmean(bias_per_session_second_mode(:,w));
    bias.sem.first_mode(w)      =sterr(bias_per_session_first_mode(:,w));
    bias.sem.second_mode(w)     =sterr(bias_per_session_second_mode(:,w));
    
    bias.per_session.all(:,w)                                            = bias_per_session(:,w) ;
    bias.mean.all(w)                                                   = nanmean(bias_per_session(:,w)) ;
    bias.sem.all(w)                                                    = sterr(bias_per_session(:,w)) ;
    bias.mean.all_diff(w)                                              = nanmean(bias_per_session_diff(:,w));
    bias.mean.first_mode_diff(w)                                       = nanmean(bias_per_session_first_mode_diff(:,w));
    bias.mean.second_mode_diff(w)                                      = nanmean(bias_per_session_second_mode_diff(:,w));
    
    
    
    bias.sem.all_diff(w)                                               = sterr(bias_per_session_diff(:,w));
    bias.sem.first_mode_diff(w)                                        = sterr(bias_per_session_first_mode_diff(:,w));
    bias.sem.second_mode_diff(w)                                       = sterr(bias_per_session_second_mode_diff(:,w));
    
    Accuracy.mean.all(w)                                               = nanmean(bias_per_session(:,w)) ;
    Accuracy.sem.all(w)                                                = nanmean(bias_per_session(:,w)) ;
    
    for c=1:n_conditions
        RT_per_session.mean.all(w,c)                = nanmean(RT_per_session_all{w,c})*1000;
        RT_per_session.mean.first_mode(w,c)         = nanmean(RT_mean_first_mode{w,c})*1000;
        RT_per_session.mean.second_mode(w,c)        = nanmean(RT_mean_second_mode{w,c})*1000;
        RT_per_session.sem.all(w,c)                 = sterr(RT_per_session_all{w,c})*1000;
        RT_per_session.sem.first_mode(w,c)          = sterr(RT_mean_first_mode{w,c})*1000;
        RT_per_session.sem.second_mode(w,c)         = sterr(RT_mean_second_mode{w,c})*1000; 
        RT_per_session.mean.all_diff(w,c)           = nanmean(RT_per_session_diff{w,c})*1000;
        RT_per_session.mean.first_mode_diff(w,c)    = nanmean(RT_mean_first_mode_diff{w,c})*1000;
        RT_per_session.mean.second_mode_diff(w,c)   = nanmean(RT_mean_second_mode_diff{w,c})*1000;        
        RT_per_session.sem.all_diff(w,c)            = sterr(RT_per_session_diff{w,c})*1000; 
        RT_per_session.sem.first_mode_diff(w,c)     = sterr(RT_mean_first_mode_diff{w,c})*1000; 
        RT_per_session.sem.second_mode_diff(w,c)    = sterr(RT_mean_second_mode_diff{w,c})*1000; 
        
        RT_raw_mean(w,c)                        = nanmean(RT_raw{w,c})*1000;
        RT_raw_std(w,c)                         = nanstd(RT_raw{w,c})*1000;        
        RT_raw_mean_first_mode(w,c)             = nanmean(RT_raw_first_mode{w,c})*1000;
        RT_raw_std_first_mode(w,c)              = nanstd(RT_raw_first_mode{w,c})*1000;        
        RT_raw_mean_second_mode(w,c)            = nanmean(RT_raw_second_mode{w,c})*1000;
        RT_raw_std_second_mode(w,c)             = nanstd(RT_raw_second_mode{w,c})*1000;
        
        Vel_per_session.mean.all(w,c)            = nanmean(Vel_per_session_all{w,c});
        Vel_per_session.mean.first_mode(w,c)     = nanmean(Vel_mean_first_mode{w,c});
        Vel_per_session.mean.second_mode(w,c)    = nanmean(Vel_mean_second_mode{w,c});
        Vel_per_session.sem.all(w,c)             = sterr(Vel_per_session_all{w,c});
        Vel_per_session.sem.first_mode(w,c)      = sterr(Vel_mean_first_mode{w,c});
        Vel_per_session.sem.second_mode(w,c)     = sterr(Vel_mean_second_mode{w,c}); 
        
        Fraction_multistep_per_session.mean.all(w,c)=nanmean( Fraction_multistep{w,c});
        Fraction_multistep_per_session.sem.all(w,c)=sterr( Fraction_multistep{w,c});
        
        accuracy_eu_per_session.mean.all(w,c)   =nanmean(Accuracy_eu_per_session_all{w,c});
        accuracy_eu_per_session.sem.all(w,c)    =sterr(Accuracy_eu_per_session_all{w,c});
        
        Bins.RThist=[0.08:0.01:0.4];
        RThist{c}(:,w)                          = hist(RT_raw{w,c},Bins.RThist);
        RThist{c}(:,w)                          = RThist{c}(:,w)./sum(RThist{c}(:,w))*100;
        
        Bins.Velhist=[0:100:1500];
        Velhist{c}(:,w)                         = hist(Vel_raw{w,c},Bins.Velhist);
        Velhist{c}(:,w)                         = Velhist{c}(:,w)./sum(Velhist{c}(:,w))*100;
                
        Bins.Acchist=[0:0.4:8];
        Acchist{c}(:,w)                         = hist(abs(Accuracy_raw{w,c}),Bins.Acchist);
        Acchist{c}(:,w)                         = Acchist{c}(:,w)./sum(Acchist{c}(:,w))*100;
        
%         accuracy_mean_S_TB(w,c)                 = nanmean(accuracy_S_TB{w,c});
%         accuracy_sem_S_TB(w,c)                  = sterr(accuracy_S_TB{w,c});
    end
    for p=1:n_positions
        Hit_rate_IN_per_target_mean(p,w)     =nanmean(Hit_rate_IN_per_target{p,w});
        Hit_rate_CH_per_target_mean(p,w)     =nanmean(Hit_rate_CH_per_target{p,w});
        Choice_rate_per_target_mean(p,w)     =nanmean(Choice_rate_per_target{p,w});
    end
    
    %Conditions={'Choice_Left','Choice_Right','Instructed_Left','Instructed_Right'};
    condition_difference_labels                                           = {'R_IN_L_IN','R_CH_L_CH','R_CH_R_IN','L_CH_L_IN'};
    condition_difference_indexes                                          = [4,3;2,1;2,4;1,3];
    for d=1:numel(condition_difference_labels)
        
        RT_per_session_condition_differences_all{w,d}                     =RT_per_session_all{w,condition_difference_indexes(d,1)}        - RT_per_session_all{w,condition_difference_indexes(d,2)};
        RT_per_session_condition_differences_first_mode{w,d}          =RT_mean_first_mode{w,condition_difference_indexes(d,1)}    - RT_mean_first_mode{w,condition_difference_indexes(d,2)};
        RT_per_session_condition_differences_second_mode{w,d}         =RT_mean_second_mode{w,condition_difference_indexes(d,1)}   - RT_mean_second_mode{w,condition_difference_indexes(d,2)};
        
        RT_per_session_condition_differences.mean.all(w,d)                =nanmean(RT_per_session_all{w,condition_difference_indexes(d,1)}        - RT_per_session_all{w,condition_difference_indexes(d,2)})*1000;
        RT_per_session_condition_differences.mean.first_mode(w,d)     =nanmean(RT_mean_first_mode{w,condition_difference_indexes(d,1)}    - RT_mean_first_mode{w,condition_difference_indexes(d,2)})*1000;
        RT_per_session_condition_differences.mean.second_mode(w,d)    =nanmean(RT_mean_second_mode{w,condition_difference_indexes(d,1)}   - RT_mean_second_mode{w,condition_difference_indexes(d,2)})*1000;
        
        RT_per_session_condition_differences.sem.all(w,d)                 =sterr(RT_per_session_all{w,condition_difference_indexes(d,1)}         - RT_per_session_all{w,condition_difference_indexes(d,2)})*1000;
        RT_per_session_condition_differences.sem.first_mode(w,d)      =sterr(RT_mean_first_mode{w,condition_difference_indexes(d,1)}     - RT_mean_first_mode{w,condition_difference_indexes(d,2)})*1000;
        RT_per_session_condition_differences.sem.second_mode(w,d)     =sterr(RT_mean_second_mode{w,condition_difference_indexes(d,1)}    - RT_mean_second_mode{w,condition_difference_indexes(d,2)})*1000;
    end
end

for condition=1:4
    temp_RT_in_condition=[RT_per_session_all{3:end,condition}];
    temp_RT_diff_in_condition=[RT_per_session_diff{3:end,condition}]; 
    RT_per_session_early{condition}                 = temp_RT_in_condition(:,GLO.idx_early);
    RT_per_session_late{condition}                  = temp_RT_in_condition(:,GLO.idx_late);
    RT_per_session_early_diff{condition}            = temp_RT_diff_in_condition(:,GLO.idx_early);
    RT_per_session_late_diff{condition}             = temp_RT_diff_in_condition(:,GLO.idx_late);
    
    RT_mean_early_diff(:,condition)                 = nanmean(RT_per_session_early_diff{condition},2);
    RT_mean_late_diff(:,condition)                  = nanmean(RT_per_session_late_diff{condition},2);    
    RT_mean_early(:,condition)                      = nanmean(RT_per_session_early{condition},2);
    RT_mean_late(:,condition)                       = nanmean(RT_per_session_late{condition},2);
end

after_stim_diff_matrix=repmat(windows'*1000+200,1,1);
R_idx=[4];
L_idx=[3];
for w=1:size(RT_raw_second_mode,1)
    RT_mean_second_mode_conditions_combined_after_stim{w}=vertcat(RT_mean_second_mode{w,:}).*1000 - after_stim_diff_matrix(w);
    RT_mean_second_mode_conditions_combined_after_stim_L{w}=vertcat(RT_mean_second_mode{w,L_idx}).*1000 - after_stim_diff_matrix(w);
    RT_mean_second_mode_conditions_combined_after_stim_R{w}=vertcat(RT_mean_second_mode{w,R_idx}).*1000 - after_stim_diff_matrix(w);
    
    RT_mean_first_mode_conditions_combined_after_stim{w}=vertcat(RT_mean_first_mode{w,:}).*1000 - vertcat(RT_mean_first_mode{2,:}).*1000;
    RT_mean_first_mode_conditions_combined_after_stim_L{w}=vertcat(RT_mean_first_mode{w,L_idx}).*1000 - vertcat(RT_mean_first_mode{2,L_idx}).*1000;
    RT_mean_first_mode_conditions_combined_after_stim_R{w}=vertcat(RT_mean_first_mode{w,R_idx}).*1000 - vertcat(RT_mean_first_mode{2,R_idx}).*1000;
    
    
%     RT_mean_second_mode_conditions_combined_after_stim_mean(w)=nanmean(vertcat(RT_mean_second_mode{w,:}).*1000 - after_stim_diff_matrix(w));
%     RT_mean_second_mode_conditions_combined_after_stim_sem(w)=sterr(vertcat(RT_mean_second_mode{w,:}).*1000 - after_stim_diff_matrix(w));
%     RT_mean_second_mode_conditions_combined_after_stim_mean_L(w)=nanmean(vertcat(RT_mean_second_mode{w,[1,3]}).*1000 - after_stim_diff_matrix(w));
%     RT_mean_second_mode_conditions_combined_after_stim_sem_L(w)=sterr(vertcat(RT_mean_second_mode{w,[1,3]}).*1000 - after_stim_diff_matrix(w));
%     RT_mean_second_mode_conditions_combined_after_stim_mean_R(w)=nanmean(vertcat(RT_mean_second_mode{w,[2,4]}).*1000 - after_stim_diff_matrix(w));
%     RT_mean_second_mode_conditions_combined_after_stim_sem_L(w)=sterr(vertcat(RT_mean_second_mode{w,[2,4]}).*1000 - after_stim_diff_matrix(w));
%     
    
    RT_raw_second_mode_conditions_combined{w}=[RT_raw_second_mode{w,:}];
    RT_raw_second_mode_conditions_combined_after_stim{w}=[RT_raw_second_mode{w,:}].*1000 - after_stim_diff_matrix(w);
    RT_raw_second_mode_conditions_combined_after_stim_L{w}=[RT_raw_second_mode{w,L_idx}].*1000 - after_stim_diff_matrix(w);
    RT_raw_second_mode_conditions_combined_after_stim_R{w}=[RT_raw_second_mode{w,R_idx}].*1000 - after_stim_diff_matrix(w);
    RT_raw_first_mode_conditions_combined{w}=[RT_raw_first_mode{w,:}].*1000;
    RT_raw_first_mode_conditions_combined_L{w}=[RT_raw_first_mode{w,L_idx}].*1000;
    RT_raw_first_mode_conditions_combined_R{w}=[RT_raw_first_mode{w,R_idx}].*1000;
    RT_raw_conditions_combined{w}=[RT_raw{w,:}];
end
raw_p_second_mode=round((cellfun(@(x)(numel(x)),RT_raw_second_mode_conditions_combined)-4)./cellfun(@(x)(numel(x)),RT_raw_conditions_combined)*100);
mean_p_second_mode=(cellfun(@(x,y) x./(x+y).*100,N_second_mode,N_first_mode,'UniformOutput',0));
if numel(L_idx)==2
mean_mean_p_second_mode_L=nanmean([mean_p_second_mode{end-3:end,1};mean_p_second_mode{end-3:end,3}],1);
mean_mean_p_second_mode_R=nanmean([mean_p_second_mode{end-3:end,2};mean_p_second_mode{end-3:end,4}],1);
sem_mean_p_second_mode_L=sterr([mean_p_second_mode{end-3:end,1};mean_p_second_mode{end-3:end,3}],1);
sem_mean_p_second_mode_R=sterr([mean_p_second_mode{end-3:end,2};mean_p_second_mode{end-3:end,4}],1);
else
mean_mean_p_second_mode_L=nanmean([mean_p_second_mode{end-3:end,3}],1);
mean_mean_p_second_mode_R=nanmean([mean_p_second_mode{end-3:end,4}],1);
sem_mean_p_second_mode_L=sterr([mean_p_second_mode{end-3:end,3}],1);
sem_mean_p_second_mode_R=sterr([mean_p_second_mode{end-3:end,4}],1);
    
    
end
raw_mean_second_mode_delay=round(nanmean([RT_raw_second_mode_conditions_combined_after_stim{end-3:end}]));
raw_mean_second_mode_delay_L=round(nanmean([RT_raw_second_mode_conditions_combined_after_stim_L{end-3:end}]));
raw_mean_second_mode_delay_R=round(nanmean([RT_raw_second_mode_conditions_combined_after_stim_R{end-3:end}]));
 
raw_sem_second_mode_delay=(sterr([RT_raw_second_mode_conditions_combined_after_stim{end-3:end}])); 
raw_sem_second_mode_delay_L=(sterr([RT_raw_second_mode_conditions_combined_after_stim_L{end-3:end}])); 
raw_sem_second_mode_delay_R=(sterr([RT_raw_second_mode_conditions_combined_after_stim_R{end-3:end}]));
 
raw_mean_first_mode_delay=round(cellfun(@(x) nanmean(x)-nanmean(RT_raw_first_mode_conditions_combined{1}),RT_raw_first_mode_conditions_combined));
raw_mean_first_mode_delay_L=round(cellfun(@(x) nanmean(x)-nanmean(RT_raw_first_mode_conditions_combined_L{1}),RT_raw_first_mode_conditions_combined_L));
raw_mean_first_mode_delay_R=round(cellfun(@(x) nanmean(x)-nanmean(RT_raw_first_mode_conditions_combined_R{1}),RT_raw_first_mode_conditions_combined_R));

raw_sem_first_mode_delay=round(cellfun(@(x) sterr(x),RT_raw_first_mode_conditions_combined).*10)./10;
raw_sem_first_mode_delay_L=round(cellfun(@(x) sterr(x),RT_raw_first_mode_conditions_combined_L).*10)./10;
raw_sem_first_mode_delay_R=round(cellfun(@(x) sterr(x),RT_raw_first_mode_conditions_combined_R).*10)./10;


raw_mean_first_mode_delay=raw_mean_first_mode_delay(end-3:end);
raw_mean_first_mode_delay_L=raw_mean_first_mode_delay_L(end-3:end);
raw_mean_first_mode_delay_R=raw_mean_first_mode_delay_R(end-3:end);
raw_sem_first_mode_delay=raw_sem_first_mode_delay(end-3:end);
raw_sem_first_mode_delay_L=raw_sem_first_mode_delay_L(end-3:end);
raw_sem_first_mode_delay_R=raw_sem_first_mode_delay_R(end-3:end);

mean_mean_second_mode_delay=round(nanmean(vertcat(RT_mean_second_mode_conditions_combined_after_stim{end-3:end})));
mean_mean_second_mode_delay_L=round(nanmean(vertcat(RT_mean_second_mode_conditions_combined_after_stim_L{end-3:end})));
mean_mean_second_mode_delay_R=round(nanmean(vertcat(RT_mean_second_mode_conditions_combined_after_stim_R{end-3:end})));    

mean_sem_second_mode_delay=round(sterr(vertcat(RT_mean_second_mode_conditions_combined_after_stim{end-3:end})).*10)./10;
mean_sem_second_mode_delay_L=round(sterr(vertcat(RT_mean_second_mode_conditions_combined_after_stim_L{end-3:end})).*10)./10;
mean_sem_second_mode_delay_R=round(sterr(vertcat(RT_mean_second_mode_conditions_combined_after_stim_R{end-3:end})).*10)./10;

mean_mean_first_mode_delay=round(nanmean([RT_mean_first_mode_conditions_combined_after_stim{end-3:end}],1));
mean_mean_first_mode_delay_L=round(nanmean([RT_mean_first_mode_conditions_combined_after_stim_L{end-3:end}],1));
mean_mean_first_mode_delay_R=round(nanmean([RT_mean_first_mode_conditions_combined_after_stim_R{end-3:end}],1));   

mean_sem_first_mode_delay=round(sterr([RT_mean_first_mode_conditions_combined_after_stim{end-3:end}]).*10)./10;
mean_sem_first_mode_delay_L=round(sterr([RT_mean_first_mode_conditions_combined_after_stim_L{end-3:end}]).*10)./10;
mean_sem_first_mode_delay_R=round(sterr([RT_mean_first_mode_conditions_combined_after_stim_R{end-3:end}]).*10)./10; 

after_stim_mean_delay_second_mode=RT_raw_mean_second_mode-repmat(after_stim_diff_matrix,1,4);


%% statistics
switch GLO.run_baseline_included_in_ttests
    case 0
        idx_to_compare_ttest_to                             = 1;
        start_idx                                           = 2;
    case 1
        idx_to_compare_ttest_to                             = 2;
        start_idx                                           = 1;
end

sta.bias_per_session_firstmode                                  = bias_per_session_first_mode(:,start_idx:end);
sta.bias_per_session_secondmode                                 = bias_per_session_second_mode(:,start_idx:end);
sta.bias_per_session_diff                                       = bias_per_session_diff(:,start_idx:end);
sta.bias_per_session_firstmode_diff                             = bias_per_session_first_mode_diff(:,start_idx:end);
sta.bias_per_session_secondmode_diff                            = bias_per_session_second_mode_diff(:,start_idx:end);
sta.bias_per_session                                            = bias_per_session(:,start_idx:end);
% sta.hitrate_IN_R                                                = Hits_IN_R(:,start_idx:end)./Total_trials_IN_R(:,start_idx:end);
% sta.hitrate_IN_L                                                = Hits_IN_L(:,start_idx:end)./Total_trials_IN_L(:,start_idx:end);
% sta.bias.mean.all                                                  = bias.mean.all(:,start_idx:end);

Conditions={'Choice_Left','Choice_Right','Instructed_Left','Instructed_Right'};
Conditions2={'L_CH','R_CH','L_IN','R_IN'};
for con=1:numel(Conditions) 
    sta.(['RT_' Conditions{con}])                       = [RT_per_session_all{start_idx:end,con}];
    sta.(['AccuracyEu_' Conditions{con}])               = [Accuracy_eu_per_session_all{start_idx:end,con}];
    sta.(['RT_difference_' Conditions{con}])            = [RT_per_session_diff{start_idx:end,con}];
    sta.(['RT_firstmode_' Conditions{con}])             = [RT_mean_first_mode{start_idx:end,con}];
    sta.(['RT_secondmode_' Conditions{con}])            = [RT_mean_second_mode{start_idx:end,con}];
    sta.(['RT_firstmode_difference_' Conditions{con}])   = [RT_mean_first_mode_diff{start_idx:end,con}];
    sta.(['RT_secondmode_difference_' Conditions{con}])  = [RT_mean_second_mode_diff{start_idx:end,con}];
    sta.(['Hit_rate_' Conditions{con}])                 = Hit_rate.per_session.(Conditions2{con})(:,start_idx:end);    
    sta.(['Vel_' Conditions{con}])                       = [Vel_per_session_all{start_idx:end,con}];
    sta.(['Vel_firstmode_' Conditions{con}])             = [Vel_mean_first_mode{start_idx:end,con}];
    sta.(['Vel_secondmode_' Conditions{con}])            = [Vel_mean_second_mode{start_idx:end,con}];
    
    sta_cell.(['RT_' Conditions{con}])                  = RT_raw(start_idx:end,con);
    sta_cell.(['RT_firstmode_' Conditions{con}])        = RT_raw_first_mode(start_idx:end,con);
    sta_cell.(['RT_secondmode_' Conditions{con}])       = RT_raw_second_mode(start_idx:end,con);
    
end
for pos=1:n_positions
        sta.(['AccuracyRadIN_pos_' num2str(pos)])                 = vertcat(Accuracy_IN_per_target_rad{pos,start_idx:end}).';  
        sta.(['AccuracyRadCH_pos_' num2str(pos)])                 = vertcat(Accuracy_CH_per_target_rad{pos,start_idx:end}).';  
        sta.(['PrecisionRadIN_pos_' num2str(pos)])                 = vertcat(Precision_IN_per_target_rad{pos,start_idx:end}).'; 
        sta.(['PrecisionRadCH_pos_' num2str(pos)])                 = vertcat(Precision_CH_per_target_rad{pos,start_idx:end}).'; 
        
        
        sta_cell.(['AccuracyRadIN_pos_' num2str(pos)])            = vertcat(Accuracy_IN_per_target_rad_raw(pos,start_idx:end)).';  
        sta_cell.(['AccuracyRadCH_pos_' num2str(pos)])            = vertcat(Accuracy_CH_per_target_rad_raw(pos,start_idx:end)).';  
end
    
condition_difference_labels                                           = {'R_IN_L_IN','R_CH_L_CH','R_CH_R_IN','L_CH_L_IN'};
for con=1:numel(condition_difference_labels)
    sta.(['RT_firstmode_' condition_difference_labels{con}])             = [RT_per_session_condition_differences_first_mode{start_idx:end,con}];
    sta.(['RT_secondmode_' condition_difference_labels{con}])            = [RT_per_session_condition_differences_second_mode{start_idx:end,con}];
    sta.(['RT_bothmodes_' condition_difference_labels{con}])             = [RT_per_session_condition_differences_all{start_idx:end,con}];
end

sta_fieldnames                                                  = fieldnames(sta);
sta_cell_fieldnames                                             = fieldnames(sta_cell);


% statistics on means
for SF=1:numel(sta_fieldnames)
    to_test                                                     = sta.(sta_fieldnames{SF});
    if ~isreal(to_test)
        to_test=real(to_test);
    end
    n_repeated_measurements                                     = size(to_test,1);
    n_windows_subjects                                          = size(to_test,2);
    subject                                                     = repmat(1:n_windows_subjects,[n_repeated_measurements 1]);
    subject                                                     = subject(:);
    measurement                                                 = repmat(1:n_repeated_measurements,[1 n_windows_subjects]);
    measurement                                                 = measurement(:);
    inputs_for_rm_anova                                         = [to_test(:), measurement, subject];
    
    %     % using paired ttest (non paired is ttest2)
    if strcmp(GLO.stat_to_use,'signed rank')
       if any(strfind(sta_fieldnames{SF},'mode')) || any(strfind(sta_fieldnames{SF},'pos')) ||  any(strfind(sta_fieldnames{SF},'Choice'))
            P1  = kruskalwallis(to_test,[],'off');
            stat_to_use                                           = 'ranksum';
        else
            if any(any(isnan(to_test))) || size(to_test,1)==1 || size(to_test,2)==1
               a=1; 
            end
%             try
%             P1  = friedman(to_test,1,'off');
%             catch eee
%             a=1;     
%             end
            P1  = kruskalwallis(to_test,[],'off');
            %stat_to_use                                           = GLO.stat_to_use;
            stat_to_use                                           = 'ranksum';
        end
    else
        if any(strfind(sta_fieldnames{SF},'mode')) || any(strfind(sta_fieldnames{SF},'pos'))
            P1                                                    = anova1(to_test,[],'off');
            %stat_to_use                                           = GLO.stat_to_use; %%%?
            %stat_to_use                                           = 'ttest paired'; %%%?
            stat_to_use                                           = 'ttest unpaired';
        else
            [F1,P1] = RMAOV1(inputs_for_rm_anova,0.05);
            stat_to_use                                           = GLO.stat_to_use;
        end
    end
    
    if P1<0.05
        pairs=[repmat(idx_to_compare_ttest_to,n_windows_subjects,1) [1:n_windows_subjects]'];
        pairs(idx_to_compare_ttest_to,:)=[];
        [h_bonf.(sta_fieldnames{SF}),p_bonf.(sta_fieldnames{SF}),sigPairs_bonf.(sta_fieldnames{SF})] = ttest_bonf(to_test,pairs);
        if idx_to_compare_ttest_to==1
            p_bonf.(sta_fieldnames{SF})=[NaN NaN p_bonf.(sta_fieldnames{SF})];
        elseif idx_to_compare_ttest_to==2
            p_bonf.(sta_fieldnames{SF})=[p_bonf.(sta_fieldnames{SF})(1) NaN p_bonf.(sta_fieldnames{SF})(2:end)];
        end
        for i=1:n_windows_subjects
            [~, ps_paired.(sta_fieldnames{SF})(i)]      = ttest(to_test(:,idx_to_compare_ttest_to),to_test(:,i));
            [~, ps_unpaired.(sta_fieldnames{SF})(i)]    = ttest2(to_test(:,idx_to_compare_ttest_to),to_test(:,i));
            if all(isnan(to_test(:,idx_to_compare_ttest_to)) | isnan(to_test(:,i)))
                p_signrank.(sta_fieldnames{SF})(i)     = NaN;
                p_ranksum.(sta_fieldnames{SF})(i)      = NaN;
            else
                p_signrank.(sta_fieldnames{SF})(i)     = signrank(to_test(:,idx_to_compare_ttest_to),to_test(:,i));
                p_ranksum.(sta_fieldnames{SF})(i)      = ranksum(to_test(:,idx_to_compare_ttest_to),to_test(:,i));
            end
            %[rank_ps.(sta_fieldnames{SF})(i) rank_Hs.(sta_fieldnames{SF})(i),  rank_ooo.(sta_fieldnames{SF})(i)]    = ranksum(to_test(:,idx_to_compare_ttest_to),to_test(:,i));
        end
        %[p_ano,anovatab,stats]              =anova1(to_test,[],'off');
        %pdunn.(sta_fieldnames{SF})          = dunnett(stats, [2:length(stats.means)], 1);
        %% Multiple Comparsion
         ps_paired.(sta_fieldnames{SF})     =ps_paired.(sta_fieldnames{SF})*(n_windows_subjects-1); 
        ps_unpaired.(sta_fieldnames{SF})    =ps_unpaired.(sta_fieldnames{SF})*(n_windows_subjects-1);
        p_bonf.(sta_fieldnames{SF})         =p_bonf.(sta_fieldnames{SF})*(n_windows_subjects-1);
         p_signrank.(sta_fieldnames{SF})     =p_signrank.(sta_fieldnames{SF})*(n_windows_subjects-1);
         p_ranksum.(sta_fieldnames{SF})      =p_ranksum.(sta_fieldnames{SF})*(n_windows_subjects-1);
    else
        ps_paired.(sta_fieldnames{SF})      =NaN(1,n_windows_subjects);
        ps_unpaired.(sta_fieldnames{SF})    =NaN(1,n_windows_subjects);
        p_bonf.(sta_fieldnames{SF})         =NaN(1,n_windows);
        p_signrank.(sta_fieldnames{SF})     =NaN(1,n_windows);
        p_ranksum.(sta_fieldnames{SF})      =NaN(1,n_windows);
    end
    if idx_to_compare_ttest_to ==1
        ps_paired.(sta_fieldnames{SF})      =[NaN  ps_paired.(sta_fieldnames{SF})];
        ps_unpaired.(sta_fieldnames{SF})    =[NaN  ps_unpaired.(sta_fieldnames{SF})];
        p_signrank.(sta_fieldnames{SF})     =[NaN  p_signrank.(sta_fieldnames{SF})];
        p_ranksum.(sta_fieldnames{SF})      =[NaN  p_ranksum.(sta_fieldnames{SF})];
    end
    switch stat_to_use
        case 'signed rank'
            asterisk.(sta_fieldnames{SF})=cell(size( p_signrank.(sta_fieldnames{SF})));
            asterisk.(sta_fieldnames{SF})(p_signrank.(sta_fieldnames{SF})>=0.05 | isnan(p_signrank.(sta_fieldnames{SF})))=deal({''});
            asterisk.(sta_fieldnames{SF})(p_signrank.(sta_fieldnames{SF})<0.05 & p_signrank.(sta_fieldnames{SF})>=0.01)=deal({'*'});
            asterisk.(sta_fieldnames{SF})(p_signrank.(sta_fieldnames{SF})<0.01)=deal({'**'});
            
        case 'ranksum'
            asterisk.(sta_fieldnames{SF})=cell(size( p_ranksum.(sta_fieldnames{SF})));
            asterisk.(sta_fieldnames{SF})(p_ranksum.(sta_fieldnames{SF})>=0.05 | isnan(p_ranksum.(sta_fieldnames{SF})))=deal({''});
            asterisk.(sta_fieldnames{SF})(p_ranksum.(sta_fieldnames{SF})<0.05 & p_ranksum.(sta_fieldnames{SF})>=0.01)=deal({'*'});
            asterisk.(sta_fieldnames{SF})(p_ranksum.(sta_fieldnames{SF})<0.01)=deal({'**'});
             
        case 'ttest paired'
            asterisk.(sta_fieldnames{SF})=cell(size( ps_paired.(sta_fieldnames{SF})));
            asterisk.(sta_fieldnames{SF})(ps_paired.(sta_fieldnames{SF})>=0.05 | isnan(ps_paired.(sta_fieldnames{SF})))=deal({''});
            asterisk.(sta_fieldnames{SF})(ps_paired.(sta_fieldnames{SF})<0.05 & ps_paired.(sta_fieldnames{SF})>=0.01)=deal({'*'});
            asterisk.(sta_fieldnames{SF})(ps_paired.(sta_fieldnames{SF})<0.01)=deal({'**'});
            
        case 'ttest unpaired'
            asterisk.(sta_fieldnames{SF})=cell(size( ps_unpaired.(sta_fieldnames{SF})));
            asterisk.(sta_fieldnames{SF})(ps_unpaired.(sta_fieldnames{SF})>=0.05 | isnan(ps_unpaired.(sta_fieldnames{SF})))=deal({''});
            asterisk.(sta_fieldnames{SF})(ps_unpaired.(sta_fieldnames{SF})<0.05 & ps_unpaired.(sta_fieldnames{SF})>=0.01)=deal({'*'});
            asterisk.(sta_fieldnames{SF})(ps_unpaired.(sta_fieldnames{SF})<0.01)=deal({'**'});
            
        case 'ttest_bonf'
            asterisk.(sta_fieldnames{SF})=cell(size( p_bonf.(sta_fieldnames{SF})));
            asterisk.(sta_fieldnames{SF})(p_bonf.(sta_fieldnames{SF})>=0.05 | isnan(p_bonf.(sta_fieldnames{SF})))=deal({''});
            asterisk.(sta_fieldnames{SF})(p_bonf.(sta_fieldnames{SF})<0.05 & p_bonf.(sta_fieldnames{SF})>=0.01)=deal({'*'});
            asterisk.(sta_fieldnames{SF})(p_bonf.(sta_fieldnames{SF})<0.01)=deal({'**'});
    end
%     if strfind(sta_fieldnames{SF},'mode')
%         asterisk.(sta_fieldnames{SF})=cell(size( ps_paired.(sta_fieldnames{SF})));
%         asterisk.(sta_fieldnames{SF})(ps_paired.(sta_fieldnames{SF})>=0.05 | isnan(ps_paired.(sta_fieldnames{SF})))=deal({''});
%         asterisk.(sta_fieldnames{SF})(ps_paired.(sta_fieldnames{SF})<0.05 & ps_paired.(sta_fieldnames{SF})>=0.01)=deal({'*'});
%         asterisk.(sta_fieldnames{SF})(ps_paired.(sta_fieldnames{SF})<0.01)=deal({'**'});
%     end
end

%statistics on raw data
for SF=1:numel(sta_cell_fieldnames)
    to_test                                                     = sta_cell.(sta_cell_fieldnames{SF});
    
    for i=1:numel(to_test)
        max_to_test(:,i)=numel(to_test{i,:}');
    end
    
    to_test_raw=NaN(numel(max_to_test),max(max_to_test));
    
    for i=1:numel(to_test)
        if ~isreal(to_test{i})
            to_test{i}=real(to_test{i});
        end
        to_test_raw(i,1:numel(to_test{i}))=to_test{i};        
    end
    to_test_raw=to_test_raw';
    if all(isempty(to_test_raw))
       to_test_raw=NaN(1,size(to_test_raw,2));
    end
    p_anova_raw_rt                                                    = anova1(to_test_raw,[],'off');
    
    % non paired ttest for concatinated raw data    
    if p_anova_raw_rt<0.05
        for w=1:n_windows_subjects
            [~, ps_raw_unpaired.(sta_cell_fieldnames{SF})(w)]    = ttest2(to_test_raw(:,idx_to_compare_ttest_to),to_test_raw(:,w));  
            hominput=[to_test_raw(:,idx_to_compare_ttest_to),ones(numel(to_test_raw(:,idx_to_compare_ttest_to)),1);...
            to_test_raw(:,w),ones(numel(to_test_raw(:,idx_to_compare_ttest_to)),1)*2];
            [ps_raw_levene.(sta_cell_fieldnames{SF})(w)]    = NaN; % Levenetest(hominput(~isnan(hominput(:,1)),:)); %
        end
        ps_raw_unpaired.(sta_cell_fieldnames{SF})   =ps_raw_unpaired.(sta_cell_fieldnames{SF})*(n_windows_subjects-1);
        ps_raw_levene.(sta_cell_fieldnames{SF})   =ps_raw_levene.(sta_cell_fieldnames{SF})*(n_windows_subjects-1);
    else
        ps_raw_unpaired.(sta_cell_fieldnames{SF})=NaN(1,n_windows_subjects);  
        ps_raw_levene.(sta_cell_fieldnames{SF})=NaN(1,n_windows_subjects);       
    end
    if idx_to_compare_ttest_to ==1
        ps_raw_unpaired.(sta_cell_fieldnames{SF})=[NaN  ps_raw_unpaired.(sta_cell_fieldnames{SF})];
        ps_raw_levene.(sta_cell_fieldnames{SF})=[NaN  ps_raw_levene.(sta_cell_fieldnames{SF})];
    end
    asterisks.(sta_cell_fieldnames{SF})=cell(size( ps_raw_unpaired.(sta_cell_fieldnames{SF})));
    asterisks.(sta_cell_fieldnames{SF})(ps_raw_unpaired.(sta_cell_fieldnames{SF})>=0.05 | isnan(ps_raw_unpaired.(sta_cell_fieldnames{SF})))=deal({''});
    asterisks.(sta_cell_fieldnames{SF})(ps_raw_unpaired.(sta_cell_fieldnames{SF})<0.05 & ps_raw_unpaired.(sta_cell_fieldnames{SF})>=0.01)=deal({'*'});
    asterisks.(sta_cell_fieldnames{SF})(ps_raw_unpaired.(sta_cell_fieldnames{SF})<0.01)=deal({'**'});
    
    asterisk_pluses.(sta_cell_fieldnames{SF})=cell(size( ps_raw_levene.(sta_cell_fieldnames{SF})));
    asterisk_pluses.(sta_cell_fieldnames{SF})(ps_raw_levene.(sta_cell_fieldnames{SF})>=0.05 | isnan(ps_raw_levene.(sta_cell_fieldnames{SF})))=deal({''});
    asterisk_pluses.(sta_cell_fieldnames{SF})(ps_raw_levene.(sta_cell_fieldnames{SF})<0.05 & ps_raw_levene.(sta_cell_fieldnames{SF})>=0.01)=deal({'+'});
    asterisk_pluses.(sta_cell_fieldnames{SF})(ps_raw_levene.(sta_cell_fieldnames{SF})<0.01)=deal({'++'});
end

condition_difference_labels                                           = {'R_IN_L_IN','R_CH_L_CH','R_CH_R_IN','L_CH_L_IN'};
condition_difference_indexes                                          = [4,3;2,1;2,4;1,3];
Conditions={'Choice_Left','Choice_Right','Instructed_Left','Instructed_Right'};
Modes={'RT_','RT_firstmode_','RT_secondmode_'};
%statistics on difrfeences in raw data
for d=1:numel(condition_difference_labels)
    for m=1:numel(Modes)
        to_test1                                                     = sta_cell.([Modes{m} Conditions{condition_difference_indexes(d,1)}]);
        to_test2                                                     = sta_cell.([Modes{m} Conditions{condition_difference_indexes(d,2)}]);
        
        for i=1:numel(to_test1)
            max_to_test(:,i)=max(numel(to_test1{i,:}'),numel(to_test1{i,:}'));
        end
        
        to_test_rt_raw1=NaN(numel(max_to_test),max(max_to_test));
        to_test_rt_raw2=NaN(numel(max_to_test),max(max_to_test));
        
        for i=1:numel(to_test1)
            to_test_rt_raw1(i,1:numel(to_test1{i}))=to_test1{i};
            to_test_rt_raw2(i,1:numel(to_test2{i}))=to_test2{i};
        end
        to_test_rt_raw1=to_test_rt_raw1';
        to_test_rt_raw2=to_test_rt_raw2';
        %p_anova_raw_rt                                                    = anova1(to_test_rt_raw,[],'off');
        p_anova_raw_rt                                                    = 0;
        
        % non paired ttest for concatinated raw data
        if p_anova_raw_rt<0.05
            
            
            for i=1+1:n_windows_subjects
                %[~, ps_raw_unpaired.([Modes{m} condition_difference_labels{d}])(i)]    = ttest2(to_test_rt_raw1(:,i) - nanmean(to_test_rt_raw1(:,idx_to_compare_ttest_to)),to_test_rt_raw2(:,i)- nanmean(to_test_rt_raw2(:,idx_to_compare_ttest_to)));
            ps_raw_unpaired.([Modes{m} condition_difference_labels{d}])(i)    = 1;
            end
            [~, ps_raw_unpaired.([Modes{m} condition_difference_labels{d}])(idx_to_compare_ttest_to)]    = ttest2(to_test_rt_raw1(:,idx_to_compare_ttest_to),to_test_rt_raw2(:,idx_to_compare_ttest_to));
            ps_raw_unpaired.([Modes{m} condition_difference_labels{d}])   =ps_raw_unpaired.([Modes{m} condition_difference_labels{d}])*(n_windows_subjects-1);
        else
            ps_raw_unpaired.([Modes{m} condition_difference_labels{d}])=NaN(1,n_windows_subjects);
        end
        if idx_to_compare_ttest_to ==1
            ps_raw_unpaired.([Modes{m} condition_difference_labels{d}])=[NaN  ps_raw_unpaired.([Modes{m} condition_difference_labels{d}])];
        end
        asterisks.([Modes{m} condition_difference_labels{d}])=cell(size( ps_raw_unpaired.([Modes{m} condition_difference_labels{d}])));
        asterisks.([Modes{m} condition_difference_labels{d}])(ps_raw_unpaired.([Modes{m} condition_difference_labels{d}])>=0.05 | isnan(ps_raw_unpaired.([Modes{m} condition_difference_labels{d}])))=deal({''});
        asterisks.([Modes{m} condition_difference_labels{d}])(ps_raw_unpaired.([Modes{m} condition_difference_labels{d}])<0.05 & ps_raw_unpaired.([Modes{m} condition_difference_labels{d}])>=0.01)=deal({'*'});
        asterisks.([Modes{m} condition_difference_labels{d}])(ps_raw_unpaired.([Modes{m} condition_difference_labels{d}])<0.01)=deal({'**'});
    end
end

% rt_all_indexes=[idx_to_compare_ttest_to:size(RT_per_session,1)];
% % correlations
% for rows=1:size(bias_per_session,2)
%     for columns=1:size(RT_per_session,2)
%         [r p]=corrcoef(RT_per_session_all{rt_all_indexes(rows),columns},bias_per_session(:,rows));
%         r_rt_bias(rows,columns)=r(1,2);
%         p_rt_bias(rows,columns)=p(1,2);
%     end
% end


if GLO.run_baseline_included_in_ttests
    start_label                                                 =1;
    start_label_significance_plots                              =1;
else
    start_label                                                 =2;
    start_label_significance_plots                              =3;
end

 
    
if GLO.table
    table_batch.title_ui_table       = [monkey ' Bias and RT for each window per session'];
    table_batch.labels_RT            = {'Choice left','Choice right','Instructed left','Instructed right'};
    
    table_batch.labels              = GLO.Labels(start_label:end);
    table_batch.labels_bias         = strcat(table_batch.labels,'_dBS');
    table_batch.labels_RT_L         = strcat(table_batch.labels,'_dRT_L');
    table_batch.labels_RT_R         = strcat(table_batch.labels,'_dRT_R');
    table_batch.labels_RT_L_CH      = strcat(table_batch.labels,'_dRT_L_CH');
    table_batch.labels_RT_R_CH      = strcat(table_batch.labels,'_dRT_R_CH');
    table_batch.labels_sig_bias     = strcat(table_batch.labels,'_*BS');
    table_batch.labels_sig_RT_L     = strcat(table_batch.labels,'_*RT_L');
    table_batch.labels_sig_RT_R     = strcat(table_batch.labels,'_*RT_R');    
    
    table_batch.labels_1st          = strcat(table_batch.labels,'_1st');
    table_batch.labels_1st_sem      = strcat(table_batch.labels,'_1st_sem');
    
    table_batch.labels_p2nd_L       = strcat(table_batch.labels,'_p2nd_L');
    table_batch.labels_p2nd_R       = strcat(table_batch.labels,'_p2nd_R');
    table_batch.labels_p2nd_L_sem   = strcat(table_batch.labels,'_p2nd_L_sem');
    table_batch.labels_p2nd_R_sem   = strcat(table_batch.labels,'_p2nd_R_sem');
    
    table_batch.labels_2nd_L        = strcat('2nd_L');
    table_batch.labels_2nd_R        = strcat('2nd_R');
    table_batch.labels_2nd_L_sem    = strcat('2nd_L_sem');
    table_batch.labels_2nd_R_sem    = strcat('2nd_R_sem');    
        
    table_batch.labels_MS_L         = strcat(table_batch.labels,'_%multi_L');
    table_batch.labels_MS_R         = strcat(table_batch.labels,'_%multi_R');
    table_batch.labels_SM_L         = strcat(table_batch.labels,'_%2ndmode_L');
    table_batch.labels_SM_R         = strcat(table_batch.labels,'_%2ndmode_R');
    
    for FN={'labels_1st','labels_1st_sem','labels_p2nd_L','labels_p2nd_R','labels_p2nd_L_sem','labels_p2nd_R_sem'}
        table_batch.(FN{:})=table_batch.(FN{:})(end-3:end);
    end
    
    table_batch.titles              = [table_batch.labels_bias table_batch.labels_RT_L table_batch.labels_RT_R table_batch.labels_RT_L_CH table_batch.labels_RT_R_CH...
                                        table_batch.labels_sig_bias table_batch.labels_sig_RT_L table_batch.labels_sig_RT_R ...
                                        table_batch.labels_MS_L table_batch.labels_MS_R table_batch.labels_SM_L  table_batch.labels_SM_R ...
                                        table_batch.labels_1st table_batch.labels_1st_sem table_batch.labels_p2nd_L table_batch.labels_p2nd_R table_batch.labels_p2nd_L_sem table_batch.labels_p2nd_R_sem ...
                                        table_batch.labels_2nd_L table_batch.labels_2nd_R table_batch.labels_2nd_L_sem table_batch.labels_2nd_R_sem];
    
    
            Fraction_multistep{w,c}(s,:)            =sum(batch_out.N{s}.multistep{w,c}) /numel(batch_out.N{s}.multistep{w,c});
            %Fraction_second_mode{w,c}(s,:)          =sum(idx_second_mode{w,c})/sum(idx_first_mode{w,c});
    
            Fraction_second_mode{w,c}(s,:)          =sum(idx_second_mode{w,c})/(sum(idx_first_mode{w,c}) + sum(idx_second_mode{w,c}));
            
    table_batch.MS_L                = round(cell2mat(Fraction_multistep(start_label:end,3)')*100);  %% 3== IN_L
    table_batch.MS_R                = round(cell2mat(Fraction_multistep(start_label:end,4)')*100);  %% 4== IN_R
    table_batch.SM_L                = round(cell2mat(Fraction_second_mode(start_label:end,3)')*100);  %% 3== IN_L
    table_batch.SM_R                = round(cell2mat(Fraction_second_mode(start_label:end,4)')*100);  %% 4== IN_R
    table_batch.bias                = round(bias_per_session_diff(:,start_label:end));
    table_batch.RT_L                = round(cell2mat(RT_per_session_diff(start_label:end,3)')*1000);  %% 3== IN_L
    table_batch.RT_R                = round(cell2mat(RT_per_session_diff(start_label:end,4)')*1000);  %% 4== IN_R   
    table_batch.RT_L_CH             = round(cell2mat(RT_per_session_diff(start_label:end,1)')*1000);  %% 1== CH_L
    table_batch.RT_R_CH             = round(cell2mat(RT_per_session_diff(start_label:end,2)')*1000);  %% 2== CH_R   
    table_batch.sig_bias            = batch_out.n_sig_bias.significant_per_session(:,start_label:end);
    table_batch.sig_RT_L            = batch_out.n_sig_RTs.significant_per_session.L_IN(:,start_label:end);  %% 3== IN_L
    table_batch.sig_RT_R            = batch_out.n_sig_RTs.significant_per_session.R_IN(:,start_label:end);  %% 4== IN_R  
    
    table_batch.first               = repmat(mean_mean_first_mode_delay,n_sessions,1);
    table_batch.first_sem           = repmat(mean_sem_first_mode_delay,n_sessions,1);
    table_batch.p2nd_L              = repmat(mean_mean_p_second_mode_L,n_sessions,1);
    table_batch.p2nd_R              = repmat(mean_mean_p_second_mode_R,n_sessions,1);
    table_batch.p2nd_L_sem          = repmat(sem_mean_p_second_mode_L,n_sessions,1);
    table_batch.p2nd_R_sem          = repmat(sem_mean_p_second_mode_R,n_sessions,1);    
    table_batch.second_L            = repmat(mean_mean_second_mode_delay_L,n_sessions,1);
    table_batch.second_R            = repmat(mean_mean_second_mode_delay_R,n_sessions,1);
    table_batch.second_L_sem        = repmat(mean_sem_second_mode_delay_L,n_sessions,1);
    table_batch.second_R_sem        = repmat(mean_sem_second_mode_delay_R,n_sessions,1);   
    
    table_batch.data                = [table_batch.bias table_batch.RT_L table_batch.RT_R table_batch.RT_L_CH table_batch.RT_R_CH...
                                        table_batch.sig_bias table_batch.sig_RT_L table_batch.sig_RT_R...
                                        table_batch.MS_L table_batch.MS_R table_batch.SM_L table_batch.SM_R ...
                                        table_batch.first table_batch.first_sem table_batch.p2nd_L table_batch.p2nd_R table_batch.p2nd_L_sem table_batch.p2nd_R_sem ...         = sem_mean_p_second_mode_R;    
                                        table_batch.second_L table_batch.second_R table_batch.second_L_sem table_batch.second_R_sem];
    
    GLO.table_per_batch.data=table_batch.data;
    GLO.table_per_batch.titles=table_batch.titles;
    
    cell2uitable(table_batch.data,table_batch.titles, table_batch.title_ui_table);    
end

%% PLOTING



% PLOTTING PARAMETERS
col                                                         = jet(numel(GLO.Labels)); 

col(1,:)                                                    = [0.2 0.2 0.2];
col(2,:)                                                    = [0.5 0.5 0.5];
alfa                                                        = 0.5;
if GLO.type_to_use                 == 3
 % elegant memory saccade colors
col                                                    = [0.2 0.2 0.2;0.5 0.5 0.5;0.859 0.275 0.6; 0.65 0.247 0.6; 0.35 0.3412 0.6471; 0.122 0.255 0.604];   
end
    


GLO.fontsize_titles_big             = 24;
GLO.fontsize_labels_big             = 20;
GLO.fontsize_ticks_big              = 18;
GLO.fontsize_legends_big             = 18;
GLO.fontsize_text_big               = 18;

GLO.fontsize_titles_small           = 20;
GLO.fontsize_labels_small           = 16;
GLO.fontsize_ticks_small            = 14;
GLO.fontsize_legends_small           = 14;
GLO.fontsize_text_small             = 14;

GLO.fontsize_titles_extra_small     = 16;
GLO.fontsize_labels_extra_small     = 12;
GLO.fontsize_ticks_extra_small      = 10;
GLO.fontsize_legends_extra_small     = 10;
GLO.fontsize_text_extra_small       = 10;

GLO.fontsize_temp                   = 14;

GLO.linewidth                       = 1;
GLO.col                             = col;
GLO.markersize                      = 3;
GLO.alfa                            = alfa;
GLO.n_windows                       = n_windows;


% PLOT limits
switch GLO.effector_to_use
    
    case 0
        y_lim_RT                                                    = [100 450];
        y_lim_RT_firstmode                                          = [100 250];
        y_lim_RT_secondmode                                         = [100 450];
        y_tick_RT                                                   = [100:50:450];
        y_lim_RT_errorbars                                          = [125 325];
        y_lim_RT_scatterplot                                        = [80 400];
        y_tick_RT_errorbars                                         = [100:50:350];
        y_lim_RT_diff_errorbars                                     = [-30 150];
        y_tick_RT_diff_errorbars                                    = [-50:50:150];
        y_lim_RT_diff                                               = [-80 230];
        y_lim_RT_con_diff                                           = [-100 200];
        y_tick_RT_diff                                              = [-60:50:250];
        y_lim_RT_indexes                                            = [125 350];
        y_tick_RT_indexes                                           = [100:50:350];
        
        y_lim_RT_comp_errorbars                                     = [100 400];
        y_tick_RT_comp_errorbars                                    = [100:50:400];
        
    case 4
        y_lim_RT                                                    = [100 400];
        y_lim_RT_firstmode                                          = [100 400];
        y_lim_RT_secondmode                                         = [100 400];
        y_tick_RT                                                   = [100:50:400];
        y_lim_RT_errorbars                                          = [100 450];
        y_tick_RT_errorbars                                         = [100:50:450];
        y_lim_RT_scatterplot                                        = [200 600];
        y_lim_RT_diff_errorbars                                     = [-30 13];
        y_tick_RT_diff_errorbars                                    = [-50:50:150];
        y_lim_RT_diff                                               = [-60 150];
        y_lim_RT_con_diff                                           = [-100 200];
        y_tick_RT_diff                                              = [-60:50:150];
        y_lim_RT_indexes                                            = [150 450];
        y_tick_RT_indexes                                           = [50:50:550];
        y_lim_RT_comp_errorbars                                     = [100 400];
        y_tick_RT_comp_errorbars                                    = [100:50:400];
end

y_lim_bias                                                  = [-10 100];
y_lim_bias_diff                                             = [-70 90];
% y_lim_bias_errorbars                                        = [0 70];
y_lim_bias_errorbars                                        = [0 100];
% y_lim_bias_diff_errorbars                                   = [-30 50];
y_lim_bias_diff_errorbars                                   = [-50 100];
y_lim_bias_indexes                                          = [-10 90];
y_lim_bias_diff_indexes                                     = [-40 50];
y_lim_N_sessions                                            = [-15 15];


print_out=horzcat(' ', monkey, ', ', batch_title, ', ', num2str(n_sessions), ' sites, ', GLO.excentricity, ' targets', ', ', 'saccade defined as the ', GLO.saccade_type);


if any(ismember(GLO.plot_to_show,[1,-1]))
        % FIGURE 1, PLOT 1: BIAS (Spatial preference modulation)
    plot_1_title='Summary 1- Bias_RT';
    summary_1                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_1_title);
    subplot(4,4,[09, 13])
    hold on
    title('Spatial preference modulation','fontsize',GLO.fontsize_titles_small);
    % shadedErrorBar(1:numel(bias.mean.all),bias.mean.all,bias.sem.all,{'r-o','markerfacecolor','r'})
    for idx_error_bar                                           = start_label:numel(GLO.Labels)
        e_bias(idx_error_bar)                                   = errorbar(idx_error_bar,bias.mean.all(idx_error_bar),bias.sem.all(idx_error_bar));
        if GLO.ttest_text
            text(idx_error_bar-0.2,bias.mean.all(idx_error_bar)+bias.sem.all(idx_error_bar)+diff(y_lim_bias_errorbars)/50,asterisk.bias_per_session(idx_error_bar));
        end
        set(e_bias(idx_error_bar),'Marker','s','MarkerFaceColor',col(idx_error_bar,:),'MarkerEdgeColor',col(idx_error_bar,:),'Color',col(idx_error_bar,:),'Linewidth',GLO.linewidth,'MarkerSize',GLO.markersize);
        
        h_pre_children                                          = get(e_bias(idx_error_bar),'children');
        Xdata                                                   = get(h_pre_children(2),'Xdata');
        tempo                                                   = 4:3:length(Xdata);
        tempo(3:3:end) = [];
        xleft                                                   = tempo; xright = tempo+1;
        xmean                                                   = (xleft+xright)/2;
        Xdata(xleft)                                            = [idx_error_bar-0.3 idx_error_bar-0.3];
        Xdata(xright)                                           = [idx_error_bar+0.3 idx_error_bar+0.3];
        set(h_pre_children(2),'Xdata',Xdata)
    end
    axis('square')
    set(gca,'ylim',y_lim_bias_errorbars,'xlim',[start_label-1 numel(GLO.Labels)+1],'Xtick',start_label:numel(GLO.Labels),'xticklabel',GLO.Labels(start_label:end),'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
    xlabel('Stimulation onset to target presentation [ms]','fontsize',GLO.fontsize_labels_small)
    ylabel('Contraversive target selection [%]','fontsize',GLO.fontsize_labels_small)
    box on
    
    
    % FIGURE 1, PLOT 2: BIAS DIFFERENCE TO TRIAL BASELINE (Spatial preference difference to trial baseline)
    subplot(4,4,[10, 14])
    hold on
    title({'Spatial preference '; 'difference to trial baseline'},'fontsize',GLO.fontsize_titles_small);
    % shadedErrorBar(1:numel(bias.mean.all),bias.mean.all,bias.sem.all,{'r-o','markerfacecolor','r'})
    for idx_error_bar                                           = start_label:numel(GLO.Labels)
        e_bias2(idx_error_bar)                                  = errorbar(idx_error_bar,bias.mean.all_diff(idx_error_bar),bias.sem.all_diff(idx_error_bar));
        if GLO.ttest_text
            text(idx_error_bar-0.2,bias.mean.all_diff(idx_error_bar)+bias.sem.all_diff(idx_error_bar)+diff(y_lim_bias_diff_errorbars)/50,asterisk.bias_per_session_diff(idx_error_bar));
        end
        set(e_bias2(idx_error_bar),'Marker','s','MarkerFaceColor',col(idx_error_bar,:),'MarkerEdgeColor',col(idx_error_bar,:),'Color',col(idx_error_bar,:),'Linewidth',GLO.linewidth,'MarkerSize',GLO.markersize);
        
        h_pre_children                                          = get(e_bias2(idx_error_bar),'children');
        Xdata                                                   = get(h_pre_children(2),'Xdata');
        tempo                                                   = 4:3:length(Xdata);
        tempo(3:3:end)                                          = [];
        xleft                                                   = tempo;                                    xright = tempo+1;
        xmean                                                   = (xleft+xright)/2;
        Xdata(xleft)                                            = [idx_error_bar-0.3 idx_error_bar-0.3];
        Xdata(xright)                                           = [idx_error_bar+0.3 idx_error_bar+0.3];
        set(h_pre_children(2),'Xdata',Xdata)
    end
    axis('square')
    set(gca,'ylim',y_lim_bias_diff_errorbars,'xlim',[start_label-1 numel(GLO.Labels)+1],'Xtick',start_label:numel(GLO.Labels),'xticklabel',GLO.Labels(start_label:end),'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
    xlabel('Stimulation onset to target presentation [ms]','fontsize',GLO.fontsize_labels_small)
    ylabel('Contraversive target selection difference [%]','fontsize',GLO.fontsize_labels_small)
    line([0,numel(GLO.Labels)+1],[0,0],'Color',[.8 .8 .8],'LineStyle',':')
    box on
    
    
    % FIGURE 1, PLOT 3: REACTION TIMES (4 PLOTS, 'Choice left','Choice right','Instructed left','Instructed right')
    Labels_RT_comparison                                        = {'Choice contra','Choice ipsi','Instructed contra','Instructed ipsi'};
    subplot_assignment                                          = [3 4 1 2];
    % plot_matrix={idx_L_CH_S, idx_R_CH_S,  idx_L_IN_S, idx_R_IN_S; idx_L_CH_B,  idx_R_CH_B, idx_L_IN_B, idx_R_IN_B};
    for k                                                       = 1:size(RT_per_session.mean.all,2)
        subplot_number(k)                                       = subplot(4,8,subplot_assignment(k));
        hold on
        for idx_error_bar = start_label:numel(GLO.Labels)
            
            e_rt(idx_error_bar,k)                               = errorbar(idx_error_bar,RT_per_session.mean.all(idx_error_bar,k),RT_per_session.sem.all(idx_error_bar,k));
            if GLO.ttest_text
                text(idx_error_bar-0.2,RT_per_session.mean.all(idx_error_bar,k)+RT_per_session.sem.all(idx_error_bar,k)+0.01,asterisk.(['RT_' Conditions{k}])(idx_error_bar));
            end
            set(e_rt(idx_error_bar,k),'Marker','s','MarkerFaceColor',col(idx_error_bar,:),'MarkerEdgeColor',col(idx_error_bar,:),'Color',col(idx_error_bar,:),'Linewidth',GLO.linewidth,'MarkerSize',GLO.markersize);
            
            h_pre_children                                      = get(e_rt(idx_error_bar,k),'children');
            Xdata                                               = get(h_pre_children(2),'Xdata');
            tempo                                               = 4:3:length(Xdata);
            tempo(3:3:end)                                      = [];
            xleft                                               = tempo; xright = tempo+1;
            xmean                                               = (xleft+xright)/2;
            Xdata(xleft)                                        = [idx_error_bar-0.3 idx_error_bar-0.3];
            Xdata(xright)                                       = [idx_error_bar+0.3 idx_error_bar+0.3];
            set(h_pre_children(2),'Xdata',Xdata)
        end
        set(gca,'ylim',y_lim_RT_errorbars,'xlim',[start_label-1 numel(GLO.Labels)+1],'Xtick',start_label:numel(GLO.Labels),'Ytick',y_tick_RT,'xticklabel',GLO.Labels(start_label:end),'fontsize',GLO.fontsize_ticks_extra_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
        title(Labels_RT_comparison{k},'interpreter','none','fontsize',GLO.fontsize_titles_extra_small)
        xlabel('Stimulation onset [ms]','fontsize',GLO.fontsize_labels_extra_small)
        ylabel('Reaction time [ms]','fontsize',GLO.fontsize_labels_extra_small)
        box on
        
    end
    linkaxes([subplot_number(1) subplot_number(2) subplot_number(3) subplot_number(4)],'xy');
    
    
    
    %%
    
    % FIGURE 1, PLOT 4: REACTION TIME DIFFERENCE (4 PLOTS, 'Choice left','Choice right','Instructed left','Instructed right')
    %Labels_RT_comparison                                        = {'Choice contra','Choice ipsi','Instructed contra','Instructed ipsi'};
    subplot_assignment                                          = [11 12 9 10];
    % plot_matrix={idx_L_CH_S, idx_R_CH_S,  idx_L_IN_S, idx_R_IN_S; idx_L_CH_B,  idx_R_CH_B, idx_L_IN_B, idx_R_IN_B};
    %     ERROR BAR
    for k                                                       = 1:size(RT_per_session.mean.all,2)
        subplot_number(k)                                   = subplot(4,8,subplot_assignment(k));
        hold on
        for idx_error_bar                                       = start_label:numel(GLO.Labels)
            e_rt(idx_error_bar,k)                               = errorbar(idx_error_bar,RT_per_session.mean.all_diff(idx_error_bar,k),RT_per_session.sem.all_diff(idx_error_bar,k));
            if GLO.ttest_text
                text(idx_error_bar-0.2,RT_per_session.mean.all_diff(idx_error_bar,k)+RT_per_session.sem.all_diff(idx_error_bar,k)+0.01,asterisk.(['RT_difference_' Conditions{k}])(idx_error_bar));
            end
            set(e_rt(idx_error_bar,k),'Marker','s','MarkerFaceColor',col(idx_error_bar,:),'MarkerEdgeColor',col(idx_error_bar,:),'Color',col(idx_error_bar,:),'Linewidth',GLO.linewidth,'MarkerSize',GLO.markersize);
            
            h_pre_children                                      = get(e_rt(idx_error_bar,k),'children');
            Xdata                                               = get(h_pre_children(2),'Xdata');
            tempo                                               = 4:3:length(Xdata);
            tempo(3:3:end)                                      = [];
            xleft                                               = tempo; xright = tempo+1;
            xmean                                               =(xleft+xright)/2;
            Xdata(xleft)                                        = [idx_error_bar-0.3 idx_error_bar-0.3];
            Xdata(xright)                                       = [idx_error_bar+0.3 idx_error_bar+0.3];
            set(h_pre_children(2),'Xdata',Xdata)
        end
        set(gca,'ylim',y_lim_RT_diff_errorbars,'xlim',[start_label-1 numel(GLO.Labels)+1],'Xtick',start_label:numel(GLO.Labels),'Ytick',y_tick_RT_diff,'xticklabel',GLO.Labels(start_label:end),'fontsize',GLO.fontsize_ticks_extra_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
        
        %     title(Labels_RT_comparison{k},'interpreter','none','fontsize',GLO.fontsize_titles_extra_small)
        xlabel('Stimulation onset [ms]','fontsize',GLO.fontsize_labels_extra_small)
        ylabel('RT difference [ms]','fontsize',GLO.fontsize_labels_extra_small)
        line([0,numel(GLO.Labels)+1],[0,0],'Color',[.8 .8 .8],'LineStyle',':')
        set(gca, 'FontName', 'Arial','Linewidth',GLO.linewidth);
        box on
    end
    linkaxes([subplot_number(1) subplot_number(2) subplot_number(3) subplot_number(4)],'xy');
    
    % FIGURE 1, PLOT 5: EARLY VS LATE DIFFERENCE INDEX (Spatial preference difference to trial baseline early vs late per session)
    subplot(4,4,[3, 7])
    hold on
    title({'Spatial preference difference to trial baseline'; 'early vs late per session'},'fontsize',GLO.fontsize_titles_small);
    line([0,0],[-110,110])
    line([-110,110],[0,0])
    plot(bias_diff_stimulated_early,bias_diff_stimulated_late,'o','MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',GLO.markersize)
    if numel(bias_diff_stimulated_early) > 3 & all(~isnan(bias_diff_stimulated_early))
        u=convhull(bias_diff_stimulated_early,bias_diff_stimulated_late);
        fill (bias_diff_stimulated_early(u), bias_diff_stimulated_late(u), 'm', 'facealpha', alfa ,'LineStyle','none');
        plot (bias_diff_stimulated_early, bias_diff_stimulated_late,'o','MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',GLO.markersize)
    else
        plot (bias_diff_stimulated_early, bias_diff_stimulated_late, 'o','MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',GLO.markersize)
    end
    box on
    set(gca,'ylim',y_lim_bias_diff_indexes,'xlim',y_lim_bias_diff_indexes,'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
    axis('equal')
    axis([y_lim_bias_diff_indexes,y_lim_bias_diff_indexes])
    xlabel('Early contra target selection difference [%]','fontsize',GLO.fontsize_labels_small,'interpreter','none')
    ylabel('Late contra target selection difference [%]','fontsize',GLO.fontsize_labels_small,'interpreter','none')
    
    
    % FIGURE 1, PLOT 6: BASELINE VS STIMULATION INDEX (Spatial preference stimulated vs trial baseline per session)
    subplot(4,4,[4, 8])
    hold on
    title({'Spatial preference'; 'stimulated vs trial baseline per session'},'fontsize',GLO.fontsize_titles_small);
    line([0,0],[-10,110])
    line([-10,110],[0,0])
    if numel(bias_per_session_trial_baseline) > 3 & all(~isnan(bias_per_session_stimulated_early))
        p=convhull(bias_per_session_trial_baseline,bias_per_session_stimulated_early);
        fill (bias_per_session_trial_baseline(p), bias_per_session_stimulated_early(p), 'b', 'facealpha', alfa ,'LineStyle','none');
        g(1)=plot (bias_per_session_trial_baseline, bias_per_session_stimulated_early,'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',GLO.markersize);
        w=convhull(bias_per_session_trial_baseline,bias_per_session_stimulated_late);
        fill (bias_per_session_trial_baseline(w), bias_per_session_stimulated_late(w), 'r', 'facealpha', alfa ,'LineStyle','none');
        g(2)=plot (bias_per_session_trial_baseline, bias_per_session_stimulated_late,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',GLO.markersize);
    else
        g(1)=plot(bias_per_session_trial_baseline,bias_per_session_stimulated_early,'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',GLO.markersize);
        g(2)=plot(bias_per_session_trial_baseline,bias_per_session_stimulated_late,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',GLO.markersize);
    end
    box on
    line(y_lim_bias_indexes,y_lim_bias_indexes,'LineStyle',':')
    set(gca,'ylim',y_lim_bias_indexes,'xlim',y_lim_bias_indexes,'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
    axis('equal')
    axis([y_lim_bias_indexes,y_lim_bias_indexes])
    xlabel('Contraversive target selection baseline [%]','fontsize',GLO.fontsize_labels_small,'interpreter','none')
    ylabel('Contraversive target selection stimulated [%]','fontsize',GLO.fontsize_labels_small,'interpreter','none')
    L1=legend(g,{'Early','Late'});
    set(L1,'Position',[0.93 0.8 0.04 0.04]);
    
    
    % FIGURE 1, PLOT 7: REACTION TIME DIFFERENCE VS BIAS DIFFERENCE (Reaction time vs spatial preference difference to trial baseline)
    Labels_RT_comparison                                        = {'Choice contra','Choice ipsi','Instructed contra','Instructed ipsi'};
    %                       1               2               3               4
    
    %% !!!!!!!!!!!!!!!
    idx_for_RT_vs_bias=find(strcmp(Labels_RT_comparison,GLO.RT_vs_bias_condition));
    temp_RT_vs_bias_condition_diff                                          = [RT_per_session_diff{:,idx_for_RT_vs_bias}];
    temp_RT_vs_bias_condition                                               = [RT_per_session_all{:,idx_for_RT_vs_bias}];
    temp_RT_vs_bias_condition_mean                                          = nanmean(temp_RT_vs_bias_condition(:,2));
    
    subplot(4,4,[11,15])
    hold on
    title({'Reaction time vs spatial preference'; 'difference to trial baseline'},'fontsize',GLO.fontsize_titles_small);
    for i=start_label:size(bias_per_session_diff,2)
        plot(temp_RT_vs_bias_condition_diff(:,i)*1000,bias_per_session_diff(:,i),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'LineStyle','none','MarkerSize',GLO.markersize)
        counter=0;
        temp_x=[];temp_y=[];
        for t=1:numel(temp_RT_vs_bias_condition_diff(:,1))
            if isfinite(temp_RT_vs_bias_condition_diff(t,i)) && isfinite(bias_per_session_diff(t,i))
                counter=counter+1;
                temp_x(counter)=temp_RT_vs_bias_condition_diff(t,i)*1000; temp_y(counter)=bias_per_session_diff(t,i);
            end
        end
        if sum(temp_x)~=0;
            if numel(temp_x) > 3
                ch_idx=convhull(temp_x,temp_y);
                fill (temp_x(ch_idx), temp_y(ch_idx), col(i,:), 'facealpha', alfa ,'LineStyle','none');
                plot (temp_x, temp_y, 'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize)
            else
                plot (temp_x, temp_y, 'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize)
            end
        end
    end
    set(gca,'ylim',y_lim_bias_diff,'xlim',y_lim_RT_diff,'Xtick',y_tick_RT_diff,'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
    box on
    axis('square')
    line([0,0],y_lim_RT_diff,'Color',[.8 .8 .8],'LineStyle',':')
    line(y_lim_bias_diff,[0,0],'Color',[.8 .8 .8],'LineStyle',':')
    axis([y_lim_RT_diff, y_lim_bias_diff])
    xlabel('Reaction time difference [ms]','fontsize',GLO.fontsize_labels_small,'interpreter','none')
    ylabel('Contraversive target selection difference [%]','fontsize',GLO.fontsize_labels_small,'interpreter','none')
    
    % FIGURE 1, PLOT 8: REACTION TIME VS BIAS (Reaction time vs spatial preference)
    subplot(4,4,[12,16])
    hold on
    title({'Reaction time vs spatial preference';GLO.RT_vs_bias_condition},'fontsize',GLO.fontsize_titles_small);
    for i=start_label:size(bias_per_session,2)
        plot(temp_RT_vs_bias_condition(:,i)*1000,bias_per_session(:,i),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'LineStyle','none','MarkerSize',GLO.markersize)
        counter=0;
        temp_x=[];temp_y=[];
        for t=1:numel(temp_RT_vs_bias_condition(:,1))
            if isfinite(temp_RT_vs_bias_condition(t,i)) && isfinite(bias_per_session(t,i))
                counter=counter+1;
                temp_x(counter)=temp_RT_vs_bias_condition(t,i)*1000; temp_y(counter)=bias_per_session(t,i);
            end
        end
        if numel(temp_x) > 3 && numel(unique(temp_x))~=1 && numel(unique(temp_y))~=1
            k=convhull(temp_x,temp_y);
            fill (temp_x(k), temp_y(k), col(i,:), 'facealpha', alfa ,'LineStyle','none');
            h(i)=plot (temp_x, temp_y, 'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize);
        elseif isempty(temp_x) & isempty(temp_y)
            h(i)=plot(NaN);
        else
            h(i)=plot (temp_x, temp_y, 'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize);
        end
    end
    box on
    
    set(gca,'ylim',y_lim_bias,'xlim',y_lim_RT,'Xtick',y_tick_RT,'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
    axis('square')
    axis([y_lim_RT, y_lim_bias])
    
    xlabel('Reaction time [ms]','fontsize',GLO.fontsize_labels_small,'interpreter','none')
    ylabel('Contraversive target selection [%]','fontsize',GLO.fontsize_labels_small,'interpreter','none')
    line(y_lim_RT,[bias.mean.all(2) bias.mean.all(2)],'LineStyle','-')
    line([temp_RT_vs_bias_condition_mean temp_RT_vs_bias_condition_mean]*1000,y_lim_bias,'LineStyle','-')
    L2=legend(h(start_label:end),GLO.Labels(start_label:end));%,'Location','NorthEastOutside');
    set(L2,'Position',[0.93 0.3 0.04 0.04]);
    
    title_and_save(summary_1,plot_1_title,print_out)
end


 %% FIGURE 2 Bias and RT in two modes
 if any(ismember(GLO.plot_to_show,[2,-1]))
     plot_2_title='Summary 2- bias and RT in two modes';
     summary_2                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_2_title);
     
     struct_to_plot=bias;
     condition_subfields={'all','first_mode','second_mode'};
     condition_titles={'Bias (all)','During stimulation','After stimulation'};
     indexes_to_plot=start_label:numel(GLO.Labels);
     asterisk_labels={'bias_per_session','bias_per_session_firstmode','bias_per_session_secondmode'};
     xlim=[start_label-1 numel(GLO.Labels)+1];
     error_bar_plots(struct_to_plot, condition_subfields, condition_titles, indexes_to_plot, asterisk, asterisk_labels, [3,5],[1,6,11], xlim, y_lim_bias_errorbars, [0:10:100], 'Contraversive target selection [%]','small')
     
     struct_to_plot=RT_per_session;
     Labels_RT_comparison                                        = {'Choice contra','Choice ipsi','Instructed contra','Instructed ipsi'};
     Labels_RT_asterisks                                        = {'Choice_Left','Choice_Right','Instructed_Left','Instructed_Right'};
     subplot_number=[];
     subplot_assignments=[2,7,12;3,8,13;4,9,14;5,10,15];
     y_lims_RT=[y_lim_RT-[0 100];y_lim_RT_firstmode;y_lim_RT_secondmode];
     for condition=1:4
         
         asterisk_labels={['RT_' Labels_RT_asterisks{condition}],['RT_firstmode_' Labels_RT_asterisks{condition}],['RT_secondmode_' Labels_RT_asterisks{condition}],};
         condition_titles={[Labels_RT_comparison{condition} ' (all)'],['During stimulation'],['After stimulation']};
         sph=error_bar_plots(struct_to_plot, condition_subfields, condition_titles, indexes_to_plot, asterisk, asterisk_labels, [3,5],subplot_assignments(condition,:), xlim, y_lims_RT, [0:50:500], 'Reaction time [ms]','small',condition);
         subplot_number=[subplot_number sph];
     end
     linkaxes(subplot_number,'x')
     title_and_save(summary_2,plot_2_title,print_out);
 end
 
  %% FIGURE 2b Bias and RT in two modes
 if any(ismember(GLO.plot_to_show,[2,-1]))
     plot_2_title='Summary 2b- bias and RT differences to baseline in two modes';
     summary_2                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_2_title);
     
     struct_to_plot=bias;
     condition_subfields={'all_diff','first_mode_diff','second_mode_diff'};
     condition_titles={'Bias (all)','During stimulation','After stimulation'};
     indexes_to_plot=start_label:numel(GLO.Labels);
     asterisk_labels={'bias_per_session','bias_per_session_firstmode','bias_per_session_secondmode'};
     xlim=[start_label-1 numel(GLO.Labels)+1];
     error_bar_plots(struct_to_plot, condition_subfields, condition_titles, indexes_to_plot, asterisk, asterisk_labels, [3,5],[1,6,11], xlim, y_lim_bias_diff_errorbars, [-100:10:100], 'Contraversive target selection difference [%]','small')
     
     struct_to_plot=RT_per_session;
     Labels_RT_comparison                                        = {'Choice contra','Choice ipsi','Instructed contra','Instructed ipsi'};
     Labels_RT_asterisks                                        = {'Choice_Left','Choice_Right','Instructed_Left','Instructed_Right'};
     subplot_number=[];
     subplot_assignments=[2,7,12;3,8,13;4,9,14;5,10,15];
     y_lims_RT=[y_lim_RT;y_lim_RT_firstmode;y_lim_RT_secondmode];
     for condition=1:4
         
         asterisk_labels={['RT_difference_' Labels_RT_asterisks{condition}],['RT_firstmode_difference_' Labels_RT_asterisks{condition}],['RT_secondmode_difference_' Labels_RT_asterisks{condition}],};
         condition_titles={[Labels_RT_comparison{condition} ' (all)'],['During stimulation'],['After stimulation']};
         sph=error_bar_plots(struct_to_plot, condition_subfields, condition_titles, indexes_to_plot, asterisk, asterisk_labels, [3,5],subplot_assignments(condition,:), xlim, y_lim_RT_diff_errorbars, y_tick_RT_diff_errorbars, 'Reaction time difference [ms]','small',condition);
         subplot_number=[subplot_number sph];
     end
     linkaxes(subplot_number,'x')
     title_and_save(summary_2,plot_2_title,print_out);
 end
 
 %% FIGURE 3: SIGNIFICANT SESSIONS (Spatial preference modulation per session)     
if any(ismember(GLO.plot_to_show,[3,-1]))
    plot_3_title='Summary 3- Significance across sessions';
    summary_3                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_3_title);

    % FIGURE 3, PLOT 1: BIAS SIGNIFICANT SESSIONS (Spatial preference modulation per session)     
    struct_to_plot=batch_out.n_sig_bias;
    condition_subfields={'all'};
    condition_titles={'Spatial preference modulation per session'};
    indexes_to_plot=start_label_significance_plots:numel(GLO.Labels);
    xlim=[start_label_significance_plots-1 numel(GLO.Labels)+1];
    significance_plots(struct_to_plot, condition_subfields, condition_titles, indexes_to_plot, [1,2],1, xlim, y_lim_N_sessions, [-20:2:20], 'inc_dec', 'big')

    % FIGURE 3, PLOT 2: RTs SIGNIFICANT SESSIONS (Spatial preference modulation per session)  
    struct_to_plot=batch_out.n_sig_RTs;
    condition_subfields={'L_CH','R_CH','L_IN','R_IN'};
    condition_titles={'Choice contra RT','Choice ipsi RT','Instructed contra RT','Instructed ipsi RT'};
    indexes_to_plot=start_label_significance_plots:numel(GLO.Labels);
    xlim=[start_label_significance_plots-1 numel(GLO.Labels)+1];
    significance_plots(struct_to_plot, condition_subfields, condition_titles, indexes_to_plot, [2,4],[3 4 7 8], xlim, y_lim_N_sessions, [-20:2:20], 'inc_dec_RT', 'small')    
    title_and_save(summary_3,plot_3_title,print_out);
end



% if any(ismember(GLO.plot_to_show,[2,-1]))
%         % FIGURE 2, PLOT 1: BIAS SIGNIFICANT SESSIONS (Spatial preference modulation per session)
%     plot_2_title='Summary 2- Significance across sessions';
%     summary_2                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_2_title);
% 
%     subplot(1,2,1)
%     hold on
%     title('Spatial preference modulation per session','fontsize',GLO.fontsize_titles_small);
%     % shadedErrorBar(1:numel(bias.mean.all),bias.mean.all,bias.sem.all,{'r-o','markerfacecolor','r'})
%     for idx_error_bar                                           = start_label_significance_plots:numel(GLO.Labels)
%         bar(idx_error_bar,batch_out.n_sig_bias.plus.all(idx_error_bar),'EdgeColor',col(idx_error_bar,:),'FaceColor','w','Linewidth',GLO.linewidth)
%         bar(idx_error_bar,batch_out.n_sig_bias.minus.all(idx_error_bar)*-1,'EdgeColor',col(idx_error_bar,:),'FaceColor','w','Linewidth',GLO.linewidth)
%         bar(idx_error_bar,batch_out.n_sig_bias.plus_significant.all(idx_error_bar),'EdgeColor',col(idx_error_bar,:),'FaceColor',col(idx_error_bar,:),'Linewidth',GLO.linewidth)
%         bar(idx_error_bar,batch_out.n_sig_bias.minus_significant.all(idx_error_bar)*-1,'EdgeColor',col(idx_error_bar,:),'FaceColor',col(idx_error_bar,:),'Linewidth',GLO.linewidth)
%     end
%     fake_bar_edge=bar(1,NaN,'EdgeColor','k','FaceColor','w');
%     fake_bar_face=bar(1,NaN,'EdgeColor','k','FaceColor','k');
%     legend([fake_bar_edge fake_bar_face],{'non significant','significant'},'Location','NorthWest','fontsize',GLO.fontsize_legends_big-4)
%     
%     axis('square')
%     set(gca,'ylim',y_lim_N_sessions,'xlim',[start_label_significance_plots-1 numel(GLO.Labels)+1],'Xtick',start_label_significance_plots:numel(GLO.Labels),'xticklabel',GLO.Labels(start_label_significance_plots:end),'Ytick',-14:2:14,'yticklabel',[14:-2:0,2:2:14],'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
%     xlabel('Stimulation onset to target presentation [ms]','fontsize',GLO.fontsize_legends_big)
%     ylabel('Number of sessions','fontsize',GLO.fontsize_legends_big)
%     text(4.6,y_lim_N_sessions(2)-1,'Increased contraversive selection','fontsize',GLO.fontsize_text_big)
%     text(4.5,y_lim_N_sessions(1)+1,'Decreased contraversive selection','fontsize',GLO.fontsize_text_big)
%     line([0,numel(GLO.Labels)+1],[0,0],'Color',[.8 .8 .8],'LineStyle',':')
%     box on
%     
%     % FIGURE 2, PLOT 2: REACTION TIMES SIGNIFICANT SESSIONS (4 PLOTS, 'Choice left RT','Choice right RT','Instructed left RT','Instructed right RT')
%     Labels_RT_comparison                                        = {'Choice left RT','Choice right RT','Instructed left RT','Instructed right RT'};
%     RT_fields                                                   = {'L_CH','R_CH','L_IN','R_IN'};
%     
%     subplot_assignment                                          = [3 4 7 8];
%     % plot_matrix={idx_L_CH_S, idx_R_CH_S,  idx_L_IN_S, idx_R_IN_S; idx_L_CH_B,  idx_R_CH_B, idx_L_IN_B, idx_R_IN_B};
%     %     ERROR BAR
%     for k                                                       = 1: numel(RT_fields)
%         subplot_number(k)                                   = subplot(2,4,subplot_assignment(k));
%         hold on
%         for idx_error_bar                                           = start_label_significance_plots:numel(GLO.Labels)
%             bar(idx_error_bar,batch_out.n_sig_RTs.plus.(RT_fields{k})(idx_error_bar),'EdgeColor',col(idx_error_bar,:),'FaceColor','w','Linewidth',GLO.linewidth)
%             bar(idx_error_bar,batch_out.n_sig_RTs.minus.(RT_fields{k})(idx_error_bar)*-1,'EdgeColor',col(idx_error_bar,:),'FaceColor','w','Linewidth',GLO.linewidth)
%             bar(idx_error_bar,batch_out.n_sig_RTs.plus_significant.(RT_fields{k})(idx_error_bar),'EdgeColor',col(idx_error_bar,:),'FaceColor',col(idx_error_bar,:),'Linewidth',GLO.linewidth)
%             bar(idx_error_bar,batch_out.n_sig_RTs.minus_significant.(RT_fields{k})(idx_error_bar)*-1,'EdgeColor',col(idx_error_bar,:),'FaceColor',col(idx_error_bar,:),'Linewidth',GLO.linewidth)
%         end
%         switch GLO.effector_to_use
%             case 0
%                 set(gca,'ylim',y_lim_N_sessions,'xlim',[start_label_significance_plots-1 numel(GLO.Labels)+1],'Xtick',start_label_significance_plots:numel(GLO.Labels),'xticklabel',GLO.Labels(start_label_significance_plots:end),'Ytick',-14:2:14,'yticklabel',[14:-2:0,2:2:14],'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
%             case 4
%                 set(gca,'ylim',y_lim_N_sessions,'xlim',[start_label_significance_plots-1 numel(GLO.Labels)+1],'Xtick',start_label_significance_plots:numel(GLO.Labels),'xticklabel',GLO.Labels(start_label_significance_plots:end),'Ytick',-14:2:14,'yticklabel',[14:-2:0,2:2:14],'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
%         end
%         axis('square')
%         xlabel('Stimulation onset to target presentation [ms]','fontsize',GLO.fontsize_labels_extra_small)
%         ylabel('Number of sessions','fontsize',GLO.fontsize_labels_extra_small)
%         line([0,numel(GLO.Labels)+1],[0,0],'Color',[.8 .8 .8],'LineStyle',':')
%         box on
%         title(Labels_RT_comparison(k),'fontsize',GLO.fontsize_titles_small)
%         text(3.9,y_lim_N_sessions(2)-2,'Increased reaction time','fontsize',GLO.fontsize_labels_extra_small)
%         text(3.5,y_lim_N_sessions(1)+2,'Decreased reaction time','fontsize',GLO.fontsize_labels_extra_small)
%     end
%     linkaxes([subplot_number(1) subplot_number(2) subplot_number(3) subplot_number(4)],'xy');
%     title_and_save(summary_2,plot_2_title,print_out)    
% end


%% FIGURE 4 compiled RT Histograms
if any(ismember(GLO.plot_to_show,[4,-1]))
    
    plot_4_title='Summary 4- RT Histograms';
    summary_4                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_4_title);
    
    Labels_RT_comparison                                        = {'Choice left','Choice right','Instructed left','Instructed right'};
    Labels_RT_comparison_labels                                 = {'Choice contraversive','Choice ipsiversive','Instructed contraversive','Instructed ipsiversive'};
    subplot_assignment                                          = [3 4 1 2];
    area3options=struct('barwidth',0.0001,'ColormapWall',flipud(col(start_label:end,:)),...
        'ColormapRoof',flipud(col(start_label:end,:)),'FacealphaWall',1,...0.31,...
        'FacealphaRoof',1,'FacealphaFloor',0,...
        'SurfRoof',0,'LegendSym','none','LegendText',{fliplr(GLO.Labels)},... %'wall'
        'TScaling',[0.5 0.5; 0.5 0.5; 1 1],...%[0.2 0.8; 0.2 0.8; 0.2 0.8],
        'Edgecolor','k','Linestyle','-');%
    azimuth=0;%30;
    elevation=40;%25;
    X=(repmat(Bins.RThist',1,n_windows-start_label+1));
    for k=1:numel(subplot_assignment)
        subplot_handle(k)=subplot(2,2,subplot_assignment(k));
        hold on
        a3=area3(X',[start_label:n_windows],flipud(RThist{k}(:,start_label:end)'),area3options);
        
        %a3=area3(X',[1:n_windows],RThist{1}',area3options);
        set(gca,'ytick',[start_label:n_windows],'yticklabel',fliplr(GLO.Labels(start_label:end)),...
            'xtick',Bins.RThist(1:4:end),'xticklabel',Bins.RThist(1:4:end)*1000,'ztick',[0:10:100],'zticklabel',[0:10:100],'zlim',[0 50],'fontsize',GLO.fontsize_ticks_small)
        view(azimuth,elevation)
        
        xlabel('Reaction time [ms]','fontsize',GLO.fontsize_labels_small) %'Position',xlabelpos,
        ylabel('Stimulation onset [ms]','fontsize',GLO.fontsize_labels_small,'rotation',90)
        zlabel('            Saccade density [%]','fontsize',GLO.fontsize_labels_small)
        title(Labels_RT_comparison_labels{k},'interpreter','none','fontsize',GLO.fontsize_titles_small)
        %grid on
    end
    title_and_save(summary_4,plot_4_title,print_out)
end


%% FIGURE 5 Compiled RT errorbar (standard deviation) plots
if any(ismember(GLO.plot_to_show,[5,-1]))
    
    plot_5_title='Summary 5- Compiled rt error bars';
    summary_5                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_5_title);
    
    % FIGURE 5 plot 2: compiled first mode RTs
    Labels_RT_comparison                                        = {'Choice left','Choice right','Instructed left','Instructed right'};
    Labels_RT_comparison_labels                                 = {'Choice contra','Choice ipsi','Instructed contra','Instructed ipsi'};
    subplot_assignment                                          = [3 4 1 2];
    % plot_matrix={idx_L_CH_S, idx_R_CH_S,  idx_L_IN_S, idx_R_IN_S; idx_L_CH_B,  idx_R_CH_B, idx_L_IN_B, idx_R_IN_B};
    for k                                                       = 1:size(RT_per_session.mean.all,2)
        subplot_number(k)                                       = subplot(3,4,subplot_assignment(k));
        hold on
        for idx_error_bar = start_label:numel(GLO.Labels)
            
            e_rt(idx_error_bar,k)                               = errorbar(idx_error_bar,RT_raw_mean(idx_error_bar,k),RT_raw_std(idx_error_bar,k));
            if GLO.ttest_text
                text(idx_error_bar-0.2,RT_raw_mean(idx_error_bar,k)+RT_raw_std(idx_error_bar,k)+0.01,asterisks.(['RT_' Conditions{k}])(idx_error_bar));
            end
            set(e_rt(idx_error_bar,k),'Marker','s','MarkerFaceColor',col(idx_error_bar,:),'MarkerEdgeColor',col(idx_error_bar,:),'Color',col(idx_error_bar,:),'Linewidth',GLO.linewidth,'MarkerSize',GLO.markersize);
            
            h_pre_children                                      = get(e_rt(idx_error_bar,k),'children');
            Xdata                                               = get(h_pre_children(2),'Xdata');
            tempo                                               = 4:3:length(Xdata);
            tempo(3:3:end)                                      = [];
            xleft                                               = tempo; xright = tempo+1;
            xmean                                               = (xleft+xright)/2;
            Xdata(xleft)                                        = [idx_error_bar-0.3 idx_error_bar-0.3];
            Xdata(xright)                                       = [idx_error_bar+0.3 idx_error_bar+0.3];
            set(h_pre_children(2),'Xdata',Xdata)
        end
        set(gca,'ylim',y_lim_RT_comp_errorbars,'xlim',[start_label-1 numel(GLO.Labels)+1],'Xtick',start_label:numel(GLO.Labels),'Ytick',y_tick_RT_comp_errorbars,'xticklabel',GLO.Labels(start_label:end),'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
        title(Labels_RT_comparison_labels{k},'interpreter','none','fontsize',GLO.fontsize_titles_small)
        xlabel('Stimulation onset [ms]','fontsize',GLO.fontsize_labels_small)
        ylabel('Reaction time [ms]','fontsize',GLO.fontsize_labels_small)
        box on
        
    end
    linkaxes([subplot_number(1) subplot_number(2) subplot_number(3) subplot_number(4)],'xy');
    
    % FIGURE 5 plot 2: compiled first mode RTs
    subplot_assignment                                          = [7 8 5 6];
    % plot_matrix={idx_L_CH_S, idx_R_CH_S,  idx_L_IN_S, idx_R_IN_S; idx_L_CH_B,  idx_R_CH_B, idx_L_IN_B, idx_R_IN_B};
    for k                                                       = 1:size(RT_per_session.mean.all,2)
        subplot_number(k)                                       = subplot(3,4,subplot_assignment(k));
        hold on
        for idx_error_bar = start_label:numel(GLO.Labels)
            
            e_rt(idx_error_bar,k)                               = errorbar(idx_error_bar,RT_raw_mean_first_mode(idx_error_bar,k),RT_raw_std_first_mode(idx_error_bar,k));
            if GLO.ttest_text
                text(idx_error_bar-0.2,RT_raw_mean_first_mode(idx_error_bar,k)+RT_raw_std_first_mode(idx_error_bar,k)+0.01,asterisks.(['RT_firstmode_' Conditions{k}])(idx_error_bar));
            end
            set(e_rt(idx_error_bar,k),'Marker','s','MarkerFaceColor',col(idx_error_bar,:),'MarkerEdgeColor',col(idx_error_bar,:),'Color',col(idx_error_bar,:),'Linewidth',GLO.linewidth,'MarkerSize',GLO.markersize);
            
            h_pre_children                                      = get(e_rt(idx_error_bar,k),'children');
            Xdata                                               = get(h_pre_children(2),'Xdata');
            tempo                                               = 4:3:length(Xdata);
            tempo(3:3:end)                                      = [];
            xleft                                               = tempo; xright = tempo+1;
            xmean                                               = (xleft+xright)/2;
            Xdata(xleft)                                        = [idx_error_bar-0.3 idx_error_bar-0.3];
            Xdata(xright)                                       = [idx_error_bar+0.3 idx_error_bar+0.3];
            set(h_pre_children(2),'Xdata',Xdata)
        end
        set(gca,'ylim',y_lim_RT_comp_errorbars,'xlim',[start_label-1 numel(GLO.Labels)+1],'Xtick',start_label:numel(GLO.Labels),'Ytick',y_tick_RT_comp_errorbars,'xticklabel',GLO.Labels(start_label:end),'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
        title([Labels_RT_comparison_labels{k} ' during stimulation'],'interpreter','none','fontsize',GLO.fontsize_titles_small)
        xlabel('Stimulation onset [ms]','fontsize',GLO.fontsize_labels_small)
        ylabel('Reaction time [ms]','fontsize',GLO.fontsize_labels_small)
        box on
        
    end
    linkaxes([subplot_number(1) subplot_number(2) subplot_number(3) subplot_number(4)],'xy');
    
    % FIGURE 5 plot 3: compiled second mode RTs
    subplot_assignment                                          = [11 12 9 10];
    % plot_matrix={idx_L_CH_S, idx_R_CH_S,  idx_L_IN_S, idx_R_IN_S; idx_L_CH_B,  idx_R_CH_B, idx_L_IN_B, idx_R_IN_B};
    for k                                                       = 1:size(RT_per_session.mean.all,2)
        subplot_number(k)                                       = subplot(3,4,subplot_assignment(k));
        hold on
        for idx_error_bar = start_label:numel(GLO.Labels)
            
            e_rt(idx_error_bar,k)                               = errorbar(idx_error_bar,RT_raw_mean_second_mode(idx_error_bar,k),RT_raw_std_second_mode(idx_error_bar,k));
            if GLO.ttest_text
                text(idx_error_bar-0.2,RT_raw_mean_second_mode(idx_error_bar,k)+RT_raw_std_second_mode(idx_error_bar,k)+0.01,asterisks.(['RT_secondmode_' Conditions{k}])(idx_error_bar));
            end
            set(e_rt(idx_error_bar,k),'Marker','s','MarkerFaceColor',col(idx_error_bar,:),'MarkerEdgeColor',col(idx_error_bar,:),'Color',col(idx_error_bar,:),'Linewidth',GLO.linewidth,'MarkerSize',GLO.markersize);
            
            h_pre_children                                      = get(e_rt(idx_error_bar,k),'children');
            Xdata                                               = get(h_pre_children(2),'Xdata');
            tempo                                               = 4:3:length(Xdata);
            tempo(3:3:end)                                      = [];
            xleft                                               = tempo; xright = tempo+1;
            xmean                                               = (xleft+xright)/2;
            Xdata(xleft)                                        = [idx_error_bar-0.3 idx_error_bar-0.3];
            Xdata(xright)                                       = [idx_error_bar+0.3 idx_error_bar+0.3];
            set(h_pre_children(2),'Xdata',Xdata)
        end
        set(gca,'ylim',y_lim_RT_comp_errorbars,'xlim',[start_label-1 numel(GLO.Labels)+1],'Xtick',start_label:numel(GLO.Labels),'Ytick',y_tick_RT_comp_errorbars,'xticklabel',GLO.Labels(start_label:end),'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
        title([Labels_RT_comparison_labels{k} ' after stimulation'],'interpreter','none','fontsize',GLO.fontsize_titles_small)
        xlabel('Stimulation onset [ms]','fontsize',GLO.fontsize_labels_small)
        ylabel('Reaction time [ms]','fontsize',GLO.fontsize_labels_small)
        box on
    end
    linkaxes([subplot_number(1) subplot_number(2) subplot_number(3) subplot_number(4)],'xy');
    title_and_save(summary_5,plot_5_title,print_out)
end


%% FIGURE 6, Hitrates and Selections per position (4 PLOTS, 'Choice left RT','Choice right RT','Instructed left RT','Instructed right RT')
    
if any(ismember(GLO.plot_to_show,[6,-1]))
    plot_6_title='Summary 6- Selection per position and hitrate tables';
    summary_6                                                  = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_6_title);
    
    n_windows=numel(windows);
    subplot(3,3,1)
    scatter(real(positions),imag(positions),2000,'r','o')
    for pos=1:n_positions
        text(real(double(positions(pos))),imag(double(positions(pos))),num2str(pos),'fontsize',GLO.fontsize_text_big)
    end
  
    
    subplot(3,3,2)
    
    for pos=1:n_positions
        hold on
        current_hit_rate=sum(Choice_rate_per_target_mean(pos,:));
        for w=1:n_windows
            bar(pos,current_hit_rate,'FaceColor',col(w,:))
            current_hit_rate=current_hit_rate-Choice_rate_per_target_mean(pos,w);
        end
    end
    ylabel('Choice ratio','fontsize',GLO.fontsize_labels_small,'interpreter','none')
    
    
    subplot(3,3,3)
    
    for pos=1:n_positions
        hold on
        for w=1:n_windows
            bar(pos,Choice_rate_per_target_mean(pos,w)+n_windows-w,'FaceColor',col(w,:))
            bar(pos,n_windows-w,'FaceColor','w') ;
        end
    end
    ylabel('Choice ratio','fontsize',GLO.fontsize_labels_small,'interpreter','none')
    
    subplot(3,3,4)
    
    text(0.05,n_windows+2,'Stim onset', 'fontsize', GLO.fontsize_text_extra_small) ;
    text(0.2, n_windows+2,'mean left', 'fontsize', GLO.fontsize_text_extra_small) ;
    text(0.4, n_windows+2,'sem left',  'fontsize', GLO.fontsize_text_extra_small) ;
    text(0.6, n_windows+2,'mean right', 'fontsize', GLO.fontsize_text_extra_small) ;
    text(0.8, n_windows+2,'sem right',  'fontsize', GLO.fontsize_text_extra_small) ;
    for w=1:n_windows
        text(0.05,w,num2str(GLO.Labels{w}), 'fontsize', GLO.fontsize_text_extra_small) ;
        text(0.2, w,num2str((Hit_rate.mean.L_IN(w))), 'fontsize', GLO.fontsize_text_extra_small) ;
        text(0.4, w,num2str((Hit_rate.sem.L_IN(w))),  'fontsize', GLO.fontsize_text_extra_small) ;
        %text(0.2, w,num2str(round(Hit_rate.mean.L_IN(w))), 'fontsize', GLO.fontsize_text_extra_small) ;
        %text(0.4, w,num2str(round(Hit_rate.sem.L_IN(w))),  'fontsize', GLO.fontsize_text_extra_small) ;
        text(0.6, w,num2str((Hit_rate.mean.R_IN(w))), 'fontsize', GLO.fontsize_text_extra_small) ;
        text(0.8, w,num2str((Hit_rate.sem.R_IN(w))),  'fontsize', GLO.fontsize_text_extra_small) ;
        %text(0.6, w,num2str(round(Hit_rate.mean.R_IN(w))), 'fontsize', GLO.fontsize_text_extra_small) ;
        %text(0.8, w,num2str(round(Hit_rate.sem.R_IN(w))),  'fontsize', GLO.fontsize_text_extra_small) ;
    end
    title('Hitrate for instructed trials','fontsize', GLO.fontsize_titles_small)
    set(gca,'ylim',[0 n_windows+4],'xlim',[0 1],'Linewidth',GLO.linewidth);
    
    subplot(3,3,7)
    
    text(0.05,n_windows+2,'Stim onset', 'fontsize', GLO.fontsize_text_extra_small) ;
    text(0.2, n_windows+2,'mean left', 'fontsize', GLO.fontsize_text_extra_small) ;
    text(0.4, n_windows+2,'sem left',  'fontsize', GLO.fontsize_text_extra_small) ;
    text(0.6, n_windows+2,'mean right', 'fontsize', GLO.fontsize_text_extra_small) ;
    text(0.8, n_windows+2,'sem right',  'fontsize', GLO.fontsize_text_extra_small) ;
    for w=1:n_windows
        text(0.05,w,num2str(GLO.Labels{w}), 'fontsize', GLO.fontsize_text_extra_small) ;
        text(0.2, w,num2str(round(Hit_rate.mean.L_CH(w))), 'fontsize', GLO.fontsize_text_extra_small) ;
        text(0.4, w,num2str(round(Hit_rate.sem.L_CH(w))),  'fontsize', GLO.fontsize_text_extra_small) ;
        text(0.6, w,num2str(round(Hit_rate.mean.R_CH(w))), 'fontsize', GLO.fontsize_text_extra_small) ;
        text(0.8, w,num2str(round(Hit_rate.sem.R_CH(w))),  'fontsize', GLO.fontsize_text_extra_small) ;
    end
    title('Hitrate for choice trials','fontsize', GLO.fontsize_titles_small)
    set(gca,'ylim',[0 n_windows+4],'xlim',[0 1],'Linewidth',GLO.linewidth);
    
    
    subplot(3,3,5)
    for pos=1:n_positions
        hold on
        current_hit_rate=sum(Hit_rate_IN_per_target_mean(pos,:));
        for w=1:n_windows
            
            bar(pos,current_hit_rate,'FaceColor',col(w,:))
            %         bar(pos,current_hit_rate,'FaceColor',col(w,:))
            %         current_hit_rate=current_hit_rate-Hit_rate_IN_per_target_mean(pos,w);
            %        bar(pos,Hit_rate_IN_per_target_mean(pos,w)+n_windows-w,'FaceColor',col(w,:))
            %        bar(pos,n_windows-w,'FaceColor','w') ;
            current_hit_rate=current_hit_rate-Hit_rate_IN_per_target_mean(pos,w);
        end
    end
    ylabel('Hitrate instructed 1','fontsize',GLO.fontsize_labels_small,'interpreter','none')
    
    
    subplot(3,3,6)
    for pos=1:n_positions
        hold on
        for w=1:n_windows
            bar(pos,Hit_rate_IN_per_target_mean(pos,w)+n_windows-w,'FaceColor',col(w,:))
            bar(pos,n_windows-w,'FaceColor','w') ;
        end
    end
    ylabel('Hitrate instructed 2','fontsize',GLO.fontsize_labels_small,'interpreter','none')
    
    subplot(3,3,8)
    for pos=1:n_positions
        hold on 
        current_hit_rate=sum(Hit_rate_CH_per_target_mean(pos,~isnan(Hit_rate_CH_per_target_mean(pos,:)))) + sum(isnan(Hit_rate_CH_per_target_mean(pos,:)));
        for w=1:n_windows
            %         bar(pos,current_hit_rate,'FaceColor',col(w,:))
            %         current_hit_rate=current_hit_rate-Hit_rate_CH_per_target_mean(pos,w);
            %         bar(pos,Hit_rate_CH_per_target_mean(pos,w)+n_windows-w,'FaceColor',col(w,:))
            %         bar(pos,n_windows-w,'FaceColor','w') ;
            if ~isnan(Hit_rate_CH_per_target_mean(pos,w))
                bar(pos,current_hit_rate,'FaceColor',col(w,:))
                current_hit_rate=current_hit_rate-Hit_rate_CH_per_target_mean(pos,w);
            else
                bar(pos,current_hit_rate,'FaceColor','m')
                current_hit_rate=current_hit_rate-1;
            end
        end
    end
    ylabel('Hitrate choice 1','fontsize',GLO.fontsize_labels_small,'interpreter','none')
    
    
    subplot(3,3,9)
    for pos=1:n_positions
        hold on
        for w=1:n_windows
            if ~isnan(Hit_rate_CH_per_target_mean(pos,w))
            bar(pos,Hit_rate_CH_per_target_mean(pos,w)+n_windows-w,'FaceColor',col(w,:))
            else
             bar(pos,1+n_windows-w,'FaceColor','m')   
            end
            bar(pos,n_windows-w,'FaceColor','w') ;
        end
    end
    ylabel('Hitrate choice 2','fontsize',GLO.fontsize_labels_small,'interpreter','none')
    title_and_save(summary_6,plot_6_title,print_out)    
end



 %% FIGURE 7 Hitrates per condition (redundant?)
 if any(ismember(GLO.plot_to_show,[7,-1]))
     plot_7_title='Summary 7- Hitrates per condition';
     summary_7                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_7_title);
     
     struct_to_plot=batch_out.n_sig_hits;
     condition_subfields={'L_CH','R_CH','L_IN','R_IN'};
     condition_titles={'Choice contra Hitrate','Choice ipsi Hitrate','Instructed contra Hitrate','Instructed ipsi Hitrate'};
     indexes_to_plot=start_label_significance_plots:numel(GLO.Labels);
     xlim=[start_label_significance_plots-1 numel(GLO.Labels)+1];
     significance_plots(struct_to_plot, condition_subfields, condition_titles, indexes_to_plot, [2,4],[1 2 3 4], xlim, y_lim_N_sessions, [-20:2:20], 'inc_dec_hitrate', 'small')
     
     asterisk_labels={'Hit_rate_Choice_Left','Hit_rate_Choice_Right','Hit_rate_Instructed_Right','Hit_rate_Instructed_Right'};
     %condition_subfields={'L_CH','R_CH','L_IN','R_IN'};
     %condition_titles={'Choice left Hits','Choice right Hits','Instructed left Hits','Instructed right Hits'};
     indexes_to_plot=start_label:numel(GLO.Labels);
     xlim=[start_label-1 numel(GLO.Labels)+1];
     error_bar_plots(Hit_rate, condition_subfields, condition_titles, indexes_to_plot, asterisk, asterisk_labels, [2,4],[5 6 7 8], xlim, [80 101], [80:2:100], 'Hitrate [%]', 'small' );
     title_and_save(summary_7,plot_7_title,print_out);
 end
 

 
 %% FIGURE 8, Accuracy in Positions, 2D INSTRUCTED
 if any(ismember(GLO.plot_to_show,[8,-1]))
    plot_8_title='Summary 8a- Accuracy in Positions 2D Instructed';
    summary_8                                                  = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_8_title);
   
    angles=[0:pi/100:2*pi];
    circle_x=cos(angles);
    circle_y=sin(angles);
    stepsize_quadrants=pi/2;
    all_phis=[-pi:stepsize_quadrants:pi];
    stepsize_plot=pi/60;
    n_windows=numel(windows);
    
    for p=1:n_positions    
        center=[real(positions(p)) imag(positions(p))]; 
        
        hold on
        cl=line(5*circle_x+center(1),5*circle_y+center(2));
        set(cl,'Color','r','LineWidth',GLO.linewidth+2)
        cross=plot(center(1),center(2),'MarkerSize',10,'Marker','+','MarkerEdgeColor','r','MarkerFaceColor','r');
        deltaphi=atan2(center(2),center(1));
        
        %Accuracy_IN_per_target_rad{p,w}(s)
        
        for w=start_label:n_windows     
            current_color=GLO.col(w,:);
%             complex_mean=nanmean(Accuracy_IN_per_target_xy_raw{p,w});
%             real_std=std(real(Accuracy_IN_per_target_rad_raw{p,w})); 
%             imag_std=std(imag(Accuracy_IN_per_target_rad_raw{p,w})); 
            
%             
            complex_mean=nanmean(Accuracy_IN_per_target_xy{p,w});
            real_std=nanmean(real(Precision_IN_per_target_rad{p,w})); 
            imag_std=nanmean(imag(Precision_IN_per_target_rad{p,w})); 
            ellipse_r=real_std*imag_std./sqrt(imag_std.^2*cos(angles - deltaphi).^2+real_std.^2*sin(angles  - deltaphi).^2);
            ellipse_x = circle_x.*ellipse_r;
            ellipse_y = circle_y.*ellipse_r;
%             text(center(1)-sign(center(1))*6,center(2)+3-w/2,asterisks.(['AccuracyRadIN_pos_' num2str(p)]){w},...
%                 'fontsize',GLO.fontsize_text_small,'color',current_color);
            text(center(1)-sign(center(1))*6,center(2)+3-w/2,asterisk.(['AccuracyRadIN_pos_' num2str(p)]){w},...
                'fontsize',GLO.fontsize_text_small,'color',current_color);
            text(center(1)-sign(center(1))*9,center(2)+3-w/2,asterisk.(['PrecisionRadIN_pos_' num2str(p)]){w},...
                'fontsize',GLO.fontsize_text_small,'color',current_color);
            hold on
            %scatter(real(batch_out.per_tar_pos{s}.accuracy_xy_IN_raw{p,w})+center(1),imag(batch_out.per_tar_pos{s}.accuracy_xy_IN_raw{p,w})+center(2),20,'Marker','o','MarkerFaceColor',current_color)     
            sch=plot(real(complex_mean)+center(1),imag(complex_mean)+center(2),'MarkerSize',5,'Marker','.','MarkerEdgeColor',current_color,'Color',current_color);
            cl=line(ellipse_x+real(complex_mean)+center(1),ellipse_y+imag(complex_mean)+center(2));
            set(cl,'Color',current_color,'LineWidth',GLO.linewidth+1);%             
            
        end
    end   
    set(gca,'fontsize', GLO.fontsize_ticks_big);
    xlabel('Horizontal eye position [deg]','fontsize', GLO.fontsize_labels_big);
    ylabel('Vertical eye position [deg]', 'fontsize',GLO.fontsize_labels_big);    
    title('Saccade accuracy','fontsize', GLO.fontsize_titles_big);
    axis equal
    title_and_save(summary_8,plot_8_title,print_out);
 end

 %% FIGURE 8b, Accuracy in Positions, 2D CHOICE
 if any(ismember(GLO.plot_to_show,[8,-1]))
    plot_8_title='Summary 8b- Accuracy in Positions 2D Choice';
    summary_8                                                  = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_8_title);
   
    angles=[0:pi/100:2*pi];
    circle_x=cos(angles);
    circle_y=sin(angles);
    stepsize_quadrants=pi/2;
    all_phis=[-pi:stepsize_quadrants:pi];
    stepsize_plot=pi/60;
    n_windows=numel(windows);
    
    for p=1:n_positions    
        center=[real(positions(p)) imag(positions(p))]; 
        
        hold on
        cl=line(5*circle_x+center(1),5*circle_y+center(2));
        set(cl,'Color','r','LineWidth',GLO.linewidth+2)
        cross=plot(center(1),center(2),'MarkerSize',10,'Marker','+','MarkerEdgeColor','r','MarkerFaceColor','r');
        deltaphi=atan2(center(2),center(1))-pi/2;
        
        %Accuracy_IN_per_target_rad{p,w}(s)
        
        for w=start_label:n_windows     
            current_color=GLO.col(w,:);
%             complex_mean=nanmean(Accuracy_CH_per_target_xy_raw{p,w});
%             real_std=std(real(Accuracy_CH_per_target_rad_raw{p,w})); 
%             imag_std=std(imag(Accuracy_CH_per_target_rad_raw{p,w})); 
            
%             
            complex_mean=nanmean(Accuracy_CH_per_target_xy{p,w});
            real_std=nanmean(real(Precision_CH_per_target_rad{p,w})); 
            imag_std=nanmean(imag(Precision_CH_per_target_rad{p,w})); 
            ellipse_r=real_std*imag_std./sqrt(imag_std.^2*sin(angles - deltaphi).^2+real_std.^2*cos(angles  - deltaphi).^2);
            ellipse_x = circle_x.*ellipse_r;
            ellipse_y = circle_y.*ellipse_r;
%             text(center(1)-sign(center(1))*6,center(2)+3-w/2,asterisks.(['AccuracyRadIN_pos_' num2str(p)]){w},...
%                 'fontsize',GLO.fontsize_text_small,'color',current_color);
            text(center(1)-sign(center(1))*6,center(2)+3-w/2,asterisk.(['AccuracyRadCH_pos_' num2str(p)]){w},...
                'fontsize',GLO.fontsize_text_small,'color',current_color);
            text(center(1)-sign(center(1))*9,center(2)+3-w/2,asterisk.(['PrecisionRadCH_pos_' num2str(p)]){w},...
                'fontsize',GLO.fontsize_text_small,'color',current_color);
            hold on
                 
            sch=plot(real(complex_mean)+center(1),imag(complex_mean)+center(2),'MarkerSize',5,'Marker','.','MarkerEdgeColor',current_color,'Color',current_color);
            cl=line(ellipse_x+real(complex_mean)+center(1),ellipse_y+imag(complex_mean)+center(2));
            set(cl,'Color',current_color,'LineWidth',GLO.linewidth+1);%             
            
        end
    end   
    set(gca,'fontsize', GLO.fontsize_ticks_big);
    xlabel('Horizontal eye position [deg]','fontsize', GLO.fontsize_labels_big);
    ylabel('Vertical eye position [deg]', 'fontsize',GLO.fontsize_labels_big);    
    title('Saccade accuracy','fontsize', GLO.fontsize_titles_big);
    axis equal
    title_and_save(summary_8,plot_8_title,print_out);
 end

%  %% FIGURE 8c, Accuracy in Positions, 2D INSTRUCTED pooled
%  if any(ismember(GLO.plot_to_show,[8,-1]))
%     plot_8_title='Summary 8c- Accuracy in Positions 2D Instructed, pooled';
%     summary_8                                                  = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_8_title);
%    
%     angles=[0:pi/100:2*pi];
%     circle_x=cos(angles);
%     circle_y=sin(angles);
%     stepsize_quadrants=pi/2;
%     all_phis=[-pi:stepsize_quadrants:pi];
%     stepsize_plot=pi/60;
%     n_windows=numel(windows);
%     
%     for p=1:n_positions    
%         center=[real(positions(p)) imag(positions(p))]; 
%         
%         hold on
%         cl=line(5*circle_x+center(1),5*circle_y+center(2));
%         set(cl,'Color','r','LineWidth',GLO.linewidth+2)
%         cross=plot(center(1),center(2),'MarkerSize',10,'Marker','+','MarkerEdgeColor','r','MarkerFaceColor','r');
%         deltaphi=atan2(center(2),center(1))-pi/2;
%         
%         %Accuracy_IN_per_target_rad{p,w}(s)
%         
%         for w=start_label:n_windows     
%             current_color=GLO.col(w,:);
%             complex_mean=nanmean(Accuracy_IN_per_target_xy_raw{p,w});
%             real_std=std(real(Accuracy_IN_per_target_rad_raw{p,w})); 
%             imag_std=std(imag(Accuracy_IN_per_target_rad_raw{p,w})); 
%             
% %             
% %             complex_mean=nanmean(Accuracy_IN_per_target_xy{p,w});
% %             real_std=nanmean(real(Precision_IN_per_target_rad{p,w})); 
% %             imag_std=nanmean(imag(Precision_IN_per_target_rad{p,w})); 
%             ellipse_r=real_std*imag_std./sqrt(imag_std.^2*sin(angles - deltaphi).^2+real_std.^2*cos(angles  - deltaphi).^2);
%             ellipse_x = circle_x.*ellipse_r;
%             ellipse_y = circle_y.*ellipse_r;
%             text(center(1)-sign(center(1))*6,center(2)+3-w/2,asterisks.(['AccuracyRadIN_pos_' num2str(p)]){w},...
%                 'fontsize',GLO.fontsize_text_small,'color',current_color);
%             text(center(1)-sign(center(1))*9,center(2)+3-w/2,asterisk_pluses.(['AccuracyRadIN_pos_' num2str(p)]){w},...
%                 'fontsize',GLO.fontsize_text_small,'color',current_color);
% %             text(center(1)-sign(center(1))*6,center(2)+3-w/2,asterisk.(['AccuracyRadIN_pos_' num2str(p)]){w},...
% %                 'fontsize',GLO.fontsize_text_small,'color',current_color);
% %             text(center(1)-sign(center(1))*9,center(2)+3-w/2,asterisk.(['PrecisionIN_pos_' num2str(p)]){w},...
% %                 'fontsize',GLO.fontsize_text_small,'color',current_color);
%             hold on
%                  
%             sch=plot(real(complex_mean)+center(1),imag(complex_mean)+center(2),'MarkerSize',5,'Marker','.','MarkerEdgeColor',current_color,'Color',current_color);
%             cl=line(ellipse_x+real(complex_mean)+center(1),ellipse_y+imag(complex_mean)+center(2));
%             set(cl,'Color',current_color,'LineWidth',GLO.linewidth+1);%             
%             
%         end
%     end   
%     set(gca,'fontsize', GLO.fontsize_ticks_big);
%     xlabel('Horizontal eye position [deg]','fontsize', GLO.fontsize_labels_big);
%     ylabel('Vertical eye position [deg]', 'fontsize',GLO.fontsize_labels_big);    
%     title('Saccade accuracy','fontsize', GLO.fontsize_titles_big);
%     axis equal
%     title_and_save(summary_8,plot_8_title,print_out);
%  end
%  
%   %% FIGURE 8d, Accuracy in Positions, 2D CHOICE pooled
%  if any(ismember(GLO.plot_to_show,[8,-1]))
%     plot_8_title='Summary 8d- Accuracy in Positions 2D Choice, pooled';
%     summary_8                                                  = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_8_title);
%    
%     angles=[0:pi/100:2*pi];
%     circle_x=cos(angles);
%     circle_y=sin(angles);
%     stepsize_quadrants=pi/2;
%     all_phis=[-pi:stepsize_quadrants:pi];
%     stepsize_plot=pi/60;
%     n_windows=numel(windows);
%     
%     for p=1:n_positions    
%         center=[real(positions(p)) imag(positions(p))]; 
%         
%         hold on
%         cl=line(5*circle_x+center(1),5*circle_y+center(2));
%         set(cl,'Color','r','LineWidth',GLO.linewidth+2)
%         cross=plot(center(1),center(2),'MarkerSize',10,'Marker','+','MarkerEdgeColor','r','MarkerFaceColor','r');
%         deltaphi=atan2(center(2),center(1))-pi/2;
%         
%         %Accuracy_IN_per_target_rad{p,w}(s)
%         
%         for w=start_label:n_windows     
%             current_color=GLO.col(w,:);
%             complex_mean=nanmean(Accuracy_CH_per_target_xy_raw{p,w});
%             real_std=std(real(Accuracy_CH_per_target_rad_raw{p,w})); 
%             imag_std=std(imag(Accuracy_CH_per_target_rad_raw{p,w})); 
%             
% %             
% %             complex_mean=nanmean(Accuracy_CH_per_target_xy{p,w});
% %             real_std=nanmean(real(Precision_CH_per_target_rad{p,w})); 
% %             imag_std=nanmean(imag(Precision_CH_per_target_rad{p,w})); 
%             ellipse_r=real_std*imag_std./sqrt(imag_std.^2*sin(angles - deltaphi).^2+real_std.^2*cos(angles  - deltaphi).^2);
%             ellipse_x = circle_x.*ellipse_r;
%             ellipse_y = circle_y.*ellipse_r;
%             text(center(1)-sign(center(1))*6,center(2)+3-w/2,asterisks.(['AccuracyRadIN_pos_' num2str(p)]){w},...
%                 'fontsize',GLO.fontsize_text_small,'color',current_color);
%             text(center(1)-sign(center(1))*9,center(2)+3-w/2,asterisk_pluses.(['AccuracyRadIN_pos_' num2str(p)]){w},...
%                 'fontsize',GLO.fontsize_text_small,'color',current_color);
% %             text(center(1)-sign(center(1))*6,center(2)+3-w/2,asterisk.(['AccuracyRadCH_pos_' num2str(p)]){w},...
% %                 'fontsize',GLO.fontsize_text_small,'color',current_color);
% %             text(center(1)-sign(center(1))*9,center(2)+3-w/2,asterisk.(['PrecisionCH_pos_' num2str(p)]){w},...
% %                 'fontsize',GLO.fontsize_text_small,'color',current_color);
%             hold on
%                  
%             sch=plot(real(complex_mean)+center(1),imag(complex_mean)+center(2),'MarkerSize',5,'Marker','.','MarkerEdgeColor',current_color,'Color',current_color);
%             cl=line(ellipse_x+real(complex_mean)+center(1),ellipse_y+imag(complex_mean)+center(2));
%             set(cl,'Color',current_color,'LineWidth',GLO.linewidth+1);%             
%             
%         end
%     end   
%     set(gca,'fontsize', GLO.fontsize_ticks_big);
%     xlabel('Horizontal eye position [deg]','fontsize', GLO.fontsize_labels_big);
%     ylabel('Vertical eye position [deg]', 'fontsize',GLO.fontsize_labels_big);    
%     title('Saccade accuracy','fontsize', GLO.fontsize_titles_big);
%     axis equal
%     title_and_save(summary_8,plot_8_title,print_out);
%  end
 
 %% FIGURE 9 RT differences in two modes
 if any(ismember(GLO.plot_to_show,[9,-1]))
     plot_9_title='Summary 9- RT condition differences in two modes';
     summary_9                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_9_title);
     
     struct_to_plot=RT_per_session_condition_differences;
     indexes_to_plot=start_label:numel(GLO.Labels);
     condition_subfields={'all','first_mode','second_mode'};
     condition_titles={'bias (all)','bias in first mode','bias in second mode'};
     
     %condition_difference_labels                                           = {'I_IN_C_IN','I_CH_C_CH','I_CH_I_IN','C_CH_C_IN'};
     Labels_RT_comparison                                        = {'IP - CO (instructed)','IP - CO (choice)','CH - IN (ipsi)','CH - IN (contra)'};
     Labels_RT_asterisks                                         = {'R_IN_L_IN','R_CH_L_CH','R_CH_R_IN','L_CH_L_IN'};
     xlim=[start_label-1 numel(GLO.Labels)+1];
     subplot_number=[];
     subplot_assignments=[1,5,9;2,6,10;3,7,11;4,8,12];
     for condition=1:4
         
         asterisk_labels={['RT_bothmodes_' Labels_RT_asterisks{condition}],['RT_firstmode_' Labels_RT_asterisks{condition}],['RT_secondmode_' Labels_RT_asterisks{condition}],};
         condition_titles={[Labels_RT_comparison{condition} ' (all)'],[Labels_RT_comparison{condition} ' during stimulation'],[Labels_RT_comparison{condition} ' after stimulation']};
         sph=error_bar_plots(struct_to_plot, condition_subfields, condition_titles, indexes_to_plot, asterisk, asterisk_labels, [3,4],subplot_assignments(condition,:), xlim, [-30 30], [-50:10:50], 'Reaction Time difference [ms]',  'small',condition);
         subplot_number=[subplot_number sph];
     end
     linkaxes(subplot_number,'xy')
     title_and_save(summary_9,plot_9_title,print_out);
 end
 
 
 %% FIGURE 10: REACTION TIME DIFFERENCES VS LEFT TARGET SELECTION DIFFERENCES SEPERATED BY WINDOWS(4 PLOTS, 'Choice left RT','Choice right RT','Instructed left RT','Instructed right RT')  
 if any(ismember(GLO.plot_to_show,[10,-1]))
    start_label_window_by_window=3;
    plot_10_title='Summary 10- Bias vs rt seperated by windows';
    summary_10                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_10_title);
    condition_labels                                        = {'Choice left','Choice right','Instructed left','Instructed right'};
    %RT_conditions                                               = {'L_CH','R_CH','L_IN','R_IN'};
    
    plot_correlations_per_window(RT_per_session_diff,bias_per_session_diff,y_lim_RT_diff,y_lim_bias_diff,'Left selection difference [%]',condition_labels,batch_out,start_label_window_by_window,'extra_small')
    title_and_save(summary_10,plot_10_title,print_out)
 end
 
 %% FIGURE 10b: REACTION TIME DIFFERENCES VS LEFT TARGET SELECTION DIFFERENCES SEPERATED BY WINDOWS(4 PLOTS, 'Choice left RT','Choice right RT','Instructed left RT','Instructed right RT')    
 if any(ismember(GLO.plot_to_show,[10,-1]))
    start_label_window_by_window=3;
    plot_10_title='Summary 10b- Bias vs rt differences in conditions seperated by windows';
    summary_10                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_10_title);
    %Labels_RT_comparison                                        = {'Choice left','Choice right','Instructed left','Instructed right'};
    %RT_conditions                                                = {'L_CH','R_CH','L_IN','R_IN'};
    condition_labels                                             = {'R_IN_L_IN','R_CH_L_CH','R_CH_R_IN','L_CH_L_IN'};
    plot_correlations_per_window(RT_per_session_condition_differences_all,bias_per_session_diff,y_lim_RT_con_diff,y_lim_bias_diff,'Left selection difference [%]',condition_labels,batch_out,start_label_window_by_window,'extra_small')

    title_and_save(summary_10,plot_10_title,print_out)
 end
 
  
 %% FIGURE 100: REACTION TIME DIFFERENCES VS LEFT TARGET SELECTION DIFFERENCES SEPERATED BY WINDOWS(4 PLOTS, 'Choice left RT','Choice right RT','Instructed left RT','Instructed right RT')  
 if any(ismember(GLO.plot_to_show,[10,-1]))
    start_label_window_by_window=3;
    plot_10_title='Summary 100- Reaction time ipsi versus contra';
    summary_10                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_10_title);
    condition_labels                                        = {'Choice left','Choice right','Instructed left','Instructed right'};
    %RT_conditions                                               = {'L_CH','R_CH','L_IN','R_IN'};
    
    plot_correlations_per_window(RT_per_session_diff,RT_per_session_diff,y_lim_RT_diff,y_lim_RT_diff,'RT difference [ms]',condition_labels,batch_out,start_label_window_by_window,'extra_small')
    title_and_save(summary_10,plot_10_title,print_out)
 end
 
%% FIGURE 11: REACTION TIME VS BIAS (Reaction time vs spatial preference) (4 PLOTS, 'Choice left RT','Choice right RT','Instructed left RT','Instructed right RT')
 if any(ismember(GLO.plot_to_show,[11,-1]))
    plot_11_title='Summary 11- Correlations';
    summary_11                                                  = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_11_title);
    
    Labels_RT_comparison                                        = {'Choice left','Choice right','Instructed left','Instructed right'};
    subplot_assignment                                          = [7 8 3 4];
    for c=1:numel(subplot_assignment)
        subplot_handle(c)=subplot(2,4,subplot_assignment(c));
        temp_RT_per_session_mean                                                = RT_per_session.mean.all(:,c);
        temp_RT_per_session_baseline_mean                                       = RT_per_session.mean.all(2,c);
        %title({'Reaction time vs spatial preference';GLO.RT_vs_bias_condition},'fontsize',GLO.fontsize_titles_small);
        title(Labels_RT_comparison{c},'interpreter','none','fontsize',GLO.fontsize_titles_small)
        hold on
        
        for i=start_label:size(bias.mean.all,2)
            h(i)=plot(temp_RT_per_session_mean(i,1),bias.mean.all(1,i),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'LineStyle','none','MarkerSize',GLO.markersize);
        end
        pol = polyfit(temp_RT_per_session_mean(2:end),bias.mean.all(2:end)',1);
        refline(pol(1),pol(2))
        box on
        [r p]=corr([temp_RT_per_session_mean(3:end) bias.mean.all(3:end)'],'type','Spearman');
        set(gca,'ylim',y_lim_bias,'xlim',y_lim_RT,'Xtick',y_tick_RT,'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
        axis('square')
        axis([y_lim_RT, y_lim_bias])
        
        xlabel('Reaction time [ms]','fontsize',GLO.fontsize_labels_small,'interpreter','none')
        ylabel('Left target selection [%]','fontsize',GLO.fontsize_labels_small,'interpreter','none')
        line(y_lim_RT,[bias.mean.all(2) bias.mean.all(2)],'LineStyle','-')
        line([temp_RT_per_session_baseline_mean temp_RT_per_session_baseline_mean],y_lim_bias,'LineStyle','-')
        text(200,80,['r ' num2str(r(1,2)) '  p ' num2str(p(1,2))])
    end
    linkaxes([subplot_handle(1), subplot_handle(2), subplot_handle(3), subplot_handle(4)],'xy');
    L2=legend(h(start_label:end),GLO.Labels(start_label:end));%,'Location','NorthEastOutside');
    set(L2,'Position',[0.93 0.3 0.04 0.04]);
    title_and_save(summary_11,plot_11_title,print_out)    
end
 
%% FIGURE 12 : separation into early and late windows, correlation between RTs / RTs and bias 
if any(ismember(GLO.plot_to_show,[12,-1]))
    plot_12_title='Summary 12- Bias vs rt in early and late windows';
    summary_12                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_12_title);
    % PLOT 1: REACTION TIMES STIMULATED EARLY & LATE VS LEFT TARGET SEELCTION STIMULATED (4 PLOTS, 'Choice left RT','Choice right RT','Instructed left RT','Instructed right RT')
    hold on
    Labels_RT_comparison                                        = {'Choice contra','Choice ipsi','Instructed contra','Instructed ipsi'};
    subplot_assignment                                          = [5 6 1 2];
    % bias_per_session_stimulated_early
    % bias_per_session_trial_baseline
    for k=1:numel(subplot_assignment)
        temp_RT_vs_bias_condition=[RT_per_session_all{:,k}];
        temp_RT_vs_bias_condition_mean=nanmean(temp_RT_vs_bias_condition(:,2));
        subplot_handle(k)=subplot(2,4,subplot_assignment(k));
        hold on
        p_index=~isnan(bias_per_session_stimulated_early) & ~isnan(RT_mean_early(:,k));
            bias_plotted=bias_per_session_stimulated_early(p_index);
            RT_plotted=RT_mean_early(p_index,k)*1000;
        if sum(p_index) > 3 && numel(unique(bias_plotted))~=1 && numel(unique(RT_plotted))~=1
            p=convhull(RT_plotted,bias_plotted);
            fill (RT_plotted(p),bias_plotted(p), 'b', 'facealpha', alfa ,'LineStyle','none');
        end
        g(1)=plot(RT_mean_early(:,k)*1000,bias_per_session_stimulated_early,'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',GLO.markersize);
        w_index=~isnan(bias_per_session_stimulated_late) & ~isnan(RT_mean_late(:,k));
        
            bias_plotted=bias_per_session_stimulated_late(w_index);
            RT_plotted=RT_mean_late(w_index,k)*1000;
        if sum(w_index) > 3 && numel(unique(bias_plotted))~=1 && numel(unique(RT_plotted))~=1
            
            w=convhull(RT_plotted,bias_plotted);
            fill (RT_plotted(w),bias_plotted(w),  'r', 'facealpha', alfa ,'LineStyle','none');
        end
        g(2)=plot(RT_mean_late(:,k)*1000,bias_per_session_stimulated_late,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',GLO.markersize);
        box on
        
        set(gca,'xlim',y_lim_RT_indexes,'ylim',y_lim_bias_indexes,'Xtick',y_tick_RT_indexes,'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
        axis('square')
        axis([y_lim_RT_indexes,y_lim_bias_indexes])
        title(Labels_RT_comparison{k},'interpreter','none','fontsize',GLO.fontsize_titles_small)
        
        L1=legend(g,{'Early','Late'},'Location','SouthEast');
        % set(L1,'Position',[0.2 0.5 0.04 0.04]);
        
        ylabel('Contraversive target selection [%]','fontsize',GLO.fontsize_labels_small,'interpreter','none')
        xlabel('Reaction time [ms]','fontsize',GLO.fontsize_labels_small,'interpreter','none')
        
         line(y_lim_RT,[bias.mean.all(2) bias.mean.all(2)],'LineStyle','-')
         line([temp_RT_vs_bias_condition_mean temp_RT_vs_bias_condition_mean]*1000,y_lim_bias,'LineStyle','-')
    end
    linkaxes([subplot_handle(1), subplot_handle(2), subplot_handle(3), subplot_handle(4)],'xy');
    
    
    % FIGURE 12, PLOT 2: REACTION TIME VS BIAS (Reaction time vs spatial preference) (4 PLOTS, 'Choice left RT','Choice right RT','Instructed left RT','Instructed right RT')
    hold on
    Labels_RT_comparison                                        = {'Choice contra','Choice ipsi','Instructed contra','Instructed ipsi'};
    subplot_assignment                                          = [7 8 3 4];
    for c=1:numel(subplot_assignment)
        subplot_handle(c)=subplot(2,4,subplot_assignment(c));
        temp_RT_vs_bias_condition                                                = [RT_per_session_all{:,c}];
        temp_RT_vs_bias_condition_mean                                            = nanmean(temp_RT_vs_bias_condition(:,2));
        hold on
        %title({'Reaction time vs spatial preference';GLO.RT_vs_bias_condition},'fontsize',GLO.fontsize_titles_small);
        title(Labels_RT_comparison{c},'interpreter','none','fontsize',GLO.fontsize_titles_small)
        
        for i=start_label:size(bias_per_session,2)
            plot(temp_RT_vs_bias_condition(:,i)*1000,bias_per_session(:,i),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'LineStyle','none','MarkerSize',GLO.markersize)
            counter=0;
            temp_x=[];temp_y=[];
            for s=1:numel(temp_RT_vs_bias_condition(:,1))
                if isfinite(temp_RT_vs_bias_condition(s,i)) && isfinite(bias_per_session(s,i))
                    counter=counter+1;
                    temp_x(counter)=temp_RT_vs_bias_condition(s,i)*1000; temp_y(counter)=bias_per_session(s,i);
                    
                end
            end
            if numel(temp_x) > 3 && numel(unique(temp_x))~=1 && numel(unique(temp_y))~=1
                k=convhull(temp_x,temp_y);
                fill (temp_x(k), temp_y(k), col(i,:), 'facealpha', alfa ,'LineStyle','none');
                h(i)=plot (temp_x, temp_y, 'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize);
            elseif isempty(temp_x) & isempty(temp_y)
                h(i)=plot(NaN);
            else
                
                h(i)=plot (temp_x, temp_y, 'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize);
            end
        end
        box on
        
        set(gca,'ylim',y_lim_bias,'xlim',y_lim_RT,'Xtick',y_tick_RT,'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
        axis('square')
        axis([y_lim_RT, y_lim_bias])
        
        xlabel('Reaction time [ms]','fontsize',GLO.fontsize_labels_small,'interpreter','none')
        ylabel('Contraversive target selection [%]','fontsize',GLO.fontsize_labels_small,'interpreter','none')
        line(y_lim_RT,[bias.mean.all(2) bias.mean.all(2)],'LineStyle','-')
        line([temp_RT_vs_bias_condition_mean temp_RT_vs_bias_condition_mean]*1000,y_lim_bias,'LineStyle','-')
    end
    linkaxes([subplot_handle(1), subplot_handle(2), subplot_handle(3), subplot_handle(4)],'xy');
    L2=legend(h(start_label:end),GLO.Labels(start_label:end));%,'Location','NorthEastOutside');
    set(L2,'Position',[0.93 0.3 0.04 0.04]);
    title_and_save(summary_12,plot_12_title,print_out)
    
end

 %% FIGURE 13: REACTION TIME EARLY VS REACTION TIME LATE SEPARATED BY WINDOWS(4 PLOTS, 'Choice left RT','Choice right RT','Instructed left RT','Instructed right RT')  
 if any(ismember(GLO.plot_to_show,[13,-1]))
    start_label_window_by_window=3;
    plot_13_title='Summary 13- RT vs RT seperated by windows';
    summary_13                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_13_title);
    Labels_RT_comparison                                        = {'Choice contra','Choice ipsi','Instructed contra','Instructed ipsi'};
    condition_labels                                            = {'Choice','Choice','Instructed','Instructed'};
    side_labels                                                 = {'contra','ipsi','contra','ipsi'};
    RT_conditions                                               = {'L_CH','R_CH','L_IN','R_IN'};
    window_to_get_correlated                                    = '80';
    windows_to_correlate_to                                     = {'-120','-80'};
    %bias_and_RT_correlation_plot(RT_conditions,Labels_RT_comparison);
    
    for wc=1:numel(windows_to_correlate_to)
        subplot(numel(Labels_RT_comparison),numel(windows_to_correlate_to)+1,(wc-1)*(numel(windows_to_correlate_to)+1)+1);
        hold on
        i=find(strcmp(GLO.Labels,windows_to_correlate_to{wc}));
        y_idx=find(strcmp(GLO.Labels,window_to_get_correlated));
        if isempty(i) || isempty(y_idx)
            continue;
        end
        %bias_per_session_diff
        counter=0;
        temp_x=[];temp_y=[];significant_index_late=logical([]);significant_index_early=logical([]);significant_index_early_and_late=logical([]);non_significant_index=logical([]);
        for s=1:numel(bias_per_session_diff(:,1))
            if  isfinite(bias_per_session_diff(s,i)) && isfinite(bias_per_session_diff(s,y_idx))
                counter=counter+1;
                temp_x(counter)=bias_per_session_diff(s,i); temp_y(counter)=bias_per_session_diff(s,y_idx);
                tmp_significant_index_bias_Y                 =  batch_out.n_sig_bias.significant_per_session(s,y_idx);
                tmp_significant_index_bias_X                 =  batch_out.n_sig_bias.significant_per_session(s,i);
                significant_index_early_and_late(counter)   =  tmp_significant_index_bias_Y && tmp_significant_index_bias_X;
                significant_index_early(counter)            =  tmp_significant_index_bias_X && ~tmp_significant_index_bias_Y;
                significant_index_late(counter)             = ~tmp_significant_index_bias_X &&  tmp_significant_index_bias_Y;
                non_significant_index(counter)              = ~tmp_significant_index_bias_X && ~tmp_significant_index_bias_Y;
            end
        end
        
        
        if isempty(temp_x) && isempty(temp_y)
            h(wc)=plot(NaN,NaN,  'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize);
        else
            [RHO,PVAL] = corr(temp_x(:),temp_y(:),'type','Spearman');RHO=round(RHO*100)/100;PVAL=round(PVAL*100)/100;
            if sum(significant_index_early_and_late)~=0
                h(wc)=plot (temp_x(significant_index_early_and_late),  temp_y(significant_index_early_and_late),  'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize+4);
            end
            if sum(significant_index_early)~=0
                h(wc)=plot (temp_x(significant_index_early),  temp_y(significant_index_early),  '>','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize+4);
            end
            if sum(significant_index_late)~=0
                h(wc)=plot (temp_x(significant_index_late),  temp_y(significant_index_late),  '^','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize+4);
            end
            if sum(non_significant_index)~=0
                hns(wc)=plot (temp_x(non_significant_index), temp_y(non_significant_index), 'o','MarkerFaceColor','None',  'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize+4);
            end
            
        end
        x_and_y_nonnan = ~isnan(temp_x) & ~isnan(temp_y);
        if ~isempty(x_and_y_nonnan)
            fitcoeffs = polyfit(temp_x(x_and_y_nonnan), temp_y(x_and_y_nonnan), 1);
            fittedX = linspace(min(temp_x),max(temp_x), 200);
            fittedY = polyval(fitcoeffs, fittedX);
            plot(fittedX, fittedY,'color',col(i,:), 'LineWidth', GLO.linewidth +1);
        end
        box on
        
        set(gca,'fontsize',GLO.fontsize_ticks_extra_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
        axis('square')
        y_lim=[-30 90];
        x_lim=[-70 10];
        axis([x_lim, y_lim])
        %axis([-50,150, -40, 60])
        if c==1
            title(['Bias late versus early'],'interpreter','none','fontsize',GLO.fontsize_titles_extra_small)
        end
        ylabel(['Contraversive selection in ' GLO.Labels{y_idx} ' [%]'],'fontsize',GLO.fontsize_labels_extra_small,'interpreter','none')
        xlabel(['Contraversive selection in ' GLO.Labels{i} ' [%]'],'fontsize',GLO.fontsize_labels_extra_small,'interpreter','none')
        
        text(x_lim(1) + diff(x_lim)*0.03,y_lim(2)*0.9,['R = ' num2str(RHO)]);
        text(x_lim(1) + diff(x_lim)*0.54,y_lim(2)*0.9,['P = ' num2str(PVAL)]);
        
        line(x_lim,[0 0],'LineStyle','-')
        line([0 0],y_lim,'LineStyle','-')
    end
    
    for c=1:numel(Labels_RT_comparison)
        %cc=(mod(c,2)-0.5)*2+c;
        %temp_RT_vs_bias_condition_mean                                           = nanmean(temp_RT_vs_bias_condition(:,2));
        
        for k=1:numel(windows_to_correlate_to)
            temp_RT_vs_bias_condition_y=[RT_per_session_diff{:,c}];
            subplot_handle(c)=subplot(numel(Labels_RT_comparison),numel(windows_to_correlate_to)+1,(c-1)*(numel(windows_to_correlate_to)+1)+k+1);
            i=find(strcmp(GLO.Labels,windows_to_correlate_to{k}));
            %label_windows(k)=1;
            y_idx=find(strcmp(GLO.Labels,window_to_get_correlated));            
            if isempty(i) || isempty(y_idx)
               continue; 
            end
            
%             if mod(k+c,2)   
%                 temp_RT_vs_bias_condition                                               = [RT_per_session_diff{:,cc}];
%                 s_label=side_labels{cc};
%                 x_label={['dRT [ms]']};
%                 %x_label={Labels_RT_comparison{cc};['dRT in ' GLO.Labels{i} ' [ms]']};
%             else
                temp_RT_vs_bias_condition                                                = [RT_per_session_diff{:,c}];
                %s_label=side_labels{c};
                %x_label={['dRT [ms]']};
            %end
            hold on
            counter=0;
            temp_x=[];temp_y=[];significant_index_late=logical([]);significant_index_early=logical([]);significant_index_early_and_late=logical([]);non_significant_index=logical([]);
            for s=1:numel(temp_RT_vs_bias_condition(:,1))
                if isfinite(temp_RT_vs_bias_condition(s,i)) && isfinite(temp_RT_vs_bias_condition_y(s,y_idx))
                    counter=counter+1;
                    temp_x(counter)=temp_RT_vs_bias_condition(s,i)*1000; temp_y(counter)=temp_RT_vs_bias_condition_y(s,y_idx)*1000;
                    tmp_significant_index_RTs_Y                 =  batch_out.n_sig_RTs.significant_per_session.(RT_conditions{c})(s,y_idx);
                    tmp_significant_index_RTs_X                 =  batch_out.n_sig_RTs.significant_per_session.(RT_conditions{c})(s,i);
                    significant_index_early_and_late(counter)   =  tmp_significant_index_RTs_Y && tmp_significant_index_RTs_X;
                    significant_index_early(counter)            =  tmp_significant_index_RTs_X && ~tmp_significant_index_RTs_Y;
                    significant_index_late(counter)             = ~tmp_significant_index_RTs_X &&  tmp_significant_index_RTs_Y;
                    non_significant_index(counter)              = ~tmp_significant_index_RTs_X && ~tmp_significant_index_RTs_Y;
                end
            end
            
            if isempty(temp_x) && isempty(temp_y)
                h(k)=plot(NaN,NaN,  'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize);
            else
                [RHO,PVAL] = corr(temp_x(:),temp_y(:),'type','Spearman');RHO=round(RHO*100)/100;PVAL=round(PVAL*100)/100;
                if sum(significant_index_early_and_late)~=0
                    h(k)=plot (temp_x(significant_index_early_and_late),  temp_y(significant_index_early_and_late),  'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize+4);
                end
                if sum(significant_index_early)~=0
                    h(k)=plot (temp_x(significant_index_early),  temp_y(significant_index_early),  '>','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize+4);
                end
                if sum(significant_index_late)~=0
                    h(k)=plot (temp_x(significant_index_late),  temp_y(significant_index_late),  '^','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize+4);
                end
                if sum(non_significant_index)~=0
                    hns(k)=plot (temp_x(non_significant_index), temp_y(non_significant_index), 'o','MarkerFaceColor','None',  'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize+4);
                end
                                
            end
            x_and_y_nonnan = ~isnan(temp_x) & ~isnan(temp_y);
            if ~isempty(x_and_y_nonnan)
            fitcoeffs = polyfit(temp_x(x_and_y_nonnan), temp_y(x_and_y_nonnan), 1);
            fittedX = linspace(min(temp_x),max(temp_x), 200);
            fittedY = polyval(fitcoeffs, fittedX);
            plot(fittedX, fittedY,'color',col(i,:), 'LineWidth', GLO.linewidth +1);
            end
            box on
            
            y_lim=[-20 200];
            x_lim=[-70 10];
        
            set(gca,'fontsize',GLO.fontsize_ticks_extra_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
            axis('square')
            axis([x_lim, y_lim])
            %axis([-50,150, -40, 60])
            if c==1
                title([GLO.Labels{y_idx} ' vs ' GLO.Labels{i}],'interpreter','none','fontsize',GLO.fontsize_titles_extra_small)
            end
            if k==1
                ylabel({Labels_RT_comparison{c};['dRT in ' GLO.Labels{y_idx} ' [ms]']},'fontsize',GLO.fontsize_labels_extra_small,'interpreter','none')
            end
            if c==numel(Labels_RT_comparison)
%                 if mod(k,2)                    
%                     xlabel(['dRT same side in ' GLO.Labels{i} ' [ms]'],'fontsize',GLO.fontsize_labels_extra_small,'interpreter','none')
%                 else
%                     xlabel(['dRT opp side in ' GLO.Labels{i} ' [ms]'],'fontsize',GLO.fontsize_labels_extra_small,'interpreter','none')
%                 end
                    xlabel(['dRT ' GLO.Labels{i} ' [ms]'],'fontsize',GLO.fontsize_labels_extra_small,'interpreter','none');
            end
            
            text(x_lim(1) + diff(x_lim)*0.03,y_lim(2)*0.9,['R = ' num2str(RHO)]);
            text(x_lim(1) + diff(x_lim)*0.54,y_lim(2)*0.9,['P = ' num2str(PVAL)]);
                        
            line(x_lim,[0 0],'LineStyle','-')
            line([0 0],y_lim,'LineStyle','-')
        end
    end
    linkaxes([subplot_handle(1), subplot_handle(2), subplot_handle(3), subplot_handle(4)],'xy');
%     L2=legend(h(1:numel(windows_to_correlate_to)),GLO.Labels(label_windows));%,'Location','NorthEastOutside');
%     set(L2,'Position',[0.93 0.3 0.04 0.04]);
    title_and_save(summary_13,plot_13_title,print_out)
 end

 
  %% FIGURE 13b: REACTION TIME EARLY VS REACTION TIME LATE SEPARATED BY WINDOWS(4 PLOTS, 'Choice left RT','Choice right RT','Instructed left RT','Instructed right RT')  
 if any(ismember(GLO.plot_to_show,[13,-1]))
    start_label_window_by_window=3;
    plot_13_title='Summary 13b- RT vs RT condition differences seperated by windows';
    summary_13                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_13_title);
    
    condition_labels                                             = {'R_IN_L_IN','R_CH_L_CH','R_CH_R_IN','L_CH_L_IN'};   
    Labels_RT_comparison                                        = {'Choice contra','Choice ipsi','Instructed contra','Instructed ipsi'};
    RT_conditions                                               = {'L_CH','R_CH','L_IN','R_IN'};
    window_to_get_correlated                                    = '80';
    windows_to_correlate_to                                     = {'-120','-80'};
    for c=1:numel(condition_labels)
        for k=1:numel(windows_to_correlate_to)
            temp_RT_vs_bias_condition                                                = [RT_per_session_condition_differences_all{:,c}];
            subplot_handle(c)=subplot(numel(condition_labels),numel(windows_to_correlate_to),(c-1)*(numel(windows_to_correlate_to))+k);
            i=find(strcmp(GLO.Labels,windows_to_correlate_to{k}));
            label_windows(k)=1;
            y_idx=find(strcmp(GLO.Labels,window_to_get_correlated));
            if isempty(i) || isempty(y_idx)
               continue; 
            end
            hold on
            counter=0;
            temp_x=[];temp_y=[];significant_index_late=logical([]);significant_index_early=logical([]);significant_index_early_and_late=logical([]);non_significant_index=logical([]);
            for s=1:numel(temp_RT_vs_bias_condition(:,1))
                if isfinite(temp_RT_vs_bias_condition(s,i)) && isfinite(temp_RT_vs_bias_condition(s,y_idx))
                    counter=counter+1;
                    temp_x(counter)=temp_RT_vs_bias_condition(s,i)*1000; temp_y(counter)=temp_RT_vs_bias_condition(s,y_idx)*1000;
                    tmp_significant_index_RTs_Y                 =  batch_out.n_sig_RTs.significant_per_session.(RT_conditions{c})(s,y_idx);
                    tmp_significant_index_RTs_X                 =  batch_out.n_sig_RTs.significant_per_session.(RT_conditions{c})(s,i);
                    significant_index_early_and_late(counter)   =  tmp_significant_index_RTs_Y && tmp_significant_index_RTs_X;
                    significant_index_early(counter)            =  tmp_significant_index_RTs_X && ~tmp_significant_index_RTs_Y;
                    significant_index_late(counter)             = ~tmp_significant_index_RTs_X &&  tmp_significant_index_RTs_Y;
                    non_significant_index(counter)              = ~tmp_significant_index_RTs_X && ~tmp_significant_index_RTs_Y;
                end
            end
            
            if isempty(temp_x) && isempty(temp_y)
                h(k)=plot(NaN,NaN,  'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize);
            else
                [RHO,PVAL] = corr(temp_x(:),temp_y(:),'type','Spearman');RHO=round(RHO*100)/100;PVAL=round(PVAL*100)/100;
                if sum(significant_index_early_and_late)~=0
                    h(k)=plot (temp_x(significant_index_early_and_late),  temp_y(significant_index_early_and_late),  'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize+4);
                end
                if sum(significant_index_early)~=0
                    h(k)=plot (temp_x(significant_index_early),  temp_y(significant_index_early),  '>','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize+4);
                end
                if sum(significant_index_late)~=0
                    h(k)=plot (temp_x(significant_index_late),  temp_y(significant_index_late),  '^','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize+4);
                end
                if sum(non_significant_index)~=0
                    hns(k)=plot (temp_x(non_significant_index), temp_y(non_significant_index), 'o','MarkerFaceColor','None',  'MarkerEdgeColor',col(i,:),'MarkerSize',GLO.markersize+4);
                end
                                
            end
            x_and_y_nonnan = ~isnan(temp_x) & ~isnan(temp_y);
            if ~isempty(x_and_y_nonnan)
            fitcoeffs = polyfit(temp_x(x_and_y_nonnan), temp_y(x_and_y_nonnan), 1);
            fittedX = linspace(min(temp_x),max(temp_x), 200);
            fittedY = polyval(fitcoeffs, fittedX);
            plot(fittedX, fittedY,'color',col(i,:), 'LineWidth', GLO.linewidth +1);
            end
            box on
            
            set(gca,'fontsize',GLO.fontsize_ticks_extra_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
            axis('square')
            axis([y_lim_RT_con_diff, y_lim_RT_con_diff])
            %axis([-50,150, -40, 60])
            if c==1
                title([GLO.Labels{y_idx} ' vs ' GLO.Labels{i}],'interpreter','none','fontsize',GLO.fontsize_titles_extra_small)
            end
            if k==1
                ylabel({condition_labels{c};['dRT in ' GLO.Labels{y_idx} ' [ms]']},'fontsize',GLO.fontsize_labels_extra_small,'interpreter','none')
            end
            if c==numel(condition_labels)
                xlabel([GLO.Labels{i} ' [ms]'],'fontsize',GLO.fontsize_labels_extra_small,'interpreter','none')
            end
            
            text(y_lim_RT_con_diff(1) + diff(y_lim_RT_con_diff)*0.03,y_lim_RT_con_diff(2)*0.9,['R = ' num2str(RHO)]);
            text(y_lim_RT_con_diff(1) + diff(y_lim_RT_con_diff)*0.54,y_lim_RT_con_diff(2)*0.9,['P = ' num2str(PVAL)]);
                        
            line(y_lim_RT_con_diff,[0 0],'LineStyle','-')
            line([0 0],y_lim_RT_con_diff,'LineStyle','-')
        end
    end
    linkaxes([subplot_handle(1), subplot_handle(2), subplot_handle(3), subplot_handle(4)],'xy');
%     L2=legend(h(1:numel(windows_to_correlate_to)),GLO.Labels(label_windows));%,'Location','NorthEastOutside');
%     set(L2,'Position',[0.93 0.3 0.04 0.04]);
    title_and_save(summary_13,plot_13_title,print_out)
 end

%% FIGURE 14 compiled Velocity Histograms
if any(ismember(GLO.plot_to_show,[14,-1]))
    plot_14_title='Summary 14- Velocity Histograms';
    summary_14                                                  = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_14_title);
    
    Labels_RT_comparison                                        = {'Choice contra','Choice ipsi','Instructed contra','Instructed ipsi'};
    subplot_assignment                                          = [3 4 1 2];
    area3options=struct('barwidth',0.0001,'ColormapWall',flipud(col(start_label:end,:)),...
        'ColormapRoof',flipud(col(start_label:end,:)),'FacealphaWall',1,...0.31,...
        'FacealphaRoof',1,'FacealphaFloor',0,...
        'SurfRoof',0,'LegendSym','none','LegendText',{fliplr(GLO.Labels)},... %'wall'
        'TScaling',[0.5 0.5; 0.5 0.5; 1 1],...%[0.2 0.8; 0.2 0.8; 0.2 0.8],
        'Edgecolor','k','Linestyle','-');%
    azimuth=0;%30;
    elevation=40;%25;
    X=(repmat(Bins.Velhist',1,n_windows-start_label+1));
    for k=1:numel(subplot_assignment)
        subplot_handle(k)=subplot(2,2,subplot_assignment(k));
        hold on
        a3=area3(X',[start_label:n_windows],flipud(Velhist{k}(:,start_label:end)'),area3options);
        
        %a3=area3(X',[1:n_windows],RThist{1}',area3options);
        set(gca,'ytick',[start_label:n_windows],'yticklabel',fliplr(GLO.Labels(start_label:end)),...
            'xtick',Bins.Velhist,'xticklabel',Bins.Velhist,'ztick',[0:20:100],'zticklabel',[0:20:100],'zlim',[0 70],'fontsize',GLO.fontsize_ticks_small)
        view(azimuth,elevation)
        
        xlabel('Velocity [deg/sec]','fontsize',GLO.fontsize_labels_small) %'Position',xlabelpos,
        ylabel('Stimulation onset [ms]','fontsize',GLO.fontsize_labels_small,'rotation',90)
        zlabel('                        Velocity density [%]','fontsize',GLO.fontsize_labels_small)
        title(Labels_RT_comparison{k},'interpreter','none','fontsize',GLO.fontsize_titles_small)
        %grid on
    end
    title_and_save(summary_14,plot_14_title,print_out)    
end

 %% FIGURE 14b Velocity Errorbars
 if any(ismember(GLO.plot_to_show,[14,-1]))
     plot_14b_title='Summary 14b- Velocity errorbars in two modes';
     summary_14b                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_14b_title);
     
     %struct_to_plot=bias;
     condition_subfields={'all','first_mode','second_mode'};
     %condition_titles={'Bias (all)','During stimulation','After stimulation'};
     indexes_to_plot=start_label:numel(GLO.Labels);
     %asterisk_labels={'bias_per_session','bias_per_session_firstmode','bias_per_session_secondmode'};
     xlim=[start_label-1 numel(GLO.Labels)+1];
     %error_bar_plots(struct_to_plot, condition_subfields, condition_titles, indexes_to_plot, asterisk, asterisk_labels, [3,5],[1,6,11], xlim, y_lim_bias_errorbars, [0:10:100], 'Contraversive target selection [%]','small')
     
     struct_to_plot=Vel_per_session;
     Labels_RT_comparison                                        = {'Choice contra','Choice ipsi','Instructed contra','Instructed ipsi'};
     Labels_RT_asterisks                                        = {'Choice_Left','Choice_Right','Instructed_Left','Instructed_Right'};
     subplot_number=[];
     subplot_assignments=[1,5,9;2,6,10;3,7,11;4,8,12];
     y_lims_RT=[y_lim_RT;y_lim_RT_firstmode;y_lim_RT_secondmode];
     for condition=1:4
         
         asterisk_labels={['Vel_' Labels_RT_asterisks{condition}],['Vel_firstmode_' Labels_RT_asterisks{condition}],['Vel_secondmode_' Labels_RT_asterisks{condition}],};
         condition_titles={[Labels_RT_comparison{condition} ' (all)'],['During stimulation'],['After stimulation']};
         sph=error_bar_plots(struct_to_plot, condition_subfields, condition_titles, indexes_to_plot, asterisk, asterisk_labels, [3,4],subplot_assignments(condition,:), xlim, [600 1000], [0:100:2000], 'Velocity [deg/sec]','small',condition);
         subplot_number=[subplot_number sph];
     end
     linkaxes(subplot_number,'y')
     title_and_save(summary_14b,plot_14b_title,print_out);
 end


%% FIGURE 15 RT per condition over time
if any(ismember(GLO.plot_to_show,[15,-1]))
    plot_15_title='Summary 15- Scatterplot RTs over time';
    summary_15                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_15_title);
    Labels_RT_comparison                                        = {'Choice contra','Choice ipsi','Instructed contra','Instructed ipsi'};
    subplot_assignment                                          = [3 4 1 2];
    for k=1:numel(subplot_assignment)
        subplot_handle(k)=subplot(2,2,subplot_assignment(k));
         hold on
        for n_windows      = start_label:numel(GLO.Labels)
            h10=plot3(1:length(RT_raw{n_windows,k}),RT_raw{n_windows,k}*1000,repmat(n_windows,1,length(RT_raw{n_windows,k})),'o','MarkerFaceColor',col(n_windows,:),'MarkerEdgeColor',col(n_windows,:),'Color',col(n_windows,:));
        end
        title(Labels_RT_comparison{k},'interpreter','none','fontsize',GLO.fontsize_titles_small)
        set(gca,'ylim',y_lim_RT_scatterplot,'xlim',[0 max(max(cellfun(@(x)numel(x),RT_raw(start_label:end,:))))+1],'Ytick',y_tick_RT,'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
        xlabel('Trial','fontsize',GLO.fontsize_labels_small)
        ylabel('Reaction time [ms]','fontsize',GLO.fontsize_labels_small)
        box on
    end
    
    L10=legend(h10(start_label:end),GLO.Labels(start_label:end));%,'Location','NorthEastOutside');
    set(L10,'Position',[0.5 0.5 0.04 0.04]);
    linkaxes([subplot_handle(1), subplot_handle(2), subplot_handle(3), subplot_handle(4)],'xy');
    title_and_save(summary_15,plot_15_title,print_out)
end


% if any(ismember(GLO.plot_to_show,[11,-1]))
%     %% FIGURE 11 choices over time per batch
%     plot_11_title='Summary 11- Choices over time,';
%     summary_11                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_11_title);
%     idx_run=1;
%     for runniou= 1:length(batch_out.moving_average_t)
%         
%         
%         
%         Moving_average_L_CH_current_run                         = Moving_average_L_CH(:, runniou);
%         Moving_average_all_CH_current_run                       = Moving_average_all_CH(:, runniou);
%         
% %     batch_out.moving_average_L=[batch_out.moving_average_r{runniou}.L_baseline_choice batch_out.moving_average_t{runniou}.L_baseline_choice batch_out.moving_average_t{runniou}.L_stim_choice];
% %     batch_out.moving_average_all=[batch_out.moving_average_r{runniou}.all_baseline_choice batch_out.moving_average_t{runniou}.all_baseline_choice batch_out.moving_average_t{runniou}.all_stim_choice];
% 
%     window_size=100;
%     color_brown=[255 0 255]./255;
%     average_window=ones(1,window_size)/window_size;
% %     batch_out.moving_average_L_sum      =sum(vertcat(batch_out.moving_average_L{1,2:end}),1);
% %     batch_out.moving_average_all_sum    =sum(vertcat(batch_out.moving_average_all{1,2:end}),1);
% %     
%     
%     
%     batch_out.moving_average_L_sum      =sum(vertcat(Moving_average_L_CH_current_run{2:end}),1);
%     batch_out.moving_average_all_sum    =sum(vertcat(Moving_average_all_CH_current_run{2:end}),1);
% %    idx_non_nan_choices=[find(batch_out.moving_average_all_sum==1) numel(batch_out.moving_average_all_sum)];
% %     x_values_filter=[idx_non_nan_choices(1):idx_non_nan_choices(end)]';
%     filtered_choice_fraction=filter(average_window,1,batch_out.moving_average_L_sum)./filter(average_window,1,batch_out.moving_average_all_sum);       
%         x_values_fit=find(~isnan(filtered_choice_fraction))';
%     smooth_choices=smooth(filtered_choice_fraction,window_size/5);
% 
% 
%     fit_tryout=fit(x_values_fit,filtered_choice_fraction(x_values_fit)','poly1');
%     y_values_fit=feval(fit_tryout,x_values_fit);
%     subplot(3,length(batch_out.moving_average_t),idx_run)
%     hold on
%     h10a=plot3(1:length(filtered_choice_fraction),filtered_choice_fraction,repmat(numel(GLO.Labels),1,length(filtered_choice_fraction)),'Linestyle','--','MarkerFaceColor',color_brown,'MarkerEdgeColor',color_brown,'Color',color_brown,'Linewidth',GLO.linewidth+1);
%     subplot(3,length(batch_out.moving_average_t),length(batch_out.moving_average_t)+idx_run)
%     hold on
%     h10b=plot3(1:length(smooth_choices),smooth_choices,repmat(numel(GLO.Labels),1,length(smooth_choices)),'Linestyle','--','MarkerFaceColor',color_brown,'MarkerEdgeColor',color_brown,'Color',color_brown,'Linewidth',GLO.linewidth+1);
%     subplot(3,length(batch_out.moving_average_t),length(batch_out.moving_average_t)*2+idx_run)
%     hold on
%     h10c=plot3(x_values_fit,y_values_fit,repmat(numel(GLO.Labels),1,length(x_values_fit)),'MarkerFaceColor',color_brown,'MarkerEdgeColor',color_brown,'Linestyle','--','Color',color_brown,'Linewidth',GLO.linewidth+1);
%     
%     for n_windows      = start_label:numel(GLO.Labels)
%         
% %         L_choices=cumsum(batch_out.moving_average_L{1,n_windows})./cumsum(batch_out.moving_average_all{1,n_windows});
% %         L_choices=cumsum(batch_out.moving_average_L{1,n_windows})./cumsum(batch_out.moving_average_all{1,n_windows});
% %         idx_non_nan_choices=find(~isnan(L_choices));
% %        idx_non_nan_choices=[find(batch_out.moving_average_all{1,n_windows}) numel(batch_out.moving_average_all{1,n_windows})];
% %         x_values_filter=[idx_non_nan_choices(1):idx_non_nan_choices(end)]';
%         filtered_choice_fraction=filter(average_window,1,Moving_average_L_CH_current_run{n_windows})./filter(average_window,1,Moving_average_all_CH_current_run{n_windows});  
%         N_values_filter=length(Moving_average_L_CH_current_run{n_windows});
%         x_values_fit=find(~isnan(filtered_choice_fraction))';
%         smooth_choices=smooth(filtered_choice_fraction,window_size/5);
%         fit_tryout=fit(x_values_fit,filtered_choice_fraction(x_values_fit)','poly1');
%         y_values_fit=feval(fit_tryout,x_values_fit);
%         subplot(3,length(batch_out.moving_average_t),idx_run)
%         h10a=plot3(1:N_values_filter,filtered_choice_fraction,repmat(n_windows,1,N_values_filter),'Linestyle','-','MarkerFaceColor',col(n_windows,:),'MarkerEdgeColor',col(n_windows,:),'Color',col(n_windows,:),'Linewidth',GLO.linewidth+1);
%         xlabel('Trial','fontsize',GLO.fontsize_labels_extra_small)
%         ylabel('Filtered/smoothened left/all cumulative sum','fontsize',GLO.fontsize_labels_extra_small)
%         set(gca,'xlim',[0 length(batch_out.moving_average_all_sum)+10],'ylim',[-0.2 1.2],'Ytick',[-0.1:0.2:1.1],'Xtick',[0:25:N_values_filter+10],'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
%         subplot(3,length(batch_out.moving_average_t),length(batch_out.moving_average_t)+idx_run)
%         h10b=plot3(1:N_values_filter,smooth_choices,repmat(n_windows,1,N_values_filter),'Linestyle','-','MarkerFaceColor',col(n_windows,:),'MarkerEdgeColor',col(n_windows,:),'Color',col(n_windows,:),'Linewidth',GLO.linewidth+1);
%         xlabel('Trial','fontsize',GLO.fontsize_labels_extra_small)
%         ylabel('Filtered/smoothened left/all cumulative sum','fontsize',GLO.fontsize_labels_extra_small)
%         set(gca,'xlim',[0 length(batch_out.moving_average_all_sum)+10],'ylim',[-0.2 1.2],'Ytick',[-0.1:0.2:1.1],'Xtick',[0:25:N_values_filter+10],'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
%         subplot(3,length(batch_out.moving_average_t),length(batch_out.moving_average_t)*2+idx_run)
%         h10c=plot3(x_values_fit,y_values_fit,repmat(n_windows,1,length(x_values_fit)),'MarkerFaceColor',col(n_windows,:),'MarkerEdgeColor',col(n_windows,:),'Linestyle','-','Color',col(n_windows,:),'Linewidth',GLO.linewidth+1);
%         xlabel('Trial','fontsize',GLO.fontsize_labels_extra_small)
%         ylabel('9th degree polinomial regression ','fontsize',GLO.fontsize_labels_extra_small)
%         set(gca,'xlim',[0 length(batch_out.moving_average_all_sum)+10],'ylim',[-0.2 1.2],'Ytick',[-0.1:0.2:1.1],'Xtick',[0:25:N_values_filter+10],'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
%     end
% idx_run=idx_run+1;
% end
%     L10a=legend(h10a(start_label:end),['All choices' GLO.Labels(start_label:end)],'Location','NorthEast');
%     title_and_save(summary_11,plot_11_title,print_out)
% end


    %% FIGURE 16 Choices over time average
if any(ismember(GLO.plot_to_show,[16,-1]))
%     
%     Moving_average_L_CH=                 Moving_average_L_CH_per_run;
%     Moving_average_all_CH=                Moving_average_all_CH_per_run;
    %
    plot_16_title='Summary 16- Choices over time, average';
    summary_16                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_16_title);
    idx_run=1;
    
    max_N_trials=min(min(cellfun(@(x) numel(x),Moving_average_all_CH(2:end,:))));
    window_size=30;
    color_brown=[255 0 255]./255;
    average_window=ones(1,window_size)/window_size;
    
    for w = 1:numel(GLO.Labels)
        for runniou= 1:size(Moving_average_all_CH,2)
            current_N_trials=numel(Moving_average_all_CH{w,runniou});
            if current_N_trials<max_N_trials
                Moving_average_all_CH{w,runniou}(current_N_trials+1:max_N_trials)=false;
                Moving_average_L_CH{w,runniou}(current_N_trials+1:max_N_trials)=false;
            else
                Moving_average_all_CH{w,runniou}=Moving_average_all_CH{w,runniou}(1:max_N_trials);
                Moving_average_L_CH{w,runniou}=Moving_average_L_CH{w,runniou}(1:max_N_trials);
            end
            Filtered_L_moving_average{w,runniou}=filter(average_window,1,Moving_average_L_CH{w,runniou});
            Filtered_all_moving_average{w,runniou}=filter(average_window,1,Moving_average_all_CH{w,runniou});
        end
        Moving_average_across_sessions_L_CH{w}=nansum(vertcat(Moving_average_L_CH{w,:}),1);
        Moving_average_across_sessions_all_CH{w}=nansum(vertcat(Moving_average_all_CH{w,:}),1);
        Moving_average_across_sessions_mean_choice_fraction{w}=nanmean(double(vertcat(Moving_average_L_CH{w,:}))./double(vertcat(Moving_average_all_CH{w,:})),1);
        Moving_average_across_sessions_sem_choice_fraction{w}=sterr(double(vertcat(Filtered_L_moving_average{w,:}))./double(vertcat(Filtered_all_moving_average{w,:})),1);
    end
    for runniou= 1:length(batch_out.moving_average_t)
        Moving_average_all_windows_L_CH(runniou,:)=nansum(vertcat(Moving_average_L_CH{:,runniou}),1);
        Moving_average_all_windows_all_CH(runniou,:)=nansum(vertcat(Moving_average_all_CH{:,runniou}),1);
        Filtered_all_windows_L_moving_average(runniou,:)=filter(average_window,1,Moving_average_all_windows_L_CH(runniou,:));
        Filtered_all_windows_all_moving_average(runniou,:)=filter(average_window,1,Moving_average_all_windows_all_CH(runniou,:));
    end
    
    Moving_average_across_sessions_sem_all_windows_choice_fraction=sterr(Filtered_all_windows_L_moving_average./Filtered_all_windows_all_moving_average,1);
    %
    Moving_average_accross_sessions_L_sum      =nansum(vertcat(Moving_average_across_sessions_L_CH{2:end}),1);
    Moving_average_accross_sessions_all_sum    =nansum(vertcat(Moving_average_across_sessions_all_CH{2:end}),1);
    filtered_choices=filter(average_window,1,Moving_average_accross_sessions_all_sum);
    filtered_choice_fraction=filter(average_window,1,Moving_average_accross_sessions_L_sum)./filter(average_window,1,Moving_average_accross_sessions_all_sum);
    x_values_fit=find(~isnan(filtered_choice_fraction))';x_values_fit=x_values_fit(x_values_fit>window_size);
    smooth_choices=smooth(filtered_choice_fraction,window_size/2);
    fit_tryout=fit(x_values_fit,filtered_choice_fraction(x_values_fit)','poly1');
    y_values_fit=feval(fit_tryout,x_values_fit);
    subplot(3,2,1)
    hold on
    % h10a=plot3(1:length(filtered_choice_fraction),filtered_choice_fraction,repmat(numel(GLO.Labels),1,length(filtered_choice_fraction)),'Linestyle','--','MarkerFaceColor',color_brown,'MarkerEdgeColor',color_brown,'Color',color_brown,'Linewidth',GLO.linewidth+1);
    %h10a2=shadedErrorBar(1:length(filtered_choice_fraction),y,errBar,lineProps,transparent)
    h10a=shadedErrorBar(1:length(filtered_choice_fraction),filtered_choice_fraction,Moving_average_across_sessions_sem_all_windows_choice_fraction,{'Color',color_brown,'Linewidth',GLO.linewidth+1},1);
    subplot(3,2,3)
    hold on
    h10b=plot3(1:length(filtered_choices),filtered_choices,repmat(numel(GLO.Labels),1,length(filtered_choice_fraction)),'Linestyle','--','MarkerFaceColor',color_brown,'MarkerEdgeColor',color_brown,'Color',color_brown,'Linewidth',GLO.linewidth+1);
    subplot(3,2,5)
    hold on
    h10c=plot3(1:length(smooth_choices),smooth_choices,repmat(numel(GLO.Labels),1,length(smooth_choices)),'Linestyle','--','MarkerFaceColor',color_brown,'MarkerEdgeColor',color_brown,'Color',color_brown,'Linewidth',GLO.linewidth+1);
    subplot(1,2,2)
    hold on
    h10d=plot3(x_values_fit,y_values_fit,repmat(numel(GLO.Labels),1,length(x_values_fit)),'MarkerFaceColor',color_brown,'MarkerEdgeColor',color_brown,'Linestyle','--','Color',color_brown,'Linewidth',GLO.linewidth+1);
    
    for n_windows      = start_label:numel(GLO.Labels)
        filtered_choices=filter(average_window,1,Moving_average_across_sessions_all_CH{n_windows});
        filtered_choice_fraction=filter(average_window,1,Moving_average_across_sessions_L_CH{n_windows})./filter(average_window,1,Moving_average_across_sessions_all_CH{n_windows});
        N_values_filter=length(Moving_average_across_sessions_L_CH{n_windows});
        x_values_fit=find(~isnan(filtered_choice_fraction))';x_values_fit=x_values_fit(x_values_fit>window_size);
        smooth_choices=smooth(filtered_choice_fraction,window_size/2);
        fit_tryout=fit(x_values_fit,filtered_choice_fraction(x_values_fit)','poly1');
        y_values_fit=feval(fit_tryout,x_values_fit);
        
        subplot(3,2,1)
        %h10a=plot3(1:N_values_filter,filtered_choice_fraction,repmat(n_windows,1,N_values_filter),'Linestyle','-','MarkerFaceColor',col(n_windows,:),'MarkerEdgeColor',col(n_windows,:),'Color',col(n_windows,:),'Linewidth',GLO.linewidth+1);
        h10a=shadedErrorBar(1:N_values_filter,filtered_choice_fraction,Moving_average_across_sessions_sem_choice_fraction{n_windows},{'Color',col(n_windows,:),'Linewidth',GLO.linewidth+1},1);
        xlabel('Trial','fontsize',GLO.fontsize_labels_small)
        ylabel('Left/all choices','fontsize',GLO.fontsize_labels_small)
        title(['Filtered with window = ' num2str(window_size) ' trials'],'fontsize',GLO.fontsize_titles_extra_small)
        set(gca,'xlim',[window_size max_N_trials],'ylim',[-0.2 1.2],'Ytick',[-0.1:0.2:1.1],'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
        
        subplot(3,2,3)
        h10b=plot3(1:N_values_filter,filtered_choices,repmat(n_windows,1,N_values_filter),'Linestyle','-','MarkerFaceColor',col(n_windows,:),'MarkerEdgeColor',col(n_windows,:),'Color',col(n_windows,:),'Linewidth',GLO.linewidth+1);
        xlabel('Trial','fontsize',GLO.fontsize_labels_small)
        ylabel('All choices','fontsize',GLO.fontsize_labels_small)
        title(['Filtered with window = ' num2str(window_size) ' trials'],'fontsize',GLO.fontsize_titles_extra_small)
        set(gca,'xlim',[window_size max_N_trials],'ylim',[-0.2 1.2],'Ytick',[-0.1:0.2:1.1],'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
        
        subplot(3,2,5)
        h10c=plot3(1:N_values_filter,smooth_choices,repmat(n_windows,1,N_values_filter),'Linestyle','-','MarkerFaceColor',col(n_windows,:),'MarkerEdgeColor',col(n_windows,:),'Color',col(n_windows,:),'Linewidth',GLO.linewidth+1);
        xlabel('Trial','fontsize',GLO.fontsize_labels_small)
        ylabel('Left/all choices','fontsize',GLO.fontsize_labels_small)
        title(['Filtered with window = ' num2str(window_size) ', then smoothed with window = ' num2str(window_size/2) ' trials'],'fontsize',GLO.fontsize_titles_extra_small)
        set(gca,'xlim',[window_size max_N_trials],'ylim',[-0.2 1.2],'Ytick',[-0.1:0.2:1.1],'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
        
        subplot(1,2,2)
        h10d=plot3(x_values_fit,y_values_fit,repmat(n_windows,1,length(x_values_fit)),'MarkerFaceColor',col(n_windows,:),'MarkerEdgeColor',col(n_windows,:),'Linestyle','-','Color',col(n_windows,:),'Linewidth',GLO.linewidth+1);
        xlabel('Trial','fontsize',GLO.fontsize_labels_small)
        ylabel('Left/all choices','fontsize',GLO.fontsize_labels_small)
        title({'Linear regression'},'fontsize',GLO.fontsize_titles_small)
        set(gca,'xlim',[window_size max_N_trials],'ylim',[-0.2 1.2],'Ytick',[-0.1:0.2:1.1],'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
    end
    
    %     L10a=legend(h10a(start_label:end),['All choices' GLO.Labels(start_label:end)],'Location','NorthEast');
    title_and_save(summary_16,plot_16_title,print_out)
end

 %% FIGURE 17 Bias and RT in two modes
 if any(ismember(GLO.plot_to_show,[17,-1]))
     plot_17_title='Summary 17- Fraction of multiple saccades per session';
     summary_17                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_17_title);
%      
     struct_to_plot=Fraction_multistep_per_session;
     condition_subfields={'all'};
     Labels_RT_comparison                                        = {'Choice contra','Choice ipsi','Instructed contra','Instructed ipsi'};
     Labels_RT_asterisks                                        = {'Choice_Left','Choice_Right','Instructed_Left','Instructed_Right'};
     subplot_number=[];
     subplot_assignments=[1;2;3;4];
     xlim=[start_label-1 numel(GLO.Labels)+1];
     indexes_to_plot=start_label:numel(GLO.Labels);
     for condition=1:4
         
         asterisk_labels={['RT_' Labels_RT_asterisks{condition}]}; %% incorrect asterisks!!!
         condition_titles={[Labels_RT_comparison{condition}]};
         sph=error_bar_plots(struct_to_plot, condition_subfields, condition_titles, indexes_to_plot, asterisk, asterisk_labels, [2,2],subplot_assignments(condition,:), xlim, [0 1], [0:0.1:1], 'Fraction multistep','small',condition);
         subplot_number=[subplot_number sph];
         for w=start_label:numel(GLO.Labels)
             scatter(repmat(w,n_sessions,1),Fraction_multistep{w,condition},GLO.markersize,'MarkerEdgeColor',col(w,:),'marker','o')
         end
     end
     linkaxes(subplot_number,'x')
     title_and_save(summary_17,plot_17_title,print_out);
 end
 
 
 
 
 
 
 
 
 

if any(ismember(GLO.plot_to_show,[101,-1]))
        % FIGURE 101, PLOT 1: BIAS (Spatial preference modulation)
    plot_101_title='Summary 101- Bias_RT_session_per_session';
    summary_101                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_101_title);
    subplot(4,4,[09, 13])
    hold on
    title('Spatial preference modulation','fontsize',GLO.fontsize_titles_small);
%     shadedErrorBar(1:numel(bias.mean.all),bias.mean.all,bias.sem.all,{'r-o','markerfacecolor','r'})
copper_colors=copper(size(bias.per_session.all,1));
[~,colorsort]=sort(bias.per_session.all(:,2));
bias_to_plot=bias.per_session.all(colorsort,:);
%copper_colors=copper_colors(colorsort,:);
for idx=1:n_sessions
     plot(start_label:numel(GLO.Labels),bias_to_plot(idx,start_label:numel(GLO.Labels))','-','Color',copper_colors(idx,:))
end
 %    plot(start_label:numel(GLO.Labels),bias.per_session.all(:,start_label:numel(GLO.Labels))','-','Color',[0.5 0.5 0.5])
    
    for idx_error_bar                                           = start_label:numel(GLO.Labels)
%         e_bias(idx_error_bar)                                   = errorbar(idx_error_bar,bias.mean.all(idx_error_bar),bias.sem.all(idx_error_bar));

        e_bias_per_session(:,idx_error_bar)                                   = plot(idx_error_bar,bias.per_session.all(:,idx_error_bar),'-o');
%         if GLO.ttest_text
%             text(idx_error_bar-0.2,bias.mean.all(idx_error_bar)+bias.sem.all(idx_error_bar)+diff(y_lim_bias_errorbars)/50,asterisk.bias_per_session(idx_error_bar));
%         end
%         set(e_bias(idx_error_bar),'Marker','s','MarkerFaceColor',col(idx_error_bar,:),'MarkerEdgeColor',col(idx_error_bar,:),'Color',col(idx_error_bar,:),'Linewidth',GLO.linewidth,'MarkerSize',GLO.markersize);
         set(e_bias_per_session(:,idx_error_bar),'Marker','s','MarkerFaceColor',col(idx_error_bar,:),'MarkerEdgeColor',col(idx_error_bar,:),'Color',col(idx_error_bar,:),'Linewidth',GLO.linewidth,'MarkerSize',GLO.markersize);
        
%         h_pre_children                                          = get(e_bias(idx_error_bar),'children');
%         Xdata                                                   = get(h_pre_children(2),'Xdata');
%         tempo                                                   = 4:3:length(Xdata);
%         tempo(3:3:end) = [];
%         xleft                                                   = tempo; xright = tempo+1;
%         xmean                                                   = (xleft+xright)/2;
%         Xdata(xleft)                                            = [idx_error_bar-0.3 idx_error_bar-0.3];
%         Xdata(xright)                                           = [idx_error_bar+0.3 idx_error_bar+0.3];
%         set(h_pre_children(2),'Xdata',Xdata)
    end
    axis('square')
    set(gca,'ylim',[-10 110],'xlim',[start_label-1 numel(GLO.Labels)+1],'Xtick',start_label:numel(GLO.Labels),'xticklabel',GLO.Labels(start_label:end),'fontsize',GLO.fontsize_ticks_small,'FontName', 'Arial','Linewidth',GLO.linewidth);
    xlabel('Stimulation onset to target presentation [ms]','fontsize',GLO.fontsize_labels_small)
    ylabel('Contraversive target selection [%]','fontsize',GLO.fontsize_labels_small)
    box on
    
    title_and_save(summary_101,plot_101_title,print_out)
end

 
 
 
 
 
 
huha=1;

function subplot_number=error_bar_plots(struct_to_plot, condition_subfields, condition_titles, indexes_to_plot, asterisk, asterisk_labels, all_subplots,subplot_assignment, xlim, ylims, ytick, y_label, fontsizes, condition)
global GLO

for k=1:numel(condition_subfields)
    if numel(ylims)>3
        ylim=ylims(k,:);
    else
        ylim=ylims;
    end
    subplot_number(k) = subplot(all_subplots(1), all_subplots(2), subplot_assignment(k));
    hold on
    for w=indexes_to_plot
        if exist('condition')
            errorbar_handle(w)  = errorbar(w,struct_to_plot.mean.(condition_subfields{k})(w,condition),struct_to_plot.sem.(condition_subfields{k})(w,condition));
            
            if GLO.ttest_text
                text(w-diff(xlim)/30,struct_to_plot.mean.(condition_subfields{k})(w,condition)+struct_to_plot.sem.(condition_subfields{k})(w,condition)+diff(ylim)/50,asterisk.(asterisk_labels{k}){w},'fontsize',eval(['GLO.fontsize_text_' fontsizes]));
            end
        else
            errorbar_handle(w)  = errorbar(w,struct_to_plot.mean.(condition_subfields{k})(w),struct_to_plot.sem.(condition_subfields{k})(w));
            
            if GLO.ttest_text
                text(w-diff(xlim)/30,struct_to_plot.mean.(condition_subfields{k})(w)+struct_to_plot.sem.(condition_subfields{k})(w)+diff(ylim)/50,asterisk.(asterisk_labels{k}){w},'fontsize',eval(['GLO.fontsize_text_' fontsizes]));
            end
        end
        set(errorbar_handle(w),'Marker','s','MarkerFaceColor',GLO.col(w,:),'MarkerEdgeColor',GLO.col(w,:),'Color',GLO.col(w,:),'Linewidth',GLO.linewidth,'MarkerSize',GLO.markersize);
        h_pre_children                                          = get(errorbar_handle(w),'children');
        Xdata                                                   = get(h_pre_children(2),'Xdata');
        tempo                                                   = 4:3:length(Xdata);
        tempo(3:3:end) = [];
        xleft                                                   = tempo; xright = tempo+1;
        Xdata(xleft)                                            = [w-0.3 w-0.3];
        Xdata(xright)                                           = [w+0.3 w+0.3];
        set(h_pre_children(2),'Xdata',Xdata)
    end
    if exist('condition')
    plot(indexes_to_plot,struct_to_plot.mean.(condition_subfields{k})(indexes_to_plot,condition));
    else
    plot(indexes_to_plot,struct_to_plot.mean.(condition_subfields{k})(indexes_to_plot));        
    end
    %title([condition_titles{k}])
    
    set(gca,'ylim',ylim,'xlim',xlim,'Xtick',indexes_to_plot,'Ytick',ytick,'xticklabel',GLO.Labels(indexes_to_plot),'fontsize',eval(['GLO.fontsize_ticks_' fontsizes]),'FontName', 'Arial','Linewidth',GLO.linewidth);
    title(condition_titles{k},'interpreter','none','fontsize',eval(['GLO.fontsize_titles_' fontsizes]))
    xlabel('Stimulation onset [ms]','fontsize',eval(['GLO.fontsize_labels_' fontsizes]))
    ylabel(y_label,'fontsize',eval(['GLO.fontsize_labels_' fontsizes]))
    box on
end

linkaxes(subplot_number,'x');

function subplot_number=significance_plots(struct_to_plot, condition_subfields, condition_titles, indexes_to_plot, all_subplots,subplot_assignment, xlim, ylim, ytick, text_option, fontsizes)
global GLO

total_max=[];
for k                                                       = 1: numel(condition_subfields)
    total_max=max([total_max,diff([struct_to_plot.minus.(condition_subfields{k})*-1;struct_to_plot.plus.(condition_subfields{k})])]);
end
for k                                                       = 1: numel(condition_subfields)
    subplot_number(k)                                       = subplot(all_subplots(1),all_subplots(2),subplot_assignment(k));
    hold on
    fake_bar_edge=bar(-1,NaN,'EdgeColor','k','FaceColor','w');
    fake_bar_face=bar(-1,NaN,'EdgeColor','k','FaceColor','k');
    for w = indexes_to_plot
        hold on
        bar(w,struct_to_plot.plus.(condition_subfields{k})(w),'EdgeColor',GLO.col(w,:),'FaceColor','w','Linewidth',GLO.linewidth)
        bar(w,struct_to_plot.minus.(condition_subfields{k})(w)*-1,'EdgeColor',GLO.col(w,:),'FaceColor','w','Linewidth',GLO.linewidth)
        bar(w,struct_to_plot.plus_significant.(condition_subfields{k})(w),'EdgeColor',GLO.col(w,:),'FaceColor',GLO.col(w,:),'Linewidth',GLO.linewidth)
        bar(w,struct_to_plot.minus_significant.(condition_subfields{k})(w)*-1,'EdgeColor',GLO.col(w,:),'FaceColor',GLO.col(w,:),'Linewidth',GLO.linewidth)
    end

    ylim=[-total_max total_max]*1.2;
    set(gca,'ylim',ylim*1.1);
    ytick=get(gca,'ytick');
    ytick=unique(round(ytick));
    ytick=ytick(ytick>=-total_max & ytick<=total_max);    
    %ytick=-4:2:4;
    yticklabel=abs(ytick);
%     set(gca,'ylim',ylim*1.1,'xlim',xlim,'Xtick',indexes_to_plot,'xticklabel',GLO.Labels(indexes_to_plot),...
%         'Ytick',ytick,'yticklabel',[max(ytick):-2:0,2:2:max(ytick)],'fontsize',fontsize,'FontName', 'Arial','Linewidth',GLO.linewidth);
    set(gca,'ylim',ylim*1.1,'xlim',xlim,'Xtick',indexes_to_plot,'xticklabel',GLO.Labels(indexes_to_plot),...
        'Ytick',ytick,'yticklabel',yticklabel,'fontsize',eval(['GLO.fontsize_ticks_' fontsizes]),'FontName', 'Arial','Linewidth',GLO.linewidth);
     axis('square')
    xlabel('Stimulation onset [ms]','fontsize',eval(['GLO.fontsize_labels_' fontsizes]))
    ylabel('Number of sessions','fontsize',eval(['GLO.fontsize_labels_' fontsizes]))
    line([0,numel(GLO.Labels)+1],[0,0],'Color',[.8 .8 .8],'LineStyle',':')
    box on
    title(condition_titles{k},'fontsize',eval(['GLO.fontsize_titles_' fontsizes]),'Interpreter','None')
    switch text_option
        case 'inc_dec'
            text(indexes_to_plot(1) - 0.5,ylim(2),'Increased contraversive selection','fontsize',eval(['GLO.fontsize_text_' fontsizes]))
            text(indexes_to_plot(1) - 0.3,ylim(1),'Decreased contraversive selection','fontsize',eval(['GLO.fontsize_text_' fontsizes]))
            set(gca,'fontsize',GLO.fontsize_ticks_big)
            xlabel('Stimulation onset [ms]','fontsize',eval(['GLO.fontsize_labels_' fontsizes]))
            ylabel('Number of sessions','fontsize',eval(['GLO.fontsize_labels_' fontsizes]))
        case 'inc_dec_RT'
            text(indexes_to_plot(1) - 0.5,ylim(2),'Increased reaction time','fontsize',eval(['GLO.fontsize_text_' fontsizes]))
            text(indexes_to_plot(1) - 0.3,ylim(1),'Decreased reaction time','fontsize',eval(['GLO.fontsize_text_' fontsizes]))
        case 'inc_dec_hitrate'
            text(indexes_to_plot(1) - 0.5,ylim(2),'Increased hitrate','fontsize',eval(['GLO.fontsize_text_' fontsizes]))
            text(indexes_to_plot(1) - 0.3,ylim(1),'Decreased hitrate','fontsize',eval(['GLO.fontsize_text_' fontsizes]))
        case 'inc_dec_accuracy'
            text(indexes_to_plot(1) - 0.5,ylim(2),'Increased accuracy','fontsize',eval(['GLO.fontsize_text_' fontsizes]))
            text(indexes_to_plot(1) - 0.3,ylim(1),'Decreased accuracy','fontsize',eval(['GLO.fontsize_text_' fontsizes]))
        case 'bigger_smaller'
            text(indexes_to_plot(1) - 0.5,ylim(2),[condition_subfields{k}(1:4) ' > ' condition_subfields{k}(end-4:end)],'fontsize',eval(['GLO.fontsize_text_' fontsizes]),'Interpreter','none')
            text(indexes_to_plot(1) - 0.3,ylim(1),[condition_subfields{k}(1:4) ' < ' condition_subfields{k}(end-4:end)],'fontsize',eval(['GLO.fontsize_text_' fontsizes]),'Interpreter','none')
    end
    current_figure_position=get(gca,'position');
    legend_position=[current_figure_position(1) + current_figure_position(3)*0.48, current_figure_position(2) + current_figure_position(4)*0.25, current_figure_position(3)*0.4 current_figure_position(4)*0.1];
    legend([fake_bar_edge fake_bar_face],{'non significant','significant'},'Location','None','position',legend_position,'fontsize',eval(['GLO.fontsize_legends_' fontsizes]),'box','off','Xcolor',[1 1 1],'Ycolor',[1 1 1])
end
linkaxes(subplot_number,'xy');

function plot_correlations_per_window(X,Y,X_lim,Y_lim,Y_label,condition_labels,batch_out,start_label_window_by_window,fontsizes)
global GLO
RT_conditions                                                = {'L_CH','R_CH','L_IN','R_IN'};

for c=1:numel(condition_labels)
    temp_X_per_condition                                                = [X{:,c}];
    %temp_RT_vs_bias_condition_mean                                           = nanmean(temp_X_per_condition(:,2));
    cc=(mod(c,2)-0.5)*2+c;
    if iscell(Y) % not so cool
        temp_Y_per_condition                                                = [Y{:,c}].*1000;
        temp_X_per_condition                                                = [X{:,cc}];
        opposite_or_not='opposite side ';
    else
        temp_Y_per_condition                                                = Y;
        opposite_or_not='';
    end
    for i=start_label_window_by_window:GLO.n_windows
        subplot_handle(c)=subplot(numel(condition_labels),GLO.n_windows-2,(c-1)*(GLO.n_windows-2)+i-start_label_window_by_window+1);
        
        
        hold on
        counter=0;
        temp_x=[];temp_y=[];significant_index_bias=logical([]);significant_index_RTs=logical([]);significant_index_RTs_and_bias=logical([]);non_significant_index=logical([]);
        for s=1:numel(temp_X_per_condition(:,1))
            if isfinite(temp_X_per_condition(s,i)) && isfinite(temp_Y_per_condition(s,i))
                counter=counter+1;
                temp_x(counter)=temp_X_per_condition(s,i)*1000; temp_y(counter)=temp_Y_per_condition(s,i);
                if iscell(Y)
                    tmp_significant_index_Y      =batch_out.n_sig_RTs.significant_per_session.(RT_conditions{c})(s,i);
                    tmp_significant_index_X       =batch_out.n_sig_RTs.significant_per_session.(RT_conditions{cc})(s,i);%% not yet correct
                else
                    tmp_significant_index_Y      =batch_out.n_sig_bias.significant_per_session(s,i);
                    tmp_significant_index_X       =batch_out.n_sig_RTs.significant_per_session.(RT_conditions{c})(s,i);%% not yet correct                    
                end
                significant_index_RTs_and_bias(counter)=tmp_significant_index_Y && tmp_significant_index_X;
                significant_index_RTs(counter)  =    tmp_significant_index_X && ~tmp_significant_index_Y;
                significant_index_bias(counter) =   ~tmp_significant_index_X &&  tmp_significant_index_Y;
                non_significant_index(counter)  =   ~tmp_significant_index_X && ~tmp_significant_index_Y;
            end
        end
        
        if isempty(temp_x) && isempty(temp_y)
            h(i)=plot(NaN,NaN,  'o','MarkerFaceColor',GLO.col(i,:),'MarkerEdgeColor',GLO.col(i,:),'markersize',GLO.markersize);
        else
            [RHO,PVAL] = corr(temp_x(:),temp_y(:),'type','Spearman');RHO=round(RHO*100)/100;PVAL=round(PVAL*100)/100;
            if sum(significant_index_RTs_and_bias)~=0
                h(i)=plot (temp_x(significant_index_RTs_and_bias),  temp_y(significant_index_RTs_and_bias),  'o','MarkerFaceColor',GLO.col(i,:),'MarkerEdgeColor',GLO.col(i,:),'markersize',GLO.markersize+4);
            end
            if sum(significant_index_RTs)~=0
                h(i)=plot (temp_x(significant_index_RTs),  temp_y(significant_index_RTs),  '>','MarkerFaceColor',GLO.col(i,:),'MarkerEdgeColor',GLO.col(i,:),'markersize',GLO.markersize+4);
            end
            if sum(significant_index_bias)~=0
                h(i)=plot (temp_x(significant_index_bias),  temp_y(significant_index_bias),  '^','MarkerFaceColor',GLO.col(i,:),'MarkerEdgeColor',GLO.col(i,:),'markersize',GLO.markersize+4);
            end
            if sum(non_significant_index)~=0
                hns(i)=plot (temp_x(non_significant_index), temp_y(non_significant_index), 'o','MarkerFaceColor','None',  'MarkerEdgeColor',GLO.col(i,:),'markersize',GLO.markersize+4);
            end
            
            text(X_lim(1) + diff(X_lim)*0.03,Y_lim(2)*0.9,['R = ' num2str(RHO)]);
            text(X_lim(1) + diff(X_lim)*0.54,Y_lim(2)*0.9,['P = ' num2str(PVAL)]);
            
        end
        fitcoeffs = polyfit(temp_x, temp_y, 1);
        fittedX = linspace(min(temp_x),max(temp_x), 200);
        fittedY = polyval(fitcoeffs, fittedX);
        plot(fittedX, fittedY,'color',GLO.col(i,:), 'Linewidth', GLO.linewidth +1);
        box on
        
        set(gca,'fontsize',eval(['GLO.fontsize_ticks_' fontsizes]),'FontName', 'Arial','Linewidth',GLO.linewidth);
        axis('square')
        axis([X_lim, Y_lim])
        %axis([-50,150, -40, 60])
        if c==1
            title(GLO.Labels{i},'interpreter','none','fontsize',eval(['GLO.fontsize_titles_' fontsizes]))
        end
        if i==start_label_window_by_window
            ylabel({condition_labels{c};Y_label},'fontsize',eval(['GLO.fontsize_labels_' fontsizes]),'interpreter','none')
        end
        if c==numel(condition_labels)
            xlabel(['RT difference ' opposite_or_not '[ms]'],'fontsize',eval(['GLO.fontsize_labels_' fontsizes]),'interpreter','none')
        end
        
        line(X_lim,[0 0],'LineStyle','-')
        line([0 0],Y_lim,'LineStyle','-')
        
    end
end
linkaxes([subplot_handle(1), subplot_handle(2), subplot_handle(3), subplot_handle(4)],'xy');
%L2=legend(h(start_label_window_by_window:end),GLO.Labels(start_label_window_by_window:end));%,'Location','NorthEastOutside');
%set(L2,'Position',[0.93 0.3 0.04 0.04]);

function title_and_save(figure_handle,plot_title,print_out)
global GLO   
    
    mtit(figure_handle,  [plot_title, print_out], 'xoff', -0.0, 'yoff', 0.04, 'color', [0 0 0], 'fontsize', GLO.fontsize_titles_big,'Interpreter', 'none');
    stampit;
    %saveas(gcf,[GLO.folder_to_save plot_title print_out]);
    if GLO.create_pdf
        switch GLO.append_pdfs
            case 0
                export_fig([GLO.folder_to_save plot_title print_out], '-pdf','-transparent') % pdf by run
            case 1
                export_fig([GLO.folder_to_save plot_title, 'appended batches'], '-pdf', '-append','-transparent') % pdf by run
        end
    end