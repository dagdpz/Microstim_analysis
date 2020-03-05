%% NASTY part for getting percentage of ipsiversive fixation breaks for ipsi!! Fix -80cue
idx_error_temp=idx_stim_in_fix{k} & idx_error{k} & idx_instructed{k};
trial_ends=cellfun(@(x) x(:,end),{out_comp{k,1}.raw(idx_error_temp).time_axis},'UniformOutput',false);
aborted_after_stim=([out_comp{k,1}.task(idx_error_temp).stim_start]+ [out_comp{k,1}.states(idx_error_temp).start_mis])< [trial_ends{:}];
abort_states=vertcat(out_comp{k,1}.states(idx_error_temp).state_abo);
all_saccade_endpos_temp=vertcat(out_comp{k,1}.raw(idx_error_temp));
all_saccade_endpos=arrayfun(@(x) double(x.x_eye(end)),all_saccade_endpos_temp)- real(vertcat(out_comp{k,1}.saccades(idx_error_temp).fix_pos));
target_positions=vertcat(out_comp{k,1}.saccades(idx_error_temp).tar_pos) - repmat(vertcat(out_comp{k,1}.saccades(idx_error_temp).fix_pos),1,1);
all_saccade_endpos_first=all_saccade_endpos;

n_contra_minus80cue=sum(real(all_saccade_endpos_first(aborted_after_stim)<0) & ismember(abort_states(aborted_after_stim),[6]))
n_ipsi_minus80cue=sum(real(all_saccade_endpos_first(aborted_after_stim)>0) & ismember(abort_states(aborted_after_stim),[6]))
n_errors_stim_minus80cue=sum(aborted_after_stim)
all_errors_minus80cue = sum(idx_error_temp)


%% NASTY part for getting percentage of ipsiversive fixation breaks for ipsi!! Cue +80
idx_error_temp=idx_stim_in_cue{k} & idx_error{k} & idx_instructed{k};
trial_ends=cellfun(@(x) x(:,end),{out_comp{k,1}.raw(idx_error_temp).time_axis},'UniformOutput',false);
aborted_after_stim=([out_comp{k,1}.task(idx_error_temp).stim_start]+ [out_comp{k,1}.states(idx_error_temp).start_mis])< [trial_ends{:}];
abort_states=vertcat(out_comp{k,1}.states(idx_error_temp).state_abo);
all_saccade_endpos_temp=vertcat(out_comp{k,1}.raw(idx_error_temp));
all_saccade_endpos=arrayfun(@(x) double(x.x_eye(end)),all_saccade_endpos_temp)- real(vertcat(out_comp{k,1}.saccades(idx_error_temp).fix_pos));
target_positions=vertcat(out_comp{k,1}.saccades(idx_error_temp).tar_pos) - repmat(vertcat(out_comp{k,1}.saccades(idx_error_temp).fix_pos),1,1);
all_saccade_endpos_first=all_saccade_endpos;

n_contra_plus80cue=sum(real(all_saccade_endpos_first(aborted_after_stim)<0) & ismember(abort_states(aborted_after_stim),[6]))
n_ipsi_plus80cue=sum(real(all_saccade_endpos_first(aborted_after_stim)>0) & ismember(abort_states(aborted_after_stim),[6]))
n_errors_stim_plus80cue=sum(aborted_after_stim)
all_errors_plus80cue = sum(idx_error_temp)



%% NASTY part for getting percentage of contraversice undershooting errors!! -80GO
idx_error_80_after_go=idx_stim_in_mem{k} & idx_error{k} & idx_instructed{k};
trial_ends=cellfun(@(x) x(:,end),{out_comp{k,1}.raw(idx_error_80_after_go).time_axis},'UniformOutput',false);
aborted_after_stim=([out_comp{k,1}.task(idx_error_80_after_go).stim_start]+ [out_comp{k,1}.states(idx_error_80_after_go).start_mis])< [trial_ends{:}];
abort_states_80_after_go=vertcat(out_comp{k,1}.states(idx_error_80_after_go).state_abo);
all_saccade_endpos_after_go=vertcat(out_comp{k,1}.saccades(idx_error_80_after_go).endpos_obs) - repmat(vertcat(out_comp{k,1}.saccades(idx_error_80_after_go).fix_pos),1,5);
target_positions_after_go=vertcat(out_comp{k,1}.saccades(idx_error_80_after_go).tar_pos) - repmat(vertcat(out_comp{k,1}.saccades(idx_error_80_after_go).fix_pos),1,1);
all_saccade_endpos_after_go_first=all_saccade_endpos_after_go(:,1);

n_contra_minus80go=sum(real(all_saccade_endpos_after_go_first(aborted_after_stim)<0) & ismember(abort_states_80_after_go(aborted_after_stim),[9,10]))
n_ipsi_minus80go=sum(real(all_saccade_endpos_after_go_first(aborted_after_stim)>0) & ismember(abort_states_80_after_go(aborted_after_stim),[9,10]))
n_errors_stim_minus80go=sum(aborted_after_stim)
all_errors_minus80go = sum(idx_error_80_after_go)



%% NASTY part for getting percentage of contraversice undershooting errors!! +80 GO
idx_error_80_after_go=[out_comp{k,1}.task.stim_start]==0.08 & idx_stim_in_acq_inv{k} & idx_error{k} & idx_instructed{k};
trial_ends=cellfun(@(x) x(:,end),{out_comp{k,1}.raw(idx_error_80_after_go).time_axis},'UniformOutput',false);
aborted_after_stim=([out_comp{k,1}.task(idx_error_80_after_go).stim_start]+ [out_comp{k,1}.states(idx_error_80_after_go).start_mis])< [trial_ends{:}];
abort_states_80_after_go=vertcat(out_comp{k,1}.states(idx_error_80_after_go).state_abo);
all_saccade_endpos_after_go=vertcat(out_comp{k,1}.saccades(idx_error_80_after_go).endpos_obs) - repmat(vertcat(out_comp{k,1}.saccades(idx_error_80_after_go).fix_pos),1,5);
target_positions_after_go=vertcat(out_comp{k,1}.saccades(idx_error_80_after_go).tar_pos) - repmat(vertcat(out_comp{k,1}.saccades(idx_error_80_after_go).fix_pos),1,1);
all_saccade_endpos_after_go_first=all_saccade_endpos_after_go(:,1);

n_contra_plus80go=sum(real(all_saccade_endpos_after_go_first(aborted_after_stim)<0) & ismember(abort_states_80_after_go(aborted_after_stim),[9,10]))
n_ipsi_plus80go=sum(real(all_saccade_endpos_after_go_first(aborted_after_stim)>0) & ismember(abort_states_80_after_go(aborted_after_stim),[9,10]))
n_errors_stim_plus80go=sum(aborted_after_stim)
all_errors_plus80go = sum(idx_error_80_after_go)






