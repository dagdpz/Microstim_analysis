%Sel_all                         = {'display',1,'summary',0,'keep_raw_data',1,'success',0,'aborted_state',9,'summary',0,'correct_offset',0,'counting_field','left_chosen_successful','max_sac_dist',[8,8],'min_sac_amplitude',8,'sac_ini_t',300,'saccade_definition',4,'inferential_on',0};
Sel_all                         = {'display',0,'summary',0,'keep_raw_data',1,'correct_offset',0,'counting_field','left_chosen_successful','max_sac_dist',[8,8],'min_sac_amplitude',8,'sac_ini_t',300,'saccade_definition',4,'inferential_on',0};

global GLO
% %direct
% file={'Y:\Data\Curius_microstim_with_parameters\20140719',8};
% GLO.windows_to_use=[3 4 5 6  7 8 11];
%GLO.type_to_use                     = 2;
% memorz
%file={'Y:\Data\Linus_microstim_with_parameters\20150227',14};
%file={'Y:\Data\Linus_microstim_with_parameters\20150305',13};
file={'Y:\Data\Curius_microstim_with_parameters\20150306',16};
%file={'Y:\Data\Curius_microstim_with_parameters\20150417',13};
GLO.windows_to_use=[17 19 20 22];
GLO.type_to_use                     = 3;

out_comp=monkeypsych_analyze_working(file,Sel_all);
abort_codes={out_comp{1}.task.abort_code}';
target_positions=unique([out_comp{1}.saccades.tar_pos] - [out_comp{1}.saccades.fix_pos] );
target_positions=target_positions(1:6);
batch_title='whatever';
errors=1;

GLO.All_stim_states                 = {'fix','fix','fix','fix','fix','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','fix','cue','cue','mem','tar_acq_inv','tar_acq_inv'};
GLO.All_stim_windows                = [-0.200 -0.160 -0.120 -0.080 -0.040 0 0.040 0.080 0.100 0.110 0.120 0.130 0.140 0.150 0.160 0.240 -0.080 0 0.080 -0.080 0 0.080];


eye_traces(out_comp,target_positions,batch_title,errors);