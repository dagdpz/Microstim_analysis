function [bias_B_trial,bias_S,mean_rt,mean_x_acc,bias_B_run,mean_rt_B_run,mean_x_acc_B_run] = bias_and_rt_microstim(varargin)
% AUDV
% LS
% look at [out_comp variable_to_test counter]=
% monkeypsych_analyze_working({'L:\Data\Linus\20150227',14},{'summary',0,'display',1,'success',0}),
% errors in cue are going to the ipsi side


% close all
batch_title                         = [];
stimulated                          = {};
baselines                           = {};
bias_B_trial                        = {};
bias_S                              = {};
mean_rt                             = {};
mean_x_acc                          = {};
bias_B_run                          = {};
mean_rt_B_run                       = {};
mean_x_acc_B_run                    = {};
saccade_definition                  = 4; %1: closest, 2: biggest, 3:last, 4:first (closest for accuracy, first for all the rest)


global GLO

GLO.specific_windows                = 1;
if nargin>0
    GLO.mon_to_test=varargin{:};
else
    GLO.mon_to_test                     = 'Linus_memory_saccades_dorsal';%'Curius_direct_saccades_dorsal';
end

switch GLO.mon_to_test
    case 'Combined_memory_dorsal',                      monkey='Combined';                  GLO.type_to_use = 3;
    case 'Combined_memory_dorsal_incomplete',           monkey='Combined';                  GLO.type_to_use = 3;
        
    case 'Combined_memory_ventral_shallow',             monkey='Combined';                  GLO.type_to_use = 3;
        %   case 'Combined_memory_ventral_deep',                monkey='Combined';                  GLO.type_to_use = 3;
    case 'Combined_memory_ventral_medium_deep',         monkey='Combined';                  GLO.type_to_use = 3;
    case 'Combined_memory_ventral_very_deep',            monkey='Combined';                 GLO.type_to_use = 3;
        
    case 'Combined_direct_dorsal',                      monkey='Combined';                  GLO.type_to_use = 2;
    case 'Combined_direct_ventral_shallow',             monkey='Combined';                  GLO.type_to_use = 2;
        %  case 'Combined_direct_ventral_deep',                monkey='Combined';                  GLO.type_to_use = 2;
    case 'Combined_direct_ventral_medium_deep',         monkey='Combined';                  GLO.type_to_use = 2;
    case 'Combined_direct_ventral_very_deep',           monkey='Combined';                  GLO.type_to_use = 2;
        
    case 'Combined_direct_dorsal_earlier_windows',      monkey='Combined';                  GLO.type_to_use = 2; GLO.specific_windows=2;
    case 'Combined_direct_outside',                     monkey='Combined';                  GLO.type_to_use = 2;
        
    case 'Combined_direct_dorsal_memory_comparison',    monkey='Combined';                  GLO.type_to_use = 2;
        
        
        
    case 'Linus_direct_saccades_dorsal',                monkey='Linus';                     GLO.type_to_use = 2;
    case 'Linus_direct_saccades_ventral_shallow',       monkey='Linus';                     GLO.type_to_use = 2;
        %   case 'Linus_direct_saccades_ventral_deep',          monkey='Linus';                     GLO.type_to_use = 2;
    case 'Linus_direct_saccades_ventral_medium_deep',   monkey='Linus';                     GLO.type_to_use = 2;
    case 'Linus_direct_saccades_ventral_very_deep',     monkey='Linus';                     GLO.type_to_use = 2;
    case 'Linus_direct_outside',                        monkey='Linus';                     GLO.type_to_use = 2;
        
    case 'Linus_memory_saccades_dorsal',                monkey='Linus';                     GLO.type_to_use = 3;
    case 'Linus_memory_saccades_dorsal_incomplete',                monkey='Linus';                     GLO.type_to_use = 3;
    case 'Linus_memory_saccades_ventral_shallow',       monkey='Linus';                     GLO.type_to_use = 3;
    case 'Linus_memory_saccades_ventral_deep',          monkey='Linus';                     GLO.type_to_use = 3;
    case 'Linus_memory_saccades_ventral_very_deep',     monkey='Linus';                     GLO.type_to_use = 3;
    case 'Linus_memory_saccades_ventral_medium_deep',   monkey='Linus';                     GLO.type_to_use = 3;
        %%!!!!!!!!!
        
    case 'Curius_direct_saccades_dorsal',               monkey='Curius';                    GLO.type_to_use = 2;
    case 'Curius_direct_saccades_ventral_shallow',      monkey='Curius';                    GLO.type_to_use = 2;
    case 'Curius_direct_saccades_ventral_medium_deep',  monkey='Curius';                    GLO.type_to_use = 2;
    case 'Curius_direct_saccades_ventral_very_deep',    monkey='Curius';                    GLO.type_to_use = 2;
        %   case 'Curius_direct_saccades_ventral_deep',         monkey='Curius';                    GLO.type_to_use = 2;
    case 'Curius_direct_outside',                       monkey='Curius';                     GLO.type_to_use = 2;
        
    case 'Curius_memory_saccades_dorsal',               monkey='Curius';                    GLO.type_to_use = 3;
    case 'Curius_memory_saccades_dorsal_incomplete',    monkey='Curius';                    GLO.type_to_use = 3;
    case 'Curius_inc',                                  monkey='Curius';                    GLO.type_to_use = 3;
    case 'Curius_dir_for_inc',                          monkey='Curius';                    GLO.type_to_use = 2;
    case 'Curius_memory_saccades_ventral_shallow',      monkey='Curius';                    GLO.type_to_use = 3;
    case 'Curius_memory_saccades_ventral_medium_deep',  monkey='Curius';                    GLO.type_to_use = 3;
    case 'Curius_memory_saccades_ventral_very_deep',    monkey='Curius';                    GLO.type_to_use = 3;
        
        
    
    case 'Shorter_duration_100ms',                      monkey='Curius';                    GLO.type_to_use = 2;
    case 'Shorter_duration_150ms',                      monkey='Curius';                    GLO.type_to_use = 2;
    case 'Shorter_duration_control',                    monkey='Curius';                    GLO.type_to_use = 2;    
    case 'Test2',                                       monkey='Linus';                    GLO.type_to_use = 2;
    case 'TestCornelius',                               monkey='Cornelius';                    GLO.type_to_use = 2;
        %   case 'Curius_memory_saccades_ventral_deep',         monkey='Curius';                    GLO.type_to_use = 3;
        
        %     case 'Curius',                                      monkey='Curius';                    GLO.type_to_use = 2;
        %     case 'Test',                                        monkey='Test';                      GLO.type_to_use = 3;
        %     case 'Test2',                                       monkey='Test2';                     GLO.type_to_use = 3;
        %     case 'Test3',                                       monkey='Test3';                     GLO.type_to_use = 3;
        %     case 'Test4',                                       monkey='Test4';                     GLO.type_to_use = 2;
end




user_drive=getUserName;
GLO.drive                           = 'Y';
switch user_drive
    case 'a.doming'
        GLO.drive                           = 'K';
    case 'lschneider'
        GLO.drive                           = 'W';
    case 'dagadmin'
        GLO.drive                           = 'Y';
end

GLO.start_date                      = 20130101; %20140101
GLO.end_date                        = 20160101;
GLO.effector_to_use                 = 0;
GLO.reach_hand                      = 2;
GLO.run_analyze                     = 1;
GLO.save_analyze_output             = 1;
GLO.run_internal_calculation        = 1;
GLO.save_batch_for_later            = 1;
GLO.create_pdf                      = 1;
GLO.append_pdfs                     = 0;
GLO.plot_batches_individually       = 0;
GLO.run_baseline_included_in_ttests = 0;
GLO.stat_to_use                     = 'signed rank';%'signed rank';%'ttest_bonf'; % ttest paired, ttest unpaired, signed rank 'signed rank'; %
GLO.ttest_text                      = 1;
GLO.RT_vs_bias_condition            = 'Instructed ipsi'; %'Choice contra','Choice ipsi','Instructed contra','Instructed ipsi'
GLO.low_excentricity_threshold      = 15;
GLO.Labels                          = {'RB','B','-200', '-160', '-120','-80','-40','Go','40','80','100','110','120','130','140','150','160','240','-80 Cue', 'Cue', '80 Cue','-80 Go', 'Go','80 Go'};
GLO.All_stim_states                 = {'fix','fix','fix','fix','fix','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','fix','cue','cue','mem','tar_acq_inv','tar_acq_inv'};
GLO.All_stim_windows                = [-0.200 -0.160 -0.120 -0.080 -0.040 0 0.040 0.080 0.100 0.110 0.120 0.130 0.140 0.150 0.160 0.240 -0.080 0 0.080 -0.080 0 0.080];
GLO.windows_to_use                  = 1:numel(GLO.All_stim_windows);
GLO.folder_to_save                  = strcat(GLO.drive,':', filesep, 'Projects', filesep, 'Pulv_microstim_behavior', filesep, 'behavior', filesep, monkey, '_summaries', filesep);
GLO.tests                           = {'Bias','farther'}; % Bias || Endpoints %  all || farther || closer
GLO.excentricity                    = GLO.tests{1,2};

% for summaries %,
 GLO.Sel_all                         = {'display',0,'summary',0,'correct_offset',1,'counting_field','left_chosen_successful','type',GLO.type_to_use,'max_sac_dist',[8,8],'min_sac_amplitude',8,'sac_ini_t',300,'saccade_definition',saccade_definition,'inferential_on',0};
GLO.plot_to_show                    = [1,2,3]; %[2];%[1:7,9:17]; 
GLO.table                           = 1;
GLO.return                          = 0;
GLO.target_pos_precision            = 2;


GLO.plot_eye_traces                 = 0;
errors_only_for_eye_trace_plot      = 1;
if GLO.plot_eye_traces
GLO.Sel_all                         = {'display',0,'summary',0,'saccade_2bo',1,'correct_offset',0,'counting_field','left_chosen_successful','type',GLO.type_to_use,'max_sac_dist',[8,8],'min_sac_amplitude',8,'sac_ini_t',300,'saccade_definition',saccade_definition,'inferential_on',0};
    GLO.Sel_all = [GLO.Sel_all {'keep_raw_data',1}];
end

% 1 - 'Bias RT'
% 2 - 'Significance across sessions'
% 3 - 'bias vs rt all'
% 4 - 'RT Histograms'
% 5 - 'Compiled rt error bars'
% 6 - 'Velocity histograms'
% 7 - 'Precision histograms'
% 8 - 'correlations'
% 9 - 'positions'
% 10 - 'scatterplot rts'
% 11 - 'choice over time per session'
% 12 - 'choice over time, averaged'
% 13 - 'Hitrates'
% 14 - 'bias and RT in two modes'
% 15 - 'RT condition differences'
% 16 - 'Bias vs RT separated by windows'
% 17 - 'Accuracy Euclidean'
% 18 - 'Accuracy per Position'
% 19 - 'Accuracy 2D'

% -1  - 'all'

if GLO.specific_windows
    switch GLO.type_to_use
        case 2
            switch GLO.effector_to_use
%                 case 0,                     GLO.windows_to_use=[]; % for saccade dataset from initial window
%                 case 4,                     GLO.windows_to_use=[]; % for reaches dataset from initial window
                case 0,                     GLO.windows_to_use=[3 4 5 6  7 8 11]; % for saccade dataset from initial window
                case 4,                     GLO.windows_to_use=[4 6 8 16]; % for reaches dataset from initial window
            end
        case 3
            switch GLO.effector_to_use
                %                 case 0,                     GLO.windows_to_use=[17 18 19 20 21 22]; % for saccade dataset from initial window
                %                 case 4,                     GLO.windows_to_use=[4 6 8 16]; % for reaches dataset from initial window
                case 0,                     GLO.windows_to_use=[17 19 20 22]; % for saccade dataset from initial window
                case 4,                     GLO.windows_to_use=[4 6 8 16]; % for reaches dataset from initial window
                    
            end
    end
    if GLO.specific_windows==2
        GLO.windows_to_use=[1 2 3 4 5 6 7 8 11 15]; % for saccade dataset from initial window
        
    end
    
    GLO.Labels=GLO.Labels([1 2 GLO.windows_to_use+2]);
end

if  any(strcmp('saccade_definition',GLO.Sel_all))
    sac_def=GLO.Sel_all{find(strcmp('saccade_definition',GLO.Sel_all))+1};
    switch sac_def
        case 1,                     GLO.saccade_type='closest and big enough';
        case 2,                     GLO.saccade_type='biggest and close enough';
        case 3,                     GLO.saccade_type='last saccade in the state';
        case 4,                     GLO.saccade_type='first saccade in the state';
    end
end

GLO.idx_GO = find(GLO.All_stim_windows(GLO.windows_to_use)==0);
GLO.idx_early = 1:GLO.idx_GO-1; % -1 to take all before GO, -2 to not take the two baselines
GLO.idx_late = GLO.idx_GO:numel(GLO.windows_to_use);




switch GLO.mon_to_test
    case 'Linus_direct_saccades_dorsal'
        monkey='Linus';
        batch_title='dorsal pulvinar';
        Aiminput={'use_dir_sti_at_win'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',0;'gridhole_y',3;'task_effector',GLO.effector_to_use};
                Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
            case 4,                     Inputsequal={'gridhole_x',0;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
                Inputslist{1}={'microstim_start',  -80,   0,  80, 240};
        end
        Inputsrange={['hits_t_2_e_' num2str(GLO.effector_to_use)],50,3000;'Current_strength',35,300;'Electrode_depth',46,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',-2;'gridhole_y',4;'task_effector',GLO.effector_to_use};
                Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
            case 4,                     Inputsequal={'gridhole_x',-2;'gridhole_y',4;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
                Inputslist{1}={'microstim_start',  -80,   0,  80, 240};
        end
        Inputsrange={['hits_t_2_e_' num2str(GLO.effector_to_use)],50,3000;'Current_strength',35,300;'Electrode_depth',45,51;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',-1;'gridhole_y',3;'task_effector',GLO.effector_to_use};
                Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
            case 4,                     Inputsequal={'gridhole_x',-1;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
                Inputslist{1}={'microstim_start',  -80,   0,  80, 240};
        end
        Inputsrange={['hits_t_2_e_' num2str(GLO.effector_to_use)],50,3000;'Current_strength',35,300;'Electrode_depth',45,50;'Session',GLO.start_date,GLO.end_date};
        filelist_hole3 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2;filelist_hole3];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        Aiminput={'use_dir_bas'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits_t_2_e_' num2str(GLO.effector_to_use)],50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
    case 'Linus_direct_saccades_ventral_shallow'
        monkey='Linus';
        batch_title='shallow ventral pulvinar';
        Aiminput={'use_dir_sti_at_win'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
                Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
            case 4,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
                Inputslist{1}={'microstim_start',  -80,   0,  80, 240};
        end
        Inputsrange={['hits_t_2_e_' num2str(GLO.effector_to_use)],50,3000;'Current_strength',35,300;'Electrode_depth',48,50.5;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use};
                Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
            case 4,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
                Inputslist{1}={'microstim_start',  -80,   0,  80, 240};
        end
        Inputsrange={['hits_t_2_e_' num2str(GLO.effector_to_use)],50,3000;'Current_strength',35,300;'Electrode_depth',48,50.5;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        Aiminput={'use_dir_bas'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits_t_2_e_' num2str(GLO.effector_to_use)],50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
    case 'Linus_direct_saccades_ventral_deep'
        monkey='Linus';
        batch_title='deep ventral pulvinar';
        Aiminput={'use_dir_sti_at_win'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
                Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
            case 4,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
                Inputslist{1}={'microstim_start',  -80,   0,  80, 240};
        end
        Inputsrange={['hits_t_2_e_' num2str(GLO.effector_to_use)],50,3000;'Current_strength',35,300;'Electrode_depth',50.75,55;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use};
                Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
            case 4,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
                Inputslist{1}={'microstim_start',  -80,   0,  80, 240};
        end
        Inputsrange={['hits_t_2_e_' num2str(GLO.effector_to_use)],50,3000;'Current_strength',35,300;'Electrode_depth',50.75,55;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        Aiminput={'use_dir_bas'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits_t_2_e_' num2str(GLO.effector_to_use)],50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
    case 'Linus_direct_saccades_ventral_medium_deep'
        monkey='Linus';
        batch_title='medium deep ventral pulvinar';
        Aiminput={'use_dir_sti_at_win'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
                Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
            case 4,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
                Inputslist{1}={'microstim_start',  -80,   0,  80, 240};
        end
        Inputsrange={['hits_t_2_e_' num2str(GLO.effector_to_use)],50,3000;'Current_strength',35,300;'Electrode_depth',50.75,51;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use};
                Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
            case 4,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
                Inputslist{1}={'microstim_start',  -80,   0,  80, 240};
        end
        Inputsrange={['hits_t_2_e_' num2str(GLO.effector_to_use)],50,3000;'Current_strength',35,300;'Electrode_depth',50.75,51;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        Aiminput={'use_dir_bas'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits_t_2_e_' num2str(GLO.effector_to_use)],50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
    case 'Linus_direct_saccades_ventral_very_deep'
        monkey='Linus';
        batch_title='very deep ventral pulvinar';
        Aiminput={'use_dir_sti_at_win'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
                Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
            case 4,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
                Inputslist{1}={'microstim_start',  -80,   0,  80, 240};
        end
        Inputsrange={['hits_t_2_e_' num2str(GLO.effector_to_use)],50,3000;'Current_strength',35,300;'Electrode_depth',51.25,55;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use};
                Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
            case 4,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
                Inputslist{1}={'microstim_start',  -80,   0,  80, 240};
        end
        Inputsrange={['hits_t_2_e_' num2str(GLO.effector_to_use)],50,3000;'Current_strength',35,300;'Electrode_depth',51.25,55;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        Aiminput={'use_dir_bas'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits_t_2_e_' num2str(GLO.effector_to_use)],50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
    case 'Linus_memory_saccades_dorsal'
        monkey='Linus';
        batch_title='dorsal pulvinar';
        Aiminput={'use_mem_sti_at_win'};
        Inputslist={};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',-1;'gridhole_y',3;'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'gridhole_x',-1;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits'],50,3000;'Current_strength',35,300;'Electrode_depth',46,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        filelist_stimulated=[filelist_hole1];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        Aiminput={'use_mem_bas'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits' ],50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
            case 'Linus_memory_saccades_dorsal_incomplete'
        monkey='Linus';
        batch_title='dorsal pulvinar';
        Aiminput={'inc_mem_sti_at_win'};
        Inputslist={};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',-1;'gridhole_y',3;'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'gridhole_x',-1;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits'],50,3000;'Current_strength',35,300;'Electrode_depth',46,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        filelist_stimulated=[filelist_hole1];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        Aiminput={'use_mem_bas'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits' ],50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
    case 'Linus_memory_saccades_ventral_shallow'
        monkey='Linus';
        batch_title='shallow ventral pulvinar';
        Aiminput={'use_mem_sti_at_win'};
        Inputslist={};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',48,50.5;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',48,50.5;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        Aiminput={'use_mem_bas'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits' ],50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
    case 'Linus_memory_saccades_ventral_deep'
        monkey='Linus';
        batch_title='deep ventral pulvinar';
        Aiminput={'use_mem_sti_at_win'};
        Inputslist={};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',50.75,55;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',50.75,55;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        Aiminput={'use_mem_bas'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits' ],50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
        
    case 'Linus_memory_saccades_ventral_very_deep'
        monkey='Linus';
        batch_title='very deep ventral pulvinar';
        Aiminput={'use_mem_sti_at_win'};
        Inputslist={};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',51.25,55;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',51.25,55;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        Aiminput={'use_mem_bas'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits' ],50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
    case 'Linus_memory_saccades_ventral_medium_deep'
        monkey='Linus';
        batch_title='medium deep ventral pulvinar';
        Aiminput={'use_mem_sti_at_win'};
        Inputslist={};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',50.75,51;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',50.75,51;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        Aiminput={'use_mem_bas'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits' ],50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
    case 'Curius_direct_saccades_dorsal'
        %GLO.end_date                        = 20150210;
        monkey='Curius';
        batch_title='dorsal pulvinar';
        Aiminput={'use_dir_sti_at_win'};
        GLO.effector_to_use=0;
        
        Inputsequal={'gridhole_x',4;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60;'Session',GLO.start_date,GLO.end_date};
        %                   Inputslist{1}={'x_distance_to_center', -23-8i, 23-8i, 23+8i,-23+8i, 24+0i,-24+0i};
        %                   Inputslist={};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        
        Inputsequal={'gridhole_x',3;'gridhole_y',5;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole3 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2;filelist_hole3];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        Aiminput={'use_dir_bas'};
        
        
        Inputsequal={'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        %                   Inputslist{1}={'x_distance_to_center', -23-8i, 23-8i, 23+8i,-23+8i, 24+0i,-24+0i};
        %                   Inputslist={'microstim_start',''};
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
        
        %         parameterfile=[GLO.drive,':', filesep,'microstim_behavior', filesep,  monkey, '_summaries', filesep, monkey, '_updated_parameters.xls'];
        %         reduce_stimulation_parameters([filelist_stimulated;filelist_baselines],parameterfile)
        
    case 'Curius_memory_saccades_dorsal'
        monkey='Curius';
        batch_title='dorsal pulvinar';
        Aiminput={'use_mem_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
        Inputsequal     ={'gridhole_x',4;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',3;'gridhole_y',5;'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole3 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2,filelist_hole3];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_mem_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
        
        %         parameterfile=[GLO.drive,':', filesep,'microstim_behavior', filesep,  monkey, '_summaries', filesep, monkey, '_updated_parameters.xls'];
        %         reduce_stimulation_parameters([filelist_stimulated;filelist_baselines],parameterfile)
        
            case 'Curius_memory_saccades_dorsal_incomplete'
        monkey='Curius';
        batch_title='dorsal pulvinar';
        Aiminput={'inc_mem_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
        Inputsequal     ={'gridhole_x',4;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',3;'gridhole_y',5;'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole3 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2,filelist_hole3];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_mem_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
        
        %         parameterfile=[GLO.drive,':', filesep,'microstim_behavior', filesep,  monkey, '_summaries', filesep, monkey, '_updated_parameters.xls'];
        %         reduce_stimulation_parameters([filelist_stimulated;filelist_baselines],parameterfile)
        
        
        
    case 'Curius_direct_saccades_ventral_shallow'
        batch_title='shallow ventral pulvinar';
        Aiminput={'use_dir_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
        Inputsequal     ={'gridhole_x',8;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',44,47.4;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal     ={'gridhole_x',9;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',44,47.4;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_dir_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
    case 'Curius_memory_saccades_ventral_shallow'
        monkey='Curius';
        batch_title='shallow ventral pulvinar';
        Aiminput={'use_mem_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
        Inputsequal     ={'gridhole_x',8;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',44,47.4;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal     ={'gridhole_x',9;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',44,47.4;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_mem_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
        
    case 'Curius_direct_saccades_ventral_medium_deep'
        monkey='Curius';
        batch_title='medium deep ventral pulvinar';
        Aiminput={'use_dir_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
        Inputsequal     ={'gridhole_x',8;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',47.4,47.9;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal     ={'gridhole_x',9;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',47.4,47.9;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_dir_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
    case 'Curius_memory_saccades_ventral_medium_deep'
        monkey='Curius';
        batch_title='medium deep ventral pulvinar';
        Aiminput={'use_mem_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
        Inputsequal     ={'gridhole_x',8;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',47.4,47.9;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal     ={'gridhole_x',9;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',47.4,47.9;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_mem_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
    case 'Curius_direct_saccades_ventral_very_deep'
        monkey='Curius';
        batch_title='very deep ventral pulvinar';
        Aiminput={'use_dir_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
        Inputsequal     ={'gridhole_x',8;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',48,49;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal     ={'gridhole_x',9;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',48,49;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_dir_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
        
        
    case 'Curius_memory_saccades_ventral_very_deep'
        monkey='Curius';
        batch_title='very deep ventral pulvinar';
        Aiminput={'use_mem_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
        Inputsequal     ={'gridhole_x',8;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',48,49;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal     ={'gridhole_x',9;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',48,49;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_mem_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
        
    case 'Combined_memory_dorsal'
        % Curius
        monkey='Curius';
        batch_title='dorsal pulvinar';
        Aiminput={'use_mem_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
        Inputsequal     ={'gridhole_x',4;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        
        Inputsequal={'gridhole_x',3;'gridhole_y',5;'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole3 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2,filelist_hole3];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_C = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_mem_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_C = DAG_arrange_baselines(baselines_in,stimulated_C);
        
        
        monkey='Linus';  Aiminput={'use_mem_sti_at_win'}; Inputslist={};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',-1;'gridhole_y',3;'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'gridhole_x',-1;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits'],50,3000;'Current_strength',35,300;'Electrode_depth',46,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        filelist_stimulated=[filelist_hole1];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_L = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        Aiminput={'use_mem_bas'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits' ],50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_L = DAG_arrange_baselines(baselines_in,stimulated_L);
        
        stimulated= [stimulated_C stimulated_L];
        baselines= [baselines_C baselines_L];
        monkey='Combined';
        
        
        
        case 'Combined_memory_dorsal_incomplete'
        % Curius
        monkey='Curius';
        batch_title='dorsal pulvinar';
        Aiminput={'inc_mem_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
        Inputsequal     ={'gridhole_x',4;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        
        Inputsequal={'gridhole_x',3;'gridhole_y',5;'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole3 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2,filelist_hole3];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_C = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'mem_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_C = DAG_arrange_baselines(baselines_in,stimulated_C);
        
        
        monkey='Linus';  Aiminput={'inc_mem_sti_at_win'}; Inputslist={};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',-1;'gridhole_y',3;'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'gridhole_x',-1;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits'],50,3000;'Current_strength',35,300;'Electrode_depth',46,60;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        filelist_stimulated=[filelist_hole1];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_L = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        Aiminput={'mem_bas'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits' ],50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_L = DAG_arrange_baselines(baselines_in,stimulated_L);
        
        stimulated= [stimulated_C stimulated_L];
        baselines= [baselines_C baselines_L];
        monkey='Combined';
        
        
        
        
        
    case 'Combined_memory_ventral_shallow'
        monkey='Curius';
        batch_title='shallow ventral pulvinar';
        Aiminput={'use_mem_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
        Inputsequal     ={'gridhole_x',8;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',44,47.4;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal     ={'gridhole_x',9;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',44,47.4;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_C = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_mem_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_C = DAG_arrange_baselines(baselines_in,stimulated_C);
        
        
        monkey='Linus';
        Aiminput={'use_mem_sti_at_win'};
        Inputslist={};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',48,50.5;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',48,50.5;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_L = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        Aiminput={'use_mem_bas'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_L = DAG_arrange_baselines(baselines_in,stimulated_L);
        
        stimulated= [stimulated_C stimulated_L];
        baselines= [baselines_C baselines_L];
        monkey='Combined';
        
        
    case 'Combined_memory_ventral_medium_deep'
        monkey='Curius';
        batch_title='medium deep ventral pulvinar';
        Aiminput={'use_mem_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
        Inputsequal     ={'gridhole_x',8;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',47.4,47.9;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal     ={'gridhole_x',9;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',47.4,47.9;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_C = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_mem_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_C = DAG_arrange_baselines(baselines_in,stimulated_C);
        
        monkey='Linus';
        Aiminput={'use_mem_sti_at_win'};
        Inputslist={};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',50.75,51;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',50.75,51;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_L = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        Aiminput={'use_mem_bas'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits' ],50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_L = DAG_arrange_baselines(baselines_in,stimulated_L);
        
        stimulated= [stimulated_C stimulated_L];
        baselines= [baselines_C baselines_L];
        monkey='Combined';
        
    case 'Combined_direct_ventral_medium_deep'
        
        monkey='Curius';
        batch_title='medium deep ventral pulvinar';
        Aiminput={'use_dir_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
        Inputsequal     ={'gridhole_x',8;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',47.4,47.9;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal     ={'gridhole_x',9;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',47.4,47.9;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_C = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_dir_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_C = DAG_arrange_baselines(baselines_in,stimulated_C);
        
        monkey='Linus';
        Aiminput={'use_dir_sti_at_win'};
        Inputslist={};
        Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',50.75,51;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',50.75,51;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_L = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_dir_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_L = DAG_arrange_baselines(baselines_in,stimulated_L);
        
        stimulated= [stimulated_C stimulated_L];
        baselines= [baselines_C baselines_L];
        monkey='Combined';
        
    case 'Combined_memory_ventral_very_deep'
        monkey='Curius';
        batch_title='very deep ventral pulvinar';
        Aiminput={'use_mem_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
        Inputsequal     ={'gridhole_x',8;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',48,49;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal     ={'gridhole_x',9;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',48,49;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_C = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_mem_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_C = DAG_arrange_baselines(baselines_in,stimulated_C);
        
        monkey='Linus';
        Aiminput={'use_mem_sti_at_win'};
        Inputslist={};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',51.25,55;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',51.25,55;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_L = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        Aiminput={'use_mem_bas'};
        switch GLO.effector_to_use
            case 0,                     Inputsequal={'task_effector',GLO.effector_to_use};
            case 4,                     Inputsequal={'task_effector',GLO.effector_to_use;'reach_hand',GLO.reach_hand};
        end
        Inputsrange={['hits' ],50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_L = DAG_arrange_baselines(baselines_in,stimulated_L);
        
        stimulated= [stimulated_C stimulated_L];
        baselines= [baselines_C baselines_L];
        monkey='Combined';
        
        
    case 'Combined_direct_dorsal'
        GLO.effector_to_use=0;
        
        batch_title='dorsal pulvinar';
        monkey='Curius';
        Aiminput={'use_dir_sti_at_win'};
        
        Inputsequal={'gridhole_x',4;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',3;'gridhole_y',5;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60};
        filelist_hole3 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2;filelist_hole3];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_C = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_dir_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_C = DAG_arrange_baselines(baselines_in,stimulated_C);
        
        
        monkey='Linus'; Aiminput={'use_dir_sti_at_win'};
        Inputsequal={'gridhole_x',0;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',46,60};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',-2;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',45,50};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',-1;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',45,51};
        filelist_hole3 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2;filelist_hole3];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_L = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_dir_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_L = DAG_arrange_baselines(baselines_in,stimulated_L);
        
        stimulated= [stimulated_C stimulated_L];
        baselines= [baselines_C baselines_L];
        monkey='Combined';
        
    case 'Combined_direct_ventral_shallow'
        monkey='Curius';
        batch_title='shallow ventral pulvinar';
        Aiminput={'use_dir_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
        Inputsequal     ={'gridhole_x',8;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',44,47.4;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal     ={'gridhole_x',9;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',44,47.4;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_C = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_dir_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_C = DAG_arrange_baselines(baselines_in,stimulated_C);
        
        monkey='Linus';
        Aiminput={'use_dir_sti_at_win'};
        Inputslist={};
        Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',48,50.5;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',48,50.5;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_L = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_dir_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_L = DAG_arrange_baselines(baselines_in,stimulated_L);
        
        stimulated= [stimulated_C stimulated_L];
        baselines= [baselines_C baselines_L];
        monkey='Combined';
        
    case 'Combined_direct_ventral_very_deep'
        monkey='Curius';
        batch_title='very deep ventral pulvinar';
        Aiminput={'use_dir_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
        Inputsequal     ={'gridhole_x',8;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',48,49;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal     ={'gridhole_x',9;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',48,49;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_C = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_dir_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        
        Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_C = DAG_arrange_baselines(baselines_in,stimulated_C);
        
        monkey='Linus';
        Aiminput={'use_dir_sti_at_win'};
        
        Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',51.25,55;'Session',GLO.start_date,GLO.end_date};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',51.25,55;'Session',GLO.start_date,GLO.end_date};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_L = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_dir_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        Inputsrange={['hits' ],50,3000;'Session',GLO.start_date,GLO.end_date};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_L = DAG_arrange_baselines(baselines_in,stimulated_L);
        
        stimulated= [stimulated_C stimulated_L];
        baselines= [baselines_C baselines_L];
        monkey='Combined';
        
    case 'Combined_direct_outside'
        GLO.effector_to_use=0;
        
        batch_title='outside pulvinar';
        monkey='Curius';
        Aiminput={'use_dir_sti_at_win'};
        
        clear Inputslist
        Inputslist={};
        Inputsequal={'gridhole_x',4;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        %Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43,43.6};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        %Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43,43.6};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        
        Inputsequal={'gridhole_x',3;'gridhole_y',5;'task_effector',GLO.effector_to_use};
        %Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',40,43.6};
        filelist_hole3 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2;filelist_hole3];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_C = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_dir_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_C = DAG_arrange_baselines(baselines_in,stimulated_C);
        
        monkey='Linus';
        Aiminput={'use_dir_sti_at_win'};
        Inputsequal={'gridhole_x',0;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        %Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,350;'Electrode_depth',40,45.9};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',-2;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        %Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,350;'Electrode_depth',40,44.9};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',-1;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        %Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,350;'Electrode_depth',40,44.9};
        filelist_hole3 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2;filelist_hole3];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_L = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_dir_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_L = DAG_arrange_baselines(baselines_in,stimulated_L);
        
        
        stimulated= [stimulated_C stimulated_L];
        baselines= [baselines_C baselines_L];
        monkey='Combined';
        
    case 'Curius_direct_outside'
        GLO.effector_to_use=0;
        
        batch_title='outside pulvinar';
        monkey='Curius';
        Aiminput={'use_dir_sti_at_win'};
        
        clear Inputslist
        Inputslist={};
        Inputsequal={'gridhole_x',4;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        %Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43,43.6};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        %Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43,43.6};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        
        Inputsequal={'gridhole_x',3;'gridhole_y',5;'task_effector',GLO.effector_to_use};
        %Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',40,43.6};
        filelist_hole3 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2;filelist_hole3];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_dir_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
    case 'Linus_direct_outside'
        batch_title='outside pulvinar';
        monkey='Linus';
        Aiminput={'use_dir_sti_at_win'};
        clear Inputslist
        Inputslist={};
        Inputsequal={'gridhole_x',0;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        %Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,350;'Electrode_depth',40,45.9};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',-2;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        %Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,350;'Electrode_depth',40,44.9};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        Inputsequal={'gridhole_x',-1;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        %Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,350;'Electrode_depth',40,44.9};
        filelist_hole3 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2;filelist_hole3];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_dir_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
    case 'Combined_direct_dorsal_earlier_windows'
        GLO.effector_to_use=0;
        
        batch_title='dorsal pulvinar earlier windows';
        monkey='Curius';
        Aiminput={'use_dir_sti_at_win'};
        GLO.effector_to_use=0;
        
        Inputsequal={'gridhole_x',4;'gridhole_y',4;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -200, -160, -120, -80, -40,   0,  40,  80, 120, 160};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',43.6,60};
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        % same hole, slightly different windows....
        Inputslist{1}={'microstim_start', -200, -160, -120, -80, -40,   0};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_C = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_dir_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_C = DAG_arrange_baselines(baselines_in,stimulated_C);
        
        monkey='Linus';
        Aiminput={'use_dir_sti_at_win'};
        Inputsequal={'gridhole_x',0;'gridhole_y',3;'task_effector',GLO.effector_to_use};
        Inputslist{1}={'microstim_start', -200, -160, -120, -80, -40,   0,  40,  80, 120};
        Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',45,60}; %% this is actually not deep enough
        filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        % same hole, slightly different windows....
        Inputslist{1}={'microstim_start', -200, -160, -120, -80, -40,   0,  40,  80, 120, 160};
        filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        
        filelist_stimulated=[filelist_hole1;filelist_hole2];
        Inputsequal_for_batch={'Session','Electrode_depth'};
        stimulated_L = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
        
        Aiminput={'use_dir_bas'};
        Inputsequal={'task_effector',GLO.effector_to_use};
        Inputsrange={'hits',50,3000};
        clear Inputslist
        Inputslist={};
        filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
        Inputsequal={'Session'};
        baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
        baselines_L = DAG_arrange_baselines(baselines_in,stimulated_L);
        
        stimulated= [stimulated_C stimulated_L];
        baselines= [baselines_C baselines_L];
        monkey='Combined';
        
    case 'Test2'
%         monkey='Curius';
%         batch_title='shallow ventral pulvinar';
%         Aiminput={'use_dir_sti_at_win'};
%         GLO.effector_to_use=0;
%         Inputslist={};
%         
%      stimulated{1} = [{[GLO.drive ':\Data\Curius_microstim\20150903']},{11}];
%      baselines{1} = [{[GLO.drive ':\Data\Curius_microstim\20150903']},{7};{[GLO.drive ':\Data\Curius_microstim\20150903']},{8}]; %% actually NO baseline
        
     
       monkey='Linus';
        batch_title='20151001';
        Aiminput={'use_mem_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
%      stimulated{1} = [{[GLO.drive ':\Data\Linus_microstim\20150930']},{9}];
%      baselines{1} = [{[GLO.drive ':\Data\Linus_microstim\20150930']},{8}]; %% actually NO baseline

%      stimulated{1} = [{[GLO.drive ':\Data\Linus_microstim\20151001']},{11}];
%      baselines{1} = [{[GLO.drive ':\Data\Linus_microstim\20151001']},{7}]; %% actually NO baseline

     stimulated{1} = [{[GLO.drive ':\Data\Linus_microstim\20151002']},{9};{[GLO.drive ':\Data\Linus_microstim\20151002']},{10}];
     baselines{1} = [{[GLO.drive ':\Data\Linus_microstim\20151002']},{7}]; %% actually NO baseline
        
%         Inputsequal     ={'gridhole_x',8;'gridhole_y',4;'task_effector',GLO.effector_to_use};
%         Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
%         Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',44,47.4;'Session',GLO.start_date,GLO.end_date};
%         filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
%         
%         Inputsequal     ={'gridhole_x',9;'gridhole_y',4;'task_effector',GLO.effector_to_use};
%         Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
%         Inputsrange     ={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',44,47.4;'Session',GLO.start_date,GLO.end_date};
%         filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
%         
%         filelist_stimulated=[filelist_hole1;filelist_hole2];
%         Inputsequal_for_batch={'Session','Electrode_depth'};
%         stimulated_C = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
%         
%         Aiminput={'use_dir_bas'};
%         Inputsequal={'task_effector',GLO.effector_to_use};
%         
%         Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
%         clear Inputslist
%         Inputslist={};
%         filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
%         Inputsequal={'Session'};
%         baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
%         baselines_C = DAG_arrange_baselines(baselines_in,stimulated_C);
%         
%         monkey='Linus';
%         Aiminput={'use_dir_sti_at_win'};
%         Inputslist={};
%         Inputsequal={'gridhole_x',5;'gridhole_y',3;'task_effector',GLO.effector_to_use};
%         Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
%         Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',48,50.5;'Session',GLO.start_date,GLO.end_date};
%         filelist_hole1 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
%         
%         Inputsequal={'gridhole_x',6;'gridhole_y',3;'task_effector',GLO.effector_to_use};
%         Inputslist{1}={'microstim_start', -120, -80, -40,   0,  40,  80, 120};
%         Inputsrange={'hits',50,3000;'Current_strength',35,300;'Electrode_depth',48,50.5;'Session',GLO.start_date,GLO.end_date};
%         filelist_hole2 = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
%         
%         filelist_stimulated=[filelist_hole1;filelist_hole2];
%         Inputsequal_for_batch={'Session','Electrode_depth'};
%         stimulated_L = DAG_get_batch_input_from_xls(filelist_stimulated,Inputsequal_for_batch,0,GLO.drive,monkey);
%         
%         Aiminput={'use_dir_bas'};
%         Inputsequal={'task_effector',GLO.effector_to_use};
%         Inputsrange={'hits',50,3000;'Session',GLO.start_date,GLO.end_date};
%         clear Inputslist
%         Inputslist={};
%         filelist_baselines = get_filelist_from_xls(Aiminput,Inputsequal,Inputsrange,Inputslist,GLO.drive,monkey);
%         Inputsequal={'Session'};
%         baselines_in = DAG_get_batch_input_from_xls(filelist_baselines,Inputsequal,0,GLO.drive,monkey);
%         baselines_L = DAG_arrange_baselines(baselines_in,stimulated_L);
%         
%         stimulated= [stimulated_C stimulated_L];
%         baselines= [baselines_C baselines_L];
%         monkey='Combined';
        
    case 'Combined_direct_dorsal_memory_comparison'    
                batch_title='memory control';
                stimulated{1} = [{[GLO.drive ':\Data\Curius_microstim_with_parameters\20150306']},{10}];
                baselines{1} = [{[GLO.drive ':\Data\Curius_microstim_with_parameters\20150306']},{8};{[GLO.drive ':\Data\Curius_microstim_with_parameters\20150306']},{9}]; %% actually NO baseline
        
                stimulated{2} = [{[GLO.drive ':\Data\Curius_microstim_with_parameters\20150417']},{11};{[GLO.drive ':\Data\Curius_microstim_with_parameters\20150417']},{15}];
                baselines{2} = [{[GLO.drive ':\Data\Curius_microstim_with_parameters\20150417']},{10}];
                
                stimulated{3} = [{[GLO.drive ':\Data\Curius_microstim_with_parameters\20150410']},{12}]; %% New one for manuscript, before skipped because memory part was incomplete 240 trials
                baselines{3} = [{[GLO.drive ':\Data\Curius_microstim_with_parameters\20150410']},{11}];
                
                stimulated{4} = [{[GLO.drive ':\Data\Linus_microstim_with_parameters\20150227']},{12}];
                baselines{4} = [{[GLO.drive ':\Data\Linus_microstim_with_parameters\20150227']},{16}];
                
                
                stimulated{5} = [{[GLO.drive ':\Data\Linus_microstim_with_parameters\20150305']},{10}];
                baselines{5} = [{[GLO.drive ':\Data\Linus_microstim_with_parameters\20150305']},{9}];
                

                
%                 stimulated_C=stimulated(1:2);
%                 stimulated_L=stimulated(3:4);
%                 baselines_C=baselines(1:2);
%                 baselines_L=baselines(3:4);

                stimulated_C=stimulated(1:3);
                stimulated_L=stimulated(4:5);
                baselines_C=baselines(1:3);
                baselines_L=baselines(4:5);

                %baselines = DAG_arrange_baselines(baselines_in,stimulated);
    case 'Test'
        monkey='Linus';
        
        baselines_in{1} = [{[GLO.drive ':\Data\Linus_microstim\20150930']},{8}]; %% actually NO baseline
        stimulated{1} = [{[GLO.drive ':\Data\Linus_microstim\20150930']},{9}];
    case 'TestCornelius'
        
        baselines_in{1} = [{[GLO.drive ':\Data\Cornelius\20171106']},{8}]; %% actually NO baseline
        stimulated{1} = [{[GLO.drive ':\Data\Cornelius\20171106']},{3,4}];
        
        %        batch_title='upper ventral pulvinar control memory';
        %
        %         stimulated{1} = [{[GLO.drive ':\Data\Curius_microstim\20150522']},{16}];
        %         baselines_in{1} = [{[GLO.drive ':\Data\Curius_microstim\20150522']},{16}]; %% actually NO baseline
        %
        %         stimulated{2} = [{[GLO.drive ':\Data\Curius_microstim\20150529']},{12};{[GLO.drive ':\Data\Curius_microstim\20150529']},{13}];
        %         baselines_in{2} = [{[GLO.drive ':\Data\Curius_microstim\20150529']},{19};{[GLO.drive ':\Data\Curius_microstim\20150529']},{20}];
        
        %
        
        
        % %%
        %         batch_title=' direct, lin 20150709, 5,3, current 250, 49.5deep';
        %         stimulated{1} = [{[GLO.drive ':\Data\Linus_microstim\20150709']},{10}];
        %         baselines_in{1} = [{[GLO.drive ':\Data\Linus_microstim\20150709']},{8};{[GLO.drive ':\Data\Linus_microstim\20150709']},{9}];
        %         baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
        %         batch_title=' direct, lin 20150710, 6,3, current 250, 51deep';
        %         stimulated{1} = [{[GLO.drive ':\Data\Linus_microstim\20150710']},{9}];
        %         baselines_in{1} = [{[GLO.drive ':\Data\Linus_microstim\20150710']},{7}];
        %         baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
        %         batch_title=' direct, lin 20150730, 6,3, current 250, 51deep';
        %         stimulated{1} = [{[GLO.drive ':\Data\Linus_microstim\20150730']},{8}];
        %         baselines_in{1} = [{[GLO.drive ':\Data\Linus_microstim\20150730']},{7}];
        %         baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
        %         batch_title=' direct, lin 20150730, 6,3, current 150, 51deep';
        %         stimulated{1} = [{[GLO.drive ':\Data\Linus_microstim\20150730']},{9}];
        %         baselines_in{1} = [{[GLO.drive ':\Data\Linus_microstim\20150730']},{7}];
        %         baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
        %         batch_title=' direct, lin 20150731, 6,3, current 250, 51deep';
        %         stimulated{1} = [{[GLO.drive ':\Data\Linus_microstim\20150731']},{8}];
        %         baselines_in{1} = [{[GLO.drive ':\Data\Linus_microstim\20150731']},{7}];
        %         baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
        %         batch_title=' direct, lin 20150731, 6,3, current 200, 50deep';
        %         stimulated{1} = [{[GLO.drive ':\Data\Linus_microstim\20150731']},{14}];
        %         baselines_in{1} = [{[GLO.drive ':\Data\Linus_microstim\20150731']},{7}];
        %         baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
        %         batch_title=' direct, lin 20150805, 5,3, current 250, 50.75deep';
        %         stimulated{1} = [{[GLO.drive ':\Data\Linus_microstim\20150805']},{8}];
        %         baselines_in{1} = [{[GLO.drive ':\Data\Linus_microstim\20150805']},{7}];
        %         baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
        %         batch_title=' direct, lin 20150805, 5,3, current 200, 51.75deep';
        %         stimulated{1} = [{[GLO.drive ':\Data\Linus_microstim\20150805']},{14}];
        %         baselines_in{1} = [{[GLO.drive ':\Data\Linus_microstim\20150805']},{7}];
        %         baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
        
        %         batch_title=' direct, lin 20150806, 5,3, current 250, 51deep';
        %         stimulated{1} = [{[GLO.drive ':\Data\Linus_microstim\20150806']},{8}];
        %         baselines_in{1} = [{[GLO.drive ':\Data\Linus_microstim\20150806']},{7}];
        %         baselines = DAG_arrange_baselines(baselines_in,stimulated);
        % %
        %         batch_title=' direct, lin 20150806, 5,3, current 250, 52deep';
        %         stimulated{1} = [{[GLO.drive ':\Data\Linus_microstim\20150806']},{14}];
        %         baselines_in{1} = [{[GLO.drive ':\Data\Linus_microstim\20150806']},{7}];
        %         baselines = DAG_arrange_baselines(baselines_in,stimulated);
        
        
        %         batch_title=' direct, lin 20150807, 6,3, current 250, 48.25deep';
        %         stimulated{1} = [{[GLO.drive ':\Data\Linus_microstim\20150807']},{8}];
        %         baselines_in{1} = [{[GLO.drive ':\Data\Linus_microstim\20150807']},{7}];
        %         baselines = DAG_arrange_baselines(baselines_in,stimulated);
        % %
        %         batch_title=' direct, lin 20150807, 6,3, current 250, 52.25deep';
        %         stimulated{1} = [{[GLO.drive ':\Data\Linus_microstim\20150807']},{14}];
        %         baselines_in{1} = [{[GLO.drive ':\Data\Linus_microstim\20150807']},{7}];
        %         baselines = DAG_arrange_baselines(baselines_in,stimulated);
        %
    
    case 'Shorter_duration_100ms'
        
             
       monkey='Curius';
        batch_title='20131114';
        Aiminput={'fix_pol_dir'};
        GLO.effector_to_use=0;
        Inputslist={};

     stimulated{1} = [{[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{21};
                      {[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{22};
                      {[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{23};
                      {[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{24};
                      {[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{25};
                      {[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{26}];
     baselines{1} = [{[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{21}]; %% actually NO baseline
    
    case 'Curius_inc'
        
        monkey='Curius';
        batch_title='20161017';
        Aiminput={'inc_mem_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
        stimulated{1} = [{[GLO.drive ':\Data\Curius_microstim_with_parameters\20150410']},{15}];
         baselines{1} = [{[GLO.drive ':\Data\Curius_microstim_with_parameters\20150806']},{17}];
            
     case 'Curius_dir_for_inc'
        
        monkey='Curius';
        batch_title='20161017';
        Aiminput={'use_dir_sti_at_win'};
        GLO.effector_to_use=0;
        Inputslist={};
        
        stimulated{1} = [{[GLO.drive ':\Data\Curius_microstim_with_parameters\20150410']},{12}];
         baselines{1} = [{[GLO.drive ':\Data\Curius_microstim_with_parameters\20150410']},{11}];
            
                       
     
    case 'Shorter_duration_150ms'
        
       monkey='Curius';
        batch_title='20131114';
        Aiminput={'fix_pol_dir'};
        GLO.effector_to_use=0;
        Inputslist={};

     stimulated{1} = [{[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{27};
                      {[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{28};
                      {[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{29};
                      {[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{30};
                      {[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{31};
                      {[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{32}];
     baselines{1} = [{[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{21}]; %% actually NO baseline
     
    case 'Shorter_duration_control'
        
       monkey='Curius';
        batch_title='20131114';
        Aiminput={'fix_pol_dir'};
        GLO.effector_to_use=0;
        Inputslist={};

     stimulated{1} = [{[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{11};
                      {[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{12};
                      {[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{13};
                      {[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{14};
                      {[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{15};
                      {[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{16};
                      {[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{17};
                      {[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{18}];
     baselines{1} = [{[GLO.drive ':\Data\Curius_microstim_with_parameters\20131114']},{21}]; %% actually NO baseline
end



if GLO.table
    all_batches_stim=[];
    all_batches_baseline=[];
    for batch=1:numel(stimulated)
        all_batches_stim= [all_batches_stim;ones(size(stimulated{batch},1),1)*batch];
    end
    for batch=1:numel(baselines)
        all_batches_baseline= [all_batches_baseline;ones(size(baselines{batch},1),1)*batch];
    end
    if strcmp(monkey,'Combined')
        all_monkeys={'Curius','Linus'};
        N_batches_C_stim=sum(cellfun(@(x)  size(x,1),stimulated_C));
        N_batches_L_stim=sum(cellfun(@(x)  size(x,1),stimulated_L));
        N_batches_C_baseline=sum(cellfun(@(x)  size(x,1),baselines_C));
        N_batches_L_baseline=sum(cellfun(@(x)  size(x,1),baselines_L));
        temp_batch_assignment_stim={all_batches_stim(1:N_batches_C_stim),all_batches_stim(N_batches_C_stim+1:N_batches_C_stim+N_batches_L_stim)};
        temp_batch_assignment_baseline={all_batches_baseline(1:N_batches_C_baseline),all_batches_baseline(N_batches_C_baseline+1:N_batches_C_baseline+N_batches_L_baseline)};
        
        temp_stim={vertcat(stimulated_C{:}),vertcat(stimulated_L{:})};
        temp_no_stim={vertcat(baselines_C{:}),vertcat(baselines_L{:})};
    else
        all_monkeys={monkey};
        temp_stim={vertcat(stimulated{:})};
        temp_batch_assignment_stim={all_batches_stim};
        temp_batch_assignment_baseline={all_batches_baseline};
        temp_no_stim={vertcat(baselines{:})};
    end
    data_ui_stim={};
    data_ui_nostim={};
    for m=1:numel(all_monkeys)
        temp_table_title='Stimulated runs';
        [table_out_stim] = table_session_by_session(temp_stim{m},temp_table_title,all_monkeys{m},temp_batch_assignment_stim{m});
        data_ui_stim=[data_ui_stim;table_out_stim.data_ui];
        
        temp_table_title='Baseline runs';
        [table_out_no_stim] = table_session_by_session(temp_no_stim{m},temp_table_title,all_monkeys{m},temp_batch_assignment_baseline{m});
        data_ui_nostim=[data_ui_nostim;table_out_no_stim.data_ui];
    end
    table_out_stim.data_ui  =data_ui_stim;
    table_out_no_stim.data_ui   =data_ui_nostim;
    table_out_stim.monkey       =monkey;
    table_out_no_stim.monkey    =monkey;
end

a=1;


%% REACHES + SACCADES DATASET
switch GLO.type_to_use
    case 2
        batch_title=[batch_title, ' direct'];
    case 3
        batch_title=[batch_title, ' memory'];
end
switch GLO.effector_to_use
    case 0
        batch_title=[batch_title, ' saccades'];
    case 4
        batch_title=[batch_title, ' reaches'];
end


%%
switch monkey
    case 'Linus'
        temp_dir=[GLO.drive ':\Projects\Pulv_microstim_behavior\behavior\Linus_summaries'];
    case 'Curius'
        temp_dir=[GLO.drive ':\Projects\Pulv_microstim_behavior\behavior\Curius_summaries'];
    case 'Combined'
        temp_dir=[GLO.drive ':\Projects\Pulv_microstim_behavior\behavior\Combined_summaries'];
    case 'Test'
        temp_dir=[GLO.drive ':\Projects\Pulv_microstim_behavior\behavior\Test_summaries'];
    case 'Test2'
        temp_dir=[GLO.drive ':\Projects\Pulv_microstim_behavior\behavior\Test2_summaries'];
    case 'Test3'
        temp_dir=[GLO.drive ':\Projects\Pulv_microstim_behavior\behavior\Test3_summaries'];
    case 'Test4'
        temp_dir=[GLO.drive ':\Projects\Pulv_microstim_behavior\behavior\Test4_summaries'];
end
summary_folder=[temp_dir filesep];

if GLO.run_analyze==0
    load([summary_folder filesep 'MPA_out_' num2str(saccade_definition) '.mat'],[GLO.mon_to_test '_stim'],[GLO.mon_to_test '_baseline']);
    out_comp_stim=eval([GLO.mon_to_test '_stim']);
    out_comp_baseline=eval([GLO.mon_to_test '_baseline']);
end

if GLO.run_internal_calculation==0
    load([summary_folder filesep 'batch_out_' num2str(saccade_definition) '.mat'],[GLO.mon_to_test '_batch_out']);
    batch_out=eval([GLO.mon_to_test '_batch_out']);
    if GLO.plot_eye_traces
        [eye_x eye_y time] = eye_traces(out_comp_stim,batch_out.target_positions,batch_title,errors_only_for_eye_trace_plot);
    end
end

if GLO.run_internal_calculation
    GLO.plot_to_save = GLO.tests{1,1};
    GLO.excentricity =  GLO.tests{1,2};
   
    sel_array=repmat({GLO.Sel_all},1,size(stimulated,2));
    Analyze_input=[stimulated;sel_array];
    %     [out_comp_stim,~,~]= monkeypsych_analyze_20140910(Analyze_input{:});
    if GLO.run_analyze
        [out_comp_stim,~,~]= monkeypsych_analyze_working(Analyze_input{:});
    end
    all_target_positions=[];
    all_fix_positions=[];
    all_tar_positions=[];
    for k=1:numel(out_comp_stim)
        switch GLO.effector_to_use
            case 0
                all_target_positions=[all_target_positions [out_comp_stim{k}.saccades.tar_pos]-[out_comp_stim{k}.saccades.fix_pos] ];
                all_fix_positions=[all_fix_positions [out_comp_stim{k}.saccades.fix_pos] ];
                all_tar_positions=[all_tar_positions [out_comp_stim{k}.saccades.tar_pos]];
            case 4
                all_target_positions=[all_target_positions [out_comp_stim{k}.reaches.tar_pos]-[out_comp_stim{k}.reaches.fix_pos]];
                all_fix_positions=[all_fix_positions [out_comp_stim{k}.reaches.fix_pos]];
                all_tar_positions=[all_tar_positions [out_comp_stim{k}.reaches.tar_pos]];
        end
    end
    target_positions          =unique(all_target_positions(~isnan(all_target_positions)));
    n_targets=numel(target_positions);
    for t=1:n_targets
        target_positions(abs(target_positions-target_positions(t))<GLO.target_pos_precision)=target_positions(t);
    end
    batch_out.target_positions          =unique(target_positions);
    if GLO.plot_eye_traces
        [eye_x eye_y time] = eye_traces(out_comp_stim,batch_out.target_positions,batch_title,errors_only_for_eye_trace_plot);
    end
    
    
    for i=1:numel(stimulated)
        if isempty(stimulated{i})
            bias_B_trial{i}={NaN};
            mean_rt{i}=NaN(numel(GLO.Labels),4);
            mean_x_acc{i}=NaN(numel(GLO.Labels),4);
            bias_S{i}={NaN(1,numel(GLO.windows_to_use))};
            idx_stim_early{i}={NaN};
            idx_stim_late{i}={NaN};
            moving_average{i}={NaN};
        else
            %[bias_B_trial{i},bias_S{i},mean_rt{i},raw_rt{i},num_hits{i},mean_x_acc{i},early_late_index{i},idx_stim_early{i},idx_stim_late{i},per_tar_pos_stim{i},moving_average{i}]=batch_criteria(out_comp_stim(i),batch_out.target_positions,monkey);
            [bias_B_trial{i},bias_S{i},mean_rt{i},raw_rt{i},num_hits{i},mean_x_acc{i},early_late_index{i},idx_stim_early{i},idx_stim_late{i},per_tar_pos{i},moving_average{i},velocities{i},acc_xy{i},acc_rad{i},N{i}] = bias_and_rt_internal_calculation(out_comp_stim(i),batch_out.target_positions);
        end
    end
    
    batch_out.per_tar_pos           = per_tar_pos;
    batch_out.bias_B_trial          = bias_B_trial;
    batch_out.mean_rt               = mean_rt;
    batch_out.raw_rt                = raw_rt;
    batch_out.num_hits              = num_hits;
    batch_out.mean_x_acc            = mean_x_acc;
    batch_out.moving_average_t      = moving_average;
    batch_out.velocities            = velocities;
    batch_out.accuracy_rad          = acc_rad;
    batch_out.accuracy_xy           = acc_xy;
    batch_out.N                     = N;
    
    batch_out.bias_S                = bias_S;
    
    batch_out.early_late_index      = early_late_index;
    batch_out.idx_stim_early        = idx_stim_early;
    batch_out.idx_stim_late         = idx_stim_late;
    
    Analyze_input=[baselines;sel_array];
    %     [out_comp_baseline,~,~]= monkeypsych_analyze_20140910(Analyze_input{:});
    if GLO.run_analyze
        [out_comp_baseline,~,~]= monkeypsych_analyze_working(Analyze_input{:});
    end
    for j=1:numel(baselines)
        if isempty(baselines{j})
            bias_B_run{j}=NaN;
            mean_rt_B_run{j}=NaN(numel(GLO.Labels),4);
            mean_x_acc_B_run{j}=NaN(numel(GLO.Labels),4);
            num_hits_B_run{j}=zeros(4,numel(GLO.Labels));
            raw_rt_B_run{j}= num2cell(repmat({NaN},numel(GLO.Labels),4));
            bias_B_run{j}= {NaN};
            moving_average_B_run{j}.L_choice=num2cell(NaN(1,numel(GLO.Labels)));
            moving_average_B_run{j}.all_choice=num2cell(NaN(1,numel(GLO.Labels)));
            moving_average_B_run{j}.session=NaN;
            moving_average_B_run{j}.run=NaN;
            velocities_B_run{j}.mean=NaN(numel(GLO.Labels),4);
            velocities_B_run{j}.std=NaN(numel(GLO.Labels),4);
            velocities_B_run{j}.raw=num2cell(NaN(numel(GLO.Labels),4));
            acc_rad_B_run{j}=velocities_B_run{j};
            acc_rad_B_run{j}.mean_eu=acc_rad_B_run{j}.mean;
            acc_rad_B_run{j}.std_eu=acc_rad_B_run{j}.std;
            acc_rad_B_run{j}.raw_eu=acc_rad_B_run{j}.raw;
            acc_xy_B_run{j}=acc_rad_B_run{j};
            N.multistep{j}=num2cell(NaN(numel(GLO.Labels),4));
            %N.total{j}=NaN(numel(GLO.Labels),4);
        else
            %[bias_B_run{j},~,mean_rt_B_run{j},raw_rt_B_run{j},num_hits_B_run{j},mean_x_acc_B_run{j},~,~,~,per_tar_pos_base{j},moving_average{j}]=batch_criteria(out_comp_baseline(j),batch_out.target_positions,monkey);
            [bias_B_run{j},~,mean_rt_B_run{j},raw_rt_B_run{j},num_hits_B_run{j},mean_x_acc_B_run{j},~,~,~,per_tar_pos_base{j},moving_average_B_run{j},velocities_B_run{j},acc_xy_B_run{j},acc_rad_B_run{j},N{j}]=bias_and_rt_internal_calculation(out_comp_baseline(j),batch_out.target_positions);
        end
    end
    
    batch_out.per_tar_pos_base      = per_tar_pos_base;
    batch_out.bias_B_run            = bias_B_run;
    batch_out.mean_rt_B_run         = mean_rt_B_run;
    batch_out.raw_rt_B_run          = raw_rt_B_run;
    batch_out.num_hits_B_run        = num_hits_B_run;
    batch_out.mean_x_acc_B_run      = mean_x_acc_B_run;
    batch_out.moving_average_r      = moving_average_B_run;
    batch_out.velocities_B_run      = velocities_B_run;
    batch_out.accuracy_rad_B_run    = acc_rad_B_run;
    batch_out.accuracy_xy_B_run     = acc_xy_B_run;    
    batch_out.N_B_run               = N;
    
    if GLO.save_analyze_output
        to_save.([GLO.mon_to_test '_stim'])=out_comp_stim;
        to_save.([GLO.mon_to_test '_baseline'])=out_comp_baseline;
        if exist([summary_folder 'MPA_out_' num2str(saccade_definition) '.mat'],'file')
            save([summary_folder 'MPA_out_' num2str(saccade_definition)],'-struct','to_save','-append')
        else
            save([summary_folder 'MPA_out_' num2str(saccade_definition)],'-struct','to_save')
        end
    end
    
    if GLO.save_batch_for_later
        to_save.([GLO.mon_to_test '_batch_out'])=batch_out;
        if exist([summary_folder 'batch_out_' num2str(saccade_definition) '.mat'],'file')
            save([summary_folder 'batch_out_' num2str(saccade_definition)],'-struct','to_save','-append');
        else
            save([summary_folder 'batch_out_' num2str(saccade_definition)],'-struct','to_save');
        end
    end
end

if GLO.table
    save([summary_folder  GLO.mon_to_test], 'table_out_stim');
end


batch_out = per_session_significance(batch_out);
batch_calculations(batch_out,monkey,batch_title);
if  GLO.plot_batches_individually
    for i=1:numel(stimulated)
        %close all
        batch_out_per_batch.target_positions        = batch_out.target_positions;
        batch_out_per_batch.per_tar_pos             = batch_out.per_tar_pos(i);
        batch_out_per_batch.bias_B_trial            = batch_out.bias_B_trial(i);
        batch_out_per_batch.bias_S                  = batch_out.bias_S(i);
        batch_out_per_batch.mean_rt                 = batch_out.mean_rt(i);
        batch_out_per_batch.raw_rt                  = batch_out.raw_rt(i);
        batch_out_per_batch.num_hits                = batch_out.num_hits(i);
        batch_out_per_batch.mean_x_acc              = batch_out.mean_x_acc(i);
        batch_out_per_batch.idx_stim_early          = batch_out.idx_stim_early(i);
        batch_out_per_batch.idx_stim_late           = batch_out.idx_stim_late(i);
        batch_out_per_batch.per_tar_pos_base        = batch_out.per_tar_pos_base(i);
        batch_out_per_batch.bias_B_run              = batch_out.bias_B_run(i);
        batch_out_per_batch.mean_rt_B_run           = batch_out.mean_rt_B_run(i);
        batch_out_per_batch.mean_x_acc_B_run        = batch_out.mean_x_acc_B_run(i);
        batch_out_per_batch.raw_rt_B_run            = batch_out.raw_rt_B_run(i);
        batch_out_per_batch.num_hits_B_run          = batch_out.num_hits_B_run(i);
        batch_out_per_batch.moving_average_r        = batch_out.moving_average_r(i);
        batch_out_per_batch.moving_average_t        = batch_out.moving_average_t(i);
        batch_out_per_batch.early_late_index        = batch_out.early_late_index(i);
        batch_out_per_batch.velocities              = batch_out.velocities(i);
        batch_out_per_batch.velocities_B_run        = batch_out.velocities_B_run(i);
        batch_out_per_batch.accuracy_rad            = batch_out.accuracy_rad(i);
        batch_out_per_batch.accuracy_xy             = batch_out.accuracy_xy(i);
        batch_out_per_batch.accuracy_rad_B_run      = batch_out.accuracy_rad_B_run(i);
        batch_out_per_batch.accuracy_xy_B_run       = batch_out.accuracy_xy_B_run(i);
        
        batch_out_per_batch.N_B_run       = batch_out.N_B_run(i);
        batch_out_per_batch.N       = batch_out.N(i);
        batch_out_per_batch.accuracy_xy_B_run       = batch_out.accuracy_xy_B_run(i);
        batch_out_per_batch.accuracy_xy_B_run       = batch_out.accuracy_xy_B_run(i);
        batch_out_per_batch.accuracy_xy_B_run       = batch_out.accuracy_xy_B_run(i);
        
        batch_title_per_batch                       = [batch_title ', batch ' num2str(i), ', ' stimulated{i}{1}(end-7:end)];
        batch_out_per_batch                         = per_session_significance(batch_out_per_batch);
        batch_calculations(batch_out_per_batch,monkey,batch_title_per_batch);
    end
end

if GLO.table
    
    idx_batch=DAG_find_column_index(table_out_stim.titles_ui,'batch');
    idx_hits=DAG_find_column_index(table_out_stim.titles_ui,'hits');
    idx_offset=DAG_find_column_index(table_out_stim.titles_ui,'offset_right');
    idx_current=DAG_find_column_index(table_out_stim.titles_ui,'current');
    batches=[table_out_stim.data_ui{:,idx_batch}]';
    for b=1:max(batches)
        hits(b,1)=sum([table_out_stim.data_ui{batches==b,idx_hits}]);
        all_offsets{b,1}=unique([table_out_stim.data_ui{batches==b,idx_offset}]);
        all_currents{b,1}=unique([table_out_stim.data_ui{batches==b,idx_current}]);
    end
        
    unique_batch_index=[true; diff(batches)~=0];
    reduced_stim_table=table_out_stim.data_ui(unique_batch_index,:);
    reduced_stim_table(:,idx_hits)=num2cell(hits);
    reduced_stim_table(:,idx_offset)=cellfun(@(x) nanmean(x),all_offsets,'UniformOutput',false);
    reduced_stim_table(:,idx_current)=cellfun(@(x) nanmean(x),all_currents,'UniformOutput',false);
    
    
    reduced_stim_table(:,end+1)=cellfun(@(x) num2str(x),all_offsets,'UniformOutput',false);
    reduced_stim_table(:,end+1)=cellfun(@(x) num2str(x),all_currents,'UniformOutput',false);
    
    table_per_batch.data=[reduced_stim_table, num2cell(GLO.table_per_batch.data)];
    table_per_batch.titles_ui=[table_out_stim.titles_ui 'all_offsets' 'all_currents' GLO.table_per_batch.titles];
    table_per_batch.titles_and_data=vertcat(table_per_batch.titles_ui, table_per_batch.data);
    
    save([summary_folder  GLO.mon_to_test], 'table_per_batch','-append');
    xlswrite([summary_folder  monkey '_effect_sizes_and_locations'],table_per_batch.titles_and_data, batch_title(1:min(31,numel(batch_title))));
end

end

function [choice_bias_trial_baseline,choice_bias_stim,mean_rt,raw_rt,num_hits,mean_x_acc,early_late_index,idx_stim_early,idx_stim_lat,per_tar_pos,moving_average,vel] = batch_criteria(out_comp,target_positions,monkey)
global GLO

excentricity=[];
switch monkey
    case 'Linus'
        temp_dir=[GLO.drive ':\microstim_behavior\Linus_summaries'];
    case 'Curius'
        temp_dir=[GLO.drive ':\microstim_behavior\Curius_summaries'];
    case 'Test'
        temp_dir=[GLO.drive ':\microstim_behavior\Test_summaries'];
end

id_2=1;
if ~isempty( out_comp{1}.emptyflag)
    for idx_tests =  1:size(GLO.tests,2)/2
        GLO.plot_to_save = GLO.tests{1,id_2};
        id_2=id_2+1;
        GLO.excentricity =  GLO.tests{1,id_2};
        id_2=id_2+1;
        [choice_bias_trial_baseline,choice_bias_stim,mean_rt,raw_rt,num_hits,mean_x_acc,early_late_index,idx_stim_early,idx_stim_lat,per_tar_pos,moving_average,vel] = bias_and_rt_internal_calculation(out_comp,target_positions);
    end
    id_2=1;
    cd(temp_dir)
end
end

function baselines_out = DAG_arrange_baselines(baselines_in,stimulated)
for k=1:numel(stimulated)
    baselines_out{k}={};
    for m=1:numel(baselines_in)
        if strcmp(baselines_in{m}{1,1},stimulated{k}{1,1})
            baselines_out{k}=baselines_in{m};
        end
    end
end
end

