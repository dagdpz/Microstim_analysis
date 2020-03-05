function [table_out] = table_session_by_session(filelist,table_title,monkey,batch_assignment)
global GLO

switch monkey
    case 'Linus'
        path= strcat(GLO.drive,':', filesep,'Data', filesep, monkey, '_microstim_with_parameters', filesep);
        data=[GLO.drive,':', filesep,'Projects' filesep 'Pulv_microstim_behavior' filesep 'behavior', filesep,  monkey, '_summaries', filesep, monkey, '_updated_parameters.xls'];
    case 'Curius'
        path= strcat(GLO.drive,':', filesep,'Data', filesep, monkey, '_microstim_with_parameters', filesep);
        data=[GLO.drive,':', filesep,'Projects' filesep 'Pulv_microstim_behavior' filesep 'behavior', filesep,  monkey, '_summaries', filesep, monkey, '_updated_parameters.xls'];
 case 'Cornelius'
        path= strcat(GLO.drive,':', filesep,'Data', filesep, 'Curius', '_microstim_with_parameters', filesep);
        data=[GLO.drive,':', filesep,'Projects' filesep 'Pulv_microstim_behavior' filesep 'behavior', filesep,  'Curius', '_summaries', filesep, 'Curius', '_updated_parameters.xls'];
end

[num,masterstring,raw] = xlsread(data,'mastertable');
masternum=[NaN(1,size(num,2));num];
idx_Session=DAG_find_column_index(masterstring,'Session');
idx_Run=DAG_find_column_index(masterstring,'Run');
idx_hits=DAG_find_column_index(masterstring,'hits');
idx_grid_x=DAG_find_column_index(masterstring,'gridhole_x');
idx_grid_y=DAG_find_column_index(masterstring,'gridhole_y');
idx_Electrode_depth=DAG_find_column_index(masterstring,'Electrode_depth');
idx_target=DAG_find_column_index(masterstring,'target');
idx_Aim=DAG_find_column_index(masterstring,'Aim');
idx_Offset_right=DAG_find_column_index(masterstring,'Offset_right');
idx_percentage_left_chosen_baseline=DAG_find_column_index(masterstring,'percentage_left_chosen_baseline');
Current_strength=DAG_find_column_index(masterstring,'Current_strength');
idx_Task_notes=DAG_find_column_index(masterstring,'notes');
Impedance_start_kilo_ohms=DAG_find_column_index(masterstring,'Impedance_start_kilo_ohms');
log_sel_equ_idx=ones(1,size(masternum,1));

for idx_session = 1:size(filelist,1)
    session_date_temp{idx_session,:}=filelist{idx_session,1};
    session_date(idx_session,:)=str2num(session_date_temp{idx_session}(end-7:end));
    run_number(idx_session,:)=filelist{idx_session,2};
    Inputsequal={'Session',session_date(idx_session,:);'Run',run_number(idx_session,:)};
    
    for equalidx=1:size(Inputsequal,1)
        idx.(Inputsequal{equalidx,1})=DAG_find_column_index(masterstring,Inputsequal{equalidx,1});
        if ischar(Inputsequal{equalidx,2})
            log_sel_equ_idx_cell=strfind([masterstring(:,idx.(Inputsequal{equalidx,1}))],Inputsequal{equalidx,2});
            log_sel_equ_idx(equalidx,:) = ~cellfun('isempty',log_sel_equ_idx_cell);
        else
            log_sel_equ_idx(equalidx,:)=masternum(:,idx.(Inputsequal{equalidx,1}))==Inputsequal{equalidx,2};
        end
    end
    
    Sel_all=all(log_sel_equ_idx,1);
    for master_idx = 1: size(masternum,1)
        if  Sel_all(master_idx)
           
            penetration_date{idx_session,1}                       = masternum(master_idx,idx_Session);
            cell_run{idx_session,1}                               = masternum(master_idx,idx_Run);
            aim{idx_session,1}                                    = raw{master_idx,idx_Aim};
            target{idx_session,1}                                 = raw{master_idx,idx_target};
            current{idx_session,1}                                = masternum(master_idx,Current_strength);
            impedance{idx_session,1}                              = masternum(master_idx,Impedance_start_kilo_ohms);
            x{idx_session,1}                                      = masternum(master_idx,idx_grid_x);
            y{idx_session,1}                                      = masternum(master_idx,idx_grid_y);
            z{idx_session,1}                                      = masternum(master_idx,idx_Electrode_depth);
            offset_right{idx_session,1}                           = masternum(master_idx,idx_Offset_right);
            offset_right{idx_session,1}                           = masternum(master_idx,idx_Offset_right);
            batch{idx_session,1}                                  = batch_assignment(idx_session);            
            hits{idx_session,1}                                   = masternum(master_idx,idx_hits);
            percentage_left_chosen_baseline{idx_session,1}        = masternum(master_idx,idx_percentage_left_chosen_baseline);
            notes{idx_session,1}                                  = raw{master_idx,idx_Task_notes};
            
            idx_session                                           = idx_session + 1;
        end
    end
end


table_out.titles_ui =[{'penetration_date'}, {'cell_run'}, {'aim'}, {'target'}, {'current'}, {'impedance'}, {'x'},{'y'},{'z'}, {'offset_right'},{'batch'}, {'hits'}, {'percentage_left_chosen_baseline'}, {'notes'}];
table_out.data_ui =[cellfun(@(var) eval(var),table_out.titles_ui,'UniformOutput',0)];
table_out.data_ui=[table_out.data_ui{:}];
%table_out.data_ui =[penetration_date, cell_run, aim, target, current, impedance, x,y,z, offset_right, percentage_left_chosen_baseline, notes];
table_out.monkey=monkey;
h=cell2uitable(table_out.data_ui, table_out.titles_ui, table_title);

    function column_index=DAG_find_column_index(inputcell,title)
        for m=1:size(inputcell,2)
            if strcmp(inputcell{1,m},title)
                column_index=m;
            end
        end
    end
end