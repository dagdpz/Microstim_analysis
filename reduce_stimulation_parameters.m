function reduce_stimulation_parameters(filelist,parameterfile)
%filelist=analysis_parameters.filelist_formatted;

reduced_parameters_file='reduced_parameters';

date_folders=vertcat(filelist{:,1});
if strcmp(date_folders(1,end),filesep)
    sessions=str2num(date_folders(:,end-8:end-1));
else    
    sessions=str2num(date_folders(:,end-7:end));
end
runs=vertcat(filelist{:,2});
[~,parametersheets] = xlsfinfo(parameterfile);

for k=1:numel(parametersheets)
    if isnan(str2double(parametersheets{k}))
        continue
    elseif any(sessions==str2double(parametersheets{k}))
        [parameter_table_num,~,parameter_table] = xlsread(parameterfile,parametersheets{k});
        current_runs=runs(sessions==str2double(parametersheets{k}));
        idx_runs=DAG_find_column_index(parameter_table,'Run');
        current_runs_xls=[false; ismember(vertcat(parameter_table_num(:,idx_runs)),current_runs)];
        reduced_parameter_entries=[parameter_table(1,:); parameter_table(current_runs_xls,:)];        
        xlswrite(reduced_parameters_file,reduced_parameter_entries,parametersheets{k});
    end
end


clear complete_mastertable
sheet_counter=0;
complete_mastertable={};
[~,parametersheets] = xlsfinfo(reduced_parameters_file);
for k=1:numel(parametersheets)
    if isempty (str2num(parametersheets{k}))
            continue
    else
        sheet_counter=sheet_counter+1;   
    end    
    [~,~,parameter_table] =xlsread(reduced_parameters_file,parametersheets{k});
    row_counter=size(complete_mastertable,1);
    n_rows=size(parameter_table,1)-1;
    if sheet_counter==1
        complete_mastertable=parameter_table;        
    else        
        table_for_updating=[parameter_table(1,:); cell(row_counter-1,size(parameter_table,2)); parameter_table(2:end,:)];
        complete_mastertable = DAG_update_mastertable_cell(complete_mastertable,table_for_updating,[row_counter+1:row_counter+n_rows]);
    end
    
end       
        xlswrite(reduced_parameters_file,complete_mastertable,'mastertable');
end