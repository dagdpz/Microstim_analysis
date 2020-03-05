function update_stimulation_parameters(mastertable_type)
global analysis_parameters
monkey=analysis_parameters.monkey;
working_folder=analysis_parameters.folders.output;
dates=analysis_parameters.dates;
original_parameters_excel_file=strcat(working_folder, analysis_parameters.files.original_parameters);
updated_parameters_excel_file=strcat(working_folder,analysis_parameters.files.updated_parameters);

if exist(updated_parameters_excel_file, 'file') && strcmp(mastertable_type,'trialinfo')
    movefile(updated_parameters_excel_file, strcat(analysis_parameters.folders.output, 'Savecopy' ,filesep, analysis_parameters.files.updated_parameters));
elseif exist(updated_parameters_excel_file, 'file') 
    copyfile(updated_parameters_excel_file, strcat(analysis_parameters.folders.output, 'Savecopy' ,filesep, analysis_parameters.files.updated_parameters));
end

mastertable_to_update=DAG_most_recent_version(working_folder,strcat(monkey, '_', mastertable_type, '_mastertable'));
[~,parametersheets] = xlsfinfo(original_parameters_excel_file);
load(mastertable_to_update);
for k=1:numel(parametersheets)
    if isnan(str2double(parametersheets{k}))
        continue
    elseif str2double(parametersheets{k}) >= dates(1) && str2double(parametersheets{k}) <= dates(2)
        stimulationparametersfile=original_parameters_excel_file;
    else
        stimulationparametersfile=updated_parameters_excel_file;
    end
    [~,~,parameter_table] = xlsread(stimulationparametersfile,parametersheets{k});
    titles_update=[mastertable(1,:)];
    Session_column_idx=ismember(titles_update,'Session');
    Session_row_idx=[false; cell2mat([mastertable(2:end,Session_column_idx)])==str2double(parametersheets{k})];
    completed_table = DAG_update_mastertable_cell(parameter_table,[mastertable(1,:); mastertable(Session_row_idx,:)],2:sum(Session_row_idx)+1);
    
    completed_table_xls=completed_table;
    for col=1:numel(completed_table)
        if isnan(completed_table{col})
            completed_table_xls{col}=[];
        elseif ~ischar(completed_table{col}) && numel(completed_table{col})>1
            kommas=repmat({','},1,numel(completed_table{col}));
            restructured_cell=vertcat(cellstr(num2str(completed_table{col}(:)))',kommas);
            restructured_string=[restructured_cell{:}];
            restructured_string(end)=[];
            completed_table_xls{col}=restructured_string;
        end
    end
    if size(completed_table_xls,2)>255
        completed_table_xls=completed_table_xls(:,1:255);
    end
    xlswrite(updated_parameters_excel_file,completed_table_xls,parametersheets{k});
end

clear complete_mastertable
sheet_counter=0;
complete_mastertable={};
for k=1:numel(parametersheets)
    if isempty (str2num(parametersheets{k}))
            continue
    else
        sheet_counter=sheet_counter+1;   
    end    
    [~,~,parameter_table] = xlsread(updated_parameters_excel_file,parametersheets{k});
    row_counter=size(complete_mastertable,1);
    n_rows=size(parameter_table,1)-1;
    if sheet_counter==1
        complete_mastertable=parameter_table;        
    else        
        table_for_updating=[parameter_table(1,:); cell(row_counter-1,size(parameter_table,2)); parameter_table(2:end,:)];
        complete_mastertable = DAG_update_mastertable_cell(complete_mastertable,table_for_updating,[row_counter+1:row_counter+n_rows]);
    end
    
end

if size(complete_mastertable,2)>255
   complete_mastertable=complete_mastertable(:,1:255);
   disp('mastertable has too many columns, size has been reduced to 255 columns')
end
xlswrite(updated_parameters_excel_file,complete_mastertable,'mastertable');
end


function completed_table = DAG_update_mastertable_cell(original_table,table_for_updating,rows_to_update)
% This function updates all rows specified with raw_data_rows of a original_table (supposedly a cell, with titles in the first row)
% With the corresponding rows of table_for_updating (supposedly a cell, with titles in the first row). 
% dimensions of both tables don't need to fit.
% Order of columns does not matter, since this function compares the title
% names in order to update the correct column (So title names are relevant and should match)
% all rows_to_update should be larger than 1, first row should be titles!


titles_update=table_for_updating(1,:);
titles_original_table=original_table(1,:);
all_rows=1:max(size(original_table,1),size(table_for_updating,1)); % 
logidx_columns=ismember(all_rows,rows_to_update);

completed_table=original_table;
for n_update_column=1:numel(titles_update)
    titles_completed_table=completed_table(1,:);
    all_updated_columns=1:size(titles_completed_table,2);
    if ismember(titles_update{n_update_column},titles_original_table)
        logidx_completed_column=ismember(titles_completed_table,titles_update{n_update_column});
        n_completed_column=all_updated_columns(logidx_completed_column);
    else
        n_completed_column=all_updated_columns(end)+1;
    end
    completed_table(1,n_completed_column)=titles_update(n_update_column);
    completed_table(logidx_columns,n_completed_column)=table_for_updating(logidx_columns,n_update_column);
end
end
