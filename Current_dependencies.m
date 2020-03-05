%monkey='Linus';
monkey='Combined';
%[num,str,all]=xlsread(['C:\Users\lschneider\Dropbox\DAG\microstim_behavior\' monkey '_microstim_dpz\' monkey '_effect_sizes_and_locations.xls'],'dorsal pulvinar direct saccades');
[num,str,all]=xlsread(['\\172.16.9.172\dag\microstim_behavior\' monkey '_summaries\' monkey '_effect_sizes_and_locations.xls'],'dorsal pulvinar direct saccades');


effects_to_look_at=all(1,:);
x_parameters={'current','offset_right'};
for x_val= 1:numel(x_parameters)
    x_to_correlate=num(:,DAG_find_column_index(all,x_parameters{x_val}));
    for eff=1:numel(effects_to_look_at)
        idx.eff=DAG_find_column_index(all,effects_to_look_at{eff});
        y_to_correlate=num(:,idx.eff);
        no_nan_idx=~isnan(x_to_correlate) & ~isnan(y_to_correlate);
        if any(no_nan_idx)
        [RHO,PVAL] = corr([x_to_correlate(no_nan_idx),y_to_correlate(no_nan_idx)],'type','spearman');
        else
          RHO=NaN(2);
          PVAL=NaN(2);          
        end
        figure
        scatter(x_to_correlate(no_nan_idx),y_to_correlate(no_nan_idx));
        xlabel(x_parameters{x_val})
        ylabel(effects_to_look_at{eff}, 'interpreter', 'none')
        title([effects_to_look_at{eff} ' R= ' num2str(RHO(1,2)) ', P= ' num2str(PVAL(1,2))],'interpreter','none')
        if exist(['\\172.16.9.172\dag\microstim_behavior\' monkey '_summaries\AACurrent_dependencies.pdf'],'file')
            export_fig(['\\172.16.9.172\dag\microstim_behavior\' monkey '_summaries\AACurrent_dependencies'],'-pdf','-transparent','-append');
        else
            export_fig(['\\172.16.9.172\dag\microstim_behavior\' monkey '_summaries\AACurrent_dependencies'],'-pdf','-transparent');
        end
        close(gcf)
    end
end