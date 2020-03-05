
load('W:\Projects\Pulv_microstim_behavior\behavior\Combined_summaries\Combined_direct_dorsal.mat')
complete_table=table_per_batch.titles_and_data;


idx_z=DAG_find_column_index(complete_table,'z');
Curius_rows=2:16;
Linus_rows=17:31;


%figure('units','normalized','outerposition',[0 0 1 1],'name','Fraction second mode across sessions');
figure
sub_titles={'Contra','Ipsi'};
bins=0:10:100;

windows={'40','80','120'};

for w=1:numel(windows)
    window=windows{w};
idx_SM(1)=DAG_find_column_index(complete_table,[window '_%2ndmode_L']);
idx_SM(2)=DAG_find_column_index(complete_table,[window '_%2ndmode_R']);
idx_bias=DAG_find_column_index(complete_table,[window '_dBS']);

for n=1:2
subplot(numel(windows),2,(w-1)*2+n)
hold on
title(['window: ' window ' ' sub_titles{n}]);

Curius_SM=[complete_table{Curius_rows,idx_SM(n)}]';
Curius_depth=[complete_table{Curius_rows,idx_z}]';
Curius_bias=[complete_table{Curius_rows,idx_bias}]';

Linus_SM=[complete_table{Linus_rows,idx_SM(n)}]';
Linus_depth=[complete_table{Linus_rows,idx_z}]';
Linus_bias=[complete_table{Linus_rows,idx_bias}]';

[R_lin_depth,p_lin_depth]=corr(Linus_depth,Linus_SM,'type','Spearman');
[R_cur_depth,p_cur_depth]=corr(Curius_depth,Curius_SM,'type','Spearman');
[R_lin_bias,p_lin_bias]=corr(Linus_bias,Linus_SM,'type','Spearman');
[R_cur_bias,p_cur_bias]=corr(Curius_bias,Curius_SM,'type','Spearman');


plot(bins,hist(Curius_SM,bins),'r','linewidth',3);
plot(bins,hist(Linus_SM,bins),'b','linewidth',3);
y_lim=get(gca,'ylim');
text(10,y_lim(2)-diff(y_lim)*3/20,sprintf('depth: C: R: %.3f p: %.3f ',R_cur_depth,p_cur_depth),'color','r');
text(10,y_lim(2)-diff(y_lim)*7/20,sprintf('depth: L: R: %.3f p: %.3f ',R_lin_depth,p_lin_depth),'color','b');
text(10,y_lim(2)-diff(y_lim)*11/20,sprintf('bias: C: R: %.3f p: %.3f ',R_cur_bias,p_cur_bias),'color','r');
text(10,y_lim(2)-diff(y_lim)*15/20,sprintf('bias: L: R: %.3f p: %.3f ',R_lin_bias,p_lin_bias),'color','b');
xlabel('fraction second mode [%]');
ylabel('N sessions');
set(gca,'ylim',y_lim);
end
end
export_fig('W:\Projects\Pulv_microstim_behavior\behavior\Combined_summaries\Fraction_second_mode', '-pdf','-transparent') % pdf by run