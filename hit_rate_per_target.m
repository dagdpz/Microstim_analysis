%% hit rate per target
global GLO

GLO.drive = 'L';
GLO.monkey = 'Curius';
GLO.date = '20150123';
if strcmp(GLO.monkey,'Curius');
GLO.datapath = [GLO.drive ':' filesep 'Data' filesep GLO.monkey filesep 'setup2_microstim' filesep GLO.date];
else 
 GLO.datapath = [GLO.drive ':' filesep 'Data' filesep GLO.monkey filesep GLO.date];   
end
GLO.folder_to_save=strcat(GLO.drive,':', filesep, 'microstim_behavior', filesep, GLO.monkey, '_summaries', filesep, 'detection_task', filesep);

dir_datapath=dir(GLO.datapath);
k=0;
for l=1: length(dir_datapath)
    if ismember('mat', dir_datapath(l).name)
        k=k+1;
        run_full{k} = [GLO.datapath filesep dir_datapath(l).name];
        run(k,:) = run_full{k}(end-5:end-4)
        
        [out_comp variable_to_test counter]= monkeypsych_analyze_20140910({GLO.datapath,str2num(run(k,:))},{'saccade_definition',4,'correct_offset',0,'summary',1,'display',0,'saccade_1bo',1,'show_only_one_sac_per_trial',0,'microstim',1,'inferential_on',0})
        
        
        
        if out_comp{1}.emptyflag==1 || all(~isfinite([out_comp{1}.saccades.tar_pos]))
            continue
        end
        
        tar_in_run= unique([out_comp{1}.saccades.tar_pos]);
        t_tar_in_run = ~isnan(tar_in_run);
        tar_in_run = tar_in_run(t_tar_in_run);
        
        
        for i = 1:numel(tar_in_run)
            
            temp_pos_all(i,:)=[out_comp{1}.saccades.tar_pos]==tar_in_run(i);
            temp_pos_suc(i,:)=[out_comp{1}.saccades.tar_pos]==tar_in_run(i) & [out_comp{1}.binary.success];
            all(i)=sum(temp_pos_all(i,:));
            suc(i)=sum(temp_pos_suc(i,:));
            labels_positions{1,i}=[sprintf('%8.2f\n',real(tar_in_run(i))) sprintf('%8.2f\n',imag(tar_in_run(i)))];
        end
        
        
        
        percentage_explored=suc./all*100;
        min_hit_rate = min(percentage_explored);
        min_position=percentage_explored==min(percentage_explored);
        
        figure('units','normalized','outerposition',[0 0 1 1],'name','Significance across sessions')
        subplot(2,2,1)
        bar(percentage_explored)
        title('Hit rate per position')
        ylabel('Percenatge')
        xlabel('Position')
        set(gca,'ylim',[0 100],'xticklabel',labels_positions)
        text(5,-15,['lowest hit rating position' labels_positions(min_position)])
        
        subplot(2,2,2)
        hold on
        scatter(real(tar_in_run),imag(tar_in_run),2000,'r','o')
        scatter(real(tar_in_run(min_position)),imag(tar_in_run(min_position)),2000,'b','o')
        
        for pos=1:numel(tar_in_run)
            text(real(tar_in_run(pos))-1,imag(tar_in_run(pos)),[num2str(percentage_explored(pos)) '%'],'fontsize',14)
        end
        
        if  sum(percentage_explored) ~= 0
            subplot(2,2,3)
            pie(percentage_explored/100,labels_positions)
            title('Normalized hit rate')
            
        end
        color_trial_1 = num2str(out_comp{1}.saccades(1).col_dim(1));
        mtit([GLO.datapath, ' run ' num2str(run(k,:)) ' RGB_R value ' color_trial_1], 'xoff', -0.0, 'yoff', 0.04, 'color', [0 0 0],'Interpreter', 'none');
%         export_fig([GLO.folder_to_save, GLO.datapath(end-7:end), ' run ', num2str(run(k,:)), ' RGB_R value ', color_trial_1], '-pdf', '-transparent') % pdf by run
%         export_fig([GLO.folder_to_save, GLO.datapath(end-7:end)], '-pdf', '-transparent', '-append') % pdf by run
        %      export_fig([GLO.folder_to_save, GLO.datapath(end-7:end) 'sum'], '-pdf', '-transparent', '-append') % pdf by run
        clear temp_pos_all temp_pos_suc all suc labels_positions percentage_explored color_trial_1 min_position
    end
end