function evoked_saccade_plotting_RL
run set_settings
global analysis_parameters

monkey                      = analysis_parameters.monkey;
batch_name                  = analysis_parameters.batches.name;
batch_title                 = analysis_parameters.batches.title;
plot_figures                = analysis_parameters.plot_figures;
batch_processing            = analysis_parameters.batch_processing;
filelist_formatted          = analysis_parameters.filelist_formatted;
temp_dir=pwd;
STATE=load_state;

if batch_processing
    load(DAG_most_recent_version(pwd,'histograms_batched'));
else
    load(DAG_most_recent_version(pwd,'saccade_histograms'));
end

 valid_batches=saccade_histograms(1).valid_batches;

type_excentricities=fieldnames(saccade_histograms);
for p = valid_batches;    
    for T=1:numel(type_excentricities)
        current=saccade_histograms(p).(type_excentricities{T});
        if ~isstruct(current)
            continue;
        end
        excentricity        =current.excentricity;
        first_3=type_excentricities{T}(1:3);
        
        if strcmp(excentricity,'_all')
            far_targets         =NaN;
        elseif strcmp(excentricity,'_close')            
            far_targets         =0;
        else strcmp(excentricity,'_far')            
            far_targets         =1;
        end
        
        if strcmp(first_3,'fix')
            trial_type         =1;
            y_lim              =30;
        elseif strcmp(first_3,'dir')
            trial_type         =2;
            y_lim              =10;
        else strcmp(first_3,'mem')
            trial_type         =3;
            y_lim              =10;
        end
                     
        stimulation_onsets      = current.stimulation_onset;
        stimulation_states      = current.stimulation_state;     
        amplitudes              = current.amplitudes;
        RT_bins                 = current.RT_bins;
        Evo_amp_TH              = current.Evo_amp_TH;        
        skipped_amps            = sum(amplitudes<Evo_amp_TH(1));
        zero_rep                = repmat([0 0 0],skipped_amps,1);
        COLORS                  = [zero_rep;hsv(numel(amplitudes)-skipped_amps-1)];
        newfilecounter          = 1;
        
        for num_stim_window=1:numel(stimulation_onsets)                       
            n.trials                = current.n_trials(num_stim_window) ;            
            n.evoked_l              = current.n.evoked_l(num_stim_window);             
            n.evoked_r              = current.n.evoked_r(num_stim_window);         
            n.return_l              = current.n.return_l(num_stim_window);         
            n.return_r              = current.n.return_r(num_stim_window);        
            n.evoked_and_return     = current.n.evoked_and_return(num_stim_window);
            RTs_histograms_L        = current.RTs_histograms_L(num_stim_window) ;
            RTs_histograms_R        = current.RTs_histograms_R(num_stim_window) ;
            Dir_all                 = current.Dir_all{num_stim_window} ;
            Amp_all                 = current.Amp_all{num_stim_window} ;
            RTs_all                 = current.RTs_all{num_stim_window} ;
            Indexes_all             = current.Indexes_all{num_stim_window} ;
            
            if stimulation_onsets(num_stim_window)>=0 % positive microstimulation start
                additional_string = [' starting from stimulation onset at', '_', STATE.ALL_NAMES{STATE.ALL==stimulation_states(num_stim_window)}, '_+', num2str(stimulation_onsets(num_stim_window).*1000), 'ms '];
            elseif  stimulation_onsets(num_stim_window)<0 % negative microstimulation start
                additional_string = [' starting from stimulation onset at', '_', STATE.ALL_NAMES{STATE.ALL==stimulation_states(num_stim_window)+1}, '_', num2str(stimulation_onsets(num_stim_window).*1000), 'ms '];
            else    % run baseline
                additional_string = ', baseline ';
            end
            
            
            %% Individual plots for each stimulation window
            if ismember(plot_figures,2)
                % Title, Filename and fodler name assignment
                if ~batch_processing
                    current_folder=[temp_dir filesep filelist_formatted{p,1}(4:7) '_' filelist_formatted{p,1}(9:11) '_' filelist_formatted{p,1}(end-7:end)];
                    pdffilename=['run ' num2str(filelist_formatted{p,2}) ' stimlocked histograms' ', ' excentricity(2:end)];
                    current_title=[monkey, '_', filelist_formatted{p,1}(end-7:end), '_', num2str(filelist_formatted{p,2}), ' Saccades/trial vs. direction and time ', additional_string, ', ', excentricity(2:end)];
                else
                    current_title=[batch_title, additional_string, ', ', excentricity(2:end)];
                    pdffilename=[batch_name ', ' excentricity(2:end)];
                    current_folder=[pwd filesep 'histograms'];
                end
                if ~exist(current_folder,'dir')
                    mkdir(current_folder);
                end
                                
                % Plotting
                plot_saccade_prob(RT_bins,n,RTs_histograms_L,RTs_histograms_R,Dir_all,Amp_all,RTs_all,Indexes_all,COLORS,current_title,trial_type,y_lim,skipped_amps,Evo_amp_TH(1),far_targets)
                
                % Appending if file already exists
                if newfilecounter==1
                    export_fig([current_folder filesep pdffilename], '-pdf', '-transparent') % pdf by run
                    newfilecounter=0;
                else
                    export_fig([current_folder filesep pdffilename], '-pdf', '-transparent', '-append') % append to existing pdf
                end
                close gcf
            end
        end
        
        %% summary plot
        if ismember(plot_figures,[1,2])
            % Title and folder name assignment
            if ~batch_processing
                current_title=[monkey, '_', filelist_formatted{p,1}(end-7:end), '_', num2str(filelist_formatted{p,2}), ' Saccades/trial vs. time relative to stimulation onset, ', excentricity(2:end)];
                current_folder=[temp_dir filesep filelist_formatted{p,1}(4:7) '_' filelist_formatted{p,1}(9:11) '_' filelist_formatted{p,1}(end-7:end)];
            else
                current_title=[batch_title, ' all windows relative to stimulation onset, ', excentricity(2:end)];
                current_folder=[pwd filesep 'histograms'];
            end
            if ~exist(current_folder,'dir')
                mkdir(current_folder);
            end
            
            y_lim=y_lim/2;
            % Plotting
            plot_saccade_prob(RT_bins,current.summed_n,current.summed_RTs_histograms_L,current.summed_RTs_histograms_R,current.summed_dir_all,current.summed_amp_all,current.summed_RTs_all,...
                Indexes_all,COLORS,current_title,trial_type,y_lim,skipped_amps,Evo_amp_TH(1),far_targets)
            
            % Save plots
            if ~batch_processing
                export_fig([current_folder filesep 'run ' num2str(filelist_formatted{p,2}) ' stimlocked histograms summary - ' excentricity(2:end)], '-pdf', '-transparent') % pdf by run
                close gcf
            else
                export_fig([current_folder filesep batch_name ' summary - ' excentricity(2:end)], '-pdf', '-transparent') % pdf by run
                saveas(gcf,[current_folder filesep batch_name ' summary - ' excentricity(2:end)]);
                close gcf
            end
        end      
    end
end
end

function plot_saccade_prob(RT_bins,n,RTs_histograms_L,RTs_histograms_R,Dir_all, Amp_all,RTs_all,Indexes_all,COLORS,current_title,trial_type,y_lim,skipped_amps,min_amp,far_targets)
%n_e,n_tot removed
set(gcf,'renderer','painters');
fontsize=24;
linewidth=1;
markersize=8;
amplitudes=[0.5,1,2,4,6,8,10,15,100];
condition_tags=fieldnames(n.trials);
fonttype='Arial';
set(0,'DefaultAxesFontName', fonttype)
n_amps=numel(amplitudes)-1;

% x_lim=max(RT_bins);
if trial_type==1 %fix
    x_lim=600;
    Conditions_to_plot_S=[1 1];
    Conditions_to_plot_B=[2 2];
    Conditions_to_plot_S_dir=1;
    Conditions_to_plot_B_dir=2;
    
    ks_to_plot_L=1;
    ks_to_plot_R=2;
    ks_to_title=[1,2];
    ks_to_lable=2;
    leg1Position = [0.41 0.13 0.06 0.02]; %x_pos y_pos width height % all in proportion to figure size
    leg2Position = [0.85 0.13 0.06 0.02]; %x_pos y_pos width height % all in proportion to figure size
    
    %Labels_RT_comparison={'Contraversive saccades ','Ipsiversive saccades '};
    Labels_RT_comparison={'Saccade probability distribution ','Saccade probability distribution '};
    Labels_polar={{'Saccade endpoints in fixation task'}};
    
    dir_subplot_positions={[2,3]};
    RT_plot_positions={1,1};
    %RT_plot_positions={1,4};
    
    %n_subplots_y=2;
    n_subplots_y=1;
    n_subplots_x=3;
    polar_ticks=[0,30,60,90,120,150,180,210,240,270,300,330,2,4,6,8,10];
    max_amp=6;
    polar_ticks_deg=[2:2:max_amp];
else %dir,mem
    x_lim=400;
    Conditions_to_plot_S=[1 1 3 3 2 2 4 4];
    Conditions_to_plot_B=[5 5 7 7 6 6 8 8];
    Conditions_to_plot_S_dir=[2 1 3 4];
    Conditions_to_plot_B_dir=[6 5 7 8];
    
    ks_to_plot_L=[1 3 5 7];
    ks_to_plot_R=[2 4 6 8];
    ks_to_title=1:4;
    ks_to_lable=5:8;
    leg1Position = [0.93 0.4 0.06 0.35]; %x_pos y_pos width height % all in proportion to figure size
    leg2Position = [0.93 0.15 0.06 0.05]; %x_pos y_pos width height % all in proportion to figure size
    
    Labels_RT_comparison={{'';'Left target acquired ';'(leftward saccades)'}, {'';'Left target acquired ';'(rightward saccades)'}, {'';'Right target acquired ';'(leftward saccades)'},{'';'Right target acquired ';'(rightward saccades)'}};
    Labels_polar={'Left INSTRUCTED ','Left CHOICE', 'Right CHOICE', 'Right INSTRUCTED'};
    
    dir_subplot_positions={9,10,11,12};
    RT_plot_positions={1,2,3,4,5,6,7,8};
    
    n_subplots_y=3;
    n_subplots_x=4;
    polar_ticks=[0,30,60,90,120,150,180,210,240,270,300,330,10,20,30,40];
    if far_targets==0
        max_amp=20;
    else
        max_amp=40;
    end
end

evoked_saccades_figure=figure('units','normalized','outerposition',[0 0 1 1]);
for k = 1:numel(Conditions_to_plot_S)
    subplot(n_subplots_y,n_subplots_x,RT_plot_positions{k})
    if ismember(k,ks_to_plot_L)
        fill([0, 200, 200, 0],[-y_lim -y_lim y_lim y_lim],[1 0.8 0.8],'edgecolor','none');
        fill([200, 400, 400, 200],[-y_lim -y_lim y_lim y_lim],[0.8 1 0.8],'edgecolor','none');
        leftward_RT_B=zeros(size(RTs_histograms_L.(condition_tags{Conditions_to_plot_B(k)})(1,:)));        
        leftward_RT_S=zeros(size(RTs_histograms_L.(condition_tags{Conditions_to_plot_S(k)})(1,:))); 
        for m=(skipped_amps+1):n_amps
            hold on
            %plot(RT_bins,RTs_histograms_L.(condition_tags{Conditions_to_plot_S(k)})(m,:),'-','Color',COLORS(m,:),'LineWidth',linewidth);
            %plot(RT_bins,RTs_histograms_L.(condition_tags{Conditions_to_plot_B(k)})(m,:),'--','Color',COLORS(m,:),'LineWidth',linewidth);
            leftward_RT_B=leftward_RT_B+RTs_histograms_L.(condition_tags{Conditions_to_plot_B(k)})(m,:);
            leftward_RT_S=leftward_RT_S+RTs_histograms_L.(condition_tags{Conditions_to_plot_S(k)})(m,:);
        end
%         leftward_RT_B=leftward_RT_B/sum(leftward_RT_B);
%         leftward_RT_S=leftward_RT_S/sum(leftward_RT_S);
        plot(RT_bins(RT_bins<=0.2)*1000,leftward_RT_B(RT_bins<=0.2)*100,'-','Color',[0.5 0.5 0.5],'LineWidth',linewidth);
        h(1)=plot(RT_bins(RT_bins<=0.2)*1000,leftward_RT_S(RT_bins<=0.2)*100,'-','Color','r','LineWidth',linewidth);
        h(3)=plot(RT_bins(RT_bins>=0.2)*1000,leftward_RT_B(RT_bins>=0.2)*100,'-','Color',[0.5 0.5 0.5],'LineWidth',linewidth);
        h(2)=plot(RT_bins(RT_bins>=0.2)*1000,leftward_RT_S(RT_bins>=0.2)*100,'-','Color','g','LineWidth',linewidth);
    elseif ismember(k,ks_to_plot_R)
        rightward_RT_B=zeros(size(RTs_histograms_R.(condition_tags{Conditions_to_plot_B(k)})(1,:)));   
        rightward_RT_S=zeros(size(RTs_histograms_R.(condition_tags{Conditions_to_plot_S(k)})(1,:)));       
        for m=(skipped_amps+1):n_amps
            hold on
            %plot(RT_bins,RTs_histograms_R.(condition_tags{Conditions_to_plot_S(k)})(m,:),'-','Color',COLORS(m,:),'LineWidth',linewidth);
            %plot(RT_bins,RTs_histograms_R.(condition_tags{Conditions_to_plot_B(k)})(m,:),'--','Color',COLORS(m,:),'LineWidth',linewidth);
            
            rightward_RT_B=rightward_RT_B+RTs_histograms_R.(condition_tags{Conditions_to_plot_B(k)})(m,:);
            rightward_RT_S=rightward_RT_S+RTs_histograms_R.(condition_tags{Conditions_to_plot_S(k)})(m,:);
        end
%         rightward_RT_B=rightward_RT_B/sum(rightward_RT_B);
%         rightward_RT_S=rightward_RT_S/sum(rightward_RT_S);
        plot(RT_bins(RT_bins<=0.2)*1000,-rightward_RT_B(RT_bins<=0.2)*100,'-','Color',[0.5 0.5 0.5],'LineWidth',linewidth);
        plot(RT_bins(RT_bins<=0.2)*1000,-rightward_RT_S(RT_bins<=0.2)*100,'-','Color','r','LineWidth',linewidth);
        plot(RT_bins(RT_bins>=0.2)*1000,-rightward_RT_B(RT_bins>=0.2)*100,'-','Color',[0.5 0.5 0.5],'LineWidth',linewidth);
        plot(RT_bins(RT_bins>=0.2)*1000,-rightward_RT_S(RT_bins>=0.2)*100,'-','Color','g','LineWidth',linewidth);
        
        
        
%         all_RT_B=rightward_RT_B+leftward_RT_B;
%         all_RT_S=rightward_RT_S+leftward_RT_S;        
%         plot(RT_bins(RT_bins<=0.2),all_RT_B(RT_bins<=0.2),'--','Color','r','LineWidth',linewidth);
%         plot(RT_bins(RT_bins<=0.2),all_RT_S(RT_bins<=0.2),'-','Color','r','LineWidth',linewidth);
%         plot(RT_bins(RT_bins>=0.2),all_RT_B(RT_bins>=0.2),'--','Color','b','LineWidth',linewidth);
%         plot(RT_bins(RT_bins>=0.2),all_RT_S(RT_bins>=0.2),'-','Color','b','LineWidth',linewidth);
    end
    if trial_type==1
        set(gca,'xtick',[0,200,400,600]);
        xlabel('Time after stimulation onset(ms)','fontsize',fontsize,'FontName',fonttype)
        ylabel({'Ipsiversive               Contraversive','             Saccade probability (%)'},'FontSize',fontsize,'FontName',fonttype)
    else
        if k==1
            ylabel({'Saccades / Trials ';' (CHOICE)'},'FontSize',fontsize,'FontName',fonttype)
        elseif k==5
            ylabel({'Saccades / Trials ';'(INSTRUCTED)'},'FontSize',fontsize,'FontName',fonttype)
        end
    end
    if ismember(k,ks_to_title)
        title(Labels_RT_comparison{k},'interpreter','none','FontSize',fontsize+2,'FontName',fonttype)
    end
    if ismember(k,ks_to_lable)
        xlabel('Time after stimulation [ms]','fontsize',fontsize,'FontName',fonttype)
    end
    set(gca,'xlim',[0 x_lim],'ylim',[-y_lim, y_lim],'FontSize',fontsize-2,'FontName',fonttype,'LineWidth',linewidth,'ygrid','on')
    if ismember(k,ks_to_plot_R)
        if trial_type~=1
            set(gca,'ytick',[],'FontName',fonttype);
            set(gca,'ycolor',[1 1 1],'FontName',fonttype)
        end
    end
    if ismember(k,ks_to_plot_L)
        text(x_lim/100,y_lim-(y_lim/10),'Stimulated trials = ', 'FontSize',fontsize-2,'FontName',fonttype)
        text(x_lim*0.8,y_lim-(y_lim/10),num2str(n.trials.(condition_tags{Conditions_to_plot_S(k)})), 'FontSize',fontsize-2,'FontName',fonttype);
        text(x_lim/100,y_lim-(2*y_lim/10),'Not stimulated trials = ', 'FontSize',fontsize-2,'FontName',fonttype)
        text(x_lim*0.8,y_lim-(2*y_lim/10),num2str(n.trials.(condition_tags{Conditions_to_plot_B(k)})), 'FontSize',fontsize-2,'FontName',fonttype);
    end
%    NE      =n_e.(condition_tags{Conditions_to_plot_S(k)});
%    Ntot    =n_tot.(condition_tags{Conditions_to_plot_S(k)});
%     NnE     =n_e.(condition_tags{Conditions_to_plot_B(k)});
%     Nntot   =n_tot.(condition_tags{Conditions_to_plot_B(k)});
%    p_evoked_S      =round(NE*100./Ntot)./100;
    
%     if ismember(k,ks_to_plot_L)
%         text(0.01,y_lim-(2*y_lim/10),'p_e_v_o_k_e_d =', 'FontSize',fontsize-2,'FontName',fonttype)
%         text(0.25,y_lim-(2*y_lim/10),num2str(p_evoked_S), 'FontSize',fontsize-2,'FontName',fonttype)
%     end

    text(x_lim/30,-y_lim+(2*y_lim/10),{'Stimulation'; '   period   '}, 'FontSize',fontsize-8,'FontName',fonttype);    
    %set(stim_area,'facealpha',.2,'edgealpha',0);
    
    %line([0.2,0.2]'*1000,[0,y_lim]', 'LineWidth',linewidth, 'Color','k');
    box on
end
Whole_legend={'0.5-1° S','0.5-1° B','1-2° S','1-2° B','2-4° S','2-4° B','4-6° S','4-6° B','6-8° S','6-8° B','8-10° S','8-10° B','10-15° S','10-15° B','>15° S','>15° B'};
if trial_type~=1
    hl=legend(Whole_legend(skipped_amps*2+1:end));
    set(hl,'Position', leg1Position,'Units', 'normalized','FontName',fonttype);
end

for k = 1:numel(Conditions_to_plot_S_dir)
    subplot(n_subplots_y,n_subplots_x,dir_subplot_positions{k})
    
    plot_dir_S=[Dir_all.(condition_tags{Conditions_to_plot_S_dir(k)})]';
    plot_dir_B=[Dir_all.(condition_tags{Conditions_to_plot_B_dir(k)})]';
    plot_amp_S=[Amp_all.(condition_tags{Conditions_to_plot_S_dir(k)})]';
    plot_amp_B=[Amp_all.(condition_tags{Conditions_to_plot_B_dir(k)})]';
    idx_dir_S=plot_amp_S>min_amp;
    idx_dir_B=plot_amp_B>min_amp;
    
    inv=polar(0,max_amp);
    set(inv,'linestyle','none','markerfacecolor','w','color','w');
    hold on
    
    if trial_type==1
        plot_dir_S_early_L=NaN;
        plot_dir_S_late_L=NaN;
        plot_dir_B_early_L=NaN;
        plot_dir_B_late_L=NaN;
        
        plot_dir_S_early_R=NaN;
        plot_dir_S_late_R=NaN;
        plot_dir_B_early_R=NaN;
        plot_dir_B_late_R=NaN;
                
        plot_amp_S_early_L=NaN;
        plot_amp_S_late_L=NaN;
        plot_amp_B_early_L=NaN;
        plot_amp_B_late_L=NaN;
        
        plot_amp_S_early_R=NaN;
        plot_amp_S_late_R=NaN;
        plot_amp_B_early_R=NaN;
        plot_amp_B_late_R=NaN;
        
        N_S_early_L_30 = 0;
        N_S_late_R_30  = 0;
        N_B_early_L_30 = 0;
        N_B_late_R_30  = 0;
        
        N_S_early = 0;
        N_S_late  = 0;
        N_B_early = 0;
        N_B_late  = 0;
        
        
        Amp_evo = NaN;
        RT_evo  = NaN;
        Amp_evo_30    = NaN;
        RT_evo_30    = NaN;
        for m=(skipped_amps+1):numel(amplitudes)-1
            idx_dir_S=plot_amp_S>amplitudes(m) & plot_amp_S<amplitudes(m+1);
            idx_dir_B=plot_amp_B>amplitudes(m) & plot_amp_B<amplitudes(m+1);
            
            index_return_S_L=Indexes_all.L_with_ret.(condition_tags{Conditions_to_plot_S_dir(k)})';
            index_return_S_R=Indexes_all.R_with_ret.(condition_tags{Conditions_to_plot_S_dir(k)})';
            index_return_B_L=Indexes_all.L_with_ret.(condition_tags{Conditions_to_plot_B_dir(k)})';
            index_return_B_R=Indexes_all.R_with_ret.(condition_tags{Conditions_to_plot_B_dir(k)})';
            
            idx_dir_S_L     =(plot_dir_S<-pi/2 | plot_dir_S>pi/2) & idx_dir_S;            
            idx_dir_S_L_30  =(plot_dir_S<-pi*5/6 | plot_dir_S>pi*5/6) & idx_dir_S;
            idx_dir_S_R     =(plot_dir_S>-pi/2 & plot_dir_S<pi/2) & idx_dir_S;
            idx_dir_S_R_30  =(plot_dir_S>-pi*1/6 & plot_dir_S<pi*1/6) & idx_dir_S;
            idx_dir_B_L     =(plot_dir_B<-pi/2 | plot_dir_B>pi/2) & idx_dir_B;
            idx_dir_B_L_30  =(plot_dir_B<-pi*5/6 | plot_dir_B>pi*5/6) & idx_dir_B;
            idx_dir_B_R     =(plot_dir_B>-pi/2 & plot_dir_B<pi/2) & idx_dir_B;
            idx_dir_B_R_30  =(plot_dir_B>-pi*1/6 & plot_dir_B<pi*1/6) & idx_dir_B;
            
            idx_dir_S_early_L     =idx_dir_S_L & RTs_all.(condition_tags{Conditions_to_plot_S_dir(k)})'<0.2;
            idx_dir_S_late_L      =idx_dir_S_L & RTs_all.(condition_tags{Conditions_to_plot_S_dir(k)})'>=0.2;
            idx_dir_B_early_L     =idx_dir_B_L & RTs_all.(condition_tags{Conditions_to_plot_B_dir(k)})'<0.2;
            idx_dir_B_late_L      =idx_dir_B_L & RTs_all.(condition_tags{Conditions_to_plot_B_dir(k)})'>=0.2;
            
            idx_dir_S_early_R     =idx_dir_S_R & RTs_all.(condition_tags{Conditions_to_plot_S_dir(k)})'<0.2;
            idx_dir_S_late_R      =idx_dir_S_R & RTs_all.(condition_tags{Conditions_to_plot_S_dir(k)})'>=0.2;
            idx_dir_B_early_R     =idx_dir_B_R & RTs_all.(condition_tags{Conditions_to_plot_B_dir(k)})'<0.2;
            idx_dir_B_late_R      =idx_dir_B_R & RTs_all.(condition_tags{Conditions_to_plot_B_dir(k)})'>=0.2;
            
            
            idx_dir_S_early_L_30     =idx_dir_S_L_30 & RTs_all.(condition_tags{Conditions_to_plot_S_dir(k)})'<0.2 & index_return_S_L;
            idx_dir_S_late_R_30      =idx_dir_S_R_30 & RTs_all.(condition_tags{Conditions_to_plot_S_dir(k)})'>=0.2 & index_return_S_L;
            idx_dir_B_early_L_30     =idx_dir_B_L_30 & RTs_all.(condition_tags{Conditions_to_plot_B_dir(k)})'<0.2 & index_return_B_L;
            idx_dir_B_late_R_30      =idx_dir_B_R_30 & RTs_all.(condition_tags{Conditions_to_plot_B_dir(k)})'>=0.2 & index_return_B_L;
            
            N_S_early_L_30 = N_S_early_L_30 + sum(idx_dir_S_early_L_30);
            N_S_late_R_30  = N_S_late_R_30 + sum(idx_dir_S_late_R_30);
            N_B_early_L_30 = N_B_early_L_30 + sum(idx_dir_B_early_L_30);
            N_B_late_R_30  = N_B_late_R_30 + sum(idx_dir_B_late_R_30);            
            
            
%             N_S_early_L_30_with_return = N_S_early_L_30_with_return + sum(idx_dir_S_early_L_30);
%             N_S_late_R_30_with_return  = N_S_late_R_30_with_return + sum(idx_dir_S_late_R_30);
%             N_B_early_L_30_with_return = N_B_early_L_30_with_return + sum(idx_dir_B_early_L_30);
%             N_B_late_R_30_with_return  = N_B_late_R_30_with_return + sum(idx_dir_B_late_R_30);      
%             


%% Previous definition !!

% 
            N_S_early = N_S_early + sum((idx_dir_S_early_L | idx_dir_S_early_R) & (index_return_S_L | index_return_S_R));
            N_S_late  = N_S_late + sum((idx_dir_S_late_L | idx_dir_S_late_R) & (index_return_S_L | index_return_S_R));
            N_B_early = N_B_early + sum((idx_dir_B_early_L | idx_dir_B_early_R) & (index_return_B_L | index_return_B_R));
            N_B_late  = N_B_late + sum((idx_dir_B_late_L | idx_dir_B_late_R) & (index_return_B_L | index_return_B_R));
            
            
%             N_S_early = N_S_early + sum((idx_dir_S_early_L ) & (index_return_S_L ));
%             N_S_late  = N_S_late + sum((idx_dir_S_late_R) & (index_return_S_L ));
%             N_B_early = N_B_early + sum((idx_dir_B_early_L ) & (index_return_B_L ));
%             N_B_late  = N_B_late + sum((idx_dir_B_late_R ) & (index_return_B_L ));
            
            
            
            
            
            %             if ~isempty(idx_dir_S) && ~all(idx_dir_S==false)
            %                 p(2*m-1)=polar(plot_dir_S(idx_dir_S),plot_amp_S(idx_dir_S),'o');
            %             else
            %                 p(2*m-1)=polar(NaN,NaN,'o');
            %
            %             end
            %             hold on
            %             if ~isempty(idx_dir_B) && ~all(idx_dir_B==false)
            %                 p(2*m)=polar(plot_dir_B(idx_dir_B),plot_amp_B(idx_dir_B) ,'o');
            %             else
            %                 p(2*m)=polar(NaN,NaN ,'o');
            %             end
            %             set(p(2*m-1),'linestyle','none','markersize',markersize,'markerfacecolor',COLORS(m,:),'Color',COLORS(m,:));
            %             set(p(2*m),'linestyle','none','markersize',markersize,'Color',COLORS(m,:),'markeredgecolor',COLORS(m,:),'markerfacecolor','none');
            
            plot_dir_S_early_L=[plot_dir_S_early_L; plot_dir_S(idx_dir_S_early_L)];
            plot_dir_S_late_L=[plot_dir_S_late_L; plot_dir_S(idx_dir_S_late_L)];
            plot_dir_B_early_L=[plot_dir_B_early_L; plot_dir_B(idx_dir_B_early_L)];
            plot_dir_B_late_L=[plot_dir_B_late_L; plot_dir_B(idx_dir_B_late_L)];
            
            plot_dir_S_early_R=[plot_dir_S_early_R; plot_dir_S(idx_dir_S_early_R)];
            plot_dir_S_late_R=[plot_dir_S_late_R; plot_dir_S(idx_dir_S_late_R)];
            plot_dir_B_early_R=[plot_dir_B_early_R; plot_dir_B(idx_dir_B_early_R)];
            plot_dir_B_late_R=[plot_dir_B_late_R; plot_dir_B(idx_dir_B_late_R)];
            
            plot_amp_S_early_L=[plot_amp_S_early_L; plot_amp_S(idx_dir_S_early_L)];
            plot_amp_S_late_L=[plot_amp_S_late_L; plot_amp_S(idx_dir_S_late_L)];
            plot_amp_B_early_L=[plot_amp_B_early_L; plot_amp_B(idx_dir_B_early_L)];
            plot_amp_B_late_L=[plot_amp_B_late_L; plot_amp_B(idx_dir_B_late_L)];
            
            plot_amp_S_early_R=[plot_amp_S_early_R; plot_amp_S(idx_dir_S_early_R)];
            plot_amp_S_late_R=[plot_amp_S_late_R; plot_amp_S(idx_dir_S_late_R)];
            plot_amp_B_early_R=[plot_amp_B_early_R; plot_amp_B(idx_dir_B_early_R)];
            plot_amp_B_late_R=[plot_amp_B_late_R; plot_amp_B(idx_dir_B_late_R)];
            
            
            Amp_evo=[Amp_evo; plot_amp_S(idx_dir_S_early_L & index_return_S_L)];
            RT_evo=[RT_evo RTs_all.(condition_tags{Conditions_to_plot_S_dir(k)})(idx_dir_S_early_L & index_return_S_L).*1000];
            
            
            Amp_evo_30=[Amp_evo_30; plot_amp_S(idx_dir_S_early_L & index_return_S_L & idx_dir_S_early_L_30)];
            RT_evo_30=[RT_evo_30 RTs_all.(condition_tags{Conditions_to_plot_S_dir(k)})(idx_dir_S_early_L & index_return_S_L & idx_dir_S_early_L_30).*1000];
            
        end
         
            
            fraction_S_L_30_early      = num2str(round(N_S_early_L_30/N_S_early*100));
            fraction_S_R_30_late       = num2str(round(N_S_late_R_30/N_S_late*100));
            fraction_B_L_30_early      = num2str(round(N_B_early_L_30/N_B_early*100));
            fraction_B_R_30_late       = num2str(round(N_B_late_R_30/N_B_late*100));
            
%             fraction_S_L_30_early      = sum(idx_dir_S_early_L_30)/(sum(idx_dir_S_early_L)+sum(idx_dir_S_early_R));
%             fraction_S_R_30_late       = sum(idx_dir_S_late_R_30)/(sum(idx_dir_S_late_L)+sum(idx_dir_S_late_R));
%             fraction_B_L_30_early      = sum(idx_dir_B_early_L_30)/(sum(idx_dir_B_early_L)+sum(idx_dir_B_early_R));
%             fraction_B_R_30_late       = sum(idx_dir_B_late_R_30)/(sum(idx_dir_B_late_L)+sum(idx_dir_B_late_R));

        p(7)=polar(plot_dir_B_early_R,plot_amp_B_early_R,'o');
        p(3)=polar(plot_dir_B_early_L,plot_amp_B_early_L,'o');
        
        
        p(2)=polar(plot_dir_S_late_L,plot_amp_S_late_L,'o');
        p(4)=polar(plot_dir_B_late_L,plot_amp_B_late_L,'o');
        p(6)=polar(plot_dir_S_late_R,plot_amp_S_late_R,'o');
        p(8)=polar(plot_dir_B_late_R,plot_amp_B_late_R,'o');
        
        p(1)=polar(plot_dir_S_early_L,plot_amp_S_early_L,'o');
        p(5)=polar(plot_dir_S_early_R,plot_amp_S_early_R,'o');
        
        
        
        
        set(p(1),'linestyle','none','markersize',markersize,'markerfacecolor','none','Color','r');
        set(p(3),'linestyle','none','markersize',markersize,'markerfacecolor','none','Color',[0.5 0.5 0.5]);
        
        set(p(5),'linestyle','none','markersize',markersize,'markerfacecolor','none','Color','r');
        set(p(7),'linestyle','none','markersize',markersize,'markerfacecolor','none','Color',[0.5 0.5 0.5]);
        
        set(p(2),'linestyle','none','markersize',markersize,'markerfacecolor','none','Color','g');
       set(p(4),'linestyle','none','markersize',markersize,'markerfacecolor','none','Color',[0.5 0.5 0.5]);
        set(p(6),'linestyle','none','markersize',markersize,'markerfacecolor','none','Color','g');
        set(p(8),'linestyle','none','markersize',markersize,'markerfacecolor','none','Color',[0.5 0.5 0.5]);
        
        

        
%         text_N_L_B=num2str(round(sum(leftward_RT_B(RT_bins<=0.2))*100));
%         text_N_L_S=num2str(round(sum(leftward_RT_S(RT_bins<=0.2))*100));
%         text_N_R_B=num2str(round(sum(rightward_RT_B(RT_bins<=0.2))*100));
%         text_N_R_S=num2str(round(sum(rightward_RT_S(RT_bins<=0.2))*100));
        
        
        text_N_E_R_B=num2str(round(n.evoked_r.base/n.trials.base*100));
        text_N_E_R_S=num2str(round(n.evoked_r.stim/n.trials.stim*100));
        text_N_E_L_B=num2str(round(n.evoked_l.base/n.trials.base*100));        
        text_N_E_L_S=num2str(round(n.evoked_l.stim/n.trials.stim*100));
        
        text_N_R_R_B=num2str(round(n.return_r.base/n.trials.base*100));
        text_N_R_R_S=num2str(round(n.return_r.stim/n.trials.stim*100));
        text_N_R_L_B=num2str(round(n.return_l.base/n.trials.base*100));
        text_N_R_L_S=num2str(round(n.return_l.stim/n.trials.stim*100));
        
        
        text_N_R_R_B=num2str(round(n.return_r.base/n.trials.base*100));        
        %text_N_E_R_S=num2str(round(n.evoked_r.stim/n.trials.stim*100));
        
        text_N_E_P_R=num2str(round(n.evoked_and_return.stim/n.return_r.stim*100));
        text_N_R_A_E=num2str(round(n.evoked_and_return.stim/n.evoked_l.stim*100));
        
        
        text_N_E_P_R_B=num2str(round(n.evoked_and_return.base/n.return_r.base*100));
        text_N_R_A_E_B=num2str(round(n.evoked_and_return.base/n.evoked_l.base*100));
        
%         fraction_S_L_30_early      = num2str(round(sum(idx_dir_S_early_L_30)/(sum(idx_dir_S_early_L)+sum(idx_dir_S_early_R))*100));
%         fraction_S_R_30_late       = num2str(round(sum(idx_dir_S_late_R_30)/(sum(idx_dir_S_late_L)+sum(idx_dir_S_late_R))*100));
%         fraction_B_L_30_early      = num2str(round(sum(idx_dir_B_early_L_30)/(sum(idx_dir_B_early_L)+sum(idx_dir_B_early_R))*100));
%         fraction_B_R_30_late       = num2str(round(sum(idx_dir_B_late_R_30)/(sum(idx_dir_B_late_L)+sum(idx_dir_B_late_R))*100));
        
        
        text(-max_amp,-max_amp*0.9,['Contraversive: ' text_N_E_L_B '% of trials'],'Color',[0.5 0.5 0.5], 'FontSize',fontsize-2,'FontName',fonttype);     
        text(-max_amp,-max_amp*0.8,['Contra late: ' text_N_R_L_B '% of trials'],'Color',[0.3 0.3 0.3], 'FontSize',fontsize-2,'FontName',fonttype);       
        text(-max_amp,-max_amp*0.7,['Contraversive: ' text_N_E_L_S '% of trials'],'Color','r', 'FontSize',fontsize-2,'FontName',fonttype);                
        text(-max_amp,-max_amp*0.6,['Contraversive: ' text_N_R_L_S '% of trials'],'Color','g', 'FontSize',fontsize-2,'FontName',fonttype);          
        text(-max_amp,-max_amp*0.5,['Returns following evoked: ' text_N_R_A_E '%'],'Color','b', 'FontSize',fontsize-2,'FontName',fonttype);      
        text(-max_amp,-max_amp*0.4,['Returns following L control: ' text_N_R_A_E_B '%'],'Color',[0.5 0.5 0.5], 'FontSize',fontsize-2,'FontName',fonttype);      
        text(-max_amp,-max_amp*0.2,['Evo + ret 30 deg contra: ' fraction_S_L_30_early '%'],'Color','r', 'FontSize',fontsize-2,'FontName',fonttype);    
        text(-max_amp,-max_amp*0.1,['Control 30 deg contra: ' fraction_B_L_30_early '%'],'Color',[0.5 0.5 0.5], 'FontSize',fontsize-2,'FontName',fonttype);   
        
        
        text(-max_amp,max_amp*0.4,['Amp: ' num2str(nanmean(Amp_evo)) '+/-' num2str(nanstd(Amp_evo))  ' deg'],'Color','r', 'FontSize',fontsize-2,'FontName',fonttype);
        text(-max_amp,max_amp*0.5,['Lat: '  num2str(nanmean(RT_evo)) '+/-' num2str(nanstd(RT_evo))  ' ms'],'Color','r', 'FontSize',fontsize-2,'FontName',fonttype);     
        
        text(max_amp/6,max_amp*0.4,['Amp 30 deg: ' num2str(nanmean(Amp_evo_30)) '+/-' num2str(nanstd(Amp_evo_30))  ' deg'],'Color','g', 'FontSize',fontsize-2,'FontName',fonttype);
        text(max_amp/6,max_amp*0.5,['Lat 30 deg: '  num2str(nanmean(RT_evo_30)) '+/-' num2str(nanstd(RT_evo_30))  ' ms'],'Color','g', 'FontSize',fontsize-2,'FontName',fonttype);  
        
        text(max_amp/6,-max_amp*0.9,['Ipsiversive: ' text_N_E_R_B '% of trials'],'Color',[0.5 0.5 0.5], 'FontSize',fontsize-2,'FontName',fonttype);      
        text(max_amp/6,-max_amp*0.8,['Ipsi late: ' text_N_R_R_B '% of trials'],'Color',[0.3 0.3 0.3], 'FontSize',fontsize-2,'FontName',fonttype);       
        text(max_amp/6,-max_amp*0.7,['Ipsiversive: ' text_N_E_R_S '% of trials'],'Color','r', 'FontSize',fontsize-2,'FontName',fonttype);     
        text(max_amp/6,-max_amp*0.6,['Ipsiversive: ' text_N_R_R_S '% of trials'],'Color','g', 'FontSize',fontsize-2,'FontName',fonttype);
        text(max_amp/6,-max_amp*0.5,['Return preceded with evoked: ' text_N_E_P_R '%'],'Color','g', 'FontSize',fontsize-2,'FontName',fonttype);   
        text(max_amp/6,-max_amp*0.4,['Return preceded with L control: ' text_N_E_P_R_B '%'],'Color','b', 'FontSize',fontsize-2,'FontName',fonttype);   
        %text(max_amp/6,-max_amp*0.3,['Return in baseline: ' text_N_R_R_B '%'],'Color',[0.5 0.5 0.5], 'FontSize',fontsize-2,'FontName',fonttype);  
        text(max_amp/6,-max_amp*0.2,['Ret + Evo 30 deg ipsi: ' fraction_S_R_30_late '%'],'Color','r', 'FontSize',fontsize-2,'FontName',fonttype);    
        text(max_amp/6,-max_amp*0.1,['Control 30 deg ipsi: ' fraction_B_R_30_late '%'],'Color',[0.5 0.5 0.5], 'FontSize',fontsize-2,'FontName',fonttype);    
        
        
        
            Amp_evo=[Amp_evo; plot_amp_S(idx_dir_S_early_L)];
            RT_evo=[RT_evo RTs_all.(condition_tags{Conditions_to_plot_S_dir(k)})(idx_dir_S_early_L)];
    else
        if ~isempty(idx_dir_S) && ~all(idx_dir_S==false)
            p(1)=polar(plot_dir_S(idx_dir_S),plot_amp_S(idx_dir_S),'o');
        else
            p(1)=polar(NaN,NaN,'o');
            hold on
        end
        
        if ~isempty(idx_dir_B) && ~all(idx_dir_B==false)
            p(2)=polar(plot_dir_B(idx_dir_B),plot_amp_B(idx_dir_B) ,'o');
        else
            p(2)=polar(NaN,NaN ,'o');
        end
        set(p(1),'linestyle','none','markersize',markersize,'markerfacecolor','r','color','r');
        set(p(2),'linestyle','none','markersize',markersize,'markerfacecolor','g','color','g');
    end
    
    hHiddenText = findall(gca,'type','text');
    for radius_counter=1:numel(polar_ticks)
       hObjr = findall(hHiddenText,'string',['  ', num2str(polar_ticks(radius_counter))]);
       hObja = findall(hHiddenText,'string',num2str(polar_ticks(radius_counter)));
       oldpos_r=get(hObjr,'Position');
        if ~isempty(oldpos_r)
            %set(hObjr,'Position',[-1*max_amp/4, oldpos_r(2)-oldpos_r(2)./2.5, 0],'FontSize',fontsize-2)
            set(hObjr,'FontSize',fontsize-2)
        end
        set(hObja,'FontSize',fontsize-2)
    end
    
    for radius_counter=1:numel(polar_ticks_deg)
        %set(findall(hHiddenText,'string',num2str(polar_ticks_deg(radius_counter))),'String', [num2str(polar_ticks_deg(radius_counter)) ' deg']); 
        set(findall(hHiddenText,'string',['  ', num2str(polar_ticks_deg(radius_counter))]),'String', [num2str(polar_ticks_deg(radius_counter)) ' deg']); 
    end
    
    box on
    set(gca,'fontsize',fontsize-3,'FontName',fonttype,'LineWidth',linewidth);
    title(Labels_polar{k},'interpreter','none','FontSize',fontsize+2)
end
if trial_type~=1
    %h2=legend(p([1,6,2,3,8,4]), 'evoked L in','return R out','stim other','base L in','base R out','base other');
else
        h1=legend(h([1,2,3]), 'during stimulation','after stimulation','no stimulation');
        h2=legend(p([1,2,3]), 'during stimulation','after stimulation','no stimulation');
    %h2=legend(p([1,3]), 'during stimulation','no stimulation');
    
    %h2=legend(p(2*(skipped_amps+1)-1:2*(numel(amplitudes)-3)),Whole_legend(skipped_amps*2+1:end-4));    
    %h2=legend(p(2*(skipped_amps+1)-1:2*(numel(amplitudes)-3)),Whole_legend(skipped_amps*2+1:end-4));
end
set(h1,'Position', leg1Position,'Units', 'normalized','FontName',fonttype);
set(h2,'Position', leg2Position,'Units', 'normalized','FontName',fonttype);

%mtit(evoked_saccades_figure, current_title, 'fontsize', 12, 'color', [0 0 0], 'xoff', -0.05, 'yoff', 0.03,  'FontSize', fontsize + 5,'FontName',fonttype);
end

function par=mtit(varargin)
%MTIT		creates a major title in a figure with many axes
%
%		MTIT
%		- creates a major title above all
%		  axes in a figure
%		- preserves the stack order of
%		  the axis handles
%
%SYNTAX
%-------------------------------------------------------------------------------
%		P = MTIT(TXT,[OPT1,...,OPTn])
%		P = MTIT(FH,TXT,[OPT1,...,OPTn])
%
%INPUT
%-------------------------------------------------------------------------------
%    FH	:	a valid figure handle		[def: gcf]
%   TXT	:	title string
%
% OPT	:	argument
% -------------------------------------------
%  xoff	:	+/- displacement along X axis
%  yoff	:	+/- displacement along Y axis
%  zoff	:	+/- displacement along Z axis
%
%		title modifier pair(s)
% -------------------------------------------
%   TPx	:	TVx
%		see: get(text) for possible
%		     parameters/values
%
%OUTPUT
%-------------------------------------------------------------------------------
% par	:	parameter structure
%  .pos :	position of surrounding axis
%   .oh	:	handle of last used axis
%   .ah :	handle of invisible surrounding axis
%   .th :	handle of main title
%
%EXAMPLE
%-------------------------------------------------------------------------------
%	subplot(2,3,[1 3]);		title('PLOT 1');
%	subplot(2,3,4); 		title('PLOT 2');
%	subplot(2,3,5); 		title('PLOT 3');
%	axes('units','inches',...
%	     'color',[0 1 .5],...
%	     'position',[.5 .5 2 2]);	title('PLOT 41');
%	axes('units','inches',...
%	     'color',[0 .5 1],...
%	     'position',[3.5 .5 2 2]);	title('PLOT 42');
%	shg;
%	p=mtit('the BIG title',...
%	     'fontsize',14,'color',[1 0 0],...
%	     'xoff',-.1,'yoff',.025);
% % refine title using its handle <p.th>
%	set(p.th,'edgecolor',.5*[1 1 1]);

% created:
%	us	24-Feb-2003		/ R13
% modified:
%	us	24-Feb-2003		/ CSSM
%	us	06-Apr-2003		/ TMW
%	us	13-Nov-2009 17:38:17

defunit='normalized';
if	nargout
    par=[];
end

% check input
if	nargin < 1
    help(mfilename);
    return;
end
if	isempty(get(0,'currentfigure'))
    disp('MTIT> no figure');
    return;
end

vl=true(size(varargin));
if	ischar(varargin{1})
    vl(1)=false;
    figh=gcf;
    txt=varargin{1};
elseif	any(ishandle(varargin{1}(:)))		&&...
        ischar(varargin{2})
    vl(1:2)=false;
    figh=varargin{1};
    txt=varargin{2};
else
    error('MTIT> invalid input');
end
vin=varargin(vl);
[off,vout]=get_off(vin{:});

% find surrounding box
ah=findall(figh,'type','axes');
if	isempty(ah)
    disp('MTIT> no axis');
    return;
end
oah=ah(1);

ou=get(ah,'units');
set(ah,'units',defunit);
ap=get(ah,'position');
if	iscell(ap)
    ap=cell2mat(get(ah,'position'));
end
ap=[	min(ap(:,1)),max(ap(:,1)+ap(:,3)),...
    min(ap(:,2)),max(ap(:,2)+ap(:,4))];
ap=[	ap(1),ap(3),...
    ap(2)-ap(1),ap(4)-ap(3)];

% create axis...
xh=axes('position',ap);
% ...and title
th=title(txt,'interpreter', 'none',vout{:});
tp=get(th,'position');
set(th,'position',tp+off);
set(xh,'visible','off','hittest','on');
set(th,'visible','on');

% reset original units
ix=find(~strcmpi(ou,defunit));
if	~isempty(ix)
    for	i=ix(:).'
        set(ah(i),'units',ou{i});
    end
end

% ...and axis' order
uistack(xh,'bottom');
axes(oah);				%#ok

if	nargout
    par.pos=ap;
    par.oh=oah;
    par.ah=xh;
    par.th=th;
end

    function	[off,vout]=get_off(varargin)
        
        % search for pairs <.off>/<value>
        
        off=zeros(1,3);
        io=0;
        for	mode={'xoff','yoff','zoff'};
            ix=strcmpi(varargin,mode);
            if	any(ix)
                io=io+1;
                yx=find(ix);
                ix(yx+1)=1;
                off(1,io)=varargin{yx(end)+1};
                varargin=varargin(xor(ix,1));
            end
        end
        vout=varargin;
    end
end

 function STATE=load_state


STATE.INI_TRI           = 1;  % initialize trial
STATE.FIX_ACQ           = 2;  % fixation acquisition
STATE.FIX_HOL           = 3;  % fixation hold
STATE.TAR_ACQ           = 4;  % target acquisition
STATE.TAR_HOL           = 5;  % target hold
STATE.CUE_ON            = 6;  % cue on
STATE.MEM_PER           = 7;  % memory period
STATE.DEL_PER           = 8;  % delay period
STATE.TAR_ACQ_INV       = 9;  % target acquisition invisible
STATE.TAR_HOL_INV       = 10; % target hold invisible
STATE.ABORT             = 19;
STATE.SUCCESS           = 20;
STATE.REWARD            = 21;
STATE.ITI               = 50;
STATE.CLOSE             = 99;
STATE.SUCCESS_ABORT     = -1;

STATE.ALL               =[STATE.INI_TRI STATE.FIX_ACQ  STATE.FIX_HOL  STATE.CUE_ON  STATE.MEM_PER  STATE.DEL_PER STATE.TAR_ACQ_INV  STATE.TAR_HOL_INV...
    STATE.TAR_ACQ  STATE.TAR_HOL STATE.SUCCESS_ABORT  STATE.REWARD  STATE.ITI];
STATE.ALL_NAMES         ={'Trial Initiation', 'Fixation Acquisition', 'Fixation Hold', 'Cue', 'Memory', 'Delay', 'Target Acquisition (inv)',...
    'Target Hold (inv)','Target Acquisition', 'Target Hold', 'Success', 'Reward', 'ITI'};
STATE.ALL_NAMES_Short         ={'INI', 'FIX_ACQ', 'FIX_HOL', 'CUE', 'MEM', 'DEL', 'TAR_ACQ_INV',...
    'TAR_HOL_INV','TAR_ACQ', 'TAR_HOLd', 'SUC', 'REW', 'ITI'};

STATE.aquisition=[STATE.FIX_ACQ STATE.TAR_ACQ STATE.TAR_ACQ_INV];
end