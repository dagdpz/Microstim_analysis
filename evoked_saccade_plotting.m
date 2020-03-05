function evoked_saccade_plotting
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
            y_lim              =0.3;
        elseif strcmp(first_3,'dir')
            trial_type         =2;
            y_lim              =1;
        else strcmp(first_3,'mem')
            trial_type         =3;
            y_lim              =1;
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
            n_trials                = current.n_trials(num_stim_window) ;            
            n_evoked                = current.n_evoked(num_stim_window);
            RTs_histograms_L        = current.RTs_histograms_L(num_stim_window) ;
            RTs_histograms_R        = current.RTs_histograms_R(num_stim_window) ;
            Dir_all                 = current.Dir_all{num_stim_window} ;
            Amp_all                 = current.Amp_all{num_stim_window} ;
            RTs_all                 = current.RTs_all{num_stim_window} ;
            
            
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
                plot_saccade_prob(RT_bins,n_trials,RTs_histograms_L,RTs_histograms_R,Dir_all,Amp_all,RTs_all,COLORS,current_title,trial_type,y_lim,n_evoked,n_trials,skipped_amps,Evo_amp_TH(1),far_targets)
                
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
            plot_saccade_prob(RT_bins,current.summed_n_trials,current.summed_RTs_histograms_L,current.summed_RTs_histograms_R,current.summed_dir_all,current.summed_amp_all,current.summed_RTs_all,...
                COLORS,current_title,trial_type,y_lim,current.summed_n_evoked,current.summed_n_trials,skipped_amps,Evo_amp_TH(1),far_targets)
            
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

function plot_saccade_prob(RT_bins,n_trials,RTs_histograms_L,RTs_histograms_R,Dir_all, Amp_all,RTs_all,COLORS,current_title,trial_type,y_lim,n_e,n_tot,skipped_amps,min_amp,far_targets)
fontsize=20;
linewidth=5;
markersize=8;
amplitudes=[0.5,1,2,4,6,8,10,15,100];
condition_tags=fieldnames(n_trials);
fonttype='Arial';
set(0,'DefaultAxesFontName', fonttype)
n_amps=numel(amplitudes)-1;

% x_lim=max(RT_bins);
if trial_type==1 %fix
    x_lim=0.6;
    Conditions_to_plot_S=[1 1];
    Conditions_to_plot_B=[2 2];
    Conditions_to_plot_S_dir=1;
    Conditions_to_plot_B_dir=2;
    
    ks_to_plot_L=1;
    ks_to_plot_R=2;
    ks_to_title=[1,2];
    ks_to_lable=2;
    leg1Position = [0.9 0.16 0.06 0.3]; %x_pos y_pos width height % all in proportion to figure size
    leg2Position = [0.9 0.16 0.06 0.3]; %x_pos y_pos width height % all in proportion to figure size
    
    Labels_RT_comparison={'Leftward (evoked) saccades ','Rightward (return) saccades '};
    Labels_polar={{'saccade vector endpoints'}};
    
    dir_subplot_positions={[2,3,5,6]};
    RT_plot_positions={1,4};
    
    n_subplots_y=2;
    n_subplots_x=3;
    polar_ticks=[0,30,60,90,120,150,180,210,240,270,300,330,2,4,6,8,10];
    max_amp=10;
else %dir,mem
    x_lim=0.4;
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
        leftward_RT_B=zeros(size(RTs_histograms_L.(condition_tags{Conditions_to_plot_B(k)})(1,:)));        
        leftward_RT_S=zeros(size(RTs_histograms_L.(condition_tags{Conditions_to_plot_S(k)})(1,:))); 
        for m=(skipped_amps+1):n_amps
            hold on
            %plot(RT_bins,RTs_histograms_L.(condition_tags{Conditions_to_plot_S(k)})(m,:),'-','Color',COLORS(m,:),'LineWidth',linewidth);
            %plot(RT_bins,RTs_histograms_L.(condition_tags{Conditions_to_plot_B(k)})(m,:),'--','Color',COLORS(m,:),'LineWidth',linewidth);
            leftward_RT_B=leftward_RT_B+RTs_histograms_L.(condition_tags{Conditions_to_plot_B(k)})(m,:);
            leftward_RT_S=leftward_RT_S+RTs_histograms_L.(condition_tags{Conditions_to_plot_S(k)})(m,:);
        end
        plot(RT_bins(RT_bins<0.2),leftward_RT_B(RT_bins<0.2),'--','Color','r','LineWidth',linewidth);
        plot(RT_bins(RT_bins<0.2),leftward_RT_S(RT_bins<0.2),'-','Color','r','LineWidth',linewidth);
        plot(RT_bins(RT_bins>=0.2),leftward_RT_B(RT_bins>=0.2),'--','Color','b','LineWidth',linewidth);
        plot(RT_bins(RT_bins>=0.2),leftward_RT_S(RT_bins>=0.2),'-','Color','b','LineWidth',linewidth);
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
        plot(RT_bins(RT_bins<0.2),rightward_RT_B(RT_bins<0.2),'--','Color','r','LineWidth',linewidth);
        plot(RT_bins(RT_bins<0.2),rightward_RT_S(RT_bins<0.2),'-','Color','r','LineWidth',linewidth);
        plot(RT_bins(RT_bins>=0.2),rightward_RT_B(RT_bins>=0.2),'--','Color','b','LineWidth',linewidth);
        plot(RT_bins(RT_bins>=0.2),rightward_RT_S(RT_bins>=0.2),'-','Color','b','LineWidth',linewidth);
    end
    
    if trial_type==1
        ylabel('Saccades / Trial','FontSize',fontsize,'FontName',fonttype)
    else
        if k==1
            ylabel({'Saccades / Trial';' (CHOICE)'},'FontSize',fontsize,'FontName',fonttype)
        elseif k==5
            ylabel({'Saccades / Trial ';'(INSTRUCTED)'},'FontSize',fontsize,'FontName',fonttype)
        end
    end
    if ismember(k,ks_to_title)
        title(Labels_RT_comparison{k},'interpreter','none','FontSize',fontsize+2,'FontName',fonttype)
    end
    if ismember(k,ks_to_lable)
        xlabel('Time after stimulation [s]','fontsize',fontsize,'FontName',fonttype)
    end
    set(gca,'xlim',[0 x_lim],'ylim',[0, y_lim],'FontSize',fontsize-2,'FontName',fonttype,'LineWidth',linewidth-2)
    if ismember(k,ks_to_plot_R)
        if trial_type==1
            xlabel('Time after stimulation [s]','fontsize',fontsize,'FontName',fonttype)
        else
            set(gca,'ytick',[],'FontName',fonttype);
            set(gca,'ycolor',[1 1 1],'FontName',fonttype)
        end
    end
    if ismember(k,ks_to_plot_L)
        text(0.01,y_lim-(y_lim/10),'# B/S =', 'FontSize',fontsize-2,'FontName',fonttype)
        text(0.25,y_lim-(y_lim/10),[num2str(n_trials.(condition_tags{Conditions_to_plot_B(k)})) '/' num2str(n_trials.(condition_tags{Conditions_to_plot_S(k)}))], 'FontSize',fontsize-2,'FontName',fonttype);
    end
    NE      =n_e.(condition_tags{Conditions_to_plot_S(k)});
    Ntot    =n_tot.(condition_tags{Conditions_to_plot_S(k)});
%     NnE     =n_e.(condition_tags{Conditions_to_plot_B(k)});
%     Nntot   =n_tot.(condition_tags{Conditions_to_plot_B(k)});
    p_evoked_S      =round(NE*100./Ntot)./100;
    
    if ismember(k,ks_to_plot_L)
        text(0.01,y_lim-(2*y_lim/10),'p_e_v_o_k_e_d =', 'FontSize',fontsize-2,'FontName',fonttype)
        text(0.25,y_lim-(2*y_lim/10),num2str(p_evoked_S), 'FontSize',fontsize-2,'FontName',fonttype)
    end
    line([0.2,0.2]',[0,y_lim]', 'LineWidth',linewidth, 'Color','k');
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
        plot_dir_S_early=[];
        plot_dir_S_late=[];
        plot_dir_B_early=[];
        plot_dir_B_late=[];
                
        plot_amp_S_early=[];
        plot_amp_S_late=[];
        plot_amp_B_early=[];
        plot_amp_B_late=[];
        
        for m=(skipped_amps+1):numel(amplitudes)-1
            idx_dir_S=plot_amp_S>amplitudes(m) & plot_amp_S<amplitudes(m+1);
            idx_dir_B=plot_amp_B>amplitudes(m) & plot_amp_B<amplitudes(m+1);
            
            idx_dir_S_early     =idx_dir_S&RTs_all.(condition_tags{Conditions_to_plot_S_dir(k)})'<0.2;
            idx_dir_S_late      =idx_dir_S&RTs_all.(condition_tags{Conditions_to_plot_S_dir(k)})'>=0.2;
            idx_dir_B_early     =idx_dir_B&RTs_all.(condition_tags{Conditions_to_plot_B_dir(k)})'<0.2;
            idx_dir_B_late      =idx_dir_B&RTs_all.(condition_tags{Conditions_to_plot_B_dir(k)})'>=0.2;
            
            
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

        plot_dir_S_early=[plot_dir_S_early; plot_dir_S(idx_dir_S_early)];
        plot_dir_S_late=[plot_dir_S_late; plot_dir_S(idx_dir_S_late)];
        plot_dir_B_early=[plot_dir_B_early; plot_dir_B(idx_dir_B_early)];
        plot_dir_B_late=[plot_dir_B_late; plot_dir_B(idx_dir_B_late)];
        
        plot_amp_S_early=[plot_amp_S_early; plot_amp_S(idx_dir_S_early)];
        plot_amp_S_late=[plot_amp_S_late; plot_amp_S(idx_dir_S_late)];
        plot_amp_B_early=[plot_amp_B_early; plot_amp_B(idx_dir_B_early)];
        plot_amp_B_late=[plot_amp_B_late; plot_amp_B(idx_dir_B_late)];
        end
        p(1)=polar(plot_dir_S_early,plot_amp_S_early,'o');
        p(2)=polar(plot_dir_S_late,plot_amp_S_late,'o');
        p(3)=polar(plot_dir_B_early,plot_amp_B_early,'o');
        p(4)=polar(plot_dir_B_late,plot_amp_B_late,'o');
        set(p(1),'linestyle','none','markersize',markersize,'markerfacecolor','r','Color','r');
        set(p(2),'linestyle','none','markersize',markersize,'markerfacecolor','b','Color','b');
        set(p(3),'linestyle','none','markersize',markersize,'markerfacecolor','w','Color','r');
        set(p(4),'linestyle','none','markersize',markersize,'markerfacecolor','w','Color','b');
        
        
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
            set(hObjr,'Position',[-1*max_amp/4, oldpos_r(2)-oldpos_r(2)./2.5, 0],'FontSize',fontsize-2)
        end
        set(hObja,'FontSize',fontsize-2)
    end
    
    set(gca,'fontsize',fontsize-3,'FontName',fonttype,'LineWidth',linewidth);
    title(Labels_polar{k},'interpreter','none','FontSize',fontsize+2)
end
if trial_type~=1
    h2=legend(p, 'Stim','Base');
else
    %h2=legend(p(2*(skipped_amps+1)-1:2*(numel(amplitudes)-3)),Whole_legend(skipped_amps*2+1:end-4));
    
    %h2=legend(p(2*(skipped_amps+1)-1:2*(numel(amplitudes)-3)),Whole_legend(skipped_amps*2+1:end-4));
end
%set(h2,'Position', leg2Position,'Units', 'normalized','FontName',fonttype);

mtit(evoked_saccades_figure, current_title, 'fontsize', 12, 'color', [0 0 0], 'xoff', -0.05, 'yoff', 0.03,  'FontSize', fontsize + 5,'FontName',fonttype);
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