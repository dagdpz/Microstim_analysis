function scatterplots_microstim(monkey,X_par,Y_par,Selection,reduce_to_residuals)
COLOR=jet(16);
updated_xls=strcat(monkey, '_updated_parameters');
parameters_to_look_at={'Evoked_P', 'Delayed_P', 'Evoked_mean_RT', 'Evoked_mean_amp', 'N_stim', 'N_evoked'};
Conditions={'_L_CH','_R_CH','_L_IN','_R_IN'};
windows={'_early','_late'};
fonttype='Arial';
set(0,'DefaultAxesFontName', fonttype);


switch X_par
    case 'current'
        X_condition='current';
        %if  strcmp(Selection,'fix') || strcmp(Selection,'freegaze')
            x_lim=[25 325];
        %else
        %    x_lim=[50 300];
        %end
        X_label=strcat('Current strength [','\mu','A]');
    case 'electrode depth'
        X_condition='depth';
        switch monkey
            case 'Curius'
                x_lim=[43 47];
            case 'Linus'
                x_lim=[45.5 49];
        end
        X_label='Electrode depth [mm]';
    case 'probability of evoking saccades'
        X_condition='Evoked_P';
        x_lim=[-0.05 1.05];
        X_label='p_e_v_o_k_e_d';
    case 'probability of delaying saccades'
        X_condition='Delayed_P';
        x_lim=[-0.05 1.05];
        X_label='p_d_e_l_a_y_e_d';
    case 'Amplitudes of evoked saccades'
        X_condition='Evoked_mean_amp';
        x_lim=[0 10];
        X_label={'Mean Evoked saccade'; 'amplitude [deg]'};
end

switch Y_par
    case 'Current'
        Y_condition='current';
        y_lim=[0 350];
        Y_label=strcat('Current strength [','\mu','A]');
    case 'Electrode depth'
        Y_condition='depth';
        y_lim=[43 48];
        Y_label='Electrode depth [mm]';
    case 'Contraversive target selection'
        Y_condition='bias';
        if  strcmp(Selection,'shifting-')
            y_lim=[-0.05 0.4];
        elseif strcmp(Selection,'shifting+')
            y_lim=[-0.05 1.3];
        else
            y_lim=[-0.5 0.5];
        end
        Y_label='Contraversive target selection';
    case 'Contraversive target selection difference'
        Y_condition='bias_difference';
        if  strcmp(Selection,'shifting-')            
            y_lim=[-0.4 0.4];
        elseif strcmp(Selection,'shifting+')
            y_lim=[-0.3 1];
        else
            y_lim=[-0.5 0.5];            
        end
        Y_label={'Contraversive target selection'; 'difference (stim-baseline)'};   
    case 'Probability of evoking saccades'
        Y_condition='Evoked_P';
        y_lim=[-0.05 1.4];
        Y_label='p_e_v_o_k_e_d';
    case 'Probability of delaying saccades'
        Y_condition='Delayed_P';
        y_lim=[-0.05 1.4];
        Y_label='p_d_e_l_a_y_e_d';
    case 'Amplitudes of evoked saccades'
        Y_condition='Evoked_mean_amp';
        y_lim=[0 12];
        if strcmp(Selection,'freegaze')
            y_lim=[0 30];
        end
        Y_label={'Mean evoked saccade'; 'amplitude [deg]'};
    case 'Latencies of evoked saccades'
        Y_condition='Evoked_mean_RT';
        y_lim=[0.05 0.15];
        Y_label={'Mean evoked saccade'; 'latency [s]'};
end

switch Selection %% tricky though... ??!!
    case 'all'
        batch_spec=', all runs';
        usedcolor=nanmean(COLOR);
        selection='_all';
    case 'shifting+'
        batch_spec=', late windows';
        usedcolor=nanmean(COLOR(7:end,:));
        selection='_late';
    case 'shifting-'
        batch_spec=', early windows';
        usedcolor=nanmean(COLOR(1:6,:));
        selection='_early';
    case 'fix'
        batch_spec=', fixation runs';
        usedcolor='r';
        selection='_fix';
    case 'freegaze'
        %Aiminput{1}={'map_evo_fre'}; %!!!!!
        batch_spec=', freegaze runs';
        usedcolor='k';
        selection='_fix';
end

if ismember(X_par,{'probability of evoking saccades'}) && ismember(Y_par,{'Probability of delaying saccades'}) && strcmp(selection,'_early');
    Y_condition=[Y_condition '_late'];  %% particular condition of comparing delay in late with evoked in early
elseif ismember(Y_par,{'Probability of evoking saccades','Probability of delaying saccades','Amplitudes of evoked saccades','Latencies of evoked saccades','Contraversive target selection','Contraversive target selection difference'})
    Y_condition=[Y_condition selection];
end
if ismember(X_par,{'probability of evoking saccades'}) && ismember(Y_par,{'Probability of delaying saccades'}) && strcmp(selection,'_late');
    X_condition=[X_condition '_early'];  %% particular condition of comparing delay in late with evoked in early
elseif ismember(X_par,{'probability of evoking saccades','probability of delaying saccades','Amplitudes of evoked saccades'})
    X_condition=[X_condition selection];
end

batch_title=[Y_par ' versus ' X_par batch_spec];

[mastertable_excel_num,~,mastertable_excel_raw]=xlsread(updated_xls,'mastertable');
idx.session             =DAG_find_column_index(mastertable_excel_raw,'Session');
idx.run                 =DAG_find_column_index(mastertable_excel_raw,'Run');
idx.current             =DAG_find_column_index(mastertable_excel_raw,'Current_strength');
idx.depth               =DAG_find_column_index(mastertable_excel_raw,'Electrode_depth');
idx.bias_late           =DAG_find_column_index(mastertable_excel_raw,'bias_late');
idx.bias_early          =DAG_find_column_index(mastertable_excel_raw,'bias_early');
idx.bias_base           =DAG_find_column_index(mastertable_excel_raw,'bias_base');


comparison.session                  =mastertable_excel_num(:,idx.session);
comparison.run                      =mastertable_excel_num(:,idx.run);
comparison.current                  =mastertable_excel_num(:,idx.current);
comparison.depth                    =mastertable_excel_num(:,idx.depth);
comparison.bias_early               =mastertable_excel_num(:,idx.bias_early);
comparison.bias_late                =mastertable_excel_num(:,idx.bias_late);
comparison.bias_difference_early    =mastertable_excel_num(:,idx.bias_early) - mastertable_excel_num(:,idx.bias_base);
comparison.bias_difference_late     =mastertable_excel_num(:,idx.bias_late)  - mastertable_excel_num(:,idx.bias_base);

for par=1:numel(parameters_to_look_at)
    current_field=[parameters_to_look_at{par} '_fix'];
    idx.(current_field)=DAG_find_column_index(mastertable_excel_raw,current_field);
    if isempty(idx.(current_field))
        comparison.(current_field)=NaN(size(mastertable_excel_num,1),1);
    else
        comparison.(current_field)=mastertable_excel_num(:,idx.(current_field));
    end
    for win = 1:numel(windows)
        current_N_stim=0;
        current_N_evoked=0; 
        for con = 1:numel(Conditions)
            current_field=[parameters_to_look_at{par} windows{win} Conditions{con}];
            idx.(current_field)=DAG_find_column_index(mastertable_excel_raw,current_field);
            if isempty(idx.(current_field))
                comparison.(current_field)=NaN(size(mastertable_excel_num,1),1);
            else
                comparison.(current_field)=mastertable_excel_num(:,idx.(current_field));
            end
            if any(strfind(current_field, 'N_stim')) || any(strfind(current_field, 'N_evoked'))
                current_N_stim=current_N_stim + comparison.(current_field);
                current_N_evoked=current_N_evoked + comparison.(current_field);
            end
        end
        if any(strfind(current_field, 'N_stim'))
            comparison.(['N_stim' windows{win}])=current_N_stim;
        elseif any(strfind(current_field, 'N_evoked'))
            comparison.(['N_evoked' windows{win}])=current_N_evoked;
        end
        
    end
end

subplotmatrix=[1,1];
y_val={Y_condition};
x_val={X_condition};

markersize=25;
linewidth=4;
fontsize=28;

subplot_title{1}=batch_title;
N_stim={['N_stim' selection]};
N_evoked={['N_evoked' selection]};

Ys_for_four_plots={'Probability of evoking saccades','Probability of delaying saccades','Amplitudes of evoked saccades','Latencies of evoked saccades'};
Xs_for_four_plots={'probability of evoking saccades','probability of delaying saccades','Amplitudes of evoked saccades'};

if strcmp(Selection,'fix') || strcmp(Selection,'freegaze')    
    N_stim={'N_stim_fix'};
    N_evoked={'N_evoked_fix'};
else
    if ismember(Y_par,Ys_for_four_plots)
        subplotmatrix=[2,2];
        linewidth=3;
        fontsize=20;
        markersize=15;
        y_val=strcat(repmat({Y_condition},4,1), Conditions');
        N_stim=strcat(repmat({'N_stim'},4,1), selection, Conditions');
        N_evoked=strcat(repmat({'N_evoked'},4,1), selection, Conditions');
        subplot_title={'Choice Left','Choice Right','Instructed Left','Instructed Right'}; % Conditions={'_L_CH','_R_CH','_L_IN','_R_IN'};
        %ks_for_ylabel=[1,3];
        %ks_for_xlabel=[3,4];
        if ~ ismember(X_par,Xs_for_four_plots)
            x_val={X_condition,X_condition,X_condition,X_condition};
        end
    end
    if ismember(X_par,Xs_for_four_plots)
        subplotmatrix=[2,2];
        linewidth=3;
        fontsize=20;
        markersize=15;
        x_val=strcat(repmat({X_condition},4,1), Conditions');
        N_stim=strcat(repmat({'N_stim'},4,1), selection, Conditions');
        N_evoked=strcat(repmat({'N_evoked'},4,1), selection, Conditions');
        subplot_title={'Choice Left','Choice Right','Instructed Left','Instructed Right'}; % Conditions={'_L_CH','_R_CH','_L_IN','_R_IN'};
        %ks_for_ylabel=[1,3];
        %ks_for_xlabel=[3,4];
        if ~ ismember(Y_par,Ys_for_four_plots)
            y_val={Y_condition,Y_condition,Y_condition,Y_condition};
        end
    end
end


current_figure=figure('units','normalized','outerposition',[0 0 1 1]);

    all_sites=[comparison.session,comparison.depth];
    site_index=zeros(size(all_sites,1),1);
    unique_sites=unique(all_sites,'rows');
    for s=1:size(unique_sites,1)
        site_index(all_sites(:,1)== unique_sites(s,1) & all_sites(:,2)==unique_sites(s,2))=s;
    end
for k=1:max(numel(x_val),numel(y_val))
    if reduce_to_residuals
       subplot_title{k}=[subplot_title{k} ' residuals'];
       y_lim=y_lim-nanmean(y_lim);
       for s=1:size(unique_sites,1)
       comparison.(y_val{k})(site_index==s)=comparison.(y_val{k})(site_index==s)-nanmean(comparison.(y_val{k})(site_index==s)); 
       end
    end
    
    y_max=max(y_lim);
    x_min=min(x_lim);
    y_step=(max(y_lim)-min(y_lim))./20;
    x_step=(max(x_lim)-min(x_lim))./10;
    subplot(subplotmatrix(1),subplotmatrix(2),k)
    hs=scatter(comparison.(x_val{k}),comparison.(y_val{k}));
    
    title(strcat(subplot_title{k}), 'fontsize',fontsize+2, 'interpreter', 'none','FontName',fonttype,'FontWeight','bold')
    %if ismember(k,ks_for_ylabel)
    ylabel(Y_label,'fontsize',fontsize,'FontName',fonttype)
    %end
    %if ismember(k,ks_for_xlabel)
    xlabel(X_label,'fontsize',fontsize,'FontName',fonttype)
    %end
    
    hp=get(hs(1),'children'); % handle for plot inside scatterplot axes
    set(hp,'Marker','o','MarkerSize',markersize,'MarkerEdgeColor','k','LineWidth',linewidth);%,'MarkerFaceColor',COLOR(k+1,:));
    set(hs,'MarkerFaceColor',usedcolor,'Marker','o','MarkerEdgeColor','k','LineWidth',linewidth);
    rl=refline;
    set(rl,'Color','k','LineWidth',linewidth + 1)
    hold on
    idx_comparison=~isnan(comparison.(x_val{k})) & ~isnan(comparison.(y_val{k}));
    if any(idx_comparison)
        [RHO,PVAL] = corr(comparison.(x_val{k})(idx_comparison),comparison.(y_val{k})(idx_comparison),'type','Spearman');
    else
        [RHO,PVAL] = deal(NaN);
    end
    text(x_min+x_step/10, y_max-1.5*y_step,     ['R = ',            num2str(round(RHO*1000)/1000)]                      ,'fontsize',fontsize-2,'FontName',fonttype);
    text(x_min+x_step/10, y_max-3.5*y_step,     ['P = ',            num2str(round(PVAL*1000)/1000)]                     ,'fontsize',fontsize-2,'FontName',fonttype);
    text(x_min+3*x_step, y_max-1.5*y_step,      ['Mean of runs = ', num2str(round(nanmean(comparison.(y_val{k}))*100)/100)],'fontsize',fontsize-2,'FontName',fonttype);
    text(x_min+3*x_step, y_max-3.5*y_step,      ['SEM of runs = ',  num2str(round(sterr(comparison.(y_val{k}))*100)/100)],'fontsize',fontsize-2,'FontName',fonttype);
    text(x_min+7*x_step, y_max-1.5*y_step,      ['N_s_t_i_m = ',    num2str(nansum([comparison.(N_stim{k})]))]          ,'fontsize',fontsize-2,'FontName',fonttype);
    text(x_min+7*x_step, y_max-3.5*y_step,      ['N_e_v_o_  = ',     num2str(nansum([comparison.(N_evoked{k})]))]        ,'fontsize',fontsize-2,'FontName',fonttype);
    
    set(gca,'ylim',y_lim, 'xlim',x_lim,'fontsize',fontsize-2,'FontName',fonttype,'LineWidth',linewidth)
    box on
end
if (ismember(Y_par,Ys_for_four_plots) || ismember(X_par,Xs_for_four_plots)) && ~strcmp(Selection,'fix') && ~strcmp(Selection,'freegaze')  
    mtit(current_figure, batch_title, 'fontsize', fontsize+4,'FontWeight','bold', 'color', [0 0 0], 'xoff', -0.05, 'yoff', 0.045)
end



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

 
