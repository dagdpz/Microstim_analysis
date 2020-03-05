function batch_out = per_session_significance(batch_out)
global GLO

n_windows=numel(GLO.windows_to_use)+2;
Hits.L_CH=zeros(numel(batch_out.num_hits),n_windows);
Hits.L_IN=zeros(numel(batch_out.num_hits),n_windows);
Hits.R_CH=zeros(numel(batch_out.num_hits),n_windows);
Hits.R_IN=zeros(numel(batch_out.num_hits),n_windows);
Total.L_CH=zeros(numel(batch_out.num_hits),n_windows);
Total.L_IN=zeros(numel(batch_out.num_hits),n_windows);
Total.R_CH=zeros(numel(batch_out.num_hits),n_windows);
Total.R_IN=zeros(numel(batch_out.num_hits),n_windows);

FN_H=fieldnames(Hits);

for s=1:numel(batch_out.num_hits)
    subfieldnames=fieldnames(batch_out.per_tar_pos{s});
    for f=1:numel(subfieldnames)
        FN=subfieldnames{f};
        if isstruct(batch_out.per_tar_pos_base{s}) && isfield(batch_out.per_tar_pos_base{s},FN) % skip not existing run baselines
        batch_out.per_tar_pos{s}.(FN)(:,1)=batch_out.per_tar_pos_base{s}.(FN)(:,2);
        end
    end
    
    Hits.L_CH(s,:)=batch_out.num_hits{s}(1,:);
    Hits.R_CH(s,:)=batch_out.num_hits{s}(2,:);
    Hits.L_CH(s,1)=batch_out.num_hits_B_run{s}(1,2);
    Hits.R_CH(s,1)=batch_out.num_hits_B_run{s}(2,2);
    for p=1:numel(batch_out.target_positions)
        if real(batch_out.target_positions(p))<0
            %            Hits.L_CH(s,:)=Hits.L_CH(s,:)+batch_out.per_tar_pos{s}.hits_CH(p,:);
            Hits.L_IN(s,:)=Hits.L_IN(s,:)+batch_out.per_tar_pos{s}.hits_IN(p,:);
            Total.L_CH(s,:)=Total.L_CH(s,:)+batch_out.per_tar_pos{s}.total_CH(p,:);
            Total.L_IN(s,:)=Total.L_IN(s,:)+batch_out.per_tar_pos{s}.total_IN(p,:);
        elseif real(batch_out.target_positions(p))>0
            %            Hits.R_CH(s,:)=Hits.R_CH(s,:)+batch_out.per_tar_pos{s}.hits_CH(p,:);
            Hits.R_IN(s,:)=Hits.R_IN(s,:)+batch_out.per_tar_pos{s}.hits_IN(p,:);
            Total.R_CH(s,:)=Total.R_CH(s,:)+batch_out.per_tar_pos{s}.total_CH(p,:);
            Total.R_IN(s,:)=Total.R_IN(s,:)+batch_out.per_tar_pos{s}.total_IN(p,:);
        end
    end
    
% % % % %     For choices across all trials
%     Hits.L_CH(s,:)=batch_out.num_hits{s}(1,:)+batch_out.num_hits{s}(2,:);
%     Hits.R_CH(s,:)=batch_out.num_hits{s}(1,:)+batch_out.num_hits{s}(2,:);
%     Hits.L_CH(s,1)=batch_out.num_hits_B_run{s}(1,2)+batch_out.num_hits_B_run{s}(2,2);
%     Hits.R_CH(s,1)=batch_out.num_hits_B_run{s}(1,2)+batch_out.num_hits_B_run{s}(2,2);
%     for p=1:numel(batch_out.target_positions)
%         if real(batch_out.target_positions(p))<0
%             %            Hits.L_CH(s,:)=Hits.L_CH(s,:)+batch_out.per_tar_pos{s}.hits_CH(p,:);
%             Hits.L_IN(s,:)=Hits.L_IN(s,:)+batch_out.per_tar_pos{s}.hits_IN(p,:);
%             Total.L_IN(s,:)=Total.L_IN(s,:)+batch_out.per_tar_pos{s}.total_IN(p,:);
%         elseif real(batch_out.target_positions(p))>0
%             %            Hits.R_CH(s,:)=Hits.R_CH(s,:)+batch_out.per_tar_pos{s}.hits_CH(p,:);
%             Hits.R_IN(s,:)=Hits.R_IN(s,:)+batch_out.per_tar_pos{s}.hits_IN(p,:);
%             Total.R_IN(s,:)=Total.R_IN(s,:)+batch_out.per_tar_pos{s}.total_IN(p,:);
%         end
%             Total.R_CH(s,:)=Total.R_CH(s,:)+batch_out.per_tar_pos{s}.total_CH(p,:);
%             Total.L_CH(s,:)=Total.L_CH(s,:)+batch_out.per_tar_pos{s}.total_CH(p,:);
%     end
    
    
    % Indexes
    %
    Side_chosen_binary{2}=[ones(Hits.L_CH(s,2),1); zeros(Hits.R_CH(s,2),1)];
    Side_chosen_early=[ones(batch_out.early_late_index{s}.n_early_L_CH,1); zeros(batch_out.early_late_index{s}.n_early_R_CH,1)];
    Status=[ones(numel(Side_chosen_binary{2}),1); zeros(numel(Side_chosen_early),1)];
    Binary_data=[Side_chosen_binary{2}; Side_chosen_early];
    p_fisher_bias_early(s) = fexact(Binary_data,Status);
    if batch_out.early_late_index{s}.n_early_L_CH ==0 && batch_out.early_late_index{s}.n_early_R_CH==0
        p_fisher_bias_early(s) = NaN;
    end
    
    Side_chosen_late=[ones(batch_out.early_late_index{s}.n_late_L_CH,1); zeros(batch_out.early_late_index{s}.n_late_R_CH,1)];
    Status=[ones(numel(Side_chosen_binary{2}),1); zeros(numel(Side_chosen_late),1)];
    Binary_data=[Side_chosen_binary{2}; Side_chosen_late];
    p_fisher_bias_late(s) = fexact(Binary_data,Status);
    if batch_out.early_late_index{s}.n_late_L_CH ==0 && batch_out.early_late_index{s}.n_late_R_CH==0
        p_fisher_bias_late(s) = NaN;
    end
    
    
    Side_chosen_binary{2}=[ones(Hits.L_CH(s,2),1); zeros(Hits.R_CH(s,2),1)];
    for w=1:n_windows
        % BIAS
        Side_chosen_binary{w}=[ones(Hits.L_CH(s,w),1); zeros(Hits.R_CH(s,w),1)];
        Status=[ones(numel(Side_chosen_binary{2}),1); zeros(numel(Side_chosen_binary{w}),1)];
        Binary_data=[Side_chosen_binary{2}; Side_chosen_binary{w}];
        p_fisher_bias(s,w) = fexact(Binary_data,Status);
        if Hits.R_CH(s,w)==0 && Hits.L_CH(s,w)==0
            p_fisher_bias(s,w) = NaN;
        end
        %Bonferoni
        if w>2
            p_fisher_bias(s,w) = p_fisher_bias(s,w).*(n_windows-2);
        end
        
        
        % Hitrates
        for c=1:numel(FN_H)
            Hitrate.(FN_H{c})(s,:)=Hits.(FN_H{c})(s,:)./Total.(FN_H{c})(s,:);
            Hitrate_plus_minus.(FN_H{c})(s,w)=Hitrate.(FN_H{c})(s,w)-Hitrate.(FN_H{c})(s,2);
            Correct_binary{w}=[ones(Hits.(FN_H{c})(s,w),1); zeros(Total.(FN_H{c})(s,w)-Hits.(FN_H{c})(s,w),1)];
            Correct_binary{2}=[ones(Hits.(FN_H{c})(s,2),1); zeros(Total.(FN_H{c})(s,2)-Hits.(FN_H{c})(s,2),1)];
            Status=[ones(numel(Correct_binary{2}),1); zeros(numel(Correct_binary{w}),1)];
            Binary_data=[Correct_binary{2}; Correct_binary{w}];
            p_fisher.(FN_H{c})(s,w) = fexact(Binary_data,Status);
            if Total.(FN_H{c})(s,w)==0
                p_fisher.(FN_H{c})(s,w) = NaN;
            end
            %Bonferoni
            if w>2
                p_fisher.(FN_H{c})(s,w) = p_fisher.(FN_H{c})(s,w).*(n_windows-2);
            end
        end
        
    end
    %Bonferoni
    %p_fisher_bias(s,3:end) = p_fisher_bias(s,3:end).*sum(~isnan(p_fisher_bias(s,3:end)));
    
end

p_L_CH=Hits.L_CH./(Hits.L_CH+Hits.R_CH);
bias_plus_minus=p_L_CH-repmat(p_L_CH(:,2),1,size(p_L_CH,2));

batch_out.n_sig_bias_early.significant_per_session=p_fisher_bias_early<0.05;
batch_out.n_sig_bias_late.significant_per_session=p_fisher_bias_late<0.05;

batch_out.n_sig_bias.significant_per_session=p_fisher_bias<0.05;
batch_out.n_sig_bias.plus.all=sum(bias_plus_minus>0,1);
batch_out.n_sig_bias.minus.all=sum(bias_plus_minus<0,1);
batch_out.n_sig_bias.plus_significant.all=sum(bias_plus_minus>0  & p_fisher_bias<0.05,1);
batch_out.n_sig_bias.minus_significant.all=sum(bias_plus_minus<0 & p_fisher_bias<0.05,1);

for c=1:numel(FN_H)
    batch_out.n_sig_hits.plus.(FN_H{c})=sum(Hitrate_plus_minus.(FN_H{c})>0,1);
    batch_out.n_sig_hits.minus.(FN_H{c})=sum(Hitrate_plus_minus.(FN_H{c})<0,1);
    batch_out.n_sig_hits.plus_significant.(FN_H{c})=sum(Hitrate_plus_minus.(FN_H{c})>0  & p_fisher.(FN_H{c})<0.05,1);
    batch_out.n_sig_hits.minus_significant.(FN_H{c})=sum(Hitrate_plus_minus.(FN_H{c})<0 & p_fisher.(FN_H{c})<0.05,1);
end

for s=1:numel(batch_out.num_hits)
    % plot_matrix={idx_L_CH_S, idx_R_CH_S,  idx_L_IN_S, idx_R_IN_S; idx_L_CH_B,  idx_R_CH_B, idx_L_IN_B, idx_R_IN_B};
    Conditions={'L_CH','R_CH','L_IN','R_IN'};
    for c=1:numel(Conditions)
        
        %% RTs
       
        RTs.(Conditions{c})(s,:)=[batch_out.raw_rt{s}{:,c}]';
        RTs.(Conditions{c})(s,1)=[batch_out.raw_rt_B_run{s}{2,c}]';
        idx_empty= cellfun(@isempty,RTs.(Conditions{c})(s,:),'UniformOutput',false);
        RTs.(Conditions{c})(s,find([idx_empty{:}]))=deal({NaN});
        
%         anovainput=[RTs.(Conditions{c}){s,:}];
        
%         P_anova=anovan(anovainput,GROUP);
%         if P_anova>0.05
%             disp('anova on RTs failed')
%         end
        
        meanRTs.(Conditions{c})(s,:)=batch_out.mean_rt{s}(:,c)';
        meanRTs.(Conditions{c})(s,1)=batch_out.mean_rt_B_run{s}(2,c)';
        plus_minus_RT.(Conditions{c})(s,:)=meanRTs.(Conditions{c})(s,:)-meanRTs.(Conditions{c})(s,2);
        
        %[~,p_ttest_RT.(Conditions{c})(s,1)]=ttest2(RTs.(Conditions{c}){s,1},RTs.(Conditions{c}){s,2});
        p_ttest_RT.(Conditions{c})(s,1)=1;
        
        for w=1:n_windows
            %[~,p_ttest_RT.(Conditions{c})(s,w)]=ttest2(RTs.(Conditions{c}){s,w},RTs.(Conditions{c}){s,2});
            if strcmp(GLO.stat_to_use,'signed rank')
               X_W=RTs.(Conditions{c}){s,w};
               X_B=RTs.(Conditions{c}){s,2};
               if ~all(isnan(X_W)) &&  ~all(isnan(X_B))
                   aaa=cellfun(@(x) numel(x),RTs.(Conditions{c})(s,2:end));
                   aab=1:numel(aaa);
                   aac=arrayfun(@(x,y)  repmat(y,x,1),aaa,aab,'UniformOutput',false);
                   P1  = kruskalwallis([RTs.(Conditions{c}){s,2:end}]',vertcat(aac{:}),'off');
                   if P1<0.05
                       p_ttest_RT.(Conditions{c})(s,w)=ranksum(X_B(~isnan(X_B)),X_W(~isnan(X_W)));
                   else
                       p_ttest_RT.(Conditions{c})(s,w)=1;
                   end
               else
                   p_ttest_RT.(Conditions{c})(s,w)=1;
               end
            end
            %Bonferoni
            if w>2
                p_ttest_RT.(Conditions{c})(s,w)=p_ttest_RT.(Conditions{c})(s,w)*(n_windows-2);
            end
        end
        batch_out.n_sig_RTs.significant_per_session.(Conditions{c})=p_ttest_RT.(Conditions{c})<0.05;
        batch_out.n_sig_RTs.plus.(Conditions{c})=sum(plus_minus_RT.(Conditions{c})>0,1);
        batch_out.n_sig_RTs.minus.(Conditions{c})=sum(plus_minus_RT.(Conditions{c})<0,1);
        batch_out.n_sig_RTs.plus_significant.(Conditions{c})=sum(plus_minus_RT.(Conditions{c})>0  & p_ttest_RT.(Conditions{c})<0.05,1);
        batch_out.n_sig_RTs.minus_significant.(Conditions{c})=sum(plus_minus_RT.(Conditions{c})<0 & p_ttest_RT.(Conditions{c})<0.05,1);
        
        %% Velocities
        
        Vel.(Conditions{c})(s,:)=[batch_out.velocities{s}.raw{:,c}]';
        Vel.(Conditions{c})(s,1)=[batch_out.velocities{s}.raw{2,c}]';
        idx_empty= cellfun(@isempty,Vel.(Conditions{c})(s,:),'UniformOutput',false);
        Vel.(Conditions{c})(s,find([idx_empty{:}]))=deal({NaN});
        meanvel.(Conditions{c})(s,:)=batch_out.velocities{s}.mean(:,c)';
        meanvel.(Conditions{c})(s,1)=batch_out.velocities{s}.mean(2,c)';
        plus_minus_vel.(Conditions{c})(s,:)=meanvel.(Conditions{c})(s,:)-meanvel.(Conditions{c})(s,2);
        
        [~,p_ttest_vel.(Conditions{c})(s,1)]=ttest2(Vel.(Conditions{c}){s,1},Vel.(Conditions{c}){s,2});
        p_ttest_vel.(Conditions{c})(s,1)=1;
        
        for w=1:n_windows
            [~,p_ttest_vel.(Conditions{c})(s,w)]=ttest2(Vel.(Conditions{c}){s,w},Vel.(Conditions{c}){s,2});
            if strcmp(GLO.stat_to_use,'signed rank')          
               X_W=Vel.(Conditions{c}){s,w};
               X_B=Vel.(Conditions{c}){s,2};
               if ~all(isnan(X_W)) &&  ~all(isnan(X_B))
               p_ttest_vel.(Conditions{c})(s,w)=ranksum(X_B(~isnan(X_B)),X_W(~isnan(X_W)));
               else
                   p_ttest_vel.(Conditions{c})(s,w)=1;
               end
            end
            %Bonferoni
            if w>2
                p_ttest_vel.(Conditions{c})(s,w)=p_ttest_vel.(Conditions{c})(s,w)*(n_windows-2);
            end
        end
        batch_out.n_sig_velocities.significant_per_session.(Conditions{c})=p_ttest_vel.(Conditions{c})<0.05;
        batch_out.n_sig_velocities.plus.(Conditions{c})                 =sum(plus_minus_vel.(Conditions{c})>0,1);
        batch_out.n_sig_velocities.minus.(Conditions{c})                =sum(plus_minus_vel.(Conditions{c})<0,1);
        batch_out.n_sig_velocities.plus_significant.(Conditions{c})     =sum(plus_minus_vel.(Conditions{c})>0  & p_ttest_vel.(Conditions{c})<0.05,1);
        batch_out.n_sig_velocities.minus_significant.(Conditions{c})    =sum(plus_minus_vel.(Conditions{c})<0 & p_ttest_vel.(Conditions{c})<0.05,1);
        
                %% Accuracy
        
        Acc.(Conditions{c})(s,:)=[batch_out.accuracy_rad{s}.raw_eu{:,c}]';
        Acc.(Conditions{c})(s,1)=[batch_out.accuracy_rad{s}.raw_eu{2,c}]';
        idx_empty= cellfun(@isempty,Acc.(Conditions{c})(s,:),'UniformOutput',false);
        Acc.(Conditions{c})(s,find([idx_empty{:}]))=deal({NaN});
        meanacc.(Conditions{c})(s,:)=batch_out.accuracy_rad{s}.mean_eu(:,c)';
        meanacc.(Conditions{c})(s,1)=batch_out.accuracy_rad{s}.mean_eu(2,c)';
        plus_minus_acc.(Conditions{c})(s,:)=meanacc.(Conditions{c})(s,:)-meanacc.(Conditions{c})(s,2);
        
        [~,p_ttest_acc.(Conditions{c})(s,1)]=ttest2(Acc.(Conditions{c}){s,1},Acc.(Conditions{c}){s,2});
        p_ttest_acc.(Conditions{c})(s,1)=1;
        
        for w=1:n_windows
            [~,p_ttest_acc.(Conditions{c})(s,w)]=ttest2(Acc.(Conditions{c}){s,w},Acc.(Conditions{c}){s,2});
            if strcmp(GLO.stat_to_use,'signed rank')
                X_W=Acc.(Conditions{c}){s,w};
                X_B=Acc.(Conditions{c}){s,2};
                if ~all(isnan(X_W)) &&  ~all(isnan(X_B))
                    p_ttest_acc.(Conditions{c})(s,w)=ranksum(X_B(~isnan(X_B)),X_W(~isnan(X_W)));
                else
                    p_ttest_vel.(Conditions{c})(s,w)=1;
                end
            end
            %Bonferoni
            if w>2
                p_ttest_acc.(Conditions{c})(s,w)=p_ttest_acc.(Conditions{c})(s,w)*(n_windows-2);
            end
        end
        batch_out.n_sig_accuracies.significant_per_session.(Conditions{c})=p_ttest_acc.(Conditions{c})<0.05;
        batch_out.n_sig_accuracies.plus.(Conditions{c})                 =sum(plus_minus_acc.(Conditions{c})<0,1);
        batch_out.n_sig_accuracies.minus.(Conditions{c})                =sum(plus_minus_acc.(Conditions{c})>0,1);
        batch_out.n_sig_accuracies.plus_significant.(Conditions{c})     =sum(plus_minus_acc.(Conditions{c})<0  & p_ttest_acc.(Conditions{c})<0.05,1);
        batch_out.n_sig_accuracies.minus_significant.(Conditions{c})    =sum(plus_minus_acc.(Conditions{c})>0  & p_ttest_acc.(Conditions{c})<0.05,1);
        
    end
    
    % differences in conditions
    con_dif_labels                                           = {'R_IN_L_IN','R_CH_L_CH','R_CH_R_IN','L_CH_L_IN'};
    Conditions                                               = {'L_CH','R_CH','L_IN','R_IN'};
    con_dif_idx                                              = [4,3;2,1;2,4;1,3];
    
    for d=1:numel(con_dif_labels)
        plus_minus_RT.(con_dif_labels{d})(s,:)=meanRTs.(Conditions{con_dif_idx(d,1)})(s,:)-meanRTs.(Conditions{con_dif_idx(d,2)})(s,:);
        for w=1:n_windows
            [~,p_ttest_RT.(con_dif_labels{d})(s,w)]=ttest2(RTs.(Conditions{con_dif_idx(d,1)}){s,w},RTs.(Conditions{con_dif_idx(d,2)}){s,w});
            if strcmp(GLO.stat_to_use,'signed rank')
                X_W=RTs.(Conditions{con_dif_idx(d,1)}){s,w};
                X_B=RTs.(Conditions{con_dif_idx(d,2)}){s,w};
                if ~all(isnan(X_W)) &&  ~all(isnan(X_B))
                    p_ttest_RT.(con_dif_labels{d})(s,w)=ranksum(X_B(~isnan(X_B)),X_W(~isnan(X_W)));
                else
                    p_ttest_RT.(con_dif_labels{d})(s,w)=1;
                end
            end
            %Bonferoni
            if w>2
                p_ttest_RT.(con_dif_labels{d})(s,w)=p_ttest_RT.(con_dif_labels{d})(s,w)*(n_windows-2);
            end
        end
        batch_out.n_sig_RTs.plus.(con_dif_labels{d})                    =sum(plus_minus_RT.(con_dif_labels{d})>0,1);
        batch_out.n_sig_RTs.minus.(con_dif_labels{d})                   =sum(plus_minus_RT.(con_dif_labels{d})<0,1);
        batch_out.n_sig_RTs.plus_significant.(con_dif_labels{d})        =sum(plus_minus_RT.(con_dif_labels{d})>0  & p_ttest_RT.(con_dif_labels{d})<0.05,1);
        batch_out.n_sig_RTs.minus_significant.(con_dif_labels{d})       =sum(plus_minus_RT.(con_dif_labels{d})<0 & p_ttest_RT.(con_dif_labels{d})<0.05,1);
        
    end
end
% L_IN=cellfun(@(y) cellfun(@(x) numel(cell2mat(x)),y(2:end,3)), batch_out.raw_rt,'UniformOutput',false)
% L_IN=[L_IN{:}];
% min(min([L_IN{:}]))
% R_IN=cellfun(@(y) cellfun(@(x) numel(cell2mat(x)),y(2:end,4)), batch_out.raw_rt,'UniformOutput',false)
% R_IN=[R_IN{:}];
% min(min(R_IN))



%XX=cell2mat(batch_out.per_tar_pos)
%XX=horzcat(XX.total_IN)';
%XX(1:9:end,:)=[];


% batch_out.per_tar_pos{s}.total_IN(p,:)

a=1;

