function summary_scatterplot_frenzy(varargin)
global analysis_parameters
if nargin==0
    monkey=analysis_parameters.monkey;
    if     analysis_parameters.batch_processing
        type=analysis_parameters.batches.type;
        excentricity=analysis_parameters.batches.excentricity;
    else
        type=NaN;
        excentricity='';
    end
else
    monkey=varargin{1};
    type=varargin{2};
    excentricity='?';
end

%saving the figures and pdf export
current_folder=[pwd filesep 'scatterplots'];
if exist(current_folder) ~= 7
    mkdir(current_folder);
end
pdffilename=[current_folder filesep monkey ' all scatterplots'];
delete([pdffilename '.pdf']);


if type==1
    to_plot=1:11;
elseif type==2
    to_plot=12:37;
else
    to_plot=1:37;
end
for k=to_plot
    
reduce_to_residuals=0;
    switch k
        case 1
            X_par='electrode depth'; % 'current'; % 'probability of evoking saccades'; %
            Y_par='Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; % 'Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; %
            Selection='fix'; % 'freegaze'; % 'go';  % 'shifting'; % 'all'; %
        case 2
            X_par='current'; % 'probability of evoking saccades'; % 'electrode depth'; %
            Y_par='Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; % 'Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; %
            Selection='fix'; % 'freegaze'; % 'go';  % 'shifting'; % 'all'; %
        case 3
            X_par='electrode depth'; % 'current'; % 'probability of evoking saccades'; %
            Y_par='Amplitudes of evoked saccades'; % 'Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; %
            Selection='fix'; % 'freegaze'; % 'go';  % 'shifting'; % 'all'; %
        case 4
            X_par='current'; % 'probability of evoking saccades'; % 'electrode depth'; %
            Y_par='Amplitudes of evoked saccades'; % 'Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; %
            Selection='fix'; % 'freegaze'; % 'go';  % 'shifting'; % 'all'; %
        case 5
            X_par='electrode depth'; % 'current'; % 'probability of evoking saccades'; %
            Y_par='Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='fix'; % 'freegaze'; % 'go';  % 'shifting'; % 'all'; %
        case 6
            X_par='current'; % 'probability of evoking saccades'; % 'electrode depth'; %
            Y_par='Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='fix'; % 'freegaze'; % 'go';  % 'shifting'; % 'all'; %
        case 7
            X_par='probability of evoking saccades'; % 'electrode depth'; % 'current'; %
            Y_par='Amplitudes of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='fix';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %
        case 8
            X_par='Amplitudes of evoked saccades'; % 'electrode depth'; % 'current'; %
            Y_par='Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='fix';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %
        case 9
            X_par='probability of evoking saccades'; % 'electrode depth'; % 'current'; %
            Y_par='Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='fix';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %
            
         case 10
            X_par='current'; % 'probability of evoking saccades'; % 'electrode depth'; %
            Y_par='Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; % 'Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; %
            Selection='fix'; % 'freegaze'; % 'go';  % 'shifting'; % 'all'; %   
            reduce_to_residuals=1;
            
          case 11
            X_par='current'; % 'probability of evoking saccades'; % 'electrode depth'; %
            Y_par='Amplitudes of evoked saccades'; % 'Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; %
            Selection='fix'; % 'freegaze'; % 'go';  % 'shifting'; % 'all'; %   
            reduce_to_residuals=1;
            
        case 12
            X_par='electrode depth'; % 'current'; % 'probability of evoking saccades'; %
            Y_par='Probability of delaying saccades'; % 'Amplitudes of evoked saccades'; % 'Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; %
            Selection='shifting+'; % 'freegaze'; % 'go';  % 'shifting'; % 'all'; %
        case 13
            X_par='current'; % 'probability of evoking saccades'; % 'electrode depth'; %
            Y_par='Probability of delaying saccades'; % 'Amplitudes of evoked saccades'; % 'Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; %
            Selection='shifting+'; % 'freegaze'; % 'go';  % 'shifting'; % 'all'; %
        case 14
            X_par='probability of evoking saccades'; % 'current'; % 'probability of evoking saccades'; %
            Y_par='Probability of delaying saccades'; % 'Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; %
            Selection='shifting+'; % 'freegaze'; % 'go';  % 'shifting'; % 'all'; %        
        case 15
            X_par='probability of delaying saccades'; % 'electrode depth'; % 'current'; %
            Y_par='Contraversive target selection'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting+';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %
        case 16
            X_par='current'; % 'electrode depth'; % 'current'; %
            Y_par='Contraversive target selection'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting+';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %
        case 17
            X_par='electrode depth'; % 'electrode depth'; % 'current'; %
            Y_par='Contraversive target selection'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting+';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %
        case 18
            X_par='electrode depth'; % 'current'; % 'probability of evoking saccades'; %
            Y_par='Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; % 'Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; %
            Selection='shifting-'; % 'freegaze'; % 'go';  % 'shifting'; % 'all'; %
        case 19
            X_par='current'; % 'probability of evoking saccades'; % 'electrode depth'; %
            Y_par='Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; % 'Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; %
            Selection='shifting-'; % 'freegaze'; % 'go';  % 'shifting'; % 'all'; %
        case 20
            X_par='electrode depth'; % 'current'; % 'probability of evoking saccades'; %
            Y_par='Amplitudes of evoked saccades'; % 'Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; %
            Selection='shifting-'; % 'freegaze'; % 'go';  % 'shifting'; % 'all'; %
        case 21
            X_par='current'; % 'probability of evoking saccades'; % 'electrode depth'; %
            Y_par='Amplitudes of evoked saccades'; % 'Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; %
            Selection='shifting-'; % 'freegaze'; % 'go';  % 'shifting'; % 'all'; %
        case 22
            X_par='electrode depth'; % 'current'; % 'probability of evoking saccades'; %
            Y_par='Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting-'; % 'freegaze'; % 'go';  % 'shifting'; % 'all'; %
        case 23
            X_par='current'; % 'probability of evoking saccades'; % 'electrode depth'; %
            Y_par='Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting-'; % 'freegaze'; % 'go';  % 'shifting'; % 'all'; %
        case 24
            X_par='probability of evoking saccades'; % 'electrode depth'; % 'current'; %
            Y_par='Amplitudes of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting-';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %
        case 25
            X_par='Amplitudes of evoked saccades'; % 'electrode depth'; % 'current'; %
            Y_par='Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting-';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %
        case 26
            X_par='probability of evoking saccades'; % 'electrode depth'; % 'current'; %
            Y_par='Latencies of evoked saccades'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting-';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %
        case 27
            X_par='probability of evoking saccades'; % 'electrode depth'; % 'current'; %
            Y_par='Contraversive target selection'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting-';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %
        case 28
            X_par='current'; % 'electrode depth'; % 'current'; %
            Y_par='Contraversive target selection'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting-';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %
        case 29
            X_par='electrode depth'; % 'electrode depth'; % 'current'; %
            Y_par='Contraversive target selection'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting-';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %         
            
        case 30
            X_par='probability of evoking saccades'; % 'electrode depth'; % 'current'; %
            Y_par='Contraversive target selection difference'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting-';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %
        case 31
            X_par='current'; % 'electrode depth'; % 'current'; %
            Y_par='Contraversive target selection difference'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting-';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %
        case 32
            X_par='electrode depth'; % 'electrode depth'; % 'current'; %
            Y_par='Contraversive target selection difference'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting-';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %    
            
        case 33
            X_par='probability of delaying saccades'; % 'electrode depth'; % 'current'; %
            Y_par='Contraversive target selection difference'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting+';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %
        case 34
            X_par='current'; % 'electrode depth'; % 'current'; %
            Y_par='Contraversive target selection difference'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting+';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %
        case 35
            X_par='electrode depth'; % 'electrode depth'; % 'current'; %
            Y_par='Contraversive target selection difference'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting+';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %
        case 36
            X_par='probability of evoking saccades'; % 'electrode depth'; % 'current'; %
            Y_par='Contraversive target selection'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting+';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %
        case 37
            X_par='probability of evoking saccades'; % 'electrode depth'; % 'current'; %
            Y_par='Contraversive target selection difference'; % 'Contraversive target selection'; % 'Electrode depth'; % 'Current'; % 'Probability of evoking saccades'; % 'Amplitudes of evoked saccades'; %
            Selection='shifting+';  % 'shifting'; % 'all'; % 'fix'; % 'freegaze'; %
    end
    
    %scatterplot
    scatterplots_microstim(monkey,X_par,Y_par,Selection,reduce_to_residuals)
    
    %saving the figures and pdf export    
%    filename=[current_folder filesep monkey ' ' Y_par ' vs ' X_par ', stim in ' Selection ', ' excentricity];
%    saveas(gcf,[filename '_scatter_plot']);
    export_fig(pdffilename, '-pdf', '-transparent', '-append');
    close gcf    
    
    errorbarplots_microstim(monkey,X_par,Y_par,Selection,reduce_to_residuals)
    
    %saving the figures and pdf export    
 %   filename=[current_folder filesep monkey ' ' Y_par ' vs ' X_par ', stim in ' Selection ', ' excentricity];
%    saveas(gcf,[filename '_errorbar_plot']);
    export_fig(pdffilename, '-pdf', '-transparent', '-append');
    close gcf
end
end