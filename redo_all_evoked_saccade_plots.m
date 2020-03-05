main_dir='Y:\microstim_behavior\Curius_summaries\Evoked_saccade_plots';
temp=dir(main_dir);
Subfolders={temp([temp.isdir]).name};
for n=3:numel(Subfolders)
    folder=[main_dir filesep Subfolders{n}];
    cd(folder)
    evoked_saccade_streamline;
end