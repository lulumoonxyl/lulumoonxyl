function []=save_all_figures(outdir)

if ~isfolder(outdir)
    mkdir(outdir)
end

%% save figures
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName   = get(FigHandle, 'Name');
%      savefig(FigHandle, fullfile(outdir, append(FigName, '.fig')));
     saveas(FigHandle, fullfile(outdir, append(FigName, '.svg')));
end

end