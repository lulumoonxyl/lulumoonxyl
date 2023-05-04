function []=change_fig_prop(fig_num,varargin)
figure(fig_num)
hold on
for i=1:2:length(varargin)
    if strcmp(varargin{i},"title")
        title(varargin{i+1});
    elseif strcmp(varargin{i},"xtitle")
        xlabel(varargin{i+1});
    elseif strcmp(varargin{i},"ytitle")
        ylabel(varargin{i+1});
    elseif strcmp(varargin{i},"name")
        set(gcf,'Name',(varargin{i+1}));

    elseif strcmp(varargin{i},"ylim")
        ylim(varargin{i+1});
    elseif strcmp(varargin{i},"xlim")
        xlim(varargin{i+1});
    end
end

hold off
end