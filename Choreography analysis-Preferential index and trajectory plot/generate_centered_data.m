function [dat,dat_grouped]=generate_centered_data(dat,dat_grouped)
%move all the data points to the origin (0,0)

for i=1:length(dat_grouped.x)
    x=dat_grouped.x{i,1}(1,1);
    y=dat_grouped.y{i,1}(1,1);
    dat_grouped.xcentered{i,1}=dat_grouped.x{i,1}-x;
    dat_grouped.ycentered{i,1}=dat_grouped.y{i,1}-y;
end

dat.xcentered=vertcat(dat_grouped.xcentered{:});
dat.ycentered=vertcat(dat_grouped.ycentered{:});

end