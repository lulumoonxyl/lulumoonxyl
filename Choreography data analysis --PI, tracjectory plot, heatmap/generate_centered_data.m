function [dat,dat_grouped,fields,fields_grouped]=generate_centered_data(dat,dat_grouped,fields,fields_grouped)
%move all the data points to the origin (0,0)
x_pos=find(matches(fields,'x'));
y_pos=find(matches(fields,'y'));
a=length(fields);
x_pos_grouped=find(matches(fields_grouped,'x'));
y_pos_grouped=find(matches(fields_grouped,'y'));
l=length(fields_grouped);
for i=1:length(dat_grouped{x_pos_grouped,1})
    x=dat_grouped{x_pos_grouped,1}{i,1}(1,1);
    y=dat_grouped{y_pos_grouped,1}{i,1}(1,1);
    dat_grouped{l+1,1}{i,1}=dat_grouped{x_pos_grouped,1}{i,1}-x;
    dat_grouped{l+2,1}{i,1}=dat_grouped{y_pos_grouped,1}{i,1}-y;
end
fields_grouped{l+1}='xcentered';
fields_grouped{l+2}='ycentered';
fields{a+1}='xcentered';
fields{a+2}='ycentered';

dat{a+1,1}=vertcat(dat_grouped{l+1,1}{:});
dat{a+2,1}=vertcat(dat_grouped{l+2,1}{:});

end