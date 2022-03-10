function [dat_JB,fields_JB]=get_xy_center(dat_JB,fields_JB)
%use this function when the get_data_from_JB_file does not get x_center and
%y_center from trx_reduced.mat
xspine=find(matches(fields_JB,'xspine'));
yspine=find(matches(fields_JB,'yspine'));
l=length(dat_JB);
fields_JB{l+1,1}='x';
fields_JB{l+2,1}='y';
for i=1:length(dat_JB{xspine,1})
    dat_JB{l+1,1}{i,1}=dat_JB{xspine,1}{i,1}(:,6);
    dat_JB{l+2,1}{i,1}=dat_JB{yspine,1}{i,1}(:,6);
end 
end 