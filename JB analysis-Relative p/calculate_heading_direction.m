function [dat,fields]=calculate_heading_direction(dat,fields,vector)
%this function just calculates the heading direction at each timepoint for
%each larva based on JB data
%vector=[1 0] or [-1 0] for the different position; [-1 0] for odor in the
%left side 
xspine_pos=find(matches(fields,'xspine'));
yspine_pos=find(matches(fields,'yspine'));

l=length(dat);
fields{l+1}='heading direction';
for i=1:length(dat{xspine_pos,1})
    x1=dat{xspine_pos,1}{i,1}(:,1);
    x2=dat{xspine_pos,1}{i,1}(:,7);
    x3=x1-x2;
    
    y1=dat{yspine_pos,1}{i,1}(:,1);
    y2=dat{yspine_pos,1}{i,1}(:,7);
    y3=y1-y2;
    % change (-1,0) based on the position of odor
    degree=atan2d(y3.*(vector(1))-vector(2).*x3,x3.*(vector(1))+y3.*vector(2));

    dat{l+1,1}{i,1}=0-degree;
    clear degree x1 x2 x3 y1 y2 y3
end
end