function [tf_mean,tfl_mean,deg_mean]=cal_reorient_f(turn_x_or_pre_deg, ...
    x_or_heading_direction,et,reorient_deg_abs,bins,deg_large,x_or_theta)

%if we want to get tf based on x, we need to put turn_x and x as the first
%2 inputs
%if based on heading direction, we need to put pre_deg and heading
%direction as the first 2 inputs

%deg_large is the degree that defines large turning, ex.
%deg_large=60-->any turning event with reorient_deg_abs>=60 is considered
%as large turning event;

%x_or_theta tells us whether we would calculate the turning frequency,
%large turning freuquency, reorientation angle as a function of x or theta
%('x' or 't')
if isequal(x_or_theta,'x')

    pos=0:bins:230;
elseif isequal(x_or_theta,'t')
    pos=-180:bins:180;
end
l=length(pos)-1;
tf=cell(l,1);tfl=cell(l,1);deg=cell(l,1);%create empty cell
for j=1:l
    w=1;q=1;
    for i=1:length(x_or_heading_direction)

        if j==l
            idx=find(x_or_heading_direction{i,1}>=pos(j)&x_or_heading_direction{i,1}<=pos(j+1));
            %find the timepoint that the larva in certain xbin or heading in certian directions
        else
            idx=find(x_or_heading_direction{i,1}>=pos(j)&x_or_heading_direction{i,1}<pos(j+1));
        end
        if isempty(idx)
            continue
        end
        if idx(1)==1
            %if the first index is equal to 1, we need to make it [],
            %as we need to calculate the time using time(idx) and
            %time(idx-1)
            idx(1)=[];
        end
        if isempty(idx)
            continue
        end

        %calculate the time that this larva spend in this xbin or heading
        %in certian directions
        t=sum(et{i,1}(idx)-et{i,1}(idx-1));
        if j==l
            idx_turn=find(turn_x_or_pre_deg{i,1}>=pos(j)&turn_x_or_pre_deg{i,1}<=pos(j+1));
            %find whether there is any turning events in this xbins or this
            %heading direction
        else
            idx_turn=find(turn_x_or_pre_deg{i,1}>=pos(j)&turn_x_or_pre_deg{i,1}<pos(j+1));
        end

        %calculate the turning frequency: count of turning event/time
        tf{j,1}(w,1)=length(idx_turn)/t;

        %calculate the large turning frequency:count of large turning
        %event/time
        idx_turn_large=find(reorient_deg_abs{i,1}(idx_turn)>=deg_large);
        tfl{j,1}(w,1)=length(idx_turn_large)/t;
        w=w+1;
        %calculate the mean reorientation angle vs x if there is any
        %turning event in the xbins
        if ~isempty(idx_turn)
            deg{j,1}{q,1}=reorient_deg_abs{i,1}(idx_turn);
            q=q+1;
        end
    end
end

for i=1:length(deg)
    deg{i,1}=vertcat(deg{i,1}{:});
end
%in here, we flip the data upside down as the odor is at the left side
%and after we flip it, the odor is at the right side x=220;we do not need
%to flip for heading direction

% if isequal(x_or_theta,'x')
%     tf=flipud(tf);tfl=flipud(tfl);deg=flipud(deg);
% end
% data=table(tf,tfl,deg);
for k=1:l
    %calculate the mean value and the sem 
    tf_mean(k,1)=mean(tf{k,1}); tf_mean(k,2)=std(tf{k,1})./sqrt(length(tf{k,1}));
    tfl_mean(k,1)=mean(tfl{k,1}); tfl_mean(k,2)=std(tfl{k,1})./sqrt(length(tfl{k,1}));
    deg_mean(k,1)=mean(deg{k,1});deg_mean(k,2)=std(deg{k,1})./sqrt(length(deg{k,1}));
end
clear deg tf tfl

end