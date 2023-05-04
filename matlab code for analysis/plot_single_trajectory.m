function []=plot_single_trajectory(dat,list,t,ttl_t,tsum,series)
%the output will have three rows and five columns, so it can plot the
%example trajectories for five animals, the "list" input should include
%five colums and n row, series is the x axis of the third row of subplot--hading
%direction bin, w is the figure number

%1) single trajectories, 'Start' is the startpoint of the object, turning
%is labelled with ark dots, red dots are where turning starts, blue dots
%are where turning stop 
%2)The heading direction detected at each timepoints, red line are where
%turing occur 
%3)The histogram plot of time that each object spends in different HD bin, eliminate
%the time when bject is turning

%sum of time is the sum of the time histogram ("tsum" input), total t eliminate turn is
%using the dat.et of each object substrct the turning time ("ttl_t" input) -->they should be the same value 

%"t" input is the time spend in each HD bin for each larvae
for w=1:height(list)
figure(w)
set(gcf,'Name',append('single animal trajectory and hist of time vs HD, ',num2str(w)));
hold on
for i=1:width(list)
idx0=dat.t0_idx{list(w,i),1};idx1=dat.t1_idx{list(w,i),1};
idx3=[];
for j=1:length(idx1)
    idx3=vertcat(idx3,(idx0(j):idx1(j))');
end
subplot(3,5,i)
hold on
plot(dat.x{list(w,i),1},dat.y{list(w,i),1},'g.');
plot(dat.x{list(w,i),1}(idx3),dat.y{list(w,i),1}(idx3),'k.')
plot(dat.x{list(w,i),1}(idx0),dat.y{list(w,i),1}(idx0),'r.');
plot(dat.x{list(w,i),1}(idx1),dat.y{list(w,i),1}(idx1),'b.')
text(dat.x{list(w,i),1}(1),dat.y{list(w,i),1}(1),'start');
hold off

subplot(3,5,i+5)
hold on
plot(dat.et{list(w,i),1}(idx3),ones(length(idx3),1).*(188),'r.','MarkerSize',5);
plot(dat.et{list(w,i),1},dat.orientation{list(w,i),1},'g-');
xlabel('Time(s)');
ylabel('Heading Direction (deg)');
axis tight;
ylim([-180 190]);
hold off

subplot(3,5,i+10)
bar(series,t(list(w,i),:));
title(append('total t(exclude turn)=',num2str(ttl_t(list(w,i))),'s,','tsum=',num2str(tsum(list(w,i))),'s'));
xlabel('Heading direction (deg)');
ylabel('Time (s)')


end 
end 
end