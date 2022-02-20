function i=fix_figure_title_etc(i,j)
%i is the figure number, j is how many subplot in the figure 
figure(i)
for l=1:j
subplot(2,3,l) %this number is based on the figure 
hold on
%xlim([30 190]);
ylim([35 60]);
axis square
%ylabel('Large Turning Frequency (s^-^1)');
%axis square
end 
% 
figure(i)
subplot(1,2,1)

hold on;
ylim([-0.06 0.11])
plot([3 5],[-0.032 -0.032 ], '-k', 'LineWidth',2)
plot([3.8 4 4.2], [-0.035 -0.035 -0.035] ,'*k', 'LineWidth',1)
text([2.5], [-0.042] ,'ns', 'FontSize',12)
% figure(i)
% subplot(2,3,1)
% title('10^-^1 GA, eliminate the first 2min');
% subplot(2,3,2)
% title('10^-^2 GA, eliminate the first 2min');
% subplot(2,3,3)
% title('10^-^3 GA, eliminate the first 2min');
% subplot(2,3,4)
% title('10^-^5 EA, eliminate the first 2min');
% subplot(2,3,5)
% title('H2O, eliminate the first 2min');

%change it based on the uname 
%eliminate the first 2min

figure(i)
hold on
subplot(2,3,1)
title('10^-^1 GA');
subplot(2,3,2)
title('10^-^2 GA');
subplot(2,3,3)
title('10^-^3 GA');
subplot(2,3,4)
title('10^-^5 EA');
subplot(2,3,5)
title('H2O');
% hold off

% figure(i)
% subplot(1,1,1)
% hold on;
% title('Total PI for week3, eliminate the first 2min');
% ylim([-0.1 0.1]);
% hold off;
% 
% figure(i)
% hold on;
% subplot(1,2,2);
% title('PI based on time for week3, eliminate the first 2min');
% ylim([-0.1 0.25]);
% hold off;
end 

