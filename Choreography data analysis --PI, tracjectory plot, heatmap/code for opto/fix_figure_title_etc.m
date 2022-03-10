function i=fix_figure_title_etc(i,j)
%i is the figure number, j is how many subplot in the figure 
figure(i)
for l=1:j
subplot(2,3,l) %this number is based on the figure 
hold on
xlim([0 230]);
ylim([20 260]);
%ylabel('Large turning frequency (s^-^1), theta [-180 -135]&[135 180]');
%axis square
end 

figure(i)
hold on
subplot(2,3,1)
title('10^-^1 GA')
subplot(2,3,2)
title('10^-^2 GA')
subplot(2,3,3)
title('10^-^3 GA')
subplot(2,3,4)
title('10^-^5 EA')
subplot(2,3,5)
title('H2O')
%change it based on the uname 
%eliminate the first 2min
% figure(i)
% hold on
% subplot(3,3,1)
% title('49a@chrimson,150_3 75_2 37_2,eliminate the first 2min');
% subplot(3,3,2)
% title('49a@chrimson,255_3 127_2 63_2,eliminate the first 2min');
% subplot(3,3,3)
% title('49a@chrimson,50_3 25_2 12_2,eliminate the first 2min');
% subplot(3,3,4)
% title('42a@chrimson,150_3 75_2 37_2,eliminate the first 2min');
% subplot(3,3,5)
% title('42a@chrimson,255_3 127_2 63_2,eliminate the first 2min');
% subplot(3,3,6)
% title('42a@chrimson,50_3 25_2 12_2,eliminate the first 2min');
% subplot(3,3,7)
% title('attp2@chrimson,150_3 75_2 37_2,eliminate the first 2min');
% subplot(3,3,8)
% title('attp2@chrimson,255_3 127_2 63_2,eliminate the first 2min');
% subplot(3,3,9)
% title('attp2@chrimson, 50_3 25_2 12_2,eliminate the first 2min');

% hold off
% figure(i)
% subplot(1,2,1)
% hold on;
%  %title('Total PI, eliminate the edge and the first 2min');
% ylim([-0.06 0.03]);
% hold off;
% % 
% figure(i)
% hold on;
% subplot(1,2,2);
% %title('PI based on time, eliminate the edge and the first 2min');
% ylim([-0.1 0.1]);
% hold off;
end 

