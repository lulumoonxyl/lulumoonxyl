function i=fix_figure_title_etc(i,j)
%i is the figure number, j is how many subplot in the figure 
figure(i)
for l=1:j
subplot(2,3,l) %this number is based on the figure 
hold on
 xlim([-170 170]);
% ylim([15 255]);
ylim([0.035 0.075]);
axis square
end 

 figure(i)
subplot(2,3,1)
title('10^-^1 GA,week1');
subplot(2,3,2)
title('10^-^2 GA,week1');
subplot(2,3,3)
title('10^-^3 GA,week1');
subplot(2,3,4)
title('10^-^5 EA,week1');
subplot(2,3,5)
title('H2O,week1');
%change it based on the uname 
end 