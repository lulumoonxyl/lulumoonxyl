function t_p_mean=plot_relative_p(t_p_mean,t_p_sem,t_abs_p,x,ii,num)
%num i just for numbering the fig.#
a=(x(2)-x(1))/2;
b=x(1)+a;
c=x(end)-a;
x1=b:a*2:c;

figure(num)
hold on ;
subplot(2,3,ii)
patch([x1(1:3)';flipud(x1(1:3)')],[t_p_mean(1:3)-t_p_sem(1:3);flipud(t_p_mean(1:3)+t_p_sem(1:3))],[1 0.41 0.16],'EdgeColor',[1 1 1])
patch([x1(3:7)';flipud(x1(3:7)')],[t_p_mean(3:7)-t_p_sem(3:7);flipud(t_p_mean(3:7)+t_p_sem(3:7))],[0.6 0.7 0.8],'EdgeColor',[1 1 1])
patch([x1(7:12)';flipud(x1(7:12)')],[t_p_mean(7:12)-t_p_sem(7:12);flipud(t_p_mean(7:12)+t_p_sem(7:12))],[0.8 0.7 0.8],'EdgeColor',[1 1 1])
patch([x1(12:16)';flipud(x1(12:16)')],[t_p_mean(12:16)-t_p_sem(12:16);flipud(t_p_mean(12:16)+t_p_sem(12:16))],[0.6 1 0.8],'EdgeColor',[1 1 1])
patch([x1(16:18)';flipud(x1(16:18)')],[t_p_mean(16:18)-t_p_sem(16:18);flipud(t_p_mean(16:18)+t_p_sem(16:18))],[1 0.41 0.16],'EdgeColor',[1 1 1])
hold on;
plot(x1',t_p_mean,'r-','LineWidth',0.5);
xlabel('Heading direction(degree)')
ylabel('Relative Probability')

figure(1+num)
hold on ;
subplot(2,3,ii)
plot(x1',t_abs_p,'r-','LineWidth',2);
xlabel('Heading direction(degree)')
ylabel('Relative Probability')
hold off;
end