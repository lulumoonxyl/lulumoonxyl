for i=1:10:100
    figure(1)
    hold on 
    plot(dat_grouped.x{i,1},dat_grouped.y{i,1},'.');
    text(dat_grouped.x{i,1}(1,1),dat_grouped.y{i,1}(1,1),append(num2str(i),'-start'));
    for j=1:10:length(PI_tp{i,1})
        text(dat_grouped.x{i,1}(j,1),dat_grouped.y{i,1}(j,1),num2str(round(PI_tp{i,1}(j,1),2)));
    end 
    hold off
end 
%name=["H2O";"10^-^3 GA";"10^-^2 GA";"10^-^1 GA";"10^-^5 EA"];
%name1=["H2O";"GA";"EA"];
name1=["42a3@C13_r_e_t_i_n_o_l_(_+_)";"42a3@C13_r_e_t_i_n_o_l_(_-_)";"49a@C13_r_e_t_i_n_o_l_(_+_)";"49a@C13_r_e_t_i_n_o_l_(_-_)"];
name=["42a3@C13_retinol(+)";"42a3@C13_retinol(-)";"49a@C13_retinol(+)";"49a@C13_retinol(-)"];

figure(j)
for i =1:3
    subplot(1,3,i)
    hold on 
    axis square
%       axis normal
%        xlim([0 100]);
        ylim([0 0.1]);
%         rlim([0 0.08]);
%      ylabel('Frequency')
%       xlabel("Running duration(s)")
%        xlabel("Reorientation angle(deg)");
    %title(name1(i),'FontSize',20);
    hold off
       ax=gca;
    ax.FontSize=18;
end 
     
     axis square
    
end
%        ylabel('Frequency','FontSize',24);
%      xlabel('Tracked duration (s)','FontSize',24);
%ylabel('PI:v_x/s','FontSize',24);
%      ylabel('Y(mm)','FontSize',24);
%      xlabel('X(mm)','FontSize',24);

     %      ylabel('Number of objects tracked','FontSize',24);
     %     ylabel('Distance to odor(mm)','FontSize',24);
     %xlabel('Tracking Time (s)','FontSize',24);
     %    xlabel('Group','FontSize',24);
%            xlabel('Time (s)','FontSize',24);
     title(name1(i),'FontSize',2);
     %     legend(name1,'FontSize',20);
     for i=7:9
         figure(i)
         for j=1:2
             subplot(1,2,j)
             hold on 
             ax=gca;
             ax.FontSize=20;
             ax.FontWeight='Bold';
%              ylim([-1.5 1.5]);
%                rlim([0 8]);
%                yline(0,'r--','LineWidth',1)
               hold off
         end

      end

%% add sig stars for acceptance rate
% *:<0.05; **:<0.01; ***:<0.001
for i=1:3
figure(i)
hold on 
plot([0.85 2.15],[0.55 0.55],'k-','LineWidth',1.5);
text(1.4,0.57,'ns','FontSize',20,'FontWeight','bold');

plot([1 5],[0.03 0.03],'k-','LineWidth',1.5);
text(2.9,0.032,'ns','FontSize',20,'FontWeight','bold');
plot([2.85 4.15],[0.55 0.55],'k-','LineWidth',1.5);
text(3.4,0.57,'ns','FontSize',20,'FontWeight','bold');


plot([0.85 1.85],[0.55 0.55],'k-','LineWidth',1.5);
text(1.3,0.57,'ns','FontSize',20,'FontWeight','bold');

plot([2.85 3.85],[0.55 0.55],'k-','LineWidth',1.5);
text(3.3,0.56,'**','FontSize',40,'FontWeight','bold');

plot([1.15 2.15],[0.58 0.58],'k-','LineWidth',1.5);
text(1.5,0.6,'ns','FontSize',20,'FontWeight','bold');

plot([3.15 4.15],[0.58 0.58],'k-','LineWidth',1.5);
text(3.45,0.59,'**','FontSize',40,'FontWeight','bold');
end 








plot([0.85 1.85],[0.55 0.55],'k-','LineWidth',1.5);
text(1.4,0.57,'**','FontSize',40);

plot([1 2],[0.6 0.6],'k-','LineWidth',1.5);
text(1.4,0.63,'***','FontSize',24);

plot([3 4],[0.1 0.1],'k-','LineWidth',2);
text(3.4,0.106,'***','FontSize',30);

plot([1.25 2.25],[0.138 0.138],'k-','LineWidth',2);
text(1.8,0.1345,'*','FontSize',30);

plot([.75 1.75],[0.1 0.1],'k-','LineWidth',2);
text(1.2,0.106,'***','FontSize',30);

plot([2.75 3.75],[-0.165 -0.165],'k-','LineWidth',2);
text(3.25,-0.17,'**','FontSize',30);

plot([2.75 3.75],[-0.12 -0.12],'k-','LineWidth',2);
text(3.25,-0.126,'***','FontSize',30);
plot([3 4],[0.057 0.057],'k-','LineWidth',2);
text(2.4,0.062,'**','FontSize',30);

plot([1 2],[-0.055 -0.055],'k-','LineWidth',2);
text(1.25,-0.062,'ns','FontSize',20);

plot([3 4],[-0.055 -0.055],'k-','LineWidth',2);
text(3.25,-0.058,'***','FontSize',30);

plot([3 4],[0.1 0.1],'k-','LineWidth',2);
text(3.5,0.11,'***','FontSize',30);

plot([3 4],[0.05 0.05],'k-','LineWidth',2);
text(3.35,0.056,'***','FontSize',30);

plot([3 4],[0.1 0.1],'k-','LineWidth',2);
text(3.35,0.11,'***','FontSize',30);
ylim([-0.1 0.2])