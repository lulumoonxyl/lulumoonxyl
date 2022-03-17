function [ns,edge,both_ns,both_edge,plot_input]=get_input_from_user(ii_num)
%% get input from user
%do they want to eliminate both the seconds and edge

str_both="Do you want to eliminate the both the first several seconds and the edge 'y'/'n'?: ";
eliminate_both=eliminate_or_not(str_both);
%% input the ns and edge values
if isequal(eliminate_both,'y')
    both_ns=input_the_ns(900);
    both_edge=input_the_edge([0,230;0,255]);
else
    both_ns=[];both_edge=[];
end

%% check whether users want to delete only teh first several seconds

str_ns="Do you want to eliminate only the first several seconds? 'y'/'n': ";
eliminate_ns=eliminate_or_not(str_ns);

if isequal(eliminate_ns,'y')
    ns=input_the_ns(900);
else
    ns=[];
end
%%
str_edge="Do you want to eliminate only the edge? 'y'/'n': ";
eliminate_edge=eliminate_or_not(str_edge);

if isequal(eliminate_edge,'y')
    edge=input_the_edge([0,230;0,255]);
else
    edge=[];
end
%% Plot or not
    str(1,1)="calculate and plot the relative probability 'y'/'n': ";
    str(2,1)="calculate and plot the reorientation frequency 'y'/'n': ";
    str(3,1)="calculate and plot the acceptance rate 'y'/'n': ";
    
    plot_input='';
    for i=1:length(str)
    plot_input=plot_or_not(str(i),i,plot_input);
    end 
    
end