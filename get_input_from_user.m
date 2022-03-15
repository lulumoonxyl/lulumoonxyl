function [ns,edge,both_ns,both_edge,plot_input]=get_input_from_user(ii_num)
%%
prompt="Do you want to eliminate the both the first several seconds and the edge 'y'/'n'?: ";
eliminate_both=input(prompt);
if isequal(eliminate_both,'y') 
        prompt="Indicate how many seconds want to eliminate (default is 120s): ";
        both_ns=input(prompt);
        if isempty(both_ns)
            both_ns=120;
        end
        prompt="Indicate the range you want to eliminate in the form [xmin,xmax;ymin,ymax] (default is [5,222;25,250]): ";
        both_edge=input(prompt);
        if isempty(both_edge)
            both_edge=[5,222;25,250];
        end
   
else
    both_ns=[];both_edge=[];
end

%%
prompt="Do you want to eliminate only the first several seconds? 'y'/'n': ";
eliminate_ns=input(prompt);
if isequal(eliminate_ns,'y')
    prompt="Indicate how many seconds want to eliminate (default is 120s): ";
    ns=input(prompt);
    if isempty(ns)
        ns=120;
    end
else
    ns=[];
end
%%
prompt="Do you want to eliminate only the edge? 'y'/'n': ";
eliminate_edge=input(prompt);

if isequal(eliminate_edge,'y')
    prompt="Indicate the range you want to eliminate? [xmin,xmax;ymin,ymax] (default is [5,222;25,250]): ";
    edge=input(prompt);
    if isempty(edge)
        edge=[5,222;25,250];
    end
else
    edge=[];
end
%%
prompt="Do you plot the trajectories 'y'/'n': ";
plot_input(1)=input(prompt);

prompt="Do you plot the centered trajectories 'y'/'n': ";
plot_input(2)=input(prompt);

prompt="Do you plot the heatmap of the trajectories 'y'/'n': ";
plot_input(3)=input(prompt);

prompt="Do you plot the tracking time vs x position 'y'/'n': ";
plot_input(4)=input(prompt);

prompt="Do you plot preferential index 'y'/'n': ";
plot_input(5)=input(prompt);
end