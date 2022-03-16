function plot=plot_or_not(str,i,plot)
while (true)
    prompt=str;
    result=input(prompt);
    plot(i)=result;
    if isequal(result,'y') |isequal(result,'n')
        break
    else
        disp("Please provide the correct format of input 'y'/'n'");
    end
end
end 