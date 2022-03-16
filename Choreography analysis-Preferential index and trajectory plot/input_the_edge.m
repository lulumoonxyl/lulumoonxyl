function edge=input_the_edge(range)
while (true)
    prompt="Indicate the range you want to eliminate in the form [xmin,xmax;ymin,ymax] (if no input is received, defualt will be [5,222;25,250]): ";
    edge=input(prompt);
    if isempty(edge)
        edge=[5,222;25,250];
        break

    elseif isnumeric(edge)
        [row,col]=size(edge);
        if row==2&&col==2&&edge(1,1)>=range(1,1)&&edge(1,2)<=range(1,2)&&edge(2,1)>=range(2,1)&&edge(2,2)<=range(2,2)
            break
        else
            disp("Please provide the correct input");
        end

    else
        disp("Please provide the correct input");
    end
end
end