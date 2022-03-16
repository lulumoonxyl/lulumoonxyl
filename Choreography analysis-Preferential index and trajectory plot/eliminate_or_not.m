function eliminate=eliminate_or_not(str)
while (true)
    prompt=str;
    eliminate=input(prompt);
    if isequal(eliminate,'y') |isequal(eliminate,'n')
        break
    else
        disp("Please provide the correct format of input 'y'/'n'");
    end
end
end 