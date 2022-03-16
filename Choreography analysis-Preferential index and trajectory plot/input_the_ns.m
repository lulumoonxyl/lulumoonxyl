function ns=input_the_ns(range_of_recording)
while (true)
    prompt="Indicate how many seconds want to eliminate (if no input is received, it will be set to 120): ";
    ns=input(prompt);
    if isempty(ns)
        ns=120;
        break
    elseif isnumeric(ns)&& ns<=range_of_recording &&ns>=0
        break
    else
        disp("Please provide the correct input");
    end
end
end