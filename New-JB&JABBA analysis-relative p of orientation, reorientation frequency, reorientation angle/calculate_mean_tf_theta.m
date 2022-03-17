function data=calculate_mean_tf_theta(result,data,ii)
%this function will calculate the mean value of tf and theta
fn=fieldnames(result);
mean_fn=append('mean_sem_',fn(2:end));
for i=2:length(fn)
    for j=1:length(result.x)-1
    sem_value=std(result.(fn{i}){j,1})/sqrt(length(result.(fn{i}){j,1}));
    value=mean(result.(fn{i}){j,1});
    
    data.(mean_fn{i-1})(length(result.x)-j,2*ii)=sem_value;
    data.(mean_fn{i-1})(length(result.x)-j,2*ii-1)=value;
    end 
end

end