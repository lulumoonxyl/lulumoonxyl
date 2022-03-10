function n=get_sample_size(fields_grouped,dat_grouped,ii,n)
PI_mean=find(matches(fields_grouped,'PI mean'))
n(ii)=length(dat_grouped{PI_mean,1});
end 