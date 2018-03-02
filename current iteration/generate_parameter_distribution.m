
% load primary data 
[binding_data, sgRNA_Tarasava, host_Tarasava] = load_seq_data;
data_set_1 = binding_data(1:30);
data_set_2 = binding_data(31:60);
data_set_3 = binding_data(61:90);
% find MLEs for each data point given assumption of distribution of data
normal_dist = @(x,sigma) normpdf(x,0,sigma);

for datapoint = 1:length(binding_data)
    data = [data_set_1(ii), data_set_2, data_set_3];
    
end
% simulate 100 new data realizations

% train distributions to .1 accuracy or 5*10^4 iterations

% look at how each distance parameter is distributed


