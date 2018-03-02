function [average_binding_data, average_sgRNA_Tarasava, average_host_Tarasava] = ...
            load_avg_data
     % load primary data 
    [binding_frequency, Data_matrix_SgRNA, Data_matrix_Host,...
                Data_vector_DataSet, Data_vector_CasType] = load_seq_data;

    % split data (Tarasava)
    binding_frequency_Tarasava = binding_frequency(find(Data_vector_DataSet == 1));
    sgRNA_Tarasava = Data_matrix_SgRNA(find(Data_vector_DataSet == 1),:);
    host_Tarasava = Data_matrix_Host(find(Data_vector_DataSet == 1),:);

    % going to use average now 
    average_binding_data = (binding_frequency_Tarasava(1:30)+binding_frequency_Tarasava(31:60)...
                            + binding_frequency_Tarasava(61:90))/3;
    average_sgRNA_Tarasava = sgRNA_Tarasava(1:30,:);
    average_host_Tarasava = host_Tarasava(1:30,:);
end