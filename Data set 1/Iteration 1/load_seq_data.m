function [binding_frequency, Data_matrix_SgRNA, Data_matrix_Host,...
            Data_vector_DataSet, Data_vector_CasType] = load_seq_data
    %%%%%========================================
    % load data
    %%%%%======================================== 
    load handel
    [~, host_site] = xlsread('Formatted_data.xlsx','Sheet1','C:C');
    [~, sgRNA] = xlsread('Formatted_data.xlsx','Sheet1','D:D');
    binding_frequency = xlsread('Formatted_data.xlsx','Sheet1','E:E');
    [~, data_set] = xlsread('Formatted_data.xlsx','Sheet1','A:A');
    [~, cas9_type] = xlsread('Formatted_data.xlsx','Sheet1','B:B');
    %%%%%========================================
    % Initialize data matrices
    %%%%%========================================
    Data_matrix_SgRNA = zeros(length(binding_frequency), length(host_site{2})); 
    Data_matrix_Host = zeros(length(binding_frequency), length(sgRNA{2})); 
    Data_vector_DataSet = zeros(length(binding_frequency), 1);
    Data_vector_CasType = zeros(length(binding_frequency), 1);
    %%%%%========================================
    % 0 = G, 1 = C, 2 = A, 3 = T
    %%%%%========================================
    % convert host sequences
    for i = 1:length(host_site)
        for j = 1:length(host_site{i})
            if host_site{i}(j) == 'A'
               Data_matrix_Host(i,j) = 2; 
            elseif host_site{i}(j) == 'T'
               Data_matrix_Host(i,j) = 3; 
            elseif host_site{i}(j) == 'G'
               Data_matrix_Host(i,j) = 0; 
            else % 'C'
               Data_matrix_Host(i,j) = 1; 
            end
        end
    end
    % convert sequences complimentary to sgRNA
    % assumes that I loaded data as a complentary site (switches base pairs)
    for i = 1:length(sgRNA)
        for j = 1:length(sgRNA{i})
            if sgRNA{i}(j) == 'A'
               Data_matrix_SgRNA(i,j) = 3; 
            elseif sgRNA{i}(j) == 'T'
               Data_matrix_SgRNA(i,j) = 2; 
            elseif sgRNA{i}(j) == 'G'
               Data_matrix_SgRNA(i,j) = 1; 
            else % 'C'
               Data_matrix_SgRNA(i,j) = 0; 
            end
        end
    end
    %convert dataset and cas9 type
    for ii = 1:length(binding_frequency)
        %cas9 type
        if cas9_type{ii}(1) == 'd'
            Data_vector_DataSet(ii) = 1; % 1 = dCas9
        else
            Data_vector_DataSet(ii) = 2; % 2 = Cas9
        end
        %dataset
        if data_set{ii}(1) == 'T'
            Data_vector_DataSet(ii) = 1; % 1 = Tarasava
        else
            Data_vector_DataSet(ii) = 2; % 2 = Hsu
        end
    end
end