function [original_sites, off_target_sites, scores] = load_off_target_data()
    %%%%%========================================
    % load data
    %%%%%======================================== 
    load handel
    [~, original_sequence] = xlsread('off_target_sites_data_set_1.xlsx','Sheet1','A:A');
    [~, off_target_sequence] = xlsread('off_target_sites_data_set_1.xlsx','Sheet1','B:B');
    scores = xlsread('off_target_sites_data_set_1.xlsx','Sheet1','C:C');
    %%%%%========================================
    % Initialize data matrices
    %%%%%========================================
    original_sites = zeros(length(scores), length(off_target_sequence{2})); 
    off_target_sites = zeros(length(scores), length(original_sequence{2})); 
    %%%%%========================================
    % 0 = G, 1 = C, 2 = A, 3 = T
    %%%%%========================================
    % convert host sequences
    for i = 1:length(off_target_sequence)
        for j = 1:length(off_target_sequence{i})
            if off_target_sequence{i}(j) == 'A'
               off_target_sites(i,j) = 2; 
            elseif off_target_sequence{i}(j) == 'T'
               off_target_sites(i,j) = 3; 
            elseif off_target_sequence{i}(j) == 'G'
               off_target_sites(i,j) = 0; 
            else % 'C'
               off_target_sites(i,j) = 1; 
            end
        end
    end
    % convert sequences complimentary to original_sequence
    % assumes that I loaded data as a complentary site (switches base pairs)
    for i = 1:length(original_sequence)
        for j = 1:length(original_sequence{i})
            if original_sequence{i}(j) == 'A'
               original_sites(i,j) = 3; 
            elseif original_sequence{i}(j) == 'T'
               original_sites(i,j) = 2; 
            elseif original_sequence{i}(j) == 'G'
               original_sites(i,j) = 1; 
            else % 'C'
               original_sites(i,j) = 0; 
            end
        end
    end
end