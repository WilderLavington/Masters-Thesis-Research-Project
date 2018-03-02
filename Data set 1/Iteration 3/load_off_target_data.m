function [original_sites, off_target_sites, scores, copy_number, pam_sites] = load_off_target_data()
    %%%%%========================================
    % load data
    %%%%%======================================== 
    load handel
    [~, off_target_sequence] = xlsread('off_target_wt_PAM.xlsx','Sheet1','A:A');
    [~, pam_sequence] = xlsread('off_target_wt_PAM.xlsx','Sheet1','B:B');
    copy_number = xlsread('off_target_wt_PAM.xlsx','Sheet1','C:C');
    scores = xlsread('off_target_wt_PAM.xlsx','Sheet1','D:D');
    [~, original_sequence] = xlsread('off_target_wt_PAM.xlsx','Sheet1','E:E');
    %%%%%========================================
    % Initialize data matrices
    %%%%%========================================
    original_sites = zeros(length(original_sequence)-1, length(original_sequence{2})); 
    off_target_sites = zeros(length(off_target_sequence)-1, length(off_target_sequence{2})); 
    pam_sites = zeros(length(pam_sequence)-1, length(pam_sequence{2})); 
    %%%%%========================================
    % 0 = G, 1 = C, 2 = A, 3 = T
    %%%%%========================================
    % convert host sequences
    for i = 2:length(off_target_sequence)
        for j = 1:length(off_target_sequence{i})
            if off_target_sequence{i}(j) == 'A'
               off_target_sites(i-1,j) = 2; 
            elseif off_target_sequence{i}(j) == 'T'
               off_target_sites(i-1,j) = 3; 
            elseif off_target_sequence{i}(j) == 'G'
               off_target_sites(i-1,j) = 0; 
            else % 'C'
               off_target_sites(i-1,j) = 1; 
            end
        end
    end
    % convert sequences complimentary to original_sequence
    % assumes that I loaded data as a complentary site (switches base pairs)
    for i = 2:length(off_target_sequence)
        for j = 1:length(original_sequence{2})
            if original_sequence{2}(j) == 'A'
               original_sites(i-1,j) = 3; 
            elseif original_sequence{2}(j) == 'T'
               original_sites(i-1,j) = 2; 
            elseif original_sequence{2}(j) == 'G'
               original_sites(i-1,j) = 1; 
            else % 'C'
               original_sites(i-1,j) = 0; 
            end
        end
    end
    % convert sequences complimentary to original_sequence
    for i = 2:length(pam_sequence)
        for j = 1:length(pam_sequence{i})
            if pam_sequence{i}(j) == 'A'
               pam_sites(i-1,j) = 2; 
            elseif pam_sequence{i}(j) == 'T'
               pam_sites(i-1,j) = 3; 
            elseif pam_sequence{i}(j) == 'G'
               pam_sites(i-1,j) = 0; 
            else % 'C'
               pam_sites(i-1,j) = 1; 
            end
        end
    end
end