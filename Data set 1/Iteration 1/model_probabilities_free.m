function [P_lambda, P_mu, lambda, mu] = model_probabilities_free(N_genome_i, match_mismatch_params,...
                        anchors, Number_of_Nucleotides, forward, backword, forward_bias, reverse_bias)
    
    %%%========================================================================
    % find what type of binding event each match/mismatch is in sequence then
    % pull the correct coefficients from the match_mismatch_params 
    %%%========================================================================
    
    %%%%% generate scaling function from beta distribution
    x = linspace(0,1,33);
    forward_bp_cf = forward;
    reverse_bp_cf = backword;

    %%%%% create hash between base pairs and parameters
    % 0 = G, 1 = C, 2 = A, 3 = T
    sequence_key = [0 0; 0 1; 0 2; 0 3; 1 0; 1 1; 1 2; 1 3; 2 0; 2 1; 2 2; 2 3; 3 0; 3 1; 3 2; 3 3];
    reduced_key = [1 2 3 4 2 5 6 7 3 6 8 9 4 7 9 10; 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];
    bp_key = zeros(Number_of_Nucleotides,1);
    
    %%%%% assign complementarity values for each of the 10 combinations 
    mismatches = mismatch_locations(N_genome_i);
    for ii = 1:Number_of_Nucleotides
        for jj = 1:16
            if  N_genome_i(ii,1) == sequence_key(jj,1) &&  N_genome_i(ii,2) == sequence_key(jj,2)
                bp_key(ii) = reduced_key(2,jj);
            end
        end
    end

    %%%%% initialize rate parameters for chain
    Number_of_Matches = Number_of_Nucleotides; %nnz(mismatches);
    lambda = zeros(Number_of_Matches,1);
    mu = zeros(Number_of_Matches,1);
    anchor_start = anchors(1:2);
    anchor_end = anchors(3:4); 
    states = 0;
    
    %%%%% calculate rate parameters for chain ( treating every bp as state )
    for current_bp = 1:Number_of_Nucleotides
        %%%%% new state found
        if mismatches(current_bp) == 1 || mismatches(current_bp) == 0
            states = states + 1;
            start = current_bp;
            
            if current_bp == 1
                %%%% forward information
                lamda_current = 0;
                for bp = start:Number_of_Nucleotides
                    c_disp = abs(start - bp)+1;
                    if mismatches(bp) == 1
                        lamda_current = lamda_current +...
                            forward_bp_cf(c_disp)*match_mismatch_params((bp_key(bp)));
                    else
                        lamda_current = lamda_current -...
                            forward_bp_cf(c_disp)*match_mismatch_params((bp_key(bp)));
                    end
                end
                lamda_current = lamda_current + forward_bp_cf(c_disp+1)*anchor_end(1);
                
                %%%% reverse information
                mu_current = 0;
                mu_current = mu_current - reverse_bp_cf(1)*match_mismatch_params(bp_key(bp));
                mu_current = mu_current - reverse_bp_cf(2)*anchor_start(2);
                
                %%%% add bias
                lamda_current = lamda_current + forward_bias;
                mu_current = mu_current + reverse_bias;
                
            elseif current_bp == Number_of_Nucleotides
                
                %%%% forward information
                lamda_current = 0;
                lamda_current = lamda_current + forward_bp_cf(1)*match_mismatch_params(bp_key(current_bp));
                lamda_current = lamda_current + forward_bp_cf(2)*anchor_end(1);
                
                %%%% reverse information
                mu_current = 0;
                for bp = start:-1:1
                    c_disp = abs(start - bp)+1;
                    if mismatches(bp) == 1
                        mu_current = mu_current -...
                            reverse_bp_cf(c_disp)*match_mismatch_params(bp_key(bp));
                    else
                        mu_current = mu_current +...
                            reverse_bp_cf(c_disp)*match_mismatch_params(bp_key(bp));
                    end
                end
                mu_current = mu_current + reverse_bp_cf(c_disp+1)*anchor_start(2);
                
                %%%% add bias
                lamda_current = lamda_current + forward_bias;
                mu_current = mu_current + reverse_bias;
                
            else
                %%%% forward information
                lamda_current = 0;
                for bp = start:Number_of_Nucleotides
                    c_disp = abs(start - bp)+1;
                    if mismatches(bp) == 1
                        lamda_current = lamda_current +...
                            forward_bp_cf(c_disp)*match_mismatch_params(bp_key(bp));
                    else
                        lamda_current = lamda_current -...
                            forward_bp_cf(c_disp)*match_mismatch_params(bp_key(bp));
                    end
                end
                lamda_current = lamda_current + forward_bp_cf(c_disp+1)*anchor_end(1);
                
                %%%% reverse information
                mu_current = 0;
                for bp = start:-1:1
                    c_disp = abs(start - bp)+1;
                    if mismatches(bp) == 1
                        mu_current = mu_current -...
                            reverse_bp_cf(c_disp)*match_mismatch_params(bp_key(bp));
                    else
                        mu_current = mu_current +...
                            reverse_bp_cf(c_disp)*match_mismatch_params(bp_key(bp));
                    end
                end
                mu_current = mu_current + reverse_bp_cf(c_disp+1)*anchor_start(2);
                
                %%%% add bias
                lamda_current = lamda_current + forward_bias;
                mu_current = mu_current + reverse_bias;
            end
            
            % exponentiate for stability
            lambda(states) = exp(lamda_current);
            mu(states) = exp(mu_current);
        end
    end
    %%%%% calculate forward anchor
    forward_anchor = anchor_start(1);
    for bp = 1:Number_of_Nucleotides
        c_disp = abs(1 - bp)+1;
        if mismatches(bp) == 1
            forward_anchor = forward_anchor +...
                forward_bp_cf(c_disp)*match_mismatch_params(bp_key(bp));
        else
            forward_anchor = forward_anchor -...
                forward_bp_cf(c_disp)*match_mismatch_params(bp_key(bp));
        end
    end
    forward_anchor = exp(forward_anchor + forward_bp_cf(c_disp+1)*anchor_end(1) + forward_bias);
    
    %%%%% calculate reverse anchor
    reverse_achor = anchor_end(2);
    for bp = Number_of_Nucleotides:-1:1
        c_disp = abs(start - bp)+1;
        if mismatches(bp) == 1
            reverse_achor = reverse_achor -...
                reverse_bp_cf(c_disp)*match_mismatch_params(bp_key(bp));
        else
            reverse_achor = reverse_achor +...
                reverse_bp_cf(c_disp)*match_mismatch_params(bp_key(bp));
        end
    end
    reverse_achor = exp(reverse_achor + reverse_bp_cf(c_disp+1)*anchor_start(2) + reverse_bias);
 
    %%%%% add two extra states
    lambda = [forward_anchor; lambda; exp(anchor_end(1))];
    mu = [exp(anchor_start(2)); mu; reverse_achor];
    
    %%%%% initialize transition probabilites
    P_lambda = zeros(length(lambda),1);
    P_mu = zeros(length(lambda),1);
    
    %%%%% calculate interior probability parameters for chain
    for current_state = 1:length(lambda)
       P_lambda(current_state) = lambda(current_state)...
           /(lambda(current_state) + mu(current_state));
    end
    for current_state = 1:length(lambda)
        P_mu(current_state) = mu(current_state)...
           /(lambda(current_state) + mu(current_state));
    end
    
end