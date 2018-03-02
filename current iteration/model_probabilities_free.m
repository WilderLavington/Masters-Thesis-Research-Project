function [Ratio, P_lambda, P_mu, lambda, mu] = model_probabilities_free(N_genome_i, match_mismatch_params,...
                        anchor, Number_of_Nucleotides, forward, backword, bias, pam)
    
    %%%========================================================================
    % find what type of binding event each match/mismatch is in sequence then
    % pull the correct coefficients from the match_mismatch_params 
    %%%========================================================================
    
    %%%%% generate scaling function from beta distribution
    forward_bp_cf = forward;
    reverse_bp_cf = [forward(1), backword];
    pam = pam*anchor(1); % needed scaling for this to work
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
    lambda = zeros(Number_of_Nucleotides+2,1); % forward transition
    mu = zeros(Number_of_Nucleotides+2,1); % reverse transition 
    states = 0;
    
    %%%%% calculate rate parameters for chain ( treating every bp as state )
    for current_bp = 0:1:Number_of_Nucleotides+1
        
        states = states + 1;
        % if we are at the pam state
        if current_bp == 0
            
            %%%% forward transition probability
            % initialize forward transition as scaled pam value
            lamda_current = forward_bp_cf(1)*pam;
            
            % add remaining forward values 
            start = 1;
            for bp = start:Number_of_Nucleotides
                c_disp = abs(start - bp)+2;
                if mismatches(bp) == 1
                    lamda_current = lamda_current +...
                        forward_bp_cf(c_disp)*match_mismatch_params(bp_key(bp));
                else
                    lamda_current = lamda_current -...
                        forward_bp_cf(c_disp)*match_mismatch_params((bp_key(bp)));
                end
            end
            % add the remaining anchor value
            lamda_current = lamda_current + forward_bp_cf(c_disp+1)*anchor(2);
            
            %%%% reverse transition probability
            % only consider current state
            mu_current = -forward_bp_cf(1)*pam;
            
        % if we are at the anchor state
        elseif current_bp == Number_of_Nucleotides+1
            
            %%%% forward transition probability
            % only forward transition comes from current state
            lamda_current = forward_bp_cf(1)*anchor(2);
            
            %%%% reverse transition probability
            start = Number_of_Nucleotides;
            mu_current = -reverse_bp_cf(1)*anchor(2);
            
            % add remaining values
            for bp = start:-1:1
                c_disp = abs(start - bp)+2;
                if mismatches(bp) == 1
                    mu_current = mu_current -...
                        reverse_bp_cf(c_disp)*match_mismatch_params(bp_key(bp));
                else
                    mu_current = mu_current +...
                        reverse_bp_cf(c_disp)*match_mismatch_params(bp_key(bp));
                end
            end
            mu_current = mu_current - reverse_bp_cf(c_disp+1)*pam;

        % if we are on the interior
        else
            start = current_bp;
            %%%% forward transition probability
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
            % add the remaining anchor value
            lamda_current = lamda_current + forward_bp_cf(c_disp+1)*anchor(2);
            
            %%%% reverse transition probability
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
            % add the remaining pam value
            mu_current = mu_current - reverse_bp_cf(c_disp+1)*pam;
            
        end
          
        % exponentiate for stability
        lambda(states) = exp(lamda_current);
        mu(states) = exp(mu_current);
        
    end
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
    Ratio = (P_mu./P_lambda).*exp(bias);
end