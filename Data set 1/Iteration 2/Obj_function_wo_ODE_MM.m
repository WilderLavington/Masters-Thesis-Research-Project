function [error] = Obj_function_wo_ODE_MM(binding_rate, sg_RNA, host_site, ...
                            params, Number_of_Nucleotides, ...
                            off_target_sites, integration_time, integration_accuracy,...
                            Cas_concentration, crRNA_concentration, on_target_copy_number, ...
                            off_target_copy_number, off_target_weight, ODE_parameters, match_mismatch_params)
    %%%%%========================================
    %set number of off target sites included
    %%%%%========================================
    number_of_off_target_sites = length(off_target_sites(:,1));

    %%%%%========================================
    %set individual model parameters
    %%%%%========================================
    anchor_params = params(1:4);
    bias = params(5:6);
    forward = params(7:39);
    backword = params(40:72);
    
    %%%%%========================================
    %initialization
    %%%%%========================================
    error = 0;
    batch_size = length(binding_rate);
    P_on_target = zeros(batch_size,1);

    %%%%%========================================
    %calculate probabilities and set error
    %%%%%========================================
    for kk = 1:batch_size
        
        %%%%%========================================
        % Calculate binding probability for current site
        %%%%%========================================
        N_genome_i = [host_site(kk,:); sg_RNA(kk,:)]';
        [P_lambda, P_mu] = model_probabilities_free(N_genome_i, match_mismatch_params,...
                        anchor_params, Number_of_Nucleotides, ...
                        forward, backword, bias(1), bias(2));
        P_on_target(kk) = BDC_Analytic_Soln(P_lambda,P_mu);
        
        %%%%%========================================
        % Calculate binding probability of off target sites binding
        %%%%%========================================
        batch_size_off_target = length(off_target_sites(:,1));
        P_off_target = zeros(batch_size_off_target,1); 
        for ii = 1:batch_size_off_target 
            N_genome_i = [off_target_sites(ii,:); sg_RNA(kk,:)]';
            [P_lambda, P_mu] = model_probabilities_free(N_genome_i, match_mismatch_params,...
                        anchor_params, Number_of_Nucleotides, ...
                        forward, backword, bias(1), bias(2));
            P_off_target(ii) = BDC_Analytic_Soln(P_lambda,P_mu);
        end
        
        %%%%%========================================
        % remove off target sites from integration with P = 0
        %%%%%========================================
        remove = zeros(size(P_off_target));
        for ii = 1:length(P_off_target)
            if P_off_target(ii) > 10^(-34)
                remove(ii) = 1;
            end
        end
        P_off_target = P_off_target(find(remove==1));
        off_target_copy_number_c = off_target_copy_number(find(remove==1));
        
        %%%%%========================================
        % Integrate system of equations to get reletive binding rates
        %%%%%========================================
        Initial_conditions = [Cas_concentration,crRNA_concentration, 0, 0, on_target_copy_number', off_target_copy_number_c'];
        tspan = linspace(0,integration_time,integration_time*integration_accuracy);
        N_total = [on_target_copy_number', off_target_copy_number_c'];
        P = [P_on_target(kk), P_off_target'];
        f = @(t,x) Stiff_Solve(t, x, ODE_parameters, N_total, P);
        [t,x] = ode23s(f,tspan,Initial_conditions');
        predicted_binding_rate = on_target_copy_number-x(end,5);
        
        %%%%%========================================
        %Add in scaled off target info for perfectly complimentary site
        %%%%%========================================
        if number_of_off_target_sites ~= 0 
            total_off_target_binding = (sum(off_target_copy_number_c'-x(end,6:end)));
            rescaled = off_target_weight*(total_off_target_binding/number_of_off_target_sites)^2;
            error = error + rescaled;
        end
        %%%%%========================================
        %calculate MSE 
        %%%%%========================================
        error = error + ((binding_rate(kk)-predicted_binding_rate)./((mean(binding_rate))))^2;
       
    end
    dig = 10^10;
    msg1 = sprintf('%.10f', round(error*dig)/dig);
    msg2 = sprintf('%.2f', length(backword));
    msgs1 = ': ';
    msgs2 = '\n';
    fprintf(msg2);
    fprintf(msgs1);
    fprintf(msg1);
    fprintf(msgs2);
    disp('');
end