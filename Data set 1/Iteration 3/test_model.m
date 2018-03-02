function [error, P_on_target, P_binding_rate] = test_model(host_site, sg_RNA, binding_rate,...
                            params, Number_of_Nucleotides, ...
                            off_target_sites, integration_time, integration_accuracy,...
                            Cas_concentration, crRNA_concentration, on_target_copy_number, ...
                            off_target_copy_number, off_target_weight, ODE_parameters, ...
                            off_target_pams, add_plot)
                        
    %%%%%========================================
    %set number of off target sites included
    %%%%%========================================
    number_of_off_target_sites = length(off_target_sites(:,1));
    
    %%%%%========================================
    %set individual model parameters
    %%%%%========================================
    anchor = params(1:2);
    bias = params(3);
    forward = params(4:37);
    backword = params(38:70);
    
        %%%%%========================================
    %set mismatch parameters
    %%%%%========================================
    %[GG, GC, GA, GT, CG, CC, CA, CT, AG, AC, AA, AT, TG, TC, TA, TT ];
    match_mismatch_params = -1*[-1.97, -2.70, -1.66, -2.04,...
                             -1.44, -1.97, -0.78, -1.29,...
                             -1.29, -2.04, -1.04, -1.27,...
                             -0.78, -1.66, -0.12, -1.04]';
    %%%%%========================================
    % set pam sites with corrosponding values
    %%%%%========================================
    % ['ATG', 'GAG', 'GGA', 'GAA', 'GAT'] 0 = G, 1 = C, 2 = A, 3 = T
    pam_key = [2, 3, 0; 0, 2, 0; 0, 0, 2; 0, 2, 2; 0, 2, 3];
    pam_vals = -1*[-6.6, -6.9, -9.4, -6.8, -6.9]/max(-1*[-6.6, -6.9, -9.4, -6.8, -6.9]);
    on_target_pam = 1; %(-2.70-2.70-1.27)
    off_target_pam_values = zeros(length(off_target_sites(:,1)),1);
    for ii = 1:off_target_sites(:,1)
        [~, index] = ismember(off_target_pams(ii,:),pam_key,'rows');
        off_target_pam_values(ii) = pam_vals(index);
    end
    
    %%%%%========================================
    %initialization
    %%%%%========================================
    error = 0;
    batch_size = length(binding_rate);
    P_on_target = zeros(batch_size,1);
    P_binding_rate = zeros(batch_size,1);
    %%%%%========================================
    %calculate probabilities and set error
    %%%%%========================================
    for kk = 1:batch_size
        
        %%%%%========================================
        % Calculate binding probability for current site
        %%%%%========================================
        N_genome_i = [host_site(kk,:); sg_RNA(kk,:)]';
        [Ratio, ~, ~, ~, ~] = model_probabilities_free(N_genome_i, match_mismatch_params,...
                        anchor, Number_of_Nucleotides, ...
                        forward, backword, bias, on_target_pam); % assumes cannonical PAM
        P_on_target(kk) = BDC_Analytic_Soln(Ratio);
        
        %%%%%========================================
        % Calculate binding probability of off target sites binding
        %%%%%========================================
        batch_size_off_target = length(off_target_sites(:,1));
        P_off_target = zeros(batch_size_off_target,1); 
        for ii = 1:batch_size_off_target 
            N_genome_i = [off_target_sites(ii,:); sg_RNA(kk,:)]';
            [Ratio, ~, ~, ~, ~] = model_probabilities_free(N_genome_i, match_mismatch_params,...
                        anchor, Number_of_Nucleotides, ...
                        forward, backword, bias, off_target_pams(ii));
            P_off_target(ii) =  BDC_Analytic_Soln(Ratio);
        end
        
        %%%%%========================================
        % remove off target sites from integration with P = 0
        %%%%%========================================
        remove = zeros(size(P_off_target));
        for ii = 1:length(P_off_target)
            if P_off_target(ii) > 10^(-32)
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
        P_binding_rate(kk) = predicted_binding_rate;
        
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
        
        %%%%%========================================
        %plot values
        %%%%%========================================
        if add_plot
            close all
            hold on
            % plot figures
            figure(1)
            plot(t, 1-x(:,5))
            title('Percent Binding of On-Target Promotor Site')
            xlabel('Time')
            ylabel('Percentage of Bound Sites')
            figure(2)
            plot(t, x(:,1))
            title('dCas9 Expression')
            xlabel('Time')
            ylabel('Number of Active Molecules')
            figure(3)
            plot(t, x(:,2))
            title('crRNA Expression')
            xlabel('Time')
            ylabel('Number of Active Molecules')
            figure(4)
            plot(t, x(:,3))
            title('Intermediate Complex Expression')
            xlabel('Time')
            ylabel('Number of Active Molecules')
            figure(5)
            plot(t, x(:,4))
            title('Isomerized Complex Expression')
            xlabel('Time')
            ylabel('Number of Active Molecules')
            figure(6)
            for ii = 6:length(x(1,:))
                plot(t, 1-x(:,ii))
                hold on
            end   
            title('Percent Binding of Off-Target Sites')
            xlabel('Time')
            ylabel('Number of Active Molecules')
            pause
        end
    end
    
end