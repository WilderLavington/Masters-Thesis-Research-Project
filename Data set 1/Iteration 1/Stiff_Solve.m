function x_prime = Stiff_Solve(t,x,ODE_Parameters, N_total_current, Probability_of_binding)
    %%%% coeffifents key 
    %r_crRNA Delta_crRNA r_Cas9 Delta_Cas9 Delta_crRNA:Cas9  Lambda 
    %k_f k_I D V k_d k_c mu 
    %%%% x key 
    % N_Cas9, N_cRNA, N_intermediate, N_Complex, N_target_1, N_target_2, ... 
    binding_sites = length(N_total_current);
    x_prime = zeros(length(N_total_current)+4,1);
    binding_rate = zeros(length(N_total_current),1);
    %%%% equations 1-4 N_crRNA, N_cas9, N_intermediate, N_complex
    x_prime(1) = ODE_Parameters(4)-ODE_Parameters(5).*x(1)-ODE_Parameters(3).*x(2).*(x(1));
    x_prime(2) = ODE_Parameters(1)-ODE_Parameters(2).*x(2)-ODE_Parameters(3).*x(2).*(x(1));
    x_prime(3) = ODE_Parameters(3)*x(2)*x(1)-(ODE_Parameters(6)+ODE_Parameters(7))*(x(3));
    %%%% cleavage rate calc.
    r_RW = 6*ODE_Parameters(8)*ODE_Parameters(9)*x(4)/ODE_Parameters(10);
    for p = 1:binding_sites
        %%%% include x(4+p) to account for slow decrease in prob of bindind
        %%%% as sites become unavailable
        r_binding = r_RW.*x(4+p)*Probability_of_binding(p);
        binding_rate(p) = ODE_Parameters(11)/(ODE_Parameters(11)+ODE_Parameters(12)).*r_binding;
    end
    x_prime(4) = ODE_Parameters(7)*x(3)-ODE_Parameters(6)*x(4)-sum(binding_rate(:));
    %%%% euations 5-n rate concentration of unbound sites of type j over time.
    for p = 5:(binding_sites+4)
         x_prime(p) = ODE_Parameters(13)*(N_total_current(p-4)-x(p))-binding_rate(p-4);
    end
end

