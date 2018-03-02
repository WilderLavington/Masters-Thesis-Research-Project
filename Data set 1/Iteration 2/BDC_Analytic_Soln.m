function soln = BDC_Analytic_Soln(p,q)
    a_0 = 0;
    N = length(q);
    if all(p)
        for j = 1:N
           prod = 1;
           for k = 1:j
               prod = prod*(q(k));
           end
           for k = 1:j
               prod = prod/p(k);
           end
           a_0 = a_0 + prod;
        end
        a_0 = a_0 + 1;
    elseif all(q)
        a_0 = 1;
    else
        a_0 = inf;
    end
    
    soln = 1/a_0;
end