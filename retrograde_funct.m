function yhat = retrograde_funct(beta, funct_parameters, time_exp)
Sp = beta(1);
Dp = beta(2);

K = funct_parameters(1,1);
r = funct_parameters(1,2);
ti = funct_parameters(1,3);

P_exp = K ./ (1+exp(-r .* (time_exp - ti)));
dP_exp = (r .* K .* exp(-r .* (time_exp-ti))) ./ ...
    (exp(-r .* (time_exp-ti))+1) .^2;
yhat = (dP_exp + Dp .* P_exp) ./ Sp;

end