function yhat = anterograde_funct(beta, funct_parameters, time_exp)
Sp = beta(1);
Dp = beta(2);
delta_t = funct_parameters(1,1);
mRNA = funct_parameters(2:end,1);
yhat = zeros(length(time_exp),1);
for kk = 2:length(time_exp)
    yhat(kk,1) = (mRNA(kk,1)*delta_t*Sp + yhat(kk-1,1))/(1+delta_t*Dp);
end
