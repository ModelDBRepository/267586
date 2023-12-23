function y = expfitdual(beta,x)

y = beta(1) + beta(2)*exp(beta(3)*x(:,1)) + ...
    beta(4)*exp(beta(5)*x(:,1));