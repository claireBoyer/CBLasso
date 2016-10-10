function lasso_sol = lasso_fourier_v1(y,Phix,sign_coeff,lambda)
% with computation of (Phi^H Phi)^(-1)
lasso_sol = (Phix\y - lambda*pinv(Phix'*Phix)*sign_coeff );
end