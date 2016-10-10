function lasso_sol = lasso_fourier_v2(PhixSy,PhiSPhiInv,sign_coeff,lambda)
% without computation of (Phi^H Phi)^(-1)
lasso_sol =  PhiSPhiInv*(PhixSy -lambda*sign_coeff) ; 
end