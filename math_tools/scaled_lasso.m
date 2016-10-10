function [scaled_lasso_sol ,  sigma_hat] = scaled_lasso(y,Phix,PhixSy,PhiSPhiInv,...
                                                      sign_coeff,sigma_old,N,lambda0,...
                                                      tol,maxiter)

niter = 0 ;
%lambda = sigma_old * lambda0 ;
sigma_hat = 0.1 ;

while (abs(sigma_hat(end) - sigma_old) > tol) && (niter<maxiter) 
    niter = niter + 1 ;
    sigma_old = sigma_hat(end) ;
    lambda = sigma_hat(end) * lambda0 ;
    %scaled_lasso_sol = lasso_fourier_v2(PhixSy,PhiSPhiInv,sign_coeff,lambda*N) ;
    scaled_lasso_sol = lasso_fourier_v1(y,Phix,sign_coeff,lambda*N) ;
    sigma_hat = norm(y - Phix * scaled_lasso_sol)/sqrt(2*N) ; % should divide by sqrt(2)     
end
end




%% WORKS - it was only for debugging
% function [scaled_lasso_sol ,  sigma_hat] = scaled_lasso(y,Phix,PhixSy,PhiSPhiInv,...
%                                                       sign_coeff,sigma_old,N,lambda0,...
%                                                       tol,maxiter)
% 
% niter = 0 ;
% %lambda = sigma_old * lambda0 ;
% sigma_hat = [0.1] ;
% 
% while (abs(sigma_hat(end) - sigma_old) > tol) && (niter<maxiter) 
%     niter = niter + 1 ;
%     sigma_old = sigma_hat(end) ;
%     lambda = sigma_hat(end) * lambda0 ;
%     %scaled_lasso_sol = lasso_fourier_v2(PhixSy,PhiSPhiInv,sign_coeff,lambda*N) ;
%     scaled_lasso_sol = lasso_fourier_v1(y,Phix,sign_coeff,lambda*N) ;
%     sigma_hat = [sigma_hat ,  norm(y - Phix * scaled_lasso_sol)/sqrt(N)] ;
%     
%     
% end
% end
