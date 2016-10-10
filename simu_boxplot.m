%% This script constructs boxplot on sigma_hat 



clear all
close all

addpath(genpath('~/Documents/MATLAB/cvx')) ;
addpath('tools/') ;
addpath('math_tools/') ;


rng(1);

%% Display features
ms = 20;
lw = 1;
P = 2048*8;
u = (0:P-1)'/P;


%% Generate data

fc = 80;
N = 2*fc+1;

Fourier = @(fc,x)exp(-2i*pi*(-fc:fc)'*x(:)');

Phi  = @(fc,x,a)Fourier(fc,x)*a;
PhiS = @(fc,u,p)real( Fourier(fc,u)'* p ); % CHECK real

delta = 3/fc;


n = 3 ; % nb of spikes
sigma = 1  ; % noise level
max_ntest = 100 ;
sigma_CBLasso = zeros(max_ntest,1) ;
sigma_BLasso = zeros(max_ntest,1) ;
q = ones(N,1);
options.h = q;
alpha = 2 ;
for ntest = 1:max_ntest
    
    
    %% Support drawing under (oversafe) separability condition
    x0 = sort([rand rand rand]');
    while (abs(x0(1)-x0(2))<delta) || (abs(x0(3)-x0(2))<delta) || (abs((1+x0(1))-x0(3))<delta)
        x0 = sort([rand rand rand]');
    end
    a0 = [1 1 -1]';
    Phix0 = Phi(fc,x0,a0) ;
    
    %% Noisy data
    w = randn(N,1) + 1i * randn(N,1);
    w = sigma * w;
    y = Phix0 + w;
    
    
    f0 = PhiS(fc,u,Phix0);
    f = PhiS(fc,u,y);
    
    %% CBLasso
    
    lambda_max_CBLasso = norm(PhiS(fc,u,y),inf) /(norm(y)*sqrt(N)) ; 
    
    % CVX solver for dual polynomial - CBlasso
    lambda_CBLasso = lambda_max_CBLasso/alpha ;  
    cvx_precision best
    cvx_solver sdpt3
    cvx_begin sdp quiet
    variable X(N+1,N+1) hermitian;
    variable c(N) complex;
    X >= 0;
    X(N+1,N+1) == 1;
    X(1:N,N+1) == c .* conj(q);
    trace(X) == 2;
    norm(c)<= 1/(sqrt(N)*lambda_CBLasso) ;
    for j = 1:N-1,
        sum(diag(X,j)) == X(N+1-j,N+1);
    end
    maximize(real(c'*y))
    cvx_end
    
    CBLasso_coeff = X(1:N,N+1);
    
    
    % Roots of dual polynomial
    [roots_broot_detected , roots_all_broot] = detection_roots(fc,CBLasso_coeff,1e-2) ;
    x_CBLasso = angle(roots_broot_detected)/(2*pi);
    x_CBLasso = sort(mod(x_CBLasso,1));
    Phix_CBLasso = Fourier(fc,x_CBLasso);
    sign_coeff_CBLasso = sign(real(Phix_CBLasso'*CBLasso_coeff)); % CHECK real
    
    

    % CBlasso
    PhixSy_CBLasso = Phix_CBLasso'*y ;
    PhiSPhiInv_CBLasso = pinv(Phix_CBLasso'*Phix_CBLasso) ;
    
    sigma_old = norm(y)/sqrt(N) ;
    lambda0 =  lambda_CBLasso ;
    tol = 1e-4 ;
    maxiter = 1000 ;
    [scaled_lasso_sol ,  sigma_CBLasso(ntest)] = scaled_lasso(y,Phix_CBLasso,PhixSy_CBLasso,PhiSPhiInv_CBLasso,...
                                                            sign_coeff_CBLasso,sigma_old,N,...
                                                            lambda0,tol,maxiter) ;
                                                        
    %% Blasso 
    
    
    lambda_max_blasso = norm(PhiS(fc,u,y),inf) / N ; 
    lambda_blasso = N * lambda_max_blasso/alpha ; 
    
    % computing the dual polynomial with Blasso
    cvx_solver sdpt3 % SeDuMi %
        cvx_begin sdp quiet
        cvx_precision high;
        variable X(N+1,N+1) hermitian;
        variable c(N) complex;
        X >= 0;
        X(N+1,N+1) == 1;
        X(1:N,N+1) == c .* conj(q);
        trace(X) == 2;
        for j = 1:N-1,
            sum(diag(X,j)) == X(N+1-j,N+1);
        end
        if lambda_blasso==0
            maximize(real( c' * y ))
        else
            minimize( norm(y/lambda_blasso-c, 'fro') )
        end
        cvx_end
    blasso_coeff = X(1:N,N+1);
    
    % roots-finding
    [roots_blasso_detected , roots_all_blasso] = detection_roots(fc,blasso_coeff,1e-2) ;
    
    x_blasso = angle(roots_blasso_detected)/(2*pi);
    x_blasso = sort(mod(x_blasso,1));
    Phix_blasso = Fourier(fc,x_blasso);
    sign_coeff_blasso = sign(real(Phix_blasso'*blasso_coeff));
    
    % BLasso
    PhixSy = Phix_blasso'*y ;
    PhiSPhiInv = pinv(Phix_blasso'*Phix_blasso) ;
    Blasso_sol = lasso_fourier_v1(y,Phix_blasso,sign_coeff_blasso,lambda_blasso) ;
    sigma_BLasso(ntest) = norm(y - Phix_blasso * Blasso_sol)/sqrt(2*(N-length(Blasso_sol))) ; 
    
    if mod(ntest,10)==0
       disp(ntest)
       % save sigma_CBLasso sigma_CBLasso 
       save sigma_BLasso sigma_BLasso 
    end
    
end
save sigma_CBLasso sigma_CBLasso
save sigma_BLasso sigma_BLasso


%% Display results


fig3 = figure() ;
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
boxplot([sigma_BLasso sigma_CBLasso],'Labels',{'BLasso','CBLasso'})
set(findobj(gca,'Type','text'),'FontSize',16)

% Crop pdf
set(fig3,'Units','Inches');
pos = get(fig3,'Position');
set(fig3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(strcat('boxplot_sigma',num2str(sigma)),'-dpdf','-r0')

