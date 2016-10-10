function [roots_detected,roots_all] = detection_roots(fc,method_coeff,tol)

aux =- conv(method_coeff,flipud(conj(method_coeff))) ; % poly coeff of (- |p|^2)
aux(2*fc+1)=1+aux(2*fc+1);
roots_all = roots(flipud(aux)); % roots function convention


% Isolate roots on the unit circle
% tol = 1e-2;
roots_detected = roots_all(abs(1-abs(roots_all)) < tol);

[~,I]=sort(angle(roots_detected));
roots_detected = roots_detected(I); % roots of 1-|p|^2
roots_detected = roots_detected(1:2:end); % roots of 1-|p|

end