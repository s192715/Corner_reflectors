function [n res] = FitTest(X_data,Y_data,n_max,Thres);
%function that takes in X_data, Y_data, n_max (max polynomial degree to
%try), and Thres and threshold precision.
%Rreturns n, the degree of polynomial that best fits the data, and res
%the norm of the residuals.

n_try = 1:n_max; %range of polynomial degrees that will be tried to fit

Norm_Of_Residual_FitTest = zeros(1,n_max);
for i = 1:n_max
    [pfit S] = polyfit(X_data,Y_data,n_try(i)); %fitting data
    Norm_Of_Residual_FitTest(i) = S.normr; %quality parameter for goodness of fit
    if S.normr <= Thres
        break
    end
end
Norm_Of_Residual_FitTest = Norm_Of_Residual_FitTest(Norm_Of_Residual_FitTest ~= 0);

n = find(Norm_Of_Residual_FitTest == min(Norm_Of_Residual_FitTest));
res = Norm_Of_Residual_FitTest;


