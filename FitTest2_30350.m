function [n res error] = FitTest2_30350(Y_data,n_max,i_true,Thres); %THIS function was for method 1 new case (updated case)
%it interpolates few state vectors (n+1 vectors with degree n polynomial),
%and leaves out vector 11 (or i_true), the error is the difference in the value of
%vector 11 and the interpolated polynomial evaluated at vector 11, i.e. how
%well the fitted polynomial predicts the correct value at vector 11.

%i_true = 11
%Function returns n, the degree of polynomial that best fits the data, res
%the norm of the residuals, and error which decides which is the best fit.

X_data = 1:length(Y_data);

Norm_Of_Residual_FitTest = zeros(1,n_max);
Error = zeros(1,n_max);

for n= 1:n_max
    if rem(n,2) == 1 %if n is odd
        i = [i_true-round((n+1)/2):1:i_true+round((n+1)/2)];
        i = i(i~=i_true);
    else %if n is even
        i = [i_true-(n/2):1:i_true+((n+2)/2)];
        i = i(i~=i_true);
    end

[pfit S] = polyfit(X_data(i),Y_data(i),n);%fitting polynomial degree n to n+1 vectors
pval = polyval(pfit,i_true); %evaluating polynomial at missing vector 11

Norm_Of_Residual_FitTest(n) = S.normr;
Error(n) = abs(Y_data(i_true) - pval);

if abs(Y_data(i_true) - pval) <= Thres
    break
end

end
Error = Error(Error~=0); %to remove zero values after implementing threshold
Norm_Of_Residual_FitTest = Norm_Of_Residual_FitTest(Norm_Of_Residual_FitTest~=0);

n = find(Error == min(Error));
res = Norm_Of_Residual_FitTest;
error = Error;
end