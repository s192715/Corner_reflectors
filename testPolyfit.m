function [residual, degree] = testPolyfit(n_max,i_true,res_max,DopplerEqVec)
residual=1e10;
for n= 1:n_max

i = [i_true-round((n+1)/2):1:i_true+round((n+1)/2)]; 
i = i(i~=i_true);%testing points

X = DopplerEqVec(i);%corresponding dopplereq datapoints
fit  = polyfit(i, X,n);%the fitting
test=polyval(fit,i_true);%test variable
res=DopplerEqVec(1,i_true)-test;
    if abs(residual-res)>res_max
        if abs(residual)>abs(res)&& residual~=res
            residual=res;
            degree=n;
        end
    end
end
