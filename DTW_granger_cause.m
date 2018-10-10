function [GfStatTest,follSimVal, notIndp,BICres,pval, F,c_v, threshIndp,testStatIndp] = DTW_granger_cause(x,y,alpha,max_lag,sigmaTHS)
%============
% Varible-Lag Granger Causality using DTW
% Developer: Chainarong Amornbunchornvej, UIC, 2018
% Revised code based on [*] by Chandler Lutz
% We use the HSIC Independence test in [*2] by Arthur Gretton, 2007.
%[*] https://www.mathworks.com/matlabcentral/fileexchange/25467-granger-causality-test
%[*2] http://people.kyb.tuebingen.mpg.de/arthur/indep.htm
% Acknowledgements:
%   This code cannot be developed without the based-code by Chandler Lutz.
%   I would like to thank Chandler Lutz, Arthur Gretton to provide the public code, 
%   and all comments in the matlab community.
%
% [GfStatTest,follSimVal, notIndp,BICres,pval, F,c_v, threshIndp,testStatIndp] = DTW_granger_cause(x,y,alpha,max_lag,sigmaTHS)
% INPUT: x,y is a time series of real numbers where x(t) is a real value at
% time t
% INPUT: alpha, is the significance level of both f-test and independence
% test
% INPUT: max_lag is a maximum, time lag and sigmaTHS is a following
% relation threshold
% OUTPUT: GfStatTest = 1 if f-test of VR-Granger reject the H_0 that y does
% not VR-Granger cuases x with sig-level alpha, otherwise zero.
% OUTPUT: follSimVal is between 0.1 to 1 if x follows y,  follSimVal is
% between -1 to -0.1 if y follows l , and around zero if no following
% relation between x and y.
% OUTPUT: notIndp =1 if independence test reject H_0 that x is independent
% with y past with sig-level alpha, otherwise zero. 
% OUTPUT: BICres =1 if BIC of  regression of x on x,y past is greater than
% regression of x on x past, otherwise zero.
% OUTPUT: pval is a p-value of f-test
% OUTPUT:  threshIndp is the threshold of independence test w.r.t. alpha, testStatIndp is the value from data
% ;notIndp =1 if testStatIndp>threshIndp
%============

% ======== legacy comments from Chandler Lutz
% [F,c_v] = granger_cause(x,y,alpha,max_lag)
% Granger Causality test
% Does Y Granger Cause X?
%
% User-Specified Inputs:
%   x -- A column vector of data (effect)
%   y -- A column vector of data (cause)
%   alpha -- the significance level specified by the user
%   max_lag -- the maximum number of lags to be considered
% User-requested Output:
%   F -- The value of the F-statistic
%   c_v -- The critical value from the F-distribution
%
% The lag length selection is chosen using the Bayesian information
% Criterion 
% Note that if F > c_v we reject the null hypothesis that y does not
% Granger Cause x

% Chandler Lutz, UCR 2009
% Questions/Comments: chandler.lutz@email.ucr.edu
% $Revision: 1.0.0 $  $Date: 09/30/2009 $
% $Revision: 1.0.1 $  $Date: 10/20/2009 $
% $Revision: 1.0.2 $  $Date: 03/18/2009 $

% References:
% [1] Granger, C.W.J., 1969. "Investigating causal relations by econometric
%     models and cross-spectral methods". Econometrica 37 (3), 424–438.

% Acknowledgements:
%   I would like to thank Mads Dyrholm for his helpful comments and
%   suggestions


%======= DTW comment
% follSimVal is positive if x follows y and negative if y follows x, otherwise zero
 warning('off','all')
 if (~exist('sigmaTHS','var') || isempty(sigmaTHS))
    sigmaTHS=0.1;
 end
 %=== indepedent test
 params.sigx=-1;
 params.sigy=-1;
 params.bootForce=1;
 params.shuff=500;
 params.q=4;
 %===

%Make sure x & y are the same length
if (length(x) ~= length(y))
    error('x and y must be the same length');
end

%Make sure x is a column vector
[a,b] = size(x);
if (b>a)
    %x is a row vector -- fix this
    x = x';
end

%Make sure y is a column vector
[a,b] = size(y);
if (b>a)
    %y is a row vector -- fix this
    y = y';
end



%Make sure max_lag is >= 1
if max_lag < 1
    error('max_lag must be greater than or equal to one');
end

%First find the proper model specification using the Bayesian Information
%Criterion for the number of lags of x

T = length(x);

BIC = zeros(max_lag,1);

%Specify a matrix for the restricted RSS
RSS_R = zeros(max_lag,1);

i = 1;
while i <= max_lag
    ystar = x(i+1:T,:);
    xstar = [ones(T-i,1) zeros(T-i,i)];
    %Populate the xstar matrix with the corresponding vectors of lags
    j = 1;
    while j <= i
        xstar(:,j+1) = x(i+1-j:T-j);
        j = j+1;
    end
    %Apply the regress function. b = betahat, bint corresponds to the 95%
    %confidence intervals for the regression coefficients and r = residuals
    [b,bint,r] = regress(ystar,xstar);
    
    %Find the bayesian information criterion
    BIC(i,:) = T*log(r'*r/T) + (i+1)*log(T);
    
    %Put the restricted residual sum of squares in the RSS_R vector
    RSS_R(i,:) = r'*r;
    
    i = i+1;
    
end

[dummyX,x_lag] = min(BIC);

%First find the proper model specification using the Bayesian Information
%Criterion for the number of lags of y

BIC = zeros(max_lag,1);

%Specify a matrix for the unrestricted RSS
RSS_U = zeros(max_lag,1);

i = 1;
ystarSet={};
xstarSet={};
while i <= max_lag
    
    ystar = x(i+x_lag+1:T,:);
    xstar = [ones(T-(i+x_lag),1) zeros(T-(i+x_lag),x_lag+i)];
    %Populate the xstar matrix with the corresponding vectors of lags of x
    j = 1;
    while j <= x_lag
        xstar(:,j+1) = x(i+x_lag+1-j:T-j,:);
        j = j+1;
    end
    %Populate the xstar matrix with the corresponding vectors of lags of y
    j = 1;
    while j <= i
        xstar(:,x_lag+j+1) = y(i+x_lag+1-j:T-j,:);
        j = j+1;
    end
    %Apply the regress function. b = betahat, bint corresponds to the 95%
    %confidence intervals for the regression coefficients and r = residuals
    [b,bint,r] = regress(ystar,xstar);
    
    %Find the bayesian information criterion
    BIC(i,:) = T*log(r'*r/T) + (i+1)*log(T);
    
    RSS_U(i,:) = r'*r;
    ystarSet{i}=ystar;
    xstarSet{i}=xstar;
    i = i+1;
    
end

[dummyY,y_lag] =min(BIC);

%============
[nL,LcYdelay,follSimVal] = CreateWarpingTSFunc(y,x);
flag=1;
if LcYdelay>0 && follSimVal>=sigmaTHS
    [nrow,lastInx]=size(xstar);
    %inxdiff=length(nL)-nrow+1;
    %nL=nL(inxdiff:end);
    
    j=1;
    T2=length(nL);
    while j<max_lag
        xstar(:,end+1)=nL(T2-nrow+1-j:T2-j);
        j=j+1;
    end
    [b,bint,r] = regress(ystar,xstar);
    BICsp = T*log(r'*r/T) + (LcYdelay+1)*log(T);
    
    RSS_Usp = r'*r;
    if dummyY > BICsp % if DTW warping provine better model
        [threshIndp,testStatIndp] = hsicTestBoot(ystar,xstar(:,end),alpha,params);
        notIndp=testStatIndp>threshIndp;
        if notIndp==1
            flag=0;
            y_lag=LcYdelay;
        end
    end
end

if flag==1
    %The numerator of the F-statistic
    F_num = ((RSS_R(x_lag,:) - RSS_U(y_lag,:))/y_lag);

    %The denominator of the F-statistic
    F_den = RSS_U(y_lag,:)/(T-(x_lag+y_lag+1));
    
    BICres= dummyY < dummyX;
    
    notIndp=false;
    for k=x_lag+2:size(xstarSet{y_lag},2)
        currxCol=xstarSet{y_lag}(:,k);
        [threshIndp,testStatIndp] = hsicTestBoot(ystarSet{y_lag},currxCol,alpha,params);
        if testStatIndp>threshIndp
            notIndp=true;
            break;
        end
    end
else
    %The numerator of the F-statistic
    F_num = ((RSS_R(x_lag,:) - RSS_Usp)/LcYdelay);

    %The denominator of the F-statistic
    F_den = RSS_Usp/(T-(x_lag+LcYdelay+1));
    
    BICres= BICsp < dummyX;
end

%The F-Statistic
F = F_num/F_den;

c_v = finv(1-alpha,y_lag,(T-(x_lag+y_lag+1)));

pval=1-fcdf(F,y_lag,(T-(x_lag+y_lag+1)));

GfStatTest= F>c_v;
 warning('on','all')
end
    
    
    
    


