

# DTWGrangerFramework

The framework of VL-Granger causality inference

The main file of this framework is "DTW_granger_cause.m".

Varible-Lag Granger Causality using DTW
Developer: Chainarong Amornbunchornvej, UIC, 2018
Revised code based on [*] by Chandler Lutz
We use the HSIC Independence test in [*2] by Arthur Gretton, 2007.

[*] https://www.mathworks.com/matlabcentral/fileexchange/25467-granger-causality-test

[*2] http://people.kyb.tuebingen.mpg.de/arthur/indep.htm

 Acknowledgements:
   This code cannot be developed without the based-code by Chandler Lutz.
   I would like to thank Chandler Lutz, Arthur Gretton to provide the public code, 
   and all comments in the matlab community.

 [GfStatTest,follSimVal, notIndp,BICres,pval, F,c_v, threshIndp,testStatIndp] = DTW_granger_cause(x,y,alpha,max_lag,sigmaTHS)
 
 INPUT: x,y is a time series of real numbers where x(t) is a real value at time t
 
 INPUT: alpha, is the significance level of both f-test and independence test
 
 INPUT: max_lag is a maximum, time lag and sigmaTHS is a following relation threshold
 
 OUTPUT: GfStatTest = 1 if f-test of VR-Granger reject the H_0 that y does not VR-Granger cuases x with sig-level alpha, otherwise zero.
 
 OUTPUT: follSimVal is between 0.1 to 1 if x follows y,  follSimVal is between -1 to -0.1 if y follows l , and around zero if no following
 relation between x and y.

============
