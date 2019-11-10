function [m,lower_bound,upper_bound,R2] = FindmSlopeArea(slope_array,area_array)

%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 07/11/2019 %%%%%%%%%%%%%%%%
% Function to calculate the concavity index based on the log-log slope-area 
% relations of fluvial pixels. The concavity index is referred to as m 
% because the linear inversion scheme assumes that n=1.
% Input parameters:
% slope array - vector of size n of slope [L/L] along fluvial pixels.
% area_array - vector of size n of upstream drainage area [L^2] 
%              along fluvial pixels.
% Output:
% m - best-fit, based on linear regression, of the concavity index
% lower_bound and upper_bound - are the 95% confidence interval on m
% R2 - is a regression coefficient of determination.


pos_slope_index = find(slope_array>0); % account only for positive slopes.
plot(log(area_array(pos_slope_index)),log(slope_array(pos_slope_index)),'xr')
[fitobject,gof] = fit(log(area_array(pos_slope_index)),log(slope_array(pos_slope_index)),'poly1');
hold on;
plot(fitobject,'k')

%output data
m = - fitobject.p1;
conf = confint(fitobject,0.95);
upper_bound = -conf(1,1);
lower_bound = -conf(2,1);
R2 = gof.rsquare;

xlabel('ln(Area)')
ylabel('ln(Slope)')
text(min(log(area_array(pos_slope_index)))+1,...
    min(log(slope_array(pos_slope_index)))+1,...
    strcat('m=', num2str(m), ', R^2=', num2str(R2)));
