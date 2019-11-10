function [K, U, t] = CalibrateKTotalUplift(H,t_H,Ustar,tstar,A0,m,to_plot)

%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 07/11/2019 %%%%%%%%%%%%%%%%
% function to find the erodibility coeffcient, K, based on Total Uplift
% data and to convert the non-dimensional uplift rate, U*, and the scaled
% time, t*, to uplif rate, U, as a function of time, t, with natural units.
% Assumes block uplift conditions.
% Input parameters:
% H - Total uplift [L] from the present to time t_H in the past
% t_H - age [T] of the dated uplifted feature
% Ustar - vector of length q of the non-dimensional uplift rate history. 
%         The output of the Block uplift inversion procedure. 
% tstar - vector of length q+1 of scaled time that bounds scaled time 
%         intervals over which Ustar is valid
% A0 - reference drainage area used in the production of chi over which the
%      inversion id based. Use A0 = 1
% m - area power used in the production of chi over which the
%     inversion id based.
% to_plot - use 0 to supress plotting of the dimensional uplift rate
%           history
% Output:
% K - erodibility coefficient constrained based on input parameters.
% U - vector of the length q of the natural units [L/T] uplift rate.
% t - vector of length q+1 of natural units [T] bounds of time intervals.

del_t_star = diff(tstar);
testH_vec = cumsum(Ustar.*del_t_star);
time_ind = find(testH_vec > H,1) - 1;
rem_H = H - testH_vec(time_ind);
del_scaled_time = rem_H/Ustar(time_ind+1);
t_H_star = tstar(time_ind+1) + del_scaled_time;
K = t_H_star/(t_H*A0^m);

U = Ustar*K;
t = tstar/K;


%plot the uplift history
if to_plot
    figure;
    hold on;
    t_plot = [];
    U_plot = [];
    for i = 1:length(U)
        t_plot = [t_plot t(i) t(i+1)];
        U_plot = [U_plot U(i) U(i)];
    end
    plot(t_plot,U_plot,'LineWidth',2)
    xlabel('t [yr]','FontSize',20)
    ylabel('U [m/yr]','FontSize',20)
end





