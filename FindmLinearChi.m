function [m, chi] = FindmLinearChi(x,y,z,rec_array,area_array)

%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 07/11/2019 %%%%%%%%%%%%%%%%
% Function to calculate the concavity index, m (assuming n=1) that produces 
% the most linear river profiles in the chi-z domain that goes through the 
% origin (0,0), and the assoiated chi values  
% Input parameters:
% x,y,z - vectors of size n of x,y, and z cordinate [L] along fluvial pixels.
% rec_array - vector of size n of the reciever relation. A value j in row i
%             means that the reciever of the pixel in row i 
%            (for which x, y, and drainage area are supplied) is in row j. 
%             The reciever of the outlet is the pixel itself.
% area_array - vector of size n of upstream drainage area [L^2] 
%              along fluvial pixels.
% Output:
% m - area power that produces the most linear rivers in the chi-z domain 
%     based on R^2
% chi - vector of size n with the chi values [L] of the pixels produced
%       using m



m_try = 0.1:0.05:0.95; % vactors of various m values to try. 
% If the produced R2 vec is monotonic, it is a good practice to extend the range. 
R2 = zeros(size(m_try));
ft = fittype('a*x');

  
for i = 1:length(m_try)
    chi = CalculateChi(x,y,rec_array,area_array,m_try(i));
    [~,gof] = fit(chi,z,ft);
    R2(i) = gof.rsquare;
end

%first figure shows R^2 as a function of m. Can be used o define a range of
%acceptable m values. 
figure
subplot(2,1,1)
plot(m_try,R2,'ob')
xlabel('m')
ylabel('R^2 of linear regression through \chi-z domain')

%second figure plots the rivers in the chi-z domaim, where chi is
%calculated with the m that produces the most linear rivers in the R^2
%sense together with the linear fit that goes through the origin. An
%important test is to see if the residual has a structure. 
subplot(2,1,2)
hold on;
[maxR2,maxR2_loc] = max(R2);
chi = CalculateChi(x,y,rec_array,area_array,m_try(maxR2_loc));
for i = 1:length(chi)
    j = rec_array(i);
    if j ~= i
        plot([chi(j) chi(i)],[z(j) z(i)],'b');
    end
end

[fitobject,~] = fit(chi,z,ft);
plot(chi,chi*fitobject.a,'k')
xlabel('\chi [m]')
ylabel('z [m]')
title('\chi-z relation using m that best linearize the profile')
text(min(chi)*2,min(z)*2,...
    strcat('m=',num2str(m_try(maxR2_loc)),', R^2=',num2str(maxR2)));
m = m_try(maxR2_loc);
set(gcf,'renderer','Painters')