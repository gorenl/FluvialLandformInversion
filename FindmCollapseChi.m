function [m,chi] = FindmCollapseChi(x,y,z,rec_array,area_array)

%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 07/11/2019 %%%%%%%%%%%%%%%%
% Function to calculate the concavity index, m (assuming n=1) that collapse  
% triutaries along a single curve in the chi-z domain, 
% and the assoiated chi values.
% For each try of m value the function calculates the associated chi
% vector. Then, it splits the data into 25 bins with equal data in each
% bin. It calcultes the std on z in each chi bin and finds the average of
% this std over the bins. The m values that produces the minimal average is
% considered as a natural concavity index. The function plots the average
% of the std for different m values to help define an accptable range for
% which the mean(std) is not very different from the minimal value. 
% Input parameters:
% x,y,z - vectors of size n of x,y, and z cordinate [L] along fluvial pixels.
% rec_array - vector of size n of the reciever relation. A value j in row i
%             means that the reciever of the pixel in row i 
%            (for which x, y, and drainage area are supplied) is in row j. 
%             The reciever of the outlet is the pixel itself.
% area_array - vector of size n of upstream drainage area [L^2] 
%              along fluvial pixels.
% Output:
% m - area power that maximizes collapse of tributaries in the chi-z domain 
%     based min(mean(std(z in bin)))
% chi - vector of size n with the chi values [L] of the pixels produced
%       using m




m_try = 0.1:0.05:0.95;% vactors of various m values to try. 
% If the trend of scat_mertic is monotonic, it is a good practice to extend the range. 
scat_metric = zeros(size(m_try));
n_bins = 25; % can be modified

%figure   
for i = 1:length(m_try)
    chi = CalculateChi(x,y,rec_array,area_array,m_try(i));
    %calculate mean of scatter
    val_in_bins = floor(length(chi)/n_bins);
    chi_scat = zeros(1,n_bins);
    sorted_chi_z = sortrows([chi,z],1);
    for k = 1:n_bins-1
        chi_scat(k) = std(sorted_chi_z(val_in_bins*(k-1)+1:val_in_bins*k,2));
    end
    chi_scat(end) = std(sorted_chi_z(val_in_bins*(n_bins-1)+1:end,2));
    scat_metric(i) = mean(chi_scat);
   
end
%first figure to plot the measn(std(over bins)) for each value of m 
figure
subplot(2,1,1)
plot(m_try,scat_metric,'ob')
xlabel('m')
ylabel('mean of z scatter in bins')

%second figure to plot the chi-z relation that minimizes the scatter
subplot(2,1,2)
hold on;
[min_scat_metric,min_scat_metric_loc] = min(scat_metric);
m = m_try(min_scat_metric_loc);
chi = CalculateChi(x,y,rec_array,area_array,m);
for i = 1:length(chi)
    j = rec_array(i);
    if j ~= i
        plot([chi(j) chi(i)],[z(j) z(i)],'b');
    end
end
xlabel('\chi [m]')
ylabel('z [m]')
title('\chi-z relation using m that minimizes scatter')
text(min(chi)*2,max(z)*0.7,...
    strcat('m=',num2str(m),', min of mean of scatter=',...
    num2str(min_scat_metric)));
set(gcf,'renderer','Painters')
