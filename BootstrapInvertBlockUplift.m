function [Ustar_mat, tstar_best] = BootstrapInvertBlockUplift(chi,z,Gamma,q,percent_smaple,num_iterations,K)

%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 07/11/2019 %%%%%%%%%%%%%%%%
% function to apply the block uplift linear inversion as part of a 
% bootstrap algorith
% Input parameters: 
% chi - a vector of size n with chi [L] values
% z - a vector of size n with elevation data [L]
% Gamma - the dampening coefficient
% q - number of time intervals in the inversion
% percent_smaple - the ratio of data to sample in each inversion
%                  application 0-1
% num_iterations - number of bootstrap iterations
% K - erodibility coeffcient [L^{1-2m}/T]. Used for plotting purposes. To
%     present non-dimensional inversion use K = 1, to present inversion 
%     with dimensional units use the natural K. The plt assumes that A0=1.
%Output:
% Ustar_mat - Matrix of size num_iterationsXq. Each row is the inversion
%             result of a particular bootstarp realization. If K is not 1
%             the inversion results are 
% tstar_best - vector of size q+1 with the boundaries of the time intervals
%              as derived from the inversion using all data.
%              First element is zero 




close all;

non_zero_elements = find(z~=0);
nz_chi = chi(non_zero_elements);
nz_z = z(non_zero_elements);

Ustar_mat = zeros(num_iterations,q);
tstar_mat = zeros(num_iterations,q+1);

%bootstrap iterations with part of the data
data_length= length(non_zero_elements);
sample_length = floor(data_length*percent_smaple);
for i = 1:num_iterations
    sample = sort(randperm(data_length,sample_length));
    [Ustar ,tstar, ~] = InvertBlockUplift(nz_chi(sample),nz_z(sample),Gamma,q,0);
    Ustar_mat(i,:) = Ustar;
    tstar_mat(i,:) = tstar;
end

%inversion with all data
[Ustar_best ,tstar_best, ~] = InvertBlockUplift(nz_chi,nz_z,Gamma,q,0);


%plotting the results
figure
hold on
for i = 1:num_iterations
    Ustar = Ustar_mat(i,:);
    tstar = tstar_mat(i,:);
    tstar_plot = [];
    Ustar_plot = [];
    for j = 1:q
        tstar_plot = [tstar_plot tstar(j) tstar(j+1)];
        Ustar_plot = [Ustar_plot Ustar(j) Ustar(j)]; 
    end
    plot(tstar_plot/K/1e6,Ustar_plot*K/1e-3,'Color', [0.9 0.9 0.9],'LineWidth',1)
end


%add the best fit solution
tstar_plot = [];
Ustar_plot = [];
for j = 1:q
    tstar_plot = [tstar_plot tstar_best(j) tstar_best(j+1)];
    Ustar_plot = [Ustar_plot Ustar_best(j) Ustar_best(j)];
end
plot(tstar_plot/K/1e6,Ustar_plot*K/1e-3,'Color', [0 0 0],'LineWidth',2);

% add mean solution
lightBlue = [91, 207, 244] / 255;
meanUstar = mean(Ustar_mat);
tstar_plot = [];
Ustar_plot = [];
for j = 1:q
    tstar_plot = [tstar_plot tstar_best(j) tstar_best(j+1)];
    Ustar_plot = [Ustar_plot meanUstar(j) meanUstar(j)];
end
plot(tstar_plot/K/1e6,Ustar_plot*K/1e-3,'Color','m','LineWidth',1);

%add std
stdUstar = std(Ustar_mat);
tstar_plot = [];
Ustar_up_plot = [];
Ustar_down_plot = [];
for j = 1:q
    tstar_plot = [tstar_plot tstar_best(j) tstar_best(j+1)];
    Ustar_up_plot = [Ustar_up_plot meanUstar(j)+stdUstar(j)  meanUstar(j)+stdUstar(j)];
    Ustar_down_plot = [Ustar_down_plot meanUstar(j)-stdUstar(j)  meanUstar(j)-stdUstar(j)];
end
plot(tstar_plot/K/1e6,Ustar_up_plot*K/1e-3,'Color','m','LineStyle',':','LineWidth',1);
plot(tstar_plot/K/1e6,Ustar_down_plot*K/1e-3,'Color','m','LineStyle',':','LineWidth',1);


xlabel('t [Ma]','FontSize',20)
ylabel('U [mm/yr]','FontSize',20)
set(gcf,'renderer','Painters')