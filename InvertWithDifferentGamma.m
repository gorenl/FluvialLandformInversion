function InvertWithDifferentGamma(chi,z,q)

%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 07/11/2019 %%%%%%%%%%%%%%%%
%function to plot an L curve of the misfit as a function of the damping
%coeffcient for Block Uplift inversion
%Input parameters: 
%chi - a vector of size n with chi [L] values
%z - a vector of size n with elevation data [L]
%q - number of time intervals in the inversion
%Output:
% L-curve figure.

close all
% Range of Gamma values. Depending on the data this can be modified.
Gamma_vec = logspace(-1,2,100); 
Misfit = zeros(size(Gamma_vec));
for i = 1:length(Gamma_vec)
    [~,~,Misfit(i)] = InvertBlockUplift(chi,z,Gamma_vec(i),q,0);
end
figure
plot(1./Gamma_vec,Misfit,'LineWidth',2)
axis([-1 max(1./Gamma_vec)*1.1 0 max(Misfit)*1.1])
xlabel('1/\Gamma','FontSize',20)
ylabel('Misfit [m]','FontSize',20)