function [Ustar,tstar,Misfit] = InvertBlockUplift(chi,z,Gamma,q,to_plot)

%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 07/11/2019 %%%%%%%%%%%%%%%%
%function to inverst the fluvial topography using the Block uplift linear
%inversion scheme to produce the history of uplift rate as a function of
%time. When chi is used as in input parameters, the code produces the 
%non-dimensional uplif rate as a function scaled time. If tau is used
%instead of chi the code produces the dimensional uplift rate as a function
%of time. In the latter case the plot axes titles should be modified from
%U* and t* to U and t.
%Input parameters:
%chi - a vector of size n with chi [L] values
%z - a vector of size n with elevation data [L]
%Gamma - the dampening coefficient
%q - number of time intervals in the inversion
%to_plot - if equal to 0 function doesn't plot the inferred uplift rate 
%          history. Otherwise function plots the history
%Output:
%Ustar - vector of size q with the inferred uplift rate history. First
%element is the prerset.
%tstar - vector of size q+1 with the boundaries of the time intervals. 
%        First element is zero 
%misfit - misfit between data and model topography
%Note: in this implementation the length of the time intervals is not
%constant, but instead, each time interval is resolved by the same number
%of data pixels.

close all
non_zero_elements = find(z~=0);
nz_chi = chi(non_zero_elements);
nz_z = z(non_zero_elements);
N = length(nz_chi);
chi_z_mat = [nz_chi,nz_z];
sorted_chi_z_mat = sortrows(chi_z_mat,1);
sorted_chi=sorted_chi_z_mat(:,1);
sorted_z = sorted_chi_z_mat(:,2);

%constract time intervals such that each time interval is resolved by the
%same number of pixels.

val_per_dt = floor(N/q);

%construct a row vector of scaled times, boundaries of time intervals.  
scaled_t_vec = zeros(1,q+1);
scaled_t_vec(1) = 0;
for i = 2:q
    scaled_t_vec(i) =  sorted_chi((i-1)*val_per_dt);
end
scaled_t_vec(q+1) = sorted_chi(end);

%construct a vector of scaled time intervals 
scaled_dt_vec = diff(scaled_t_vec);

%construct the forward model matrix
Astar = zeros(N,q);
for i = 1:N
    filled_full_elements = find (scaled_t_vec >= sorted_chi(i),1)-2;
    Astar(i,1:filled_full_elements) = scaled_dt_vec(1:filled_full_elements);
    Astar(i,filled_full_elements+1) = sorted_chi(i) - ...
        sum(scaled_dt_vec(1:filled_full_elements));
end

%build the prior model
U_pri_star = ones(q,1)*mean(sorted_z./sorted_chi);

%construct the least square estimate
denom = Astar'*Astar+Gamma^2*eye(q);
nom=Astar'*(sorted_z-Astar*U_pri_star);
Ustar = U_pri_star + denom\nom; %size qX1
tstar = scaled_t_vec'; %size q+1X1 

Misfit = 1/((N-q))*sqrt(sum((sorted_z - Astar*Ustar).^2));

if to_plot   
    figure;
    %plot the uplift history as a straircase plot 
    hold on;
    tstar_plot = [];
    Ustar_plot = [];
    for i = 1:q
        tstar_plot = [tstar_plot tstar(i) tstar(i+1)];
        Ustar_plot = [Ustar_plot Ustar(i) Ustar(i)];
    end
    plot(tstar_plot,Ustar_plot,'LineWidth',2)
    xlabel('t^* [m]','FontSize',20)
    ylabel('U^*','FontSize',20) 
    figure
    %plot data z vs. model z. 
    plot(sorted_z,Astar*Ustar,'x')
    xlabel('data z [m]')
    ylabel('modelled z [m]')
end
    