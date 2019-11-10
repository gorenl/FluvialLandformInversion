function [Up,tstar,Misfit] = InvertParabola(chi,z,x,rec_array,Gamma,q,K,to_plot)
%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 07/11/2019 %%%%%%%%%%%%%%%%
%function to 
%
%
%
%Input parameters:
%chi - a vector of size n with chi [L] values
%z - a vector of size n with elevation data [L]
% x - vector of size n of the x coordinate [L] of each pixel. Here, x
%     coordinates are assumed to be along strike, in the direction where U
%     changes
%rec_array - vector of size n of the reciever relation. A value j in row i
%             means that the reciever of the pixel in row i 
%            (for which x, y, and drainage area are supplied) is in row j. 
%             The reciever of the outlet is the pixel itself.
%Gamma - the dampening coefficient
%q - number of time intervals in the inversion
% K - erodibility coeffcient [L^{1-2m}/T]. Used for plotting purposes. To
%     present non-dimensional inversion use K = 1, to present inversion 
%     with dimensional units use the natural K. The plt assumes that A0=1.
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
nz_x = x(non_zero_elements);
%nz_rec = rec_array(non_zero_elements);
id = (1:length(chi))';
nz_id = id(non_zero_elements);
N = length(nz_chi);
id_chi_z_mat = [nz_id,nz_chi,nz_z];
sorted_id_chi_z_mat = sortrows(id_chi_z_mat,2);
sorted_chi=sorted_id_chi_z_mat(:,2);


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


% define the space and time variable As^* matrix
% a row for each pixel
% the columns are ordered according to 1U1 1U2 .... (1U1) is at
% location of pixel 1 at time step 1 (present day)
A = zeros(N,N*q);
for i = 1:N
  
    % global position in term of time.
    % Initialize to the upper bound of the integral for that pixel
    global_time = nz_chi(i);
    time_int = 1;
    % local time is time left for the end of this time interval
    local_time = scaled_dt_vec(time_int);
    id = nz_id(i);
    curr_ind = i;
    %id of the receiver
    id_rec = rec_array(id);
    % row location of the receiver
    next_ind = find(nz_id == id_rec,1);
    while global_time > 1e-12
        %calculate the column in A
        col = (time_int-1)*N + curr_ind;
        if isempty(next_ind)
            chi_next = 0;
        else
            chi_next = nz_chi(next_ind);
        end
        if global_time - chi_next >=local_time+1e-10 %time left is smaller than the difference between chi's
            A(i,col) = local_time;
            global_time = global_time - local_time;
            time_int = time_int+1;
            if time_int <= length(scaled_dt_vec)% just not to go over the limit
                local_time = scaled_dt_vec(time_int);
            else
                'stop here'
            end
            % no need to change curr_ind and next_ind
        else
            A(i,col) = global_time - chi_next;
            local_time = local_time - ...
                (global_time - chi_next);
            global_time = chi_next;
            curr_ind = next_ind;
            id = id_rec;
            id_rec = rec_array(id);
            next_ind = find(nz_id == id_rec,1);
            % no need to change time_int
        end
    end
    if abs(sum(A(i,:))-nz_chi(i))>1e-10
        'problem in filling A'
    end
   
end


%Building the parabila matrix Bp
Bp = zeros(q*N,3*q);


% devision by 1e3 is to convert meters to kms.
for i = 1:q
    for j = 1:N
    
        Bp((i-1)*N+j,(i-1)*3+1) = (nz_x(j)/1e3)^2;
        Bp((i-1)*N+j,(i-1)*3+2) = nz_x(j)/1e3;        
        Bp((i-1)*N+j,(i-1)*3+3) = 1;
    end
end

Ap = A*Bp;


%build the prior model
U_pri_star = ones(3*q,1)*mean(nz_z./nz_chi);

%construct the least square estimate
denom = Ap'*Ap+Gamma^2*eye(3*q);
nom=Ap'*(nz_z-Ap*U_pri_star);
Up = U_pri_star + denom\nom; %size 3qX1
tstar = scaled_t_vec'; %size q+1X1 

Misfit = 1/((N-3*q))*sqrt(sum((nz_z - Ap*Up).^2));


%presenting the results in plots

x_domain = 0:1:max(x)/1e3;
U_space_time_mat = zeros(length(scaled_dt_vec),length(x_domain));
for i = 1:length(scaled_dt_vec)
    U_space_time_mat(i,:) = Up((i-1)*3+1)*(x_domain).^2+...
        Up((i-1)*3+2).*(x_domain)+Up((i-1)*3+3);
end

tstar_plot = [];
Ustar_mat_plot = [];
for i = 1:q
    tstar_plot = [tstar_plot tstar(i) tstar(i+1)];
    Ustar_mat_plot = [Ustar_mat_plot ;U_space_time_mat(i,:) ;U_space_time_mat(i,:)];
end





if to_plot
    
    figure
    [X,T] = meshgrid(x_domain,tstar_plot);
    if K == 1 %non-dimensional inversion
        surf(X,T,Ustar_mat_plot,'edgecolor','none')
        ylabel('scaled t [m]');
        colorbar
        zlabel('Non-dimensional Rock Uplift Rate')
    else
        surf(X,T/K/1e6,Ustar_mat_plot*K/1e-3,'edgecolor','none')
        
        ylabel('time [Ma]');
        colorbar
        zlabel('Rock Uplift Rate [mm/yr]')
    end
    xlabel('x [km]');
    
end









