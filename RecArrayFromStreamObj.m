
%%%%% Liran Goren, gorenl@bgu.ac.il, December 2020 %%%%%%%%
%
%
%
%%%% The script uses topotoolbox objects to generate 
%%%% a receiver array vector 
%%%% and a matrix of X|Y|Z|A|Flow distance| 
%%%% The receiver array is needed for 
%%%% the fluvial inversion algorithms. 
%%%% So far the script was tested only once on data 
%%%% sent to me by Simone Racano.
%%%% The test data (not available) was of a single basin
%%%% if more than stream network is made of more than one
%%%% basin, then code adjustments are probably necessary.
%%%% The script uses topotoolbox funcationality

% object_data assumes GridObj: DEM and A, Streamobj: S 
load('object_data.mat')


[X,Y,Z,A,dist] = STREAMobj2XY(S,DEM,A,S.distance);
num_nodes = length(X);
rec_array = zeros(num_nodes,1);
k = 1;
TableData = zeros(num_nodes,5);
first_nan = 1;
for i = 1:num_nodes  
    if(isnan(X(i)))
        if first_nan
            k = k-1;
            first_nan = 0;
            rec_array(k) = 0; %the outlet
        else
            k = k-2;
            junction = find(TableData(:,1) == X(i-1) & TableData(:,2) ==Y(i-1));
            rec_array(k) = junction(1);
            rec_array(k+1) = 0;
            TableData(k+1,:) = [0 0 0 0 0];
        end
    else
        TableData(k,:) = [X(i) Y(i) Z(i) A(i) dist(i)];
        rec_array(k) = k+1;
    end
    k = k+1;
end
first_zero = find(TableData(:,1) == 0,1);
TableData = TableData(1:first_zero-1,:);
rec_array = rec_array(1:first_zero-1);

% For debugging.  
% If the the figure reproduces well the drainage netwrok, 
% it means that the receiver array is valid.
figure; 
hold on;
for i = 1:first_zero-1
     my_rec = rec_array(i);
     if (my_rec ~=0)
     plot([TableData(i,1) TableData(my_rec,1)],[TableData(i,2) TableData(my_rec,2)],'b');

     end
end