function chi = CalculateChi(x,y,rec_array,area_array,m)
%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 07/11/2019 %%%%%%%%%%%%%%%%
% Function to calculate the chi for river pixels of a fluvial network.
% Input parameters:
% x - vector of size n of the x coordinate [L] of each pixel.
% y - vector of size n of the y coordinate [L] of each pixel.
% rec_array - vector of size n of the reciever relation. A value j in row i
%             means that the reciever of the pixel in row i 
%            (for which x, y, and drainage area are supplied) is in row j. 
%             The reciever of the outlet is the pixel itself.
% area_array - vector of size n of upstream drainage area [L^2] 
%              along fluvial pixels.
% m - scalar. The area power in the stream power law. The analysis assumes 
%     that n = 1.   
% OUtput:
% chi - vector of size n with the chi values [L] of the pixels. 
%
% Note: row i in x,y,area_array and chi refers to the same pixel.




A0 = 1; %reference drainage area

%The easisest way to calculate chi is to sort the network from outlet to
%heads. For that all the data is grouped and sorted based on drainage area.
%such that we start with the largest area.

DataMat = [(1:length(x))',x,y,area_array];
SortedDataMat = sortrows(DataMat,-4);

chi = zeros(size(x));
for i = 1:length(x)
    my_id = SortedDataMat(i,1);
    my_rec = rec_array(my_id);
    if my_id==my_rec % outlets are defined as their own receivers
        chi(i) = 0;
    else
        j = find(SortedDataMat(:,1) == my_rec); % find the row of my receiver
        dist = sqrt((SortedDataMat(i,2)-SortedDataMat(j,2))^2 +...
            (SortedDataMat(i,3)-SortedDataMat(j,3))^2); 
        chi(i) = chi(j) + dist*(A0/SortedDataMat(i,4))^m;
    end
end

SortedDataMat = [SortedDataMat chi]; %add chi to the data
DeSortedDataMat = sortrows(SortedDataMat,1); %resort based on ids
chi = DeSortedDataMat(:,end);