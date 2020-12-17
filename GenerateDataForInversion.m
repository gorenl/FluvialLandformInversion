
%%%%% Liran Goren, gorenl@bgu.ac.il, December 2020 %%%%%%%%
%
%
%
% Script to generate data for the inversion code.
% Calls to RecArrayFromStreamObj, so assumes 
% original data from topotoolbox. 
% See RecArrayFromStreamObj for more details



pixelsize=10; %m %your pixel size
pixeltoarea=pixelsize^2; %m^2
m_of_concavity=0.45; %your coice of the m exponent. 

%calling the script GenerateData
RecArrayFromStreamObj;
% The script produces TableData with X|Y|Z|A|Flow direction matrix
% and a rec_array vector

%generate the results matrix 
%|ID|X|Y|Z|A|Flow dist|Dist to receiver|Chi|
r = length(rec_array);
DataTauMat=zeros(r,8);
DataTauMat(:,1)=(1:r);
DataTauMat(:,2:6)=TableData(:,1:5);
DataTauMat(:,5)= DataTauMat(:,5)*pixeltoarea;
%distance between donor-receiver
for i =1:r
    j=rec_array(i);
    if j~=0
        DataTauMat(i,7)= DataTauMat(i,6) - DataTauMat(j,6);
    end
end
DataTauMat = sortrows(DataTauMat,6); %sort by area
% now we can scan from bottom to top
% calculate chi
for i=2:r
    me = DataTauMat(i,1);
    my_rec=rec_array(me);
    j = find(DataTauMat(:,1) == my_rec);
    if j ~= 0
        DataTauMat(i,8) = DataTauMat(j,8) +...
            DataTauMat(i,6)/DataTauMat(i,5)^m_of_concavity;
    end
end

figure; %for debugging
scatter(DataTauMat(:,2),DataTauMat(:,3),20,DataTauMat(:,8));
colormap('jet')

DataTauMat = sortrows(DataTauMat,1);

ElevationChi(:,1)=DataTauMat(:,4);
ElevationChi(:,2)=DataTauMat(:,8);
save('ElevationChi','ElevationChi','-ascii');
save('AllData','DataTauMat','-ascii');
save('river_network','rec_array','-ascii');



 