function [ Fil, Col, Int ] = DinamicInt_GroupedPixels(RawData, RoiRadius, MovingRadius, CenterPoint)
%make a Intensity string of a movie for the RoiRadius. This Roi is Placed
%in the position that get the Maximum Int average

RawData = double(RawData);

%% 
radius1=RoiRadius; % radius of internal circle to cuantificate
radius2=MovingRadius; % radius to area to evaluate.

%% Make Mean Filter

Disk1 = fspecial('disk',radius1)>0;
Filter = Disk1./sum(sum(Disk1));

%% Make Mask

[ sizeX, sizeY, sizeZ ] = size(RawData);

if nargin < 4
    centre=[ceil(sizeX/2) ceil(sizeY/2)];
else
    centre=CenterPoint;
end

Mask=zeros(sizeX);
Disk = fspecial('disk',radius2)==0;
Mask(centre(1)-radius2:centre(1)-radius2+size(Disk,1)-1, centre(2)-radius2:centre(2)-radius2+size(Disk,2)-1)=double(~Disk);

% Apply Mask

[Int,FilTemp] = max(bsxfun(@times,imfilter(RawData,Filter),Mask));
[Int,Col] = max(Int);
Col = squeeze(Col);
Int = squeeze(Int);
FilTemp = squeeze(FilTemp);

for i = length(Col):-1:1
    
    Fil(i) = FilTemp(Col(i),i);
    
end

Fil = Fil';

end

    
