function [trans1, trans2] = getTransversal(cx,cy, img, varargin)

if isempty(varargin)
    locPeak = round(length(img)/2);
else
    locPeak = varargin{1};
end
limit= length(img);

l1 = min((cx -1), (cy - 1));
l2 = min((limit - cx), (limit - cy));
l3 = min((cx - 1), (limit - cy));
l4 = min((limit - cx), (cy - 1));


%% transversal 1
indxTrans1x = ((cx-l1):1:(cx + l2)); 
indxTrans1y = ((cy-l1):1:(cy + l2));
indxTrns1= [indxTrans1y; indxTrans1x]';
trans1 = zeros(length(indxTrns1),1);
for i = 1:length(indxTrns1)
    trans1(i,1) = img(indxTrns1(i,1),indxTrns1(i,2));
end


%% transversal 2
indxTrans2x = ((cx-l3):1:(cx + l4)); 
indxTrans2y = ((cy + l3):-1:(cy - l4));
indxTrns2= [indxTrans2y; indxTrans2x]';
trans2 = zeros(length(indxTrns2),1);
for i = 1:length(indxTrns2)
    trans2(i,1) = img(indxTrns2(i,1),indxTrns2(i,2));
end

if ~isempty(varargin)
%     tail1 = ones(61,1).*mean([trans1(1:5); trans1(length(trans1)-4: end)]);
%     tail2 = ones(61,1).*mean([trans2(1:5); trans2(length(trans2)-4: end)]);
%     lPart = varargin{1}-1;
%     tail1(31-l1:((31-l1) + length(trans1) -1)) = trans1;
%     tail2(31-l3:((31-l3)+length(trans2)-1)) = trans2;
%     trans1 = tail1((31-lPart):(31 + lPart));
%     trans2 = tail2((31-lPart):(31 + lPart));

    tail1 = ones(81,1).*mean([trans1(1:5); trans1(length(trans1)-4: end)]);
    tail2 = ones(81,1).*mean([trans2(1:5); trans2(length(trans2)-4: end)]);
    lPart = varargin{1}-1;
    tail1(41-l1:((41-l1) + length(trans1) -1)) = trans1;
    tail2(41-l3:((41-l3)+length(trans2)-1)) = trans2;
    trans1 = tail1((41-lPart):(41 + lPart));
    trans2 = tail2((41-lPart):(41 + lPart));

end

%% correct pixels
ideaLength = length(trans1);
actuaLength =  sqrt(2)* ideaLength;
x1 = ((1:length(trans1)));
x2 = (1:(ideaLength/actuaLength):length(trans1));

interTrans1a = (interp1(x1,trans1,x2, 'spline'))';
interTrans1b = flip((interp1(x1,flip(trans1),x2, 'spline'))');
% interTrans1c = (interp1(x1,trans1,x2))';
% interTrans1d = flip((interp1(x1,flip(trans1),x2))');


interTrans1 = (interTrans1a + interTrans1b)./2;

interTrans2a = (interp1(x1,trans2,x2, 'spline'))';
interTrans2b = flip((interp1(x1,flip(trans2),x2, 'spline'))');
interTrans2 = (interTrans2a + interTrans2b)./2;

pix2rmv = round((length(interTrans1) - length(trans1))/2);

interTrans1(1:pix2rmv) = [];
interTrans1((length(interTrans1)-(pix2rmv-1)):end) = [];

interTrans2(1:pix2rmv) = [];
interTrans2((length(interTrans2)-(pix2rmv-1)):end) = [];
%% give it same amplitude 
trans1 = interTrans1 ;%+ (trans1(locPeak) - interTrans1(locPeak));
trans2 = interTrans2 ;%+ (trans2(locPeak) - interTrans2(locPeak));

end