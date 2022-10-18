function syn = synapseDetection(RawData, InfoExp, sigmaRange, AmpTreshold, fixedRadius)  
%% syn with fields : rowNum, columnNum, rowNumCorr, columnNumCorr, SmallMovieXY, Sigma, Bckgrd, Amp, mask

sigma = InfoExp.sigma;

date = strsplit(InfoExp.date,'.');
year = char(date(3));
year = erase(year, '20');
baseID = [char(date(1)) char(date(2)) year erase(InfoExp.name, 'Data')];

disp 'Looking for Synapses'
SigmaMin = sigmaRange(1);
SigmaMax = sigmaRange(2);

if length(AmpTreshold) == 2
    AmpMin = AmpTreshold(1);
    AmpMax = AmpTreshold(2);
else
    AmpMin = AmpTreshold;
    AmpMax = 100; %% absur value, all accepted
    
end


%% 2D GAUSSIAN FIT
[pstruct, ~, ~, ~] = pointSourceDetectionDario(RawData, sigma, 'mode','xyASc'); %-------------------------> DETECTION SYNAPSES!!

%% depurate depurate synapses using the cut offs
LogicVect = pstruct.s'<= SigmaMin | pstruct.s'>= SigmaMax  | pstruct.A'<= AmpMin | pstruct.A'>= AmpMax | pstruct.x'< 10 | pstruct.x'>(size(RawData,1)-10) | pstruct.y'<10 | pstruct.y'>(size(RawData,2)-10);

%% column number
columnNum = round(pstruct.x'); 
columnNum(LogicVect,:) = [];

%% row number 
rowNum = round(pstruct.y');
rowNum(LogicVect,:) = [];

%% sigma
sigma = pstruct.s';
sigma(LogicVect,:) = [];

%% Background
bckgrd = pstruct.c';
bckgrd(LogicVect,:) = [];

%% Amplitude
Amp = pstruct.A';
Amp(LogicVect,:) = [];
%%
syn = struct;

%synID = round((rowNum .* columnNum)./(sigma));
for i = 1:length(columnNum)
    %% add row and column to the syn output
    syn(i).synID = str2double([baseID num2str(i)]);
    syn(i).rowNum = rowNum(i);
    syn(i).columnNum = columnNum(i);
    %% correction factor regarding the border issue
    columnCorrFact = (columnNum(i) - 21 < 0)*(21 - columnNum(i)) + (columnNum(i) + 21 > size(RawData,2))*(-columnNum(i) - 21 + size(RawData,2));
    rowCorrFact = (rowNum(i) - 21 < 0)*(21 - rowNum(i)) + (rowNum(i) + 21 > size(RawData,1))*(-rowNum(i) - 21 + size(RawData,1));
    
    syn(i).rowNumCorr = rowNum(i) + rowCorrFact;
    syn(i).columnNumCorr = columnNum(i) + columnCorrFact;  
    
    smallMovieXYTemp = [21-rowCorrFact 21-columnCorrFact];
    syn(i).smallMovieXY = smallMovieXYTemp; 


    factor = 1+exp(-3*(sigma(i)-SigmaMin));
    r = factor*sigma(i);
    if fixedRadius > 0
        r = fixedRadius;
    end
    syn(i).mask = bsxfun(@plus, ((1:41) - smallMovieXYTemp(2)).^2, (transpose(1:41) - smallMovieXYTemp(1)).^2) < r^2;
    
    %% add the other variables (one by one) to the syn structure, redundace? 
    syn(i).sigma = sigma(i);
    
   % syn(i).synID = round(([syn(i).rowNum] * [syn(i).columnNum])/([syn(i).sigma]));
    
    syn(i).Amp = Amp(i);
    syn(i).bckgrd = bckgrd(i);
    
end
disp([num2str(length(syn)) ' synapses located'])
end
