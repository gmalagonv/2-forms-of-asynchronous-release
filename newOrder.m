%% NEW ORDER : DETECT EVENTS IN EVERY FRAME
function newOrder(MovieN, info, binFrames,OnlypositDfFlag)%, selectedWindowLength)%, FlagLight,FlagTtest)
%% MovieN is ordinal 

tic;
%% SOME FLAGS & VARIABLES
selectedWindowLength = 0; %% if 0, use 120 | if 1, use ALL (meaning, total number of APs delivered)| if other num, use that number.    
memorySaveFlag = 0;
radiusMask = 0; %if 0: radius is of the syn mask is autocalculated; else: radius used = radiusMask.
minNumSynch = 5; % default value: 5. Minimun number of synchronous events to be cosnidered a valid synapse;
ignoreTail = 0;
%% LOCATION DATA
if strcmp(computer,'PCWIN64') == 1 % GM, apr 2019, if working in windows
    cd D:/Gerardo/data
    strHome = 'D:/Gerardo/nextcloud/analysis';    
else
    if "chimbita" == getComputerName()
        cd /media/Data
    else
        cd ~/data        
    end    
   strHome = '/home/gerard/nextcloud/analysis';
end

%% Get info Experiment.
names = fieldnames(info);
eval ([sprintf('Data.InfoExp = info.%s;', names{MovieN})]);
disp(names{MovieN});
date = strsplit(Data.InfoExp.date,'.');
year = char(date(3));
cd(char(fullfile(year,Data.InfoExp.date)));   %GM dec.2018

%% load the WHOLE movie ***************************************************
% delete previous flags regarding 'TAIL' and 'descartFrame'
[~, firstPreframe] = max(Data.InfoExp.descartFrame);
%Data.InfoExp.descartFrame = [];
Data.InfoExp.OnlyStimFrames = false;
Data.InfoExp.tail = false;
%eval ([sprintf('Movie = bfOpen3DVolume_GM(''%s.cxd'', Data.InfoExp);',names{MovieN})]);
%Andor iXon camera or Orca Hamamatsu

%% NOTE:
% the frames4synDetect can change for ex. if we need only stim frames etc. for
% example, we could not average them and load only the stim events etc...: 
%  if infoExp.OnlyStimFrames == true
%         if infoExp.tail == true %% remove tail from discrete selection of frames
%             numImages = length(infoExp.descartFrame)*(infoExp.NumberAP);
%             imagesIndexList = (1:numImages);
%         end        
%         logicVect = logical(repmat(infoExp.descartFrame,1,round(numImages/length(infoExp.descartFrame))));
%         imagesIndexList = imagesIndexList(logicVect);
%         numImages = length(imagesIndexList);               
%     end
 
if isfield(Data.InfoExp, 'FileExtension') && strcmp(Data.InfoExp.FileExtension,'.tif') 
    ext = '.tif';
else
    ext = '.cxd';
end


nameMov = [names{MovieN} ext];
%% check if file already concatenated
[concatpresent] = checkIfConcInFile(nameMov);
if concatpresent
    disp('file concatenated')
    %% MEMORY SAVE
    if memorySaveFlag == 1
        [~,totNumFrames] = bfOpen3DVolume_GM(nameMov,1, (1:1:2));
        frames2load = (1:1:totNumFrames);
        avgFlag = 1;
    else 
        frames2load = 1; % 1 means all
        avgFlag = 0;
    end

    %%
    [Movie,totNumFrames] = bfOpen3DVolume_GM(nameMov, avgFlag, frames2load, 0);
    %[volume, numImages] = bfOpen3DVolume_GM(filename, avgFlag, frames2load, varargin)
else
    disp('file NO concatenated, load by parts')
    [Movie,totNumFrames] = movieMergerMTLB(nameMov);
    
    
end


%%*************************************************************************
%% Get Z projection & normalize
[~,~,lMovie] = size(Movie);
if lMovie > 1
    RawData = AverageZDimension(Movie);
else
    RawData = Movie;
end
    
RawData = (RawData-min(RawData(:)))/(max(RawData(:))-min(RawData(:))); % normalization

%% synapse selection parameters for sigma: [162.5 520] in nm
Data.AnalizeSetups.LookingSyn.SigmaRange = [162.5 520]./Data.InfoExp.PixelSize; % in pixels
Data.AnalizeSetups.LookingSyn.AmplitudTreshold = 0.02;
if isfield(Data.InfoExp, 'FileExtension') && strcmp(Data.InfoExp.FileExtension,'.tif') 
    Data.AnalizeSetups.LookingSyn.AmplitudTreshold = [0.02 0.5];
end
%sigmaSyst = Data.InfoExp.sigma;

%% finding SYNAPSES in the Z-projected image
% syn with fields: rowNum, columnNum, rowNumCorr, columnNumCorr, SmallMovieXY, Sigma, Background, A, Mask
syn = synapseDetection(RawData, Data.InfoExp, Data.AnalizeSetups.LookingSyn.SigmaRange, Data.AnalizeSetups.LookingSyn.AmplitudTreshold, radiusMask);

%% remove redundant locs, min dist == 10 pix, just to avoid to show up in the same frame
logicCorr = comparecoor_corrB([syn.rowNum], [syn.columnNum] , 10);
syn = syn(logicCorr);
disp([num2str(sum(~logicCorr)) ' redundant locations, ' num2str(length(syn)), ' synapses left'])
%% ADD the goodness of fit as a criteria to choose the best
%% determinate the presence of a Tail 
%[ ~, ~, totNumFrames ] = size(Movie);
lastFrameBeforeTail = ((Data.InfoExp.NumberAP)*((1/Data.InfoExp.Frequency)/Data.InfoExp.timeframe)) + firstPreframe;
if (lastFrameBeforeTail + 5) < totNumFrames
    tailFlag(1,1) = 1;
    tailFlag(2,1) = lastFrameBeforeTail;
else
    tailFlag = 0;
end

%% if the number of synapses is not huge, Load the whole movie to make the roi extraction faster
if length(syn) < 5 && memorySaveFlag == 1
    [Movie,~] = bfOpen3DVolume_GM(nameMov,0, 1);
end
%% xtract small movies and traces froms ROIS
disp 'extracting traces and movies from ROIs'
roi = RoisAndTraces(Movie, totNumFrames,syn, nameMov ,tailFlag, ignoreTail);

%% remove unnecesary fields
syn = rmfield(syn,'mask'); %% added to roi
%roi = rmfield(roi,'movie');
disp 'done'

%% depurate 'syn' and 'roi' based on Tail if so
if ignoreTail == 0
    if tailFlag(1,1) == 1 
        %positiveTail = ([roi.df_tail] > 0)  ;
        positiveTail = ([roi.tailSDtrace] == 1)  ;
        roi = roi(positiveTail); %% remove all the negative DF.
        syn = syn(positiveTail);
        disp([num2str(sum(~positiveTail)) ' tailless synapses, ' num2str(length(syn)), ' synapses left'])

    end
else
    if tailFlag(1,1) == 1 
        Data.InfoExp.ignoreTail = true;
        disp('tail ignored')
    end
end

%% 'good' synapses have at least 50 % of their after stims frames (5) a positive df
[~, msk, ~] =  testSimpleDf(roi, lastFrameBeforeTail,0);
roi = roi(msk);
syn = syn (msk);
disp([num2str(sum(~msk)) ' synapses removed by DF, ' num2str(length(syn)), ' synapses left'])
%% events Localization **************************************** EV LOC ************************************************                 <----------------------
disp ('events Localization')
%locEv = locEvs(syn, roi,1, 0,Data.InfoExp);
locEv = locEvs(syn, roi,1, 0,Data.InfoExp);
%locEvs(syn, roi,fixBckFlag,olDetectprogFlag, InfoExp)

%% clasify Events
if ~isempty(fields(locEv))
    locEv  = clasifyEvents(locEv, Data.InfoExp);
else
    disp('nothing found, no saved')
    return
end

%% depurate syn and roi
syn = syn(ismember([syn.synID], unique([locEv.synID])));
roi = roi(ismember([roi.synID], unique([locEv.synID])));

%% dfCriteria
%% testing shit. nov 2021 
% locEv = traceSmallmaskPre(Data.InfoExp,roi, locEv ,5, 0, 0);
% locEv = transTraceIntPre(roi,locEv,0,0, 0);
%locEv = dfCriteriaB (locEv);
%%

locEv = dfCriteria (locEv, roi, Data.InfoExp, OnlypositDfFlag, 1);
%locEv = dfCriteria (locEv, roi, InfoExp,flagDepurate, flagCrit);

% apply df criteria to the non Synch events Only
%% following 3 lines commented out by nov 2020 
% if OnlypositDfFlag ~= 1
%     %locEv = locEv(logical([locEv.synch]+([locEv.df] > 0)));
%     locEv = locEv(logical([locEv.synch]+([locEv.dfCrit])));
% end
    
%% SDcriterias
locEv = SDcriteria (locEv, roi, 0);
locEv = SDcriteriaB (locEv, roi, 0);
locEv = SDcriteriaC (locEv, roi, 0, Data.InfoExp);
locEv = SDcriteriaBpreF (locEv, roi, 0);
% apply the SD criteria to the non synch events
%% following line commented out by nov 2020
%locEv = locEv(logical([locEv.synch]+([locEv.SDa])));

%% Synch events perSynapse
if selectedWindowLength ~= 0
    if selectedWindowLength == 1
        locEv = bestWin2 (locEv, minNumSynch, Data.InfoExp.NumberAP, Data.InfoExp, 1);
        Data.InfoExp.selectedWindowLength = Data.InfoExp.NumberAP;
    else
        locEv = bestWin2 (locEv, minNumSynch, selectedWindowLength, Data.InfoExp, 1);
        Data.InfoExp.selectedWindowLength = selectedWindowLength;
    end
    %% DEFAULT OPTION:
else 
    locEv = bestWin2 (locEv, minNumSynch, 120, Data.InfoExp, 1);
    Data.InfoExp.selectedWindowLength = 120;
end
Data.InfoExp.minNumSynch = minNumSynch;
%winlength = 120;
%% depurate syn and roi
syn = syn(ismember([syn.synID], unique([locEv.synID])));
roi = roi(ismember([roi.synID], unique([locEv.synID])));

% %% Add frameRelative & stimNum
% fr = [locEv.frame]' - firstPreframe;
% frR = (fr-(floor(fr/20))*20);
% frR(frR == 0) = 20;
% st = (floor((fr-1)./20)) + 1;
% for ff = 1:length(locEv)
%     locEv(ff).frameRelat = frR(ff);
%     locEv(ff).stimNum = st(ff);
% end
%% Add dt trace to each event
[locEv, syn] = pieceOftrace2(locEv,roi, syn, Data.InfoExp);

%% clusters
[syn, locEv] = allTypeOfClusters(syn,locEv, Data, 1);

%% distance to the center
[locEv, syn] = distances2centers(locEv, syn, Data.InfoExp);

%%
Data.perSyn = syn;
Data.roi = roi;
Data.allEv = locEv;
%% distances and time  from prev Ev
Data = distFromLastEv(Data);
Data = timeFromLastEv(Data);

%% multifunctional function with some useful stuff (2021)
Data = normByRadius2(Data, 1, Data.InfoExp.PixelSize);

%% Data = transTrace(DataName,rmvBckgFlag, mskFlag,recenterFlag, transvFlag, textInput)
Data = transTrace(Data,1,0,0,0,1);
%% small mask
Data = traceSmallmask(Data, 5, 0, 0);
%% int Trace
Data = transTraceInt(Data,0,0,0);

%% remove repeated IDs (if so..)
Data = corrDuplisynIDS(Data);
%% march 2022, add some fields
Data = normByRadius3(Data);
Data = addFieldSyn2dEvs(Data, 'areaSynch');
Data = addSomeFitVariables(Data);
if length(Data.allEv) >= 1
    Data = detectDoublets(Data,'traceDt');
end

%% ASSIGN TO THE WORKSPACE
%assignin('base', 'Data',Data);
%%

cd (strHome)
cd (year)

 disp([num2str(length(Data.perSyn)) ' synapse(s) with ' num2str(length(Data.allEv)) ' events located'])
if ~isempty(Data.allEv)
    if binFrames == 1
        %fname = sprintf('spH_%s_%s', Data.InfoExp.date, names{MovieN});
        fname = sprintf('NOr_%s_%s', Data.InfoExp.date, names{MovieN});
    else
        %fname = sprintf('spH_%s_%s_%s', Data.InfoExp.date, names{MovieN}, ['_binFrame' num2str(binFrames)]);
        fname = sprintf('NOr_%s_%s_%s', Data.InfoExp.date, names{MovieN}, ['_binFrame' num2str(binFrames)]);
    end
    if tailFlag(1,1) == 1 && ignoreTail == 1
        fname = [fname, '_IgTail'];
    end
    
    if Data.InfoExp.selectedWindowLength ~= 120
        fname = [fname, '_win' num2str(Data.InfoExp.selectedWindowLength)];
    end
    

    if isfolder(Data.InfoExp.date) == 0    %GM dec.2018
        mkdir (Data.InfoExp.date);
    end    
    eval ([sprintf('save %s/%s.mat Data ''-v7.3''', Data.InfoExp.date, fname)]);
    disp 'saved'
else
    disp('nothing found, no saved')

   
end
disp('Done! clap your hands and say...')
soundEnd;
clearvars

toc;




end
