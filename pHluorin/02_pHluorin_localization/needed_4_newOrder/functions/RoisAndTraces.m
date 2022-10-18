function roi = RoisAndTraces(movie, totNumFrames, syn, nameMov, tailFlag, ignoreTail)
%% nov 2021
%% ignoreTail = tail present (needed to calculate lastFrameBeforeTail), but ignored

% rowNumCorr, columnNumCorr, movie, mask, SigmaFitSyn, smallMovieXY,

%function roi = RoisAndTraces(movie, syn, tailFlag)
% movie, trace, tracesF, traceDt, snr, df_tail

 
if tailFlag(1,1) == 1
    lastFrameBeforeTail = tailFlag(2,1);
    tailFlag = true;
else
    tailFlag = false;    
end



roi = struct;
[~,~, numFramesLoaded] = size(movie);

%%
synID = [syn.synID];%round(([syn.rowNum] .* [syn.columnNum])./([syn.sigma]));
%%
% TailBased depuration: 
if tailFlag == true && ignoreTail == 0
    maskSyn = testTailRoi (syn, totNumFrames, movie, lastFrameBeforeTail);
    syn = syn(maskSyn);
end 



%% if BIG movie
if numFramesLoaded == 1
    numParts = 10;
    lengthPart = floor(totNumFrames / numParts);
    
    for part = 1:numParts
        fstFramePart = lengthPart * (part-1) + 1;
        if part ~= numParts
            lastFramePart = fstFramePart + lengthPart -1;
        else
            lastFramePart = totNumFrames;
        end
        framesPart = (fstFramePart : 1 : lastFramePart);
        [movie, ~ ]= bfOpen3DVolume_GM(nameMov, 0, framesPart);
        finalLengthpart =  length(framesPart);
        for  k = 1:length(syn)
            if part == 1
                roi(k).synID = synID(k);
                roi(k).movie = reshape(movie((syn(k).rowNumCorr -20 : syn(k).rowNumCorr + 20),(syn(k).columnNumCorr - 20 : syn(k).columnNumCorr + 20),:),[41,41, finalLengthpart]);
            else
                 roi(k).movie = cat(3, roi(k).movie, reshape(movie((syn(k).rowNumCorr -20 : syn(k).rowNumCorr + 20),(syn(k).columnNumCorr - 20 : syn(k).columnNumCorr + 20),:),[41,41, finalLengthpart]));
            end
            
            
            
        end
        
        
    end
    
    
end
    

for i = 1:length(syn)
    if numFramesLoaded ~= 1
        %% add ID
        roi(i).synID = synID(i);
        %% movie
        roi(i).movie = reshape(movie((syn(i).rowNumCorr -20 : syn(i).rowNumCorr + 20),(syn(i).columnNumCorr - 20 : syn(i).columnNumCorr + 20),:),[41,41, numFramesLoaded]);
    end
    
    %% trace
    temp = double(syn(i).mask).*double(roi(i).movie(:,:,:));
    area = sum(sum(double(syn(i).mask)));
    roi(i).mask = syn(i).mask;
    roi(i).trace = squeeze(double(sum(sum(temp))/area)); 
    %% traceBckg maybe not needed
    
%     maskBckg = double(~(syn(i).mask));
%     temp = maskBckg.* double(roi(i).movie(:,:,:));
%     area = sum(sum(maskBckg));
%     roi(i).traceBckg = squeeze(double(sum(sum(temp))/area)); 
    
    %% 'filtered' trace
    RoiRadius = ceil(1.5*syn(i).sigma);
    MovingRadius = 3;
    [ ~, ~, roi(i).traceF] = DinamicInt_GroupedPixels(roi(i).movie,  RoiRadius,  MovingRadius, syn(i).smallMovieXY);
     %% detrend an 'un-tailed' trace
     if tailFlag
        %trace2dt = roi(i).trace(1 : lastFrameBeforeTail);
        [roi(i).movie, roi(i).traceDt, ~, ~, ~] = detrendMovie(roi(i).movie, roi(i).mask, lastFrameBeforeTail,0, 0);
        
        %[movie_dt, trace_dtIn, trace_dtOut, trace_dt, templateFactor] = detrendMovie(movie, mask, lastFrameBeforeTail,eqPreFrFlag,outMaskFlag, varargin)

        %% OPTION A
        %% move tail to the level of the Detrended
        tail = roi(i).trace((lastFrameBeforeTail + 1):end);
        % denoise Tail
        tail = denoiseForwrdBack_g(tail,20,10);
        %remove 0 val
        tail(tail == 0) = mean(tail);
        f0 = mean(roi(i).traceDt((lastFrameBeforeTail-10):lastFrameBeforeTail));% last 10 points
        dif = tail(1) - f0;
        tail = tail - dif;
        [~, locMax] = max(tail);
        if locMax < 4
            locMax = 4;
        end
        
        if (locMax + 3) > length(tail)
            locMax = length(tail) - 3;
        end
        f = mean(tail((locMax - 3):(locMax + 3))); %% 6 pnts
        df = (f-f0)/f0;
        roi(i).df_tail = df;
        SDtrace =  std((roi(i).traceDt((lastFrameBeforeTail-10):lastFrameBeforeTail)));
        roi(i).locMaxtail = locMax + lastFrameBeforeTail;
        roi(i).tailSDtrace = f > (f0 + 2*SDtrace);
        
        if roi(i).tailSDtrace == false && ignoreTail == 0
            roi(i).movie = 0;
            roi(i).trace = 0;
            roi(i).traceBckg = 0;
            roi(i).traceF = 0;
            roi(i).traceDt = 0;
            roi(i).traceBckgDt = 0;
            
        
        end
        
        %% OPTION B
        
%         f0 = mean(roi(i).trace((lastFrameBeforeTail - 10):lastFrameBeforeTail));
%         f = mean(roi(i).trace((lastFrameBeforeTail + 8):(lastFrameBeforeTail + 18 )));
%         df = (f-f0)/f0;
%         roi(i).df_tail = df;
%         
     else
         %trace2dt = roi(i).trace;
         [roi(i).movieDt, roi(i).traceDt, roi(i).tracetraceBckgDt, ~, ~] = detrendMovie(roi(i).movie, roi(i).mask, length(roi(i).trace), 0,1);

         %[roi(i).movieDt, roi(i).traceDt, roi(i).tracetraceBckgDt, ~, ~] = detrendMovie(roi(i).movie, roi(i).mask, length(roi(i).trace), 0);
% [movie_dt, trace_dtIn, trace_dtOut, trace_dt, templateFactor] = detrendMovie(movie, mask, lastFrameBeforeTail,eqPreFrFlag,outMaskFlag, varargin)
     end
     %[roi(i).traceDt,~,~] = nt_detrend(trace2dt,50,[],'sinusoids',[],4,[]);
    
    % snr calculated over the detrended trace
    traceUp = ([roi(i).traceDt]-mean([roi(i).traceDt])) + 100;
    roi(i).snr =  mean(traceUp)/std(traceUp);
    
    
    
end
% if isfield(roi,'tailSDtrace')
%     disp([num2str(sum([roi.tailSDtrace])) ' syn. with tail bigger than SD'])

end
