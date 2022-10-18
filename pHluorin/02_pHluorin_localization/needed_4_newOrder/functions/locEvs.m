function locEv = locEvs(syn, roi,fixBckFlag,olDetectprogFlag, InfoExp)
tic;
sigma = InfoExp.sigma;

%%
locEv = struct;
[~, firstPreframe] = max(InfoExp.descartFrame);
%numAP = InfoExp.NumberAP;

%% actually totNumFrames = lastFrameBeforeTail
totNumFrames = ((InfoExp.NumberAP)*((1/InfoExp.Frequency)/InfoExp.timeframe)) + firstPreframe;
totNumSyn = length(syn);

%% fit parameters
A = cell(totNumFrames,totNumSyn);
x = A;
y = A;
s = A;
c = A;
val = A;
%% standard deviations, estimated by error propagation
A_pstd = A;
x_pstd = A;
y_pstd = A;
c_pstd = A;

%% synapse info
%synNums = A;
synID = A;
synSigma = A;
%% frames
%frames = A;
frames = A;
%% logical Mask
selector = false(totNumFrames,totNumSyn);

for synNum = 1:totNumSyn    
    movie = roi(synNum).movie;    
    mask =  roi(synNum).mask;
    ID = syn(synNum).synID;%round((syn(synNum).rowNum .* syn(synNum).columnNum)./(syn(synNum).sigma));
    ssigma = syn(synNum).sigma;
    %%
        %[maskDT, ~] = calculateMask(movie, mask, sigma);
    %[movie, ~, ~] = detrendMovie(movie, mask, 1, 1);
    
    if olDetectprogFlag == 1 
        premovie = cat(3, movie(:,:,1), movie);
        %[premovie,~ ]= avgPreframes(movie, 5,0);
        %[bck, ~] = localBck(5,200,movie);

        %preframe = avgPreframes(movie, numPrevFr, meanFlag)
        threshEv = 2.5;
        
    end
    
    if fixBckFlag > 0
        if fixBckFlag == 1
            %[bck, ~] = localBck(5,movie, InfoExp);
            [bck, ~] = localBck2(5,1,movie, InfoExp);
            %[bck, ~] = localBck3(movie,InfoExp,1);
% % %             
%             for fr = 1:length(movie)
%                 bck(:,:,fr) = medfilt2(bck(:,:,fr));
%             end
            
        elseif fixBckFlag == 2
            %[bck,~ ]= avgPreframes(movie, 5,0);
            %% other option, just preprame filtered
            premovie = cat(3, movie(:,:,1), movie);
            for fr = 1:length(movie)
                bck(:,:,fr) = medfilt2(premovie(:,:,fr));
            end
                
   
            
        elseif fixBckFlag == 3
             premovie = cat(3, movie(:,:,1), movie);
             %[premovie,~ ]= avgPreframes(movie, 2,1);
             [bck, ~] = localBck(5,movie, InfoExp);
             
             for fr = 1:length(movie)
                 bck(:,:,fr) = bck(:,:,fr).*(1/3) + premovie(:,:,fr).*(2/3);
             end
        else
            [bck, ~] = avgPreframesB(movie, 0, 1);
%             [bck, ~] = localBck(5,movie, InfoExp);
%              
%              for fr = 1:length(movie)
%                  bck(:,:,fr) = double(bck(:,:,fr)).*(1/3) + premovie(:,:,fr).*(2/3);
%              end

            
        end
    else
        bck = 0;
    end
    
  %%    test  
    %% select a subset of frames, not needeed to localize in all
    [~, ~, framesS] =  testSimpleDf(roi(synNum), totNumFrames, 1);
    %[~, ~, framesS] =  testSimpleDfB(roi(synNum),firstPreframe, totNumFrames);
    framesS = cell2mat(framesS);
    framesS = SDcriteriaBpreF_trace (roi(synNum).traceDt, framesS);
    %framesS = (1:totNumFrames);
  
    %% compare bckg with trace (trace NEEDS to be bigger than bckg), another criteria to avoid to check unnecesary frames.
    if  bck ~= 0
        temp = (mask).* double(bck);
        area =  sum(sum(mask));
        bck_tr =  squeeze(double(sum(sum(temp))/area));
        trace = roi(synNum).traceDt; 
        mskPos = (trace(framesS)-bck_tr(framesS)) > 0;
        framesS = framesS(mskPos);
    end
    %% apply selection
    movie = movie(:,:,framesS);
    if bck ~= 0
        bck = bck(:,:,framesS);
    end
    if olDetectprogFlag == 1
        premovie = premovie(:,:,framesS);
    end
    
    
    %% ************************* PAR FOR LOOP HERE
    parfor frame = 1:length(framesS)
%         
%          [pstruct, ~, ~, ~] = FusionLoc(movie(:,:,frame), premovie(:,:,frame), sigma,'mask', double(mask), 'DeltaAcTresh', threshEv,'OnlyOne', true, ...
%                                         'imposedBackground',bck, 'mode','xycA');
        [pstruct, ~, ~, ~, ~ ] = FusionLoc_oneFrame_GM(movie(:,:,frame), sigma,'mask', double(mask),'OnlyOne', true, 'mode','xyA', 'imposedBackground', bck(:,:,frame));



        if ~isempty(pstruct) %&& pstruct.A_pstd <= 100 && pstruct.s >= 1 && pstruct.s <= 2
            
            A{frame, synNum} = [pstruct.A];
            x{frame, synNum} = [pstruct.x];
            y{frame, synNum} = [pstruct.y];
            s{frame, synNum} = [pstruct.s];
            c{frame, synNum} = [pstruct.c];
            val{frame, synNum} = [pstruct.valA];
             
            
            A_pstd{frame, synNum} = [pstruct.A_pstd];
            x_pstd{frame, synNum} = [pstruct.x_pstd];
            y_pstd{frame, synNum} = [pstruct.y_pstd];            
            c_pstd{frame, synNum} = [pstruct.c_pstd];
            
            
            % sigma_r : standard deviation of the background (residual)
%         SE_sigma_r : standard error of sigma_r
%            pval_Ar : p-value of an amplitude vs. background noise test (p < 0.05 -> significant amplitude)
            
            
            
            
            synID{frame, synNum} = ones(size([pstruct.A])).*ID;
            synSigma{frame, synNum} = ones(size([pstruct.A])).*ssigma;
            %synNums{frame, synNum} = ones(size([pstruct.A])).*synNum;            
            %frames{frame, synNum} = ones(size([pstruct.A])).*frame;
            frames{frame, synNum} = framesS(frame);
            
            selector(frame, synNum) = true;
    
        end    
    end
    
    
end

A = A(selector);
maskLongOnes = locateLongs (A);
A = [A{:}]; 

x = [x{selector}];
y = [y{selector}];
s = [s{selector}];
c = [c{selector}];
val = [val{selector}];

A_pstd = [A_pstd{selector}];
x_pstd = [x_pstd{selector}];
y_pstd = [y_pstd{selector}];
c_pstd = [c_pstd{selector}];

%synNums = [synNums{selector}];
synID = [synID{selector}];
synSigma = [synSigma{selector}]; 
%frames = [frames{selector}];
frames = [frames{selector}];
%% separate events where more than one was found


%% FILL UP STRUCTURE
for k = 1:length(A)
    locEv(k).synID = synID(k);
    %locEv(k).synNum = synNums(k);
    locEv(k).synSigma = synSigma(k);
    
    %locEv(k).frame = frames(k);
    locEv(k).frame = frames(k);
    
    locEv(k).A = A(k);
    locEv(k).x = x(k);
    locEv(k).y = y(k);
    locEv(k).s = s(k);
    locEv(k).c = c(k);
    locEv(k).val = val(k);
    
    
    locEv(k).A_pstd = A_pstd(k);
    locEv(k).x_pstd = x_pstd(k);
    locEv(k).y_pstd = y_pstd(k);
    locEv(k).c_pstd = c_pstd(k);
    locEv(k).long = maskLongOnes(k);
end
%% remove some multiple locs per frame by 'continuity'
if sum(maskLongOnes) > 1
    locEv = checKContinuityLongsB(locEv);
end


%% Add frameRelative & stimNum
if isfield(locEv,'frame')
    fr = [locEv.frame]' - firstPreframe;
    frR = (fr-(floor(fr/20))*20);
    frR(frR == 0) = 20;
    st = (floor((fr-1)./20)) + 1;
    for ff = 1:length(locEv)
        locEv(ff).frameRelat = frR(ff);
        locEv(ff).stimNum = st(ff);
    end
    %% remove the rest of multiple locs selecting the best fit
    disp ('bestFitLongs started')
    locEv = bestFitLongs(locEv);

end


toc;

end
