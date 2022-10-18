function trace = traceInt(movie, x, y,frameResponse,rmvBckgFlag,foldFlag, varargin)

locPeak = 16; %% 16 in transTrace,  where the peak in the trace should be
trimmerFlag = 1; 

if rmvBckgFlag == 1
    NumFrame2Remv = 6;% 6
    if (frameResponse - NumFrame2Remv ) <= 0
        NumFrame2Remv = frameResponse -1;
    end
    bckg = double(median(movie(:, :, (frameResponse - NumFrame2Remv):(frameResponse - 2)), 3));
    movie = double(movie) - double(bckg);    
    
end

%% square frames
[Lframe1,~, Lmovie]= size(movie);

if isempty(varargin)
    numFramBef = frameResponse -1;
    numFramAft = Lmovie -frameResponse;
else
    numFramBef = varargin{1};
    numFramAft = varargin{2};
    
end

indxMov = (frameResponse - numFramBef:frameResponse + numFramAft);
indxMov(indxMov < 1) = 1;
indxMov(indxMov > Lmovie) = Lmovie;

limsx = (x - (locPeak - 1)):1:(x + (locPeak - 1));
limsy = (y - (locPeak - 1)):1:(y + (locPeak - 1));

limsx(limsx < 1) = 1;
limsx(limsx > Lframe1) = Lframe1;
limsy(limsy < 1) = 1;
limsy(limsy > Lframe1) = Lframe1;

%% x and y traces at the center
trace = zeros(length(indxMov), 1);
pnt = 1;
for fr = indxMov
    
    frame = movie(:,:,fr);
    xtrc = double(frame(limsy,x));
    ytrc = double((frame(y,limsx))');
    [t1trc, t2trc] = getTransversal(x,y, frame, locPeak);
    if foldFlag == 1
        xtrc = [(xtrc(1:15) + flip(xtrc(17:end)))./2; xtrc(16)];%;flip((TotTrace(1:15) + flip(TotTrace(17:end)))./2)];
        ytrc = [(ytrc(1:15) + flip(ytrc(17:end)))./2; ytrc(16)];
        t1trc = [(t1trc(1:15) + flip(t1trc(17:end)))./2; t1trc(16)];
        t2trc = [(t2trc(1:15) + flip(t2trc(17:end)))./2; t2trc(16)];
        
    end
    
    
    TotTrace = (xtrc + ytrc + t1trc + t2trc) ./ 4;
    
    if trimmerFlag ==1
        borderL = 11;%
        trimmer = true(size(xtrc));
        trimmer(1:borderL) = false;
        trimmer(length(trimmer)-borderL+1:end) = false;
    else
        trimmer = true(size(xtrc));

    end
    
    if foldFlag == 1
        %TotTrace = [(TotTrace(1:15) + flip(TotTrace(17:end)))./2; TotTrace(16)];%;flip((TotTrace(1:15) + flip(TotTrace(17:end)))./2)];
       trimmer = trimmer(1:16);
    end

    
    trace(pnt,1)= trapz(TotTrace(trimmer));
    pnt = pnt +1;
    
    
end
