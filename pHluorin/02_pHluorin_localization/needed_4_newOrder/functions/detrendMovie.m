function [movie_dt, trace_dtIn, trace_dtOut, trace_dt, templateFactor] = detrendMovie(movie, mask, lastFrameBeforeTail,eqPreFrFlag,outMaskFlag, varargin)

%% varargin{1} = mask claculated by [maskMovie, NUmcalcMask] = calculateMask(movie, excluMask, sigma)
%% if mask input is set to 1, the varargin{1} will also be used to get the trace, if not, it will use the mask one (added to mask2)
%% additonal mask 2 plot, used to generate the traces
w = 0;

%% 
if ~isempty(varargin)
    mask2 = mask;
    mask = varargin{1};
end


%% nov 2020
%% calculate lastFrameBeforeTail
% numStims = 200;
% freqStims = 1;
% timeframe = 0.05;
% firstPreframe = 4;
% lastFrameBeforeTail = ((numStims)*((1/freqStims)/timeframe)) + firstPreframe; %lastFrameBeforeTail = ((Data.InfoExp.NumberAP)*((1/Data.InfoExp.Frequency)/Data.InfoExp.timeframe)) + firstPreframe;
% 

movie = movie(:,:,(1:lastFrameBeforeTail));
[~, ~ ,framesMsk]= size(mask);

if framesMsk > 1
    mask = mask(:,:,(1:lastFrameBeforeTail));
end


maskIn = mask;
maskOut = ~mask;


if outMaskFlag == 1
    mask2DT = maskOut;
    numMsk = 1;
elseif outMaskFlag == 0
    mask2DT =  maskIn;
    numMsk = 1;
else
    mask2DT = maskOut;
    numMsk = 2;
end

for msk = 1:numMsk
    if msk ~= 1
        mask2DT = ~mask2DT;
    end
    temp2DT = double(mask2DT).*double(movie);
    area2DT = squeeze(sum(sum(double(mask2DT))));  
    trace2DT = squeeze(double(sum(sum(temp2DT))))./area2DT;
    

%% detrend trace
    [trace_dt,~,~] = nt_detrend(trace2DT,50,[],'sinusoids',[],4,[]);
    %[trace_dt,wi,~] = nt_detrend(trace2DT,50,[],'polynomials',[],4,[]);
    
    template = trace2DT - trace_dt;
    
    if msk == 1
        [~, locmin] = min(template);
        templateFactor = template./(template(locmin));
    else
        templateFactorB = template./(template(locmin));
        templateFactor = (templateFactor*(1-w) + templateFactorB*(w));
    end
    
    
end


%% multiply by movie by factor (this is kind of stupid and wordy way but whatever)
 movie_dt = uint16(zeros (size(movie)));
 for i = 1: length(templateFactor)
     movie_dt(:,:,i) = (movie(:,:,i))./templateFactor(i);
 end
 
%xtract roi from dt_movie
% temp = double(Data.synapsis(indexData(i)).Mask).*double(roi.all(i).movPerSyn(:,:,:));
% area = sum(sum(double(Data.synapsis(indexData(i)).Mask)));    
% roi.all(i).roi_trace = squeeze(double(sum(sum(temp))/area));

if ~isempty(varargin) && mask ~= 1
    maskIn = mask2;
    maskOut = ~mask2;
end

temp = double(maskIn).*double(movie_dt);
area = squeeze(sum(sum(double(maskIn))));  
trace_dtIn = squeeze(double(sum(sum(temp))))./area;

%% jan 2021
if eqPreFrFlag == 1
   eqtrace = eqpreframes(trace_dtIn,2,200);
   
    templateFactor = eqtrace ./ trace_dtIn;
   
    for i = 1: length(templateFactor)
        movie_dt(:,:,i) = double(movie_dt(:,:,i)).*templateFactor(i);
    end
    
    temp = double(maskIn).*double(movie_dt);
    area = squeeze(sum(sum(double(maskIn))));  
    trace_dtIn = squeeze(double(sum(sum(temp))))./area;
 
end

temp = double(maskOut).*double(movie_dt);
area = squeeze(sum(sum(double(maskOut))));  
trace_dtOut = squeeze(double(sum(sum(temp))))./area;

 


