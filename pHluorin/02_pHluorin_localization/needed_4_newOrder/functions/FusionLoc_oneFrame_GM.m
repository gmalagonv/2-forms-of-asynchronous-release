
%[pstruct, mask, imgLM, imgLoG] = pointSourceDetection(img, sigma, mode)
%
% Inputs :   
%                 img : input image
%               sigma : standard deviation of the Gaussian PSF
%
% Options (as 'specifier'-value pairs): 
%
%              'mode' : parameters to estimate. Default: 'xyAc'.
%             'alpha' : alpha value used in the statistical tests. Default: 0.05.
%              'mask' : mask of pixels (i.e., cell mask) to include in the detection. Default: all.
%       'FitMixtures' : true|{false}. Toggles mixture-model fitting.
%       'MaxMixtures' : maximum number of mixtures to fit. Default: 5.
%   'RemoveRedundant' : {true}|false. Discard localizations that coincide within 'RedundancyRadius'.
%  'RedundancyRadius' : Radius for filtering out redundant localizatios. Default: 0.25
%         'Prefilter' : {true}|false. Prefilter to calculate mask of significant pixels.
%     'RefineMaskLoG' : {true}|false. Apply threshold to LoG-filtered img to refine mask of significant pixels.
%   'RefineMaskValid' : {true}|false. Return only mask regions where a significant signal was localized.
%        'ConfRadius' : Confidence radius for positions, beyond which the fit is rejected. Default: 2*sigma
%        'WindowSize' : Window size for the fit. Default: 4*sigma, i.e., [-4*sigma ... 4*sigma]^2
%
% Outputs:  
%             pstruct : output structure with Gaussian parameters, standard deviations, p-values
%                mask : mask of significant (in amplitude) pixels
%               imgLM : image of local maxima
%              imgLoG : Laplacian of Gaussian-filtered image

% Francois Aguet, April 2011 (last modified: 09/30/2013)

function [pstruct, mask, imgLM, imgLoG, numLm] = FusionLoc_oneFrame_GM(img, sigma, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
%ip.addRequired('frame', @isnumeric);
ip.addRequired('img', @isnumeric);
ip.addRequired('sigma', @isscalar);
%%
%ip.addParamValue('DeltaAcTresh', 2, @isnumeric);%commented out,GM July 18, 2020
ip.addParamValue('OnlyOne', true ,@islogical);
ip.addParamValue('BestAfterFit', false, @islogical); %% GM, july 2020
ip.addParamValue('simpleBackground', false, @islogical); %% GM, Nov 2020
ip.addParamValue('imposedBackground', 0, @isnumeric);%% GM, Nov 2020
ip.addParamValue('imposedA', false, @islogical);%% GM, Nov 2020

%%

ip.addParamValue('Mode', 'xyAc', @ischar);
ip.addParamValue('Alpha', 0.05, @isscalar);
ip.addParamValue('Mask', [], @(x) isnumeric(x) || islogical(x));
%%
ip.addParamValue('FitMixtures', false, @islogical); %% GM, somehow here, loc of events by this way have bigger STD, and more multiple
%%
ip.addParamValue('MaxMixtures', 5, @isposint);
ip.addParamValue('RemoveRedundant', true, @islogical);
ip.addParamValue('RedundancyRadius', 0.25, @isscalar);
ip.addParamValue('Prefilter', true, @islogical);
ip.addParamValue('RefineMaskLoG', true, @islogical);
ip.addParamValue('RefineMaskValid', true, @islogical);
ip.addParamValue('ConfRadius', []); % Default: 2*sigma, see fitGaussians2D.
ip.addParamValue('WindowSize', []); % Default: 4*sigma, see fitGaussians2D.
ip.KeepUnmatched = true;
ip.parse(img, sigma, varargin{:});
mode = ip.Results.Mode;
alpha = ip.Results.Alpha;









% postFrame = double(postFrames(:,:,379));
% preFrame = double(preFrames(:,:,379));
% sigma = 1.5;
% Mask = double(Masks(:,:,379));
% alpha = 0.05;


%% commented out, GM July 18 2020, 
%DeltaAcTresh = ip.Results.DeltaAcTresh;
%TreshFactor = 1;



    
%     if i == 1
%         img = preFrame;
%     else
%         img = postFrame;
%     end
numLm = 0;
%img = frame;
i = 1; %% maybe i is not needed, but I'm lazy. GM july 2020
if ~isa(img, 'double')
    img = double(img);
end

% Gaussian kernel
w = ceil(4*sigma); %w = ceil((4*0.5)*sigma);
x = -w:w;
g = exp(-x.^2/(2*sigma^2));
u = ones(1,length(x));

% convolutions
imgXT = padarrayXT(img, [w w], 'symmetric');
fg = conv2(g', g, imgXT, 'valid');
fu = conv2(u', u, imgXT, 'valid');
fu2 = conv2(u', u, imgXT.^2, 'valid');

% Laplacian of Gaussian
gx2 = g.*x.^2;
imgLoG = 2*fg/sigma^2 - (conv2(g, gx2, imgXT, 'valid')+conv2(gx2, g, imgXT, 'valid'))/sigma^4;
imgLoG = imgLoG / (2*pi*sigma^2);

% 2-D kernel
g = g'*g;
n = numel(g);
gsum = sum(g(:));
g2sum = sum(g(:).^2);

% solution to linear system
A_est(:,:,i) = (fg - gsum*fu/n) / (g2sum - gsum^2/n);
c_est(:,:,i) = (fu - A_est(:,:,i)*gsum)/n;


%% mask relevant pixels over Noise.

if ip.Results.Prefilter

    J = [g(:) ones(n,1)]; % g_dA g_dc
    C = inv(J'*J);

    f_c = fu2 - 2*c_est(:,:,i).*fu + n*c_est(:,:,i).^2; % f-c
    RSS = A_est(:,:,i).^2*g2sum - 2*A_est(:,:,i).*(fg - c_est(:,:,i)*gsum) + f_c;
    RSS(RSS<0) = 0; % negative numbers may result from machine epsilon/roundoff precision
    sigma_e2 = RSS/(n-3);

    sigma_A = sqrt(sigma_e2*C(1,1));

    % standard deviation of residuals
    sigma_res = sqrt(RSS/(n-1));

    kLevel = norminv(1-alpha/2.0, 0, 1);

    SE_sigma_c = sigma_res/sqrt(2*(n-1)) * kLevel;
    df2 = (n-1) * (sigma_A.^2 + SE_sigma_c.^2).^2 ./ (sigma_A.^4 + SE_sigma_c.^4);
    scomb = sqrt((sigma_A.^2 + SE_sigma_c.^2)/n);
    T = (A_est(:,:,i) - sigma_res*kLevel) ./ scomb;
    pval = tcdf(-T, df2);

    % mask of admissible positions for local maxima
    mask = pval < 0.05;
else
    mask = true(size(img));
end
% all local max
allMax = locmax2d(imgLoG, 2*ceil(sigma)+1);

% local maxima above threshold in image domain
imgLM = allMax .* mask;


if sum(imgLM(:))~=0 % no local maxima found, likely a background image

    if ip.Results.RefineMaskLoG
        % -> set threshold in LoG domain
        logThreshold = min(imgLoG(imgLM~=0));
        logMask = imgLoG >= logThreshold;

        % combine masks
        mask = mask | logMask;
    end

    % re-select local maxima
    imgLM = allMax .* mask;

    % apply exclusion mask
    if ~isempty(ip.Results.Mask)
        imgLM(ip.Results.Mask==0) = 0;
    end
end
%% *** commented out, GM July 21, 2020 
% if sum(imgLM(:))~=0
% %             TreshFactor = 1+(A_est(imgLM~=0)-1000)/3000; %this generate a factor that increase proportional if preframe have something to localize.
%     TreshFactor = 1+(max(A_est(imgLM~=0))-1000)/4000; %this generate a factor that increase proportional if preframe have something to localize.
% 
% end
    
%%   ********************************************


pstruct = [];
if sum(imgLM(:))~=0 % no local maxima found, likely a background image
    if ip.Results.RefineMaskLoG
        % -> set threshold in LoG domain
        logThreshold = min(imgLoG(imgLM~=0));
        logMask = imgLoG >= logThreshold;
        
        % combine masks
        mask = mask | logMask;
    end
    
    % re-select local maxima
    imgLM = allMax .* mask;
    
    % apply exclusion mask
    if ~isempty(ip.Results.Mask)
        imgLM(ip.Results.Mask == 0) = 0;
    end
      
    
    % Treshold x DeltaInt between pre and post Frame.
    %% July 18, GM. comented out, only one frame checked *************************************
    
%     Ac_est = A_est + c_est;
%     
%     h = fspecial('average', round(2*sigma));
%     e = imdilate(Ac_est, ceil(h))-imfilter(Ac_est, h);
%     DeltaAc_est_Treshold = median(e(:))*DeltaAcTresh*TreshFactor;
%     
%     DeltaAc_est = Ac_est(:,:,2)-imdilate(Ac_est(:,:,1), strel('disk', ceil(sigma)));
%     DeltaAc_est_Mask = DeltaAc_est; 
%     DeltaAc_est_Mask(DeltaAc_est_Mask<DeltaAc_est_Treshold) = 0; 
%     DeltaAc_est_Mask(DeltaAc_est_Mask>=DeltaAc_est_Treshold) = 1;
% 
%     imgLM = imgLM .* DeltaAc_est_Mask;
%     
%     %     
%     A_est = A_est(:,:,2);
%     c_est = c_est(:,:,2);
%% ************************************************************************

 [lmy, lmx] = find(imgLM~=0);     
 lmIdx = sub2ind(size(img), lmy, lmx);  

%% July 18, GM. comented out, only one frame checked     
%     if ip.Results.OnlyOne  % para que solo se quede con el mejor maximo?
%         [~,I] = max(DeltaAc_est(lmIdx));
%         lmIdx =  lmIdx(I);
%         [lmy,lmx] = ind2sub(size(img),lmIdx);
%     end
%% **************************************************************** 
%% July 18, GM

% if ip.Results.OnlyOne  % select only one max.... here the diference with FusionLoc
%     Ac_est = A_est + c_est;
%      [~,I] = max(Ac_est(lmIdx));
%      lmIdx =  lmIdx(I);
%      [lmy,lmx] = ind2sub(size(img),lmIdx);          
% end

if ip.Results.OnlyOne  % select only one max.... here the diference with FusionLoc
    %Ac_est = A_est + c_est;
     [~,I] = max(imgLM(lmIdx));
     lmIdx =  lmIdx(I);
     [lmy,lmx] = ind2sub(size(img),lmIdx);          
end


%%
    
    if ~isempty(lmIdx)
        numLm = numLm + length(lmx); %% GM 
        % run localization on local maxima
        
        %% GM, nov 2020. overwrite the background value simply by the AVG of the contour outside the mask
        if ip.Results.simpleBackground
            c_est(lmIdx) = mean(double(img) .* double(~(mask)),'all');
        end
        %% GM, nov 2020. overwrite the background value, plug a value coming from the avg of the expected unresponsive frames in the detrended movie ()
        %% it would be used ip.Results.simpleBackground or  ip.Results.imposedBackground
%         if length(ip.Results.imposedBackground) > 1
%             c_est = ip.Results.imposedBackground;
%         end
        
           %% GM nov 2020
   
        
        if sum(sum(ip.Results.imposedBackground)) > 1
            c_estOld = c_est;
            %c_est = ip.Results.imposedBackground;
    %        
            %c_est = (ip.Results.imposedBackground + c_est)./2.07;
            c_calc = ip.Results.imposedBackground;
            %c_est = (median(cat(3, c_calc, c_est),3)).*0.95;
            c_est = c_calc;
            %c_est = mean(ip.Results.imposedBackground, c_est,3);


    % %        
           if ~isempty(lmIdx)
                %% smooth the calculated background
               %c_est = medfilt2(c_est); % 
               
                %% factor1 1 (over A)
                % original:
                
%                 %% this was used succesfully 8888888888888
%                   fact = (c_est - c_estOld);
%                   A_est = A_est - fact;
% %                  %**********************************
% %                  
% %                     %% testB:  (used in 141, xy, local)
%                  res = img(lmIdx)-(c_est(lmIdx) + A_est(lmIdx));
%                  if  res > 0
%                      %A_est = A_est + res;
%                      c_est = c_est + res/2;
%                  else
%                      A_est = A_est + res;
%                      %c_est = c_est + res;
% 
%                  end
%                  
%% 88888888888888888888888888888888888888888888888888888888888888888888888888888888888
                  %% factor2 (over c)
                  %% original:
                  %factor2 = (img(lmIdx)-c_est(lmIdx))-A_est(lmIdx);
                  %c_est = c_est + factor2;
                  %***************************
                  
                   %% factor2 (140 in excel)
                  %% testA: decent, an option 
%                     idealA = img(lmIdx)-c_est(lmIdx);
%                     factor2 = idealA - A_est(lmIdx);
%                     A_est = A_est + factor2*(1/4);
%                     c_est = c_est - factor2*(3/4); %-
% %***********************


                   %% factor2 (140 in excel)
                  %% testA2: decent, an option 
%                     idealA = img(lmIdx)-c_est(lmIdx);
%                     factor2 = idealA - A_est(lmIdx);
%                     if factor2 > 0
%                         A_est = A_est + factor2*(3/4);
%                         c_est = c_est - factor2*(1/4); %-
%                     else
%                         A_est = A_est - factor2*(3/4);
%                         c_est = c_est + factor2*(1/4); %-
%                     end
                    %% testA3: decent, an option 
%                     idealA = img(lmIdx)-c_est(lmIdx);
%                     factor2 = idealA - A_est(lmIdx);
%                     if factor2 > 0
%                         A_est = A_est + factor2*(4/5);
%                         c_est = c_est - factor2*(1/5); %-
%                     else
%                         A_est = A_est - factor2*(4/5);
%                         c_est = c_est + factor2*(1/5); %-
%                     end
%                 
% %***********************


                      
% 
% %********* 
               
                    %% testC:  ( xy, local)
%                     res = img(lmIdx)-(c_est(lmIdx) + A_est(lmIdx));
%                      if res > 0
%                          %A_est = A_est + res;
%                          c_est = c_est - res*(3/4);
%                          A_est = A_est + res*(1/4);
%                      else
%                          A_est = A_est - res*(1/4);
%                          c_est = c_est + res*(3/4);
%                          
%                      end
%                      
             




                  %% test:
%                     idealA = img-c_est;
%                     
%                     factor2 = idealA - A_est;
%                     factor2 = medfilt2(factor2);
%                     A_est = A_est + factor2*(1/4);
%                     c_est = c_est - factor2*(3/4);




%                
                  %***************************
                
                 
           end
        end

        
        %% GM, nov 2020. overwrite the A value, plug simply the value of the img
        if ip.Results.imposedA
            A_est = (img-c_est);
        end
        
        
        if ~ip.Results.FitMixtures            
            pstruct = fitGaussians2D(img, lmx, lmy, A_est(lmIdx), sigma*ones(1,length(lmIdx)),...
                double(c_est(lmIdx)), mode, 'mask', mask, 'alpha', alpha,...
                'ConfRadius', ip.Results.ConfRadius, 'WindowSize', ip.Results.WindowSize);
        else
            pstruct = fitGaussianMixtures2D(img, lmx, lmy, A_est(lmIdx), sigma*ones(1,length(lmIdx)),...
                c_est(lmIdx), 'mask', mask, 'alpha', alpha, 'maxM', ip.Results.MaxMixtures);
        end
        
        % remove NaN values
        idx = ~isnan([pstruct.x]);
        if sum(idx)~=0
            fnames = fieldnames(pstruct);
            for k = 1:length(fnames)
                pstruct.(fnames{k}) = pstruct.(fnames{k})(idx);
            end
            
            % significant amplitudes
            idx = [pstruct.hval_Ar] == 1;
            
            % eliminate duplicate positions (resulting from localization)
            if ip.Results.RemoveRedundant
                pM = [pstruct.x' pstruct.y'];
                idxKD = KDTreeBallQuery(pM, pM, ip.Results.RedundancyRadius*ones(numel(pstruct.x),1));
                idxKD = idxKD(cellfun(@numel, idxKD)>1);
                for k = 1:length(idxKD)
                    RSS = pstruct.RSS(idxKD{k});
                    idx(idxKD{k}(RSS ~= min(RSS))) = 0;
                end
            end
            %% July 23,2020. GM
            if ip.Results.BestAfterFit && sum(idx) > 1
               %factor = [pstruct.x_pstd] + [pstruct.y_pstd] + [pstruct.A_pstd] + 1./[pstruct.A];
               factor = [pstruct.x_pstd] + [pstruct.y_pstd] + ([pstruct.A_pstd]./[pstruct.A]);
               idxidx  =  (1:1:length(idx));
               idxidx = idxidx(idx);
               [~, locMin] = min(factor(idx));
               locMin =  idxidx(locMin);
               idx(locMin) = 0;
                 
            end           
            %%
            
            
            if sum(idx)>0
                fnames = fieldnames(pstruct);
                for k = 1:length(fnames)
                    pstruct.(fnames{k}) = pstruct.(fnames{k})(idx);
                end
                pstruct.hval_Ar = logical(pstruct.hval_Ar);
                pstruct.hval_AD = logical(pstruct.hval_AD);
                pstruct.isPSF = ~pstruct.hval_AD;
                
                % adjust mixture index if only a single component remains
                if ip.Results.FitMixtures
                    mv = 0:max(pstruct.mixtureIndex);
                    multiplicity = getMultiplicity(pstruct.mixtureIndex);
                    pstruct.mixtureIndex(ismember(pstruct.mixtureIndex, mv(multiplicity==1))) = 0;
                end
            else
                pstruct = [];
            end
        else
            pstruct = [];
        end
  
    end
end

if ~isempty(pstruct) && ip.Results.RefineMaskValid
    CC = bwconncomp(mask);
    labels = labelmatrix(CC);
    loclabels = labels(sub2ind(size(img), pstruct.y_init, pstruct.x_init));
    idx = setdiff(1:CC.NumObjects, loclabels);
    CC.PixelIdxList(idx) = [];
    CC.NumObjects = length(CC.PixelIdxList);
    mask = labelmatrix(CC)~=0;
end  
  

%% commented out by GM july 2020. Not needed for this purporses
 if ~isempty(pstruct)
     pstruct.initA = A_est(lmIdx);
     pstruct.initc = c_est(lmIdx);
     pstruct.valA = img(lmIdx);
     pstruct.A_est = A_est;
     pstruct.c_est = c_est;
 end
%     
%     pstruct.DeltaInt_est = DeltaAc_est(lmIdx);
%     
%     AvegFilter = fspecial('gaussian',7);
% 
%     pstruct.LocFrame.Frame = postFrame;
%     pstruct.PreFrame.Frame = preFrame;
% 
%     Maximo = double(max([max(pstruct.LocFrame.Frame(:)) max(pstruct.PreFrame.Frame(:))]));
%     Minimo = double(min([min(pstruct.LocFrame.Frame(:)) min(pstruct.PreFrame.Frame(:))]));
% 
%     pstruct.LocFrame.FrameNormalizado = (double(pstruct.LocFrame.Frame)-Minimo)/(Maximo-Minimo);
%     pstruct.PreFrame.FrameNormalizado = (double(pstruct.PreFrame.Frame)-Minimo)/(Maximo-Minimo);
% 
%     pstruct.LocFrame.FrameNormalizado2 = (double(pstruct.LocFrame.Frame)-pstruct.c)/(Maximo-pstruct.c);
%     pstruct.LocFrame.FrameNormalizado2(pstruct.LocFrame.FrameNormalizado2<0) = 0;
%     pstruct.PreFrame.FrameNormalizado2 = (double(pstruct.PreFrame.Frame)-pstruct.c)/(Maximo-pstruct.c);
%     pstruct.PreFrame.FrameNormalizado2(pstruct.PreFrame.FrameNormalizado2<0) = 0;
%     pstruct.LocFrame.FrameFilter = imfilter(pstruct.LocFrame.FrameNormalizado2, AvegFilter, 'replicate');
%     pstruct.LocFrame.FrameFilter = pstruct.LocFrame.FrameFilter/max(pstruct.LocFrame.FrameFilter(:));
%     pstruct.PreFrame.FrameFilter = imfilter(pstruct.PreFrame.FrameNormalizado2, AvegFilter, 'replicate');
%     pstruct.PreFrame.FrameFilter = pstruct.PreFrame.FrameFilter/max(pstruct.LocFrame.FrameFilter(:));
%     
%     [ FF, CC ] = size(pstruct.PreFrame.FrameNormalizado);
%     pstruct.EventPicture = zeros(CC*2+4,FF*2+4);
%     pstruct.EventPicture(2:FF+1,2:CC+1) = pstruct.PreFrame.FrameNormalizado;
%     pstruct.EventPicture(2:FF+1,CC+4:end-1) = pstruct.LocFrame.FrameNormalizado;
%     pstruct.EventPicture(FF+4:end-1,2:CC+1) = pstruct.PreFrame.FrameFilter;
%     pstruct.EventPicture(FF+4:end-1,CC+4:end-1) = pstruct.LocFrame.FrameFilter;
%     
%     pstruct.EventPictureLoc(:,:,3) = pstruct.EventPicture;
%     pstruct.EventPictureLoc(:,:,2) = pstruct.EventPicture;
%     pstruct.EventPictureLoc(:,:,1) = pstruct.EventPicture;
%     XX = floor(pstruct.x);YY = floor(pstruct.y);
%     pstruct.EventPictureLoc(1+YY,CC+3+XX,:) = [ 1,0,0 ];
%     pstruct.EventPictureLoc(FF+3+YY,CC+3+XX,:) = [ 1,0,0 ];
% 
 end
    








