function averages = AVGtracesFromALL(ALL, criteria, addtionalField,  binSize, numExp,normFlag, andorFlag,varargin)
%% numExp = number of exponnetials for fit
%% varargin{1} number of trace to avg (number is the order it shows up in the structure fomr left to right)
%% synchFlag = to know if we are analysing 1st frame Data. 
%% 06.13.2020. highPassFlag + if traces are going to be filtered using a high pass filter. Cutoff freq = 0.1, conservative. 
%% 09.17.2020. varargin added as input, works as a : Plot Flag, being varargin{1} the number of the figure used 
%% or  if varargin{1}= 1, use the next figure availabe. varargin{2} = number of bin to plot. 
%% 11.05.2020. add the possibility to xtract traces from NEWORDER function
%**************************************************************************
%% FLAGS
highPassFlag = 0;
synchFlag = 0; %  for old program


flagRandReplacement = 0; 
%proximalLIM = 50;
%distalLIM = 250;
minNumEv = 5;
percentageFlag = 0; % 1 or 2, 1= Analysis per syn, but only the ones with a min number of events. 2 = all events in the % claculated per syn. 04.08.2020, GM. 
%percentageNolimFlag = 1;
percentage = 2;
preStim = 20;%20; %% nomally 20

%%


%%

if nargin >= 7 %&& EvClass == 1
    allFields = fieldnames(ALL.allEv);
    trace2use = varargin{2};
else
    trace2use = 1;

end
%% tau
if length(numExp) > 1
    numExp = 3;
    flagLocalTau = true;
%     meanTrace = forAVGtracesFromALL(ALL,preStim, highPassFlag,normFlag, andorFlag,trace2use);
%     [~,params,~]= expFitDecayPLusLineTrace(meanTrace, 20, 19, 5, 0, 1);    
%     Tau = params(2);


else
    %numExp = 3;
    flagLocalTau = false;
    Tau = 0.0439;
end



%% default criteria : 'dist2centFSynch'
if ~ischar(criteria) && criteria == 1
    criteria = 'dist2centFSynch';
end
    


if isfield(ALL.allEv, 'synch') 
    separated = struct;
    masks = cell(3,2);
    %% SYNCH
    masks{1,1} = logical([ALL.allEv.synch]);
    masks{1,2} = 'Synch';
    %% ASYNCH 2nd
    masks{2,1} = logical([ALL.allEv.asynch2nd]);
    masks{2,2} = 'Asynch2nd';
    %% ASYNCH ALLOC
    masks{3,1} = logical([ALL.allEv.asynchAlloc]);
    masks{3,2} = 'AsynchAll';
%     %% SYNCH + 2nd ASYNCH
%     masks{2,1} = logical(logical([ALL.allEv.synch]) + logical([ALL.allEv.asynch2nd]));
%     masks{2,2} = 'SynchPlus2nd';

%     %% ALL
%     masks{5,1} = true(1, length(ALL.allEv));
%     masks{5,2} = 'All';
    %%
else
    
    masks = cell(1,1);
    %% ALL
    masks{1,1} = true(1, length(ALL.allEv));
    masks{1,2} = 'All';
    %% Sept 2020, to add a df
    if ~isfield(ALL.allEv, "df")
        ALL = addDfNorm2ALL(ALL, synchFlag);
    end
    
    
end
[numEvClasses, ~ ]= size(masks);
ALLorig = ALL;

if nargin >= 6
    numEvClasses = varargin{1};
end
    



for EvClass = 1:numEvClasses
    averages = struct;
    mskClass = masks{EvClass,1};
    if sum(mskClass) >= 1 
        ALL.allEv = ALLorig.allEv(mskClass);

        %% here traces are binned by dist2centF, so binSize is a distance in nm

        all_traces = struct2cell(ALL.allEv);

        long_length = cell2mat(cellfun(@size,(all_traces(:,:,1)),'UniformOutput',false));
        long_loc = find((long_length(:,1) > 1) == true);

        %% 03.26.2020, GM
        if nargin >= 7 %&& EvClass == 1
            %allFields = fieldnames(ALL.allEv);
            disp(['TRACES USED: ' (allFields{long_loc(trace2use)})])
            %trace2use = varargin{2};
           
%         
%         else
%             trace2use = 1;
            
        end
        
        all_traces_subs = all_traces(long_loc(trace2use),:,:); %%
        nanMask = (~any( isnan(squeeze(cell2mat( all_traces_subs)))))';

        %f1s = cellfun(@(x) (x(preStim + 1)-x(preStim)),all_traces_subs,'UniformOutput',false);
        f1s = cellfun(@(x) (x(preStim + 1)-mean(x(preStim-4:preStim))),all_traces_subs,'UniformOutput',false);

        f1s = squeeze(cell2mat(f1s));

        f1msk  = f1s > 0;
      %% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        nanMask = nanMask & f1msk;
      %% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        all_traces_subs = all_traces_subs(nanMask);
        %% test oct 2020, used to scale downs data taken with Andor camera
        if andorFlag >= 1 
            all_traces_subs =  cellfun(@(x) (x ./19.135), all_traces_subs,'UniformOutput',false);
            if andorFlag == 2 
                all_traces_subs =  cellfun(@(x) (x ./19.0370), all_traces_subs,'UniformOutput',false);
            end
        end
        
            
        
        
        
        
        %% 06.13.2020, GM. highPassFlag
        if highPassFlag == 1
            all_traces_subs = cellfun(@(x) highpass(x,0.1,20), all_traces_subs,'UniformOutput',false);
            %all_traces_subs = cellfun(@highpass, all_traces_subs, 0.1, mat2cell(20), 'UniformOutput',false); 
        end
     
      %  all_traces_subs = cellfun(@(x) ((x-mean(x(preStim-4:preStim)))), all_traces_subs,'UniformOutput',false); %% first point(before stim) set at 500, just to not to have negative values;

        %% 10.02.2020 GM, normFlag NORMALIZATION, commented out by DEC 2021
        if normFlag > 0
            for tr = 1:length(all_traces_subs)
                tra_ce = all_traces_subs{tr};
                %tra_ce = (tra_ce - min(tra_ce));
                %tra_ce = (tra_ce - mins{tr})./(maxs{tr}-mins{tr}) ;
                if normFlag == 1 
                    f0 = mean(tra_ce((preStim-4):preStim));
                    %std0 = std(tra_ce((preStim-4):preStim));
                    all_traces_subs{tr} = (tra_ce -f0 )./f0;
                elseif normFlag == 2 
%                   
                    all_traces_subs{tr} = (tra_ce - min(tra_ce))./(max(tra_ce)- min(tra_ce));
                elseif normFlag == 3
                    
                    %val = mean(tra_ce((preStim-4):preStim));
                    val = min(tra_ce);
                    tra_ce = tra_ce - val;
                    tra_ce = tra_ce ./(tra_ce(preStim +1));
%                     if tra_ce(preStim +1) ~= 1
%                         disp ('here')
%                     end
                    all_traces_subs{tr} = tra_ce;
                    %disp('here')
                elseif normFlag == 4
                    %val = mean(tra_ce((preStim-4):preStim));
                    val = tra_ce(preStim +1);%mean(tra_ce(preStim-4:preStim));
                    %tra_ce = tra_ce - val;
                    tra_ce = tra_ce ./val;
                    all_traces_subs{tr} = tra_ce;



                end

                %all_traces_subs{tr} = (all_traces_subs{tr}-f0{tr})/f0{tr};
            end

        end

%
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>
   %% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        all_traces_subs = cellfun(@(x) ((x-mean(x(preStim-4:preStim)))), all_traces_subs,'UniformOutput',false); %% first point(before stim) set at 500, just to not to have negative values;
        %% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       if normFlag == 5
           avgPeak = (mean(cell2mat(all_traces_subs),3));
           avgPeak = avgPeak(preStim + 1);
           all_traces_subs = cellfun(@(x) (x ./ avgPeak), all_traces_subs,'UniformOutput',false);
           
       end
%        
     
        
        
        all_traces_dns  =  cellfun(@(x) denoiseForwrdBack_g(x,20,10),all_traces_subs,'UniformOutput',false); 
        all_traces_dns  =  cellfun(@(x1,x2) ((x1+x2)/2), all_traces_subs,all_traces_dns, 'UniformOutput',false); 
        
        %all_traces_subs = cellfun(@(x) (x-(x(preStim))), all_traces_subs,'UniformOutput',false);
        %all_traces_subs = cellfun(@(x) ((x-x(preStim))), all_traces_subs,'UniformOutput',false); %% first point(before stim) set at 500, just to not to have negative values;





        %% DENOISED TRACE %% OBSOLETE BY NOV 2021
%         if length(long_loc) 
%             all_traces_dns = all_traces(long_loc(2),:,:);
%             all_traces_dns = cellfun(@(x) (x-x(preStim)) ,all_traces_dns,'UniformOutput',false); %% first point(before stim) set at 500, just to not to have negative values;
%             all_traces_dns = all_traces_dns(nanMask);
%         end

        %% 04.09.2020, GM.remove NaNs comming from traces that do not have enough pnts before the stim . 
%         if numEvClasses == 1
%                
%         else
%             all_dist2centF = [ALL.allEv.(criteria)];
%             %all_dist2centF = [ALL.allEv.dist2centFSynch];
%             %all_dist2centF = [ALL.allEv.distFromLast];
%             %all_dist2centF = [ALL.allEv.dist2centFAll];
%             %all_dist2centF = [ALL.allEv.timeFromLast];
%             %all_dist2centF = [ALL.allEv.distFromLast];
% 
%         end
        all_dist2centF = [ALL.allEv.(criteria)];
        all_dist2centF = all_dist2centF(nanMask);
        all_frameRelat = [ALL.allEv.frameRelat];
        all_frameRelat = all_frameRelat(nanMask);
        if ischar(addtionalField) 
            all_addtionalField = [ALL.allEv.(addtionalField)];
            all_addtionalField = all_addtionalField(nanMask);
        end

    %     if synchFlag == 1 %|| ~isfield(ALL.perSyn, 'forGraphFrom2')
    %         %% synch
    %         all_traces_slope = cellfun(@(x) (x(preStim + 2)-x(preStim + 1)),all_traces_subs,'UniformOutput',false);
    %         all_traces_ratio = cellfun(@(x) ((x(preStim + 1)-x(preStim + 3))/(x(preStim + 1))) ,all_traces_subs,'UniformOutput',false);
    %         %% 04.16.2020 slopes of the slow component
    %         xfit  = (0:1:17)';
    %         all_traces_slopeSlow = cellfun(@(x) polyfit(xfit, (x((preStim + 3):(preStim + 20))),1) ,all_traces_subs,'UniformOutput',false);
    %         all_traces_slopeSlow = cellfun(@(x) (x(1)),all_traces_slopeSlow,'UniformOutput',false);
    % 
    %     else
        %% SLOPES
        all_traces_slope = cellfun(@(x) (x(preStim + 2)-x(preStim + 1)),all_traces_subs,'UniformOutput',false);
         %% RATIOS
    %     all_traces_ratio = cellfun(@(x) (abs((x(preStim + 1))-(mean(x(preStim + 2: preStim +3)) - mean(x(preStim-4:preStim)))) ...
    %         /(x(preStim + 1)- mean(x(preStim-4:preStim)))) ,all_traces_subs,'UniformOutput',false);

    %     all_traces_ratio = cellfun(@(x) (abs((x(preStim + 1)- mean(x(preStim-4:preStim)))-(mean(x(preStim + 2: preStim +3)) - mean(x(preStim-4:preStim)))) ...
    %         /(x(preStim + 1)- mean(x(preStim-4:preStim)))) ,all_traces_subs,'UniformOutput',false);
    %% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    option A: double exp fit of the decay FIT   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       % all_traces_ratio = cellfun(@(x) (expFitDecayTrace(x, preStim)),all_traces_subs,'UniformOutput',false);

        all_traces_ratio = cell (length(all_traces_subs), 1);
        all_traces_fit = cell (length(all_traces_subs), 1);
        all_traces_onlyFast = cell (length(all_traces_subs), 1);
        all_traces_onlySlow = cell (length(all_traces_subs), 1);
        
        fitParams = cell (length(all_traces_subs), 1);
        
  %% single or double exponential fit
        if numExp<=2
  
  
            [~,params,~]= expFitDecayTrace(mean(cell2mat(squeeze(all_traces_subs)')')', 20, numExp);

            if numExp == 2
                taus = [params(2),params(4)];
             %factor = 4;
                minTau = min(taus); %/ 2;
                maxTau = max(taus);% * 4;
                disp('double exponential Fit')
            else
                minTau = params(2);
                maxTau = 0;
                disp('single exponential Fit')
            end

            parfor xx = 1:length(all_traces_subs)
                trc = all_traces_subs{xx};
                %trc = all_traces_dns{xx};
                [all_traces_ratio{xx}, fitParams{xx}, ~]= expFitDecayTrace(trc, preStim, numExp, minTau, maxTau);

            end
        elseif numExp == 3
            
            if flagLocalTau
                [~,params,~]= expFitDecayPLusLineTrace(mean(cell2mat(squeeze(all_traces_subs)')')', 20, (21 - mean(all_frameRelat)), 5,0,1);
                Tau = params(2);
            end
%       
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  
            
        
            disp('single exponential & Line Fit')
            %slope = params(4);

    %         
            parfor xx = 1:length(all_traces_subs)
                trc = all_traces_subs{xx};
                lTrace = 21 - all_frameRelat(xx);
                %trc = all_traces_dns{xx};
                [all_traces_ratio{xx}, fitParams{xx}, all_traces_fit{xx}, all_traces_onlyFast{xx},all_traces_onlySlow{xx}]= expFitDecayPLusLineTrace(trc, preStim, lTrace,5 ,0,1,Tau);%, slope);
                %expFitDecayPLusLineTrace(trace, prestim, lTrace, filterFlag,varargin)

            end
            %% GM, Dec 2021. repeat the worst fits, this time low-pass filtering the linear component
            disp('remove the bad fits');
            fitParamsM = cell2mat(fitParams);
            [~, locResidual] = size(fitParamsM);
            locResidual = locResidual -1;
            mskfitParamsM = fitParamsM(:,locResidual) > (mean(fitParamsM(:,locResidual)) + 2*std(fitParamsM(:,locResidual)));
            %fitParams(mskfitParamsM,:) = NaN;
            %% NAN
            indxFitParamsM = (1:length(mskfitParamsM));
            indxFitParamsM  = indxFitParamsM(mskfitParamsM);
            for kk = 1:length(indxFitParamsM)
                fitParams{indxFitParamsM(kk)} = nan(1, (locResidual+1));
                
                %% remove traces too...
                all_traces_subs{indxFitParamsM(kk)} = nan(size(all_traces_subs{1}));
                all_traces_fit{indxFitParamsM(kk)} = nan(size(all_traces_fit{1}));
                all_traces_onlyFast{indxFitParamsM(kk)} = nan(size(all_traces_onlyFast{1}));
                all_traces_onlySlow{indxFitParamsM(kk)} = nan(size(all_traces_onlySlow{1}));
                
            end
%             
            
            
        end



    %% frequency domain: 
        all_freqSpect = cell (length(all_traces_subs), 1);    
        all_pwSpect = cell (length(all_traces_subs), 1);

        parfor yy = 1:length(all_traces_subs)
            trc = all_traces_subs{yy};
            [all_freqSpect{yy}, all_pwSpect{yy}] = freqAnalys(trc, preStim, 18);
        end





    %% option B
    % 
    %     all_traces_ratio = cell (length(all_traces_subs), 1);
    %     %f1s = zeros (length(all_traces_subs), 1);
    %     for xx = 1:length(all_traces_subs)
    %         trc = all_traces_subs{xx};
    %         %f1 = median(trc((preStim-4):preStim));
    %         f1 = trc(preStim + 1) ;%- f0;
    %         fm = median(trc(preStim + 2: preStim + 3));% - f0;
    %         %all_traces_ratio{xx} = (f1 - fm)/ f1;
    %         %f1s(xx) = f1;
    %         if f1 > fm
    %             rat = (f1 - fm)/ f1;
    %             if rat < 1
    %                 all_traces_ratio{xx} = rat;
    %             else
    %                 all_traces_ratio{xx} = 1;
    %             end
    %             
    %         else
    %             all_traces_ratio{xx} = 0;
    %             %f1s(xx) = NaN;
    %         end
    %         
    %         
    %         
    %     end
    %     %%%
        %% 04.16.2020 slopes of the slow component
        xfit  = (0:1:17)';
       
       

   
         
         
%         all_traces_slopeSlow = cellfun(@(x) polyfit(xfit, (x((preStim + 3):(preStim + 20))),1) ,all_traces_dns,'UniformOutput',false);
%% new  
        all_traces_slopeSlow = cellfun(@(x) polyfit(xfit, (x((preStim + 3):(preStim + 20))),1) ,all_traces_subs,'UniformOutput',false);
        %all_traces_slopeSlow = cellfun(@(x) polyfit(xfit, (x((preStim + 3):(preStim + 20))),1) ,all_traces_dns,'UniformOutput',false);
%%        
        
        all_traces_slopeSlow = cellfun(@(x) (x(1)*20),all_traces_slopeSlow,'UniformOutput',false);
     


        %% average allEv traces inside a given bin
       if length(binSize) > 1
           if binSize(1) == 0
               [BINS,~] = discretize(all_dist2centF, binSize(2));
           else
               BINS  = eqFreqBins(all_dist2centF, binSize(2));
           end
           
       else
           Maxbin = (ceil(max(all_dist2centF)/binSize))*binSize;
           edges = 0:binSize:Maxbin;
           BINS = discretize(all_dist2centF,edges);
           
       end
        
        indexesRand = (1:1:length(all_dist2centF));
        
        %% test: rmeove bad fits. DEC 20201
        
%         all_traces_subs(mskfitParamsM) = NaN; 
%         all_traces_slope(mskfitParamsM) = NaN;
%         all_traces_ratio(mskfitParamsM) = NaN;
%         fitParams(mskfitParamsM) = NaN;
%         all_traces_slopeSlow(mskfitParamsM) = NaN;
%         all_freqSpect(mskfitParamsM) = NaN; 
%         all_pwSpect(mskfitParamsM) = NaN;
        
        
        

        for i = 1:max(BINS)
            mask = (BINS == i);
            %% apply BIN's mask
            if length(long_loc) > 1
                traces_dns = all_traces_dns(mask);
            end
            traces_subs = all_traces_subs(mask); 
            traces_slope = all_traces_slope(mask);
            traces_ratio = all_traces_ratio(mask);
            traces_fit = all_traces_fit(mask);
            traces_onlyFast = all_traces_onlyFast(mask);
            traces_onlySlow  = all_traces_onlySlow(mask);
            
            traces_fitParams = fitParams(mask);
            traces_slopeSlow = all_traces_slopeSlow(mask);
            traces_freqSpect = all_freqSpect(mask); 
            traces_pwSpect = all_pwSpect(mask);

            %% 
            averages.bin(i).bin = [min(all_dist2centF(mask)) max(all_dist2centF(mask))];%[edges(i),edges(i+1)];
            averages.bin(i).numTraces = sum(mask);
            %% dist2cent
            averages.bin(i).dist2centF = (all_dist2centF(mask))';
            averages.bin(i).meandist2centF = mean(all_dist2centF(mask));
            averages.bin(i).semdist2centF = sem_gm(all_dist2centF(mask));
            %% addtionalField
            if ischar(addtionalField)
                averages.bin(i).(addtionalField) = all_addtionalField(mask)';
                averages.bin(i).(['mean' addtionalField]) = mean(all_addtionalField(mask));
                averages.bin(i).(['sem' addtionalField]) = sem_gm(all_addtionalField(mask));
            end
            
            %% Traces
            % original
            averages.bin(i).Traces = squeeze(cell2mat(traces_subs));
            averages.bin(i).meanTrace = (mean(cell2mat(traces_subs),3, 'omitnan'));
            
            averages.bin(i).varTrace = (var(cell2mat(traces_subs),0,3, 'omitnan'));
            averages.bin(i).stdTrace = (std(cell2mat(traces_subs),0,3, 'omitnan'));
            averages.bin(i).semTrace = [averages.bin(i).stdTrace] ./ sqrt(length(traces_subs));
            averages.bin(i).cvTrace = [averages.bin(i).stdTrace] ./ [averages.bin(i).meanTrace]; 
            
            % denoised
            if length(long_loc) > 1
                averages.bin(i).Traces_dns = squeeze(cell2mat(traces_dns));
                averages.bin(i).meanTrace_dns = (mean(cell2mat(traces_dns),3));
            end
            

            %% Slopes
            %fast
            averages.bin(i).Slopes = -1*( squeeze(cell2mat(traces_slope)));
            averages.bin(i).meanSlope = -1*( mean(cell2mat(traces_slope)));    
            averages.bin(i).semSlope = sem_gm(cell2mat(traces_slope));    
            %slow
            averages.bin(i).SlopesSlw = -1*(squeeze(cell2mat(traces_slopeSlow)));
            averages.bin(i).meanSlopeSlw = -1*(mean(cell2mat(traces_slopeSlow)));    
            averages.bin(i).semSlopeSlw = sem_gm(cell2mat(traces_slopeSlow));    

            %% ratios    

            selecratios = squeeze(cell2mat(traces_ratio));
            %badRatios = any( isnan( selecratios ) | isinf( selecratios ));% | selecratios < 0 | selecratios > 3, 2 );
            %selecratios = selecratios( ~any( isnan( selecratios ) | isinf( selecratios ) | selecratios < -10 | selecratios > 10, 2 ),: );%% 03.26.2020, GM
            %selecratios(badRatios) = NaN;
            averages.bin(i).ratios = selecratios;
            averages.bin(i).meanRatio = mean(selecratios, 'omitnan');
            averages.bin(i).semRatio = sem_gm(selecratios);
            averages.bin(i).medianRatio = median(selecratios, 'omitnan');
            fit = squeeze(cell2mat(traces_fitParams));
            averages.bin(i).fit = fit;
            if ~isempty(fit)
                if numExp >= 2
                    

                    averages.bin(i).meanFstau = mean(fit(:,2), 'omitnan');
                    averages.bin(i).medianFstau = median(fit(:,2), 'omitnan');
                    averages.bin(i).semFstau = sem_gm(fit(:,2));
                    
                    averages.bin(i).meanpeakFast = mean(fit(:,1), 'omitnan');
                    averages.bin(i).medianpeakFast = median(fit(:,1), 'omitnan');
                    averages.bin(i).sempeakFast = sem_gm(fit(:,1));

                    averages.bin(i).meanSltau = mean(fit(:,4), 'omitnan');
                    averages.bin(i).medianSltau = median(fit(:,4), 'omitnan');
                    averages.bin(i).semSltau = sem_gm(fit(:,4));
                    
                    averages.bin(i).meanIntercept = mean(fit(:,5), 'omitnan');
                    averages.bin(i).medianIntercept = median(fit(:,5), 'omitnan');
                    averages.bin(i).semIntercept = sem_gm(fit(:,5));
                    
              
                    
                    
                    averages.bin(i).meanFitTrace = (mean(cell2mat(traces_fit'),2, 'omitnan'));
                    averages.bin(i).stdFitTrace = (std(cell2mat(traces_fit'),0,2, 'omitnan'));
                    averages.bin(i).semFitTrace = [averages.bin(i).stdFitTrace] ./ sqrt(length(traces_fit));
                    
                    averages.bin(i).meanFastTrace = (mean(cell2mat(traces_onlyFast'),2, 'omitnan'));
                    averages.bin(i).stdFastTrace = (std(cell2mat(traces_onlyFast'),0,2, 'omitnan'));
                    averages.bin(i).semFastTrace = [averages.bin(i).stdFastTrace] ./ sqrt(length(traces_onlyFast));
                    
                    averages.bin(i).meanSlowTrace = (mean(cell2mat(traces_onlySlow'),2, 'omitnan'));
                    averages.bin(i).stdSlowTrace = (std(cell2mat(traces_onlySlow'),0,2, 'omitnan'));
                    averages.bin(i).semSlowTrace = [averages.bin(i).stdSlowTrace] ./ sqrt(length(traces_onlySlow));
                    
                    
                    

                    
                else
                    averages.bin(i).meanTau = mean(fit(:,2), 'omitnan');
                    averages.bin(i).medianTau = median(fit(:,2), 'omitnan');
                    averages.bin(i).semTau = sem_gm(fit(:,2));                   
                    
                end
            else
                if numExp >= 2
                    averages.bin(i).meanFastContrib = [];
                    averages.bin(i).medianFastContrib = [];
                    averages.bin(i).semFastContrib = [];

                    averages.bin(i).meanFstau = [];
                    averages.bin(i).medianFstau = [];
                    averages.bin(i).semFstau = [];

                    averages.bin(i).meanSltau = [];
                    averages.bin(i).medianSltau = [];
                    averages.bin(i).semSltau = [];
                else
                    averages.bin(i).meanTau = [];
                    averages.bin(i).medianTau = [];
                    averages.bin(i).semTau = [];                   
                    
                    
                end
             
            end


            %% freq stuff
            averages.bin(i).freqSpect = squeeze(cell2mat(traces_freqSpect));
            averages.bin(i).meanfreqSpect = mean(squeeze(cell2mat(traces_freqSpect)), 'omitnan');
            averages.bin(i).pwSpect = squeeze(cell2mat(traces_pwSpect));
            averages.bin(i).meanpwSpect = mean(squeeze(cell2mat(traces_pwSpect)), 'omitnan');
            averages.bin(i).sempwSpect = sem_gm(squeeze(cell2mat(traces_pwSpect)),2);





            %% random selection of same number of traces***********************************************************************
            if flagRandReplacement == 1
                X = randperm (length(mask),sum(mask)); %with replacement
            else
                X1 = randperm (length(indexesRand),sum(mask));
                X = indexesRand(X1); 
            end 

            maskRnd = false(length(mask),1);
            maskRnd(X)= true;

            if length(long_loc) > 1
                traces_dnsRnd = all_traces_dns(maskRnd);
            end

            traces_subsRnd = all_traces_subs(maskRnd);
            averages.binRand(i).numTraces = sum(maskRnd);
            averages.binRand(i).trace = (mean(cell2mat(traces_subsRnd),3));
            if length(long_loc) > 1
                averages.binRand(i).trace_dns = (mean(cell2mat(traces_dnsRnd),3));
            end

            averages.binRand(i).indexes = X;
            averages.binRand(i).meandist2centF = mean(all_dist2centF(maskRnd));
            if flagRandReplacement ~= 1
                 indexesRand(X1) = [];
            end
        end
        %% PERsYN! not yet set for NEW ORDER
        %%  % of max d2centF, percentageFlag must be = 1, a min amount of traces needed to average.
        if percentageFlag > 0
            averages.perSyn.traceProx = [];
            averages.perSyn.traceDist = [];
            %% for perSyn flag =2; 
            tracesProx =  {};
            tracesDist = {};
            dist2centProx = [];
            dist2centDist = [];    
            slopeProx = {};
            slopeDist = {};
            ratioProx = {};
            ratioDist = {};

            %%
            %dist2centF_all = [ALL.bin.dist2centF];
            ID_all = [ALL.allEv.synNum].*[ALL.allEv.area];
            ID_all = ID_all(nanMask); % 04.09.2020, GM, remove events whose trace is NaNs
            ID_syn = [ALL.perSyn.synNum].*[ALL.perSyn.area];
            synNum_syn = [ALL.perSyn.synNum];

            for j = 1:length(synNum_syn)
                maskSyn = (ID_all == ID_syn(j));
                tracesSyn = all_traces_subs(maskSyn);
                tracesSyn_dns = all_traces_dns(maskSyn);
                dist2centSyn = all_dist2centF(maskSyn);
                traces_slopeSyn = all_traces_slope(maskSyn);
                traces_ratioSyn = all_traces_ratio(maskSyn);


                proximalLIM = (percentage*(max(dist2centSyn)-min(dist2centSyn))/100) + min(dist2centSyn);       
                distalLIM = ((100 -(percentage + 0))*(max(dist2centSyn)-min(dist2centSyn))/100) + min(dist2centSyn); % percentage +10 .. just to include more events         
                maskProximal = (dist2centSyn <= proximalLIM);
                maskDistal = (dist2centSyn >= distalLIM);
              %% alternative way to get distal and prox mask: same thing... 04.09.2020 GM
        %         dist2centSynPerCent = 100*((dist2centSyn - min(dist2centSyn))./(max(dist2centSyn)-min(dist2centSyn)));        
        %         maskProximal = (dist2centSynPerCent <= percentage);
        %         maskDistal = (dist2centSynPerCent >= (100 -(percentage+0)));
        %         
                if percentageFlag == 1
                    if sum(maskProximal) >= minNumEv && sum(maskDistal) >= minNumEv

                        if isempty([averages.perSyn.traceDist])
                            loc = 1;
                        else
                            loc = length(averages.perSyn) + 1;
                        end

                        averages.perSyn(loc).synNum = synNum_syn(j);
                        averages.perSyn(loc).synIndex = j;

                        averages.perSyn(loc).proximalLIM = proximalLIM;
                        averages.perSyn(loc).numTracesProx = sum(maskProximal);
                        averages.perSyn(loc).meandist2centFProx = mean(dist2centSyn(maskProximal));
                        averages.perSyn(loc).traceProx = (mean(cell2mat(tracesSyn(maskProximal)),3));
                        averages.perSyn(loc).traceProx_dns = (mean(cell2mat(tracesSyn_dns(maskProximal)),3));

                        averages.perSyn(loc).distalLIM = distalLIM;
                        averages.perSyn(loc).numTracesDist = sum(maskDistal);
                        averages.perSyn(loc).meandist2centFDist = mean(dist2centSyn(maskDistal));
                        averages.perSyn(loc).traceDist = (mean(cell2mat(tracesSyn(maskDistal)),3));
                        averages.perSyn(loc).traceDist_dns = (mean(cell2mat(tracesSyn_dns(maskDistal)),3));

                        %averages.perSyn(loc).dist2centMIN = min(dist2centSyn);
                        %averages.perSyn(loc).dist2centMAX = max(dist2centSyn);        
                        %averages.perSyn(loc).dist2centMEAN = mean(dist2centSyn);
                        %averages.perSyn(loc).dist2centCV = std(dist2centSyn)/mean(dist2centSyn);


                        %% plot

                        %figure(loc), plot(averages.perSyn(loc).proximal), hold on, plot(averages.perSyn(loc).distal), ...
                        %plot(averages.perSyn(loc).proximal_dns), hold on, plot(averages.perSyn(loc).distal_dns), hold off;
                    end
                else %04.08.2020, GM. 
                    tracesProx= cat(3,tracesProx, tracesSyn(maskProximal));
                    tracesDist = cat(3,tracesDist, tracesSyn(maskDistal));

                    dist2centProx = [dist2centProx (dist2centSyn(maskProximal))];
                    dist2centDist = [dist2centDist (dist2centSyn(maskDistal))];

                    slopeProx = cat(3, slopeProx, traces_slopeSyn(maskProximal));
                    slopeDist = cat(3, slopeDist, traces_slopeSyn(maskDistal));

                    ratioProx = cat(3, ratioProx,traces_ratioSyn(maskProximal));
                    ratioDist = cat(3, ratioDist,traces_ratioSyn(maskDistal));




                    %proximal = tracesSyn(maskProximal);%proximal = [proximal; squeeze(tracesSyn(maskProximal))];%04.08.2020, GM.
                    %proximal = [proximal; tracesSyn(maskProximal)];
                end

            end
            if percentageFlag > 1
                averages.perSyn(1).traceProx = mean(cell2mat(tracesProx),3);%proximal;%(mean(cell2mat(proximal),3));
                averages.perSyn(1).traceDist = mean(cell2mat(tracesDist),3);

                averages.perSyn(1).dist2centFProx = dist2centProx;
                averages.perSyn(1).meandist2centFProx = mean(dist2centProx);
                averages.perSyn(1).dist2centFDist = dist2centDist;
                averages.perSyn(1).meandist2centFDist = mean(dist2centDist);

                averages.perSyn(1).SlopesProx = squeeze(cell2mat(slopeProx));    
                averages.perSyn(1).meanSlopeProx = mean(cell2mat(slopeProx));
                averages.perSyn(1).semSlopeProx = sem_gm(cell2mat(slopeProx));        

                averages.perSyn(1).SlopesDist = squeeze(cell2mat(slopeDist));
                averages.perSyn(1).meanSlopeDist = mean(cell2mat(slopeDist));
                averages.perSyn(1).semSlopeDist = sem_gm(cell2mat(slopeDist));

                selecRatiosProx = squeeze(cell2mat(ratioProx));
                selecRatiosProx = selecRatiosProx( ~any( isnan( selecRatiosProx ) | isinf( selecRatiosProx ) | selecRatiosProx < -10 | selecRatiosProx > 10, 2 ),: );%% 03.26.2020, GM
                averages.perSyn(1).RatiosProx = selecRatiosProx;
                averages.perSyn(1).meanRatioProx = mean(selecRatiosProx);
                averages.perSyn(1).semRatioProx = sem_gm(selecRatiosProx);

                selecRatiosDist = squeeze(cell2mat(ratioDist));
                selecRatiosDist = selecRatiosDist( ~any( isnan( selecRatiosDist ) | isinf( selecRatiosDist ) | selecRatiosDist < -10 | selecRatiosDist > 10, 2 ),: );%% 03.26.2020, GM
                averages.perSyn(1).RatiosDist = selecRatiosDist;
                averages.perSyn(1).meanRatioDist = mean(selecRatiosDist);
                averages.perSyn(1).semRatioDist = sem_gm(selecRatiosDist);       



            end
        end


        %% oct 2020 baseline removal of the avg trace
        for bn = 1:length(averages.bin)
            if ~isempty(averages.bin(bn).meanTrace)
                %averages.bin(bn).meanTraceS = averages.bin(bn).meanTrace(20:39)-(averages.bin(bn).meanTrace(40:end)- averages.bin(bn).meanTrace(40));
                averages.bin(bn).meanTraceS = averages.bin(bn).meanTrace(preStim:(preStim+19))-(averages.bin(bn).meanTrace((preStim+20):(preStim+39))- averages.bin(bn).meanTrace(preStim+20));
            end
        end
        %%
%         if nargin > 8
%             if varargin{3} == 1
%                 figNum =  length(findobj('type','figure')) + 1;
%             else
%                 figNum = varargin{3};
%             end
% 
%             for bn = 1 : varargin{4}
%                 figure(figNum), plot(averages.bin(bn).meanTrace), hold on
%                 figure(figNum + 1), plot(averages.bin(bn).meanTraceS), hold on
%             end
% 
%             if length(varargin) > 4
%                 figure(figNum), title(varargin{5});
%                 figure(figNum + 1), title([varargin{5} '-sub']);
%             else
%                 figure(figNum), title([inputname(1) masks{EvClass,2}]);
%                 figure(figNum + 1), title([inputname(1) masks{EvClass,2} '-sub']);
%             end
%             figure(figNum), hold off;
%             figure(figNum + 1), hold off;
%         end

        if numEvClasses > 1
            separated.(masks{EvClass,2}) = averages;
        end

    
    if numEvClasses > 1
        averages = separated;
    end
    end
end

end
