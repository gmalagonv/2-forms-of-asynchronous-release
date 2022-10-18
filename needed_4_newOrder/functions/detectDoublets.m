function Data = detectDoublets(Data,typeTrace)
warning('off','all')
%% function 2 find possible MVR evs. TypeTrace is the criteria used.
%% Normally is a trace (For ex. traceDt that includes the whole synapse)
%% but also can be the Amplitude of the fit (A) while the events were localized.
%% the program uses the mean or median of ONLY synchronous and asynchronous events.
%% in the same way, the field added to perSyn is median or mean of only synchronous and asynchronous


% typeTrace : 'traceDt', 'traceLoc', 'A'


% if  isfield (Data, 'ROI')
%     traces = Data.ROI.trace;
% else
%     traces = Data.perSyn.trace;
%     
% end

%% for detection doublets & triplets
level = 1.5;
%type = 2; % 1 = median, 2= mean, 
%results = struct; % zeros (length(Data.perSyn), 1);
minumasynch = 4;
%%

indxALLEvs = (1:length(Data.allEv));
% if evsforcriteria == 1
%     mskSynchAsynch = logical([Data.allEv.synch]);
% else

% end

mskSynch = [Data.allEv.synch];
mskAsynch = [Data.allEv.asynch2nd];
mskBoth = logical(mskSynch + mskAsynch);

if isfield(Data.allEv,'mvr')
    Data.allEv = rmfield(Data.allEv,'mvr');
end

if isfield(Data.allEv,'normAmpli')
    Data.allEv = rmfield(Data.allEv,'normAmpli');
end

[lType, ~] = size( [Data.allEv.(typeTrace)]);
nameField = ['normAmp'  typeTrace];
nameMVR = ['mvr'  typeTrace];

% if type == 1
%     operation = 'median';
% else
%     operation = 'mean';
% end
    
mameSynField = ['mean' nameField];
mameSynField2 = ['testNorm' typeTrace];

    
for syn = 1:length(Data.perSyn)
    synID = Data.perSyn(syn).synID;
   %% BOTH, synch and Asynch
    indxEvs = indxALLEvs(logical(([Data.allEv.synID] == synID) .* mskBoth));
    indxAsynch = indxALLEvs(logical(([Data.allEv.synID] == synID) .* mskAsynch));
    mskAsynchsyn = ismember(indxEvs, indxAsynch);
    %ampliEvs = zeros(length(indxEvs),1);
    
    if lType > 1
        pre = [Data.allEv(indxEvs).(typeTrace)];
        pre = mean(pre((15:20), :), 1);

        val = [Data.allEv(indxEvs).(typeTrace)];
        val =  val(21, :);
        ampliEvs = (val - pre)';
        ampliEvs = (ampliEvs - min(ampliEvs)) + 10;
        %disp('here')
    else
        %ampliEvs = [Data.allEv(indxEvs).(typeTrace)];
        %ampliEvs = ([Data.allEv(indxEvs).(typeTrace)] + [Data.allEv(indxEvs).c])./[Data.allEv(indxEvs).c];
        ampliEvs = [Data.allEv(indxEvs).(typeTrace)]./[Data.allEv(indxEvs).c];
    end
     %% test normality
    [Data.perSyn(syn).(mameSynField2),~] = adtest(ampliEvs);
    corrected = false;

    if Data.perSyn(syn).(mameSynField2) == 1
        ampliEvsC = ampliEvs;        
        while length(ampliEvsC) > 4
           [~, lmax] = max(ampliEvsC);
           ampliEvsC(lmax) = [];
           if ~adtest(ampliEvsC)
               corrected = true;
               if sum(mskAsynchsyn) >= minumasynch
                   [~, ~, ci, ~] = ttest(ampliEvs, mean(ampliEvs, 'omitnan'));
                   meanAsynch = mean(ampliEvs(mskAsynchsyn));
                   if meanAsynch > ci(1) && meanAsynch < ci(2)
                       valSingle = meanAsynch;
                   else
                       valSingle = mean(ampliEvs, 'omitnan');

                   end

               else
                   valSingle = mean(ampliEvs, 'omitnan');
               end



%                valSingle = mean(ampliEvsC, 'omitnan');  

%                if type == 1
%                    valSingle = median(ampliEvsC, 'omitnan');            
%                else
%                    valSingle = mean(ampliEvsC, 'omitnan');  
%                end
               break
           
           end

        end
        if ~corrected
            numBins = round((max(ampliEvs)- min(ampliEvs))/ ((2 * iqr(ampliEvs))/ (length(ampliEvs))^(1/3))) + 1;
            [y, x]= hist(ampliEvs, numBins);
            binWidth = x(2)-x(1);
            x = [(x(1) - binWidth), x, (x(end) + binWidth)];
            y = [0, y, 0];
            [~, lm] = max(y(1:round((length(y))/2)));
    %         if lm > 
    %             disp('here')
    %         end
            valSingle = x(lm);
            %disp('here')

        end




    else
        if sum(mskAsynchsyn) >= minumasynch
            [~, ~, ci, ~] = ttest(ampliEvs, mean(ampliEvs, 'omitnan'));
            meanAsynch = mean(ampliEvs(mskAsynchsyn));
            if meanAsynch > ci(1) && meanAsynch < ci(2)
                valSingle = meanAsynch;
            else
                valSingle = mean(ampliEvs, 'omitnan');
            end

        else
            valSingle = mean(ampliEvs, 'omitnan');
        end
%         if type == 1
%             valSingle = median(ampliEvs, 'omitnan');            
%         else
%             valSingle = mean(ampliEvs, 'omitnan');  
%         end
    end

    Data.perSyn(syn).(mameSynField) = mean(ampliEvs, 'omitnan');    
    normAmpliEvs = ampliEvs ./valSingle;    
    %disp('here')
    normAmpliEvsLogic = normAmpliEvs  >= level;
    
 
    for ev = 1:length(indxEvs)  
        
        Data.allEv(indxEvs(ev)).(nameMVR)= normAmpliEvsLogic(ev); 
        Data.allEv(indxEvs(ev)).(nameField) = normAmpliEvs(ev);
%         pre =  mean(Data.allEv(indxEvs(ev)).(typeTrace)(15:20));
%         ampliEvs(ev) = Data.allEv(indxEvs(ev)).(typeTrace)(21) - pre;
    end
    
     
    
    
    
    
end

for ev = 1:length(Data.allEv)
    if isempty(Data.allEv(ev).(nameMVR))
        Data.allEv(ev).(nameMVR) = NaN;
        Data.allEv(ev).(nameField) = NaN;
        
    end

end
warning('on','all')
end
