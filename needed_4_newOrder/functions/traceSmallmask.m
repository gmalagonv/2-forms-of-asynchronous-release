function Data = traceSmallmask(DataName, radius, rmvBckgFlag, noiseFlag)

%%radius = 5;

if ischar(DataName)
    Data = evalin('base', DataName);
else
    Data = DataName;
end


lengthInterv = length(Data.InfoExp.descartFrame); 
numFramBef = lengthInterv;
numFramAft = (lengthInterv*2) - 2;

%numFramTot = numFramBef + numFramAft + 1;
area = sum(sum(double(bsxfun(@plus, ((1:41) - 21).^2, (transpose(1:41) - 21).^2) < radius^2)));


if length([Data.perSyn.synID])== length(unique([Data.perSyn.synID])) && length([Data.perSyn.synID]) == length(unique([Data.allEv.synID]))
else
    Data = corrDuplisynIDS(Data);
end



for synNum = 1:length(Data.roi)
    synID = Data.roi(synNum).synID;
    mskSyn = [Data.allEv.synID] == synID;
    

    if synNum == 1
        locEv = 0;

    end
    
    for ev = 1:sum(mskSyn)
        
        frameN = Data.allEv(locEv + ev).frame;
        x = round(Data.allEv(locEv + ev).x);
        y = round(Data.allEv(locEv + ev).y);
        
        if noiseFlag == 1
            x = x + 10;
            y = y + 10;
            if x > 41
                x = x - 20;
            end
            if y > 41
                y = y - 20;
            end
            
            
        end
             
        maskEv = bsxfun(@plus, ((1:41) - x).^2, (transpose(1:41) - y).^2) < radius^2;
        
        indxMov = (frameN - numFramBef:frameN + numFramAft);
        indxMov(indxMov < 1) = 1;
        indxMov(indxMov > length(Data.roi(synNum).movie)) = length(Data.roi(synNum).movie);
        temp = double(maskEv).*double(Data.roi(synNum).movie(:,:,indxMov));
        
        
        if rmvBckgFlag == 1
            NumFrame2Remv = 6;
            if (frameN - NumFrame2Remv ) <= 0
                NumFrame2Remv = frameN -1;
            end
            bckg = double(maskEv) .* mean(Data.roi(synNum).movie(:, :, (frameN - NumFrame2Remv):(frameN - 2)), 3);
            temp = temp- bckg;

        end
        if noiseFlag == 0
            Data.allEv(locEv + ev).traceLoc = squeeze(double(sum(sum(temp))/area));
            locPeak = 21;
            Data.allEv(locEv + ev).traceLocPeak = Data.allEv(locEv + ev).traceLoc(locPeak) - mean(Data.allEv(locEv + ev).traceLoc((locPeak-6):(locPeak-1)));
            
            Data.allEv(locEv + ev).traceLocrat1 = Data.allEv(locEv + ev).traceLoc(locPeak)/Data.allEv(locEv + ev).traceLoc(locPeak + 1);
            Data.allEv(locEv + ev).traceLocrat2 = Data.allEv(locEv + ev).traceLoc(locPeak + 1)/Data.allEv(locEv + ev).traceLoc(locPeak + 2);
           
        else
            Data.allEv(locEv + ev).traceLocN = squeeze(double(sum(sum(temp))/area)); 
        end

        
        
        
  
        %locEv(evIndxSyn(ev)).traceDt = (traceDt(frameN(ev)-numFramBef:frameN(ev)+numFramAft));
        
        
        
        
        
        %%
        
        
     
%         frame = Data.roi(synNum).movie(:,:,frameN);
%         Befframe = Data.roi(synNum).movie(:,:,frameN - 1);
%         
%         if (frameN + 1) < length(Data.roi(synNum).movie)
%             Aftframe1 = Data.roi(synNum).movie(:,:,frameN + 1);
%             Aftframe2 = Data.roi(synNum).movie(:,:,frameN + 2);
%         else
%             Aftframe1 = NaN(size(frame));
%             Aftframe2 = NaN(size(frame));
%         end
%         if mskFlag == 1
%             area =  sum(sum(Data.roi(synNum).mask));
%             frame = double(frame) .* Data.roi(synNum).mask;
%             Befframe = double(Befframe) .* Data.roi(synNum).mask;
%             Aftframe1 = double(Aftframe1) .* Data.roi(synNum).mask;
%             Aftframe2 = double(Aftframe2) .* Data.roi(synNum).mask;
% 
%         end


        
        
        
    end
    locEv = locEv + ev;
        

        
       
        
end
   
 

%% 
% msk = isnan([Data.allEv.tTrcxyAft2]');
% msk = msk(:,1);
% indx = (1:length(Data.allEv))';
% indx = indx(msk);
% for kk = 1:length(indx)
%     Data.allEv(indx(kk)).tTrcxyAft1 = [];
%     Data.allEv(indx(kk)).tTrcxyAft2 = [];
% end

end