function Data = transTrace(DataName,rmvBckgFlag, mskFlag,recenterFlag, foldFlag,  transvFlag)
%% function to get trace from frame with ev.
%% additionally, it gives an estimate of the asymetry of such a trces (aysmetryTraces) only implemented when transvFlag == 1; 

%% recenterFlag is kind of useless, too slow
%% corrType : startegy 2 avoid negative integral. 1= substract smalles pnt of the 4 traces.
%% maskFlag. if 1, mask added. if 2, a field called areaMask added.  
corrType = 0; %2
integrals = 0; %% if 1, calculates the integral and the asymetry of the traces etc... maybe useless
asymetryByFit = 0;
trimmerFlag = 1; %% little mask applied to the trace when intergral are masured (thta is, when corrType = 2 & for all int and ratios etc...)



%%
if ischar(DataName)
    Data = evalin('base', DataName);
else
    Data = DataName;
end

if length([Data.perSyn.synID])== length(unique([Data.perSyn.synID])) && length([Data.perSyn.synID]) == length(unique([Data.allEv.synID]))
else
    Data = corrDuplisynIDS(Data);
end



for synNum = 1:length(Data.roi)
    synID = Data.roi(synNum).synID;
    mskSyn = [Data.allEv.synID] == synID;
    if mskFlag ~= 0
        mask = Data.roi(synNum).mask;
        areaMask =  sum(sum(mask));
    end

    if synNum == 1
        locEv = 0;

    end
    
    for ev = 1:sum(mskSyn)
%         if locEv + ev > 7305
%             disp ('here')
%         end
        frameN = Data.allEv(locEv + ev).frame;
        
        
     
        frame = Data.roi(synNum).movie(:,:,frameN);
        Befframe = Data.roi(synNum).movie(:,:,frameN - 1);
        if (frameN + 1) < length(Data.roi(synNum).movie)
            Aftframe1 = Data.roi(synNum).movie(:,:,frameN + 1);
            Aftframe2 = Data.roi(synNum).movie(:,:,frameN + 2);
        else
            Aftframe1 = NaN(size(frame));
            Aftframe2 = NaN(size(frame));
        end
      
        if mskFlag == 1
            frame = double(frame) .* mask;
            Befframe = double(Befframe) .* mask;
            Aftframe1 = double(Aftframe1) .* mask;
            Aftframe2 = double(Aftframe2) .* mask;
        end

        
        
            
        if rmvBckgFlag ~= 0
            NumFrame2Remv = 12;% 6
            if (frameN - NumFrame2Remv ) <= 0
                NumFrame2Remv = frameN -1;
            end
            bckg =double(median(Data.roi(synNum).movie(:, :, (frameN - NumFrame2Remv):(frameN - 2)), 3));
            
            %bckg = double(median(Data.roi(synNum).movie(:, :, (frameN - NumFrame2Remv):(frameN - 2)), 3));
            if mskFlag == 1
                bckg = double(bckg) .* mask;
            end
            
            if rmvBckgFlag == 1
                frame = double(frame) - bckg;
                Befframe = double(Befframe) - bckg;
                Aftframe1 = double(Aftframe1) - bckg;
                Aftframe2 = double(Aftframe2) - bckg;
            end


            
        end

        x = round(Data.allEv(locEv + ev).x);
        y = round(Data.allEv(locEv + ev).y);


        limsx = (x - 15):1:(x + 15);
        limsy = (y - 15):1:(y + 15);
        
        limsx(limsx < 1) = 1;
        limsx(limsx > 41) = 41;
        limsy(limsy < 1) = 1;
        limsy(limsy > 41) = 41;
        %% x and y traces at the center
        
        xtrc = double(frame(limsy,x));
        ytrc = double((frame(y,limsx))');
        
        if recenterFlag == 1
            [~,fitParams1x]= gaussFit3CoCenter(xtrc);
            if round(fitParams1x(2)) ~= 16 && round(fitParams1x(2)) < 19 && round(fitParams1x(2)) > 13
                ycorr = round(fitParams1x(2))-16;
                y = y + ycorr;
                limsy = (y - 15):1:(y + 15);
                limsy(limsy < 1) = 1;
                limsy(limsy > 41) = 41;
                %disp('here')
%                 xtrcBef = recenterTrace(xtrcBef, round(fitParams1x(2)), 16);
%                 xtrc = recenterTrace(xtrc, round(fitParams1x(2)), 16);
%                 xtrcAft1 = recenterTrace(xtrcAft1, round(fitParams1x(2)), 16);
%                 xtrcAft2 = recenterTrace(xtrcAft2, round(fitParams1x(2)), 16);
                
            end
            
             [~,fitParams1y]= gaussFit3CoCenter(ytrc);
            if round(fitParams1y(2)) ~= 16 && round(fitParams1y(2)) < 19 && round(fitParams1y(2)) > 13
                xcorr = round(fitParams1y(2))-16;
                x = x + xcorr;
                limsx = (x - 15):1:(x + 15);
                limsx(limsx < 1) = 1;
                limsx(limsx > 41) = 41;
%                 ytrcBef = recenterTrace(ytrcBef, round(fitParams1y(2)), 16);
%                 ytrc = recenterTrace(ytrc, round(fitParams1y(2)), 16);
%                 ytrcAft1 = recenterTrace(ytrcAft1, round(fitParams1y(2)), 16);
%                 ytrcAft2 = recenterTrace(ytrcAft2, round(fitParams1y(2)), 16);
                
            end
        xtrc = double(frame(limsy,x));
        ytrc = double((frame(y,limsx))');
        end
        
        
        
        
        %%
        
        xtrcBef = double(Befframe(limsy,x));
        %xtrc = double(frame(limsy,x));
        xtrcAft1 = double(Aftframe1(limsy,x));
        xtrcAft2 = double(Aftframe2(limsy,x));
        
        ytrcBef = double((Befframe(y,limsx))');
        %ytrc = double((frame(y,limsx))');
        ytrcAft1 = double((Aftframe1(y,limsx))');
        ytrcAft2 = double((Aftframe2(y,limsx))');
        
        if transvFlag == 1
            [t1trcBef, t2trcBef] = getTransversal(x,y, Befframe, 16);
            [t1trc, t2trc] = getTransversal(x,y, frame, 16);
            
            
            
            
            if asymetryByFit == 1
                [maxAsym, medAsym, sigma0, sigma1] = aysmetryTraces(xtrc, ytrc, t1trc, t2trc, 1,asymetryByFit);
            else
                [maxAsym, medAsym] = aysmetryTraces(xtrc, ytrc, t1trc, t2trc, 1,asymetryByFit);
            end
            
            if sum(sum(isnan(Aftframe1))) == 0
                [t1trcAft1, t2trcAft1] = getTransversal(x,y, Aftframe1, 16);
                [t1trcAft2, t2trcAft2] = getTransversal(x,y, Aftframe2, 16);
                
            else
                t1trcAft1 = xtrcAft1;
                t2trcAft1 = xtrcAft1;
                t1trcAft2 = xtrcAft1;
                t2trcAft2 = xtrcAft1;
            end
     
            
        end
        
        
        
        
        if recenterFlag == 2
            %(locEv + ev)
            %if (locEv + ev) ~= 429 && (locEv + ev) ~= 4358  && (locEv + ev) ~= 5264
            
            [fitParams1x, ~]= gaussFit3CoCenter(xtrcAft1);
            if round(fitParams1x(2)) ~= 16 && round(fitParams1x(2)) < 20 && round(fitParams1x(2)) > 12
                [fitParams2x, ~]= gaussFit3CoCenter(xtrcAft2);
                if round(fitParams2x(2)) ~= 16 && abs(round(fitParams2x(2)) - round(fitParams1x(2))) < 2
                    xtrcAft1 = recenterTrace(xtrcAft1, round(fitParams1x(2)), 16);
                    xtrcAft2 = recenterTrace(xtrcAft2, round(fitParams2x(2)), 16);
                end
            end
            
            [fitParams1y, ~]= gaussFit3CoCenter(ytrcAft1);
            if round(fitParams1y(2)) ~= 16 && round(fitParams1y(2)) < 20 && round(fitParams1y(2)) > 12
                [fitParams2y, ~]= gaussFit3CoCenter(ytrcAft2);
                if round(fitParams2y(2)) ~= 16 && abs(round(fitParams2y(2)) - round(fitParams1y(2))) < 2
                    ytrcAft1 = recenterTrace(ytrcAft1, round(fitParams1y(2)), 16);
                    ytrcAft2 = recenterTrace(ytrcAft2, round(fitParams2y(2)), 16);
                end
            end
            %end
            
            
            
        end
        
        if rmvBckgFlag == 2
            [bckx, ~]= gaussFit3CoCenter(bckg(limsy,x), [0 16 0 0]);
            [bcky, ~]= gaussFit3CoCenter(bckg(y,limsx)', [0 16 0 0]);
            
            
            xtrcBef = xtrcBef - bckx;%Befframe(limsy,x);
            xtrc = xtrc - bckx;%frame(limsy,x);
            xtrcAft1 = xtrcAft1 - bckx;%Aftframe1(limsy,x);
            xtrcAft2 = xtrcAft2 - bckx;%Aftframe2(limsy,x);
            
            ytrcBef = ytrcBef - bcky;%(Befframe(y,limsx))';
            ytrc = ytrc - bcky;%(frame(y,limsx))';
            ytrcAft1 = ytrcAft1 - bcky;%(Aftframe1(y,limsx))';
            ytrcAft2 = ytrcAft2 - bcky;%(Aftframe2(y,limsx))';
            
            
            if transvFlag == 1
                [bckt1trc, bckt2trc] = getTransversal(x,y, bckg, 16);
                [bckt1trc, ~]= gaussFit3CoCenter(bckt1trc, [0 16 0 0]);
                [bckt2trc, ~]= gaussFit3CoCenter(bckt2trc, [0 16 0 0]);
                
                t1trcBef = t1trcBef - bckt1trc;
                t2trcBef =  t2trcBef - bckt2trc;
                
                t1trc = t1trc - bckt1trc;
                t2trc = t2trc - bckt2trc;
                
                if sum(sum(isnan(Aftframe1))) == 0
                    t1trcAft1 = t1trcAft1 - bckt1trc;
                    t1trcAft2 = t1trcAft2 - bckt1trc;
                    
                    t2trcAft1 = t2trcAft1 - bckt2trc;
                    t2trcAft2 = t2trcAft2 - bckt1trc;
                    
                   
                else
                    t1trcAft1 = xtrcAft1;
                    t2trcAft1 = xtrcAft1;
                    t1trcAft2 = xtrcAft1;
                    t2trcAft2 = xtrcAft1;
                end
                
                
                
                
            
                
            end

            
        end
        if foldFlag == 1
            xtrcBef = [(xtrcBef(1:15) + flip(xtrcBef(17:end)))./2; xtrcBef(16);flip((xtrcBef(1:15) + flip(xtrcBef(17:end)))./2)];
            xtrc = [(xtrc(1:15) + flip(xtrc(17:end)))./2; xtrc(16);flip((xtrc(1:15) + flip(xtrc(17:end)))./2)];
            xtrcAft1 = [(xtrcAft1(1:15) + flip(xtrcAft1(17:end)))./2; xtrcAft1(16);flip((xtrcAft1(1:15) + flip(xtrcAft1(17:end)))./2)];
            xtrcAft2 = [(xtrcAft2(1:15) + flip(xtrcAft2(17:end)))./2; xtrcAft2(16);flip((xtrcAft2(1:15) + flip(xtrcAft2(17:end)))./2)];
            
            ytrcBef = [(ytrcBef(1:15) + flip(ytrcBef(17:end)))./2; ytrcBef(16);flip((ytrcBef(1:15) + flip(ytrcBef(17:end)))./2)];
            ytrc = [(ytrc(1:15) + flip(ytrc(17:end)))./2; ytrc(16);flip((ytrc(1:15) + flip(ytrc(17:end)))./2)];
            ytrcAft1 = [(ytrcAft1(1:15) + flip(ytrcAft1(17:end)))./2; ytrcAft1(16);flip((ytrcAft1(1:15) + flip(ytrcAft1(17:end)))./2)];
            ytrcAft2 = [(ytrcAft2(1:15) + flip(ytrcAft2(17:end)))./2; ytrcAft2(16);flip((ytrcAft2(1:15) + flip(ytrcAft2(17:end)))./2)];
            
            if transvFlag == 1
                t1trcBef = [(t1trcBef(1:15) + flip(t1trcBef(17:end)))./2; t1trcBef(16);flip((t1trcBef(1:15) + flip(t1trcBef(17:end)))./2)];
                t1trc = [(t1trc(1:15) + flip(t1trc(17:end)))./2; t1trc(16);flip((t1trc(1:15) + flip(t1trc(17:end)))./2)];
                t1trcAft1 = [(t1trcAft1(1:15) + flip(t1trcAft1(17:end)))./2; t1trcAft1(16);flip((t1trcAft1(1:15) + flip(t1trcAft1(17:end)))./2)];
                t1trcAft2 = [(t1trcAft2(1:15) + flip(t1trcAft2(17:end)))./2; t1trcAft2(16);flip((t1trcAft2(1:15) + flip(t1trcAft2(17:end)))./2)];
                
                t2trcBef = [(t2trcBef(1:15) + flip(t2trcBef(17:end)))./2; t2trcBef(16);flip((t2trcBef(1:15) + flip(t2trcBef(17:end)))./2)];
                t2trc = [(t2trc(1:15) + flip(t2trc(17:end)))./2; t2trc(16);flip((t2trc(1:15) + flip(t2trc(17:end)))./2)];
                t2trcAft1 = [(t2trcAft1(1:15) + flip(t2trcAft1(17:end)))./2; t2trcAft1(16);flip((t2trcAft1(1:15) + flip(t2trcAft1(17:end)))./2)];
                t2trcAft2 = [(t2trcAft2(1:15) + flip(t2trcAft2(17:end)))./2; t2trcAft2(16);flip((t2trcAft2(1:15) + flip(t2trcAft2(17:end)))./2)];
            end
            
            
            
            
        end
        Data.allEv(locEv + ev).tTrcx = xtrc;
        Data.allEv(locEv + ev).tTrcy = ytrc;

        if transvFlag ~= 1
            Data.allEv(locEv + ev).tTrcxyBef = ((xtrcBef + ytrcBef) ./ 2);%(Befframe(limsy,x) + (Befframe(y,limsx))')./2;
            Data.allEv(locEv + ev).tTrcxy = ((xtrc + ytrc) ./ 2);%(frame(limsy,x) + (frame(y,limsx))')./2;
            Data.allEv(locEv + ev).tTrcxyAft1 = ((xtrcAft1 + ytrcAft1) ./ 2);%(Aftframe1(limsy,x) + (Aftframe1(y,limsx))')./2;
            Data.allEv(locEv + ev).tTrcxyAft2 = ((xtrcAft2 + ytrcAft2) ./ 2);%(Aftframe2(limsy,x) + (Aftframe2(y,limsx))')./2;
        else
            Data.allEv(locEv + ev).tTrcxyBef = ((xtrcBef + ytrcBef + t1trcBef + t2trcBef) ./ 4);%(Befframe(limsy,x) + (Befframe(y,limsx))')./2;
            Data.allEv(locEv + ev).tTrcxy = ((xtrc + ytrc + t1trc + t2trc) ./ 4);%(frame(limsy,x) + (frame(y,limsx))')./2;
            Data.allEv(locEv + ev).tTrcxyAft1 = ((xtrcAft1 + ytrcAft1 + t1trcAft1 + t2trcAft1) ./ 4);%(Aftframe1(limsy,x) + (Aftframe1(y,limsx))')./2;
            Data.allEv(locEv + ev).tTrcxyAft2 = ((xtrcAft2 + ytrcAft2 + t1trcAft2 + t2trcAft2) ./ 4);%(Aftframe2(limsy,x) + (Aftframe2(y,limsx))')./2;
        end
        if mskFlag == 1
            Data.allEv(locEv + ev).tTrcxyBef = Data.allEv(locEv + ev).tTrcxyBef./ areaMask;%(xtrcBef + ytrcBef + t1trcBef + t2trcBef) ./ 4;%(Befframe(limsy,x) + (Befframe(y,limsx))')./2;
            Data.allEv(locEv + ev).tTrcxy = Data.allEv(locEv + ev).tTrcxy./ areaMask;%(xtrc + ytrc + t1trc + t2trc) ./ 4;%(frame(limsy,x) + (frame(y,limsx))')./2;
            Data.allEv(locEv + ev).tTrcxyAft1 = Data.allEv(locEv + ev).tTrcxyAft1./ areaMask;%(xtrcAft1 + ytrcAft1 + t1trcAft1 + t2trcAft1) ./ 4;%(Aftframe1(limsy,x) + (Aftframe1(y,limsx))')./2;
            Data.allEv(locEv + ev).tTrcxyAft2 = Data.allEv(locEv + ev).tTrcxyAft2./ areaMask;%v(xtrcAft2 + ytrcAft2 + t1trcAft2 + t2trcAft2) ./ 4;%
            
            
        end
        if mskFlag == 2
            Data.allEv(locEv + ev).areaMask = areaMask;
            
        end
        
        
        %% Add short mask
        if trimmerFlag ==1
            borderL = 12;
            trimmer = true(size(xtrc));
            trimmer(1:borderL) = false;
            trimmer(length(trimmer)-borderL+1:end) = false;
        else
            trimmer = true(size(xtrc));
%             Data.allEv(locEv + ev).tTrcxyBef = Data.allEv(locEv + ev).tTrcxyBef(trimmer);%(xtrcBef + ytrcBef + t1trcBef + t2trcBef) ./ 4;%(Befframe(limsy,x) + (Befframe(y,limsx))')./2;
%             Data.allEv(locEv + ev).tTrcxy = Data.allEv(locEv + ev).tTrcxy(trimmer);%(xtrc + ytrc + t1trc + t2trc) ./ 4;%(frame(limsy,x) + (frame(y,limsx))')./2;
%             Data.allEv(locEv + ev).tTrcxyAft1 = Data.allEv(locEv + ev).tTrcxyAft1(trimmer);%(xtrcAft1 + ytrcAft1 + t1trcAft1 + t2trcAft1) ./ 4;%(Aftframe1(limsy,x) + (Aftframe1(y,limsx))')./2;
%             Data.allEv(locEv + ev).tTrcxyAft2 = Data.allEv(locEv + ev).tTrcxyAft2(trimmer);%v(xtrcAft2 + ytrcAft2 + t1trcAft2 + t2trcAft2) ./ 4;%
        %
        end
        
        %% integrals, correct for not having negative integrals. 

        if corrType ~= 0 && integrals ~= 0
            
            
            if corrType == 1 
                corr = min([min(Data.allEv(locEv + ev).tTrcxyBef), min(Data.allEv(locEv + ev).tTrcxy),...
                min(Data.allEv(locEv + ev).tTrcxyAft1), min(Data.allEv(locEv + ev).tTrcxyAft2)]);
                if corr < 0
                    Data.allEv(locEv + ev).tTrcxyBef = Data.allEv(locEv + ev).tTrcxyBef - corr;
                    Data.allEv(locEv + ev).tTrcxy = Data.allEv(locEv + ev).tTrcxy - corr;
                    Data.allEv(locEv + ev).tTrcxyAft1 = Data.allEv(locEv + ev).tTrcxyAft1 - corr;
                    Data.allEv(locEv + ev).tTrcxyAft2 = Data.allEv(locEv + ev).tTrcxyAft2 - corr;
                end
            else
                ints(1,1) = trapz(Data.allEv(locEv + ev).tTrcxyBef(trimmer) );
                ints(2,1) = trapz(Data.allEv(locEv + ev).tTrcxy(trimmer));
                ints(3,1) = trapz(Data.allEv(locEv + ev).tTrcxyAft1(trimmer));
                ints(4,1) = trapz(Data.allEv(locEv + ev).tTrcxyAft2(trimmer));
                mskInts = ints < 0;
                
                if sum(mskInts) ~= 0
                    corr = ones(size(Data.allEv(locEv + ev).tTrcxyBef)).*((abs(min(ints)) + 1)/(length(Data.allEv(locEv + ev).tTrcxyBef(trimmer))-1));
                    Data.allEv(locEv + ev).tTrcxyBef = (Data.allEv(locEv + ev).tTrcxyBef + corr);
                    Data.allEv(locEv + ev).tTrcxy = (Data.allEv(locEv + ev).tTrcxy + corr);
                    Data.allEv(locEv + ev).tTrcxyAft1 = (Data.allEv(locEv + ev).tTrcxyAft1 + corr);
                    Data.allEv(locEv + ev).tTrcxyAft2 = (Data.allEv(locEv + ev).tTrcxyAft2 + corr);
                end

                
            end
        end
        %%
       if transvFlag == 1 && integrals == 1
            Data.allEv(locEv + ev).maxAsym = maxAsym;
            Data.allEv(locEv + ev).medAsym = medAsym;
            if asymetryByFit == 1
                Data.allEv(locEv + ev).sigma0 = sigma0;
                Data.allEv(locEv + ev).sigma1 = sigma1;
            end
       end
%        
%        bef = trapz(Data.allEv(locEv + ev).tTrcxyBef(trimmer));
%        one = trapz(Data.allEv(locEv + ev).tTrcxy(trimmer));
%        two = trapz(Data.allEv(locEv + ev).tTrcxyAft1(trimmer));
%        three = trapz(Data.allEv(locEv + ev).tTrcxyAft2(trimmer));
%        
       %% integrals....
       if integrals == 1 
           Data.allEv(locEv + ev).intTrcBef = trapz(Data.allEv(locEv + ev).tTrcxyBef(trimmer));
           Data.allEv(locEv + ev).intTrcxy = trapz(Data.allEv(locEv + ev).tTrcxy(trimmer));%(one-bef)/bef;%
           Data.allEv(locEv + ev).intTrcxyAft1 = trapz(Data.allEv(locEv + ev).tTrcxyAft1(trimmer));%(two-bef)/bef;%
           Data.allEv(locEv + ev).intTrcxyAft2 = trapz(Data.allEv(locEv + ev).tTrcxyAft2(trimmer));%(three-bef)/bef;%



           Data.allEv(locEv + ev).ratioInt = Data.allEv(locEv + ev).intTrcxy/(Data.allEv(locEv + ev).intTrcxyAft1);
           Data.allEv(locEv + ev).ratioInt2 = Data.allEv(locEv + ev).intTrcxyAft1/(Data.allEv(locEv + ev).intTrcxyAft2);

            %Data.allEv(locEv + ev).difussed = Data.allEv(locEv + ev).ratioInt < 2 && Data.allEv(locEv + ev).ratioInt > 1;
           Data.allEv(locEv + ev).as= Data.allEv(locEv + ev).ratioInt < 1.5 && Data.allEv(locEv + ev).ratioInt > 0.5;
           Data.allEv(locEv + ev).as2= Data.allEv(locEv + ev).ratioInt < 1.5 && Data.allEv(locEv + ev).ratioInt > 0.5 &&...
                                       Data.allEv(locEv + ev).ratioInt2 < 1.5 && Data.allEv(locEv + ev).ratioInt2 > 0.5;
       end
        

     
        
        
    end
    locEv = locEv + ev;
        

        
       
        
        %Data.traces(ev).yAft19 = (Data.roi(synNum).movie(y,limsx,frame + 19))';
end
   
 

%% 
if ~isempty(Data.allEv)
    msk = isnan([Data.allEv.tTrcxyAft2]);
    msk = msk(:,1);
    indx = (1:length(Data.allEv))';
    indx = indx(msk);
    for kk = 1:length(indx)
        Data.allEv(indx(kk)).tTrcxyAft1 = [];
        Data.allEv(indx(kk)).tTrcxyAft2 = [];
    end
end

end