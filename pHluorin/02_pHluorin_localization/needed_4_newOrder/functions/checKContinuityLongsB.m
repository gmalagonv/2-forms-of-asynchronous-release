  function locEv = checKContinuityLongsB(locEv)

initNumMult = sum([locEv.long]);  
disp ([num2str(initNumMult) ' multiple locs. per frame found'])

varRepet = 1;
while varRepet ~= 0
    totNumEv = length(locEv);
    % maskPre = true (totNumEv,1);
    % maskPost = true (totNumEv,1);
    maskTot = true(totNumEv,1);


    [ONduration, ~, indexesON]  = intervalsOnOff([locEv.long]);
    ONduration = ONduration';
    xy_pstdAll = (sqrt((([locEv.x_pstd]).^2)+(([locEv.y_pstd]).^2) ))';

    frames = [locEv.frame]';
    %xy_pstd = xy_pstd.*([locEv.long])';

    for i = 1:length(ONduration)
        preIndx = min (indexesON{i}) - 1;
        postIndx = max (indexesON{i}) + 1;
        %framesInvolved = frames(indexesON{i});
        numFrames = length(unique(frames(indexesON{i}))); 
        indx = indexesON{i};
        frames2check = [min(unique(frames(indexesON{i}))) max(unique(frames(indexesON{i})))]; % only the 1st and the last
        for f = 1:2 %% pre and post
            indx2check = indx(frames(indx) ==  frames2check(f));
            %% pre
            if f == 1
                localMaskPre = false(size(indx2check));
                %kIndx = preIndx;

                if (frames(preIndx) + 1) == frames(min(indx2check))
                    for m = 1:length(indx2check)
                        distance = pdist_several(locEv(preIndx).x, locEv(preIndx).y, locEv(indx2check(m)).x, locEv(indx2check(m)).y);
                        %dist(indx2check(m),1) = distance(3);
                        localMaskPre(m) = (distance(3) < xy_pstdAll(preIndx));
                    end

                if numFrames > 1 && sum(localMaskPre) == 1
                    maskTot(indx2check(~localMaskPre)) = false;
                    locEv(indx2check(localMaskPre)).long = false;
                end
                end            

            else
            %% post
                localMaskPost = false(size(indx2check));
                if (frames(postIndx) - 1) == frames(max(indx2check))
                    for m = 1:length(indx2check)
                        distance = pdist_several(locEv(postIndx).x, locEv(postIndx).y, locEv(indx2check(m)).x, locEv(indx2check(m)).y);
                        %dist(indx2check(m),1) = distance(3);
                        localMaskPost(m) = (distance(3) < xy_pstdAll(postIndx));
                    end
                if numFrames > 1 && sum(localMaskPost) == 1
                    maskTot(indx2check(~localMaskPost)) = false;
                    locEv(indx2check(localMaskPost)).long = false;
                end
                end            
            end
        end

        if numFrames == 1
            localMaskTot = logical(localMaskPre + localMaskPost);
            if sum(localMaskTot) == 1
                maskTot(indx(~localMaskTot)) = false;
                locEv(indx(localMaskTot)).long = false;
            end
        else

        end   

    end
    varRepet = sum(~maskTot);    
    locEv = locEv(maskTot);
end
finNumMult = sum([locEv.long]);
disp ([num2str(round((initNumMult-finNumMult)*100/initNumMult)) '% removed by continuity'])

end
