function all = normByRadius2(all, extlPolygon, PixelSize)%,field)
%% if extlPolygon == 1; limit AZ = Hull Polygon + 25 nm extra. better use it

%% april 2021. function with multiple functions (yes I know, not ideal): 
% 1) calculates eccentricity of the synapse
% 2) recalculates the area (useful when used over simulated data)
% 3) calculates angle in radians of every event (assuming center in the hull envelope center)
% 4) calculates the normalized distance to the center

if length([all.perSyn.synID])== length(unique([all.perSyn.synID])) && length([all.perSyn.synID]) == length(unique([all.allEv.synID]))
else
    all = corrDuplisynIDS(all);
end

%% syn

synIDs = unique([all.allEv.synID], 'stable');
numSyns = length(synIDs);
%synchMask = [all.allEv.synch];
indxAll = (1:1:length(all.allEv));
% all.allEv.norm1 = zeros(length(all.allEv),1);
% all.allEv.angle = zeros(length(all.allEv),1);
for syn = 1 : numSyns
    mskEvs = ismember([all.allEv.synID], synIDs(syn));
    %synID = all.perSyn(syn).synID;

    %% Hull Synch
    xall = [all.allEv(mskEvs).x];
    yall = [all.allEv(mskEvs).y];
    xCentEv = [min(xall)-5, max(xall)+5];
    indxEvSyn = indxAll(mskEvs);
    
    synchMask = [all.allEv(indxEvSyn).synch];
    xsynch = xall(synchMask);
    ysynch = yall(synchMask);
%     indxSynch = (1:1:length(indxEvSyn));
%     indxSynch = indxSynch (synchMask);
    
    if sum(synchMask) > 1
        indxHull = convhull(xsynch,ysynch);
        xHull = xsynch(indxHull);
        yHull = ysynch(indxHull);
        [centerF, all.perSyn(syn).areaSynch] = centroidHull(xHull,yHull);
        
  %% polygon plus 25 nm
       
        p = polybuffer(polyshape(xHull, yHull), (25/PixelSize));
      %%  
        xallCenter = xall - centerF(1);
        yallCenter = yall - centerF(2);
        
        
        [angle,dist2cent] = cart2pol(xallCenter,yallCenter);
        %% add ratio minRadius/maxRadius per synapse
        %1) test a number of angles:
        xHullCenter = xHull - centerF(1);
        yHullCenter = yHull - centerF(2);
        angles = linspace(0,1, 1000);
        dst = sqrt((max(abs(xHullCenter)))^2 + (max(abs(yHullCenter)))^2)+1;
        allRadius = struct;
        %figure(200), hold on, plot(xHullCenter, yHullCenter);
        for ang = 1:length(angles)
            [x1,y1] = pol2cart(angles(ang)*pi,dst);
            [x2,y2] = pol2cart((angles(ang)-1)*pi,dst);
            intrc = InterX([xHullCenter;yHullCenter],[[x1,x2];[y1,y2]]);
            allRadius(ang).length = pdist ([[intrc(1,1), intrc(2,1)];[intrc(1,2), intrc(2,2)]]);
            allRadius(ang).x = intrc(1,:);
            allRadius(ang).y = intrc(2,:);
            allRadius(ang).ang = angles(ang);
            %figure(200), plot(allRadius.x, allRadius.y);
        end
       
        [maxLength, lcMx] = max([allRadius.length]);
        %[~, lcmn] = min([allRadius.length]);
        % find the  line perpendicular to the max
        perpAng = allRadius(lcMx).ang + 0.5;
        [x1,y1] = pol2cart(perpAng*pi,dst);
        [x2,y2] = pol2cart((perpAng-1)*pi,dst);
        intrc = InterX([xHullCenter;yHullCenter],[[x1,x2];[y1,y2]]);
        perpLength = pdist ([[intrc(1,1), intrc(2,1)];[intrc(1,2), intrc(2,2)]]);
%         perpx = intrc(1,:);
%         perpy = intrc(2,:);
        all.perSyn(syn).eccentricity = perpLength/maxLength;
        
        
        %plot:
%         figure(200), clf;
%         figure(200), hold on, plot(xHullCenter, yHullCenter), plot(allRadius(lcMx).x, allRadius(lcMx).y), plot(allRadius(lcmn).x, allRadius(lcmn).y);
%         figure(200), plot(perpx, perpy), hold off;
%         maxAll = max(max(abs(xHullCenter)), max(abs(yHullCenter)));
%         figure(200), xlim([-maxAll maxAll]), ylim([-maxAll maxAll]), hold off;
        %%
        
        
   
        for ev = 1:length(indxEvSyn)
            
            all.allEv(indxEvSyn(ev)).tetha = angle(ev)/pi;
            %% in : events in are maximun 25 nm away from hull envelope (described in polygon p)
            all.allEv(indxEvSyn(ev)).in = inpolygon(xall(ev),yall(ev),p.Vertices(:,1),p.Vertices(:,2));
            
            if ismember(xall(ev), xHull) &&  extlPolygon == 0
                all.allEv(indxEvSyn(ev)).norm1 = 1;
            else

   
                %% function center-event line
                m = (yall(ev)-centerF(2))/(xall(ev)-centerF(1));
                b = yall(ev)- m*(xall(ev));
                yCentEv = (xCentEv.*m) + b;
                %% find all the sides that center-event line crosses; pick up the one with the shortest P-ev distance:
                if extlPolygon == 0
                    P = InterX([xHull;yHull],[xCentEv;yCentEv]);
                else
                    P =  InterX([p.Vertices(:,1)';p.Vertices(:,2)'],[xCentEv;yCentEv]);
                end
                
                
                distev2P = squareform(pdist([[xall(ev) yall(ev)]; P'])); 
                distev2P = distev2P(:,1); distev2P(distev2P == 0) = [];
                [~, minP] = min(distev2P);
                P = (P(:,minP))';
                all.allEv(indxEvSyn(ev)).norm1 = dist2cent(ev) / pdist([[centerF(1) centerF(2)]; [P(1) P(2)]]);

            end

        end
    else
        for ev = 1:length(indxEvSyn)
            all.allEv(indxEvSyn(ev)).norm1 = NaN;
            
        end
    end
    
end
end
