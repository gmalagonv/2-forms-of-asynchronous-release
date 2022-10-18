function [ Cluster ] = ClusterVesicFusion( Data,SynN,MaxDist,figure,indicesAccpEv,varargin)
%% if events do not come from a Data file, use varargin. Put in Data a any file with the pixelSize info. 
% SynN = 1;
MaxDist = MaxDist/Data.InfoExp.PixelSize; %in pixels 
if nargin == 4 %GM Nov 2018
    Cluster.NEvents = Data.synapsis(SynN).NlocEvents;
    XY = [ [Data.synapsis(SynN).events_localization.Event.x]' [Data.synapsis(SynN).events_localization.Event.y]' ];
else
    if isempty(varargin)
        Cluster.NEvents = length(indicesAccpEv); %GM Nov 2018
        XY = [ [Data.synapsis(SynN).events_localization.Event(indicesAccpEv).x]' [Data.synapsis(SynN).events_localization.Event(indicesAccpEv).y]' ]; %GM Nov 2018
    else
        XY = varargin{1};
        Cluster.NEvents = length(varargin{1});
    end
end
    

if Cluster.NEvents >= 2
    %% Hierarchical Clustering

    Y = pdist(XY);
    %Y = squareform(Y); %% corrected in May 11, 2021 (late :()
    if length(Y) < 3
        disp('here')
    end
        
    Z = linkage(Y,'complete');

    T = cluster(Z,'cutoff',MaxDist,'criterion','distance');

    for i = 1:max(T)

        Cluster.cluster(i).loc = XY(T==i, :);

        if size(Cluster.cluster(i).loc,1)~=1
            [~,ctrs] = kmeans(Cluster.cluster(i).loc,1);
            Cluster.cluster(i).centroid = ctrs;

            c = Cluster.cluster(i).centroid;
            c(2:1+length(Cluster.cluster(i).loc),:) = Cluster.cluster(i).loc;
            dist = pdist(c);
            Cluster.cluster(i).maxRad = max(dist(1:size(Cluster.cluster(i).loc,1)));

        else

            Cluster.cluster(i).centroid = Cluster.cluster(i).loc;
            Cluster.cluster(i).maxRad = 0;
        end
    end
    
else
    
    if Cluster.NEvents == 1
        
        Cluster.cluster.loc = XY;
        Cluster.cluster.centroid = XY;
        Cluster.cluster.maxRad = 0;
        
    else
        Cluster = [];
    end
end

if figure == 1
    H = dendrogram(Z,'ColorThreshold','default');
    set(H,'LineWidth',2)

    %% Legend and Labels.

    % hTitle  = title ('Clusters');
    hTitle  = title (sprintf('%d Fusion/%d Clusters', Cluster.NEvents,length(Cluster.cluster)));
    hXLabel = xlabel('\bfEvent', 'FontName', 'AvantGarde', 'FontSize', 17);
    hYLabel = ylabel('\bflength \rm\it(pixels)', 'FontName', 'AvantGarde', 'FontSize', 17);
    
    %% Adjust Font and Axes Properties.


    set( gca                       , ...
        'FontName'   , 'Helvetica' );
    set([hTitle, hXLabel, hYLabel ], ...
        'FontName'   , 'AvantGarde');

    set([hXLabel, hYLabel ]  , ...
        'FontSize'   , 12          );
    set( hTitle                    , ...
        'FontSize'   , 16          , ...
        'FontWeight' , 'bold'      );

    set(gca, ...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'off'      , ...
      'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off'      , ...
      'XGrid'       , 'off'      , ...
      'XColor'      , [.3 .3 .3], ...
      'YColor'      , [.3 .3 .3]);
    
    
    
    
    %% Draw Cut Line.
    
    line([-100000 100000],[MaxDist MaxDist],'Color','r','LineStyle','--','LineWidth',2);

    %% Export to EPS


    set(gcf, 'PaperPositionMode', 'auto');
    eval ([sprintf('print -depsc2 Dendogram%d.eps;', length(dir( fullfile('Dendogram*.eps') ))+1)]);
%     eval ([sprintf('print -depsc2 DendogramSynN%d.eps;', SynN)]);
        
    
end
