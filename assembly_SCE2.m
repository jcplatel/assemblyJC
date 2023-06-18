%default values
%MinPeakDistancesce=5 frames= 500ms
% synchronous_frames=2
%MinPeakDistance=5 frames= 500ms
%sce_n_cells_threshold=5 or 5% ???   seem best results between 10 to 15
%cells  should we simulate ???

%to do:
% nice raster, sorted raster, nice SCE on sorted raster , assembly surrounded, SCE organized by assembly 
clear
close all
%disp(daytime)
%for kmeans_surrogate=[10 10 100 100 200 200 1000 1000]
tic
%for sce_n_cells_threshold=10:20%[5,10,15,20]%5:30%
%for z=1:5
openingnwb


%% Load settings
MinPeakDistancesce=3 ;
%synchronous_frames=2;
MinPeakDistance=3;

%imaging_sampling_rate=15.7;
synchronous_frames=round(0.2*sampling_rate,0); %200ms *sampling rate

% PathSave='/Users/platel/Desktop/exp/analysis/';
PathSave='/Users/platel/Desktop/exp/analysis/';

 daytime = datestr(now,'yy_mm_dd_HH_MM_SS');
 namefull=[PathSave daytime name '/'];
 mkdir (namefull)    % make folder for saving analysis
disp(['make new folder ' namefull])

%% Load current data
% 
% plane0=double(plane0(iscell0(:,1)>0,:));
% plane1=double(plane1(iscell1(:,1)>0,:));
% 
Tr1b=F;

% % Tr1b=plane1;

%% preprocessing

%smooth speed
%speedraw=speed;
speedraw=speed;
 speed =smoothdata(speed,'gaussian',50);
%speed already smoothed


[NCell,Nz] = size(Tr1b);

% median normalize
Tr1b=Tr1b./median(Tr1b,2);
disp('median normalization')

%bleaching correction

traces=Tr1b;
for k = 1:NCell
    Tr1b(k, :) = detrend(Tr1b(k,:),'Continuous',false);
 end
disp('bleaching correction')

% Savitzky-Golay filter
 Tr1b = sgolayfilt(Tr1b',3,5)';
disp('sgolayfilter')

[NCell,Nz] = size(Tr1b);

% %refine only for speed<1
WinRest=find(speed<=1);
WinActive=find(speed>1);
% Tr1b = Tr1b(:,WinRest);
[NCell,Nz] = size(Tr1b);
disp (['Ncells= ' num2str(NCell)])
MovT=transpose(1:Nz);  %put real time

%disp (['remove running frame: from ' num2str(length(traces))  ' to '  num2str(Nz)  ])
%% Detect small calcium transients

%f = figure('visible','off');
figure
for i = 1:NCell
    plot(MovT,Tr1b(i,:)+i-1)
    hold on
end

% Detect Calcium Transients using a sliding window
Raster = zeros(NCell,Nz);
Acttmp2 = cell(1,NCell);

for i=1:NCell    
    % th(i)=3*iqr(Tr1b(i,:));
    th(i)=max ([  3*iqr(Tr1b(i,:)),  3*std(Tr1b(i,:)) ,0.2]) ;
    % th(i)=3*std(Tr1b(i,:));
    [amplitude,locs] = findpeaks(Tr1b(i,:),'MinPeakProminence',th(i),'MinPeakDistance',MinPeakDistance);
    % valeurs_identiques = intersect (locs,WinActive);
    % locs_sans_ide=setdiff(locs(:), valeurs_identiques);
    % Acttmp2{i}=locs_sans_ide;%%%%%%%%findchangepts(y,MaxNumChanges=10,Statistic="rms")
    Acttmp2{i}=locs;
    ampli{i}=amplitude;
end

%maybe we should get rid of interneurons=too high activity like 5% highest
% %activity
% for i=1:NCell    
%     length(cell2mat(Acttmp2(1,i)))
% end

%f = figure('visible','off');
for i = 1:NCell
    Raster(i,Acttmp2{i}) = 1;           %Raster = real raster of cell activity
    %plot(MovT,Tr1b(i,:)+i-1)
    hold on
    plot(MovT(Acttmp2{i}),Tr1b(i,Acttmp2{i})+i-1,'.r')
end



% Sum activity over n (synchronous_frames ) consecutive frames

MAct = zeros(1,Nz-synchronous_frames);          %MAct= Sum active cells 
for i=1:Nz-synchronous_frames
    MAct(i) = sum(max(Raster(:,i:i+synchronous_frames),[],2));
end
disp(['Sum transient: ' num2str(sum(MAct))])

    %%%%shuffling to find threshold for number of cell for sce detection
    NShfl=100;
    for n=1:NShfl
        for c=1:NCell
            k = randi(Nz-length(WinActive));
            Rastersh(c,:)= circshift(Raster(c,:),k,2);
        end
        for i=1:Nz-synchronous_frames
            MActsh(i) = sum(max(Rastersh(:,i:i+synchronous_frames),[],2));
        end
    
        Sumactsh(:,n)=MActsh;
    end
        percentile = 99; % Calculate the 5% highest point or 99
        sce_n_cells_threshold = prctile(Sumactsh, percentile,"all");
        disp(['sce_n_cells_threshold: ' num2str(sce_n_cells_threshold)])
        %sce_n_cells_threshold = 5;
        %sce_n_cells_threshold =median(Sumactsh,'all');
        %sce_n_cells_threshold =3*iqr(Sumactsh,'all');    
    
    % Select synchronies (RACE)         % TRace=localisation SCE 
    %MinPeakDistancesce=5;
    [pks,TRace] = findpeaks(MAct,'MinPeakHeight',sce_n_cells_threshold,'MinPeakDistance',MinPeakDistancesce);
    % sumpeaks=sum(peaks)
    NRace = length(TRace);
    disp(['nSCE: '  num2str(NRace)])

    % Create RasterPlots%%%%%%very weird here increase RACE from n-1 :n+2 
    Race = zeros(NCell,NRace);      %Race=cells that participate in SCE 
    RasterRace = zeros(NCell,Nz);
    for i = 1:NRace
        Race(:,i) = max(Raster(:,TRace(i)-1:TRace(i)+2),[],2);    %maybe this window can be optimized ???
        RasterRace(Race(:,i)==1,TRace(i)) = 1;                          %Raster with only activity in SCE
    end
    
    % % Display race
    % for i = 1:length(TRace)
    %     line(MovT(TRace(i))*[1 1],[0 NCell+1],'Color','g');
    % end

     %exportgraphics(gcf,[namefull 'activityall.tif'],'Resolution',1200)
    % exportgraphics(gcf,[namefull 'activityall.pdf'],'ContentType','vector')

%display raster
    % % Display race

    %©
    % figure
    % for i = 1:length(TRace)
    % %     line(MovT(TRace(i))*[1 1],[0 NCell+1],'Color','g');
    % % end
    % hold on 
    % imagesc(Raster)
    % colormap jet
    % axis image
    % xlabel('sorted SCE #')     %was RACE
    % ylabel('sorted SCE #')  
    %plot MAct under
    %plot speed under

    % exportgraphics(gcf,[namefull 'raster.png'],'Resolution',300)


    %return
    % break
    %% Save
      save([namefull,'Acttmp2.mat'],'Acttmp2')      %cell activity (time of peaks) 
      save([namefull,'Raster.mat'],'Race')      %Raster = real raster of cell activity
      save([namefull,'Race.mat'],'Race')        %Race=cells X SCE 
      save([namefull,'RasterRace.mat'],'Race')      %Raster with only activity in SCE over time
      save([namefull,'TRace.mat'],'TRace')      %time of significant SCE
    
    %% Clustering
    
    [NCell,NRace] = size(Race);
    %[IDX2,sCl,M,S] = kmeansopt(Race,10,'var');
    %M=covariance matrix  S= Silhouette but of all clustering X10 times
    %sCl=median des silhouettes des clusters du meilleur clustering      IDX2=
    [IDX2,sCl,M,S] = kmeansopt(Race,200,'var');%increase to 50,  100 ?    
    % M = CovarM(Race);
    % IDX2 = kmedoids(M,NCl);
    NCl = max(IDX2);
    [~,x2] = sort(IDX2);
    MSort = M(x2,x2);
    disp(['nClusters: '  num2str(NCl)])

    %Race clusters
    R = cell(0);
    CellScore = zeros(NCell,NCl);
    CellScoreN = zeros(NCell,NCl);
    for i = 1:NCl
        R{i} = find(IDX2==i);
        CellScore(:,i) = sum(Race(:,R{i}),2);
        CellScoreN(:,i) = CellScore(:,i)/length(R{i});
    end
    %Assign cells to cluster with which it most likely spikes
    [~,CellCl] = max(CellScoreN,[],2);
    %Remove cells with less than 2 spikes in a given cluster   ###should
    %export this raw CellCl  
    CellCl(max(CellScore,[],2)<2) = 0;
    [X1,x1] = sort(CellCl);
    
%write list of cells in each cluster

    assemblyraw= cell(0);
    k = 0;
    for i = 1:NCl
        k = k+1;
        assemblyraw{k} = transpose(find(CellCl==i));
    end

    save([namefull,'assemblyraw.mat'],'assemblyraw')  
    
    % figure
     f = figure('visible','off');
    subplot(1,2,1)
    imagesc(MSort)
    colormap jet
    axis image
    xlabel('sorted SCE #')     %was RACE
    ylabel('sorted SCE #')     %was RACE
    
    subplot(1,2,2)
    imagesc(Race(x1,x2),[-1 1.2])
    axis image
    xlabel('sorted SCE #')     %was RACE
    ylabel('sorted Cell #')
    
    exportgraphics(gcf,[namefull 'clusters.png'],'Resolution',300)
    
    %% Save Clusters
     %save([namefull,'Clusters.mat'],'IDX2')   % liste de tous les SCE et dans quel cluster ils sont
        
    %% Remove cluster non-statistically significant: basically do some sce permutation 
    % and some clustering with the same number of cluster and if the best silhouette
    % of one cluster is higher than real data then remove this cluster.
    % since they are sorted remove the worst one 
    
   
    kmeans_surrogate=100; %or more...
    sClrnd = zeros(1,kmeans_surrogate);
    parfor i = 1:kmeans_surrogate  
        sClrnd(i) = kmeansoptrnd(Race,10,NCl); 
    end
    %NClOK = sum(sCl>max(sClrnd));
    NClOK =sum(sCl>prctile(sClrnd,95));
    sClOK = sCl(1:NClOK)';
    disp(['nClustersOK: ' num2str(NClOK)])

    %new list of cells %JC

    assemblystat= cell(0);
    k = 0;
    for i = 1:NClOK
        k = k+1;
        assemblystat{k} = transpose(find(CellCl==i));
    end

    save([namefull,'assemblystat.mat'],'assemblystat')
    save([namefull,'NClustersOK.mat'],'NClOK')
    
    RaceOK = Race(:,IDX2<=NClOK);
    NRaceOK = size(RaceOK,2);
    disp(['nSCEOK: ' num2str(NRaceOK)])    

    % we should plot Race OK before othogonalisation
    %%%%
    %
    %pause (1)
    if NClOK>1
         RACE_Ortho
    else
        NCl=NClOK;
         assemblyortho= cell(0);
         assemblystat= cell(0);
       % assemblyraw=[]
    end

    ncluster(sce_n_cells_threshold)=NCl;
    test=strcat ('percentile sce threshold=', num2str(percentile), ' sce_n_cells_threshold=', num2str( sce_n_cells_threshold), ' ncluster= ' , num2str(NCl) , (  ' NRace= ') , num2str(NRace), (  ' NRaceOK= ') ,  num2str(NRaceOK) ) ;
    %disp(strcat ('synchronous_frames= ', num2str(synchronous_frames),  ' minpeak distactivity', num2str(MinPeakDistance ),  ' ncluster= ' , num2str(NCl) , (  ' NRace= ') , num2str(NRace), (  ' NRaceOK= ') ,  num2str(NRaceOK) )) 
    % disp(daytime)
    % disp(test)

    %open file identifier
     fid=fopen([namefull,'settings.txt'],'w');
     fprintf(fid, test);
     fclose(fid);
    brainbowassemblies2023_06_09
    %write all data in excel
    exportdata
    toc

%end
