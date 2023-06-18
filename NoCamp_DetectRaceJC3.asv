%default values
%MinPeakDistancesce=5 frames= 500ms
% synchronous_frames=2
%MinPeakDistance=5 frames= 500ms
%sce_n_cells_threshold=5 or 5% ???   seem best results between 10 to 15
%cells  should we simulate ???

clear
close all

%for sce_n_cells_threshold=10:20%[5,10,15,20]%5:30%
openingnwb


tic
%% Load settings
MinPeakDistancesce=5 ;
synchronous_frames=2;
MinPeakDistance=5;

%imaging_sampling_rate=15.7;
synchronous_frames=round(0.2*sampling_rate,0); %200ms *sampling rate

% PathSave='/Users/platel/Desktop/exp/analysis/';
PathSave='/Users/platel/Desktop/exp/analysis/';

 daytime = datestr(now,'yy_mm_dd_HH_MM_SS');
 namefull=[PathSave daytime name '/'];
 mkdir (namefull)    % make folder for saving analysis

%% Load current data
% 
% plane0=double(plane0(iscell0(:,1)>0,:));
% plane1=double(plane1(iscell1(:,1)>0,:));
% 
Tr1b=F;

% % Tr1b=plane1;

%% preprocessing

%smooth speed
speed =smoothdata(speed,'gaussian',50);

% median normalize
Tr1b=Tr1b./median(Tr1b,2);

[NCell,Nz] = size(Tr1b);

%bleaching correction

% traces=Tr1b;
% 
% for k = 1:NCell
%     p0 = polyfit(1:Nz, traces(k, :), 1);
%     traces(k, :) = traces(k, :) ./ polyval(p0, 1:Nz);
% end

% Savitzky-Golay filter
 Tr1b = sgolayfilt(Tr1b',2,3)';

% 
% %refine only for speed<1
WinRest=find(speed<1);
Tr1b = Tr1b(:,WinRest);
[NCell,Nz] = size(Tr1b);

MovT=transpose(1:Nz);  %put real time

%% Detect small calcium transients

% figure
% for i = 1:NCell
%     plot(MovT,Tr1b(i,:)+i-1)
%     hold on
% end

% Detect Calcium Transients using a sliding window
Raster = zeros(NCell,Nz);
Acttmp2 = cell(1,NCell);

for i=1:NCell    
    th=3*iqr(Tr1b(i,:));
    %th(i)=3*std(Tr1b(i,:));
    [~,locs] = findpeaks(Tr1b(i,:),'MinPeakProminence',th,'MinPeakDistance',MinPeakDistance);
    Acttmp2{i}=locs;%%%%%%%%findchangepts(y,MaxNumChanges=10,Statistic="rms")
end

%figure
for i = 1:NCell
    Raster(i,Acttmp2{i}) = 1;           %Raster = real raster of cell activity
    % plot(MovT,Tr1b(i,:)+i-1)
    % hold on
    % plot(MovT(Acttmp2{i}),Tr1b(i,Acttmp2{i})+i-1,'.r')
end

% Sum activity over n (synchronous_frames ) consecutive frames

MAct = zeros(1,Nz-synchronous_frames);          %MAct= Sum active cells 
for i=1:Nz-synchronous_frames
    MAct(i) = sum(max(Raster(:,i:i+synchronous_frames),[],2));
end

%for z=1:2
    %%%%shuffling to find threshold for number of cell for sce detection
    NShfl=100;
    for n=1:NShfl
        for c=1:NCell
            k = randi(Nz);
            Rastersh(c,:)= circshift(Raster(c,:),k,2);
        end
        for i=1:Nz-synchronous_frames
            MActsh(i) = sum(max(Rastersh(:,i:i+synchronous_frames),[],2));
        end
    
        Sumactsh(:,n)=MActsh;
    end
        percentile = 95; % Calculate the 5% highest point or 99
        sce_n_cells_threshold = prctile(Sumactsh, percentile,"all");
        %sce_n_cells_threshold = 15;
        %sce_n_cells_threshold =median(Sumactsh,'all');
        %sce_n_cells_threshold =3*iqr(Sumactsh,'all');    
    
    % Select synchronies (RACE)         % TRace=localisation SCE 
    %MinPeakDistancesce=5;
    [pks,TRace] = findpeaks(MAct,'MinPeakHeight',sce_n_cells_threshold,'MinPeakDistance',MinPeakDistancesce);
    % sumpeaks=sum(peaks)
    NRace = length(TRace);
    
    % Create RasterPlots%%%%%%very weird here increase RACE from n-1 :n+2 
    Race = zeros(NCell,NRace);      %Race=cells that participate in SCE 
    RasterRace = zeros(NCell,Nz);
    for i = 1:NRace
        Race(:,i) = max(Raster(:,TRace(i)-1:TRace(i)+2),[],2);    %maybe this window can be optimized ???
        RasterRace(Race(:,i)==1,TRace(i)) = 1;                          %Raster with only acticity in SCE
    end
    
    % % Display race
    % for i = 1:length(TRace)
    %     line(MovT(TRace(i))*[1 1],[0 NCell+1],'Color','g');
    % end
    %toc
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
    [IDX2,sCl,M,S] = kmeansopt(Race,100,'var');%increase to 50,  100 ?    
    % M = CovarM(Race);
    % IDX2 = kmedoids(M,NCl);
    NCl = max(IDX2);
    [~,x2] = sort(IDX2);
    MSort = M(x2,x2);
    
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
    %Remove cells with less than 2 spikes in a given cluster
    CellCl(max(CellScore,[],2)<2) = 0;
    [X1,x1] = sort(CellCl);
    
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
        
    %% Remove cluster non-statistically significant
    
    sClrnd = zeros(1,100);
    for i = 1:20
        sClrnd(i) = kmeansoptrnd(Race,10,NCl);  %increase to 50, 100 ???
    end
    NClOK = sum(sCl>max(sClrnd));
    sClOK = sCl(1:NClOK)';
    
    save([namefull,'NClustersOK.mat'],'NClOK')
    
    RaceOK = Race(:,IDX2<=NClOK);
    NRaceOK = size(RaceOK,2);
    
    pause (1)
    if NClOK>1
         RACE_Ortho
    else
        NCl=NClOK;
    end
    toc
    ncluster(sce_n_cells_threshold)=NCl;
    test=strcat ('percentile sce threshold=', num2str(percentile), ' sce_n_cells_threshold=', num2str( sce_n_cells_threshold), ' ncluster= ' , num2str(NCl) , (  ' NRace= ') , num2str(NRace), (  ' NRaceOK= ') ,  num2str(NRaceOK) ) ;
    %disp(strcat ('synchronous_frames= ', num2str(synchronous_frames),  ' minpeak distactivity', num2str(MinPeakDistance ),  ' ncluster= ' , num2str(NCl) , (  ' NRace= ') , num2str(NRace), (  ' NRaceOK= ') ,  num2str(NRaceOK) )) 
    disp(daytime)
    disp(test)

    %open file identifier
     fid=fopen([namefull,'settings.txt'],'w');
     fprintf(fid, test);
     fclose(fid);

%end
%end
