clear
close all

% % delete(gcp('nocreate'))
% parpool ('processes',4)

tic

%% Load settings

% filename='/Users/platel/Desktop/exp/Sofia/mc_tiff_downs_6.5_10_um.tif';     %to fill up 
filename='E:\backup_various\sofia\best_data_mc_downsampled\P2_P3_P4\emx1\ani80_2023-02-23\tiff\downs_20_80_b_p3.tif';

%filename='E:\backup_various\sofia\best_data_mc_downsampled\P0_P1\emx1\ani4mc_tiff_downs_6.5_10_um.tif';

PathSave='E:\backup_various\sofia\analysisJC\';      %to fill up 
[pathold,name,ext] = fileparts(filename);

daytime = datestr(now,'yy_mm_dd_HH_MM_SS');
namefull=[PathSave name daytime  '/'];
mkdir (namefull)    % make folder for saving analysis

MinPeakDistancesce = 3 ;
MinPeakDistance = 3 ;
kmean_iter=100;
sampling_rate=15;           % to fill up
synchronous_frames=round(0.2*sampling_rate,0); %200ms *sampling rate
size_pixel=20;%in um
% sce_n_cells_threshold = 10000/size_pixel^2; %10000=100um^2 basically size threshold of 100um 100 
% sce_n_cells_threshold = 5000/size_pixel^2; %10000=100um^2 basically size threshold of 100um =50
sce_n_cells_threshold = 150;

%% Load  data

%F=double(F(iscell(:,1)>0,:));

% name='/Users/platel/Desktop/exp/Sofia/mc_tiff_downs_6.5_10_um.tif';
% % tifInfo = imfinfo('/Users/platel/Desktop/exp/Sofia/ani55/mc_concat_reca_940_92um-001-_grid_13_factor.tif');
tifInfo = imfinfo(filename);
% 
numPages = numel(tifInfo);

sizex=tifInfo(1).Width;
sizey=tifInfo(1).Height;
% %tifStack = cell(1, numPages);
F=zeros(sizex*sizey, numPages);
% %55*55
% 
for i = 1:numPages
    F(:,i) = reshape(imread(filename, i),sizex*sizey,1);
end

% %% preprocessing

% % median normalize
[NCell,Nz] = size(F);
% sce_n_cells_threshold =50;%5/100*NCell; %to adjust
% F=F./median(F,2);
F=normalize(F,2,"center","median");
% ws = warning('off','all');
% for i=1:NCell
%     p0=polyfit(1:Nz,F(i,:),3);
%     F(i,:)=F(i,:)./polyval(p0,1:Nz);
% end
% warning(ws)%% preprocessing
% 
% disp('median normalization')
% 
% %bleaching correction
% 
% % traces=Tr1b;
% % for k = 1:NCell
% %     Tr1b(k, :) = detrend(Tr1b(k,:),'Continuous',false);
% %  end
% % disp('bleaching correction')
% 
% % Savitzky-Golay filter
F = sgolayfilt(F',3,5)';
% disp('sgolayfilter')



% %refine only for speed<1
% WinRest=find(speed<=1);
WinActive=[];%find(speed>1);
% Tr1b = Tr1b(:,WinRest);

disp (['Ncells= ' num2str(NCell)])
MovT=transpose(1:Nz);  %put real time

%disp (['remove running frame: from ' num2str(length(traces))  ' to '  num2str(Nz)  ])
%% Detect small calcium transients

% f = figure('visible','off');
% figure
% for i = 1:NCell
%     plot(MovT,F(i,:)+i-1)
%     hold on
% end

% Detect Calcium Transients using a sliding window
Raster = zeros(NCell,Nz);
Acttmp2 = cell(1,NCell);
ampli = cell(1,NCell);
minithreshold=0.1;

for i=1:NCell    
    th=max ([  3*iqr(F(i,:)),  3*std(F(i,:)) ,minithreshold]) ;
    [amplitude,locs] = findpeaks(F(i,:),'MinPeakProminence',th,'MinPeakDistance',MinPeakDistance);
    Acttmp2{i}=locs;%%%%%%%%findchangepts(y,MaxNumChanges=10,Statistic="rms")
end


%f = figure('visible','off');
for i = 1:NCell
    Raster(i,Acttmp2{i}) = 1;           %Raster = real raster of cell activity
    %plot(MovT,Tr1b(i,:)+i-1)
    % hold on
    % plot(MovT(Acttmp2{i}),Tr1b(i,Acttmp2{i})+i-1,'.r')
end

% Sum activity over n (synchronous_frames ) consecutive frames

MAct = zeros(1,Nz-synchronous_frames);          %MAct= Sum active cells 
for i=1:Nz-synchronous_frames
    MAct(i) = sum(max(Raster(:,i:i+synchronous_frames),[],2));
end
disp(['Sum transient: ' num2str(sum(MAct))])

    %%%%shuffling to find threshold for number of cell for sce detection
% MActsh = zeros(1,Nz-synchronous_frames);   
% Rastersh=zeros(NCell,Nz);   
% NShfl=100;
% Sumactsh=zeros(Nz-synchronous_frames,NShfl);   
% for n=1:NShfl
% 
%     for c=1:NCell
%         k = randi(Nz-length(WinActive));
%         Rastersh(c,:)= circshift(Raster(c,:),k,2);
%     end
% 
%     for i=1:Nz-synchronous_frames   %need to use WinRest???
%         MActsh(i) = sum(max(Rastersh(:,i:i+synchronous_frames),[],2));
%     end
% 
%     Sumactsh(:,n)=MActsh;
% 
% end
% 
% percentile = 95; % Calculate the 5% highest point or 99
%  sce_n_cells_threshold = prctile(Sumactsh, percentile,"all");
% % sce_n_cells_threshold =10;
% 
% disp(['sce_n_cells_threshold: ' num2str(sce_n_cells_threshold)])
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

% figure
% for i = 1:length(TRace)
%      line(MovT(TRace(i))*[1 1],[0 NCell+1],'Color','g');
% end
% hold on 
% imagesc(Raster)
% colormap jet
% axis image
% xlabel('sorted SCE #')     %was RACE
% ylabel('sorted SCE #')  
% plot MAct under
% plot speed under

% exportgraphics(gcf,[namefull 'raster.png'],'Resolution',300)


%return
% break
%% Save
  % save([namefull,'Acttmp2.mat'],'Acttmp2')      %cell activity (time of peaks) 
  % save([namefull,'Raster.mat'],'Race')      %Raster = real raster of cell activity
  % save([namefull,'Race.mat'],'Race')        %Race=cells X SCE 
  % save([namefull,'RasterRace.mat'],'Race')      %Raster with only activity in SCE over time
  % save([namefull,'TRace.mat'],'TRace')      %time of significant SCE

%% Clustering

[NCell,NRace] = size(Race);
%[IDX2,sCl,M,S] = kmeansopt(Race,10,'var');
%M=covariance matrix  S= Silhouette but of all clustering X10 times
%sCl=median des silhouettes des clusters du meilleur clustering      IDX2=
% [IDX2,sCl,M,S] = kmeansopt(Race,kmean_iter,'var');%increase to 50,  100 ?  

[IDX2,sCl,M,S] = kmeansopt(Race,100,'var');% find best silhouette from 2:18

NCl = max(IDX2);
[IDX2,sCl,M,S] = kmeansoptbarrel(Race,kmean_iter,'var',NCl);% 

% M = CovarM(Race);
% IDX2 = kmedoids(M,NCl);
NCl = max(IDX2);
[~,x2] = sort(IDX2);
MSort = M(x2,x2);
disp(['nClusters: '  num2str(NCl)])

% save([namefull,'assemblyraw.mat'],'assemblyraw')  

% figure
%  f = figure('visible','off');
%  figure
% subplot(1,2,1)
% imagesc(MSort)
% colormap jet
% axis image
% xlabel('sorted SCE #')     %was RACE
% ylabel('sorted SCE #')     %was RACE
% 
% subplot(1,2,2)
% imagesc(Race(x1,x2),[-1 1.2])
% % axis image
% xlabel('sorted SCE #')     %was RACE
% ylabel('sorted Cell #')
% axis tight
% 
% exportgraphics(gcf,[namefull 'clusters.png'],'Resolution',300)

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
% [~,CellCl] = max(CellScoreN,[],2);  removed 2024-03-21 to prevent
% orthogonalisation
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


%new list of cells %JC

assemblystat= cell(0);
k = 0;
for i = 1:NClOK
    k = k+1;
    assemblystat{k} = transpose(find(CellCl==i));
end

% save([namefull,'assemblystat.mat'],'assemblystat')
% save([namefull,'NClustersOK.mat'],'NClOK')

RaceOK = Race(:,IDX2<=NClOK);
NRaceOK = size(RaceOK,2);
disp(['nSCEOK: ' num2str(NRaceOK)])    

f = figure('visible','off');
 % figure
subplot(1,2,1)
imagesc(MSort)
colormap jet
axis image
xlabel('sorted SCE #')     %was RACE
ylabel('sorted SCE #')     %was RACE

subplot(1,2,2)
imagesc(Race(x1,x2),[-1 1.2])
% axis image
xlabel('sorted SCE #')     %was RACE
ylabel('sorted Cell #')
axis tight
namegraph='RACE.png';
exportgraphics(gcf,[namefull namegraph],'Resolution',300);
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
% test=strcat ('percentile sce threshold=', num2str(percentile), ' sce_n_cells_threshold=', num2str( sce_n_cells_threshold), ' ncluster= ' , num2str(NCl) , (  ' NRace= ') , num2str(NRace), (  ' NRaceOK= ') ,  num2str(NRaceOK) ) ;
% %disp(strcat ('synchronous_frames= ', num2str(synchronous_frames),  ' minpeak distactivity', num2str(MinPeakDistance ),  ' ncluster= ' , num2str(NCl) , (  ' NRace= ') , num2str(NRace), (  ' NRaceOK= ') ,  num2str(NRaceOK) )) 
% % disp(daytime)
% % disp(test)
% 
% %open file identifier
% 
%  fid=fopen([namefull,'settings.txt'],'w');
%  fprintf(fid, test);
%  fclose(fid);

%brainbowassemblies2023_06_09
%write all data in excel

save([namefull,'results.mat'])  
plot_grid_20240321

% exportdata

toc
