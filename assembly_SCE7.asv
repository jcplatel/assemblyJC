%% Load settings
MinPeakDistancesce=5 ;% 5 default
MinPeakDistance=3;% 3 default
% sampling_rate=10;
synchronous_frames=round(0.2*sampling_rate,0); %200ms *sampling rate
% synchronous_frames=2; %2 default
% kmean_iter=100;
% kmeans_surrogate=100;
percentile=NaN;

%% Load current data

% Tr1b=double(F(iscell(:,1)>0,:));
Tr1b=double(F);
% clear F
% Tr1b=F(colorcell<7,:);
% if mean(Tr1b,"all")>10000
%     disp('pb intensity DF')
%     return
% end
speedsm =smoothdata(speed,'gaussian',50);
[NCell,Nz] = size(Tr1b);

% median normalize

% disp('median normalization')

%bleaching correction
Tr1b = sgolayfilt(Tr1b',3,9)';
% ws = warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
% for i=1:NCell
%     p0=polyfit(1:Nz,Tr1b(i,:),3);
%     Tr1b(i,:)=Tr1b(i,:)./polyval(p0,1:Nz);
% end
% warning(ws)%% preprocessing

window_size = 2000; % Largeur de la fenêtre en points temporels
percentile_value=5;
num_blocks = ceil(Nz / window_size);
for n=1:NCell
    trace=Tr1b(n,:);
    F0 = nan(Nz, 1);
    % Calcul du percentile bloc par bloc
    for i = 1:num_blocks
        start_idx = (i-1) * window_size + 1;
        end_idx = min(i * window_size, Nz);
        F0(start_idx:end_idx) = prctile(trace(start_idx:end_idx), percentile_value);
    end
    F0 = movmedian(F0, window_size, 'omitnan');
    F0 = smoothdata (F0,1,"gaussian",window_size/2);
    Tr1b(n,:)=(trace-F0')./F0';
end

% %refine only for speed<2
WinRest=find(speedsm<=2);
WinActive=find(speedsm>2);

% Tr1b = Tr1b(:,WinRest);
[NCell,Nz] = size(Tr1b);
disp (['Ncells= ' num2str(NCell)])
MovT=transpose(1:Nz);  %put real time

% Detect Calcium Transients using a sliding window
Raster = zeros(NCell,Nz);
Acttmp2 = cell(1,NCell);
ampli = cell(1,NCell);
minithreshold=0.1; 

for i=1:NCell    
    
    th(i)=2.576*std(Tr1b(i,WinRest));
    % th(i)=max ([2*iqr(Tr1b(i,WinRest)) ,2.576*std(Tr1b(i,WinRest))]) ;
    % [amplitude,locs] = findpeaks(Tr1b(i,:),'MinPeakProminence',th(i),'MinPeakDistance',MinPeakDistance);
    [amplitude,locs] = findpeaks(Tr1b(i,:),'MinPeakProminence',th(i) ,'MinPeakDistance',MinPeakDistance);%2.576= 99%, 3.291=99.9
    valeurs_identiques = intersect (locs,WinActive);
    % ampALL{i}=amplitude;
    [locs_sans_ide , idx ]=setdiff(locs(:), valeurs_identiques);
    ampli_sans_ide=amplitude(idx);
    % ampALL{i}=ampli_sans_ide;
    % Acttmp2{i}=locs_sans_ide(ampli_sans_ide>0.05 & ampli_sans_ide<1);
    Acttmp2{i}=locs_sans_ide;%%%%%%%%findchangepts(y,MaxNumChanges=10,Statistic="rms")
    %Acttmp2{i}=locs;
    % ampli{i}=amplitude;

end

% %f = figure('visible','off');
for i = 1:NCell
    if Acttmp2{i}>0
        Raster(i,Acttmp2{i}) = 1;           %Raster = real raster of cell activity
    end
    % plot(MovT,Tr1b(i,:)+i-1)
    % hold on
    % plot(MovT(Acttmp2{i}),Tr1b(i,Acttmp2{i})+i-1,'.r')
end

% Sum activity over n (synchronous_frames ) consecutive frames

MAct = zeros(1,Nz-synchronous_frames);          %MAct= Sum active cells 
for i=1:Nz-synchronous_frames
    MAct(i) = sum(max(Raster(:,i:i+synchronous_frames),[],2));
end
% sum(MAct>10)
% MAct (MAct>100)=0;
 disp(['Sum transient: ' num2str(sum(MAct))])

%%%%shuffling to find threshold for number of cell for sce detection
% MActsh = zeros(1,Nz-synchronous_frames);   
% Rastersh=zeros(NCell,Nz);   
% NShfl=100;
% Sumactsh=zeros(Nz-synchronous_frames,NShfl);   
% 
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
% end
% 
% percentile = 95; % Calculate the 5% highest point or 99
% sce_n_cells_threshold = prctile(Sumactsh, percentile,"all");
% % 
% % disp(['sce_n_cells_threshold: ' num2str(sce_n_cells_threshold)])
sce_n_cells_threshold = 10;
% if sce_n_cells_threshold>20 ; sce_n_cells_threshold=20; end
% if sce_n_cells_threshold<5 ; sce_n_cells_threshold=5; end
% Select synchronies (RACE)         % TRace=localisation SCE 

[pks,TRace] = findpeaks(MAct,'MinPeakHeight',sce_n_cells_threshold,'MinPeakDistance',MinPeakDistancesce);
% idx=find(pks<70);
%TRace=TRace(idx);
NRace = length(TRace);

disp(['nSCE: '  num2str(NRace)])

% Create RasterPlots%%%%%%very weird here increase RACE from n-1 :n+2 
Race = zeros(NCell,NRace);      %Race=cells that participate in SCE 
RasterRace = zeros(NCell,Nz);
for i = 1:NRace
    Race(:,i) = max(Raster(:,TRace(i)-1:TRace(i)+2),[],2);    %maybe this window can be optimized ???
    RasterRace(Race(:,i)==1,TRace(i)) = 1;                          %Raster with only activity in SCE time: all frames
end

%% Clustering

[NCell,NRace] = size(Race);
% NClini=10;
if ~exist ('NClini','var')
    NClini=0;
end
%Race=gpuArray(Race);
[IDX2,sCl,M,S,NClini] = kmeansopttest(Race,kmean_iter,'var',NClini);%increase to 50,  100 ?    
% M = CovarM(Race);
% IDX2 = kmedoids(M,NCl);
NCl = max(IDX2);
[~,x2] = sort(IDX2);
MSort = M(x2,x2);
disp(['nClusters: '  num2str(NCl)])


%% Save Clusters
 %save([namefull,'Clusters.mat'],'IDX2')   % liste de tous les SCE et dans quel cluster ils sont
    
%% Remove cluster non-statistically significant: basically do some sce permutation 
% and some clustering with the same number of cluster and if the best silhouette
% of one cluster is higher than real data then remove this cluster.
% since they are sorted remove the worst one 

% kmeans_surrogate=500; %or more...
sClrnd = zeros(1,kmeans_surrogate);
cRace = parallel.pool.Constant(Race);
parfor i = 1:kmeans_surrogate  
    sClrnd(i) = kmeansoptrnd(cRace.Value,10,NCl); 
end


%NClOK = sum(sCl>max(sClrnd));
NClOK =sum(sCl>prctile(sClrnd,95)); %use max, use 99%  ?
sClOK = sCl(1:NClOK)';
disp(['nClustersOK: ' num2str(NClOK)])


%%%new add
%Race clusters
R = cell(0);        % which sce for each cluster
CellScore = zeros(NCell,NClOK);  % number of time a cell participate in a sce of a given cluster
CellScoreN = zeros(NCell,NClOK);
for i = 1:NClOK
    R{i} = find(IDX2==i);
    CellScore(:,i) = sum(Race(:,R{i}),2);
    CellScoreN(:,i) = CellScore(:,i)/length(R{i});
end
%Assign cells to cluster with which it most likely spikes in term of
%percentage

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

%%%%%%%

RaceOK = Race(:,IDX2<=NClOK);
NRaceOK = size(RaceOK,2);
% disp(['nSCEOK: ' num2str(NRaceOK)])    

figure('visible','off');
% figure
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

exportgraphics(gcf,strcat(namefull ,'clusters.png'),'Resolution',300)
close gcf

if NClOK>1
     RACE_Orthojc
     %%call rasterplot here
else
    NCl=NClOK;
     assemblyortho= cell(0);
     assemblystat= cell(0);
   % assemblyraw=[]
end

if NCl ==0
    NClOK=0;
     assemblyortho= cell(0);
     assemblystat= cell(0);
end
disp(['NCl: ' num2str(NCl)])  
% ncluster(sce_n_cells_threshold)=NCl;

%write all data in excel

% exportdata
% save([namefull,'results.mat'])  

