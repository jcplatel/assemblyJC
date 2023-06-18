clear
close all

openingnwb
%% Load current data

% load('WinRest')  %period of rest
% load('Tr1b')    %fluo
% load('Cells')       %movie
% load('MovT')        % times of ????   images ???
% load('Speed')      % instantaneous speed
imaging_sampling_rate=8;
synchronous_frames=round(0.2*imaging_sampling_rate,0); %200ms *sampling rate
PathSave='/Users/platel/Desktop/exp/analysis/';

 daytime = datestr(now,'yy_mm_dd_HH_MM_SS');
 namefull=[PathSave daytime name '/'];
 mkdir (namefull)    % make folder for saving analysis
disp(['make new folder ' namefull])


 %Tr1b=double(plane0(iscell(:,1)>0,:));
% plane1=double(plane1(iscell1(:,1)>0,:));

% Tr1b=vertcat(plane0, plane1);

% MovT=transpose(1:12378);  %put real time
%  speed = readNPY('/Users/platel/Desktop/exp/speed.npy');
%  sce=readNPY('/Users/platel/Desktop/exp/sce.npy');
% allcells=readNPY('/Users/platel/Desktop/exp/aurelie/suite2p_444175_221125_plane1/F.npy');
% iscell=readNPY('/Users/platel/Desktop/exp/aurelie/suite2p_444175_221125_plane1/iscell.npy');
% Tr1b=double(allcells(iscell(:,1)>0,:));
Tr1b=Tr1b./median(Tr1b,2);

%bleaching correction

%traces=Tr1b;
for k = 1:NCell
    Tr1b(k, :) = detrend(Tr1b(k,:),'Continuous',false);
end

%Tr1b=readNPY('/Users/platel/Desktop/exp/sce.npy');
WinRest=find(speed<1);

%% Detect small calcium transients
[NCell,Nz] = size(Tr1b);

sce_n_cells_threshold=round(0.05*NCell);


% Savitzky-Golay filter
Tr1b = sgolayfilt(Tr1b',3,5)';

% figure
% for i = 1:NCell
%     plot(MovT,Tr1b(i,:)+i-1)
%     hold on
% end

% Detect Calcium Transients using a sliding window
%TrRest = Tr1b(:,WinRest);
TrRest = Tr1b;
Raster = zeros(NCell,Nz);
%WinSize = 40;
WinSize = round(4 * imaging_sampling_rate);

parfor i=1:NCell    
    Acttmp = zeros(1,Nz);
    Sigtmp = zeros(1,Nz);
    Trtmp = Tr1b(i,:);
    %Remove points with high baseline
    ThBurst = median(Trtmp) + iqr(Trtmp)/2;
    for j = WinSize+1:Nz-WinSize
        if speed(j,1)<1
            Wintmp = j-WinSize:j+WinSize;
            Mediantmp = median(Trtmp(Wintmp));
            %Not active in 10 last frames and not within burst activity
            if sum(Acttmp(j-10:j-1)) == 0 && Mediantmp < ThBurst 
                Acttmp(j) = Trtmp(j) - Mediantmp > 3*iqr(Trtmp(Wintmp));
                Sigtmp(j) = (Trtmp(j) - Mediantmp) / iqr(Trtmp(Wintmp));
            end
        end
    end
    Acttmp2{i} = find(Acttmp);
    Sigtmp2{i} = Sigtmp(Acttmp2{i});

end
for i = 1:NCell
    Raster(i,Acttmp2{i}) = 1;
    %plot(MovT(Acttmp2{i}),Tr1b(i,Acttmp2{i})+i-1,'.r')
end

% Sum activity over n synchronous_frames consecutive frames

MAct = zeros(1,Nz-synchronous_frames);
for i=1:Nz-synchronous_frames
    MAct(i) = sum(max(Raster(:,i:i+synchronous_frames),[],2));
end

% Select synchronies (RACE)
%Th = 5;
[pks,TRace] = findpeaks(MAct,'MinPeakHeight',sce_n_cells_threshold,'MinPeakDistance',5);

NRace = length(TRace);

% Create RasterPlots
Race = zeros(NCell,NRace);
RasterRace = zeros(NCell,Nz);
for i = 1:NRace
    Race(:,i) = max(Raster(:,TRace(i)-1:TRace(i)+2),[],2);
    RasterRace(Race(:,i)==1,TRace(i)) = 1;
end

% Display race
% for i = 1:length(TRace)
%    % line(MovT(TRace(i))*[1 1],[0 NCell+1],'Color','g');
% end
%toc
%return
% break
% %% Save
% save([PathSave,'Acttmp2.mat'],'Acttmp2')
% save([PathSave,'Race.mat'],'Race')
save([namefull,'TRace.mat'],'TRace')   
% 
% 

%% Clustering
[NCell,NRace] = size(Race);
[IDX2,sCl,M,S] = kmeansopt(Race,200,'var');
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

figure
subplot(1,2,1)
imagesc(MSort)
colormap jet
axis image
xlabel('RACE #')
ylabel('RACE #')

subplot(1,2,2)
imagesc(Race(x1,x2),[-1 1.2])
axis image
xlabel('RACE #')
ylabel('Cell #')
exportgraphics(gcf,[name 'clusters.png'],'Resolution',300)

%% Save Clusters
% save([PathSave,'Clusters.mat'],'IDX2')


%% Remove cluster non-statistically significant

sClrnd = zeros(1,20);
for i = 1:100
    sClrnd(i) = kmeansoptrnd(Race,10,NCl);  %was at 10
end
%NClOK = sum(sCl>max(sClrnd));
NClOK =sum(sCl>prctile(sClrnd,95));
sClOK = sCl(1:NClOK)';
% 
% save([PathSave,'NClustersOK.mat'],'NClOK')

RaceOK = Race(:,IDX2<=NClOK);
NRaceOK = size(RaceOK,2);
if NClOK>1
    RACE_Ortho
end
toc
pause(1)
ncluster(sce_n_cells_threshold)=NCl
nracemax(sce_n_cells_threshold)=NRaceOK
test=strcat ('MinPeakDistancesce= ', num2str(5), ' sce_n_cells_threshold=', num2str( sce_n_cells_threshold), ' synchronous_frames = ', num2str(synchronous_frames),    ' ncluster= ' , num2str(NCl) , (  ' NRace= ') , num2str(NRace), (  ' NRaceOK= ') ,  num2str(NRaceOK) ) ;
%disp(strcat ('synchronous_frames= ', num2str(synchronous_frames),  ' minpeak distactivity', num2str(MinPeakDistance ),  ' ncluster= ' , num2str(NCl) , (  ' NRace= ') , num2str(NRace), (  ' NRaceOK= ') ,  num2str(NRaceOK) )) 
disp(test)
%save([PathSave,'NClustersOK.mat'],'NClOK')

%fileName={'new1.txt', 'new2.txt', 'new3.txt'};
%open file identifier
 fid=fopen([name,'settings.txt'],'w');
 fprintf(fid, test);
 fclose(fid);


