%default values
%MinPeakDistancesce=5 frames= 500ms
% synchronous_frames=2
%MinPeakDistance=5 frames= 500ms
%sce_n_cells_threshold=5 or 5% ???   seem best results between 10 to 15
%cells  should we simulate ???

clear
close all

for sce_n_cells_threshold=8:12%[5,10,15,20]%5:30%

for z=1:3

tic
%% Load settings
MinPeakDistancesce=5 ;
synchronous_frames=2;
MinPeakDistance=5;
%sce_n_cells_threshold = 18;
imaging_sampling_rate=15.7;
%synchronous_frames=round(0.2*imaging_sampling_rate,0); %200ms *sampling rate

% PathSave='/Users/platel/Desktop/exp/analysis/';
PathSave='/Users/platel/Desktop/exp/erwan/analysis/';

 dt = datestr(now,'yy_mm_dd_HH_MM');
 name=[PathSave dt 'assembly_ErwanP23D_230201_230216_17_230224/'];
 mkdir (name)    % make folder for saving analysis

%% Load current data
% 
% speed = readNPY('/Users/platel/Desktop/exp/speed175plane0.npy');
% %speed=readNPY('/Users/platel/Desktop/exp/aurelie/152/speed152plane0.npy');
%speed=readNPY('/Volumes/Crucial X6/444112/speed444112.npy');
%speed=readNPY('/Volumes/Crucial X6/444118/speed444118_220912.npy');

 % plane0=readNPY('/Users/platel/Desktop/exp/aurelie/suite2p_444175_221125_plane0/F.npy');
 % plane1=readNPY('/Users/platel/Desktop/exp/aurelie/suite2p_444175_221125_plane1/F.npy');
% plane0=readNPY('/Users/platel/Desktop/exp/aurelie/152/Fplane0.npy');
% plane1=readNPY('/Users/platel/Desktop/exp/aurelie/152/Fplane1.npy');
% plane0=readNPY('/Volumes/Crucial X6/444112/suite2p_444112_220927_plane0/F.npy');
% plane1=readNPY('/Volumes/Crucial X6/444112/suite2p_444112_220927_plane1/F.npy');
% plane0=readNPY('/Volumes/Crucial X6/444118/suite2p_444118_220912_plane0/F.npy');
% plane1=readNPY('/Volumes/Crucial X6/444118/suite2p_444118_220912_plane1/F.npy');

% iscell0=readNPY('/Users/platel/Desktop/exp/aurelie/suite2p_444175_221125_plane0/iscell.npy');
% iscell1=readNPY('/Users/platel/Desktop/exp/aurelie/suite2p_444175_221125_plane1/iscell.npy');
% iscell0=readNPY('/Users/platel/Desktop/exp/aurelie/152/iscellplane0.npy');
% iscell1=readNPY('/Users/platel/Desktop/exp/aurelie/152/iscellplane1.npy');
% iscell0=readNPY('/Volumes/Crucial X6/444112/suite2p_444112_220927_plane0/iscell.npy');
% iscell1=readNPY('/Volumes/Crucial X6/444112/suite2p_444112_220927_plane1/iscell.npy');
% % iscell0=readNPY('/Volumes/Crucial X6/444118/suite2p_444118_220912_plane0/iscell.npy');
% % iscell1=readNPY('/Volumes/Crucial X6/444118/suite2p_444118_220912_plane1/iscell.npy');

% speed=readNPY('/Users/platel/Desktop/exp/erwan/speedsuite2p_230215_230302_12_230307.npy');
% plane0=readNPY('/Users/platel/Desktop/exp/erwan/suite2p_230215_230302_12_230307_plane0/F.npy');
% iscell=readNPY('/Users/platel/Desktop/exp/erwan/suite2p_230215_230302_12_230307_plane0/iscell.npy');
speed=readNPY('/Users/platel/Desktop/exp/erwan/speedP23D_230201_230216_17_230224.npy');
plane0=readNPY('/Users/platel/Desktop/exp/erwan/suite2p_230201_230216_17_230224_plane0/F.npy');
iscell=readNPY('/Users/platel/Desktop/exp/erwan/suite2p_230201_230216_17_230224_plane0/iscell.npy');


%suite2p_230215_230302_12_230307_plane0
% tic
% k=18000;
% for n=1:k
%     tmp=imread('/Users/platel/Desktop/exp/Sofia/brainbow ani51/TSeries-03282023-1355_940_135umdeep-grid.tif','index',n);
%     Flin(:,n) = reshape(tmp, [], 1);
%     Fraw(:,:,n)=tmp;
% end


Tr1b=double(plane0(iscell(:,1)>0,:));
% plane1=double(plane1(iscell1(:,1)>0,:));
% 
% % Tr1b=vertcat(plane0, plane1);
% Tr1b=double(Flin(:,:));
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
% [NCell,Nz] = size(Tr1b);
% Tr1b = Tr1b(:,WinRest);
%TrRest = Tr1b;

figure
for i = 1:NCell
    plot(MovT,Tr1b(i,:)+i-1)
    hold on
end

% Detect Calcium Transients using a sliding window
%Tr1b = Tr1b(:,WinRest);
%TrRest = Tr1b;
Raster = zeros(NCell,Nz);
%WinSize = 40;
%WinSize = round(4 * imaging_sampling_rate);

% Acttmp = zeros(1,Nz);
% Sigtmp = zeros(1,Nz);

Acttmp2 = cell(1,NCell);

parfor i=1:NCell    

    %MinPeakDistance=5;  % for cell activities they can be close
    th=3*iqr(Tr1b(i,:));
    %th(i)=3*std(Tr1b(i,:));
    [~,locs] = findpeaks(Tr1b(i,:),'MinPeakProminence',th,'MinPeakDistance',MinPeakDistance);
    Acttmp2{i}=locs;%%%%%%%%findchangepts(y,MaxNumChanges=10,Statistic="rms")

end
%figure
for i = 1:NCell
    Raster(i,Acttmp2{i}) = 1;
        
    plot(MovT,Tr1b(i,:)+i-1)
    hold on
    plot(MovT(Acttmp2{i}),Tr1b(i,Acttmp2{i})+i-1,'.r')
end

% Sum activity over n (synchronous_frames ) consecutive frames

MAct = zeros(1,Nz-synchronous_frames);
for i=1:Nz-synchronous_frames
    MAct(i) = sum(max(Raster(:,i:i+synchronous_frames),[],2));
end

% Select synchronies (RACE)
%MinPeakDistancesce=5;
[pks,TRace] = findpeaks(MAct,'MinPeakHeight',sce_n_cells_threshold,'MinPeakDistance',MinPeakDistancesce);
% sumpeaks=sum(peaks)
NRace = length(TRace);

% Create RasterPlots%%%%%%very weird here increase RACE from n-1 :n+2 
Race = zeros(NCell,NRace);
RasterRace = zeros(NCell,Nz);
for i = 1:NRace
    Race(:,i) = max(Raster(:,TRace(i)-1:TRace(i)+2),[],2);
    RasterRace(Race(:,i)==1,TRace(i)) = 1;
end

% % Display race
% for i = 1:length(TRace)
%     line(MovT(TRace(i))*[1 1],[0 NCell+1],'Color','g');
% end
%toc
%return
% break
%% Save
  save([name,'Acttmp2.mat'],'Acttmp2')
  save([name,'Race.mat'],'Race')
  save([name,'TRace.mat'],'TRace')



%% Clustering

[NCell,NRace] = size(Race);
%[IDX2,sCl,M,S] = kmeansopt(Race,10,'var');
[IDX2,sCl,M,S] = kmeansopt(Race,10,'var');
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
 save([name,'Clusters.mat'],'IDX2')


%% Remove cluster non-statistically significant

sClrnd = zeros(1,100);

for i = 1:20
    sClrnd(i) = kmeansoptrnd(Race,10,NCl);  %was at 10 but too variable from one time to the other
end
NClOK = sum(sCl>max(sClrnd));
sClOK = sCl(1:NClOK)';

save([name,'NClustersOK.mat'],'NClOK')

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
test=strcat (' sce_n_cells_threshold=', num2str( sce_n_cells_threshold),  ' nassemblies= ' , num2str(NCl) , (' NRace= ') , num2str(NRace), (  ' NRaceOK= ') ,  num2str(NRaceOK) ) ;
%disp(strcat ('synchronous_frames= ', num2str(synchronous_frames),  ' minpeak distactivity', num2str(MinPeakDistance ),  ' ncluster= ' , num2str(NCl) , (  ' NRace= ') , num2str(NRace), (  ' NRaceOK= ') ,  num2str(NRaceOK) )) 
disp(test)
%save([PathSave,'NClustersOK.mat'],'NClOK')

%fileName={'new1.txt', 'new2.txt', 'new3.txt'};
%open file identifier
 fid=fopen([name,'settings.txt'],'w');
 fprintf(fid, test);
 fclose(fid);

end
end
