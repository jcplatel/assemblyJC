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
    [amplitude,locs] = findpeaks(Tr1b(i,:),'MinPeakProminence',th(i) ,'MinPeakDistance',MinPeakDistance);%2.576= 99%, 3.291=99.9
    valeurs_identiques = intersect (locs,WinActive);
    [locs_sans_ide , idx ]=setdiff(locs(:), valeurs_identiques);
    ampli_sans_ide=amplitude(idx);

    Acttmp2{i}=locs_sans_ide;%%%%%%%%findchangepts(y,MaxNumChanges=10,Statistic="rms")
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


sce_n_cells_threshold = 10;

[pks,TRace] = findpeaks(MAct,'MinPeakHeight',sce_n_cells_threshold,'MinPeakDistance',MinPeakDistancesce);

NRace = length(TRace);

disp(['nSCE: '  num2str(NRace)])

% Create RasterPlots%%%%%%very weird here increase RACE from n-1 :n+2 
Race = zeros(NCell,NRace);      %Race=cells that participate in SCE 
RasterRace = zeros(NCell,Nz);
for i = 1:NRace
    Race(:,i) = max(Raster(:,TRace(i)-1:TRace(i)+2),[],2);    %maybe this window can be optimized ???
    RasterRace(Race(:,i)==1,TRace(i)) = 1;                          %Raster with only activity in SCE time: all frames
end
