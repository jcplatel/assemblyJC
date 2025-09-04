%New method for statistically finding assemblies and RACE recruitment
%with Bonferroni correction, from clustering output 

%load('Race.mat')
%load('Clusters.mat') 
%load('NClustersOK.mat')
%load('TRace')

[NCell,NRace] = size(Race);
NCl = NClOK;
NShuf = 1000;
Nt = 14000; % For multiple movies

%% Statistiscal definition of cell assemblies

%Count number of participation to each cluster
CellP = zeros(NCell,NCl); CellR = zeros(NCell,NCl);
%here we use NClOK so we stick only to the statistically significant clusters
for i = 1:NCl
    CellP(:,i) = sum(Race(:,IDX2 == i),2);  %number of participation to each cluster
    CellR(:,i) = CellP(:,i)/sum(IDX2 == i);    %divided per nb of cell in this cluster
end

%Test for statistical significance (for highly active cell and or clusters with a lot of SCE)
CellCl = zeros(NCl,NCell, 'uint16'); %Binary matrix of cell associated to clusters
parfor j = 1:NCell
    %Random distribution among Clusters
    RClr = zeros(NCl,NShuf);
    Nrnd = sum(Race(j,:));
    for l = 1:NShuf
        Random = randperm(NRace);   % permutation des SCE
        Random = Random(1:Nrnd);    % on ne garde que le nombre de SCE de la cellule en cours
        Racer = zeros(1,NRace);
        Racer(Random) = 1;  
        for i = 1:NCl
            RClr(i,l) = sum(Racer(:,IDX2 == i),2);      %nombre de SCE qui sont dans le cluster i
        end
    end
    RClr = sort(RClr,2);
    %         ThMin = mean(Random) - 2*std(Random);
    %Proba above 95th percentile
    ThMax = RClr(:,round(NShuf*(1-0.05/NCl)));  %calcul proba: nombre de fois à 95% que cette cellule est dans le cluster i
    for i = 1:NCl
        CellCl(i,j) = double(CellP(j,i)>ThMax(i));% - double(RCl(:,j)<ThMin);       on compare nombre reel à seuil permutation 
    end
end

%save([namefull,'CellClstatraw.mat'],'CellCl')  
A0 = find(sum(CellCl) == 0); %Cells not in any cluster
A1 = find(sum(CellCl) == 1); %Cells in one cluster
A2 = find(sum(CellCl) >= 2); %Cells in several clusters
%maybe export CellCl from here before ortho  
CellClraw=CellCl;%THIS IS THE GUY WE CAN LOOK AT

%Keep cluster where they participate the most THIS IS THE PART TO CHANGE
for i = A2
    [~,idx] = max(CellR(i,:));
    CellCl(:,i) = 0;
    CellCl(idx,i) = 1;
end
C0 = cell(0);
k = 0;
for i = 1:NCl
    if length(find(CellCl(i,:)))>9 %  limit of cell per assembly , was at 5
        k = k+1;
        C0{k} = find(CellCl(i,:));
    end
end
assemblyortho=C0;
%save([namefull 'assemblyortho'],'assemblyortho')%

%Participation rate to its own cluster
CellRCl = max(CellR([A1 A2],:),[],2);
% save([namefull 'CellRCl'],'CellRCl')
%% Assign RACE to groups of cells      #####didn't get that....

NCl = length(C0);
[NCell,NRace] = size(Race);
%Cell count in each cluster
RCl = zeros(NCl,NRace);
PCl = zeros(NCl,NRace);

for i = 1:NCl
    RCl(i,:) = sum(Race(C0{i},:));   %sum sce from assembly assigned to cluster
end

RCln = zeros(NCl,NRace, 'uint16');
parfor j = 1:NRace  %permutation of cells 
    %Random distribution among Clusters
    RClr = zeros(NCl,NShuf);
    Nrnd = sum(Race(:,j));
    for l = 1:NShuf
        Random = randperm(NCell);
        Random = Random(1:Nrnd);
        Racer = zeros(NCell,1);
        Racer(Random) = 1;
        for i = 1:NCl
            RClr(i,l) = sum(Racer(C0{i}));
        end
    end
    %         ThMin = mean(Random) - 2*std(Random);
    RClr = sort(RClr,2);
    %         ThMin = mean(Random) - 2*std(Random);
    %Proba above 95th percentile
    ThMax = RClr(:,round(NShuf*(1-0.05/NCl)));
    for i = 1:NCl
        PCl(i,j) = double(RCl(i,j)>ThMax(i));% - double(RCl(:,j)<ThMin);
    end
    %Normalize (probability)
    RCln(:,j) = RCl(:,j)/sum(Race(:,j));
end
% save([namefull,'RaceCl2'],'PCl')

%% Show Sorted Rasterplot

%load('RaceCl2')
%load('CellClList2.mat')
%load('Race.mat')
[NCell,NRace] = size(Race);
NCl = length(C0);

%Recreate CellCl (equivalent of RCl for the cells)
CellCl = zeros(NCl,NCell);
for i = 1:NCl
    CellCl(i,C0{i}) = 1;
end

NCl = length(C0);

Cl0 = find(sum(PCl,1) == 0);%statistically no sce
Cl1 = find(sum(PCl,1) == 1);%only 1 assembly
Cl2 = find(sum(PCl,1) == 2);%2 assemblies
Cl3 = find(sum(PCl,1) == 3);%3 assemblies
Cl4 = find(sum(PCl,1) == 4);%4 assemblies

Bin = 2.^(0:NCl-1);

%Sort Cl1
[~,x01] = sort(Bin*PCl(:,Cl1));
Cl1 = Cl1(x01);

%Sort Cl2
[~,x02] = sort(Bin*PCl(:,Cl2));
Cl2 = Cl2(x02);

%Sort Cl3
[~,x03] = sort(Bin*PCl(:,Cl3));
Cl3 = Cl3(x03);

RList = [Cl0 Cl1 Cl2 Cl3 Cl4];
%x1 from DetectRace

[X1,x1] = sort(Bin*CellCl);  %resort assembly from this CellCl 

%Mov2 = find(TRace>Nt,1);
%Race2d = Race;
%Race2d(:, Mov2:end) = Race(:,Mov2:end)/2;

% figure
% % f = figure('visible','off');f
% imagesc(Race(x1,RList))

%Show with replay
% load('TestRACESeq')
% load('CorrRACESeq')
% figure
% RaceReplay = Race*0;RaceFwdReplay = Race*0;RaceBwdReplay = Race*0;
% RaceReplay(:,p<0.05) = Race(:,p<0.05);
% RaceFwdReplay(:,p<0.05 & r>0) = Race(:,p<0.05 & r>0);
% RaceBwdReplay(:,p<0.05 & r<0) = Race(:,p<0.05 & r<0);
% imagesc(Race(x1,RList) - RaceReplay(x1,RList)/2)
% imagesc(Race(x1,RList) - RaceFwdReplay(x1,RList)/3 - 2*RaceBwdReplay(x1,RList)/3)
% colormap hot
% axis tight
% exportgraphics(gcf,[namefull 'cell_assemblies.png'],'Resolution',300)

