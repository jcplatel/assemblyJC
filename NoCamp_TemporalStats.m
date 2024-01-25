% load('TrDC')
% load('DC')
[NDC,Nt] = size(TrDC);
% load('Speed')
% load('MovT')
% Speed = circshift(speed,[0 -14]); %Delay btw ephy and image
Sp = Speed;
SpBlur=GaussBlur1d(Sp,Nt/10,2);
DCind=1:NDC;
d = cumsum(Sp)*0.1;
[~,xDel]=sort(ShiftDC);
TrB=TrDC(xDel,:);
TrBd=diff(GaussBlur1d(TrB,Nt/25,2),[],2); %Nt/20
figure
imagesc(TrBd,[0 max(TrBd(:))/3])
axis([1 1600 1 NDC])
hold on
plot(NDC-SpBlur/3,'g')
colormap hot
t = (1:Nt)/10;

%% Manually detect sequences
j=1;k=1;
while true
    [x,~,b]=ginput(2);
    disp('ok')
    if b==1
        a(:,j)=x;% *10 if added
        j=j+1;
    end
    if b(1)==3
        axis([1+k*1500 1600+k*1500 1 NDC])
        pause(0.1)
        k=k+1;
    end
    if b(1)==2
        break
    end
end

%%
Raster=zeros(NDC,Nt);
Cor1d2=[];
Shift1d2=[];
Max1d2=[];
Dist1d2=[];
NSeq=length(a);
% NSeq=NEp;
% a=round([max(S1/10-15,1);E1/10]);
for j=1:NSeq
%     Win=xLim(1,j):xLim(2,j);
    Win=round(a(1,j)):round(a(2,j));
    %re-order cells
    SigRef2=exp(-(Win-mean(Win)).^2/4^2);
    dt2=length(Win);
    dSeq{j}=d(Win)-d(Win(1));
    for i=1:NDC
        tmp=covnorm(TrBd(i,Win),SigRef2,dt2);
        Cor1d2(j,i)=max(tmp);
        if Cor1d2(j,i)>0.3
            Shift1d2(j,i)=find(tmp==Cor1d2(j,i))-dt2-1+mean(Win)-min(Win);
            Raster(i,round(Shift1d2(j,i)+Win(1)))=1;
            Dist1d2(j,i)=dSeq{j}(max(min(round(Shift1d2(j,i))+2,length(Win)),1));
            Max1d2(j,i) = TrBd(i,round(Shift1d2(j,i))+Win(1)-1);
            Speed1d2(j,i) = mean(Sp((round(Shift1d2(j,i))+Win(1)-1)+(-2:2)));
            if Max1d2(j,i)<0.002 %0.01
                Shift1d2(j,i)=nan;
                Dist1d2(j,i)=nan;
                Max1d2(j,i) = nan;
            end
        else
            Shift1d2(j,i)=nan;
            Dist1d2(j,i)=nan;
            Max1d2(j,i) = nan;
        end
    end
end

%% Normalisation
% SeqOK=nanmean(Dist1d2,2)>0.035;
Shift1d2OK=[];
Dist1d2OK=[];
for i=1:NSeq
    Shift1d2OK(i,:)=Shift1d2(i,:)-nanmedian(Shift1d2(i,:));
    Dist1d2OK(i,:)=Dist1d2(i,:)-nanmedian(Dist1d2(i,:));
end
% !!! Better with dist for H63 !!!
[~,xDelM]=sort(nanmedian(Shift1d2OK,1));
% xDelM = xDelM(DCRace);
% xDelM = xDelM(~ismember(1:NDC,[DC1]));
Shift1d2OK=Shift1d2OK(:,xDelM);
Shift1d2=Shift1d2(:,xDelM);
Dist1d2OK=Dist1d2OK(:,xDelM);
Dist1d2=Dist1d2(:,xDelM);
figure
plot(Shift1d2OK'/10,'.')
hold on
errorbar(nanmedian(Shift1d2OK,1)/10,naniqr(Shift1d2OK,1)/10/2,'.k','MarkerSize',18)
xlabel('Cell #')
ylabel('Time (s)')
figure
plot(Dist1d2OK','.')
hold on
errorbar(nanmedian(Dist1d2OK,1),naniqr(Dist1d2OK,1)/2,'.k','MarkerSize',18)
xlabel('Cell #')
ylabel('Distance (cm)')

%%
NDC = length(xDelM);
figure
% subplot(4,1,1)
imagesc(t,1:NDC,TrBd(xDelM,:),[0 max(TrBd(:))/2])
colormap hot
hold on
plot(t,NDC-SpBlur/2,'g')
for i=1:NSeq
    tmp=~isnan(Shift1d2(i,:));
    plot(t(round(a(1,i)+Shift1d2(i,tmp))),DCind(tmp),'og')
%     plot(t(round(S1(i)/10+Shift1d2(i,tmp))),DCind(tmp),'*g')
end
xlabel('Time (in s)')
axis ij