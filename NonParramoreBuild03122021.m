%inputs needed to enter: 
%       IslandSize (whatever estimate you got) 
%       Desired Beach slope (using 0.02 from Davidson Arnott 'Conceptual model for...'

%name the island you wish to work on, Matlab will load up the jpg with that name
Iname='Smith3';
data=imread('Smith3.jpg');      %data reads in as an 
data2=double(data(:,:,2));      

ddata=im2double(data); %convert image to a double array with all values between 0,1

figure
colormap jet
imagesc(ddata)
colorbar
Values=unique(ddata);

A=ddata;

NdexNum=8;
[NdexMap,map6]=rgb2ind(A,NdexNum);
if 1==mod(size(NdexMap,1),2)
    NdexMap=NdexMap(1:size(NdexMap,1)-1,:);
end
if 1==mod(size(NdexMap,2),2)
    NdexMap=NdexMap(:,1:size(NdexMap,2)-1);
end

figure
imagesc(NdexMap)
title(sprintf('%s Island reduced to %d color indices',Iname,NdexNum))
colorbar

%finding minimal index - will likely be the black lines
% for i=1:6
% k=i-1;
% totals(i,1)=k;
% totals(i,2)=sum(sum(IND6==k));
% end
% KillColor=totals(find(totals(:,2)==min(totals(:,2))),1);
% disp(sprintf('removing color %d',KillColor))
% [n1,n2]=size(IND6);
% for nn=1:4
% %     nn
%     for ii=2:n1-1
%         for jj=2:n2-1
%             if IND6(ii,jj)==KillColor
%                 vec1=double([IND6(ii-1,jj) IND6(ii+1,jj) IND6(ii,jj-1) IND6(ii,jj+1)]);
%                 vec1(vec1==KillColor)=NaN;
%                 IND6(ii,jj)=mode(vec1);
%             end
%         end
%     end
% 

Problem=7; %index associated with outlines in the colormap
disp(sprintf('Removing color %d',Problem));
[n1,n2]=size(NdexMap);
for nn=1:5
%     nn
    for ii=2:n1-1
        for jj=2:n2-1
            if NdexMap(ii,jj)==Problem
                vec1=double([NdexMap(ii-1,jj) NdexMap(ii+1,jj) NdexMap(ii,jj-1) NdexMap(ii,jj+1)]);
                vec1(vec1==Problem)=NaN;
                NdexMap(ii,jj)=mode(vec1);
            end
        end
    end
end

Problem=4; %index associated with outlines in the colormap
disp(sprintf('Removing color %d',Problem));
[n1,n2]=size(NdexMap);
for nn=1:5
%     nn
    for ii=2:n1-1
        for jj=2:n2-1
            if NdexMap(ii,jj)==Problem
                vec1=double([NdexMap(ii-1,jj) NdexMap(ii+1,jj) NdexMap(ii,jj-1) NdexMap(ii,jj+1)]);
                vec1(vec1==Problem)=NaN;
                NdexMap(ii,jj)=mode(vec1);
            end
        end
    end
end


figure
imagesc(NdexMap)
title(sprintf('%s Island with index %d (outlines) removed',Iname,Problem))
colorbar

% Reduce the number of pixels by factor of 2 in each direction (4 overall)
% % % % disp('Reducing number of pixels');
% % % % A2=zeros(size(NdexMap)/2);
% % % % [m2,n2]=size(A2);
% % % % for ii=1:m2
% % % %     for jj=1:n2
% % % %         A2(ii,jj)=mode(NdexMap(2*ii-1:2*ii,2*jj-1:2*jj),'all');
% % % %     end
% % % % end
% % % % figure
% % % % imagesc(A2)
% % % % title(sprintf('%s Island, 6 colors, 4x fewer pixels',Iname))
% % % % colorbar
% % % % 
% % % % NdexMap=A2;

%rescale the image;
[n1,n2]=size(NdexMap);
IslandSize=9065000; %island surface area in square meters
% ScalingFactor=IslandSize/(n1*n2);
H2O=mode(mode(NdexMap));                %most common number on map will be water
NonH2O=sum(sum(NdexMap~=H2O));          %just marsh and island area
ScalingFactor=sqrt(IslandSize/NonH2O);  %scale estimated area of the island by current area of island in image

NdexMap_temp=imresize(NdexMap,ScalingFactor);
NdexMap=NdexMap_temp;

% Replace the color 2 with 5  
% IND6(IND6==2)=5;

% % Replace the color 4 with whatever occurs the most around it
% disp('Replacing color 4');
% for nn=1:1
% %     nn
%     for ii=2:n1-1
% %         ii
%         for jj=2:n2-1
%             vec1=double([IND6(ii-1,jj) IND6(ii+1,jj) IND6(ii,jj-1) IND6(ii,jj+1)]);
%             vec1(vec1==4)=NaN;
%             [m1,f]=mode(vec1);
%             if f==3 || f==4
%                 IND6(ii,jj)=m1;
%             end
%         end
%     end
% end

% This section smooths out pixels that don't match the stuff around it in
% any of four directions - run it a couple times
nnend=2;
dd=6; % (How many cells do you want to look at in each direction?)
for nn=1:nnend
%     nn
    IND6C=NdexMap;
    disp(strcat(['N,S,E,W check ',int2str(nn),'/',int2str(nnend)]));
    % First check north, south, east, west...
    for ii=dd:n1-dd+1
%         ii
        for jj=dd:n2-dd+1
            m1=mode(IND6C(ii-dd+1:ii-1,jj)); m2=mode(IND6C(ii+1:ii+dd-1,jj));
            m3=mode(IND6C(ii,jj-dd+1:jj-1)); m4=mode(IND6C(ii,jj+1:jj+dd-1));
            if IND6C(ii,jj)~=m1 && IND6C(ii,jj)~=m2 && IND6C(ii,jj)~=m3 && IND6C(ii,jj)~=m4
                NdexMap(ii,jj)=mode([m1,m2,m3,m4]);
            end
        end
    end
    IND6C=NdexMap;
    disp(strcat(['NE,NW,SE,SW check ',int2str(nn),'/',int2str(nnend)]));
    % Then check northeast, northwest, southeast, and southwest
    for ii=dd:n1-dd+1
%         ii
        for jj=dd:n2-dd+1
            m1=mode(diag(IND6C(ii-dd+1:ii-1,jj-dd+1:jj-1))); m2=mode(diag(IND6C(ii+1:ii+dd-1,jj-dd+1:jj-1)));
            m3=mode(diag(IND6C(ii-dd+1:ii-1,jj+1:jj+dd-1))); m4=mode(diag(IND6C(ii+1:ii+dd-1,jj+1:jj+dd-1)));
            if IND6C(ii,jj)~=m1 && IND6C(ii,jj)~=m2 && IND6C(ii,jj)~=m3 && IND6C(ii,jj)~=m4
                NdexMap(ii,jj)=mode([m1,m2,m3,m4]);
            end
        end
    end
end

figure
imagesc(NdexMap)
title(sprintf('%s Island, outlines removed and smoothed',Iname))
colorbar

% for ii=1:n1
%     for jj=1:n2
%         if NdexMap(ii,jj)==6
%             NdexMap(ii,jj)=3;
%         end
%     end
% end

%%
%working with the 6 color matrix to produce an elevation map with desireable orientation
%if IND6 has already been saved, and you just need to use this subsection
%for tweaking the image:
% % % % Iname='Cobb3'; %<<--enter island name
% % % % filename=sprintf('%s_6color.mat',Iname);
% % % % NdexMap=importdata(filename);
% % % % %else, comment out above 3 lines.
% % % % 
% % % % %roate by desired number of degrees (may take some guessing) then reset
% % % % %added 0 cells as more ocean
% % % % IND6plus=NdexMap(:,:)+1;
% % % % angle=35; %enter degress to rotate counter-clockwise
% INDrot=imrotate(IND6plus,angle);
% INDrot(INDrot==0)=6;
% INDrot=INDrot(:,:)-1;
% IND6=INDrot;

X=NdexMap; %making a copy to play with

H=zeros(size(X,1),size(X,2));   %hopefully will be the new elevation matrix

H2O=mode(mode(NdexMap));                %most common number on map will be water
for ii=1:size(X,1)
    for jj=1:size(X,2)
        if X(ii,jj)==H2O %all water area (white)
            H(ii,jj)=0;
        else %X(ii,jj)~=5  % subaerial plus marsh
            H(ii,jj)=10;
        end
    end
end


Hnew=H;
%now we need the indexes of the boundary cells and remove half oz zeros
%cols and rows some number of times



WE_shore=zeros(size(Hnew,1),2);
for ii=1:size(Hnew,1)
    RowNow=Hnew(ii,:);
    check1=sum(RowNow>0);
    if check1~=0
        WE_shore(ii,1)=find(RowNow>0,1,'first');
        WE_shore(ii,2)=find(RowNow>0,1,'last');
    end
end

MaxWest=(min(WE_shore(WE_shore~=0)));
MaxEast=(max(WE_shore(WE_shore~=0)));
NS_shore=zeros(size(Hnew,2),2);
for jj= 1:size(Hnew,2)
    ColNow=Hnew(:,jj);
    check2=sum(ColNow>0);
    if check2 ~=0
        NS_shore(jj,1)=find(ColNow>0,1,'first');
        NS_shore(jj,2)=find(ColNow>0,1,'last');
    end
end

Wv=WE_shore(:,1); Ev=WE_shore(:,2); %vectors of col index for west/east shore
Nv=NS_shore(:,1); Sv=NS_shore(:,2); %vectos of row numbers for north/south shore
W0=min(Wv(Wv>0));E0=max(Ev(Ev>0)); %min max col cell(western-most and eastern-most)
N0=min(Nv(Nv>0));S0=max(Sv(Sv>0)); %max minrow cell (northern-most and southern-most)


% HO=zeros(size(H,1),size(H,2)); %produces an outline image, just to check
% for ii=1:length(WE_shore)
%     if WE_shore(ii,1)~=0 && WE_shore(ii,2)~=0
%         HO(ii, WE_shore(ii,1))=1;HO(ii, WE_shore(ii,2))=1;
%     end
% end
% for jj=1:length(NS_shore)
%     if NS_shore(jj,1)~=0 && NS_shore(jj,2)~=0
%         HO(NS_shore(jj,1),jj)==1;HO(NS_shore(jj,2),jj)==1;
%     end
% end
% figure
% imagesc(HO)
% 
% DDrop=0.1;  %elevation decrease (change in height)
% PerLength=10; %nomber of cells to skip before decreasing again (change in distance)
% dy=DDrop;
% dx=PerLength;
% 
% Hcopy=H;
% 
% for n=1:length(WE_shore)
%     RowNow=H(n,:);
%     west=WE_shore(n,1);
%     east=WE_shore(n,2);
%     wcnt=0;
%     if west~=0
%         RowNow(west-1)=-1;
%         for kk=-2:-1:-west+1
%             RowNow(west+kk)=-1+RowNow(west+kk)+dy*(kk/dx);
%         end
%     end
% %     H(n,:)=RowNow;
%     if east~=0
%         RowNow(east+1)=-1;
%         for kk=1:(length(RowNow)-east)
%             RowNow(east+kk)=-1+RowNow(east+kk)+dy*(-kk/dx);
%         end
%     end
%     H(n,:)=RowNow;
% end
% 
% %copy first and last row to get gradient change above, replace the above
% %water cells with -0.1 and use the same process as above on the whole row
% %going up and down from the north and south tip
% RstarN=find((WE_shore>0),1,'first');
% RowStarN=H(RstarN,:);
% DelInd=find(RowStarN>=0);
% Mn=mean(DelInd);
% RowStarN(DelInd)=-1;
% for n=-1:-1:-RstarN+1
%     H(RstarN+n,:)=RowStarN(:)+dy*(n/dx);
% end
% RstarS=find((WE_shore(:,1)>0),1,'last');
% RowStarS=H(RstarS,:);
% DelInd=find(RowStarS>=0);
% RowStarS(DelInd)=-1;
% for n=1:length(WE_shore)-RstarS
%     H(RstarS+n,:)=RowStarS(:)+dy*(-n/dx);
% end
% 
% figure
% imagesc(H)
DDrop=0.1;  %elevation decrease (change in height)
PerLengthEast=100; %nomber of cells to skip before decreasing again (change in distance)
PerLengthWest=500;
dy=DDrop;
dxE=PerLengthEast;
dxW=PerLengthWest;
dxNS=floor(.5*(dxE+dxW));

Hcopy=H;

for n=1:length(WE_shore)
    RowNow=H(n,:);
    west=WE_shore(n,1);
    east=WE_shore(n,2);
    wcnt=0;
    if west~=0 && west~=1
        RowNow(west-1)=-1;
        for kk=-2:-1:-west+1
            RowNow(west+kk)=-1+RowNow(west+kk)+dy*(kk/dxW);
        end
    end
%     H(n,:)=RowNow;
    if east~=0 &&east~=size(H,2)
        RowNow(east+1)=-1;
        for kk=1:(length(RowNow)-east)
            RowNow(east+kk)=-1+RowNow(east+kk)+dy*(-kk/dxE);
        end
    end
    H(n,:)=RowNow;
end

%copy first and last row to get gradient change above, replace the above
%water cells with -0.1 and use the same process as above on the whole row
%going up and down from the north and south tip
RstarN=find((WE_shore>0),1,'first');
RowStarN=H(RstarN,:);
DelInd=find(RowStarN>=0);
Mn=mean(DelInd);
RowStarN(DelInd)=-1;
for n=-1:-1:-RstarN+1
    H(RstarN+n,:)=RowStarN(:)+dy*(n/dxNS);
end
RstarS=find((WE_shore(:,1)>0),1,'last');
RowStarS=H(RstarS,:);
DelInd=find(RowStarS>=0);
RowStarS(DelInd)=-1;
for n=1:length(WE_shore)-RstarS
    H(RstarS+n,:)=RowStarS(:)+dy*(-n/dxNS);
end

figure
imagesc(H)

    
Helev=zeros(size(H,1),size(H,2));
nearshore=-1:0.1:-.01;
marsh=-.5:.1:.5;
LowE=1:.1:2.5;
HighE=2.5:.1:5;
% BchWidth=20;%randi([10,20],5,1); %dx
Slope=.02; %desired beach slope (Davidson Arnott 'Conceptual model on the effects of sea level rise...'
BchMax=1.5; %dy
BchWidth=BchMax/Slope;
for ii=1:size(Helev,1)
    for jj=1:size(Helev,2)
        if X(ii,jj)==0
            Helev(ii,jj)=marsh(randi([1,numel(marsh)]));
%             Helev(ii,jj)=-.5;
        elseif X(ii,jj)==1 || X(ii,jj)==2 || X(ii,jj)==3
            Helev(ii,jj)=LowE(randi([1,numel(LowE)]));
%             Helev(iri,jj)=1;
        elseif X(ii,jj)==6||X(ii,jj)==7 
            Helev(ii,jj)=HighE(randi([1,numel(HighE)]));
%             Helev(ii,jj)=3;
        end
        if H(ii,jj)==10
            H(ii,jj)=0;
        end
        west=WE_shore(ii,1);
        east=WE_shore(ii,2);
        Iwidth=east-west;
        if WE_shore(ii,1)~=0 && WE_shore(ii,2)~=0
            for jj = 0:BchWidth;
                Helev(ii,east-jj)=0+BchMax*(jj/BchWidth);
                H(ii,east-jj)=0;
            end
        end
    end
end

% figure
% imagesc(H)
% figure
% imagesc(Helev)

H=round((H+Helev)*10,0);
figure
imagesc(H)
colorbar
title(sprintf('Elevation Matrix H for %s Island',Iname))
% save([Iname '_elevation'],'H');
%        
%%  Polishing off teh water gradient with a Gaussian filter


% N=zeros(size(H,1),size(H,2));
% M=zeros(size(H,1),size(H,2));
% P=zeros(size(H,1),size(H,2));
% for i=1:size(H,1)
%     for j=1:size(H,2)
%         if H(i,j)>=-5 && H(i,j)<=BchMax
%             M(i,j)=H(i,j);
%         end
%     end
% end
% [rowN colN]=find(H<-5);
% [rowP colP]=find(H>BchMax);
% for i=1:length(rowN)
% R=rowN(i);C=colN(i);
% N(R,C)=H(R,C);
% end
% for i=1:length(rowP)
%     R=rowP(i); C=colP(i);
%     P(R,C)=H(R,C);
% end
% Nfilt=imgaussfilt(N,10);
% Mfilt=imgaussfilt(M,2);
% Pfilt=imgaussfilt(P,3);
% Hfilt=Nfilt+Mfilt+Pfilt;


[rowN colN]=find(H<1);
[rowP colP]=find(H>=1);
N=zeros(size(H,1),size(H,2));
P=zeros(size(H,1),size(H,2));
for i=1:length(rowN)
R=rowN(i);C=colN(i);
N(R,C)=H(R,C);
end
for i=1:length(rowP)
    R=rowP(i); C=colP(i);
    P(R,C)=H(R,C);
end
Nfilt=imgaussfilt(N,50);
Pfilt=imgaussfilt(P,3);
Hfilt=Nfilt+Pfilt;

figure
imagesc(Hfilt)
title(sprintf('Elevation map for %s Island with water gradient',Iname))
colorbar

Hround=round(H,0);
% [rowLT colLT]=find(Hfilt<=BchMax);
% [rowGT colGT]=find(Hfilt>BchMax);
% LowLand=zeros(size(H,1),size(H,2));
% HighLand=zeros(size(H,1),size(H,2));
% for i=1:length(rowLT)
%     R=rowLT(i);C=colLT(i);
%     LowLand(R,C)=Hfilt(R,C);
% end
% for i=1:length(rowGT)
%     R=rowGT(i);C=colGT(i);
%     HighLand(R,C)=Hfilt(R,C);
% end
% LLfilt=imgaussfilt(LowLand,3);
% H_HighLow=LLfilt+HighLand;
% 
% figure
% imagesc(H_HighLow)
% colorbar
% 
% 
% 
% delta=.1;
% L=1;
% Hnew=BeachBuilder(H_HighLow,delta,L,BchMax,BchWidth);

% delta=.1;
% L=1;
% Hnew=BeachBuilder(Hfilt,delta,L,BchMax,BchWidth);

figure
imagesc(Hfilt)
title(sprintf('Elevation map for %s Island with water gradient and beach',Iname))
colorbar

FinalH=round(Hfilt);
% Hstar=FinalH;
% delta=.1;
% L=1;
% PC=Hstar*0;
% t=0;
% flag=1;
% cnt=0;
%     while cnt<30
%         [Hstar,flag,CellCt]=AVALANCHEtime03062021(Hstar,delta,L,flag,PC,t);
%         if flag==0
%             cnt=cnt+1;
%         else
%             cnt=0;
%         end
%     end

% save([Iname '_elevation'],'H');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


function N1=NewNhd(R,C,H)
if R<=(size(H,1)-1)
    R=floor(R);
    Rp=R+1;
else
    Rp=floor(R);
end
if R>=2
    R=floor(R);
    Rm=R-1;
else
    Rm=floor(R);
end
if C<=(size(H,2)-1)
    C=floor(C);
    Cp=C+1;
else
    Cp=floor(C);
end
if C>=2
    C=floor(C);
    Cm=C-1;
else
    Cm=floor(C);
end
N1=zeros(1,8);
%     PB4=PA(:,:,i);
N1=[H(Rm,Cm),H(Rm,C),H(Rm,Cp),H(R,Cp),H(Rp,Cp),H(Rp,C),H(Rp,Cm),H(R,Cm)];

% for k=1:8
%     if N1(k)==-999
%         N1(k)=0;
%         if Rm==R
%     elseif N1(k)<0.01
%         N1(k)=0;
%     end
% end
end




%%%%%
% %ekevation assignments before I changed things
% for ii=1:size(X,1)
%     for jj=1:size(X,2)
%         if X(ii,jj)==5
%             H(ii,jj)=-5;
%         elseif X(ii,jj)==2 || X(ii,jj)==4
%             H(ii,jj)=nearshore(randi([1,numel(nearshore)]));
%         elseif X(ii,jj)==0
%             H(ii,jj)=marsh(randi([1,numel(marsh)]));
%         elseif X(ii,jj)==1
%             H(ii,jj)=LowE(randi([1,numel(LowE)]));
%         elseif X(ii,jj)==3
%             H(ii,jj)=HighE(randi([1,numel(HighE)]));
%         end
%     end
% end
