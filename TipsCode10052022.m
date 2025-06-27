function [Htips] = TipsCode10052022(Hstar)
%filename='Parramore03312021.mat';
%filename='testland2copy.txt';
% filename = 'SmithNew06132021.mat';
% H=importdata(filename);

H = Hstar(:,:,1);
[n1, n2]=size(H(:,:));
ESL=zeros(n1,1);
WSL=zeros(n1,1);
DC=zeros(n1,1);
dH= zeros(n1,1);
xDT=zeros(n1,1);  %%%% x coordinates of Dune Top
yDT=zeros(n1,1);  %%%% y coordinates of Dune Top
xDB=zeros(n1,1);  %%%% x coordinates of Dune Bottom
yDB=zeros(n1,1);  %%%% y coordinates of Dune Bottom

%% (Code 1) To determine the eastern shore line (ESL)
for i=1:n1
    for j=n2:-1:1
        if H(i,j)>0
            ESL(i)=j; %%%%%%%%% Eastern shore line (ESL)%%%%%%%%% gives column
            break
       end
    end
end

% To determine the Western Shore Line
for i=1:n1
    for j=1:n2
        if H(i,j)>=0
            WSL(i)=j; %%%%%%%%% Western shore line (WSL)%%%%%%%%%
            break
       end
    end
end


%% (Code 2) To determine the j-th column of first dune crest at DC(i) from Eastern Shore line 
% To determine the dune height at H(i,j) from DC(i)=j%%%%%%
for i=1:n1
    if ESL(i)>0
        shore=ESL(i);
        for j=shore:-1:1
            if H(i,j)>H(i,j-1)
                DC(i)=j; %%%%%%%%%% DC=dune crest %%%%%%%%%%
               dH(i)= H(i,j); %%%%%%% dH= dune heights %%%%%%%
                
                break
            end
        end
    end
end

%% (Code 3) To find the maximum height of dune crest from dH 
%To find the x coordinate of dune that start from top of the island(Ytop)
%To find the x coordinate of dune that ends at the bottom (Ybot)


Ymax= max(dH); %%%% maximum height of dune crest
Ytop = find((dH>0),1,'first'); %%%%% x coordinate where the dune start from top
Ybot=  find ((dH>0),1, 'last'); %%%%% x coordinate where the dune ends at the bottom 

%%% (x,y) coordinates of Ytop: 
xYtop=Ytop;  
yYtop=ESL(Ytop); %% y coordinate of ESL at x where the dune start from top of the island
XYTop= [xYtop yYtop]; 

%%% (x,y) coordinates of Ybot: 

xYbot=Ybot;
yYbot=ESL(Ybot);  %% y coordinate of ESL at x where the dune ends at the bottom of the island
XYBot= [xYbot yYbot];


 %% (Code 4) To find the x and y coordinates of Dune Crest 

for i=2:n1-1
    if dH(i)>=dH(i-1)
  
    
        if dH(i) > dH (i+1)
           
        xDT(i) = i; %%%% X coordinates of dune tops= xDT
        yDT(i) = DC(i); %%%%% y coordinates of dune top = yDT
        
        
        end 
    end
end


  
xDT(xDT==0)=[]; %%% xDT without zeros
yDT(yDT==0)=[]; %%% yDT without zeros

%% (Code 5) To find the x and y coordinates of Dune bottom
          
for i=2:n1-1
    if dH(i)<dH(i-1)
            
            if dH(i)<= dH (i+1)
            if dH(i)~=0
            
           

           
        xDB(i) = i; %%%%%% X coordinates of dune bottoms = xDB
        yDB(i) = DC(i); %%%%%%% y coordinates dune bottoms= yDB
            end
            end
            
    end
end

xDB(xDB==0)=[]; %%%%%% xDB without zeros
yDB(yDB==0)=[]; %%%%%% yDB without zeros

% figure
% imagesc(H)
% hold on
% % scatter(ESL,1:n1,'ro')
% scatter(DC,1:n1,'yo')  
% title('before')
% colorbar
% caxis([-5 26])
 

%% (Code 6) To create new matrix MNew with x and y coordinates starting from top of the island,dune tops(MAX) and dune bottoms(MIN)(alternating) and bottom of the island.

xB= xDB; %%%%% whichever length is shorter needs to be xB and needs to make sure xDT is the start of new combined vector
         
xT=xDT;  %%%%% xDT should be one more than xDB

LenxB = length(xB);
C = zeros(1, LenxB+length(xT));
C(2:2:LenxB*2) = xB;
C(1:2:LenxB*2) = xT(1:LenxB);
C(LenxB*2+1:end) = xT(LenxB+1:end);


yB= yDB; %%%%% whichever length is shorter needs to be E and needs to make sure yt is the start of new combined vector
yT= yDT;   %%% yDT should be 1 more than yDB

LenyB= length(yB);
D = zeros(1, LenyB+length(yT));
D(2:2:LenyB*2) = yB;
D(1:2:LenyB*2) = yT(1:LenyB);
D(LenyB*2+1:end) = yT(LenyB+1:end); 

M= [C; D]';


MNew = [XYTop; M; XYBot];%%% includes coordinates of top and bottom

%% (Code 7 Beth)To create the effects of the cold current running N to S
Htips = Hstar(:,:,1);
n = length(MNew);
norbandstart = MNew(1,1);
southbandstop = MNew(n,1) ;

% Defining stop point for N band and start point for S band
for i = 1:50
if MNew(i+1,1) - MNew(1,1) >= 15 %this distance can change
norbandstop = MNew(i+1,1);
break
end
end
for i = 1:50
if MNew(n,1) - MNew(n-i,1)>= 15 %this distance can change
southbandstart = MNew(n-i,1);
break
end
end

delta1 = 2 ; % amount of sand moved

%NORTH
%begin by defining the band removal width
ESL_north = ESL(norbandstart:norbandstop);
WSL_north = WSL(norbandstart:norbandstop);
ESLmax_north = max(ESL_north);
WSLmin_north = min((WSL_north(WSL_north>0)));
band = floor((ESLmax_north-WSLmin_north)/2);
%removing from norbandstart-5 to norbandstop
for i = (norbandstart - 10) : norbandstop
for j = (ESLmax_north - band) : (ESLmax_north + 10) %ESL
if Hstar(i,j) >= -4 && Hstar(i,j)<= 2
Htips(i,j) = Hstar(i,j) - delta1 ;
end
end
for j = (WSLmin_north - 10) : (WSLmin_north + band) %WSL
if Hstar(i,j) >= -4 && Hstar(i,j)<= 2
Htips(i,j) = Hstar(i,j) - delta1 ;
end
end
end

%removed_N = sum((Hstar(:,:,1) - Htips), 'all'); %calculates how much was removed

%SOUTH
%begin by defining the band deposit width
ESL_south = ESL(southbandstart:southbandstop);
WSL_south = WSL(southbandstart:southbandstop);
ESLmax_south = max(ESL_south);

WSLmin_south = min(WSL_south(WSL_south>0));
band = floor((ESLmax_south-WSLmin_south)/2);
%depositing from southbandstart to southbandstop+5
for i = southbandstart: (southbandstop+10)
for j = (ESLmax_south- band) : (ESLmax_south + 10) %ESL
if Hstar(i,j) >= -4 && Hstar(i,j)<= 2
Htips(i,j) = Hstar(i,j) + delta1 ;
end
end
for j = (WSLmin_south - 10) : (WSLmin_south + band) %WSL
if Hstar(i,j) >= -4 && Hstar(i,j)<= 2
Htips(i,j) = Hstar(i,j) + delta1 ;
end
end
end

%added_S = (sum((Hstar(:,:,1)-Htips),'all') + removed_N); %calculates how much was added

Hstar(:,:,1) = Htips;

% % before and after figure %commented out for main code
% tiledlayout(1,2);
% nexttile
% imagesc(Hstar(:,:,1))
% title('Before')
% colorbar
% clim([-5 26])
% nexttile
% imagesc(Htips)
% title('After')
% colorbar
% clim([-5 26])
end
    