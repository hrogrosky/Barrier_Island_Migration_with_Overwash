function [Hstar] = TipsCode10052022(Hstar)
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
        if H(i,j)>0
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
 


%% (Code 7 Beth) To create the effects of the cold current running N to S
delta = 1 ; % amount of sand moved

%%%%%% NORTH %%%%%%

for i = xYtop-25:xYtop+25
    % Eastern shore line
    if ESL(i) ~= 0
    for j = ESL(i)+5:-1:(ESL(i)-5)
        if H(i,j) >= -4 && H(i,j)<= 2
            H(i,j) = H(i,j) - delta;
        end
    end
    end
    
    % Western shore line
    if WSL(i) ~= 0
    for j = WSL(i)+5:-1:(WSL(i)-5)
        if H(i,j) >= -4 && H(i,j)<= 2
            H(i,j) = H(i,j) - delta;
        end
    end
    end
end

%%%%%%% SOUTH %%%%%%

for i = xYbot-25:xYbot+25
    % Eastern shore line
    if ESL(i) ~= 0
    for j = ESL(i)+5:-1:(ESL(i)-5)
        if H(i,j) >= -4 && H(i,j)<= 2
            H(i,j) = H(i,j) + delta;
        end
    end
    end
    
    % Western shore line
    if WSL(i) ~= 0
    for j = WSL(i)+5:-1:(WSL(i)-5)
        if H(i,j) >= -4 && H(i,j)<= 2
            H(i,j) = H(i,j) + delta;
        end
    end
    end
end

H = Hstar;
end
   
%just one figure
% figure
% imagesc(H(:,:))
% hold on
% colormap(jet())
% colorbar
% caxis([-20 25])