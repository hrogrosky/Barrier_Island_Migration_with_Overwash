function [H_overwash, PC1_overwash, PC2_overwash, PC3_overwash, PC4_overwash, DC, n1, surge_hist, EastSL, DC_top, DC_bottom, count_hist,OW_hist] = Overwash072324_v1(Hstar, PC,PC1, PC2, PC3, PC4, time,t, surge_hist,DC_top, DC_bottom, count_hist,window, window_top, OW_hist, k1, k2, SS_scale)

% We only used these files if we want to run the overwash code on it's own,
% not inside of the maincode. Mostly used for troubleshooting. Need to
% uncomment the island you want to use, line 11, and comment out line 18.

%filename='Smith03312021.mat'
%filename='Parramore03312021.mat';
%filename= 'testland2 copy.txt';

%H = importdata(filename);

%load in wave height data, covert to array
wv_data = readtable("wallops2022.csv",'VariableNamingRule','preserve');
wvht = wv_data(:,7);
wvht = table2array(wvht);

H=Hstar(:,:,1);
[n1, n2]=size(H(:,:));
EastSL=zeros(n1,1);
WSL=zeros(n1,1);
island_width = zeros(n1,1);
alpha = zeros(n1,1);
DC=zeros(n1,1);
dH= zeros(n1,1);
xDT=zeros(n1,1);  %%%% x coordinates of Dune Top
yDT=zeros(n1,1);  %%%% y coordinates of Dune Top
xDB=zeros(n1,1);  %%%% x coordinates of Dune Bottom
yDB=zeros(n1,1);  %%%% y coordinates of Dune Bottom
up=1;

%make copy of H for 'before' image
before_H = H;


%% To determine the eastern shore line (ESL)
for i=1:n1
    for j=n2:-1:1
        if H(i,j)>0
        %if H(i,j)>=0
            EastSL(i)=j; %%%%%%%%% Eastern shore line (ESL)%%%%%%%%%
            break
       end
    end
end


 %% To determine the western shore line (WSL)
for i=1:n1
    if EastSL(i)>0
    for j=EastSL(i):-1:1
        if H(i,j)<=0 && H(i-50,j)<=0 && H(i-100,j)<=0
            WSL(i)=j; %%%%%%%%% Western shore line (ESL)%%%%%%%%%
            island_width(i) = EastSL(i)-WSL(i);
            break
       end
    end
    end
end      



%% To determine the j-th column of first dune crest at DC(i) from Eastern Shore line 
n3 = 8;
 for i=1:n1
    if EastSL(i)>0
        shore=EastSL(i); %adjust this line to move in from ESL
        for j= shore 
            l=0;
            j1 = max(1,j-n3*l);
            j2 = max(1,j-n3*(l+1));
            j3 = max(1,j-n3*(l+2));
            avg1 = mean(H(i,j2:j1));
            avg2 = mean(H(i,j3:j2));
                while avg2 >= avg1 %&& avg2>=0
                l = l + 1;
                j1 = max(1,j-n3*l);
                j2 = max(1,j-n3*(l+1));
                j3 = max(1,j-n3*(l+2)); 
                avg1 = mean(H(i,j2:j1));
                avg2 = mean(H(i,j3:j2));
                if avg1 <= 0 && avg2 <= 0 %% sets the ESL to be DC if avg starts to go negative
                    DC(i) = EastSL(i);
                    dH(i) = H(i,EastSL(i));
                break
                end   
                end
        DC(i)=j1; %%%%%%%%%% DC=dune crest %%%%%%%%%%
        dH(i)=H(i,j1);%%%%%%% dH= dune heights %%%%%%%

        end
    end   
 end

%%% Calculating the average/min/max of DC height in a set window
DC_window = zeros(length(H),1);

 
for i = 1:(1+window) %%% for the first 5 rows
    if island_width(i) > 180
    DC_window(i) = mean(nonzeros(dH(i:i+window_top)));
    end
    if island_width(i) <= 180
    DC_window(i) = mean(nonzeros(dH(i:i+window)));
    end
end

for i = (2+window):(n1-window-1) %%% for everything in between
    if island_width(i) > 180
    DC_window(i) = mean(nonzeros(dH(i-window:i+window_top)));
    end
    if island_width(i) <= 180
    DC_window(i) = mean(nonzeros(dH(i-window:i+window)));
    end
end

for i = (n1-window):n1 %%% for the last 5 rows
    if island_width(i) > 180
    DC_window(i) = mean(nonzeros(dH(i-window_top:i)));
    end
    if island_width(i) <= 180
    DC_window(i) = mean(nonzeros(dH(i-window:i)));
    end
end


  %store DC at rows 500 and 2600 at each time step
top_transect = 1000;
bot_transect = 2600;    
DC_top(t+1,1) = t;
DC_top(t+1,2) = top_transect;
DC_top(t+1,3) = DC(top_transect);
DC_bottom(t+1,1) = t;
DC_bottom(t+1,2) = bot_transect;
DC_bottom(t+1,3) = DC(bot_transect);

 if t == 0
     EastSL_initial = EastSL;
 end
 if t ~= 0
     EastSL_final = EastSL;
 end

%% initialize left and right boundary (CURRENTLY NOT IN USE. NOT A BIG EFFECT)
% right_bnd = zeros(length(DC),1);
% left_bnd = zeros(length(DC),1);
% 
% if t ~= 0
% for i = 1:n1
%     if EastSL(i) > 0
% right_bnd(i) = EastSL(i) - 25;
% left_bnd(i) = EastSL(i) - 175;
%     end
% end
% 
% for i = 1:length(DC)
%     if DC(i) > right_bnd(i)
%         DC(i) = right_bnd(i);
%         dH(i) = H(i,DC(i));
%     end
%     if DC(i) < left_bnd(i)
%         DC(i) = left_bnd(i);
%         dH(i) = H(i,DC(i));
%     end
% end
% end

% %% To find the maximum height of dune crest from dH 
% %To find the x coordinate of dune that start from top of the island(Ytop)
% %To find the x coordinate of dune that ends at the bottom (Ybot)

% 
Ytop=find((DC>0),1,'first');
Ybot= find((DC>0),1,'last');
% 
% %%% (x,y) coordinates of Ytop: 
xYtop=Ytop;  
yYtop=EastSL(Ytop); %% y coordinate of ESL at x where the dune start from top of the island
XYTop= [xYtop yYtop]; 
% 
% %%% (x,y) coordinates of Ybot: 
% 
xYbot=Ybot;
yYbot=EastSL(Ybot);  %% y coordinate of ESL at x where the dune ends at the bottom of the island
XYBot= [xYbot yYbot];


%% Check to see if Overwash runs

n = numel(wvht);        % determines the number of elements of vector wvht
rng("shuffle");       %ensures the random numbers do not repeat everytime MATLAB restarts
i = randi(n);        % returns a random integer between 1 and n
SS = 10*wvht(i)*(SS_scale+((t/26)*0.01)); % Storm surge multiplied by 10 since depth is 0.1 and SS_scale
move = 10 ; %Distance that we are moving inland from ESL (Not using now)
% wdist=20; % Distance that we are moving to the water from ESL
Hcopy=H; % Copy of H

%PLANT LOSS DUE TO SALT WATER
loss_P1 = 0.08; % literature suggests salt tolerant but experiments differ
loss_P2 = 0.04; % relatively salt tolerent
loss_P3 = 0.16; % P3 is salt sensitive - this is Morella cerifera (woody shrub)
loss_P4 = 0.12; % not sure about this one

% Initializing vectors of storm surge history, slabs removed and deposited
surge_hist(t+1,1) = t;
surge_hist(t+1,2) = SS;
count = zeros(n1,1);
count2 = zeros(n1,1);

%defining the deposit start location, d, to the left of DC
%defining the removal start location, wdist, to the right of the ESL

% alpha = 0.15;
for i = 1:n1
    if island_width(i)>180
        alpha(i) = 0.3;
    elseif island_width(i) <= 180
        alpha(i) = 0.2;
    end
end
% 
% d = round(alpha*SS);

beta = 0.5;
wdist = round(beta*SS);

% k=0.5; %number of cells removed if DC height < SS
% k2 = 0.05; %number of cells removed if DC height > SS (0.01 for Smith, 0.05 for Parramore)


% The overwash routine runs if the randomly selected storm surge is greater
% than a specified value. Right now we use 20 (2m). Removing cells focuses
% on the entire length of the island (xYtop to xYbot) from the ESL + 20 to 
% (wdist = 20), all the way to the the location of the first dune
% crest - 5 (moving right to left). Depositing cells focuses on the entire 
% length of the island (xYtop to xYbot) from the location of the first dune
% crest - 5, all the way to max((DC(i)-5-count(i)),1). Previously the western
% boundary was the WSL, but we ran in to issues where we hadn't deposited
% enough sand before hitting the boundary- this allows us to move past it.
% Depositing stops when the amount removed = amount deposited. 

if SS >=  1
for i= xYtop:xYbot
if DC(i) ~= 0

    %   REMOVING CELLS if DC < SS
     for j = min((EastSL(i)+wdist),length(H)):-1:max((DC(i)-5),1)
        if H(i,j)<= SS && DC_window(i)<SS %if the elevation at each cell is less than SS

            % If there are plants present, we will calculate a random
            % number between 0 and 1. If that number is greater than the
            % plant covereage in that cell, we will remove k cells.

            if PC(i,j) ~= 0
               remove_rand = rand(1,1);
            if remove_rand > PC(i,j)
                H(i,j)= H(i,j) - k1;  
                count(i) = count(i) +1;
                % decrease percent cover of affected cells, with a min of 0
                PC1(i,j,1) = max(PC1(i,j,1)-loss_P1, 0);
                PC2(i,j,1) = max(PC2(i,j,1)-loss_P2, 0);
                PC3(i,j,1) = max(PC3(i,j,1)-loss_P3, 0);
                PC4(i,j,1) = max(PC4(i,j,1)-loss_P4, 0);
            end
            end

            % If there are no plants present we remove k cells
            if PC(i,j) == 0 
                H(i,j) = H(i,j) - k1; 
                count(i) = count(i) +1;
            end
            
        end
     end
 
 % REMVOING CELLS IF DC > SS
    for j = min((EastSL(i)),length(H)):-1:max(DC(i),1)
        if DC_window(i)>SS %if the elevation of DC window is greater than SS

            % If there are plants present, we will calculate a random
            % number between 0 and 1. If that number is greater than the
            % plant covereage in that cell, we will remove k cells.

            if PC(i,j) ~= 0
               remove_rand = rand(1,1);
            if remove_rand > PC(i,j)
                H(i,j)= H(i,j) - k2;  
                count(i) = count(i) +1;
                % decrease percent cover of affected cells, with a min of 0
                PC1(i,j,1) = max(PC1(i,j,1)-loss_P1, 0);
                PC2(i,j,1) = max(PC2(i,j,1)-loss_P2, 0);
                PC3(i,j,1) = max(PC3(i,j,1)-loss_P3, 0);
                PC4(i,j,1) = max(PC4(i,j,1)-loss_P4, 0);
            end
            end

            % If there are no plants present we remove k cells
            if PC(i,j) == 0 
                H(i,j) = H(i,j) - k2; 
                count(i) = count(i) +1;
            end
            
        end
     end

%   DEPOSITING CELLS
if count(i) >= 1

    OW_hist(i) = OW_hist(i) + 1; %keeping track of # of times overwash ran per row

       for  j= min(DC(i)-round(alpha(i)*SS),length(H))
        if j > 0
        if H(i,j) <= SS && DC_window(i)<SS 
        while count2(i) < count(i)
            j = max(j,1);
            % If there are plants present, we will calculate a random
            % number between 0 and 1. If that number is less than the
            % plant covereage in that cell, we will add 2k cells. If it is
            % greater than plant coverage, we add k cells. 

            if PC(i,j) ~= 0 && j>0
                   deposit_rand = rand(1,1);
                if deposit_rand < PC(i,j)
                H(i,j)=H(i,j)+ (2*k1);
                count2(i) = count2(i) + 2;
                j = min(j-1,length(H));
                % decrease percent cover of affected cells
                PC1(i,j,1) = max(PC1(i,j,1)-loss_P1, 0);
                PC2(i,j,1) = max(PC2(i,j,1)-loss_P2, 0);
                PC3(i,j,1) = max(PC3(i,j,1)-loss_P3, 0);
                PC4(i,j,1) = max(PC4(i,j,1)-loss_P4, 0);
                end
                if deposit_rand >= PC(i,j)
                H(i,j)=H(i,j)+ k1;
                count2(i) = count2(i) + 1;
                j = min(j-1,length(H));
                % decrease percent cover of affected cells
                PC1(i,j,1) = max(PC1(i,j,1)-loss_P1, 0);
                PC2(i,j,1) = max(PC2(i,j,1)-loss_P2, 0);
                PC3(i,j,1) = max(PC3(i,j,1)-loss_P3, 0);
                PC4(i,j,1) = max(PC4(i,j,1)-loss_P4, 0);
                end
            end

            %if there are no plants present, we dump k cells

            if PC(i,j) == 0  
              H(i,j)=H(i,j)+ k1;
              j = min(j-1,length(H));
              count2(i) = count2(i) + 1;  
            end

          
         end %end of while loop
        end
        end
       end
end 

end
end
end %end of overwash loop

 % this stores the number of cells removed in row 500 and 2600
            % at each time step
count_hist(t+1,1) = t;
count_hist(t+1,2) = count(500);
count_hist(t+1,3) = count(2600);

PC1_overwash = PC1(:,:,1);
PC2_overwash = PC2(:,:,1);
PC3_overwash = PC3(:,:,1);
PC4_overwash = PC4(:,:,1);

H_overwash = H;


end % end of function

       
% % % before and after figure %commented out for main code
% sz = 15;
% tiledlayout(1,2);
% nexttile
% imagesc(before_H)
% title('Before')
% colorbar
% clim([-5 26])
% nexttile
% imagesc(H)
% hold on
% scatter(DC, 1:n1, sz ,'w','filled')
% contour(H(:,:,1),[-0 -0],'color','w','linewidth',1.1)
% title('After')
% colorbar
% clim([-5 26])       



%just one figure
% climsH=[min(min(H(:,1))) max(max(H(:,1)))];
% figure
% imagesc(H, climsH)
% hold on
% scatter(EastSL,1:n1,'ro')
% scatter(WSL,1:n1,'gr')
% scatter(DC, 1:n1 ,'w', 'filled')
% hold on
% % contour(H,[0 0], 'color', 'm')
% legend('ESL', 'WSL')
% colormap(jet())
% colorbar
% caxis([-1 35])
% 
% 
% for i = 1:length(EastSL)
%     if EastSL(i)>0 && WSL(i)>0
%         island_width(i) = EastSL(i)-WSL(i);
%     end
% end
% 
% plot(EastSL)
% hold on
% plot(EastSL_smooth)
