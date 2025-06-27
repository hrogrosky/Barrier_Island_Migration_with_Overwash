%% Loading OBSERVED Parramore 1984 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%% CHANGE FILE NAME IN LINE 81 %%%%%%%%%%


I=imread('Parramore 1984.png');
[m1,m2]=size(I(:,:,1));
H_1984 = zeros(m1,m2);

% this loop isolates cells whose color is close to 0 (black), making
% those entries 1 and everything else 0
for i = 1:m1
    for j = 1:m2
        if I(i,j,1) <= 10
            H_1984(i,j) = 1;
        end
    end
end

%removing the border from the image (was black but we don't want it)
H_1984(4:8,:) = 0;
H_1984(3299:end,:) = 0;

%% Parramore 1998 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=imread('Parramore 1998.png');
[m1,m2]=size(I(:,:,1));
H_1998 = zeros(m1,m2);

% this loop isolates cells whose color is close to 0 (black), making
% those entries 1 and everything else 0
for i = 1:m1
    for j = 1:m2
        if I(i,j,1) <= 10
            H_1998(i,j) = 1;
        end
    end
end

%removing the border from the image (was black but we don't want it)
H_1998(2:6,:) = 0;
H_1998(:,2:7) = 0;
H_1998(:,2547:end) = 0;
H_1998(3294:3299,:) = 0;

%shift 82 columns to the right and down 60
H_1998 = imresize(H_1998,0.92);
H_1998(H_1998 ~= 0) = 1;
H_1998 = circshift(H_1998,[60 55]);


%% Parramore 2011 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=imread('Parramore 2011.png');
[m1,m2]=size(I(:,:,1));
H_2011 = zeros(m1,m2);

% this loop isolates cells whose color is close to 0 (black), making
% those entries 1 and everything else 0
for i = 1:m1
    for j = 1:m2
        if I(i,j,1) <= 10
            H_2011(i,j) = 1;
        end
    end
end

%removing the border from the image (was black but we don't want it)
H_2011(1:2,:) = 0;
H_2011(:,1:2) = 0;
H_2011(:,2547:end) = 0;
H_2011(3298:end,:) = 0;

%shifting 155 units up and 175 units to the left
H_2011 = imresize(H_2011,1.08);
H_2011(H_2011 ~= 0) = 1;
H_2011 = circshift(H_2011,[-155 -175]);


%% Loading MODEL Parramore year 0 (1984)
filenamestring='Pte702_03_Mar_2025_15_35_49'; % CHANGE THIS LINE EVERY RUN


yearnum=0; 
islandstring=filenamestring(1);
cd /Users/beththomas/Documents/'Documents - Beth MacBook Air'/MATLAB/'Barrier Island'
curdir=pwd;
cd(filenamestring);
filenamestringmat=strcat([filenamestring,'t',int2str(yearnum),'.mat']);
load(filenamestringmat);
cd ..

% resaving the H matrix so positive elevation is 1, anything else is 0
[n1,n2]=size(H);
H_Mod_yr0 = zeros(n1,n2);

for i = 1:n1
    for j = 1:n2
        if H(i,j) > 0
            H_Mod_yr0(i,j) = 1;
        end
    end
end


%scale by 1.1, shift 120 columns to the right (second #) and down 120
%(first #)
H_Mod_yr0 = imresize(H_Mod_yr0,1.1);
H_Mod_yr0(H_Mod_yr0 ~= 0) = 1;
H_Mod_yr0 = circshift(H_Mod_yr0,[-120 -120]);


Area_yr0 = sum(sum(H_Mod_yr0>0));

% plot the OBS 1984 and MOD 1984 contours to see how we line up
% contour(H_1984,[1 1],'Color','blue','LineWidth',2); 
% hold on
% contour(H_Mod_yr0,[1 1],'Color','black','LineWidth',2); 
% legend('Obs','Mod');
% set(gca, 'YDir','reverse')


%% Loading MODEL Parramore year 14 (1998)
yearnum=14; 
islandstring=filenamestring(1);
cd /Users/beththomas/Documents/'Documents - Beth MacBook Air'/MATLAB/'Barrier Island'
curdir=pwd;
cd(filenamestring);
filenamestringmat=strcat([filenamestring,'t',int2str(yearnum),'.mat']);
load(filenamestringmat);
cd ..

% resaving the H matrix so positive elevation is 1, anything else is 0
[n1,n2]=size(H);
H_Mod_yr14 = zeros(n1,n2);

for i = 1:n1
    for j = 1:n2
        if H(i,j) > 0
            H_Mod_yr14(i,j) = 1;
        end
    end
end


%scale by 1.1, shift 120 columns to the right and down 120
H_Mod_yr14 = imresize(H_Mod_yr14,1.1);
H_Mod_yr14(H_Mod_yr14 ~= 0) = 1;
H_Mod_yr14 = circshift(H_Mod_yr14,[-120 -120]);

Area_yr14 = sum(sum(H_Mod_yr14>0));

%% Loading MODEL Parramore year 27 (2011)
yearnum=27; 
islandstring=filenamestring(1);
cd /Users/beththomas/Documents/'Documents - Beth MacBook Air'/MATLAB/'Barrier Island'
curdir=pwd;
cd(filenamestring);
filenamestringmat=strcat([filenamestring,'t',int2str(yearnum),'.mat']);
load(filenamestringmat);
cd ..

% resaving the H matrix so positive elevation is 1, anything else is 0
[n1,n2]=size(H);
H_Mod_yr27 = zeros(n1,n2);

for i = 1:n1
    for j = 1:n2
        if H(i,j) > 0
            H_Mod_yr27(i,j) = 1;
        end
    end
end


%scale by 1.1, shift 120 columns to the right and down 120
H_Mod_yr27 = imresize(H_Mod_yr27,1.1);
H_Mod_yr27(H_Mod_yr27 ~= 0) = 1;
H_Mod_yr27 = circshift(H_Mod_yr27,[-120 -120]);

Area_yr27 = sum(sum(H_Mod_yr27>0));

%% Print Model Area Vector
% Plot all 3
% contour(H_Mod_yr0,[1 1],'Color','blue','LineWidth',2); 
% hold on
% contour(H_Mod_yr14,[1 1],'Color','green','LineWidth',2); 
% hold on
% contour(H_Mod_yr27,[1 1],'Color','red','LineWidth',2);
% set(gca, 'YDir','reverse')

Model_Area = [Area_yr0, Area_yr14, Area_yr27]

Area_prct = [(Area_yr0/sum(sum(H_1984)))*100, (Area_yr14/sum(sum(H_1998)))*100, (Area_yr27/sum(sum(H_2011)))*100]

%% CALCULATE OBSERVED AREA
sum_0 = 0;

for i = 1:3300
    for j = 1:2550
        if H_1984(i,j)>0 && H_Mod_yr0(i,j)>0
            sum_0 = sum_0 + 1;
        end
    end
end


sum_14 = 0;

for i = 1:3036
    for j = 1:2346
        if H_1998(i,j)>0 && H_Mod_yr14(i,j)>0
            sum_14 = sum_14 + 1;
        end
    end
end


sum_27 = 0;

for i = 1:3484
    for j = 1:2615
        if H_2011(i,j)>0 && H_Mod_yr27(i,j)>0
            sum_27 = sum_27 + 1;
        end
    end
end

%% Overlapping Contours
x_km = ([1:1:width(H)]*2);
y_km = ([1:1:length(H)]*2);

tlc = tiledlayout(1,3);
tlc.TileSpacing = 'tight';


%plotting 1984 contour
nexttile(1)
contour(H_Mod_yr0,[1,1],'Color','blue','LineWidth',2);
hold on
contour(H_1984,[1,1],'Color','red','LineWidth',2);
set(gca, 'YDir','reverse');
xticklabels(x_km);
yticklabels(y_km);
xlabel('x (km)');
ylabel('y (km)');
text(-100,3600,'(a)','FontSize',16);
title('Year 0 Contours');
ax=gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.XLabel.FontSize = 16;
ax.YLabel.FontSize = 16;
ax.Title.FontSize = 16;

%plotting 1998 contour
nexttile(2)
contour(H_Mod_yr14,[1,1],'Color','blue','LineWidth',2);
hold on
contour(H_1998,[1,1],'Color','red','LineWidth',2);
set(gca, 'YDir','reverse')
xticklabels(x_km);
set(gca,'YTickLabel',{' '})
text(-100,3600,'(b)','FontSize',16);
title('Year 14 Contours');
ax=gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.XLabel.FontSize = 16;
ax.YLabel.FontSize = 16;
ax.Title.FontSize = 16;

%plotting 2011 contour
nexttile(3)
contour(H_Mod_yr27,[1,1],'Color','blue','LineWidth',2);
hold on
contour(H_2011,[1,1],'Color','red','LineWidth',2);
legend('Model Contour', 'Observed Contour',fontsize=16);
set(gca, 'YDir','reverse')
xticklabels(x_km);
set(gca,'YTickLabel',{' '})
text(-100,3600,'(c)','FontSize',16);
title('Year 27 Contours');
ax=gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.XLabel.FontSize = 16;
ax.YLabel.FontSize = 16;
ax.Title.FontSize = 16;


title(tlc,'Parameterization',fontsize=24)

%%
overlap_model = [(sum_0/Area_yr0)*100 (sum_14/Area_yr14)*100 (sum_27/Area_yr27)*100]
overlap_obs = [(sum_0/sum(sum(H_1984)))*100 (sum_14/sum(sum(H_1998)))*100 (sum_27/sum(sum(H_2011)))*100]
