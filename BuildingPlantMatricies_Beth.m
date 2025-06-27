%% Using Island Matrix to create 4 plant matricies
%Updated Parramore 7/7/2024

%load in the image for updated Parramore
Iname = 'Parramore2021';
data=imread('Parramore2021.png');      
data2=double(data(:,:,2));      

ddata=im2double(data); %convert image to a double array with all values between 0,1

% figure
% colormap jet
% imagesc(ddata)
% colorbar
Values=unique(ddata);

A=ddata;

%we need to assign an index value to each color on the jpeg:

NdexNum=10;      %number of colors on the jpeg, 
%                   remember to include white and black - might need to play with this a bit, 
%                   but you at least need to have a distinct number associated with the OUTLINES
%                   once those outlines are removed, the other areas should have more distinct coloring
%                   (if area not distinct after outlines removed, come back here and choose higher NdexNum)

[NdexMap,map6]=rgb2ind(A,NdexNum);  %NdexMap is your color-coded image
% it helps us out a lot to have even number of rows/cols, so we just delete one if necessary:
if 1==mod(size(NdexMap,1),2)
    NdexMap=NdexMap(1:size(NdexMap,1)-1,:);
end
if 1==mod(size(NdexMap,2),2)
    NdexMap=NdexMap(:,1:size(NdexMap,2)-1);
end

% figure(2)
% imagesc(NdexMap)
% title(sprintf('%s Island reduced to %d color indices',Iname,NdexNum))
% colorbar

%% OUTLINES %%%%
%removing the outlines (there are multiple colors associated with outlines: 1,2,0,7,8)
Problem=1; %index associated with outlines in the colormap
[n1,n2]=size(NdexMap);
for nn=1:5      %loop will cycle through five times to try any remove all of the instances of that color index
%     nn
    for ii=2:n1-1
        for jj=2:n2-1
            if NdexMap(ii,jj)==Problem      %if we find an instance of the problem color
                vec1=double([NdexMap(ii-1,jj) NdexMap(ii+1,jj) NdexMap(ii,jj-1) NdexMap(ii,jj+1)]); %take a vector of the 4-cell neighborhood (Nhd)
                vec1(vec1==Problem)=NaN;    %ignore elements of the Nhd that are same as problem
                NdexMap(ii,jj)=mode(vec1);  %replace the current cell with the most common value surrounding it
            end
        end
    end
end

Problem=2; %index associated with outlines in the colormap
for nn=1:5      %loop will cycle through five times to try any remove all of the instances of that color index
%     nn
    for ii=2:n1-1
        for jj=2:n2-1
            if NdexMap(ii,jj)==Problem      %if we find an instance of the problem color
                vec1=double([NdexMap(ii-1,jj) NdexMap(ii+1,jj) NdexMap(ii,jj-1) NdexMap(ii,jj+1)]); %take a vector of the 4-cell neighborhood (Nhd)
                vec1(vec1==Problem)=NaN;    %ignore elements of the Nhd that are same as problem
                NdexMap(ii,jj)=mode(vec1);  %replace the current cell with the most common value surrounding it
            end
        end
    end
end

Problem=0; %index associated with outlines in the colormap
for nn=1:5      %loop will cycle through five times to try any remove all of the instances of that color index
%     nn
    for ii=2:n1-1
        for jj=2:n2-1
            if NdexMap(ii,jj)==Problem      %if we find an instance of the problem color
                vec1=double([NdexMap(ii-1,jj) NdexMap(ii+1,jj) NdexMap(ii,jj-1) NdexMap(ii,jj+1)]); %take a vector of the 4-cell neighborhood (Nhd)
                vec1(vec1==Problem)=NaN;    %ignore elements of the Nhd that are same as problem
                NdexMap(ii,jj)=mode(vec1);  %replace the current cell with the most common value surrounding it
            end
        end
    end
end

Problem=5; %index associated with outlines in the colormap
for nn=1:5      %loop will cycle through five times to try any remove all of the instances of that color index
%     nn
    for ii=2:n1-1
        for jj=2:n2-1
            if NdexMap(ii,jj)==Problem      %if we find an instance of the problem color
                vec1=double([NdexMap(ii-1,jj) NdexMap(ii+1,jj) NdexMap(ii,jj-1) NdexMap(ii,jj+1)]); %take a vector of the 4-cell neighborhood (Nhd)
                vec1(vec1==Problem)=NaN;    %ignore elements of the Nhd that are same as problem
                NdexMap(ii,jj)=mode(vec1);  %replace the current cell with the most common value surrounding it
            end
        end
    end
end

Problem=8; %index associated with outlines in the colormap
for nn=1:5      %loop will cycle through five times to try any remove all of the instances of that color index
%     nn
    for ii=2:n1-1
        for jj=2:n2-1
            if NdexMap(ii,jj)==Problem      %if we find an instance of the problem color
                vec1=double([NdexMap(ii-1,jj) NdexMap(ii+1,jj) NdexMap(ii,jj-1) NdexMap(ii,jj+1)]); %take a vector of the 4-cell neighborhood (Nhd)
                vec1(vec1==Problem)=NaN;    %ignore elements of the Nhd that are same as problem
                NdexMap(ii,jj)=mode(vec1);  %replace the current cell with the most common value surrounding it
            end
        end
    end
end

% figure(3)
% imagesc(NdexMap)
% title(sprintf('%s Island reduced to %d color indices',Iname,NdexNum))
% colorbar
%% Creating 4 plant matrices and a beach matrix
%now that we have the colors indexed, we can create plant matrices based on
%the color. We also keep track of the beach area 

%initialize matricies 
P1 = zeros(size(NdexMap));
P2 = P1; 
P3 = P1;
P4 = P1;
Beach = P1;

[n1, n2] = size(NdexMap);

for i = 1:n1
    for j = 1:n2
        if NdexMap(i,j) == 9 %index 9 is P1 and P2
            P1(i,j) = 1;
            P2(i,j) = 1;
        end
         if NdexMap(i,j) == 4 %index 4 is P3
            P3(i,j) = 1;
         end
         if NdexMap(i,j) == 6 %index 6 is P4
            P4(i,j) = 1;
         end
        if NdexMap(i,j) == 3 %make water -1
            P1(i,j) = -1;
            P2(i,j) = -1;
            P3(i,j) = -1;
            P4(i,j) = -1;
        end
        if NdexMap(i,j) == 7 %identifies the beach area 
            Beach(i,j) = 1;
        end
    end
end

%plot initial plant matricies (no percent cover incorporated) 
% can comment out later 
% 
% figure(4)
% tiledlayout(2,2);
% 
% % Tile 1
% nexttile
% imagesc(P1)
% title('P1 - Ammophila')
% 
% % Tile 2
% nexttile
% imagesc(P2)
% title('P2 - Spartina Patens')
% 
% % Tile 3
% nexttile
% imagesc(P3)
% title('P3 - Morella')
% 
% % Tile 4
% nexttile
% imagesc(P4)
% title('P4 - Spartina Alterniflora')
%% Using the plant matrices to build the elevation matrix H
% If P1 is present, assign a random value from 1m to 5m (P1) or 0.75m to 3m (P2) 
% (not sure what to do about P2 since it is the same as P1?)
% If P3 is present, assign a random value from 1.5m to 2.5m
% If P4 is present, assign a random value from -0.5m to 1m
% If beach area, assign a random value from 0m to 1m

H = ones(n1,n2)*-20;

for i = 1:n1
    for j = 1:n2
        if P1(i,j) == 1
            H(i,j) = randi([5 15],1,1);
        end
        if P3(i,j) == 1
            H(i,j) = randi([15 25],1,1);
        end
        if P4(i,j) == 1
            H(i,j) = randi([-5 5],1,1);
        end
        if Beach(i,j) == 1
            H(i,j) = randi([0 5],1,1);
        end
    end
end

%% SMOOTHING
% This section smooths out the elevation matrix twice since heights were
% randomly selected and island will look very spotty

nnend=2; %(couple of times)

for nn=1:nnend

[n1, n2]=size(H(:,:));
EastSL=zeros(n1,1);
WestSL=zeros(n1,1);
H_smooth = H;

%%%% To determine the eastern shore line and western shore line
for i=1:n1
    for j=n2:-1:1
        if H(i,j)>0
            EastSL(i)=j; 
            break
       end
    end
end

for i=1:n1
    for j=1:n2
        if H(i,j)>0
            WestSL(i)=j; 
            break
       end
    end
end

% actual smoothing routine in moving 3x3 square 
 for i = 1:n1
    if EastSL(i)>0
        Eshore=EastSL(i); 
        Wshore=WestSL(i);
        for j= (Eshore):-1:(Wshore) 
            l=1;
            j1 = max(j - (l),1);
            j2 = j + (l);
            i1 = max(i - (l),1);
            i2 = i + (l);
            avg1 = round(sum(sum(H(i1:i2,j1:j2)))/9); %mean using original matrix 
%              mode1 = mode(mode(H(i1:i2,j1:j2)));
            H_smooth(i,j)= avg1; %only change center cell
            if j1 <= WestSL(i)
                break
            end
        end
    end
 end

 H = H_smooth;
end

figure(5)
colormap jet
imagesc(H)
title('Final Smoothed Island - water not addressed yet')
colorbar