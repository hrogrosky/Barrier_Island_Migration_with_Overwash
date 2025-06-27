%trying to make this code easier to work with 05/19/2021
%inputs needed to enter: 
%       IslandSize (whatever estimate you got) 
%       Desired Beach slope (using 0.02 from Davidson Arnott 'Conceptual model for...'

%name the island you wish to work on, Matlab will load up the jpg with that name
Iname='SmithVEG2011';
% Iname='Smith3';
% Iname='Parramore'

loadprompt='Would you like to use a saved index map? (1 for yes, any other number for no) \n';
loadNdex=input(loadprompt);

if loadNdex~=1
data=imread('SmithVEG2011.png');      %data reads in as an 
data2=double(data(:,:,2));      

ddata=im2double(data); %convert image to a double array with all values between 0,1

figure
colormap jet
imagesc(ddata)
colorbar
Values=unique(ddata);

A=ddata;

%we need to assign an index value to each color on the jpeg:

NdexNum=8;      %number of colors on the jpeg, 
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

figure(2)
imagesc(NdexMap)
title(sprintf('%s Island reduced to %d color indices',Iname,NdexNum))
colorbar

% Now we need to remove the thick black lines that are in between te regions
%   you should have an output of the color-indexed image, 
%   the prompt will let you select the color index which appears to be 
%   associated with the outlines. 
%       if it isn't a distinct color, you may need to go back and try a higher
%       NdexNum value. Keep increasing until the outline is clear.

prompt1='What index appears to be associated with the island outlines? \n';
Problem=input(prompt1);
% Problem=7; %index associated with outlines in the colormap
disp(sprintf('Removing color %d',Problem));
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

figure(2)
imagesc(NdexMap)
title(sprintf('%s Island reduced to %d color indices',Iname,NdexNum-1))
colorbar

%often removing one color is not enough - depends on the island and the number of original indices used, if we need to remove another color index:
prompt2='What index now appears to be associated with the island outline? \n ****Remember: zoom in, look for numbers not clearly part of a section, or only existing between areas**** \n ****(if the image looks good already, enter -1)**** \n';
Problem=input(prompt2);
if Problem~=-1
% Problem=4; %index associated with outlines in the colormap
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


figure(2)
imagesc(NdexMap)
title(sprintf('%s Island reduced to %d color indices, index %d removed',Iname,NdexNum-2,Problem))
colorbar
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This next section will reduce the number of pixels in the image, which will make the image smaller and cut down computation time
%             this may throw off your overall island sizes, so it's not recommended unless you are very sure of what you are doing
%             OR if you are just playing around and trying to get the settings right
% 
% 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rescale the image; (I'm pretty sure we don't need to do this either, because we resize the image in the MainCode now
%         I'll leave it here just in case
[n1,n2]=size(NdexMap);
% IslandSize=9065000; %island surface area in square meters
% % ScalingFactor=IslandSize/(n1*n2);
% H2O=mode(mode(NdexMap));                %most common number on map will be water
% NonH2O=sum(sum(NdexMap~=H2O));          %just marsh and island area
% ScalingFactor=sqrt(IslandSize/NonH2O);  %scale estimated area of the island by current area of island in image
% 
% NdexMap_temp=imresize(NdexMap,ScalingFactor);
% NdexMap=NdexMap_temp;

% Replace the color 2 with 5  
% IND6(IND6==2)=5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% This section smooths out pixels that don't match the stuff around it in
% any of four directions - run it a couple times
nnend=2; %(couple of times)
dd=6; % (How many cells do you want to look at in each direction?)
for nn=1:nnend
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

figure(3)
imagesc(NdexMap)
title(sprintf('%s Island, outlines removed and smoothed',Iname))
colorbar

% The above section can take a while, if you get to a point where you like the index map you should save it:
SVprompt1='Do you want to save the Index Map? (1 for yes, any other number for no) \n';
Sv1=input(SVprompt1);
if Sv1==1
    svtxt=strcat([Iname,'IndexMap']);
    save(svtxt);
end

%then you can pick up the code from here later on by loading what you had previously:
elseif loadNdex==1
filename=strcat([Iname,'IndexMap']);
NdexStruct=load(filename);
NdexMap=NdexStruct.NdexMap;
figure(3)
imagesc(NdexMap)
title(sprintf('%s Island, outlines removed and smoothed',Iname))
colorbar
X=NdexMap;
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This first section here re-orients the island to be roughly vertical, (shoreline parallel to y-axis)
%   this can be helpful sometimes when trying to get the image "just right" since column procedures go right to left,
%   but it's better to avoid it if you can, because you have to create a lot of buffer space (extra rows/columns) to rotate, 
%   and it's tricky to get the island back into the correct position.


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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
X=NdexMap; %making a copy to play with so we don't lose original
H=zeros(size(X,1),size(X,2));   %hopefully will be the new elevation matrix

% we need to know where the water/ocean is (i.e. non-land, non-marsh):
%       We are going to make the H map all 10s where there is island, and 0 everywhere else
%       This makes it easier to find the outlines of the island

H2O=mode(mode(NdexMap));                %most common number on map will be water
for ii=1:size(X,1)
    for jj=1:size(X,2)
        if X(ii,jj)==H2O %all water area (white)
            H(ii,jj)=0; %make the water 0 in the NEW map (this line is actually redundant, since H is all 0 at this point)
        else  % make subaerial plus marsh area (actual island) 10 for now
            H(ii,jj)=10;
        end
    end
end

Hnew=H; %again, a copy in case we mess something up...

WE_shore=zeros(size(Hnew,1),2);%a n1x2 vector with the West shoreline col number in position 1, and East shoreline col number in position 2
for ii=1:size(Hnew,1)
    RowNow=Hnew(ii,:);
    check1=sum(RowNow>0);
    if check1~=0
        WE_shore(ii,1)=find(RowNow>0,1,'first');
        WE_shore(ii,2)=find(RowNow>0,1,'last');
    end
end

MaxWest=(min(WE_shore(WE_shore~=0))); %westernmost point of island
MaxEast=(max(WE_shore(WE_shore~=0))); %easternmost point of island
NS_shore=zeros(size(Hnew,2),2); %another vector, n2x2 with the NORTH shoreline row number is position 1 and South shoreline row number in position 2
for jj= 1:size(Hnew,2)
    ColNow=Hnew(:,jj);
    check2=sum(ColNow>0);
    if check2 ~=0
        NS_shore(jj,1)=find(ColNow>0,1,'first');
        NS_shore(jj,2)=find(ColNow>0,1,'last');
    end
end
%this just pulls apart the above vectors to make the next steps a little easier:

Wv=WE_shore(:,1); Ev=WE_shore(:,2); %vectors of col index for west/east shore
Nv=NS_shore(:,1); Sv=NS_shore(:,2); %vectos of row numbers for north/south shore
W0=min(Wv(Wv>0));E0=max(Ev(Ev>0)); %min max col cell(western-most and eastern-most)
N0=min(Nv(Nv>0));S0=max(Sv(Sv>0)); %max minrow cell (northern-most and southern-most)


prompt3='For the ocean slope, what is change in y (height)? \n ****(I used 0.1, or one cell height, for my thesis run maps)**** \n';
DDrop=input(prompt3);
% DDrop=0.1;  %elevation decrease (change in height)
prompt4='For the ocean slope, what is change in x (distance) going east? \n ****This should be smaller than the change in distance for the marsh (going west)**** \n ****(I used 100 for my thesis run maps)**** \n';
PerLengthEast=input(prompt4);
% PerLengthEast=100; %nomber of cells to skip before decreasing again (change in distance)
prompt5='For the marsh slope, what is change in x (distance) going west? \n ****This should be larger than the change in distance for the ocean (going east)**** \n ****(I used 500 for my thesis run maps)**** \n';
PerLengthWest=input(prompt5);
% PerLengthWest=500;
%putting in terms of dy/dx for clarity:
dy=DDrop;
dxE=PerLengthEast;
dxW=PerLengthWest;
dxNS=floor(.5*(dxE+dxW));   %North/South gradient is just the average of the east/west gradient - this was for simplicity and easy smoothing later

Hcopy=H;%keeping a copy of the land/water areas 
%               We want H to be our final elevation map, so we adjust the ocean depths in that array
%               You will need to make sure neither the west/east shorelines are in the first/last index 
%                   if they are, you will need to add some extra columns of zeros at beginning/end of map

for n=1:length(WE_shore)
    RowNow=H(n,:);      %I usually work one row at a time, MATLAB likes to work with vectors more than arrays
    west=WE_shore(n,1); %west is the index of western shoreline
    east=WE_shore(n,2); %east is the index of eastern shoreline
%     filling out the western underwater area:
    if west~=0 && west~=1 %if west shoreline exists and west shoreline isn't in the first column...
        RowNow(west-1)=-1;
        for kk=-2:-1:-west+1
            RowNow(west+kk)=-1+RowNow(west+kk)+dy*(kk/dxW); %depth is applied linearly using change in height/distance inputs 
%                                                                 (starts at -1 to ensure it stays negative, we will smooth the whole map later)
        end
    end
%     filling out the eastern underwater area
    if east~=0 &&east~=size(H,2)% if east shoreline exists and east shoreline isn't in the last column...
        RowNow(east+1)=-1;
        for kk=1:(length(RowNow)-east)
            RowNow(east+kk)=-1+RowNow(east+kk)+dy*(-kk/dxE);%depth is applied linearly using change in height/distance inputs 
%                                                                 (starts at -1 to ensure it stays negative, we will smooth the whole map later)
        end
    end
    H(n,:)=RowNow;  %replacing the vector we worked on into the corresponding row of array H
end

%the area before the island (to the north) and the area after the island (to the south) need to be done separately:
    %copy first and last row to get gradient change from above, replace the pos. elev.
    %cells with -0.1 and then use the same process as above on the whole column
    %going up and down from the north and south tip
    
RstarN=find((WE_shore>0),1,'first');%index of first row that has the island in it
RowStarN=H(RstarN,:);               %copy of that first row
DelInd=find(RowStarN>=0);%DelInd=index is where we need to replace the positive elevation with -0.1
Mn=mean(DelInd);%i don't think this is used for anything...
RowStarN(DelInd)=-1; %is negative .1 meters, as stated above
% now moving back (or north, rather) from the row we just replaced and applying the water gradient:
for n=-1:-1:-RstarN+1       
    H(RstarN+n,:)=RowStarN(:)+dy*(n/dxNS);
end
% process repeats to go south from bottom tip of island:
RstarS=find((WE_shore(:,1)>0),1,'last');
RowStarS=H(RstarS,:);
DelInd=find(RowStarS>=0);
RowStarS(DelInd)=-1;
for n=1:length(WE_shore)-RstarS
    H(RstarS+n,:)=RowStarS(:)+dy*(-n/dxNS);
end

figure
imagesc(H)

% now we need to replace all of the subaerial/ positive elevation area of island with reasonable elevations
% this is based on the elevations that we know are associated with each height--
figure(3)
ElevPromptS='Using the Elevation Map, which color index appears to be associated with the marsh?';
answer1=inputdlg(ElevPromptS);
MarshIndex=str2num(answer1{1});
ElevPromptG='Using the Elevation Map, which color index appears to be associated with the lowland/grasses? (if there are multiple, list all of them WITHOUT commas)';
answer2=inputdlg(ElevPromptG);
LowGrassIndex=str2num(answer2{1});
ElevPromptM='Using the Elevation Map, which color index appears to be associated with Morella?';
answer3=inputdlg(ElevPromptM);
MorellaIndex=str2num(answer3{1});

%Helev will be our elevation matrix while we are building it. 
%       We use the copy of the index matrix X to verify what zone we are in
%       then apply some (randomly chosen) elevation values from a range that makes sense
%       for that zone (low lying/grass plant zone, morella zone, marsh zone)
Helev=zeros(size(H,1),size(H,2)); 
nearshore=-1:0.1:-.01;
marsh=-.5:.1:.5;
LowE=1:.1:1.5;
HighE=2:.1:2.8; %maxing out a meter higher than the morella upper bound - smoothing will
%                     bring this value down close to the real upper bound

% while we are adding elevations, we also want to add the beach to the eastern shore of island:

Slope=.02; %desired beach slope (Davidson Arnott 'Conceptual model on the effects of sea level rise...'
BchMax=1.5; %dy (max beach height - been using 1.5 for a while now)
BchWidth=BchMax/Slope;

for ii=1:size(Helev,1)
    for jj=1:size(Helev,2)
        GrassCheck=ismember(LowGrassIndex,X(ii,jj)); %(checking if the current color index of X(ii,jj) is in one of the areas associated with the grasses) 
        MorellaCheck=ismember(MorellaIndex,X(ii,jj)); %(same check, but this time testing if it's one of the indices associated with morella)
        if X(ii,jj)==MarshIndex %(if X(ii,jj) is the marsh index - marsh usually only has one index associated with it, so we don't need the above check)
            Helev(ii,jj)=marsh(randi([1,numel(marsh)])); %assign a random marsh elev. value
        elseif sum(GrassCheck>0) %if X(ii,jj) is one of the grass indices...
            for mm=1:length(LowGrassIndex) %we loop through all of the grass indices
                inputnow=LowGrassIndex(mm); %establishing which index we are looking at
                    if X(ii,jj)==inputnow %if X(ii,jj) is the established index...
                        Helev(ii,jj)=LowE(randi([1,numel(LowE)])); %give it a randomly chosen elevation from the grass range
                    end
            end
        elseif sum(MorellaCheck>0) %Works same as with grasses, but this time checks for indices associated with morella
            for mm=1:length(MorellaIndex)
                inputnow=MorellaIndex(mm);
                    if X(ii,jj)==inputnow
                        Helev(ii,jj)=HighE(randi([1,numel(HighE)])); %and assigns a random value from morella range
                    end
            end
        end
        if H(ii,jj)==10 %now that we have the elevations that we want for this cell we can 0 out the place-holder in H
            H(ii,jj)=0; %(if we don't so this the H map ends up +10 meters too tall
        end
        west=WE_shore(ii,1); %now we need to build the beach
        east=WE_shore(ii,2);
        Iwidth=east-west;
        if WE_shore(ii,1)~=0 && WE_shore(ii,2)~=0 %if we have shorelines
            for jj = 0:BchWidth %we build the beach up linearly to the max height
                Helev(ii,east-jj)=0+BchMax*(jj/BchWidth); %start with 0 bc we define shoreline as first cell from the east that is >=0
                H(ii,east-jj)=0; %getting rid of place-holders in H map
            end
        end
    end
end

%we did all the elvations in terms of meters, so we want to change the entire map to slabs (10 slabs per meter)
%we also want to have slabs in integer values so we round to whole numbers
% H=round((H+Helev)*10,0);
H=(H+Helev)*10;
figure
imagesc(H)
colorbar
title(sprintf('Elevation Matrix H for %s Island',Iname))
% save([Iname '_elevation'],'H');
%        
%%  Polishing off the water gradient with a Gaussian filter
% right now we have the basic island we want, but it is likely very rough looking beause of the random numbers
% and the north-south/east-west creation of the underwater areas. We want all of this to look as
% smooth and natural as possible, so we use the gaussian filter to "blur" the image

% The underwater areas need much more smoothing to look natural, and we don't want to smooth the
% island areas too much or we lose the lowland/morella regions that we just set. 
% (Smoothing the marsh as though underwater appears to work fine)
% To do this: We split map into a "negative" (less than marsh max 1m) and positive (>1m) and
%     use different sigmas (standard deviations) for the gauss filter 
%   (larger for smoother/underwater, smaller for coarse/subaerial island)

[rowN colN]=find(H<1);
[rowP colP]=find(H>=1);
N=zeros(size(H,1),size(H,2)); %for "negative" (<1) parts of H
P=zeros(size(H,1),size(H,2)); %for "positive" (>1) parts of H
for i=1:length(rowN)
    R=rowN(i);C=colN(i);
    N(R,C)=H(R,C);
end
for i=1:length(rowP)
    R=rowP(i); C=colP(i);
    P(R,C)=H(R,C);
end
Nfilt=imgaussfilt(N,50); %large filter=smoother underwater
Pfilt=imgaussfilt(P,3); %small filter = rougher subaerial island features
Hfilt=Nfilt+Pfilt;  %put the two back together





figure
imagesc(Hfilt)
title(sprintf('Smoothed Elevation map for %s Island with water gradient and beach',Iname))
colorbar

%smoothing has created non-integer values, so we need to round one last time:
FinalH=round(Hfilt);

FinalSavePrompt='Would you like to save the final elevation map? \n (1 for yes, and other number for no) \n';
SvF=input(FinalSavePrompt);
if SvF==1
    SvFtxt=strcat([Iname,'ElevMap']);
    save(SvFtxt,'FinalH');
end
    

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%previous versions of this used this local neighborhood function from AVALANCHE
% while I don't think it is used in this version, I will leave it here in case
% it is useful for future edits.

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
