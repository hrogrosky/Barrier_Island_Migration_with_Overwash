%% MainCode092820
MPswitch=1;
% ATcntr=0;
% SLRswitch=0;
% is the updated main code for barrier island evolution.
% this version of the code is build upon the original code by Greg Robson

% off-hand key notes (NOT EVEN REMOTELY COMPREHENSIVE)
%   **ARRAYS
%       H - main elevation matrix
%           (only updated in initialization and at end of the main loop)
%       Hstar - dummy elevation array for changes in subroutines
%           *third dimension unused - was for land categorization
%       P_i~Plant perent coverage matrices,
%       P3d - tracks dead morella
%           ~Pi(:,:,1)  current cell percent cover
%           ~Pi(:,:,2)  initial elevation of current cell
%               (currently unused, will be needed to check for burial)
%       PC - Combine Pi's for a total percent cover array - (no cell >1)
%             (Pi stores value -999 in cells without that plant - i.e.
%             water, so PC is the array with only cells greater than 0
%       Ptot= Sun of all PCi
%       W - water table data (currently unused)
%       S - salinity (?) data (currently unused)
%   **General PARAMETERS
%       time = number of iterations (two week time steps)
%       delta - slab heigh (meters)
%       L - slab length and width (meters)

%   **Routines and Governing Parameters
%     Main Code
%       t - current time step
%     (2)AeolianTransport
%       Pe - prob. of erosion
%       Pd - prob. deposition
%       n - number of slabs which can move due to wind
%    (3) Plant Propogation
%       currently runs twice per year
%       first time ignores morella and is meant to simulate springtime
%           growth for grasses (P1,P2,P4 - allow them more chance to grow)
%       PiIC - Plant initial percent cover for P1,P2,P3,P4
%       PiErosionCoefficient - factors into AT
%       P-Burial - number of slabs until death

%    (4) Marine Processes
%       currently runs once every 3 months
%
%     avalanche:
%           whole domain (underwater and island) ever 5th year time step
%           once per year runs THREE times per step over only subaerial (ground)
%               after SwampProcess, AeolianTransport, and before loop end
%               otherwise runs once per time step, before loop end
%
%   This version contains NO use of Land Categorization
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                    %%%%%LOAD ELEVATION MATRIX%%%%%




%%%%%% PARAMETERS YOU MAY NEED TO CHANGE %%%%%%
%      restartrun=0 --> use image file as IC; 
%      restartrun=1 --> use previous simulation result as IC
restartrun=0;
%      Code automatically saves images as .jpg.  To save images as .fig as
%      well, set saveasdotfit=1
saveasdotfig=0;
%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%

Sv1=[]; Sv1f=[]; Sv1Cent=[]; Sv2=[]; Sv2f=[]; Sv2Cent=[]; Sv3=[]; Sv3f=[]; Sv3Cent=[]; 






%%%%%%%%%% FOR INITIAL CONDITION FROM AN IMAGE %%%%%%%%%%%
if restartrun==0
    t0=0; 
    
    %%%%%% NEW PARAMETERS YOU WILL WANT TO CHANGE %%%%%%
    %%% WHICH ISLAND %%%
    islandstring='P';     % Use 'P' for Parramore; 'S' for Smith
    
    %%% WHICH VIMS SLR ACCELERATION CURVE %%%
    SLRrun=2;  %1 - historic; 2 - low; 3 - high; 4 - highest.  Based on SLR graph from VIMS
    
    %%% BRUUN RULE COEFFICIENT
    BruunParam=2000;  % Suggested values:  500, 1000, 1500, 2000, 2500 %was 4000 on 10/3/23
    
    %%% DIRECTORY WHERE MAINCODE FILE IS KEPT
    cd  '/Users/beththomas/Documents/Documents - Beth MacBook Air/MATLAB/Barrier Island';
    
    %%% DETERMINE HOW OFTEN TO PLOT AND/OR SAVE DATA %%%
    DurImFreq=26*27;       %how often to output "during" elevation images (26~yearly, 130~every 5 years)
        MElevCont=1; %make 1 if you want contour
        SElevCont=0;
    PlantDurImFreq=26*27; %same, just for plant images %WAS 130
    DurCUImFreq=130;
    DursavedataFreq=130; %how often to save data (in 2-week intervals; 26~yearly, 78~every 3 years, 130~every 5 years)
    Dursavedatavec = [0, 26*14, 26*27];

    %%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%
    
    curdir=pwd;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if islandstring=='P'
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % PARRAMORE ISLAND DATA:%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % IslandArea=14400000;
        % 12km x 1200 m = 14400000 - Relative Influence Antecendent... paper, 12.4 m/yr
%         filename='ParramoreElevMap05252021.mat'; %reduced island
        filename='ParramoreElevMap07082024.mat'; %including all tiny islands
        IslandArea=19070000; %for parramore
        data=importdata(filename);
        H=data;                         %H will be the main elevation matrix
        TranChk=1;
    elseif islandstring=='S'
        %%%%%%%%%%%%%%%%%%%%%
        % SMITH ISLAND DATA:%
        %%%%%%%%%%%%%%%%%%%%%
        filename='SmithNew06132021.mat'
        IslandArea=9065000; %for Smith
        data=importdata(filename);
        H=data;                         %H will be the main elevation matrix
        TranChk=2;
    end


    % filename='EricsIsland.txt'
    % IslandArea=41000;%EricsIsland


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    ScaleFactor=floor(sqrt(IslandArea/(sum(sum(H(:,:)>0)))));
    % InterpFactor=.75*sqrt(ScaleFactor);   %ScaleFactor 2
    InterpFactor=0.5*sqrt(ScaleFactor);     %ScaleFactor 4

    H=interp2(H,InterpFactor);
    ScaleFactor=floor(sqrt(IslandArea/(sum(sum(H(:,:)>0)))));
    if TranChk==1
        TranRow=floor(.25*size(H,1));   %choose row for transect for parramore
    elseif TranChk==2
        TranRow=floor(0.75*size(H,1));    %choose row for transect for smith
    elseif exists(TranChk)==0
        TranRow=floor(0.5*size(H,1));
    end


    %%%%%% Reducing the size of updated Parramore (7/11/24) %%%%%%%%
    %%%% H was 6755 x 4755 , this should reduce both dim by 2 %%%%%%
 if TranChk == 1
        A2=zeros(floor(size(H)/2));
        [m2,n2]=size(A2);
        for ii=1:m2
        for jj=1:n2
            A2(ii,jj)=mode(H(2*ii-1:2*ii,2*jj-1:2*jj),'all');
        end
        end
end

H = A2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%END - LOAD ELEVATION MATRIX%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                    %%%%%    PARAMETERS     %%%%%                                                                                                                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                                                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n')
    fprintf('Initializing...')
               %cheat sheet:
    time= 26*27; %26*10; %8*26+702;  %2600~100yrs;%1300~50yrs;%780~30yrs;%520~20yrs;260~10yrs;
    delta=0.1;
    % L=1;%<----------ScaleFactor
    L=ScaleFactor;
    BchMax =1.5;   %greatest elevation that the beach spreads inland to
    OW=0; %switch for overwashing - might use when we get storms incorporated (if windspeed>16m/s OW=1 for on, etc)
    BchW=10;    %unused - maximum beach width
    MPswitch=1;     %changes marine processes:    --(0 or any ~=1,2,3) for OFF     
    %                                             --(1) for ON {new version - update using same equil. slope}       MarineProcesses03312021
    %                                             --(2)for ON {old version - update equil. slope every 12 wks}      OLDMarineProcesses03312021
    surge_hist = zeros(time+1,2); %initiates matrix to store time step and storm surge
    error_vec = zeros(time+1,4);  %initiates matrix to store cells<0, =0, and >0 at each time step
    DepositSupply=3;     %number of slabs of sediment supplied to beach each time MarinePRocesses is run (3 months
    SLRswitch=0;    %turns sea level rise ON(1)/OFF(0 or any ~=1)
    SLRyrs=[ 8 7 6 5]; %years at which to perform SLR at rate SLR and/or migration...
    Ma= 1; %1.0075;       %shoreline migration growth factor, use 1 for no accel/noSLR
    MigAccel=0;     %will need to track magnitude of migration change 
    %     MigYr=15; %initial yearly migration of the island in meters/year
    ScaleFactor=floor(sqrt(IslandArea/(sum(sum(H(:,:,1)>0)))));

    % QUADRATIC SEA LEVEL RISE (ADDED 4/29, RO)
    SLRaccelratevec=[0.0000 0.00014 0.00034 0.00051]*12/39';         % accel rate (m/yr^2), tuned to match red, yellow, green curves of VIMS plot
    SLRaccelrate=SLRaccelratevec(SLRrun);
    SLRinitrate=0.175/39;                                    %linear SLR (m/yr), tuned to match historic curve of VIMS plot, %0.175/39.3701;   
    yrmax=108; % 2100-1992
    SLRvec=SLRinitrate*(0:yrmax)+SLRaccelrate*(0:yrmax).^2; % Assume quadratic growth from 1992-2100
    SLRvec=[(-8:-1)*SLRinitrate SLRvec];
    % For 50-year run, starting at 2020, ending at 2070:  (index=1, 1984; index=37, 2020)
    SLRvec=SLRvec(37:end);
    MigYr=BruunParam*(SLRvec(2)-SLRvec(1));
    
    % R=zeros(size(H,1),1);
    MaxSwampWidth=1000; %how far the swamp should go out into the water - just have to make something up for now

    WindData=1; %0~steady windspeed/direction ,1~windspeed and direction from empirical data

    if WindData==0
        Windspeed=16;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        WindDir=1;  %wind direction            %  1=N  2=NE  3=E  4=SE  5=S  6=SW  7=W  8=NW  %                                                                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    Windmin=6;
    StormThreshold=16;
    ATcntr=0; %testing variable - to count number of times AT is called

    %%%%%%%%%%%%%%%%%%%%%%%%%   PLANT PARAMETERS    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Plant initial conditions:
    PlantRangeArray=[1 5;0.75 3;1.5 5.0;-0.5 1]; %all of the elevation ranges  for p1-p4
    P1IC=0.8;% 0.9; %was 0.8
    P2IC=0.8;% 0.9; %was 0.8
    P3IC=0.8;% 0.9; %was 0.8
    P4IC=0.8;% 0.9; %was 0.8
    P1PctMax=0.6; %was 0.6               %largest percentage we will allow any plant population on a given cell to attain
    P2PctMax=0.6; %was 0.6
    P3PctMax=1; %was 1
    P4PctMax=0.6; %was 0.6
    PctMax=[P1PctMax P2PctMax P3PctMax P4PctMax];
    MasterMax=1;  %was 1 %The most any cell can permit - 80% plant coverage
    KillSwitch=0;   %this will kill plants on the bottom half of the island if set to 1, set to 0 (or anything else) to turn off
    alpha=.01; %propagation rate for each populated cell
    DBE=0.03;     %death by elevation rate for each populated cell outside of plant's elevation range
    P3Killzone=200;
    gdrange1=[-.01:.01:.04]; %range of percent values for growth/death for plant populations at (0, 50)% cover
    gdrange2=[-.01:.01:.04];%[-.04:.01:.04]; %range of percent values for growth/death for plant pops greater than 50% cover


    %%%% CREATE FILENAME FOR THIS RUN USING PARAMETER VALUES ABOVE, ADDED 4/29, RO
    filenamestring=strcat([islandstring,'t0',int2str(t0),'te',num2str(time),'BP',int2str(BruunParam),'SLRi',num2str(SLRinitrate),'SLRa',num2str(SLRaccelrate),'MMax',num2str(MasterMax)]);
    filenamestring=strrep(filenamestring,'.','p');
    filenamestring=strrep(filenamestring,'-','n');

    % CREATE DIRECTORY FOR THIS RUN, ALL RESULTS STORED HERE, ADDED 4/29, RO
    mkdir(filenamestring);

    %Elevation Matrix dependent parameters:
    ROW=size(H,1);
    COLUMN=size(H,2);
    H(:,:,2)=zeros(ROW,COLUMN);     %Creating extra dimension for H
    Hstar=H;        %initializing dummy matrix
    MigCnt=zeros(ROW,1);   %used to store how many meters the shoreline should have receded by, shoreline moves when MigCnt>ScalingFactor
    M=max(max(H(:,:,1)));
    m=min(min(H(:,:,1)));
    W=zeros(ROW,COLUMN);  %do not comment out - needed as input
    S=zeros(ROW,COLUMN);  %do not comment out - needed as input

    %parameters unused so far (check with Greg)
    % numslabg=0;
    % q=0;
    % Salt=1;
    % ps=0.6;
    AV_ARRAY=zeros(time,1000,3);        %tracks cell movement in AV. routine
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Imaging Options %%%%%

    %  ***note that we do not currently have a running land categorization
    %     subroutine, so leave those options commented or at 0 - I think this
    %     will eventually come back into use so I am leaving the lines there***


    ElevImCU_during=0;  %close up images about the centroid - Close up elev, close up elev/contour split box, ALSO quadbox plant images
    Contour30=1;      %contours at 0, 13, and 27 years
    DurContour=0;       %images of the island outline

    TransectIm=0;       %images of 2D island transect
        TranImFreq=260; %every 9 years  
        if TransectIm==1 && exist('TranRow') == 0
            TranRow=floor(.5*size(H,1));   %choose row for transect
        end

    MnPlantCvrIm=1;      %mean percent cover for each plant, displayed at end of routine 1 to turn on, any other number to turn off
    
    ElevIm_before=1;        ElevIm_during=1;        ElevIm_after=0;
    LandCatIm_before=0;     LandCatIm_during=0;     LandCatIm_after=0;          %leave these at 0 until we have a working LandCategorization routine
    PlantPCIm_before=0;     PlantPCIm_during=1;     PlantPCIm_after=0;
    SingleStepICIm=0;   
    
    climsH=[min(min(H(:,:,1))) max(max(H(:,:,1)))];
    climsP=[-1 1];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%              %%%%%%Initializing Plant Arrays%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %GRASS #1 (Ammophila)(burial resistant)
    P1=zeros(ROW,COLUMN,2);
    P1Burial=2;         %meters of burial until death
    P1ErosionCoefficient= (2/3);     %used in formula   d=c1P1+c2P2+c3P3+c4P4

    %GRASS #2  (Spartina)
    P2=zeros(ROW,COLUMN,2);
    P2Burial=0.5;         %meters of burial until death
    P2ErosionCoefficient= (1/3);

    %SHRUB #1  (Morella)(bird-dispersed seeds)
    P3=zeros(ROW,COLUMN,2);
    P3Burial=1.5;         %meters of burial until death
    P3ErosionCoefficient= 1;

    %dead morella will track when morella reaches death and store dead
    %debris data in P3d(:,:,1) for some number of years which will be counted in
    %P3d(:,:,2) - this is only updated inside of PlantProp
    P3d=zeros(ROW,COLUMN,2);
    P3dErosionCoefficient=1;

    %SECOND TYPE OF SPARTINA!!!
    P4=zeros(ROW,COLUMN,2);
    P4Burial=0.5;
    P4ErosionCoefficient= (1/3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INITIAL PLANT CONDITIONS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % Each plant species has a preferred "altitude" that it grows best in
    % % % (min and max among all species is 0.75m and 5m, respectively).
    % % % We set each cell of "appropriate" elevation to be covered by some% of the
    % % % respective plant species. For cells which are submerged under water, we
    % % % set the proportion of plant coverage to -1.
    % % %
    % % % We have deleted the initialization of plant population since we are using
    % % % real data.
    % % %It is still important to ensure that the 2nd layer of each
    % % % plant matrix is initiated.
    fprintf('\n')
    fprintf('Seeding the island...')
    PC1=zeros(size(H,1),size(H,2));
    PC2=zeros(size(H,1),size(H,2));
    PC3=zeros(size(H,1),size(H,2));
    PC4=zeros(size(H,1),size(H,2));
    [Hstar,MeanBeachWidth,ESL,WSL,MESL,MWSL,OL]=Shoreline03312021(Hstar,delta,L,BchMax);%function to return shoreline(s)
    for i=1:size(H,1)
        for j=1:size(H,2)
            if (delta*H(i,j,1))>=PlantRangeArray(1,1) && (delta*H(i,j,1))<=PlantRangeArray(1,2)
                P1(i,j,1)=P1IC;
                P1(i,j,2)=H(i,j,1);
            elseif (H(i,j,1)-W(i,j))<0
                P1(i,j,1)=-999;
                P1(i,j,2)=0;
            end
            PC1(i,j)=P1(i,j,1)*(P1(i,j,1)>0);
            if (delta*H(i,j,1))>=PlantRangeArray(2,1) && (delta*H(i,j,1))<=PlantRangeArray(2,2)
                P2(i,j,1)=P2IC;
                P2(i,j,2)=H(i,j,1);
            elseif (H(i,j,1)-W(i,j))<0
                P2(i,j,1)=-999;
                P2(i,j,2)=0;%
            end
            PC2(i,j)=P2(i,j,1)*(P2(i,j,1)>0);
            if (delta*H(i,j,1))>=PlantRangeArray(3,1) && (delta*H(i,j,1))<=PlantRangeArray(3,2)
                P3(i,j,1)=P3IC;
                P3(i,j,2)=H(i,j,1);
            elseif (H(i,j,1)-W(i,j))<0
                P3(i,j,1)=-999;
                P3(i,j,2)=0;
            end
            PC3(i,j)=P3(i,j,1)*(P3(i,j,1)>0);
            if (delta*H(i,j,1))>=PlantRangeArray(4,1) && (delta*H(i,j,1))<=PlantRangeArray(4,2)
                if ESL(i)~=0
                    if MWSL(i)<j && MESL(i)>=j
                        if rand<0.5
                            P4(i,j,1)=P4IC;
                            P4(i,j,2)=H(i,j,1);
                        end
                    else
                        P4(i,j,1)=0;
                        P4(i,j,2)=0;
                    end
                end
    %             WesternCells=zeros(1,MaxSwampWidth);
    %             EasternCells=zeros(1,MaxSwampWidth);
    %             for mm=1:min(j-1,size(WesternCells,2))
    %                 WesternCells(mm)=delta*H(i,j-mm,1);
    %             end
    %             for mm=1:min(size(EasternCells,2),size(H,2)-j)
    %                 EasternCells(mm)=delta*H(i,j+mm,1);
    %             end
    %             CheckWesternCells=WesternCells<=-.5;%must be within MaxSwampWidth cells to the west of a cell which is outside of marsh range
    %             CheckEasternCells=EasternCells>0;
    %             if sum(CheckWesternCells)==0 || sum(CheckEasternCells)==0
    %                 P4(i,j,1)=-999;
    %                 P4(i,j,2)=0;
    %             end%
    %             if sum(CheckWesternCells)~=0 &&.5<rand
    %                 P4(i,j,1)=P4IC;
    %                 P4(i,j,2)=H(i,j,1);%
    %             end
            elseif (delta*H(i,j,1)-W(i,j))<-0.5 || (delta*H(i,j,1)-W(i,j))>1
                            if (delta*H(i,j,1)<PlantRangeArray(4,1)) %if less than min height (-0.5)
                                P4(i,j,1)=-999;
                                P4(i,j,2)=0;
                            elseif (delta*H(i,j,1)>PlantRangeArray(4,2)) %elseif greater than max height, make 0 (no death by elev. just kill)
                                P4(i,j,1)=0;
                                P4(i,j,2)=0;
                            end
            end
            PC4(i,j)=P4(i,j,1)*(P4(i,j,1)>0);
        end
    end
    t=NaN;
    [P1,P2,P3,P4,P3d]=PlantPropagation06132021(Hstar,t,P1,P2,P3,P4,W,S,delta,P3d,MaxSwampWidth,PlantRangeArray,alpha,DBE,gdrange1,gdrange2,PctMax,MasterMax,MWSL,MESL,ESL,P3Killzone);
        for i=1:size(H,1)
            for j=1:size(H,2)
                PC1(i,j)=P1(i,j,1)*(P1(i,j,1)>0);
                PC2(i,j)=P2(i,j,1)*(P2(i,j,1)>0);                                                                  %
                PC3(i,j)=P3(i,j,1)*(P3(i,j,1)>0);
                PC4(i,j)=P4(i,j,1)*(P4(i,j,1)>0);
            end
        end
    PC=(P1ErosionCoefficient.*PC1(:,:))+(P2ErosionCoefficient.*PC2(:,:))+(P3ErosionCoefficient.*PC3(:,:))+(P4ErosionCoefficient.*PC4(:,:))+(P3dErosionCoefficient.*P3d(:,:,1));
    Ptot=P1(:,:,1)+P2(:,:,1)+P3(:,:,1)+P4(:,:,1)+P3d(:,:,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%END - INITIALIZING PLANT ARRAYS%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%killing plants on half of island test%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if KillSwitch==1;
        for i=floor(0.5*(size(H,1))):size(H,1)     %killing bottom half plants
            for j=1:size(H,2)
                if P1(i,j,1)>0
                    P1(i,j,1)=0;
                    PC1(i,j)=0;
                end
                if P2(i,j,1)>0
                    P2(i,j,1)=0;
                    PC2(i,j)=0;
                end
                if P3(i,j,1)>0
                    P3(i,j,1)=0;
                    PC3(i,j)=0;
                end
                if P4(i,j,1)>0
                    P4(i,j,1)=0;
                    PC4(i,j)=0;
                end
            end
        end
        PC=(P1ErosionCoefficient.*PC1(:,:))+(P2ErosionCoefficient.*PC2(:,:))+(P3ErosionCoefficient.*PC3(:,:))+(P4ErosionCoefficient.*PC4(:,:))+(P3dErosionCoefficient.*P3d(:,:,1));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%END of half island plant test conditions%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end































%%%%%%%%%% FOR INITIAL CONDITION FROM PREVIOUS RUN %%%%%%%%%%%
if restartrun==1  
    
    
    %%%%%% NEW PARAMETERS YOU WILL WANT TO CHANGE %%%%%%
    
    % Enter folder name here!!!!
    filenamestring='St00te702BP2500SLRi0p004359SLRa3p0769en06MMax1'; 
    
    % Enter number of year you want to load data for
    yearnum=20; 
    % Enter folder where main code is stored here.  Data should be stored in
    % subdirectory of this folder
    cd  '/Users/beththomas/Documents/Documents - Beth MacBook Air/MATLAB/Barrier Island' ;

    % Ending time.  IMPORTANT:  For entering time, remember you are not starting at zero!
    % e.g., to restart a simulation after yearnum=27 years and run for another 10 years, use:                         
    time=(yearnum+10)*26;  %2600~100yrs;%1300~50yrs;%780~30yrs;%520~20yrs;260~10yrs;

    %%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%
    
    % Matlab moves to the subdirectory, loads the data, and identifies whether
    % working with Smith or Parramore, then returns to main directory
    curdir=pwd;
    cd(filenamestring);
    filenamestringmat=strcat([filenamestring,'t',int2str(yearnum),'.mat']);
    load(filenamestringmat);
    islandstring=filenamestring(1);
    cd ..

    t0=yearnum*26+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if islandstring=='P'
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % PARRAMORE ISLAND DATA:%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % IslandArea=14400000;
        % 12km x 1200 m = 14400000 - Relative Influence Antecendent... paper, 12.4 m/yr
        %     filename='Parramore03312021.mat'
        IslandArea=19070000; %for parramore
        %     data=importdata(filename);
        %     H=data;                         %H will be the main elevation matrix
        TranChk=1;
    elseif islandstring=='S'
        %%%%%%%%%%%%%%%%%%%%%
        % SMITH ISLAND DATA:%
        %%%%%%%%%%%%%%%%%%%%%
        %     filename='Smith03312021.mat'
        IslandArea=9065000; %for Smith
        %     data=importdata(filename);
        %     H=data;                         %H will be the main elevation matrix
        TranChk=2;
    end


    % filename='EricsIsland.txt'
    % IslandArea=41000;%EricsIsland


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % % % ScaleFactor=floor(sqrt(IslandArea/(sum(sum(H(:,:)>0)))));
    % InterpFactor=.75*sqrt(ScaleFactor);   %ScaleFactor 2
    % % % InterpFactor=0.5*sqrt(ScaleFactor);     %ScaleFactor 4

    % % % H=interp2(H,InterpFactor);
    % % % ScaleFactor=floor(sqrt(IslandArea/(sum(sum(H(:,:)>0)))));
    if TranChk==1
        TranRow=floor(.25*size(H,1));   %choose row for transect for parramore
    elseif TranChk==2
        TranRow=floor(0.75*size(H,1));    %choose row for transect for smith
    else exists(TranChk)==0
        TranRow=floor(0.5*size(H,1));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%END - LOAD ELEVATION MATRIX%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                    %%%%%    PARAMETERS     %%%%%                                                                                                                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                                                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n')
    fprintf('Initializing...')

    delta=0.1;
    % L=1;%<----------ScaleFactor
    L=ScaleFactor;
    BchMax =1.5;   %greatest elevation that the beach spreads inland to
    OW=0; %switch for overwashing - might use when we get storms incorporated (if windspeed>16m/s OW=1 for on, etc)
    BchW=10;    %unused - maximum beach width
    MPswitch=1;     %changes marine processes:    --(0 or any ~=1,2,3) for OFF     
    %                                             --(1) for ON {new version - update using same equil. slope}       MarineProcesses03312021
    %                                             --(2)for ON {old version - update equil. slope every 12 wks}      OLDMarineProcesses03312021
    DepositSupply=3;     %number of slabs of sediment supplied to beach each time MarinePRocesses is run (3 months
    SLRswitch=0;    %turns sea level rise ON(1)/OFF(0 or any ~=1)
    SLRyrs=[ 8 7 6 5]; %years at which to perform SLR at rate SLR and/or migration...
    Ma= 1; %1.0075;       %shoreline migration growth factor, use 1 for no accel/noSLR
    % MigAccel=0;     %will need to track magnitude of migration change 
    MigYr=BruunParam*(SLRvec(2)-SLRvec(1));
    %     MigYr=15; %initial yearly migration of the island in meters/year

    % % % 
    ScaleFactor=floor(sqrt(IslandArea/(sum(sum(H(:,:,1)>0)))));
    % % % 

    %     % QUADRATIC SEA LEVEL RISE (ADDED 4/29, RO)
    %     SLRaccelratevec=[0.0 0.00014 0.00033 0.00051]';         % accel rate (ft/yr^2), tuned to match red, yellow, green curves of VIMS plot
    %     SLRrun=4;                                               %1 - historic; 2 - low; 3 - high; 4 - highest.  Based on SLR graph from VIMS
    %     SLRaccelrate=SLRaccelratevec(SLRrun);
    %     SLRinitrate=0.25/39;                                    %linear SLR (m/yr), tuned to match historic curve of VIMS plot, %0.175/39.3701;   
    %     yrmax=108;
    %     SLRvec=SLRinitrate*(0:yrmax)+SLRaccelrate*(0:yrmax).^2; % Assume quadratic growth from 1992-2100
    %     SLRvec=[SLRinitrate*(-9:-1) SLRvec];                    % Assume linear growth from 1984-1992
    %     BruunParam=500;
        
    % R=zeros(size(H,1),1);
    MaxSwampWidth=1000; %how far the swamp should go out into the water - just have to make something up for now

    % % % WindData=1; %0~steady windspeed/direction ,1~windspeed and direction from empirical data

    if WindData==0
        Windspeed=16;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        WindDir=1;  %wind direction            %  1=N  2=NE  3=E  4=SE  5=S  6=SW  7=W  8=NW  %                                                                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    Windmin=6;
    StormThreshold=16;
    ATcntr=0; %testing variable - to count number of times AT is called

    %%%%%%%%%%%%%%%%%%%%%%%%%   PLANT PARAMETERS    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Plant initial conditions:
%     PlantRangeArray=[1 5;0.75 3;1.5 2.5;-0.5 0]; %all of the elevation ranges  for p1-p4
    % % % P1IC=0.5; 0.9;
    % % % P2IC=0.5; 0.9;
    % % % P3IC=0.5; 0.9;
    % % % P4IC=0.5; 0.9;
    % % % P1PctMax=.6;                %largest percentage we will allow any plant population on a given cell to attain
    % % % P2PctMax=.6;
    % % % P3PctMax=.8;
    % % % P4PctMax=.8;
    % % % PctMax=[P1PctMax P2PctMax P3PctMax P4PctMax];
    % % % MasterMax=0.5;% The most any cell can permit - 80% plant coverage
    KillSwitch=0;   %this will kill plants on the bottom half of the island if set to 1, set to 0 (or anything else) to turn off
    alpha=.01; %propagation rate for each populated cell
%     DBE=.3;     %death by elevation rate for each populated cell outside of plant's elevation range
%     gdrange1=[-.02:.01:.08]; %range of percent values for growth/death for plant populations at (0, 50)% cover
%     gdrange2=[-.02:.01:.08];%[-.04:.01:.04]; %range of percent values for growth/death for plant pops greater than 50% cover


    %%%% CREATE FILENAME FOR THIS RUN USING PARAMETER VALUES ABOVE, ADDED 4/29, RO
    % filenamestring=strcat([islandstring,'t',num2str(time),'SLRi',num2str(SLRinitrate),'SLRa',num2str(SLRaccelrate),'P1',num2str(P1IC),'P2',num2str(P2IC),'P3',num2str(P3IC),'P4',num2str(P4IC),'MMax',num2str(MasterMax),'t0_',num2str(t0)]);
    % filenamestring=strrep(filenamestring,'.','p');
    % filenamestring=strrep(filenamestring,'-','n');

    % CREATE DIRECTORY FOR THIS RUN, ALL RESULTS STORED HERE, ADDED 4/29, RO
    % mkdir(filenamestring);

    %Elevation Matrix dependent parameters:
    ROW=size(H,1);
    COLUMN=size(H,2);
    % % % H(:,:,2)=zeros(ROW,COLUMN);     %Creating extra dimension for H
    % % % Hstar=H;        %initializing dummy matrix
    MigCnt=zeros(ROW,1);   %used to store how many meters the shoreline should have receded by, shoreline moves when MigCnt>ScalingFactor
    M=max(max(H(:,:,1)));
    m=min(min(H(:,:,1)));
    W=zeros(ROW,COLUMN);  %do not comment out - needed as input
    S=zeros(ROW,COLUMN);  %do not comment out - needed as input

    %parameters unused so far (check with Greg)
    % numslabg=0;
    % q=0;
    % Salt=1;
    % ps=0.6;
    % % % AV_ARRAY=zeros(time,1000,3);        %tracks cell movement in AV. routine
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Imaging Options %%%%%

    %  ***note that we do not currently have a running land categorization
    %     subroutine, so leave those options commented or at 0 - I think this
    %     will eventually come back into use so I am leaving the lines there***

    DurImFreq=26*5;      %how often to output "during" elevation images (26~yearly, 130~every 5 years)
        MElevCont=0;
        SElevCont=0;
    PlantDurImFreq=130; %same, just for plant images
    DurCUImFreq=130;
    DursavedataFreq=130;

    ElevImCU_during=1;    %(moved to top)
    Contour30=1;          %contours at 0, 13, and 27 years (moved to top)
    DurContour=0;         %images of the island outline (moved to top)

    TransectIm=1;       %images of 2D island transect
        TranImFreq=260; %every 9 years  
        if TransectIm==1 && exist('TranRow') == 0
            TranRow=floor(.5*size(H,1));   %choose row for transect
        end

    MnPlantCvrIm=1;      %mean percent cover for each plant, displayed at end of routine 1 to turn on, any other number to turn off
    ElevIm_before=1;        ElevIm_during=1;        ElevIm_after=0;
    LandCatIm_before=0;     LandCatIm_during=0;     LandCatIm_after=0;          %leave these at 0 until we have a working LandCategorization routine
    PlantPCIm_before=0;     PlantPCIm_during=1;     PlantPCIm_after=0;
    SingleStepICIm=0;   %plots from initialization process :
    %includes elevation and Plaint initial conditions for
    %one pass through Swamp, PlantPropagation and
    %AeolionTransport
    %     climsH=[min(min(H(:,:,1))) max(max(H(:,:,1)))];
    %     climsP=[-1 1];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%              %%%%%%Initializing Plant Arrays%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %GRASS #1 (Ammophila)(burial resistant)
    % P1=zeros(ROW,COLUMN,2);
    P1Burial=2;         %meters of burial until death
    P1ErosionCoefficient=(2/3);     %used in formula   d=c1P1+c2P2+c3P3+c4P4

    %GRASS #2  (Spartina)
    % P2=zeros(ROW,COLUMN,2);
    P2Burial=0.5;         %meters of burial until death
    P2ErosionCoefficient=(1/3);

    %SHRUB #1  (Morella)(bird-dispersed seeds)
    % P3=zeros(ROW,COLUMN,2);
    P3Burial=1.5;         %meters of burial until death
    P3ErosionCoefficient=1;

    %dead morella will track when morella reaches death and store dead
    %debris data in P3d(:,:,1) for some number of years which will be counted in
    %P3d(:,:,2) - this is only updated inside of PlantProp

    %%%%%%%%%% REED, NEED TO COMMENT THIS LINE OUT ON NEXT TEST!!!
    % % % P3d=zeros(ROW,COLUMN,2);
    %%%%%%%%%%

    P3dErosionCoefficient=1;

    %SECOND TYPE OF SPARTINA!!!
    % P4=zeros(ROW,COLUMN,2);
    P4Burial=0.5;
    P4ErosionCoefficient=(1/3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INITIAL PLANT CONDITIONS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % Each plant species has a preferred "altitude" that it grows best in
    % % % (min and max among all species is 0.75m and 5m, respectively).
    % % % We set each cell of "appropriate" elevation to be covered by some% of the
    % % % respective plant species. For cells which are submerged under water, we
    % % % set the proportion of plant coverage to -1.
    % % %
    % % % We have deleted the initialization of plant population since we are using
    % % % real data.
    % % %It is still important to ensure that the 2nd layer of each
    % % % plant matrix is initiated.
    % % % fprintf('\n')
    % % % fprintf('Seeding the island...')
    % % % PC1=zeros(size(H,1),size(H,2));
    % % % PC2=zeros(size(H,1),size(H,2));
    % % % PC3=zeros(size(H,1),size(H,2));
    % % % PC4=zeros(size(H,1),size(H,2));
    % % % [Hstar,MeanBeachWidth,ESL,WSL,MESL,MWSL,OL]=Shoreline03312021(Hstar,delta,L,BchMax);%function to return shoreline(s)
    % % % for i=1:size(H,1)
    % % %     for j=1:size(H,2)
    % % %         if (delta*H(i,j,1))>=1 && (delta*H(i,j,1))<=5
    % % %             P1(i,j,1)=P1IC;
    % % %             P1(i,j,2)=H(i,j,1);
    % % %         elseif (H(i,j,1)-W(i,j))<0
    % % %             P1(i,j,1)=-999;
    % % %             P1(i,j,2)=0;
    % % %         end
    % % %         PC1(i,j)=P1(i,j,1)*(P1(i,j,1)>0);
    % % %         if (delta*H(i,j,1))>=0.25 && (delta*H(i,j,1))<=3
    % % %             P2(i,j,1)=P2IC;
    % % %             P2(i,j,2)=H(i,j,1);
    % % %         elseif (H(i,j,1)-W(i,j))<0
    % % %             P2(i,j,1)=-999;
    % % %             P2(i,j,2)=0;%
    % % %         end
    % % %         PC2(i,j)=P2(i,j,1)*(P2(i,j,1)>0);
    % % %         if (delta*H(i,j,1))>=1.5 && (delta*H(i,j,1))<=2.5
    % % %             P3(i,j,1)=P3IC;
    % % %             P3(i,j,2)=H(i,j,1);
    % % %         elseif (H(i,j,1)-W(i,j))<0
    % % %             P3(i,j,1)=-999;
    % % %             P3(i,j,2)=0;
    % % %         end
    % % %         PC3(i,j)=P3(i,j,1)*(P3(i,j,1)>0);
    % % %         if (delta*H(i,j,1))>=-0.5 && (delta*H(i,j,1))<=1
    % % %             if ESL(i)~=0
    % % %                 if MWSL(i)<j && MESL(i)>=j
    % % %                     if rand<0.5
    % % %                         P4(i,j,1)=P4IC;
    % % %                         P4(i,j,2)=H(i,j,1);
    % % %                     end
    % % %                 else
    % % %                     P4(i,j,1)=0;
    % % %                     P4(i,j,2)=0;
    % % %                 end
    % % %             end
    % % % %             WesternCells=zeros(1,MaxSwampWidth);
    % % % %             EasternCells=zeros(1,MaxSwampWidth);
    % % % %             for mm=1:min(j-1,size(WesternCells,2))
    % % % %                 WesternCells(mm)=delta*H(i,j-mm,1);
    % % % %             end
    % % % %             for mm=1:min(size(EasternCells,2),size(H,2)-j)
    % % % %                 EasternCells(mm)=delta*H(i,j+mm,1);
    % % % %             end
    % % % %             CheckWesternCells=WesternCells<=-.5;%must be within MaxSwampWidth cells to the west of a cell which is outside of marsh range
    % % % %             CheckEasternCells=EasternCells>0;
    % % % %             if sum(CheckWesternCells)==0 || sum(CheckEasternCells)==0
    % % % %                 P4(i,j,1)=-999;
    % % % %                 P4(i,j,2)=0;
    % % % %             end%
    % % % %             if sum(CheckWesternCells)~=0 &&.5<rand
    % % % %                 P4(i,j,1)=P4IC;
    % % % %                 P4(i,j,2)=H(i,j,1);%
    % % % %             end
    % % %         elseif (delta*H(i,j,1)-W(i,j))<-0.5 || (delta*H(i,j,1)-W(i,j))>1
    % % %                         if (delta*H(i,j,1)<PlantRangeArray(4,1)) %if less than min height (-0.5)
    % % %                             P4(i,j,1)=-999;
    % % %                             P4(i,j,2)=0;
    % % %                         elseif (delta*H(i,j,1)>PlantRangeArray(4,2)) %elseif greater than max height, make 0 (no death by elev. just kill)
    % % %                             P4(i,j,1)=0;
    % % %                             P4(i,j,2)=0;
    % % %                         end
    % % %         end
    % % %         PC4(i,j)=P4(i,j,1)*(P4(i,j,1)>0);
    % % %     end
    % % % end
    t=NaN;
    % % % [P1,P2,P3,P4,P3d]=PlantPropagation05172021(Hstar,t,P1,P2,P3,P4,W,S,delta,P3d,MaxSwampWidth,PlantRangeArray,alpha,DBE,gdrange1,gdrange2,PctMax,MasterMax,MWSL,MESL,ESL,P3Killzone)
    % % %     for i=1:size(H,1)
    % % %         for j=1:size(H,2)
    % % %             PC1(i,j)=P1(i,j,1)*(P1(i,j,1)>0);
    % % %             PC2(i,j)=P2(i,j,1)*(P2(i,j,1)>0);                                                                  %
    % % %             PC3(i,j)=P3(i,j,1)*(P3(i,j,1)>0);
    % % %             PC4(i,j)=P4(i,j,1)*(P4(i,j,1)>0);
    % % %         end
    % % %     end
    % % % PC=(P1ErosionCoefficient.*PC1(:,:))+(P2ErosionCoefficient.*PC2(:,:))+(P3ErosionCoefficient.*PC3(:,:))+(P4ErosionCoefficient.*PC4(:,:))+(P3dErosionCoefficient.*P3d(:,:,1));
    Ptot=P1(:,:,1)+P2(:,:,1)+P3(:,:,1)+P4(:,:,1)+P3d(:,:,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%END - INITIALIZING PLANT ARRAYS%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%killing plants on half of island test%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % if KillSwitch==1;
    % % %     for i=floor(0.5*(size(H,1))):size(H,1)     %killing bottom half plants
    % % %         for j=1:size(H,2)
    % % %             if P1(i,j,1)>0
    % % %                 P1(i,j,1)=0;
    % % %                 PC1(i,j)=0;
    % % %             end
    % % %             if P2(i,j,1)>0
    % % %                 P2(i,j,1)=0;
    % % %                 PC2(i,j)=0;
    % % %             end
    % % %             if P3(i,j,1)>0
    % % %                 P3(i,j,1)=0;
    % % %                 PC3(i,j)=0;
    % % %             end
    % % %             if P4(i,j,1)>0
    % % %                 P4(i,j,1)=0;
    % % %                 PC4(i,j)=0;
    % % %             end
    % % %         end
    % % %     end
    % % %     PC=(P1ErosionCoefficient.*PC1(:,:))+(P2ErosionCoefficient.*PC2(:,:))+(P3ErosionCoefficient.*PC3(:,:))+(P4ErosionCoefficient.*PC4(:,:))+(P3dErosionCoefficient.*P3d(:,:,1));
    % % % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%END of half island plant test conditions%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
















t=t0;


%%                 %%%%%%%%BEFORE IMAGES%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ElevIm_before==1
    Hfilt=round(imgaussfilt(H(:,:,1),2),0);
    %%%%%%%%Elevation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure                                          %
    colormap(jet())                                 %
    imagesc(H(:,:,1),climsH)                        %
    hold on
    if MElevCont==1
    contour(H(:,:,1),[-0 -0],'color','m','linewidth',1.1)
    end
    if SElevCont==1
    contour(Hfilt(:,:,1),[-0 -0],'color','k','linewidth',1.1)
    end
    if MElevCont==1 && SElevCont==1
        lgnd=legend('\color{white} marsh contour (-0.5m)','\color{white} sea level contour (0m)');
        set(lgnd,'color','none','location','southeast');
    elseif MElevCont==1 && SElevCont==0
        lgnd=legend('\color{white} island contour (0 m)');
        set(lgnd,'color','none','location','southeast');        
    elseif SElevCont==1 && MElevCont==0
        lgnd=legend('\color{white} sea level contour (0m)');
        set(lgnd,'color','none','location','southeast');
    end
    colorbar                                        %
    title('Before Image');                          %
    %%%%
    figureHandle = gcf;
    figfontsize=15; % sometimes matlab chooses pretty small fonts...
    set(gca,'fontsize',figfontsize);
    set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
    figfilenamestring=strcat(['FigElb_',filenamestring]);
    cd(filenamestring);
    if restartrun==0
    if saveasdotfig==1 saveas(gcf,figfilenamestring); end% save as .fig
    saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
    end
    cd ..
    %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if PlantPCIm_before==1
    %%%%%%PLANTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% AMMOPHILA
    figure
    colormap gray
    imagesc(P1(:,:,1),climsP)
    colorbar
    title('P1 - Ammophila IC')
    %%%%
    figureHandle = gcf;
    figfontsize=15; % sometimes matlab chooses pretty small fonts...
    set(gca,'fontsize',figfontsize);
    set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
    figfilenamestring=strcat(['FigP1b_',filenamestring]);
    cd(filenamestring);
    if restartrun==0
    if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
    saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
    end
    cd ..
    %%%%

    %%% SPARTINA
    figure
    colormap gray
    imagesc(P2(:,:,1),climsP)
    colorbar
    title('P2 - Spartina patens IC')
    %%%%
    figureHandle = gcf;
    figfontsize=15; % sometimes matlab chooses pretty small fonts...
    set(gca,'fontsize',figfontsize);
    set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
    figfilenamestring=strcat(['FigP2b_',filenamestring]);
    cd(filenamestring);
    if restartrun==0
    if saveasdotfig==1 saveas(gcf,figfilenamestring); end% save as .fig
    saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
    end
    cd ..
    %%%%
    
    %%% MORELLA
    figure
    colormap gray
    imagesc(P3(:,:,1),climsP)
    colorbar
    title('P3 - Morella IC')
    %%%%
    figureHandle = gcf;
    figfontsize=15; % sometimes matlab chooses pretty small fonts...
    set(gca,'fontsize',figfontsize);
    set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
    figfilenamestring=strcat(['FigP3b_',filenamestring]);
    cd(filenamestring);
    if restartrun==0
    if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
    saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
    end
    cd ..
    %%%%
    
    %%% SPARTINA (marsh)
    figure
    colormap gray
    imagesc(P4(:,:,1),climsP)
    colorbar
    title('P4 - Spartina alterniflora IC')
    %%%%
    figureHandle = gcf;
    figfontsize=15; % sometimes matlab chooses pretty small fonts...
    set(gca,'fontsize',figfontsize);
    set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
    figfilenamestring=strcat(['FigP4b_',filenamestring]);
    cd(filenamestring);
    if restartrun==0
    if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
    saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
    end
    cd ..
    %%%%
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%Single Step Subroutine Images%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if restartrun==0
fprintf('\n')


fprintf('Running startup sequences:')
fprintf('\n')
fprintf('Running Swamp Processes...')
%Initial - Run Swamp
    [Hstar,ColumnArraySwamp1,ColumnArraySwamp2,AdjascentLengthSwamp,OppositeLocationSwamp,flag,INum]=SwampProcesses03312021(Hstar,delta,L,flag,PC);
fprintf('done')
fprintf('\n')
fprintf('Running Avalanche...')
% (1.5) post swamp avalache
PassCount=1;
CellsMoved=zeros(1,1000);
[Hstar,flag,CellCt]=AVALANCHE03312021(Hstar,delta,L,flag,PC);
CellsMoved(PassCount)=CellCt;
while flag==1
    PassCount=PassCount+1;
    [Hstar,flag,CellCt]=AVALANCHE03312021(Hstar,delta,L,flag,PC);
    CellsMoved(PassCount)=CellCt;%
end
CellsMoved=CellsMoved(1,1:PassCount);
SwampAvArray=[NaN PassCount CellsMoved];    %tracking number of passes and number of cells moved per pass
AvIC_Array1=SwampAvArray;                   %AvIC_Array1 tracks avalanche data for initialization routines
H=Hstar;                                        %update H
fprintf('done')


%****************** NOT USING LAND CATEGORIZATION - BUT MAY AGAIN LATER - DO NOT DELETE!!
for i=1:size(Hstar,1)
    for Icount=1:INum
        if ColumnArraySwamp2(i,Icount)~=0
%             for j=ColumnArraySwamp1(i,Icount):ColumnArraySwamp2(i,Icount)
%                 Hstar(i,j,2)=2;
%             end
        end
    end
end
%****************** NOT USING LAND CATEGORIZATION - BUT MAY AGAIN LATER - DO NOT DELETE!!

%post swamp processes/avalanche1 images
if SingleStepICIm==1
    figure('NumberTitle','off','Name','After Swamp');
    colormap(jet())        %
    %     clims=[min(min(Hstar(:,:,1))) max(max(Hstar(:,:,1)))];          %
    imagesc(Hstar(:,:,1),climsH)%
    title('ICs After Swamp')
    colorbar
    %%%%
    figureHandle = gcf;
    figfontsize=15; % sometimes matlab chooses pretty small fonts...
    set(gca,'fontsize',figfontsize);
    set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
    figfilenamestring=strcat(['FigHstarb_',filenamestring]);
    cd(filenamestring);
    if restartrun==0
    if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
    saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
    end
    cd ..
    %%%%

end
fprintf('\n')
fprintf('Running Marine Processes...')

%Initial - Marine Processes
MeanBeachWidthVec=zeros(1,time+2);  %for tracking beach width evolution
if MPswitch==1
    % need version of PC with erosion coefficients AND with negative values (zeros only where plants COULD grow) to calculate migration
    PCmp=-999*ones(ROW,COLUMN);
    for i=1:ROW
        for j=1:COLUMN
            if P1(i,j,1)>=0 || P2(i,j,1)>=0 ||P3(i,j,1)>=0 %||P4(i,j,1)>=0
                PCmp(i,j)=(P1ErosionCoefficient*P1(i,j,1)*(P1(i,j,1)>0))+(P2ErosionCoefficient*P2(i,j,1)*(P2(i,j,1)>0))+(P3ErosionCoefficient*P3(i,j,1)*(P3(i,j,1)>0));%+(P4ErosionCoefficient*P4(i,j,1)*(P4(i,j,1)>0));
            end
        end
    end         
%     Pt1=max(P1(:,:,1),0);
%     Pt2=max(P2(:,:,1),0);
%     Pt3=max(P3(:,:,1),0);
%     Pt4=max(P4(:,:,1),0);
%     PCmp=(P1ErosionCoefficient.*Pt1)+(P2ErosionCoefficient.*Pt2)+(P3ErosionCoefficient.*Pt3)+(P4ErosionCoefficient.*Pt4);
    [Hstar,P3,MeanBeachWidth,ESL,WSL,MESL,MWSL,OppositeLocation,SLRyrs,MigCnt,MigAccel]=MarineProcesses05282021(Hstar,delta,L,P1,P2,P3,P4,BchMax,OW,t,SLRyrs,PCmp,IslandArea,MigCnt,ScaleFactor,MigYr,Ma,MigAccel,MasterMax,SLRvec,BruunParam);
    MeanBeachWidthVec(1)=MeanBeachWidth;
    H=Hstar;
elseif MPswitch==2
    [Hstar,PlantColumnArray,PlantColumnArray2,P3,flag,MeanBeachWidth]=OLDMarineProcesses03312021(Hstar,delta,L,flag,PC,P3,BchMax,BchW);
    MeanBeachWidthVec(1)=MeanBeachWidth;
    H=Hstar;
end
fprintf('done')

if SingleStepICIm==1
    figure('NumberTitle','off','Name','After Marine Processes');
    colormap(jet())        %
    imagesc(Hstar(:,:,1),climsH)%
    title('ICs after Swamp, Marine Processes')
    colorbar
    %%%%
    figureHandle = gcf;
    figfontsize=15; % sometimes matlab chooses pretty small fonts...
    set(gca,'fontsize',figfontsize);
    set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
    figfilenamestring=strcat(['FigHstarb_',filenamestring]);
    cd(filenamestring);
    if restartrun==0
    if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
    saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
    end
    cd ..
    %%%%

end

%Redeclare PC since P3 has been updated
%PC=(P1ErosionCoefficient.*P1(:,:,1))+(P2ErosionCoefficient.*P2(:,:,1))+(P3ErosionCoefficient.*P3(:,:,1))+(P4ErosionCoefficient.*P4(:,:,1));                                             %
for i=1:size(H,1)
    for j=1:size(H,2)
        PC1(i,j)=P1(i,j,1)*(P1(i,j,1)>0);
        PC2(i,j)=P2(i,j,1)*(P2(i,j,1)>0);                                                                  %
        PC3(i,j)=P3(i,j,1)*(P3(i,j,1)>0);
        PC4(i,j)=P4(i,j,1)*(P4(i,j,1)>0);
    end
end
PC=(P1ErosionCoefficient.*PC1(:,:))+(P2ErosionCoefficient.*PC2(:,:))+(P3ErosionCoefficient.*PC3(:,:))+(P4ErosionCoefficient.*PC4(:,:))+(P3dErosionCoefficient.*P3d(:,:,1));

%post marine processes avalanche
PassCount=1;
CellsMoved=zeros(1,1000);
[Hstar,flag,CellCt]=AVALANCHE03312021(Hstar,delta,L,flag,PC);
CellsMoved(PassCount)=CellCt;
while flag==1
    PassCount=PassCount+1;
    [Hstar,flag,CellCt]=AVALANCHE03312021(Hstar,delta,L,flag,PC);
    CellsMoved(PassCount)=CellCt;%
end
CellsMoved=CellsMoved(1,1:PassCount);
AvIC_Array2=[NaN PassCount CellsMoved];
H=Hstar;

%****************** NOT USING LAND CATEGORIZATION - BUT MAY AGAIN LATER - DO NOT DELETE!!
%     %Initial Declare all categories in H
PlantColumnArray=zeros([size(H,1) 1]); %P3 isn't allowed to grow on the beach - this is to check for that
PlantColumnArray2=zeros([size(H,1) 1]);
[Hstar,P1,P2,P3,P4]=LandCategorizationUNUSED03312021(Hstar,P1,P2,P3,P4,PlantColumnArray,PlantColumnArray2,ColumnArraySwamp1,ColumnArraySwamp2,INum);
%****************** NOT USING LAND CATEGORIZATION - BUT MAY AGAIN LATER - DO NOT DELETE!!

if SingleStepICIm==1
    figure
    colormap(jet())        
    imagesc(Hstar(:,:,2),climsH)
    colorbar
    ('ICs LandCats - After Swamp, Marine Processes and Avalanche');
    %%%%
    figureHandle = gcf;
    figfontsize=15; % sometimes matlab chooses pretty small fonts...
    set(gca,'fontsize',figfontsize);
    set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
    figfilenamestring=strcat(['FigHstarb2_',filenamestring]);
    cd(filenamestring);
    if restartrun==0
    if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
    saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
    end
    cd ..
    %%%%    
end

fprintf('\n')
fprintf('Initialization complete')
end

    H2=Hstar(:,:,1);
    [n1, n2]=size(H2(:,:));
    ESL_hist=zeros(n1,time+1);
    top_hist = zeros(time+1,n2);
    mid_hist = zeros(time+1,n2);
    bot_hist = zeros(time+1,n2);
    H_cell_count = zeros(3,time+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%END INITIAL PARAMETERIZATION%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
tic
T=zeros(time,1);
%#########################################################################%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!!@@##$$               %%%%%%%%MAIN LOOP%%%%%%%%                 $$##@@!!%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%##########################################################################
tStart=cputime;
for t=t0:time
    TimeStep=t
    Ptot=P1(:,:,1)+P2(:,:,1)+P3(:,:,1)+P4(:,:,1)+P3d(:,:,1);
    if 0==mod(t,26)
        if t==0
            fprintf('\n')
            fprintf('Here we go...')
        elseif t==26
            fprintf('\n')
            fprintf('%d year has passed',t/26)
        else
            fprintf('\n')
            fprintf('%d years have passed',t/26)
        end
    end
    %     TimeStep=t
    %%    %PC=(P1ErosionCoefficient.*P1(:,:,1))+(P2ErosionCoefficient.*P2(:,:,1))+(P3ErosionCoefficient.*P3(:,:,1))+(P4ErosionCoefficient.*P4(:,:,1));
    %need to work this into the inside of PlantProcesses
    for i=1:size(H,1)
        for j=1:size(H,2)
            PC1(i,j)=P1(i,j,1)*(P1(i,j,1)>0);
            PC2(i,j)=P2(i,j,1)*(P2(i,j,1)>0);                                                                  %
            PC3(i,j)=P3(i,j,1)*(P3(i,j,1)>0);
            PC4(i,j)=P4(i,j,1)*(P4(i,j,1)>0);
        end
    end
    PC=(P1ErosionCoefficient.*PC1(:,:))+(P2ErosionCoefficient.*PC2(:,:))+(P3ErosionCoefficient.*PC3(:,:))+(P4ErosionCoefficient.*PC4(:,:))+(P3dErosionCoefficient.*P3d(:,:,1));
    
% make copy of H at start
if t == 0
    H_initial = H;
end

    %%    %(1) Run Swamp

    [Hstar,ColumnArraySwamp1,ColumnArraySwamp2,AdjascentLengthSwamp,OppositeLocationSwamp,flag,INum]=SwampProcesses03312021(Hstar,delta,L,flag,PC);

    %%    %(1.AV)            THIS AVALANCHE RUNS ONCE PER YEAR (every 26 time steps) AND IS THE FIRST AVALANCHE       
    % if (0==mod(t,26)) && t~=0
    PassCount=1;
    CellsMoved=zeros(1,1000);
    [Hstar,flag,CellCt]=AVALANCHEtime03312021(Hstar,delta,L,flag,PC,t);
    CellsMoved(PassCount)=CellCt;
    while flag==1
        PassCount=PassCount+1;
        [Hstar,flag,CellCt]=AVALANCHEtime03312021(Hstar,delta,L,flag,PC,t);
        CellsMoved(PassCount)=CellCt;%
        if CellCt<=5
            flag=0;
        end
    end
    CellsMoved=CellsMoved(1,1:PassCount);
    SwampAvArray=[NaN PassCount CellsMoved];
    AvIC_Array1=SwampAvArray;
    AV_ARRAY(TimeStep+1,1:PassCount+2,1)=SwampAvArray;

    
    %****************** NOT USING LAND CATEGORIZATION - BUT MAY AGAIN LATER - DO NOT DELETE!!
    for i=1:size(Hstar,1)
        for Icount=1:INum
            if ColumnArraySwamp2(i,Icount)~=0
                for j=ColumnArraySwamp1(i,Icount):ColumnArraySwamp2(i,Icount)
                    Hstar(i,j,2)=2;
                end
            end
        end
    end
    %****************** NOT USING LAND CATEGORIZATION - BUT MAY AGAIN LATER - DO NOT DELETE!!
    
    
    %%    %(2) Aeolian Transport
    if WindData==0
        if Windspeed>=6
            WindCnt=floor(((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*24);
            for day=1:14
                %                 for hour =1:WindCnt
                    [Hstar]=AeolianTransport03212021(Hstar,delta,L,W,Windspeed,Windmin,StormThreshold,WindDir,PC);
                %                 end
            end%
        end
    end
    
    if WindData==1 && t>0
        %WindDir=randi([1,8]);  can randomize wind or use data-weighted wind directions outlined below. (comment this if you randomize)
        for day=1:14
            WindSpeedParam=randi(2183);
            if WindSpeedParam>=1 && WindSpeedParam<=2116
                Windspeed=randi([0,5]);
            elseif WindSpeedParam>2116 && WindSpeedParam<=2183
                ATcntr=ATcntr+1;
                Windspeed=randi([6,15]);
            elseif WindSpeedParam>2183
                Windspeed=16;
            end
            WindParam=randi(2107);
            if WindParam>=1 && WindParam<=250
                WindDir=1;
            elseif WindParam>=251 && WindParam<=655
                WindDir=2;
            elseif WindParam>=656 && WindParam<=835
                WindDir=3;
            elseif WindParam>=836 && WindParam<=1091
                WindDir=4;
            elseif WindParam>=1092 && WindParam<=1464
                WindDir=5;
            elseif WindParam>=1465 && WindParam<=1717
                WindDir=6;
            elseif WindParam>=1718 && WindParam<=1865
                WindDir=7;
            elseif WindParam>=1866 && WindParam<=2107
                WindDir=8;
            end

            if Windspeed>=6
%                 WindCnt=floor(((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*24);
%                 for hour=1:WindCnt
                        [Hstar]=AeolianTransport03312021(Hstar,delta,L,W,Windspeed,Windmin,StormThreshold,WindDir,PC);
%                 end
            end
        end
    end
    
    %%    %(2.AV)                       THIS AVALANCHE RUNS ONCE PER YEAR (every 26 time steps) AND IS THE SECOND AVALANCHE
    % if (0==mod(t,26)) && t~=0
    PassCount=1;
    CellsMoved=zeros(1,1000);
    [Hstar,flag,CellCt]=AVALANCHEtime03312021(Hstar,delta,L,flag,PC,t);
    CellsMoved(PassCount)=CellCt;
    while flag==1
        PassCount=PassCount+1;
        [Hstar,flag,CellCt]=AVALANCHEtime03312021(Hstar,delta,L,flag,PC,t);
        CellsMoved(PassCount)=CellCt;%
        %             if CellCt<=5
        %                 flag=0;
        %             end
    end
    CellsMoved=CellsMoved(1,1:PassCount);
    Av2_Array=[TimeStep PassCount CellsMoved];
    AV_ARRAY(TimeStep+1,1:PassCount+2,2)=Av2_Array;
    % end
   %% Effect of current moving sand south (Beth)
%  this code is set up to move 1 cell at time

if (0==mod(t,130)) && t~=0 %runs once ever 5 years
 [Hstar]=TipsCode10052022(Hstar); 
end


    %% Smoothing routine (average height across 3x3 sqaure)
    if (0==mod(t,26))
    [Hstar] = Smoothing2023(Hstar);
    end

%     if (0==mod(t,13))
%     [Hstar] = Smoothing_water(Hstar);
%     end

%% Overwash Code (Beth modification of Kiran's code)

[H_overwash, PC1_overwash, PC2_overwash, PC3_overwash, PC4_overwash, DC, n1, surge_hist, EastSL]=Overwash072324_v2(Hstar, PC, PC1, PC2, PC3, PC4, time, t, surge_hist); 

Hstar(:,:,1) = H_overwash(:,:,1); %new matrix Hstar stores this movement
surge_hist = surge_hist;
if t == 0
     EastSL_initial = EastSL;
 end
 if t ~= 0
     EastSL_final = EastSL;
     ESL_movement = EastSL_initial - EastSL_final;
 end
PC1(:,:,1) = PC1_overwash(:,:,1);
PC2(:,:,1) = PC2_overwash(:,:,1);
PC3(:,:,1) = PC3_overwash(:,:,1);
PC4(:,:,1) = PC4_overwash(:,:,1);


%% Storing Eastern Shore Line and transects of H 
% [ESL_hist] = storing_ESL(Hstar,t,ESL_hist);
% ESL_hist = ESL_hist;
% 
[top_hist, mid_hist, bot_hist, H_cell_count] = storing_H_transects(Hstar,t,top_hist, mid_hist, bot_hist, H_cell_count);
top_hist = top_hist;
mid_hist = mid_hist;
bot_hist = bot_hist;
H_cell_count = H_cell_count;

%commented out on 7/9/24 to test new Parramore island wihtout adjusting
%this routine

    %%    %(3) Plant Propogation
%land categroization used to be required for plant processes - this might be removeable
    %****************** NOT USING LAND CATEGORIZATION - BUT MAY AGAIN LATER - DO NOT DELETE!!
%     [Hstar,P1,P2,P3,P4]=LandCategorizationUNUSED03312021(Hstar,P1,P2,P3,P4,PlantColumnArray,PlantColumnArray2,ColumnArraySwamp1,ColumnArraySwamp2,INum);
    %****************** NOT USING LAND CATEGORIZATION - BUT MAY AGAIN LATER - DO NOT DELETE!!
    
    if (0==mod(t,13))  %  choose 26 for 2-week timesteps;    52 for 1-week timesteps
        [Hstar,MeanBeachWidth,ESL,WSL,MESL,MWSL,OL]=Shoreline03312021(Hstar,delta,L,BchMax);%get shoreline to use in PlantPropagation
        [P1,P2,P3,P4,P3d]=PlantPropagation06132021(Hstar,t,P1,P2,P3,P4,W,S,delta,P3d,MaxSwampWidth,PlantRangeArray,alpha,DBE,gdrange1,gdrange2,PctMax,MasterMax,MWSL,MESL,ESL,P3Killzone);
        if P3IC==0 %killing off any morella that was established due to bird activity
            P3(:,:,1)=min(P3(:,:,1),0);
        end
    end
    
    for i=1:size(H,1)
        for j=1:size(H,2)
            PC1(i,j)=P1(i,j,1)*(P1(i,j,1)>0);
            PC2(i,j)=P2(i,j,1)*(P2(i,j,1)>0);                                                                  %
            PC3(i,j)=P3(i,j,1)*(P3(i,j,1)>0);
            PC4(i,j)=P4(i,j,1)*(P4(i,j,1)>0);
        end
    end
    PC=(P1ErosionCoefficient.*PC1(:,:))+(P2ErosionCoefficient.*PC2(:,:))+(P3ErosionCoefficient.*PC3(:,:))+(P4ErosionCoefficient.*PC4(:,:))+(P3dErosionCoefficient.*P3d(:,:,1));
    Ptot=PC1(:,:,1)+PC2(:,:,1)+PC3(:,:,1)+PC4(:,:,1)+P3d(:,:,1);
    
    %%   %(4) Marine Processes - edits
    %currently every 12 months, moves 2 slabs into ocean
    if (0==mod(t,26))
        if MPswitch==1
            PCmp=-999*ones(ROW,COLUMN);
            for i=1:ROW
                for j=1:COLUMN
                    if P1(i,j,1)>=0 || P2(i,j,1)>=0 ||P3(i,j,1)>=0 %||P4(i,j,1)>=0
                        PCmp(i,j)=(P1ErosionCoefficient*P1(i,j,1)*(P1(i,j,1)>0))+(P2ErosionCoefficient*P2(i,j,1)*(P2(i,j,1)>0))+(P3ErosionCoefficient*P3(i,j,1)*(P3(i,j,1)>0));%+(P4ErosionCoefficient*P4(i,j,1)*(P4(i,j,1)>0));
                    end
                end
            end
                %             PCmp=(P1ErosionCoefficient.*P1(:,:,1))+(P2ErosionCoefficient.*P2(:,:,1))+(P3ErosionCoefficient.*P3(:,:,1))+(P4ErosionCoefficient.*P4(:,:,1))+(P3dErosionCoefficient.*P3d(:,:,1));
            [Hstar,P3,MeanBeachWidth,ESL,WSL,MESL,MWSL,OL,SLRyrs,MigCnt,MigAccel]=MarineProcesses05282021(Hstar,delta,L,P1,P2,P3,P4,BchMax,OW,t,SLRyrs,PCmp,IslandArea,MigCnt,ScaleFactor,MigYr,Ma,MigAccel,MasterMax,SLRvec,BruunParam);
        elseif MPswitch==2
            [Hstar,PlantColumnArray,PlantColumnArray2,P3,flag,MeanBeachWidth]=OLDMarineProcesses03312021(Hstar,delta,L,flag,PC,P3,BchMax,BchW);
        end
        if MPswitch~=0
            MeanBeachWidthVec(1)=MeanBeachWidth;
            H=Hstar;
        end
        
    end

    %%  %(4.AV*)update PC before avalanching again
    %need to work this into the inside of PlantProcesses
    for i=1:size(H,1)
        for j=1:size(H,2)
            %         PC1(i,j)=P1(i,j,1)*(P1(i,j,1)>0);
            %         PC2(i,j)=P2(i,j,1)*(P2(i,j,1)>0);                                                                  %
            PC3(i,j)=P3(i,j,1)*(P3(i,j,1)>0);       %p3 is only thing that changes inside of marine processes - killed on beach
            %         PC4(i,j)=P4(i,j,1)*(P4(i,j,1)>0);
        end
    end
    PC=(P1ErosionCoefficient.*PC1(:,:))+(P2ErosionCoefficient.*PC2(:,:))+(P3ErosionCoefficient.*PC3(:,:))+(P4ErosionCoefficient.*PC4(:,:))+(P3dErosionCoefficient.*P3d(:,:,1));
    %%   %(4.AV) Avalanche after Marine Processes   THIS AVALANCHE RUNS EVERY TIME STEP AND IS THE FINAL AVALANCHE
    PassCount=1;
    CellsMoved=zeros(1,1000);
    [Hstar,flag,CellCt]=AVALANCHEtime03312021(Hstar,delta,L,flag,PC,t);
    CellsMoved(PassCount)=CellCt;
    while flag==1
        PassCount=PassCount+1;
        [Hstar,flag,CellCt]=AVALANCHEtime03312021(Hstar,delta,L,flag,PC,t);
        CellsMoved(PassCount)=CellCt;
    end
    CellsMoved=CellsMoved(1,1:PassCount);
    Av3_Array=[TimeStep PassCount CellsMoved];
    AV_ARRAY(TimeStep+1,1:PassCount+2,3)=Av3_Array;

  
    %%      %(5) Data and categories in H
    %****************** NOT USING LAND CATEGORIZATION - BUT MAY AGAIN LATER - DO NOT DELETE!!
%     [Hstar,P1,P2,P3,P4]=LandCategorizationUNUSED03312021(Hstar,P1,P2,P3,P4,PlantColumnArray,PlantColumnArray2,ColumnArraySwamp1,ColumnArraySwamp2,INum);
    %****************** NOT USING LAND CATEGORIZATION - BUT MAY AGAIN LATER - DO NOT DELETE!!
    [Hstar,Centroid]=Qdata03312021(Hstar,P1,P2,P3,P4);%,ESL,WSL,MESL,MWSL,OppositeLocation);
    CentroidVec(t+1,:)=Centroid;
    
    %calculating mean percent cover for each Pi
    p1=P1(:,:,1)>0;p2=P2(:,:,1)>0;p3=P3(:,:,1)>0;p4=P4(:,:,1)>0;
    P1gt=P1(:,:,1).*p1;
    P2gt=P2(:,:,1).*p2;
    P3gt=P3(:,:,1).*p3;
    P4gt=P4(:,:,1).*p4;
    
    MnP1(t+2)=sum(P1gt,'all')/sum(p1,'all');
    MnP2(t+2)=sum(P2gt,'all')/sum(p2,'all');
    MnP3(t+2)=sum(P3gt,'all')/sum(p3,'all');
    MnP4(t+2)=sum(P4gt,'all')/sum(p4,'all');
    
    % calculating the number of cells <0 , =0, and >0 at each time step
    error_vec(t+1,1) = t;
    error_vec(t+1,2) = sum(sum(H(:,:,1)<0));
    error_vec(t+1,3) = sum(sum(H(:,:,1)==0));
    error_vec(t+1,4) = sum(sum(H(:,:,1)>0));
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%                         Sea Level Rise                          %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SLRswitch==1
        if 0==mod(t,156)%0==mod(t,156)    %156 timesteps = 312 weeks = 6 years
            %    if 0==mod(t,260)
            Hstar=Hstar-1; %working estimate for SLR in VA is 1.66cm/yr = .1m/6yrs
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%                       End Sea Level Rise                        %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    H=round(Hstar,0);        %*****THIS IS THE ONLY UPDATE FOR H IN THE MAIN LOOP - OTHERWISE ALL UPDATES ARE MADE TO HSTAR
    
    if t==1
        Sv1=H(:,:,1);
        Sv1f=round(imgaussfilt(H(:,:,1),4),0);
        Sv1Cent=Centroid;
    elseif t==26*15
        Sv2=H(:,:,1);
        Sv2f=round(imgaussfilt(H(:,:,1),4),0);
        Sv2Cent=Centroid;
    elseif t==26*30
        Sv3=H(:,:,1);
        Sv3f=round(imgaussfilt(H(:,:,1),4),0);
        Sv3Cent=Centroid;
    end
    
    tEnd=cputime-tStart;
    
    
    T(t+1)=toc;
    fprintf('\n')
    if t~=0
        Tnew=T(t+1)-T(t);
    else
        Tnew=T(1);
    end
    fprintf('t=%d time: %0.1f',TimeStep,Tnew)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                %%%%%%%%    During Images    %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % % % % %ELEVATION
    
    if ElevIm_during==1 
        Hfilt=round(imgaussfilt(H(:,:,1),4),0);
        if (0==mod(t,DurImFreq))                                                                                                                                 %
            figure                                                                                                                                    %
            colormap(jet())
            imagesc(H(:,:,1),climsH);
            %hold on
%             sz = 15;
%             scatter(DC, 1:n1,sz,'w', 'filled');
            hold on
            if MElevCont==1 %make MElevCont == 1 in line 100 or so %THIS OUTLINES SHORELINE
                contour(H(:,:,1),[-0 -0],'color','b','linewidth',1.1) %try changing to [0 0]
                hold on
                contour(H_initial(:,:,1),[-0 -0],'color','m','linewidth',1.1)
            end
            if SElevCont==1
                contour(Hfilt(:,:,1),[-0 -0],'color','k','linewidth',1.1)
            end
            if MElevCont==1 && SElevCont==1
                lgnd=legend('\color{white} marsh contour (-0.5m)','\color{white} sea level contour (0m)');
                set(lgnd,'color','none','location','southeast');
            elseif MElevCont==1 && SElevCont==0
                lgnd=legend('\color{white} dune crest');
                set(lgnd,'color','none','location','southeast');
            elseif SElevCont==1 && MElevCont==0
                lgnd=legend('\color{white} sea level contour (0m)');
                set(lgnd,'color','none','location','southeast');
            end
            colorbar
            title(sprintf('Elevation at %d years',t/26))
            hold off
            %%%%
            figureHandle = gcf;
            figfontsize=15; % sometimes matlab chooses pretty small fonts...
            set(gca,'fontsize',figfontsize);
            set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
            figfilenamestring=strcat(['FigElt',num2str(t/26),'_',filenamestring]);
            cd(filenamestring);
            if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
            saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
            cd ..
            %%%%
    
        end
    end

    %This is the unsmoothed transect migration image
% % %     if TransectIm==1
% % %         if 0==mod(t,TranImFreq)
% % %             Trow=H(TranRow,:,1);
% % %             Trow(Trow<-5)=-6;
% % %             xlim1=find(Trow>=-5,1,'first')-100; xlim2=find(Trow>=-5,1,'last')+100;
% % %             figure(1000)
% % %             hold on
% % %             plot(1:size(H,2),Trow,'Linewidth',2)
% % %             xlim([xlim1 xlim2])
% % %             ylim([-5,max(Trow)+20])
% % %             title(sprintf('Island transect at row %d',TranRow))
% % %             hold off
% % %         end
% % %     end
    %%%Close up elev:
    if ElevImCU_during==1
        if (0==mod(t,DurCUImFreq))
            Wndw=250;
            MidRow=round(Centroid(2));
            MidCol=round(Centroid(1));
            figure
            imagesc(H(MidRow-Wndw:MidRow+Wndw,MidCol-Wndw:MidCol+Wndw,1),climsH);
            hold on
            ContourVec=[0:10:max(max(Hfilt))];
            colorbar
            contour(Hfilt(MidRow-Wndw:MidRow+Wndw,MidCol-Wndw:MidCol+Wndw,1),ContourVec);
            contour(Hfilt(MidRow-Wndw:MidRow+Wndw,MidCol-Wndw:MidCol+Wndw,1),ContourVec,'LineColor','k');
            [C,h]=contour(Hfilt(MidRow-Wndw:MidRow+Wndw,MidCol-Wndw:MidCol+Wndw,1),ContourVec,'LineColor','k');
            clabel(C,h)
            title(sprintf('Elevation Close-Up about Centroid at %d years',t/26))
            %%%%
            figureHandle = gcf;
            figfontsize=15; % sometimes matlab chooses pretty small fonts...
            set(gca,'fontsize',figfontsize);
            set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
            figfilenamestring=strcat(['FigECUt',num2str(t/26),'_',filenamestring]);
            cd(filenamestring);
            if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
            saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
            cd ..
            %%%%
% % % old one (lines at 10, 20 close up):
%             figure                                                                                                                                    %
%             colormap(jet())
%             MidPt=round(CentroidVec(t+1,:));
%             MidCol=MidPt(1); MidRow=MidPt(2);
%             imagesc(H(MidRow-50:MidRow+50,MidCol-50:MidCol+50,1),climsH);
%             hold on
%             ContourVec=[0:5:max(max(Hfilt))];
%             contour(Hfilt(MidRow-50:MidRow+50,MidCol-50:MidCol+50,1),[10 10],'Color','m','linewidth',1.1);
%             contour(Hfilt(MidRow-50:MidRow+50,MidCol-50:MidCol+50,1),[20 20],'Color','k','linewidth',1.1);
%             lgnd=legend('\color{white} marsh contour (-0.5m)','\color{white} sea level contour (0m)');
%             set(lgnd,'color','none','location','southeast');
%             colorbar
%             title(sprintf('Elevation Close-Up about Centroid at %d years',t/26))
%             hold off


            hFig=figure;
            set(gcf,'PaperPositionMode','auto')
            set(hFig,'units','inches')
            hFig.Position(1)=.5;
            hFig.Position(2)=.5;
            hFig.Position(3)=2*hFig.Position(3);
            hFig.Position(4)=1.25*hFig.Position(4);
            subplot(1,2,1)
            colormap(jet())
            imagesc(H(MidRow-Wndw:MidRow+Wndw,MidCol-Wndw:MidCol+Wndw,1),climsH);
            title(sprintf('Elevation Close-Up about Centroid at %d years',t/26))
            subplot(1,2,2)
            [C,h] = contour(Hfilt(MidRow-Wndw:MidRow+Wndw,MidCol-Wndw:MidCol+Wndw,1),ContourVec);
            clabel(C,h)
            set(gca,'YDir','reverse')
            title(sprintf('Contours of Elevation Close-Up about Centroid at %d years',t/26))
            %%%%
            figureHandle = gcf;
            figfontsize=15; % sometimes matlab chooses pretty small fonts...
            set(gca,'fontsize',figfontsize);
            set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
            figfilenamestring=strcat(['FigECCt',num2str(t/26),'_',filenamestring]);
            cd(filenamestring);
            if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
            saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
            cd ..
            %%%%

            
%             QUADBOX images for Plants - controlled by ElevImCU_during
            figure
            subplot(2,2,1)
            colormap bone
            imagesc(P1(MidRow-Wndw:MidRow+Wndw,MidCol-Wndw:MidCol+Wndw,1),climsP);
            title(sprintf('P_1 about centroid at %d years',t/26))
            subplot(2,2,2)
            colormap bone
            imagesc(P2(MidRow-Wndw:MidRow+Wndw,MidCol-Wndw:MidCol+Wndw,1),climsP);
            title(sprintf('P_2 about centroid at %d years',t/26))
            subplot(2,2,3)
            colormap bone
            imagesc(P3(MidRow-Wndw:MidRow+Wndw,MidCol-Wndw:MidCol+Wndw,1),climsP);
            title(sprintf('P_3 about centroid at %d years',t/26))
            subplot(2,2,4)
            colormap bone
            imagesc(P4(MidRow-Wndw:MidRow+Wndw,MidCol-Wndw:MidCol+Wndw,1),climsP);
            title(sprintf('P_4 about centroid at %d years',t/26))
            %%%%
            figureHandle = gcf;
            figfontsize=15; % sometimes matlab chooses pretty small fonts...
            set(gca,'fontsize',figfontsize);
            set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
            figfilenamestring=strcat(['FigPCt',num2str(t/26),'_',filenamestring]);
            cd(filenamestring);
            if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
            saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
            cd ..
            %%%%

        end
    end
    
    % % % % %Land Category
    
    if LandCatIm_during==1
        if (0==mod(t,DurImFreq)) %every 5 years                                                                                                                                  %
            figure                                                                                                                                  %
            colormap(jet())
            clims=[min(min(H(:,:,2))) max(max(H(:,:,2)))];
            imagesc(H(:,:,2),clims);
            colorbar
            title(sprintf('Land Category at time step = %d',t))
            %%%%
            figureHandle = gcf;
            figfontsize=15; % sometimes matlab chooses pretty small fonts...
            set(gca,'fontsize',figfontsize);
            set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
            figfilenamestring=strcat(['FigLCt',num2str(t/26),'_',filenamestring]);
            cd(filenamestring);
            if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
            saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
            cd ..
            %%%%

        end
    end
    
    % % % % %Plant Percent Cover%
    
    if PlantPCIm_during==1
        if (0==mod(t,PlantDurImFreq))  %  2~/mo, choose 26~/year
            
            %%% AMMOPHILA
            figure
            colormap gray
            %        clims=[-.1 max(0,max(max(P1(:,:,1))))];
            imagesc(P1(:,:,1),climsP)
            colorbar
            title(sprintf('P1-Ammophila at %d years',t/26))
            %        title(sprintf('P1-Ammophila at time step = %d',t))
            %%%%
            figureHandle = gcf;
            figfontsize=15; % sometimes matlab chooses pretty small fonts...
            set(gca,'fontsize',figfontsize);
            set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
            figfilenamestring=strcat(['FigP1t',num2str(t/26),'_',filenamestring]);
            cd(filenamestring);
            if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
            saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
            cd ..
            %%%%

            %%% SPARTINA
            figure
            colormap gray
            %        clims=[-.1 max(0,max(max(P2(:,:,1))))];
            imagesc(P2(:,:,1),climsP)
            colorbar%
            title(sprintf('P2 - Spartina patens at %d years',t/26))
            %        title(sprintf('P2 - Spartina patens at time step = %d',t))
            %%%%
            figureHandle = gcf;
            figfontsize=15; % sometimes matlab chooses pretty small fonts...
            set(gca,'fontsize',figfontsize);
            set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
            figfilenamestring=strcat(['FigP2t',num2str(t/26),'_',filenamestring]);
            cd(filenamestring);
            if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
            saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
            cd ..
            %%%%

            %%% MORELLA
            figure
            colormap gray
            %        clims=[-.1 max(0,max(max(P3(:,:,1))))];
            imagesc(P3(:,:,1),climsP)
            colorbar%
            title(sprintf('P3 - Morella at %d years',t/26))
            %        title(sprintf('P3 - Morella at time step = %d',t))
            %%%%
            figureHandle = gcf;
            figfontsize=15; % sometimes matlab chooses pretty small fonts...
            set(gca,'fontsize',figfontsize);
            set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
            figfilenamestring=strcat(['FigP3t',num2str(t/26),'_',filenamestring]);
            cd(filenamestring);
            if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
            saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
            cd ..
            %%%%

            %%% DEAD MORELLA
            figure
            colormap gray
            %        clims=[-.1 max(0,max(max(P3d(:,:,1))))];
            imagesc(P3d(:,:,1),climsP)
            colorbar
            title(sprintf('P3d - DEAD Morella at %d years',t/26))
            %        title(sprintf('P3d - DEAD Morella at time step = %d',t))
            %%%%
            figureHandle = gcf;
            figfontsize=15; % sometimes matlab chooses pretty small fonts...
            set(gca,'fontsize',figfontsize);
            set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
            figfilenamestring=strcat(['FigP3dt',num2str(t/26),'_',filenamestring]);
            cd(filenamestring);
            if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
            saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
            cd ..
            %%%%

            %%% SPARTINA (marsh)
            figure
            colormap gray
            %        clims=[-.1 max(1,max(max(P4(:,:,1))))];
            imagesc(P4(:,:,1),climsP)
            colorbar
            title(sprintf('P4 - Spartina alterniflora at %d years',t/26))
            %        title(sprintf('P4 - Spartina alterniflora at time step = %d',t))
            %%%%
            figureHandle = gcf;
            figfontsize=15; % sometimes matlab chooses pretty small fonts...
            set(gca,'fontsize',figfontsize);
            set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
            figfilenamestring=strcat(['FigP4t',num2str(t/26),'_',filenamestring]);
            cd(filenamestring);
            if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
            saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
            cd ..
            %%%%

        end
        
    end
    if sum(ismember(Dursavedatavec,t))>0
        cd(filenamestring);
        if t==26 disp(strcat(['Saving Data after 1 year']));
        else disp(strcat(['Saving Data after ',int2str(t/26),' years']));
        end
        save(strcat([filenamestring,'t',int2str(t/26),'.mat']),'H','Hstar','P1','P2','P3','P4','P3d','PC','ScaleFactor','MigAccel','WindData','PctMax','MasterMax','AV_ARRAY','Ma','t','SLRvec','BruunParam','P1IC','P2IC','P3IC','P4IC','SLRinitrate','SLRaccelrate','BruunParam','MnP1','MnP2','MnP3','MnP4','Sv1','Sv2','Sv3','Sv1f','Sv2f','Sv3f','Sv1Cent','Sv2Cent','Sv3Cent','climsH','climsP','P3Killzone','DBE','gdrange1','gdrange2','PlantRangeArray');
        cd ..
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  END OF MAIN LOOP  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%         %%%%%%%%    ESL plots over time    %%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if time == 1300
 figure
plot(ESL_hist(:,1),'b.');
hold on 
plot(ESL_hist(:,651),'r.');
hold on 
plot(ESL_hist(:,1301),'g.');
legend({'Initial', '25 years', '50 years'},'FontSize',14);
xlabel('Row location (North to South)', 'FontSize',14);
ylabel('Shore line column location', 'FontSize',14);
title('Eastern Shore Line Movement Over Time','FontSize',14)
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%         %%%%%%%%    Elevation plots over time    %%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if time == 1300
% figure;
% subplot(3,1,1)
% plot(top_hist(1,:),'b.')
% hold on
% plot(top_hist(651,:),'r.')
% hold on
% plot(top_hist(1301,:),'g.')
% legend({'Initial', '25 years', '50 years'},'FontSize',14);
% xlabel('Column location (West to East)', 'FontSize',14);
% ylabel('Elevation', 'FontSize',14);
% title('Elevation change at top of island','FontSize',14)
% 
% subplot(3,1,2)
% plot(mid_hist(1,:),'b.')
% hold on
% plot(mid_hist(651,:),'r.')
% hold on
% plot(mid_hist(1301,:),'g.')
% legend({'Initial', '25 years', '50 years'},'FontSize',14);
% xlabel('Column location (West to East)', 'FontSize',14);
% ylabel('Elevation', 'FontSize',14);
% title('Elevation change at middle of island','FontSize',14)
% 
% subplot(3,1,3)
% plot(bot_hist(1,:),'b.')
% hold on
% plot(bot_hist(651,:),'r.')
% hold on
% plot(bot_hist(1301,:),'g.')
% legend({'Initial', '25 years', '50 years'},'FontSize',14);
% xlabel('Column location (West to East)', 'FontSize',14);
% ylabel('Elevation', 'FontSize',14);
% title('Elevation change at bottom of island','FontSize',14)
% end

if time == 1300
figure;
subplot(2,1,1)
plot(top_hist(1,:),'b.')
hold on
plot(top_hist(5,:),'r.')
hold on
plot(top_hist(14,:),'g.')
legend({'Initial', '25 years', '50 years'},'FontSize',14);
xlabel('Column location (West to East)', 'FontSize',14);
ylabel('Elevation', 'FontSize',14);
title('Elevation change at top of island','FontSize',14)
subplot(2,1,2)
plot(bot_hist(1,:),'b.')
hold on
plot(bot_hist(5,:),'r.')
hold on
plot(bot_hist(14,:),'g.')
legend({'Initial', '25 years', '50 years'},'FontSize',14);
xlabel('Column location (West to East)', 'FontSize',14);
ylabel('Elevation', 'FontSize',14);
title('Elevation change at bottom of island','FontSize',14)
end

if time == 1300
figure;
plot(log(H_cell_count(1,:)),'b.')
hold on
plot(log(H_cell_count(2,:)),'r.')
hold on
plot(log(H_cell_count(3,:)),'g.')
legend({'Cells above 20', 'Cells between 10 and 20', 'Cells between 0 and 10'},'FontSize',14);
xlabel('Time in two weeks', 'FontSize',14);
xlim([0 1301]);
ylim([10 15.5]);
ylabel('log(number of cells)', 'FontSize',14);
title('Number of Cells in Different Elevation Classes','FontSize',14)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    After Images    %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % %Elevation
if ElevIm_after==1
    figure
    colormap(jet())
    imagesc(H(:,:,1),climsH);
    hold on
    contour(H(:,:,1),[-0.5 -0.5],'color','black')
    colorbar
    title('After Image')
    %%%%
    figureHandle = gcf;
    figfontsize=15; % sometimes matlab chooses pretty small fonts...
    set(gca,'fontsize',figfontsize);
    set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
    figfilenamestring=strcat(['FigElat',num2str(t/26),'_',filenamestring]);
    cd(filenamestring);
    if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
    saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
    cd ..
    %%%%

end

%Mean Plant Percent Cover Images
if MnPlantCvrIm==1
    figure
    hold on
    plot(MnP1)
    plot(MnP2)
    plot(MnP3)
    plot(MnP4)
    legend('P_1','P_2','P_3','P_4')
    xlabel('time in two week steps')
    ylabel('Mean Percent Cover')
    title(sprintf('Mean Percent Cover - %.0f years',time/26));
    %%%%
    figureHandle = gcf;
    figfontsize=15; % sometimes matlab chooses pretty small fonts...
    set(gca,'fontsize',figfontsize);
    set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
    figfilenamestring=strcat(['FigPCat',num2str(t/26),'_',filenamestring]);
    cd(filenamestring);
    if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
    saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
    cd ..
    %%%%

end

% % % % %Land Category

if LandCatIm_after==1
    figure
    colormap(jet())
    clims=[min(min(H(:,:,2))) max(max(H(:,:,2)))];
    imagesc(H(:,:,2),clims);
    colorbar
    title(sprintf('Land Category - After'))
    %%%%
    figureHandle = gcf;
    figfontsize=15; % sometimes matlab chooses pretty small fonts...
    set(gca,'fontsize',figfontsize);
    set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
    figfilenamestring=strcat(['FigLCat',num2str(t/26),'_',filenamestring]);
    cd(filenamestring);
    if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
    saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
    cd ..
    %%%%

end

% % % % %PLANTS (percent cover)

if PlantPCIm_after==1
    %%% AMMOPHILA
    figure
    colormap gray
    imagesc(P1(:,:,1),climsP)
    colorbar
    title(sprintf('P1-Ammophila After %d years',time/26));
    %%%%
    figureHandle = gcf;
    figfontsize=15; % sometimes matlab chooses pretty small fonts...
    set(gca,'fontsize',figfontsize);
    set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
    figfilenamestring=strcat(['FigP1at',num2str(t/26),'_',filenamestring]);
    cd(filenamestring);
    if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
    saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
    cd ..
    %%%%
    
    %%% SPARTINA
    figure
    colormap gray
    imagesc(P2(:,:,1),climsP)
    colorbar
    title(sprintf('P2 - Spartina patens after %d years',time/26));
    %%%%
    figureHandle = gcf;
    figfontsize=15; % sometimes matlab chooses pretty small fonts...
    set(gca,'fontsize',figfontsize);
    set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
    figfilenamestring=strcat(['FigP2at',num2str(t/26),'_',filenamestring]);
    cd(filenamestring);
    if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
    saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
    cd ..
    %%%%
    
    %%% MORELLA
    figure
    colormap gray
    imagesc(P3(:,:,1),climsP)
    colorbar%
    title(sprintf('P3 - Morella after %d years',time/26));
    %%%%
    figureHandle = gcf;
    figfontsize=15; % sometimes matlab chooses pretty small fonts...
    set(gca,'fontsize',figfontsize);
    set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
    figfilenamestring=strcat(['FigP3at',num2str(t/26),'_',filenamestring]);
    cd(filenamestring);
    if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
    saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
    cd ..
    %%%%
    
    %%% SPARTINA (marsh)
    figure
    colormap gray
    imagesc(P4(:,:,1),climsP)
    colorbar
    title(sprintf('P4 - Spartina alterniflora after %d years',time/26));
    %%%%
    figureHandle = gcf;
    figfontsize=15; % sometimes matlab chooses pretty small fonts...
    set(gca,'fontsize',figfontsize);
    set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
    figfilenamestring=strcat(['FigP4at',num2str(t/26),'_',filenamestring]);
    cd(filenamestring);
    if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
    saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
    cd ..
    %%%%
    
end

%% Contours at 0, 15, and 30 years
% if Contour30==1
% %Contours of the shoreline (0m elevation)    
% figure                                                                                                                                    %
% colormap(jet())
% hold on
% contour(Sv1f(:,:,1),[1 1],'Color','blue','LineWidth',2);           %beginning contour (after initialization)
% contour(Sv2f(:,:,1),[1 1],'Color','green','LineWidth',2);          %contour at 13 years (t=338)
% contour(Sv3f(:,:,1),[1 1],'Color','red','LineWidth',2);            %contour at 27 years (t=702)
% scatter(Sv1Cent(1),Sv1Cent(2),'o','MarkerEdgeColor','k','MarkerFaceColor','b')
% scatter(Sv2Cent(1),Sv2Cent(2),'o','MarkerEdgeColor','k','MarkerFaceColor','g')
% scatter(Sv3Cent(1),Sv3Cent(2),'o','MarkerEdgeColor','k','MarkerFaceColor','r')
% scatter(Sv1Cent(1),Sv1Cent(2),'o','MarkerEdgeColor','w',...
%                                    'MarkerFaceColor','b',...
%                                    'LineWidth',1.25)
% scatter(Sv2Cent(1),Sv2Cent(2),'o','MarkerEdgeColor','w',...
%                                    'MarkerFaceColor','g',...
%                                    'LineWidth',1.25)                               
% scatter(Sv3Cent(1),Sv3Cent(2),'o','MarkerEdgeColor','w',...
%                                    'MarkerFaceColor','r',...
%                                    'LineWidth',1.25)
% set(gca, 'YDir','reverse')
% legend('Initial','15 years','30 years','Location','southeast')
% title('Island Contours at 0m Elevation');
% %%%%
% figureHandle = gcf;
% figfontsize=15; % sometimes matlab chooses pretty small fonts...
% set(gca,'fontsize',figfontsize);
% set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
% figfilenamestring=strcat(['FigCont',num2str(t/26),'_',filenamestring]);
% cd(filenamestring);
% if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
% saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
% cd ..
% % %%%%
% % 
% % 
% % %contours including the marsh area:
% % % % % figure                                                                                                                                    %
% % % % % colormap(jet())
% % % % % hold on
% % % % % contour(Sv1f(:,:,1),[-5 -5],'Color','blue','LineWidth',2);           %beginning contour (after initialization)
% % % % % contour(Sv2f(:,:,1),[-5 -5],'Color','green','LineWidth',2);          %contour at 13 years (t=338)
% % % % % contour(Sv3f(:,:,1),[-5 -5],'Color','red','LineWidth',2);            %contour at 27 years (t=702)
% % % % % scatter(Sv1Cent(1),Sv1Cent(2),'o','MarkerEdgeColor','k','MarkerFaceColor','b')
% % % % % scatter(Sv2Cent(1),Sv2Cent(2),'o','MarkerEdgeColor','k','MarkerFaceColor','g')
% % % % % scatter(Sv3Cent(1),Sv3Cent(2),'o','MarkerEdgeColor','k','MarkerFaceColor','r')
% % % % % % scatter(Sv3Cent(1),Sv3Cent(2),'*r')
% % % % % % scatter(Sv2Cent(1),Sv2Cent(2),'*g')
% % % % % % scatter(Sv1Cent(1),Sv1Cent(2),'*b')
% % % % % set(gca, 'YDir','reverse')
% % % % % legend('Initial','13 years','27 years','Location','southeast')
% % % % % title('Island Contours with Marshland (at -0.5m Elevation)')
% end
% % 
% if TransectIm==1
%     figure
%     Trow1=Sv1f(TranRow,:,1);
%     Trow2=Sv2f(TranRow,:,1);
%     Trow3=Sv3f(TranRow,:,1);
%     Trow1(Trow1<-5)=-6; Trow2(Trow2<-5)=-6; Trow3(Trow3<-5)=-6;
%     xlim1=find(Trow>=-5,1,'first')-100; xlim2=find(Trow>=-5,1,'last')+100;
%     hold on
%     plot(1:size(H,2),Trow1,'-','Linewidth',2)
%     plot(1:size(H,2),Trow2,':','Linewidth',2)
%     plot(1:size(H,2),Trow3,':','Linewidth',2)
%     xlim([xlim1 xlim2])
%     ylim([-5,max(Trow)+20])
%     title(sprintf('Island transect at row %d',TranRow))
%     legend('Initial','15 years','30 years','Location','northwest')
%     hold off;
%     %%%%
%     figureHandle = gcf;
%     figfontsize=15; % sometimes matlab chooses pretty small fonts...
%     set(gca,'fontsize',figfontsize);
%     set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
%     figfilenamestring=strcat(['FigTransat',num2str(t/26),'_',filenamestring]);
%     cd(filenamestring);
%     if saveasdotfig==1 saveas(gcf,figfilenamestring); end % save as .fig
%     saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
%     cd ..
    %%%%
    
% end

fprintf('DONE!')
toc

%!!!!!!!!!!!!!!!!!!!!!!!!!no, seriously ... it's over.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

