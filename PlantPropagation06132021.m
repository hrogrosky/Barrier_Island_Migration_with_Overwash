function[P1,P2,P3,P4,P3d]=PlantPropagation06132021(Hstar,t,P1,P2,P3,P4,W,S,delta,P3d,MaxSwampWidth,PlantRangeArray,alpha,DBE,gdrange1,gdrange2,PctMax,MasterMax,MWSL,MESL,ESL,P3Killzone)
% function[P1,P2,P3,P4,P3d]=PlantPropagation03312021(Hstar,P1,P2,P3,P4,W,S,delta,P1Burial,P2Burial,P3Burial,P4Burial,P3d,MaxSwampWidth,t)

test=zeros(size(P1,1),size(P2,2));
clims=[0 0];

P1B4=P1(:,:,1);
P2B4=P2(:,:,1);
P3B4=P3(:,:,1);
P4B4=P4(:,:,1);

% first loop does all propagating and death of all populations,
% second loop is death by competition for cells > MasterMAX


%   To run this file you will need to specify:
%       H     - elevation matrix that is continually used in main code
%       P1    - the plant matrix for Ammophila (GRASS)
%       P2    - the plant matrix for Spartina  (GRASS)
%       P3    - the plant matrix for Morella   (SHRUB)
%       W     - the elevation matrix for the water table
%       S     - the matrix which determines available salinity at each cell
%     delta   - the thickness of each slab
%
%   The routine will return the matrix for each of the plant species after
%   propagating.


% fprintf('I am running PP! ')
%recently moved to main code:
% PlantRangeArray=[1 5;0.75 3;1.5 2.5;-0.5 1]; %all of the elevation ranges  for p1-p4
% alpha=.01; %propagation rate for each populated cell
% DBE=.3;     %death by elevation rate for each populated cell outside of plant's elevation range
% gdrange1=[-.02:.01:.08]; %range of percent values for growth/death for plant populations at (0, 50)% cover
% gdrange2=[-.02:.01:.08];%[-.04:.01:.04]; %range of percent values for growth/death for plant pops greater than 50% cover
% P1PctMax=.6;                %largest percentage we will allow any plant population on a given cell to attain
% P2PctMax=.6;
% P3PctMax=.8;
% P4PctMax=.8;
% PctMax=[P1PctMax P2PctMax P3PctMax P4PctMax];
% MasterMax=0.8;% The most any cell can permit - 80% plant coverage


SwampWidth=MaxSwampWidth; % number of cells wide that the swamp should be - should find a better way of establishing, 10 'looks right' for now
WesternCells=zeros(1,SwampWidth);

dH=delta*Hstar(:,:,1);          %to make sure everything works correctly in this subroutine we work in meters and not slabs - since plant data is in meters




%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%DEATH BY WATER TABLE STUFF%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%not using-we have(?) data, but I've never seen us use water table concept


% HWcheck=H(:,:,1)-W(:,:);
% HWcheck1=(HWcheck>0);
% HWcheck2a(:,:)=(delta*HWcheck>=-0.5);
% HWcheck2b=(delta*HWcheck<=1);

% for i=1:size(H,1)
%     for j=1:size(H,2)
%         if H(i,j,2)~=-1
%         if HWcheck2a(i,j)==1
%             if HWcheck2b(i,j) ==1
%                 if H(i,j,2)==2
%                     if P4(i,j,1)==-999
%                         P4(i,j,1)=0;
%                     end
%                 end
%             end
%         end
%
%         if HWcheck1(i,j)==1
%             if P1(i,j,1)==-999
%                 P1(i,j,1)=0;
%             end
%
%             if P2(i,j,1)==-999
%                P2(i,j,1)=0;
%             end
%
%             if P3(i,j,1)==-999
%                 P3(i,j,1)=0;
%             end
%         end
%         if HWcheck1(i,j)==0
%             P1(i,j,1)=-999;
%             P1(i,j,2)=0;
%             P2(i,j,1)=-999;
%             P2(i,j,2)=0;
%             P3(i,j,1)=-999;
%             P3(i,j,2)=0;
%         end
%         end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%GROWTH/DEATH/PROP%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%splitting P_i into arrays for P_i(:,:,1) and P_i(:,:,2)
PA=zeros(size(Hstar,1),size(Hstar,2),4); %(P)lant (A)rray for tracking plant pct cover in growth/death/prop loop
PA(:,:,1)=P1(:,:,1);
PA(:,:,2)=P2(:,:,1);
PA(:,:,3)=P3(:,:,1);
PA(:,:,4)=P4(:,:,1);
PHA=zeros(size(Hstar,1),size(Hstar,2),4); %(P)lant (H)eight (A)rray for tracking plant init. elev. in growth/death/prop loop
PHA(:,:,1)=P1(:,:,2);
PHA(:,:,2)=P2(:,:,2);
PHA(:,:,3)=P3(:,:,2);
PHA(:,:,4)=P4(:,:,2);

% plants do not propagate during initialization, skip straight to beath by competition
if isnan(t)==0
    %   Beginning of the year, plants use full gdrange (death and growth)
    if mod(t,26)==0
        for i=1:size(PA,3)
            PAB4=PA(:,:,i);
            Px=PA(:,:,i);
            Ph=PHA(:,:,i);
            for R=1:size(Hstar,1)
                for C=1:size(Hstar,2)
                    if dH(R,C)<-0.5
                        Px(R,C)=-999;
                        Ph(R,C)=0;
                    elseif ESL(R)~=0 && C>ESL(R)
                        Px(R,C)=-999;
                        Ph(R,C)=0;
                    else
                        if i==3
                            P3d(R,C,2)=P3d(R,C,2)-1*(P3d(R,C,2)>0);     %if i=3 (working on P3) we remove a counter from the p3d array if there is a value stored there
                            %                 if P3d(R,C,2)==0
                            %                     P3(R,C,1)=0;
                            %                 end
                            if dH(R,C)<0 && P3d(R,C,1)>0                %if elev. goes below zero we get rid of dead morella
                                P3d(R,C,1)=0;
                                P3d(R,C,2)=0;
                            end
                        end
                        if i~=4                                         %work on all Pi arrays except p4 first
                            if (dH(R,C)<PlantRangeArray(i,1)) || (dH(R,C)>PlantRangeArray(i,2)) %if we are outside of the elevation range for this plant...
                                if Px(R,C)~=-999                                                    %but we aren't underwater...
                                    if i==3 && Px(R,C)>0 && dH(R,C)>0 &&t~=0         %if it's p3 outside of elevation range we have to create a dead morella value
                                        PxB4=Px(R,C);                           %store p3 %cover before death by elev.
                                        Px(R,C)=Px(R,C) - DBE*Px(R,C);%-999;          %changed on 9/21 so that death by elevation is not a sudden drop to 0
                                        if PxB4>0.01 && Px(R,C)<0.01             %if newly below 1% we create a p3d cell and kill off remaining p3
                                            Px(R,C)=0;
                                            Ph(R,C)=0;
                                            P3d(R,C,1)=.05;
                                            P3d(R,C,2)=10;
                                        end
                                    elseif Px(R,C) >0 && dH(R,C)>=0                                       %if it's p1 or p2 outside of elevation we remove a percentage
                                        Px(R,C)=Px(R,C)-DBE*Px(R,C);%-999;          %changed on 9/21 see above
                                        if Px(R,C)<0.01                        %if that percentage drops below .1% we kill it
                                            Px(R,C) = 0;                         %since it is outside of it's elev. range we negate the cell for future popln.
                                            Ph(R,C)=0;
                                        end
                                    else                                    %should only have cells now underwater which we will negate
                                        
                                        Ph(R,C)=0;
                                    end
                                end
                            end
                            if (dH(R,C)>=PlantRangeArray(i,1)) && (dH(R,C)<=PlantRangeArray(i,2))   %if we are inside of the plant elevation range
                                %                 resetting negated cells which were previously submerged
                                if Px(R,C)==-999
                                    Px(R,C)=0;
                                    Ph(R,C)=Hstar(R,C,1);
                                end
                                %                 find the neighborhood, calc new popln for
                                %                 growth/propgatn :
                                %                   some % * (how many neighbors are
                                %                 populated)
                                Nhd=NewNhd(R,C,PAB4);
                                if Px(R,C) <.5
                                    PxB4=Px(R,C);
                                    kk=randi([1,length(gdrange1)]);
                                    beta=gdrange1(kk);
                                    Px(R,C)=min(Px(R,C)+sum(beta*Nhd),PctMax(i));%+beta*(beta>0)*(i==3),PctMax(i));          %grow by sum of neighbors* beta, if it's morella grow a little more?
                                    if i==3 && PxB4>=0.05 && Px(R,C)<=.05                                                      % add's chance of spread by birds anywhere inside of p3 range (same below)
                                        Px(R,C)=0;
                                        Ph(R,C)=0;
                                        P3d(R,C,1)=.05;
                                        P3d(R,C,2)=5;
                                    end
                                else
                                    kk=randi([1,length(gdrange2)]);
                                    beta=gdrange2(kk);
                                    Px(R,C)=min(Px(R,C)+sum(beta*Nhd),PctMax(i));%+beta*(beta>0)*(i==3),PctMax(i));          %grow by sum of neighbors* beta, if it's morella grow a little more?
                                end
                                %                record elevations for all cells inside the range,
                                %                elevation=zero if no popln exists on cell currently
                                if Px(R,C)==0
                                    Ph(R,C)=0;
                                elseif Px(R,C)>0
                                    Ph(R,C)=Hstar(R,C,1);
                                elseif Px(R,C)<0
                                    Px(R,C)=0;
                                    Ph(R,C)=0;
                                end
                            end
                        elseif i==4
                            if (dH(R,C)<PlantRangeArray(i,1)) %if less than min height (-0.5)
                                if Px(R,C)~=-999    %negate (should already be done from first loop after R,C declared
                                    
                                    Ph(R,C)=0;
                                end
                            elseif (dH(R,C)>PlantRangeArray(i,2)) %elseif greater than max height, make 0 (no death by elev. just kill)
                                Px(R,C)=0;
                                Ph(R,C)=0;
                            end
                            if (dH(R,C)>=PlantRangeArray(i,1)) && (dH(R,C)<=PlantRangeArray(i,2))
                                if ESL(R)~=0
                                    if MWSL(R)<=C && MESL(R)>=C
                                        if Px(R,C)==-999
                                            Px(R,C)=0;
                                            Ph(R,C)=0;
                                        end   
                                        Nhd=NewNhd(R,C,Px);
                                        if Px(R,C) <.5
                                            kk=randi([1,length(gdrange1)]);
                                            beta=gdrange1(kk);
                                            Px(R,C)=min(Px(R,C)+sum(beta*Nhd)+beta*(beta>0)*(i==3),PctMax(i));   %this also kills by a
                                            
                                        else
                                            kk=randi([1,length(gdrange2)]);
                                            beta=gdrange2(kk);
                                            Px(R,C)=min(Px(R,C)+sum(beta*Nhd)+beta*(beta>0)*(i==3),PctMax(i));
                                        end
                                    else
                                        Px(R,C)=0;
                                    end
                                end
                                %                record elevations for all cells inside the range,
                                %                elevation=zero if no popln exists on cell currently
                                if Px(R,C)==0
                                    Ph(R,C)=0;
                                elseif Px(R,C)>0
                                    Ph(R,C)=Hstar(R,C,1);
                                elseif Px(R,C)<0 && Px(R,C)~=-999
                                    Px(R,C)=0;
                                    Ph(R,C)=0;
                                end
                            end
                        end
                    end
                end
            end
            PA(:,:,i)=Px;
            PHA(:,:,i)=Ph;
        end
        
        
        %   Half year, plants only grow (no death, simulate spring time growth)
    else
        Spr_gdrange1=gdrange1(gdrange1>0);
        Spr_gdrange2=gdrange2(gdrange2>0);
        for i=1:size(PA,3)
            PAB4=PA(:,:,i);
            Px=PA(:,:,i);
            Ph=PHA(:,:,i);
            for R=1:size(Hstar,1)
                for C=1:size(Hstar,2)
                    if dH(R,C)<-.5
                        Px(R,C)=-999;
                        Ph(R,C)=0;
                    elseif ESL(R)>0 && C>ESL(R)
                        Px(R,C)=-999;
                        Ph(R,C)=0;
                    else
                        if i==3
                            if dH(R,C)<0 && P3d(R,C,1)>0                %if elev. goes below zero we get rid of dead morella
                                P3d(R,C,1)=0;
                                P3d(R,C,2)=0;
                            end
                        end
                        if i~=4                                         %work on all Pi arrays except p4 first
                            if (dH(R,C)<PlantRangeArray(i,1)) || (dH(R,C)>PlantRangeArray(i,2)) %if we are outside of the elevation range for this plant...
                                if Px(R,C)~=-999                                                    %but we aren't underwater...
                                    if i==3 && Px(R,C)>0 && dH(R,C)>0 &&t~=0         %if it's p3 outside of elevation range we have to create a dead morella value
                                        PxB4=Px(R,C);                           %store p3 %cover before death by elev.
                                        Px(R,C)=Px(R,C) - DBE*Px(R,C);%-999;          %changed on 9/21 so that death by elevation is not a sudden drop to 0
                                        if PxB4>0.01 && Px(R,C)<0.01            %if newly below 1% we create a p3d cell and kill off remaining p3
                                            Px(R,C)=0;
                                            Ph(R,C)=0;
                                            P3d(R,C,1)=.05;
                                            P3d(R,C,2)=5;
                                        end
                                    elseif Px(R,C) >0 && dH(R,C)>=0                                       %if it's p1 or p2 outside of elevation we remove a percentage
                                        Px(R,C)=Px(R,C)-DBE*Px(R,C);%-999;          %changed on 9/21 see above
                                        if Px(R,C)<0.001                        %if that percentage drops below .1% we kill it
                                            Px(R,C) = 0;                         %since it is outside of it's elev. range we drop to 0, will negate after elev drops below-0.5.
                                            Ph(R,C)=0;
                                        end
                                    else                                    %should only have cells now underwater which we will negate
                                        Px(R,C)=0;
                                        Ph(R,C)=0;
                                    end
                                end
                            end
                            if i~=3 && (dH(R,C)>=PlantRangeArray(i,1)) && (dH(R,C)<=PlantRangeArray(i,2))   %if we are inside of the plant elevation range, don't do morella this time
                                %                 resetting negated cells which were previously submerged
                                if Px(R,C)==-999
                                    Px(R,C)=0;
                                    Ph(R,C)=Hstar(R,C,1);
                                end
                                %                 find the neighborhood, calc new popln for
                                %                 growth/propgatn :
                                %                   some % * (how many neighbors are
                                %                 populated)
                                Nhd=NewNhd(R,C,PAB4);
                                if Px(R,C) <.5
                                    PxB4=Px(R,C);
                                    kk=randi([1,length(Spr_gdrange1)]);
                                    beta=Spr_gdrange1(kk);
                                    Px(R,C)=min(Px(R,C)+sum(beta*Nhd),PctMax(i));%+beta*(beta>0)*(i==3),PctMax(i));          %grow by sum of neighbors* beta, if it's morella grow a little more?
                                    if i==3 && PxB4>=0.05 && Px(R,C)<=.05                                                      % add's chance of spread by birds anywhere inside of p3 range (same below)
                                        Px(R,C)=0;
                                        Ph(R,C)=0;
                                        P3d(R,C,1)=.05;
                                        P3d(R,C,2)=10;
                                    end
                                else
                                    kk=randi([1,length(Spr_gdrange2)]);
                                    beta=Spr_gdrange2(kk);
                                    Px(R,C)=min(Px(R,C)+sum(beta*Nhd),PctMax(i));%+beta*(beta>0)*(i==3),PctMax(i));          %grow by sum of neighbors* beta, if it's morella grow a little more?
                                end
                                %                record elevations for all cells inside the range,
                                %                elevation=zero if no popln exists on cell currently
                                if Px(R,C)==0
                                    Ph(R,C)=0;
                                elseif Px(R,C)>0
                                    Ph(R,C)=Hstar(R,C,1);
                                elseif Px(R,C)<0
                                    Px(R,C)=0;
                                    Ph(R,C)=0;
                                end
                            end
                        elseif i==4
                            if (dH(R,C)<PlantRangeArray(i,1)) %if less than min height (-0.5)
                                if Px(R,C)~=-999    %negate (should already be done from first loop after R,C declared
                                    Px(R,C)=-999;
                                    Ph(R,C)=0;
                                end
                            elseif (dH(R,C)>PlantRangeArray(i,2)) %elseif greater than max height, make 0 (no death by elev. just kill)
                                Px(R,C)=0;
                                Ph(R,C)=0;
                            end
                            if (dH(R,C)>=PlantRangeArray(i,1)) && (dH(R,C)<=PlantRangeArray(i,2))
                                if ESL(R)~=0
                                    if MWSL(R)<=C && MESL(R)>=C
                                        if Px(R,C)==-999
                                            Px(R,C)=0;
                                            Ph(R,C)=0;
                                        end   
                                        Nhd=NewNhd(R,C,Px);
                                        if Px(R,C) <.5
                                            kk=randi([1,length(gdrange1)]);
                                            beta=gdrange1(kk);
                                            Px(R,C)=min(Px(R,C)+sum(beta*Nhd)+beta*(beta>0)*(i==3),PctMax(i));   %this also kills by a
                                            
                                        else
                                            kk=randi([1,length(gdrange2)]);
                                            beta=gdrange2(kk);
                                            Px(R,C)=min(Px(R,C)+sum(beta*Nhd)+beta*(beta>0)*(i==3),PctMax(i));
                                        end
                                    else
                                        Px(R,C)=0;
                                    end
                                end
                                %                record elevations for all cells inside the range,
                                %                elevation=zero if no popln exists on cell currently
                                if Px(R,C)==0
                                    Ph(R,C)=0;
                                elseif Px(R,C)>0
                                    Ph(R,C)=Hstar(R,C,1);
                                elseif Px(R,C)<0 && Px(R,C)~=-999
                                    Px(R,C)=0;
                                    Ph(R,C)=0;
                                end
                            end
                        end
                    end
                end
            end
            PA(:,:,i)=Px;
            PHA(:,:,i)=Ph;
        end
    end
end


%now updating changes in original plant matrices
P1(:,:,1)=PA(:,:,1);
P2(:,:,1)=PA(:,:,2);
P3(:,:,1)=PA(:,:,3);
P4(:,:,1)=PA(:,:,4);

P1(:,:,2)=PHA(:,:,1);
P2(:,:,2)=PHA(:,:,2);
P3(:,:,2)=PHA(:,:,3);
P4(:,:,2)=PHA(:,:,4);

% Killing Morella too close to shore
Hfilt=imgaussfilt(Hstar(:,:,1),2); %%smoothing the island to get a better sense of where the shorelines are
Iwidth=zeros(size(Hstar,1),1);
WSLtemp=Iwidth;
for R=1:size(Hstar,1)
    Rnow=Hfilt(R,:,1);
    Pnow=P3(R,:,1);
    PHnow=P3(R,:,2);
    if ESL(R)>0 && t<1 || ESL(R)>0 && isnan(t)
        j=find(Rnow(1:ESL(R)-1)>=0,1,'last');   %
        while Rnow(j)>0
            j=j-1;
        end
        j1=ESL(R);
        j2=j;
        if isempty(j)
            j2=j1;
        end
        WSLtemp(R)=j2;
        Iwidth(R)=j1-j2;
        if Iwidth(R)<300%CheckWidth>5
            for m=1:P3Killzone
                C=ESL(R)-m;
                if Pnow(C)>0
                    Pnow(C)=0;
                    PHnow(C)=0;
                end
            end
        elseif Iwidth(R)>=300%CheckWidth<=5
            for m=1:floor(.5*P3Killzone)
                C=ESL(R)-m;
                if Pnow(C)>0
                    Pnow(C)=0;
                    PHnow(C)=0;
                end
            end
        end
        
    elseif ESL(R)>0 && t>1
        j=find(Rnow(1:ESL(R)-1)>=0,1,'last');   %
        while Rnow(j)>0
            j=j-1;
        end
        j1=ESL(R);
        j2=j;
        if isempty(j)
            j2=j1;
        end
        WSLtemp(R)=j2;
        Iwidth(R)=j1-j2;
        if Iwidth(R)<300
            for m=1:min(P3Killzone,ESL(R)-1)
                C=ESL(R)-m;
                if Pnow(C)>0
                    Pnow(C)=0;
                    PHnow(C)=0;
                end
            end
        elseif Iwidth(R)>=300
            for m=1:min(floor(.5*P3Killzone),ESL(R)-1)
                C=ESL(R)-m;
                if Pnow(C)>0
                    Pnow(C)=Pnow(C)-.1*Pnow(C);
                    % %                 if Pnow(C)<0.001
                    % %                     Pnow(C)=0;
                    % %                     PHnow(C)=0;
                    % %                 end
                end
            end
        end
        %     elseif ESL(R)>0 && t>=1
        %         for m=1:P3Killzone
        %             C=ESL(R)-m;
        %             if Pnow(C)>0
        %                 Pnow(C)=Pnow(C)-.1*Pnow(C);
        % %                 if Pnow(C)<0.001
        % %                     Pnow(C)=0;
        % %                     PHnow(C)=0;
        % %                 end
        %             end
        %         end
    end
    P3(R,:,1)=Pnow;
    P3(R,:,2)=PHnow;
end

%new temp plant matrices removing the -999 values to get an accurate total
%of all poplns on cells
Pt1=max(P1(:,:,1),0);
Pt2=max(P2(:,:,1),0);
Pt3=max(P3(:,:,1),0);
Pt4=max(P4(:,:,1),0);

Ptot=Pt1+Pt2+Pt3+Pt4;
%Death by Comp
fprintf('\n')
fprintf('The plants are killing each other!!!')
for i =1:size(Hstar,1)
    for j =1:size(Hstar,2)
        if Ptot(i,j,1)>MasterMax
            Px1=Pt1(i,j,1);
            Px2=Pt2(i,j,1);
            Px3=Pt3(i,j,1);
            Px4=Pt4(i,j,1);
            if Px3==MasterMax
                P1(i,j,1)=0;
                P2(i,j,1)=0;
                P4(i,j,1)=0;
            elseif (Px3 < MasterMax) && (Px3>0)
                k=(MasterMax-Px3)/(Px1+Px2+Px4);
                P1(i,j,1)=Px1*k;
                P2(i,j,1)=Px2*k;
                P4(i,j,1)=Px4*k;
            elseif Px3==0
                %if Px1>=(k/4) && Px2>=(k/4) && Px3>=(k/4) && Px4>=(k/4)
                k=(Px1+Px2+Px4)-MasterMax;
                if Px1>=(k/3) && Px2>=(k/3) && Px4>=(k/3)
                    P1(i,j,1)=Px1-(k/3);
                    P2(i,j,1)=Px2-(k/3);
                    P4(i,j,1)=Px4-(k/3);
                    %Only two are bigger than k/3
                    %%1 2%%
                elseif Px1>=(k/3) && Px2>=(k/3) && Px4<(k/3)
                    P1(i,j,1)=Px1-((k/3)+(((k/3)-Px4)/2));
                    P2(i,j,1)=Px2-((k/3)+(((k/3)-Px4)/2));
                    P4(i,j,1)=0;
                    %%1 4%%
                elseif Px1>=(k/3) && Px2<(k/3) && Px4>=(k/3)
                    P1(i,j,1)=Px1-((k/3)+(((k/3)-Px4)/2));
                    P2(i,j,1)=0;
                    P4(i,j,1)=Px4-((k/3)+(((k/3)-Px2)/2));
                    %%2 4%%
                elseif Px1<(k/3) && Px2>=(k/3) && Px4>=(k/3)
                    P1(i,j,1)=0;
                    P2(i,j,1)=Px2-((k/3)+(((k/3)-Px1)/2));
                    P4(i,j,1)=Px4-((k/3)+(((k/3)-Px1)/2));
                    %Only one is bigger than k/3%
                    %%1%%
                elseif Px1>=(k/3) && Px2<(k/3) && Px4<(k/3)
                    P1(i,j,1)=Px1-((k/3)+((k/3)-Px2)+((k/3)-Px4));
                    P2(i,j,1)=0;
                    P4(i,j,1)=0;
                    %%2%%
                elseif Px1<(k/3) && Px2>=(k/3) && Px4<(k/3)
                    P1(i,j,1)=0;
                    P2(i,j,1)=Px2-((k/3)+((k/3)-Px1)+((k/3)-Px4));
                    P4(i,j,1)=0;
                    %%4%%
                elseif Px1<(k/3) && Px2<(k/3) && Px4>=(k/3)
                    P1(i,j,1)=0;
                    P2(i,j,1)=0;
                    P4(i,j,1)=Px4-((k/3)+((k/3)-Px1)+((k/3)-Px2));
                end
            end
        end
    end
end


%currently run this outside of PlantProcesses
for i=1:size(Hstar,1)
    for j=1:size(Hstar,2)
        PC1(i,j)=P1(i,j,1)*(P1(i,j,1)>0);
        PC2(i,j)=P2(i,j,1)*(P2(i,j,1)>0);
        PC3(i,j)=P3(i,j,1)*(P3(i,j,1)>0);
        PC4(i,j)=P4(i,j,1)*(P4(i,j,1)>0);
        
    end
end
Ptot=PC1+PC2+PC3+PC4;


    function N1=NewNhd(R,C,PAB4)
        if R<=(size(Hstar,1)-1)
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
        if C<=(size(Hstar,2)-1)
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
        N1=[PAB4(Rm,Cm),PAB4(Rm,C),PAB4(Rm,Cp),PAB4(R,Cp),PAB4(Rp,Cp),PAB4(Rp,C),PAB4(Rp,Cm),PAB4(R,Cm)];
        for k=1:8
            if N1(k)==-999
                N1(k)=0;
            elseif N1(k)<0.01
                N1(k)=0;
            end
        end
    end

end
