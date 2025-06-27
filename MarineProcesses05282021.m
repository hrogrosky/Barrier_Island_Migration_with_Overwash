function [Hstar,P3,MeanBeachWidth,ESL,WSL,MESL,MWSL,OL,SLRyrs,MigCnt,MigAccel]=MarineProcesses05172021(Hstar,delta,L,P1,P2,P3,P4,BchMax,OW,t,SLRyrs,PCmp,IslandArea,MigCnt,ScaleFactor,MigYr,Ma,MigAccel,MasterMax,SLRvec,BruunParam, Plx_scale)
alpha=1;
AoR=pi/6;
t
T=t/26
% AccelCheck=MigYr*(Ma)^T
% AccelCheck=MigYr*(1+2*SLRaccelrate/SLRinitrate*T);
AccelCheck=BruunParam*(SLRvec(T+2)-SLRvec(T+1));
if floor(AccelCheck)>MigYr
    Mig=floor(AccelCheck);
elseif AccelCheck==0
    Mig=AccelCheck;
else
    Mig=MigYr;
end

nn=20;   %number of rows above and below current row to check plant density when calculating migration %was 5 on 3/13/24

if isnan(t)==1
    Tflag=1;
    t=0;
else
    Tflag=0;
end

test=0;
SLR=0.00635;    %was for Bruun Rule testing - not currently used
Htest=0; %handy in editing to just declare this
Ptest=0; %handy in editing to just declare this

[n1 n2]=size(Hstar(:,:,1));
AL=zeros(n1,1); %AdjascentLength~width of beach
OL=zeros(n1,1); %OppositeLocation~innermost reach of beach/index of first cell>BchMax
ESL=zeros(n1,1);      %EAST shoreline (column value per row;
WSL=zeros(n1,1);    %formerly "WesternShore"
MESL=zeros(n1,1);   %formerly "SlineSwamp"
MWSL=zeros(n1,1);   %formerly "WBndrySwamp"

beta=zeros([n1 1]);             %this gets called but isn't currently used, stores profile slope
Beta=zeros([n1 1]);             %bruun rule beta

R=zeros(n1,1);
DC=zeros(n1,2);


dH=delta*Hstar(:,:,1);  %making new array of Hstar values in terms of meters instead of slabs (ease of use)

for i=1:n1
    Rnow=dH(i,:);
    Pnow=P3(i,:,1);
    PHnow=P3(i,:,2);
    if Rnow(n2)<0
        for j1=n2-1:-1:2                %declaring eastern shoreline by looking for first cell with positive elevation west of a cell with negative elevation
            if Rnow(j1)>=0 && Rnow(j1+1)<=0 && ESL(i)==0
                if Rnow(j1-1)~=0
                    ESL(i)=j1;
                end
            end
            if j1==2 && ESL(i)==0
                ESL(i)=ESL(i)+1*(Rnow(1)>=0);           %special condition for first column - avoid index errors
            end
            if ESL(i)~=0 && ESL(i)~=1                   %if we found a shoreline (that wasn't in col 1) need to declare AL and OL
                DC(i,2)=max(Rnow);                      %will use max of current row if no cell satisfies being greater than BchMax
                DC(i,1)=find(Rnow==DC(i,2),1,'last');
                flag=0;
                while flag==0
                    for j2=j1-1:-1:1                                % looking for dune crest starting with shoreline and moving west
                        if Rnow(j2)>=BchMax                         % if we find a cell >= BchMax
                            OL(i)=j2;                               % then that cell is the dune crest
                            AL(i)=j1-j2;                            % width of the beach is shoreline-dunecrest (indexes)
                            flag=1;
                            break
                        elseif Rnow(j2)<0                           % if we go below water
                            m=1;
                            while j2-m>=1
                                if Rnow(j2-m)>0
                                    j2=j2-m;                        %if we get back above water, change j2 and keep looking for swamp/dune crest
                                    m=1;
                                elseif Rnow(j2-m)<=-0.5
                                    OL(i)=(j2)+find(Rnow(j2+1:j1)==max(Rnow(j2+1:j1)),1,'last');              % make the dune crest the max height of the positive elevation portion of the island
                                    AL(i)=j1-OL(i);                 %width of beach is shoreline index - the index of max height of subaerial island
                                    flag=1;
                                    break
                                elseif m==j2-1
                                    OL(i)=DC(i,1);                  %if we search the rest of the row and don't find swamp/dune crest use max(row) as dune crest
                                    AL(i)=j1-OL(i);
                                    m=j2;
                                    flag=1;
                                else
                                    m=m+1;
                                end
                                
                            end
                            if j2==1 %j2==2                 %if we
                                OL(i)=DC(i,1);              % make the dune crest the max height of the positive elevation portion of the island
                                AL(i)=j1-OL(i);                 %width of beach is shoreline - the max height of subaerial island
                                flag=1;
                            end
                        elseif j2==1
                            OL(i)=DC(i,1);              % make the dune crest the max height of the positive elevation portion of the island
                            AL(i)=j1-OL(i);                 %width of beach is shoreline - the max height of subaerial island
                            flag=1;
                        end
                        if OL(i)>0 && AL(i)>0
                            break
                        end
                    end
                end
            end
            if OL(i)>0 && AL(i)>0
                break
            end
        end
       % if ESL(i)>1                 %if the island has been found in this row
       %     for m=0:AL(i)           %removing P3 from the beach
       %         j=OL(i)+m;
       %         if Pnow(j)>0
       %             Pnow(j)=0;
       %             PHnow(j)=0;
       %         end
       %     end
       %     P3(i,:,1)=Pnow;
       %     P3(i,:,2)=PHnow;
       % elseif ESL(i)==1        %if shoreline is first cell make OL, AL first cell
       %     OL(i)=1;
       %     AL(i)=1;
       %     if Pnow(1)>0
       %         Pnow(1)=0;
       %         PHnow(1)=0;
       %     end
       % end
%         if ESL(i)>0
%             for m=1:P3Killzone
%                 j=ESL(i)-P3Killzone;
%                 if PNow(j)>0
%                     PNow(j)=0;
%                     PHnow(j)=0;
%                 end
%             end
%             P3(i,:,1)=Pnow;
%             P3(i,:,2)=PHnow;
%         end
                    
    end
end

MBWfactor=(AL~=0);  %calculating changes in beach width - does nothing unless used outside fcn in an image
MeanBeachWidth=(sum(AL(MBWfactor)))/sum(MBWfactor); %same



%now calculating the foreshore slope and resetting cells
cnt=0;
DofCi=zeros(n1,1);
xdoc=sym('xdoc');

for i=1:n1
    if  ESL(i)>0
        cnt=cnt+1;
        Rnow=dH(i,:);
        if AL(i)~=0 %if we have a shoreline and the adjascent length is not zero
            beta(i)=tan(Rnow(OL(i))/AL(i)); %calculate the slope of the shoreline
            R(i)=1/beta(i);                     %part of Bruun rule - unused
        elseif AL(i)==0
            beta(i)=.0003;    %if there is a shore that is one cell wide (special condition above) use common shore slope
            R(i)=1/beta(i);
        end
    end
end
cnt=0;
shortR=0; %for removing NaN and infinity slopes - possible with holes and ponds, but unlikely (pretty sure I debugged this issue)
for i=1:n1
    if ESL(i)>0
        cnt=cnt+1;
        shortR(cnt)=R(i);%*(isnan(R(i)==0)); %uncomment this if NaN issue (see above) comes up
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for-Vector-of-Migration-Years%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%below is for using a vector of predetermined years to trigger migration
% YrsPast=sum(SLRyrs(SLRyrs<0));                    %any negative years in vector are summed
% YrCnt=(t/26)+YrsPast;                             %current number of years since last migration
% if YrCnt==max(SLRyrs(SLRyrs>0))                   %if current #yrs is the next number of years in vector of vaues to trigger migration
%     MigChk=1;                                     %trigger migration
%     kk=find(SLRyrs==max(SLRyrs(SLRyrs>0)))        %find that yr value in vector
%     SLRyrs(kk)=-SLRyrs(kk);                       %replace with negative year so it will be summed as years past
% else
%     MigChk=0;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%For_Yearly_Migration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Plx=zeros(n1,1);
if 0==mod(t,26)% && t~=0
    MigChk=1;
    if MigChk==1
        for i=1:n1
            if i<=nn
                PlChkArea=PCmp(1:1+nn,:,1);  %plant cover for rows near current row
              Plx(i)=min(sum(PlChkArea(PlChkArea>=0))/(1),1); %divide 121000 for Smith, 4.4*10^5 for Parramore
             %   Plx(i)=mean(PlChkArea(PlChkArea>=0));
                if sum(dH(i,:)>=0)==0
                    MigCnt(i)=MigCnt(i)+Mig;
                else
                    if Plx(i) >= 0.5
                        MigCnt(i)=MigCnt(i)+round(0.00*Mig);   %reduce migration by 70% if nearby weighted plant cover exceeds 50%
                    elseif Plx(i) < 0.5 
                       MigCnt(i)=MigCnt(i)+round((1-2*Plx(i))*Mig);  %reduce migration by 50% if nearby weighted plant cover exceeds 35% but not 50%
                    end
                end
                if MigCnt(i)>ScaleFactor
                    %                     MigR=1;                                         %set MigR=1 if using scaled version of migration
                    MigR=floor(MigCnt(i)/ScaleFactor);
                    fprintf('\n')
                    fprintf('THE SHORELINE AT ROW %d IS MIGRATING EAST BY %d COLUMNS!',i,MigR)
                    RnowDummy=zeros(1,n2);
                    P1RnowDummy=RnowDummy;
                    P2RnowDummy=RnowDummy;
                    P3RnowDummy=RnowDummy;
                    P4RnowDummy=RnowDummy;
                    RnowDummy(1:n2-MigR)=dH(i,MigR+1:n2);
                    P1RnowDummy(1:n2-MigR)=P1(i,MigR+1:n2);
                    P2RnowDummy(1:n2-MigR)=P2(i,MigR+1:n2);
                    P3RnowDummy(1:n2-MigR)=P3(i,MigR+1:n2);
                    P4RnowDummy(1:n2-MigR)=P4(i,MigR+1:n2);
                    dH(i,1:n2-MigR)=RnowDummy(1:n2-MigR);
                    P1(i,1:n2-MigR)=P1RnowDummy(1:n2-MigR);
                    P2(i,1:n2-MigR)=P2RnowDummy(1:n2-MigR);
                    P3(i,1:n2-MigR)=P3RnowDummy(1:n2-MigR);
                    P4(i,1:n2-MigR)=P4RnowDummy(1:n2-MigR);
                    MigCnt(i)=MigCnt(i)-MigR*ScaleFactor;
                end
            elseif nn<i && i<n1-nn
%                 if i>2000
%                     fprintf('oo')
%                 end
                PlChkArea=PCmp(i-nn:i+nn,:,1);  %plant cover for nows near current row
                Plx(i)=min(sum(PlChkArea(PlChkArea>=0))/(Plx_scale),1); %was 242000 for Smith, 8.85*10^3 for Parramore
              %  Plx(i)=mean(PlChkArea(PlChkArea>=0));
               if sum(dH(i,:)>=0)==0
                    MigCnt(i)=MigCnt(i)+Mig;
                else
                    if Plx(i) >= 0.5
                        MigCnt(i)=MigCnt(i)+round(0.00*Mig);   %reduce migration by 70% if nearby weighted plant cover exceeds 50%
                    elseif Plx(i) < 0.5 
                        MigCnt(i)=MigCnt(i)+round((1-2*Plx(i))*Mig);   %reduce migration by 50% if nearby weighted plant cover exceeds 35% but not 50%
                    end
                end
                if MigCnt(i)>ScaleFactor
                    %                     MigR=1;                                         %set MigR=1 if using scaled version of migration
                    MigR=floor(MigCnt(i)/ScaleFactor);
                    fprintf('\n')
                    fprintf('THE SHORELINE AT ROW %d IS MIGRATING EAST BY %d COLUMNS!',i,MigR)
                    RnowDummy=zeros(1,n2);
                    P1RnowDummy=RnowDummy;
                    P2RnowDummy=RnowDummy;
                    P3RnowDummy=RnowDummy;
                    P4RnowDummy=RnowDummy;
                    RnowDummy(1:n2-MigR)=dH(i,MigR+1:n2);
                    P1RnowDummy(1:n2-MigR)=P1(i,MigR+1:n2);
                    P2RnowDummy(1:n2-MigR)=P2(i,MigR+1:n2);
                    P3RnowDummy(1:n2-MigR)=P3(i,MigR+1:n2);
                    P4RnowDummy(1:n2-MigR)=P4(i,MigR+1:n2);
                    dH(i,1:n2-MigR)=RnowDummy(1:n2-MigR);
                    P1(i,1:n2-MigR)=P1RnowDummy(1:n2-MigR);
                    P2(i,1:n2-MigR)=P2RnowDummy(1:n2-MigR);
                    P3(i,1:n2-MigR)=P3RnowDummy(1:n2-MigR);
                    P4(i,1:n2-MigR)=P4RnowDummy(1:n2-MigR);
                    MigCnt(i)=MigCnt(i)-MigR*ScaleFactor;
                end
            elseif i>=n1-nn
                PlChkArea=PCmp(n1-nn:n1,:,1);  %plant cover for nows near current row
               Plx(i)=min(sum(PlChkArea(PlChkArea>=0))/(Plx_scale),1);
              %  Plx(i)=mean(PlChkArea(PlChkArea>=0));
               if sum(dH(i,:)>=0)==0
                    MigCnt(i)=MigCnt(i)+Mig;
                else
                   if Plx(i) >= 0.5
                        MigCnt(i)=MigCnt(i)+round(0.00*Mig);   %reduce migration by 70% if nearby weighted plant cover exceeds 50%
                    elseif Plx(i) < 0.5 
                       MigCnt(i)=MigCnt(i)+round((1-2*Plx(i))*Mig);   %reduce migration by 50% if nearby weighted plant cover exceeds 35% but not 50%
                    end
                end
                if MigCnt(i)>ScaleFactor
                    %                     MigR=1;                                         %set MigR=1 if using scaled version of migration
                    MigR=floor(MigCnt(i)/ScaleFactor);
                    fprintf('\n')
                    fprintf('THE SHORELINE AT ROW %d IS MIGRATING EAST BY %d COLUMNS!',i,MigR)
                    RnowDummy=zeros(1,n2);
                    P1RnowDummy=RnowDummy;
                    P2RnowDummy=RnowDummy;
                    P3RnowDummy=RnowDummy;
                    P4RnowDummy=RnowDummy;
                    RnowDummy(1:n2-MigR)=dH(i,MigR+1:n2);
                    P1RnowDummy(1:n2-MigR)=P1(i,MigR+1:n2);
                    P2RnowDummy(1:n2-MigR)=P2(i,MigR+1:n2);
                    P3RnowDummy(1:n2-MigR)=P3(i,MigR+1:n2);
                    P4RnowDummy(1:n2-MigR)=P4(i,MigR+1:n2);
                    dH(i,1:n2-MigR)=RnowDummy(1:n2-MigR);
                    P1(i,1:n2-MigR)=P1RnowDummy(1:n2-MigR);
                    P2(i,1:n2-MigR)=P2RnowDummy(1:n2-MigR);
                    P3(i,1:n2-MigR)=P3RnowDummy(1:n2-MigR);
                    P4(i,1:n2-MigR)=P4RnowDummy(1:n2-MigR);
                    MigCnt(i)=MigCnt(i)-MigR*ScaleFactor;
                end
            end
        end
    end
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     *****RESETTING EQUILIBRIUM PROFILE OF SHORELINE HAS BEEN REMOVED FOR THIS VERSION*****
    %      ***see Fix2MarineProcesses02202021 for most recent version of resetting profile***
    
    dH=round(dH,1); %need to round to 1 dec. place if redeclared shoreline
    Hstar(:,:,1)=(1/delta)*dH;  %convert back to slabs when redeclaring Hstar
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Tflag==1
        t=NaN;  %reset t if this is initialization MP so it doesn't throw of loop in MainCode
    end
    
end