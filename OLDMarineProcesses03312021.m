function [Hstar,PlantColumnArray,PlantColumnArray2,P3,flag,MeanBeachWidth]=OLDMarineProcesses03312021(Hstar,delta,L,flag,PC,P3,BchMax,BchW,t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MARINE PROCESSES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This subroutine will consider the periodicity of oceanic waves and how
%   they can change the landscape of coastal islands. We will construct an
%   "Equilibrium profile" which is a smooth surface maintained by the runup
%   and swash of these waves.                  

%   To run this file, you will need to specify:
%      Hstar  - elevation matrix that will be updated
%      delta  - the thickness of each slab
%        L    - the length of each slab
%   The subroutine will then return Hstar after it has been updated%
%   according to the marine processes.   
SLR=0.00635;
Htest=0; %handy in editing to just declare this
Ptest=0; %handy in editing to just declare this
% BchMax=2;   %greatest elevation that the beach spreads inland to

Li=size(Hstar,1);
Wi=size(Hstar,2);%
AdjascentLength=zeros([Li 1]); %width of beach
OppositeLocation=zeros([Li 1]); %innermost reach of beach
ColumnArray=zeros([Li 1]);      %shoreline (column value per row)
NSRowArray=zeros([1 Wi]);       %noth shoreline
SSRowArray=zeros([1 Wi]);       %south shoreline

PlantColumnArray=zeros([Li 1]); %P3 isn't allowed to grow on the beach - this is to check for that
PlantColumnArray2=zeros([Li 1]);
LongTransect=zeros([Wi 1]);     %width of beach again...hmmm
beta=zeros([Li 1]);             %this gets called but isn't currently used, stores profile slope
Beta=zeros([Li 1]);
dSline=zeros([Li 1]);

R=zeros(size(dH,1),1);

% fprintf('I am running MP! ')
                                                                                                                                                       
dH=delta*Hstar;                                                                                                                                                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the slope and R2% to help determine the Equilibrium Profile %                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                   %
% First, we find all of the rows which have shoreline and record the column
% in which we find the shoreline. We store this is ColumnArray.
for i=1:size(dH,1)  
    if dH(i,size(dH,2),1)<0
        for j=size(dH,2)-1:-1:1
            if dH(i,j,1)>=0 && dH(i,j+1,1)<0 && ColumnArray(i)==0
                ColumnArray(i)=j;
                break;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% ALL OF THIS IS TO DECLARE NORTH AND SOUTH SHORES - CURRENTLY NOT ERODING THIS WAY, BUT MAY IN THE FUTURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%same, but for columns - tips of islands. Row Array will have two rows, for
%north and south shore.
% 
% for j=1:size(dH,2)
%     ColNow=dH(:,j,1);
%     ColCheck=(dH(:,j,1)>=0);
%     NS_INum=0;
%     if sum(ColCheck)~=0
%         for i=2:size(dH,1)
%             if ColCheck(i-1,1)==0 && ColCheck(i,1)==1
%                 if j==1
%                     if i~=size(dH,1)
%                         if dH(i+1,j,1)>=0 || dH(i,j+1,1)>=0 || dH(i,j+2,1)>=0
%                             NS_INum=NS_INum+1;
%                             NSRowArray(NS_INum, j)=i;
%                         end
%                     else
%                         if dH(i,j,1)>=0 || dH(i,j+1,1)>=0 || dH(i,j+2,1)>=0
%                             NS_INum=NS_INum+1;
%                             NSRowArray(NS_INum, j)=i;
%                         end
%                     end
%                 elseif j==size(dH,2)
%                     if i~=size(dH,1)
%                         if dH(i+1,j,1)>=0 || dH(i,j-1,1)>=0 || dH(i,j-2,1)>=0
%                             NS_INum=NS_INum+1;
%                             NSRowArray(NS_INum, j)=i;
%                         end
%                     else
%                         if dH(i,j,1)>=0 || dH(i,j-1,1)>=0 || dH(i,j-2,1)>=0
%                             NS_INum=NS_INum+1;
%                             NSRowArray(NS_INum, j)=i;
%                         end
%                     end
%                 else
%                     if i~=size(dH,1)
%                         if dH(i+1,j,1)>=0 || dH(i,j+1,1)>=0 || dH(i,j-1,1)>=0
%                             NS_INum=NS_INum+1;
%                             NSRowArray(NS_INum, j)=i;
%                         end
%                     else
%                         if dH(i,j,1)>=0 || dH(i,j+1,1)>=0 || dH(i,j-1,1)>=0
%                             NS_INum=NS_INum+1;
%                             NSRowArray(NS_INum, j)=i;
%                         end
%                     end
%                 end       
%             end
%         end
%         SS_INum=0;
%         for i=size(dH,1)-1:-1:1
%             if ColCheck(i+1,1)==0 && ColCheck(i,1)==1
%                 if dH(i-1,j,1)>=0 || dH(i,j+1,1)>=0 || dH(i,j-1,1)>=0
%                     SS_INum=SS_INum+1;
%                     SSRowArray(SS_INum,j)=i;
%                     break
%                 end
%             end
%         end
%     end
% end
% 
% 
% NorShoreLength=zeros([size(NSRowArray,1) Wi]);
% NS_OppLoc=zeros([size(NSRowArray,1) Wi]);
% SouShoreLength=zeros([size(SSRowArray,1) Wi]);
% SS_OppLoc=zeros([size(SSRowArray,1) Wi]);
% 
% %Declaring North Shore Beach Region
% for j=1:size(dH,2)
%     for INum=1:size(NSRowArray,1)
% 
%         if NSRowArray(INum,j)~=0
%             ColNow=dH(:,j,1);
%             i=NSRowArray(INum,j);
%         else
%             i=-1;
%         end
%         if NSRowArray(INum,j)>0 && i>0
%             if NSRowArray(INum,j)==size(dH,2)
%                 NorShoreLength(INum,j)=1;
%                 NS_OppLoc(INum,j)=i;
%             else
%                 a=i;
%                 %zz=find(ColNow>=0,1,'first');
%                 for mm=a+1:size(dH,1)
%                     yy=mm;
%                     if ColNow(mm,1)<0
%                         break
%                     end 
%                 end
%                 for a=i+1:yy-1
%                     if ColNow(a,1)*delta<=BchMax  && ColNow(a+1,1)*delta>BchMax
%                         NorShoreLength(INum,j)=(a-i)*L;
%                         NS_OppLoc(INum,j)=a;
%                         break
%                     end
%                     if a==yy-1 && NorShoreLength(INum,j)==0 
%                         NorShoreLength(INum,j)=12;                               %chose a shore width of 2 for when the BchMax is never reached (for most of my toy I's each col has two layers of 0's so we need at least two to get to positive elevation)
%                         NS_OppLoc(INum,j)=i+12;
%                     end
%                     if NorShoreLength(INum,j)>=(yy-i)+1
%                         NorShoreLength(INum,j)=0;                               %if we go past the length of the island (hit water again) we call no shore for that column (this should never happen if above works)
%                         NS_OppLoc(INum,j)=0;
%                     end
%                 end
%             end
%             if NorShoreLength(INum,j)~=0 && NS_OppLoc(INum,j)~=0
%                 for kk=NSRowArray(INum,j):NS_OppLoc(INum,j)
%                     if P3(kk,j,1)>0
%                         P3(kk,j,1)=0;
%                     end
%                 end
%             end
%         end
%     end
%     %Declaring South Shore Beach Reqions
%     for INum=1:size(SSRowArray,1) 
%         if SSRowArray(INum,j)~=0
%             ColNow=dH(:,j,1);
%             i=SSRowArray(INum,j);
%         else
%             i=-1;
%         end
%         if SSRowArray(INum,j)>0 && i>0
%             if SSRowArray(INum,j)==1
%                 SouShoreLength(INum,j)=1;
%                 SS_OppLoc(INum,j)=1;
%             else
%                 if SSRowArray(INum,j)>0 && i>0
%                     b=i;
%                     %zz=find(ColNow>=0,1,'last');
%                     for nn=b-1:-1:2
%                         yy=nn;
%                         if ColNow(nn,1)<0
%                             break
%                         end
%                     end
%                     for b=i-1:-1:yy
%                         if ColNow(b,1)*delta<=BchMax && ColNow(b-1,1)*delta>BchMax
%                             SouShoreLength(INum,j)=(i-b)*L;
%                             SS_OppLoc(INum,j)=b;
%                             break
%                         end
%                         if a==yy 
%                             SouShoreLength(INum,j)=12;                               %chose a shore width of 2 for when the BchMax is never reached (for most of my toy I's each col has two layers of 0's so we need at least two to get to positive elevation)
%                             SS_OppLoc(INum,j)=i-12;
%                         end
%                         if SouShoreLength(INum,j)>=(i-yy)+1
%                             SouShoreLength(INum,j)=0;
%                             SS_OppLoc(INum,j)=0;
%                         end
%                     end
%                     if SouShoreLength(INum,j)~=0 && SS_OppLoc(INum,j)~=0
%                         for kk=SS_OppLoc(INum,j):SSRowArray(INum,j)
%                             if P3(kk,j,1)>0
%                                 P3(kk,j,1)=0;
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%declaring marsh shore
                                                                                                                                                       %
% First, we find all of the rows which have shoreline and record the column                                                                            %
% in which we find the shoreline. We store this is SLineSwamp.  
SlineSwamp= zeros ([Li 1]);
WesternShore=zeros([Li 1]);
for i=1:size(dH,1)                                                     
    for j=1:size(dH,2)
        if dH(i,j)>=0 && WesternShore(i)==0
            WesternShore(i)=j;
        end
        if dH(i,j)>=1 && SlineSwamp(i)==0
            SlineSwamp(i)=j;
            break;
        end
    end
end
  

%Next, we find the first cell, which is to the west of
%H(i,ColumArraySwamp1(i)) for each relevant value for i. We store this in
%SlineSwamp.
WBndrySwamp=zeros ([Li 1]);%Column value for swamp shoreline location%
for i=1:size(dH,1) 
    if SlineSwamp(i)~=0%
        for j=SlineSwamp(i):-1:2                                                                                                                              %
            if (dH(i,j-1)*delta)<-.5 && WBndrySwamp(i)==0 %&& dH(i,j-1)>=0                                                                                %
                WBndrySwamp(i)=j;                                                                                                                    %
                break;                                                                                                                                     %
            end                                                                                                                                            %
        end
    end%
end  
%now have east and west boundary for swamp

% LandFlag=0;
for i=1:size(dH,1)                                                                                                                                  %
    if ColumnArray(i)~=0
        for j=ColumnArray(i)-1:-1:2
            if ColumnArray(i)-j<BchW
                if dH(i,j-1)>=BchMax && dH(i,j)>0%was dH(i,j-1)<=BchMax before 8/20 - changed so shore would be slope from 2 and not whatever height happens to be east of 2
                    OppositeLocation(i)=j;
                    AdjascentLength(i)=L*(ColumnArray(i)-j);
                    break;
                elseif dH(i,j-1)<0
                    OppositeLocation(i)=j;
                    AdjascentLength(i)=L*(ColumnArray(i)-j);
                    break;                                
                end
            elseif ColumnArray(i)-j==BchW
                OppositeLocation(i)=j;
                AdjascentLength(i)=L*(ColumnArray(i)-j);
                break
            else
                break
            end
        end
            %if there is a single cell greater than zero followed by one
            %less than zero we skip over it and start again with the next cell greater than 0
        while OppositeLocation(i)==ColumnArray(i) && OppositeLocation(i)>2
            if dH(i,ColumnArray(i)-2)>0
                ColumnArray(i)=ColumnArray(i)-2;
                for j=ColumnArray(i)-2:-1:2
                    if ColumnArray(i)-j<BchW
                        if dH(i,j-1)>=BchMax && dH(i,j)>0%was dH(i,j-1)<=BchMax before 8/20 - changed so shore would be slope from 2 and not whatever height happens to be east of 2
                            OppositeLocation(i)=j;
                            AdjascentLength(i)=L*(ColumnArray(i)-j);
                            break;
                        elseif dH(i,j-1)<0
                            OppositeLocation(i)=j;
                            AdjascentLength(i)=L*(ColumnArray(i)-j);
                            break;                                
                        end
                    elseif ColumnArray(i)-j==BchW
                        OppositeLocation(i)=j;
                        AdjascentLength(i)=L*(ColumnArray(i)-j);
                        break
                    else
                        break
                    end
                end
            else %single positive elevation in a whole row replace with mean(neighbors to east and west)
                dH(i,ColumnArray(i))=mean([dH(i,ColumnArray(i)-1) dH(i,ColumnArray(i)+1)]);
                ColumnArray(i)=0;
                OppositeLocation(i)=0;
                AdjascentLength(i)=0;
            end
                
        end
        if OppositeLocation(i)==0 && ColumnArray(i)~=0
            OppositeLocation(i)=-999;
        
        end
        if OppositeLocation(i)>0 && ColumnArray(i)~=0
            for kk=OppositeLocation(i):ColumnArray(i)
                if P3(i,kk,1)>0
                    P3(i,kk,1)=0;
                end
            end
        end
    end
end
MBWfactor=(AdjascentLength~=0);
MeanBeachWidth=(sum(AdjascentLength(MBWfactor)))/sum(MBWfactor);
    
% if t>0
for i=1:size(dH,1)
    if OppositeLocation(i)>0
        %equilibrium foreshore slope
        beta(i)=atan((dH(i,OppositeLocation(i),1))/AdjascentLength(i));
        Beta(i)=tan((dH(i,OppositeLocation(i),1))/AdjascentLength(i));
        dSline(i)=SLR/Beta(i);
    end
end
% end  
R=floor(1./tan(beta));

for i=1:size(dH,1)
    if OppositeLocation(i)>0 && ColumnArray(i)>0  
        for k=AdjascentLength(i):-1:1                                                                   %
            dH(i,ColumnArray(i)-k)=round(dH(i,OppositeLocation(i))-((dH(i,OppositeLocation(i))/AdjascentLength(i))*(AdjascentLength(i)-k)),1);
%           dH(i,ColumnArray(i)-k)=round((Hstar(i,OppositeLocation(i))/AdjascentLength(i))*(AdjascentLength(i)-k),1);
        end
        for k=-1:-1:-AdjascentLength(i) %washing back out into ocean for distance equal to the beach width of that row
            if dH(i, ColumnArray(i)+k)~=0
                dH(i,ColumnArray(i)-k)=-dH(i,ColumnArray(i)+k);
            end
        end
    elseif OppositeLocation(i)<0 && ColumnArray(i)>0
        %OVERWASH
        IslandTransect=dH(i,WesternShore(i):ColumnArray(i),1);
        DcIndex=find(IslandTransect==max(IslandTransect),1,'last');
        prob=rand;
        if DcIndex==length(IslandTransect)                          %Trying to change this so only rows touching water get overwashed 
            for kk=length(IslandTransect)-1:-1:1                    %this overwash doesn't go straight back, but to cell diagonall back from host
                if 0.5 < prob                                                                        %50% chance if whole row is 0
                    dH(i,WesternShore(i)+kk,1)=dH(i,WesternShore(i)+kk,1)-(1*delta);
                    if i>1 && kk>=2 && dH(i-1, WesternShore(i)+(kk-1),1)<0                        %if row above top row of island
                        dH(i-1, WesternShore(i)+(kk-1),1)=dH(i-1, WesternShore(i)+(kk-1),1)+(1*delta);
                    elseif i < size(dH,1) && kk>=2 && dH(i+1, WesternShore(i)+(kk-1),1)<0  &&  (WesternShore(i)+(kk-1))>0            %if bottom row
                        dH(i+1, WesternShore(i)+(kk-1),1)=dH(i+1, WesternShore(i)+(kk-1),1)+(1*delta);
                    end
                end
            end
        else
            for kk=length(IslandTransect):-1:DcIndex
                CellShift=(1*delta)*((1-PC(1,kk))<rand);     %added ((1-PC)<rand) to factor in plants to this erosion
                if CellShift>0
                    dH(i,WesternShore(i)+kk,1)=dH(i,WesternShore(i)+kk,1)-(1*delta)*((1-PC)<rand); 
                    dH(i,WesternShore(i)+(kk-1),1)=dH(i,WesternShore(i)+(kk-1),1)-(1*delta);        %trying moving each short-row overwash cell by one
                end
            end
        end          
    end
end     
Hstar=(1/delta)*dH;
end