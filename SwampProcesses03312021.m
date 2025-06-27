function [Hstar,ColumnArraySwamp1,ColumnArraySwamp2,AdjascentLengthSwamp,OppositeLocationSwamp,flag,INum]=SwampProcesses03312021(Hstar,delta,L,flag,PC)


%   To run this file, you will need to specify:
%      Hstar  - elevation matrix that will be updated
%      delta  - the thickness of each slab
%        L    - the length of each slab 
%       flag  - boolean variable to indicate the need to avalanche

[n1 n2]=size(Hstar(:,:,1));
AdjascentLengthSwamp=zeros([n1 1]);     %width of swamp
OppositeLocationSwamp=zeros([n1 1]);    %farthest away from island that is still swamp
ColumnArraySwamp1=zeros([n1 1]);        %eastern boundary (island west coast)
ColumnArraySwamp2=zeros([n1 1]);        %western boundary (out into backbarrier/bay)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declare Marsh Region & Poll Eastern Boundary Cells to be Deposited at Western Boundary %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                                                                           
% First, we find all of the rows which have a marsh "shoreline" (western boundary for marsh) and record the column
% in which we find the shoreline. We store this is ColumnArraySwamp1
for i=1:size(Hstar,1)
    ElevCheck1=delta*Hstar(i,:,1)>=-.5;  
    Icount=0;%%
    for j=2:size(Hstar,2)  
        if ElevCheck1(j)==1 && ElevCheck1(j-1)==0
            Icount=Icount+1;
            ColumnArraySwamp1(i,Icount)=j;
        end
    end
end
INum=size(ColumnArraySwamp1,2);                 %tracking number of land masses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % BELOW - PIECES NOT CURRENTLY USED BUT MAY BE USEFUL IN THE FUTURE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % FINDING THE EASTERN BORDER OF SWAMP %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Greg's notes on this piece:
%Next, we find the first cell, which is to the east of
%H(i,ColumArraySwamp1(i)) for each relevant value for i. We store this in
%ColumArraySwamp2.
% % % 
% % % ColumnArraySwamp2=zeros(size(ColumnArraySwamp1(:,:)));
% % % for i=1:size(Hstar,1)
% % %     ElevCheck2=Hstar(i,:,1)*delta>=.5;
% % %     Schk1=sum(ColumnArraySwamp1(i,1:INum));
% % %     Schk2=sum(ElevCheck2);
% % %     if Schk1~=0>0
% % %         for Icount=1:INum
% % %             if ColumnArraySwamp1(i,Icount)~=0
% % %                 for j=ColumnArraySwamp1(i,Icount):size(Hstar,2)
% % %                     if ElevCheck2(j)==1 && ColumnArraySwamp2(i,Icount)==0
% % %                         ColumnArraySwamp2(i,Icount)=j-1;
% % %                         break;
% % %                     elseif Schk2==0
% % %                     end
% % %                 end
% % %             end
% % %         end
% % %     end
% % % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % DECLARING THE REGION BETWEEN THE TWO BOUNDARIES AS 'SWAMP' STORED IN Hstar(:,:,2) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%(Greg's notes on this piece:)
% % % %Now, we have an eastern boundary and a western boundary for our marsh for
% % % %each row on the island. We must be careful in the next step! Because of
% % % %the nature of the tips of the islands, we must acknowledge that is it
% % % %possible that for some value of i, ColumnArraySwamp1(i)~=0 but
% % % %ColumnArraySwamp2(i)=0! We must be careful to avoid computation error.
% % % 
% % % %Next, we will declare which region is our "initial" swamp region.
% % % 
% % % for i=1:size(Hstar,1)
% % %     if i==52
% % %         RowCheck=Hstar(i,:,1);
% % %     end
% % %     for Icount=1:INum
% % %         if ColumnArraySwamp2(i,Icount)~=0
% % %             for j=ColumnArraySwamp1(i,Icount):ColumnArraySwamp2(i,Icount)
% % %                 Hstar(i,j,2)=2;
% % %             end
% % %         end
% % %     end
% % % end


% % %Next, we will look at Hstar(i,ColumnArraySwamp2(i)+1) (the portion of the
% % %island which is NOT in the marsh but is directly adjacent to the marsh)
% % %and give it a 20% probability to be deposited at
% % %Hstar(i,ColumnArraySwamp1(i)-1) (the first portion of the island which has
% % %a non-positive elevation that is directly adjascent to the marsh).
% % %This will produce a "flattening" of the marsh which will intentionally
% % %migrate westward. It will also potentially produce eastward migration of
% % %the marsh (due to local elevation changes).
% % 
% % if rand<=.2
% % for i=1:size(Hstar,1)
% %     for Icount=1:INum
% %         if ColumnArraySwamp2(i,Icount)~=0 && ColumnArraySwamp1(i,Icount)>1
% %             if (ColumnArraySwamp1(i,Icount)-1)>=1
% % %                 if rand<=.2  %make this a main code parameter and pull out
% %                     Hstar(i,ColumnArraySwamp1(i,Icount)-1)=Hstar(i,ColumnArraySwamp1(i,Icount)-1)+1;
% %                     Hstar(i,ColumnArraySwamp2(i,Icount)+1)=Hstar(i,ColumnArraySwamp2(i,Icount)+1)-1;
% %                     flag=1;
% %                 end
% %             end
% %         end
% %     end
% % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % END - UNUSED PIECES % % % % % % % % % % % % % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%probability of the swamp migration - the first piece of land above marsh
%elevation range is removed from current cell and moved one cell west

if rand<=0.10           %<---probability of marsh erosion
    for i=1:n1
        if ColumnArraySwamp1(i,1)~=0
            Rnow=Hstar(i,:,1)*delta;
            CAS1=ColumnArraySwamp1(i,1);            %beginning at the swamp boundary
            j=find(Rnow(CAS1+1:n2)>=1,1,'first');   %find first cell above max swamp elev.
            if isempty(j)==0                        %if such a cell ^ exists
                Rnow(CAS1+j)=Rnow(CAS1+j)-1*delta;  %remove a slab
                Rnow(CAS1-1)=Rnow(CAS1-1)+1*delta;  %deposit westward
            end
            Hstar(i,:,1)=Rnow*(1/delta);            %update Hstar
            flag=1;                                 %since cells have been moved, trigger flag for avalanching
        end
    end
end
ColumnArraySwamp2=zeros(n1,INum);                   %need to output ColumnArraySwamp 2 for UNUSEDLandCategorization fcn
end