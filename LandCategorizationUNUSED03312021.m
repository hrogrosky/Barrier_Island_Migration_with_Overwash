function [Hstar,P1,P2,P3,P4]=LandCategorizationUNUSED03312021(Hstar,P1,P2,P3,P4,PlantColumnArray,PlantColumnArray2,ColumnArraySwamp1,ColumnArraySwamp2,INum)
%FOR THE CODE COMPILED ON SEPTEMBER 28 2020 WE DO NOT USE LAND
%CATEGORIZATION

%THIS SUBROUTINE IS STILL CALLED BY THE MAIN CODE - i HAVE BEEN LEAVING IT
%THERE AS A PLACEHOLDER SINCE WE MAY HAVE USE OF IT IN THE FUTURE. IT DOES
%NOT ALTER THE ELEVATION MAP OIN ANY WAY AND THE VALUES STORED IN H(:,:,2)
%ARE NOT UTILIZED ANYWHERE

%KEEPING THE THIRD DIMENSION ON H/HSTAR MOST CERTAINLY SLOWS DOWN THE CODE-
%BUT REMOVING IT REQUIRES A SUBSTANTIAL RE-WRITE

    ROW=size(Hstar,1);
    COLUMN=size(Hstar,2);
    %%% First, we will reset all of Hstar(:,:,2)=0, since it is much more
    %%% difficult to "forcibly" identify the dune fields.
    Hstar(:,:,2)=zeros(ROW,COLUMN);
    
    %%% Next, we will use portions of MarineProcesses to identify which
    %%% portions of the island will be categorized as "beach"
    %%% Hstar(i,j,2)=1
    %%% and portions of SwampProcesses to identify which portions of the
    %%% island will be categorized as "marsh" Hstar(i,j,2)=2.
    
    
    %%%% MARINE PROCESSES %%%%
    % Here, I am making sure that Marine Processes returns PlantColumnArray
    % and PlantColumnArray2, hopefully allowing us to avoid relocating the
    % important cells.
    for i=1:ROW
        for j=PlantColumnArray(i):-1:PlantColumnArray2(i)
            if j>0
                Hstar(i,j,2)=1;
            end
        end
    end
    
    %%%% SWAMP PROCESSES %%%%
    % Here, I am making sure that SwampProcesses returns ColumnArraySwamp1
    % and ColumnArraySwamp2, hopefully allowing us to avoid relocating the
    % important cells (similar to above).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ORIGINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for i=1:ROW
%         for Icount=1:size(ColumnArraySwamp1,2)
%             if ColumnArraySwamp2(i,Icount)~=0
%                 for j=ColumnArraySwamp1(i,Icount):ColumnArraySwamp2(i,Icount)
%                     Hstar(i,j,2)=2;
%                 end
%             end
%         end
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ORIGINAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i=1:size(Hstar,1)
%     ElevCheck1=delta*Hstar(i,:,1)>-.5;  %THIS NEEDS TO BE >-0.5 if not using the smallest matrices
% 
%     Icount=0;%%
%     for j=2:size(Hstar,2)  
%         if ElevCheck1(j)==1 && ElevCheck1(j-1)==0    %if delta*Hstar(i,j)>-.5 && Hstar(i,j-1)<=-.5 %&& ColumnArraySwamp1(i,Icount)==0 
%             Icount=Icount+1;%
%             ColumnArraySwamp1(i,Icount)=j;                                                                                                                    %
%             %break;                                                                                                                                    %
%         end                                                                                                                                            %
%     end                                                                                                                                                %
% end                                                                                                                                                    %
% INum=size(ColumnArraySwamp1,2);
% 
% ColumnArraySwamp2=zeros(size(ColumnArraySwamp1(:,:)));
% for i=1:size(Hstar,1)
%     ElevCheck2=Hstar(i,:,1)*delta>=1;
%     if sum(ColumnArraySwamp1(i,1:INum)~=0)>0
%         for Icount=1:INum
%             if ColumnArraySwamp1(i,Icount)~=0
%                 for j=ColumnArraySwamp1(i,Icount):size(Hstar,2)
%                     if ElevCheck2(j)==1 && ColumnArraySwamp2(i,Icount)==0
%                         ColumnArraySwamp2(i,Icount)=j-1;
%                         break;
%                     elseif sum(ElevCheck2)==0
% 
% 
%     %                 elseif sum(ElevCheck2)==0
%     %                     ColumnArraySwamp1(i,Icount)=0;
%                     end
%                 end
%             end
%         end
%     end
% end


for i=1:size(Hstar,1)
    for MarshCount=1:size(ColumnArraySwamp1,2)
        if ColumnArraySwamp2(i,MarshCount)~=0
            for j=ColumnArraySwamp1(i,MarshCount):ColumnArraySwamp2(i,MarshCount)
                Hstar(i,j,2)=2;
            end
%             if ColumnArraySwamp2==0
%                 Hstar(i,1,2)=0;
%             end
        end
    end
end
    
    
    %%% Lastly, we redeclare which portions of the island are underwater.
   for i=1:ROW
       for j=1:COLUMN
           if Hstar(i,j,1)<0
               Hstar(i,j,2)=-1;
           end
       end
   end

end