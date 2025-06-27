function [Hstar,Centroid]=Qdata03312021(Hstar,P1,P2,P3,P4);%,ESL,WSL,MESL,MWSL,OppositeLocation)
% function [Hstar,Centroid]=Qdata03312021(Hstar,P1,P2,P3,P4,ESL,WSL,MESL,MWSL,OppositeLocation)
%quatitative data for run analysis - this will include:
%       -number of cells with positive elev. and number of slabs above
%           sea level (for island volume)
%       -area of island (m^2)
%       -current location of island centroid 
%       -mean percent cover
[n1 n2]=size(Hstar(:,:,1));
PosCells=0;
IslandArea=0;
Centroid=[0 0];

%centroid
xi=0;yi=0;M=0;
for i=1:n1
    for j=1:n2
        if Hstar(i,j)>=0
        xi=xi+j*Hstar(i,j);
        yi=yi+i*Hstar(i,j);
        M=M+Hstar(i,j);
        end
    end
end

xbar=xi/M;
ybar=yi/M;
Centroid=[xbar ybar]

%calculating surface area

PosCells=0;
for i=1:n1
    RCellCnt=0;
    Rnow=Hstar(i,:,1);
    RCellCnt=sum(Rnow>0);
    PosCells=PosCells+RCellCnt;
end

    %%% First, we will reset all of Hstar(:,:,2)=0, since it is much more
    %%% difficult to "forcibly" identify the dune fields.
%     Hstar(:,:,2)=zeros(n1,n2);
    
% Land Catergories: -1 water, 0 dune/mainland, 1 beach, 2 marsh
% uses MWSL, MESL, WSL and ESL from MarineProcesses to designate
% PtCtr=0;
%     for i=1:n1
%         Rnow=Hstar(i,:,1);
%         CatRnow=Hstar(i,:,2);
%         if ESL(i,1)~=0
%             
%             for j=1:MWSL(i)-1
%                 CatRnow(j)=-1;
%             end
%             for j=MWSL(i):MESL(i)
%                 CatRnow(j)=2;
%             end
%             for j=WSL(i):OppositeLocation(i)-1
%                 CatRnow(j)=0;
%             end
%             if OppositeLocation(i)>0
%                 for j=OppositeLocation(i):ESL(i,1)
%                     CatRnow(j)=1;
%                 end
%             else
%                 for j=WSL(i):ESL(i)
%                     CatRnow(j)=1;
%                 end
%             end
%             for j=ESL(i,1)+1:n2
%                 CatRnow(j)=-1;
%             end
%         else
%             CatRnow(:)=-1;
%         end
%         Hstar(i,:,2)=CatRnow;
%         %gathering x,y,z coords of pos. elev. cells for centroid
% 
%     end


end