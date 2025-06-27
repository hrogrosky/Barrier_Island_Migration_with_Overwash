function [Hstar,MeanBeachWidth,ESL,WSL,MESL,MWSL,OL]=Shoreline03312021(Hstar,delta,L,BchMax)
%this function is just to return the eastern shoreline for plant initialization
%P4 not permitted to grow east of ESL - (E)astern (S)hore (L)ine

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
    end
end

for ii=1:n1
    if ESL(ii)~=0
        Rnow=dH(ii,:);
        MWSL(ii)=find(Rnow>=-0.5,1,'first');
        DC(ii)=find(Rnow==max(Rnow),1,'first');
        SwmpChk=Rnow(MWSL(ii)+1:DC(ii));
        j=find(SwmpChk>1,1,'first');
        if isempty(j)==1
            MESL(ii)=DC(ii);
        else
            MESL(ii)=MWSL(ii)+j;
        end            
    end
end

MBWfactor=(AL~=0);  %calculating changes in beach width - does nothing unless used outside fcn in an image
MeanBeachWidth=(sum(AL(MBWfactor)))/sum(MBWfactor); %same

end