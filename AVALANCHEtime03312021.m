function [Hstar,flag,CellCt] = AVALANCHEtime03312021(Hstar,delta,L,flag,PC,t)
%% This version of avalanche:           
%   ~weights the probability of avalanching to favor the direction of the 
%       greatest violation of the angle of repose.
%   ~takes t as input to track if avalanching is required over whole domain
%       (once per 5 years)
%   ~most recent edit code was ATTEMPT2newAVALANCHEtime0928
% 
%AoR=pi/6;
 AoR=pi/9; %changed this from line above on 8/6/23

%%
% tic

Hstarflag=Hstar;
CellCt=0;
flag=0;
FLAG=1;
alpha=1;
% fprintf('I am running AV! ')
Nhd=4;                      %size of neighborhood check (4 or 8)
beta1=zeros(1,Nhd);   %angle of repose
beta2=zeros(1,Nhd);   %
check1=zeros(1,Nhd);
check2=zeros(1,Nhd);

%once per 5 year we check every cell in domain, both above and  below sea level

%if (0==mod(t,26)) || exist('t','var')==0 
for R=1:size(Hstar,1)
    for C=1:size(Hstar,2)
        Erosioncheck=1;
        FLAG=1;
        while FLAG==1
            RAND=rand;
            %FLAG=1;
            %if Hstar(R,C)>0
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

                Hx=Hstar(R,C);                                              %current cell value
                N=[Hstar(Rm,C),Hstar(R,Cp),Hstar(Rp,C),Hstar(R,Cm)];        %neighborhood
                
                beta1=atan(((Hx-N)*delta)/L);
                check1=(beta1>=(AoR));
                %Avdir=find(beta1>=AoR);
                beta2=alpha*((beta1/(AoR))*(1-PC(R,C)));
                check2=(RAND<beta2);
                
                if Hstar(R,C)<2 %doesn't run avalance if cell height is less than 2
                    break
                end


                if sum(check2)~=0
                    while  sum(check1)~=0 && sum(check2)~=0  %if I use a while loop, there's no point to using probability
                        %posbeta1=beta1;
                        %posbeta1(find(beta1<0))=0;
                        %Nprob=posbeta1/sum(posbeta1);
                        AVbeta1=zeros(1,4);
                        AVbeta1(check1)=beta1(check1);
                        AVbeta1(check2)=AVbeta1(check2);
                        Nprob=AVbeta1/sum(AVbeta1);
                        prob=rand;
                        Hcheck=Hx;
                        Hx=Hx-1*check1(1)*check2(1)*(prob<Nprob(1))-1*check1(2)*check2(2)*(prob<Nprob(2))-1*check1(3)*check2(3)*(prob<Nprob(3))-1*check1(4)*check2(4)*(prob<Nprob(4));
%                         if Hx-Hcheck==0
%                             FLAG=0;
%                             break
%                         end
                        N(1)=N(1)+check1(1)*check2(1)*(prob<Nprob(1));
                        N(2)=N(2)+check1(2)*check2(2)*(prob<Nprob(2));
                        N(3)=N(3)+check1(3)*check2(3)*(prob<Nprob(3));
                        N(4)=N(4)+check1(4)*check2(4)*(prob<Nprob(4));
                        beta1=newbeta1(N,Hx);
                        beta2=newbeta2(N,beta1);
                        check1=newcheck1(N,beta1);
                        check2=newcheck2(N,beta2);
                        %flag=1;
                        if Hx-Hcheck~=0
                            CellCt=CellCt+1;
                            flag=1;
                        end
                        if sum(check1)==0 
                            Hstar(R,C)=Hx;
                            Hstar(Rm,C)=N(1);
                            Hstar(R,Cp)=N(2);
                            Hstar(Rp,C)=N(3);
                            Hstar(R,Cm)=N(4);
                            %fprintf('Changes have been made')
                            FLAG=0;
                            if Hx >=0
                                break
                            end
                        end
                        if check1*check2'==0
                            break
                        end
                    end
                    if sum(check1)==0
                        FLAG=0;
                        %no avalanching needed
                    end
                elseif sum(check2)==0 && sum(check1)~=0
                   FLAG=0;
                   %fprintf('Too many plants!')
                elseif sum(check1)==0
                    FLAG=0;     %only other case should be if both sum to 0
                end

%             else%{if Hstar(R,C)<=0}
%                 Erosioncheck=0;
%                 FLAG=0;
%             end
            FLAG=0;
        end
    end
end
%for other instances of avalanche we only check subaerial portions of island
%else
%     ISLANDindex=Hstar(:,:,1)>0; %all cells with land
%     Hstarflag=cumsum(ISLANDindex,2) == 1 & ISLANDindex; %this outputs an array same size as H but with a 1 in the first cell>0 
%     IslandCheck=sum(Hstarflag');    %outputs a row vector with 0~ no col>0, 1~col>0 (ie yes/no land in that row)
%     COLindex1 = Hstarflag*(1:size(Hstar,2))';   %the index of the first positive value in a row
%     COLindex2 = zeros(size(COLindex1));
%     for ii=1:length(COLindex1)
%         if ii == 101
%             fprintf('yay')
%         end
%         
%         if COLindex1(ii)~=0
%             RowNow=ISLANDindex(ii,:);
%             Icount=0;
%             for jj = COLindex1(ii):1:size(RowNow)
%                 if RowNow(jj)==1 && RowNow(jj+1)==0
%                     Icount=Icount+1;
%                     COLindex2(ii,Icount)=jj;
%                 end
%             end
%         end
%     end
IslandColumnArray1=zeros([size(Hstar,1) 1]);

for i=1:size(Hstar,1)
    ElevCheck1=Hstar(i,:,1)>=0;  %first cell above water
    Icount=0;%%
    for j=2:size(Hstar,2)  
        if ElevCheck1(j)==1 && ElevCheck1(j-1)==0    %if delta*Hstar(i,j)>-.5 && Hstar(i,j-1)<=-.5 %&& ColumnArraySwamp1(i,Icount)==0 
            Icount=Icount+1;%
            IslandColumnArray1(i,Icount)=j;                                                                                                                    %
            %break;                                                                                                                                    %
        end                                                                                                                                            %
    end                                                                                                                                                %
end  
INum=size(IslandColumnArray1,2);
IslandColumnArray2=zeros(size(IslandColumnArray1(:,:)));
for i=1:size(Hstar,1)
    ElevCheck2=Hstar(i,:,1)<=0;
    if sum(IslandColumnArray1(i,1:INum)~=0)>0
        for Icount=1:INum
            if IslandColumnArray1(i,Icount)~=0
                for j=IslandColumnArray1(i,Icount):size(Hstar,2)
                    if ElevCheck2(j)==1 && IslandColumnArray2(i,Icount)==0
                        IslandColumnArray2(i,Icount)=j;
                        break;
    %                 elseif sum(ElevCheck2)==0
    %                     ColumnArraySwamp1(i,Icount)=0;
                    end
                end
            end
        end
    end
end

for R=1:size(Hstar,1)
    for Icount=1:INum
        if IslandColumnArray2(R,Icount)~=0
            for C=IslandColumnArray1(R,Icount):IslandColumnArray2(R,Icount)
                Erosioncheck=1;
                            FLAG=1;
                            while FLAG==1
                                RAND=rand;
                                %FLAG=1;
                                if Hstar(R,C)>0
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

                                    Hx=Hstar(R,C);                                              %current cell value
                                    N=[Hstar(Rm,C),Hstar(R,Cp),Hstar(Rp,C),Hstar(R,Cm)];        %neighborhood

                                    beta1=atan(((Hx-N)*delta)/L);
                                    check1=(beta1>=(AoR));
                                    %Avdir=find(beta1>=AoR);
                                    beta2=alpha*((beta1/(AoR))*(1-PC(R,C)));
                                    check2=(RAND<beta2);


                    %                 if sum(check1)~=0
                    %                     if sum(check2)~=0
                    %                     end
                    %                 end
                                        if R==241
                                            if C==492
%                                                 fprintf('STOP!')
                                            end                                            
                                        end

                                    if sum(check2)~=0
                                        if check1*check2'~=0
                                            while  sum(check1)~=0 && sum(check2)~=0  %if I use a while loop, there's no point to using probability
                                                %posbeta1=beta1;
                                                %posbeta1(find(beta1<0))=0;
                                                %Nprob=posbeta1/sum(posbeta1);
                                                AVbeta1=zeros(1,4);
                                                AVbeta1(check1)=beta1(check1);
                                                AVbeta1(check2)=AVbeta1(check2);
                                                Nprob=AVbeta1/sum(AVbeta1);
                                                prob=rand;
                                                Hcheck=Hx;
                                                Hx=Hx-1*check1(1)*check2(1)*(prob<Nprob(1))-1*check1(2)*check2(2)*(prob<Nprob(2))-1*check1(3)*check2(3)*(prob<Nprob(3))-1*check1(4)*check2(4)*(prob<Nprob(4));
                        %                         if Hx-Hcheck==0
                        %                             FLAG=0;
                        %                             break
                        %                         end
                                                N(1)=N(1)+check1(1)*check2(1)*(prob<Nprob(1));
                                                N(2)=N(2)+check1(2)*check2(2)*(prob<Nprob(2));
                                                N(3)=N(3)+check1(3)*check2(3)*(prob<Nprob(3));
                                                N(4)=N(4)+check1(4)*check2(4)*(prob<Nprob(4));
                                                beta1=newbeta1(N,Hx);
                                                beta2=newbeta2(N,beta1);
                                                check1=newcheck1(N,beta1);
                                                check2=newcheck2(N,beta2);
                                                %flag=1;
                                                if Hx-Hcheck~=0
                                                    CellCt=CellCt+1;
                                                    flag=1;
                                                end
                                                if sum(check1)==0 
                                                    Hstar(R,C)=Hx;
                                                    Hstar(Rm,C)=N(1);
                                                    Hstar(R,Cp)=N(2);
                                                    Hstar(Rp,C)=N(3);
                                                    Hstar(R,Cm)=N(4);
                                                    %fprintf('Changes have been made')
                                                    FLAG=0;
                                                    if Hx >=0
                                                        break
                                                    end
                                                end
                                                if check1*check2'==0
                                                    break
                                                end
                                            end
                                        end
                                        if sum(check1)==0
                                            FLAG=0;
                                            %no avalanching needed
                                        end
                                    elseif sum(check2)==0 && sum(check1)~=0
                                       FLAG=0;
                                       %fprintf('Too many plants!')
                                    elseif sum(check1)==0
                                        FLAG=0;     %only other case should be if both sum to 0
                                    end

                                else%{if Hstar(R,C)<=0}
                                    Erosioncheck=0;
                                    FLAG=0;
                                end
                                FLAG=0;
                            end
                        end
                    end
            end

        end
    %end
%end

%this portion is for elevation maps with multiple full sized islands within the domain
%we rarely need this, but it's here just in case

% for R=1:size(Hstar,1)
%     %ISLANDindex=(Hstar(R,:,1)>0);
%     if sum(ISLANDindex)~=0      %skip whole rows of water
%         for C=1:size(Hstar,2)
%             Erosioncheck=1;
%             FLAG=1;
%             while FLAG==1
%                 RAND=rand;
%                 %FLAG=1;
%                 if Hstar(R,C)>0
%                     if R<=(size(Hstar,1)-1)
%                         R=floor(R);
%                         Rp=R+1;
%                     else 
%                         Rp=floor(R);
%                     end
%                     if R>=2
%                         R=floor(R);
%                         Rm=R-1;
%                     else
%                         Rm=floor(R);
%                     end
%                     if C<=(size(Hstar,2)-1)
%                         C=floor(C);
%                         Cp=C+1;
%                     else
%                         Cp=floor(C);
%                     end
%                     if C>=2
%                         C=floor(C);
%                         Cm=C-1;
%                     else
%                         Cm=floor(C);
%                     end
% 
%                     Hx=Hstar(R,C);                                              %current cell value
%                     N=[Hstar(Rm,C),Hstar(R,Cp),Hstar(Rp,C),Hstar(R,Cm)];        %neighborhood
% 
%                     beta1=atan(((Hx-N)*delta)/L);
%                     check1=(beta1>=(AoR));
%                     %Avdir=find(beta1>=AoR);
%                     beta2=alpha*((beta1/(AoR))*(1-PC(R,C)));
%                     check2=(RAND<beta2);
% 
% 
%     %                 if sum(check1)~=0
%     %                     if sum(check2)~=0
%     %                     end
%     %                 end
%                         if R==24
%     %                         fprintf('STOP!')
%                         end
% 
%                     if sum(check2)~=0
%                         if check1*check2'~=0
%                             while  sum(check1)~=0 && sum(check2)~=0  %if I use a while loop, there's no point to using probability
%                                 %posbeta1=beta1;
%                                 %posbeta1(find(beta1<0))=0;
%                                 %Nprob=posbeta1/sum(posbeta1);
%                                 AVbeta1=zeros(1,4);
%                                 AVbeta1(check1)=beta1(check1);
%                                 AVbeta1(check2)=AVbeta1(check2);
%                                 Nprob=AVbeta1/sum(AVbeta1);
%                                 prob=rand;
%                                 Hcheck=Hx;
%                                 Hx=Hx-1*check1(1)*check2(1)*(prob<Nprob(1))-1*check1(2)*check2(2)*(prob<Nprob(2))-1*check1(3)*check2(3)*(prob<Nprob(3))-1*check1(4)*check2(4)*(prob<Nprob(4));
%         %                         if Hx-Hcheck==0
%         %                             FLAG=0;
%         %                             break
%         %                         end
%                                 N(1)=N(1)+check1(1)*check2(1)*(prob<Nprob(1));
%                                 N(2)=N(2)+check1(2)*check2(2)*(prob<Nprob(2));
%                                 N(3)=N(3)+check1(3)*check2(3)*(prob<Nprob(3));
%                                 N(4)=N(4)+check1(4)*check2(4)*(prob<Nprob(4));
%                                 beta1=newbeta1(N,Hx);
%                                 beta2=newbeta2(N,beta1);
%                                 check1=newcheck1(N,beta1);
%                                 check2=newcheck2(N,beta2);
%                                 %flag=1;
%                                 if Hx-Hcheck~=0
%                                     flag=1;
%                                 end
%                                 if sum(check1)==0 
%                                     Hstar(R,C)=Hx;
%                                     Hstar(Rm,C)=N(1);
%                                     Hstar(R,Cp)=N(2);
%                                     Hstar(Rp,C)=N(3);
%                                     Hstar(R,Cm)=N(4);
%                                     %fprintf('Changes have been made')
%                                     FLAG=0;
%                                     if Hx >=0
%                                         break
%                                     end
%                                 end
%                                 if check1*check2'==0
%                                     break
%                                 end
%                             end
%                         end
%                         if sum(check1)==0
%                             FLAG=0;
%                             %no avalanching needed
%                         end
%                     elseif sum(check2)==0 && sum(check1)~=0
%                        FLAG=0;
%                        %fprintf('Too many plants!')
%                     elseif sum(check1)==0
%                         FLAG=0;     %only other case should be if both sum to 0
%                     end
% 
%                 else%{if Hstar(R,C)<=0}
%                     Erosioncheck=0;
%                     FLAG=0;
%                 end
%                 FLAG=0;
%             end
%         end
%     end
%         
% end    
% end


%Local functions for finding angle of repose and probability WRT plant cover
function B1=newbeta1(N,Hx)
    %B1=zeros(1,length(N));
    B1=atan(((Hx-N)*delta)/L);

    
%     for k=1:length(N)
%     B1(k)=atan(((Hx-N(k))*delta)/L);
%     end
end

function B2=newbeta2(N,beta1)
    %B2=zeros(1:length(N));
    B2=alpha*((beta1/(AoR))*(1-PC(R,C)));

%     for k=1:length(N)
%     B2(k)=alpha*((beta1(k)/(AoR))*(1-PC(R,C)));
%     end
end

function C1=newcheck1(N,beta1)
    %C1=zeros(1,length(N));
    C1=(beta1>=(AoR));
%     for k=1:length(N)
%     C1(k)=(beta1(k)>=(AoR));
%     end
end

function C2=newcheck2(N,beta2)
    %C2=zeros(1,length(N));
    C2=(RAND<beta2);
%     for k=1:length(N)
%     C2(k)=(RAND<beta2(k));
%     end
end

%FLAGCHECK=Hstarflag-Hstar;
%if sum(sum(Hstarflag-Hstar))~=0
    %flag=1;
%end


%  toc               
end 