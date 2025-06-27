function [Hstar,flag,CellCt] = AVALANCHE03312021(Hstar,delta,L,flag,PC)
%%
%This version of avalanche:
%       ~weights the probability of avalanching to favor the direction of 
%           the greatest violation of the angle of repose
%       ~ Does NOT take time as an input and is only used in initialization
%           process at beginning of the MainCode - runs over whole domain
% tic
Htest=0;
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

for R=1:size(Hstar,1)
    for C=1:size(Hstar,2)
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
                check1=(beta1>=(pi/6));
                %Avdir=find(beta1>=pi/6);
                beta2=alpha*((beta1/(pi/6))*(1-PC(R,C)));
                check2=(RAND<beta2);
                
                
%                 if sum(check1)~=0
%                     if sum(check2)~=0
%                     end
%                 end

                       
                if sum(check2)~=0
                    %if check1*check2'~=0
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
                            end
%                             if Hx <=0
%                                 break
%                             end
                            if check1*check2'==0
                                break
                            end
                        end
                    %end
                
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
    B2=alpha*((beta1/(pi/6))*(1-PC(R,C)));

%     for k=1:length(N)
%     B2(k)=alpha*((beta1(k)/(pi/6))*(1-PC(R,C)));
%     end
end

function C1=newcheck1(N,beta1)
    %C1=zeros(1,length(N));
    C1=(beta1>=(pi/6));
%     for k=1:length(N)
%     C1(k)=(beta1(k)>=(pi/6));
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

if CellCt<5
    flag=0;
end
%  toc               
end 