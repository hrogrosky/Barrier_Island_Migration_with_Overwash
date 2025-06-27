function [Hstar] = Smoothing_water(Hstar)

%%%%%%% COMMENT THESE OUT WHEN NOT TESTING CODE %%%%%%%%
%  filename='Parramore03312021.mat'; 
%  H = importdata(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=Hstar(:,:,1);
H_smooth = Hstar(:,:,1);
[n1, n2]=size(H(:,:));
EastSL=zeros(n1,1);
WestSL=zeros(n1,1);

%%%% To determine the eastern shore line (ESL) (from Overwash code)
for i=1:n1
    for j=n2:-1:1
        if H(i,j)>0
            EastSL(i)=j; 
            break
       end
    end
end

for i=1:n1
    for j=1:n2
        if H(i,j)>0
            WestSL(i)=j; 
            break
       end
    end
end

%%%%% creating 3x3 square and averaging elevation in square

%%%%%% FIRST ATTEMPT %%%%%%
% for i = 1:n1
%     if EastSL(i)>0
%         Eshore=EastSL(i); 
%         Wshore=WestSL(i);
%         if (0==mod(i,3))
%         for j= Eshore:-1:Wshore 
%             l=1;
%             j1 = j - (l);
%             j2 = j + (l);
%             i1 = i - (l);
%             i2 = i + (l);
%             avg1 = ceil(sum(sum(H(i1:i2,j1:j2)))/9);
%             
%             H(i1:i2,j1:j2)= avg1;
%             if j1 <= WestSL(i)
%                 break
%             end
%         end
%         end
%     end
% end
% 
% Hstar = H;
%%%%%%%%%%%%%%%

%%%%% SECOND ATTEMPT %%%%%
for i = 1:n1
    if EastSL(i)>0
        Eshore=EastSL(i); 
        Wshore=WestSL(i);
        %if (0==mod(i,3))
        for j= (Eshore+1):-1:(Eshore+40) 
            l=1;
            j1 = max(j - (l),1);
            j2 = j + (l);
            i1 = i - (l);
            i2 = i + (l);
            avg1 = round(sum(sum(H(i1:i2,j1:j2)))/9); %using original matrix
            
            H_smooth(i,j)= avg1; %only change center cell
            if j1 <= EastSL(i)+40
                break
            end
        end
         for j= (Wshore-1):-1:(Wshore-20) 
            l=1;
            j1 = max(j - (l),1);
            j2 = j + (l);
            i1 = i - (l);
            i2 = i + (l);
            avg1 = round(sum(sum(H(i1:i2,j1:j2)))/9); %using original matrix
            
            H_smooth(i,j)= avg1; %only change center cell
            if j1 <= WestSL(i)-20
                break
            end
        end
        %end
    end
end

Hstar = H_smooth;
%%%%%%%%%%%%

%%%% after figure commented out for main code
% imagesc(H_smooth)
% colorbar
% clim([-5 26])
end