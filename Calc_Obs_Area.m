% CALCULATING THE AREA OF SMITH IN 1984, 1998, AND 2011
% STEP ONE: upload each individual image into MATLAB
% STEP TWO: fill in island with 1's and all other 0's 
% STEP THREE: remove borders if neccessary
% STEP FOUR: adjust matrix so all of the islands are the same size 

%% Smith 1984 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=imread('Smith 1984 1.png');
[m1,m2]=size(I(:,:,1));
H_1984 = zeros(m1,m2);

% this loop isolates cells whose color is close to 0 (black), making
% those entries 1 and everything else 0
for i = 1:m1
    for j = 1:m2
        if I(i,j,1) <= 10
            H_1984(i,j) = 1;
        end
    end
end

%removing the border from the image (was black but we don't want it)
H_1984(3:5,:) = 0;
H_1984(:,1:2) = 0;
H_1984(3296:3298,:) = 0;

%% Smith 1998 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=imread('Smith 1998 1.png');
[m1,m2]=size(I(:,:,1));
H_1998 = zeros(m1,m2);

% this loop isolates cells whose color is close to 0 (black), making
% those entries 1 and everything else 0
for i = 1:m1
    for j = 1:m2
        if I(i,j,1) <= 10
            H_1998(i,j) = 1;
        end
    end
end

%removing the border from the image (was black but we don't want it)
H_1998(3:5,:) = 0;
H_1998(:,1:2) = 0;
H_1998(3296:3298,:) = 0;

% %shift 30 units to the left
H_1998 = circshift(H_1998,[0 -59]);


%% Smith 2011 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=imread('Smith 2011 1.png');
[m1,m2]=size(I(:,:,1));
H_2011 = zeros(m1,m2);

% this loop isolates cells whose color is close to 0 (black), making
% those entries 1 and everything else 0
for i = 1:m1
    for j = 1:m2
        if I(i,j,1) <= 10
            H_2011(i,j) = 1;
        end
    end
end

%removing the border from the image (was black but we don't want it)
H_2011(3:5,:) = 0;
H_2011(:,1:2) = 0;
H_2011(3296:3298,:) = 0;

% %shifting 51 units to the left
H_2011 = circshift(H_2011,[0 -51]);


% % plot all 3 at the same time %%%%%%%%%%%%%%%%%
% contour(H_1984,[1 1],'Color','blue','LineWidth',2); 
% hold on
% contour(H_1998,[1 1],'Color','green','LineWidth',2); 
% hold on
% contour(H_2011,[1 1],'Color','red','LineWidth',2);
% set(gca, 'YDir','reverse')

%% Calculate Area
Area = [sum(sum(H_1984>0)), sum(sum(H_1998>0)), sum(sum(H_2011>0))]