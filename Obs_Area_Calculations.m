% CALCULATING THE AREA OF PARRAMORE IN 1984, 1998, AND 2011
% STEP ONE: upload each individual image into MATLAB
% STEP TWO: fill in island with 1's and all other 0's 
% STEP THREE: remove borders if neccessary
% STEP FOUR: adjust matrix so all of the islands are the same size 

%% Parramore 1984 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=imread('Parramore 1984.png');
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
H_1984(4:8,:) = 0;
H_1984(3299:end,:) = 0;

%% Parramore 1998 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=imread('Parramore 1998.png');
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
H_1998(2:6,:) = 0;
H_1998(:,2:7) = 0;
H_1998(:,2547:end) = 0;
H_1998(3294:3299,:) = 0;

%shift 82 columns to the right and down 60
H_1998 = imresize(H_1998,0.92);
H_1998(H_1998 ~= 0) = 1;
H_1998 = circshift(H_1998,[60 55]);


%% Parramore 2011 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=imread('Parramore 2011.png');
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
H_2011(1:2,:) = 0;
H_2011(:,1:2) = 0;
H_2011(:,2547:end) = 0;
H_2011(3298:end,:) = 0;

%shifting 155 units up and 175 units to the left
H_2011 = imresize(H_2011,1.08);
H_2011(H_2011 ~= 0) = 1;
H_2011 = circshift(H_2011,[-155 -175]);


%% plot all 3 at the same time %%%%%%%%%%%%%%%%%
% contour(H_1984,[1 1],'Color','blue','LineWidth',2); 
% hold on
% contour(H_1998,[1 1],'Color','green','LineWidth',2); 
% hold on
% contour(H_2011,[1 1],'Color','red','LineWidth',2);
% set(gca, 'YDir','reverse')

%% Calculate Area
Area = [sum(sum(H_1984>0)), sum(sum(H_1998>0)), sum(sum(H_2011>0))]
