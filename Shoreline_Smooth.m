%% Smoothing the Eastern and Western Shorelines %%
function [Hstar,EastSL_change] = Shoreline_Smooth(Hstar)
%% definining neccessary vectors
H=Hstar(:,:,1);
[n1, n2]=size(H(:,:));
EastSL=zeros(n1,1);
EastSL_smooth = zeros(n1,1);
EastSL_1=zeros(n1,1);
EastSL_1_smooth = zeros(n1,1);
EastSL_2=zeros(n1,1);
EastSL_2_smooth = zeros(n1,1);
EastSL_3=zeros(n1,1);
EastSL_3_smooth = zeros(n1,1);
EastSL_4=zeros(n1,1);
EastSL_4_smooth = zeros(n1,1);
EastSL_change = zeros(n1,1);


% SIZE OF THE WINDOW AROUND EACH POINT
N=150; 

%% To determine the eastern shore line (ESL) (from Overwash code)
for i=1:n1
    for j=n2:-1:1
        if H(i,j)>0
            EastSL(i)=j; 
            break
       end
    end
end

%%%% Smoothing the EastSL out and calculate the change in shoreline

row_top = find((EastSL>0),1,'first');
row_bot = find((EastSL>0),1,'last');

EastSL(row_top-50:row_top-1) = EastSL(row_top);
EastSL(row_bot+50:row_top+1) = EastSL(row_bot);
row_top1 = find((EastSL>0),1,'first');
row_bot1 = find((EastSL>0),1,'last');

% use medfilt1 to smooth the ESL and WSL a 150th order one-dimensional median fit
EastSL_smooth(row_top1:row_bot1) = round(medfilt1(EastSL(row_top1:row_bot1),N));
EastSL_smooth(row_top1:row_top-1) = 0;
EastSL_smooth(row_bot+1:row_bot1) = 0;

for i=row_top:row_bot
if EastSL_smooth(i)~=0 && EastSL(i)~= 0 
if EastSL_smooth(i) > EastSL(i)
    H(i, EastSL(i):EastSL_smooth(i)) = 1; %fills in the island between old and new shoreline
end
if EastSL_smooth(i) < EastSL(i)
    H(i, EastSL_smooth(i):EastSL(i)) = 0;
end
end
end

EastSL_change = EastSL_smooth - EastSL;


%% Smooth Elevation 1
for i=1:n1
    for j=n2:-1:1
        if H(i,j)>1
            EastSL_1(i)=j; 
            break
       end
    end
end

%%%% Smoothing the EastSL out and calculate the change in shoreline

row_top = find((EastSL_1>0),1,'first');
row_bot = find((EastSL_1>0),1,'last');

EastSL_1(row_top-50:row_top-1) = EastSL_1(row_top);
EastSL_1(row_bot+50:row_top+1) = EastSL_1(row_bot);
row_top1 = find((EastSL_1>0),1,'first');
row_bot1 = find((EastSL_1>0),1,'last');

% use medfilt1 to smooth the ESL and WSL a 150th order one-dimensional median fit
EastSL_1_smooth(row_top1:row_bot1) = round(medfilt1(EastSL_1(row_top1:row_bot1),N));
EastSL_1_smooth(row_top1:row_top-1) = 0;
EastSL_1_smooth(row_bot+1:row_bot1) = 0;

for i=row_top:row_bot
if EastSL_1_smooth(i)~=0 && EastSL_1(i)~= 0 
if EastSL_1_smooth(i) > EastSL_1(i)
    H(i, EastSL_1(i):EastSL_1_smooth(i)) = 2; %fills in the island between old and new shoreline
end
if EastSL_1_smooth(i) < EastSL_1(i)
    H(i, EastSL_1_smooth(i):EastSL_1(i)) = 1;
end
end
end

%% Smooth Elevation 2
for i=1:n1
    for j=n2:-1:1
        if H(i,j)>2
            EastSL_2(i)=j; 
            break
       end
    end
end

%%%% Smoothing the EastSL out and calculate the change in shoreline

row_top = find((EastSL_2>0),1,'first');
row_bot = find((EastSL_2>0),1,'last');

EastSL_2(row_top-50:row_top-1) = EastSL_2(row_top);
EastSL_2(row_bot+50:row_top+1) = EastSL_2(row_bot);
row_top1 = find((EastSL_2>0),1,'first');
row_bot1 = find((EastSL_2>0),1,'last');

% use medfilt1 to smooth the ESL and WSL a 150th order one-dimensional median fit
EastSL_2_smooth(row_top1:row_bot1) = round(medfilt1(EastSL_2(row_top1:row_bot1),N));
EastSL_2_smooth(row_top1:row_top-1) = 0;
EastSL_2_smooth(row_bot+1:row_bot1) = 0;

for i=row_top:row_bot
if EastSL_2_smooth(i)~=0 && EastSL_2(i)~= 0 
if EastSL_2_smooth(i) > EastSL_2(i)
    H(i, EastSL_2(i):EastSL_2_smooth(i)) = 3; %fills in the island between old and new shoreline
end
if EastSL_2_smooth(i) < EastSL_2(i)
    H(i, EastSL_2_smooth(i):EastSL_2(i)) = 2;
end
end
end

%% Smooth Elevation 3
for i=1:n1
    for j=n2:-1:1
        if H(i,j)>3
            EastSL_3(i)=j; 
            break
       end
    end
end

%%%% Smoothing the EastSL out and calculate the change in shoreline

row_top = find((EastSL_3>0),1,'first');
row_bot = find((EastSL_3>0),1,'last');

EastSL_3(row_top-50:row_top-1) = EastSL_3(row_top);
EastSL_3(row_bot+50:row_top+1) = EastSL_3(row_bot);
row_top1 = find((EastSL_3>0),1,'first');
row_bot1 = find((EastSL_3>0),1,'last');

% use medfilt1 to smooth the ESL and WSL a 150th order one-dimensional median fit
EastSL_3_smooth(row_top1:row_bot1) = round(medfilt1(EastSL_3(row_top1:row_bot1),N));
EastSL_3_smooth(row_top1:row_top-1) = 0;
EastSL_3_smooth(row_bot+1:row_bot1) = 0;

for i=row_top:row_bot
if EastSL_3_smooth(i)~=0 && EastSL_3(i)~= 0 
if EastSL_3_smooth(i) > EastSL_3(i)
    H(i, EastSL_3(i):EastSL_3_smooth(i)) = 4; %fills in the island between old and new shoreline
end
if EastSL_3_smooth(i) < EastSL_3(i)
    H(i, EastSL_3_smooth(i):EastSL_3(i)) = 3;
end
end
end

%% Smooth Elevation 4
for i=1:n1
    for j=n2:-1:1
        if H(i,j)>4
            EastSL_4(i)=j; 
            break
       end
    end
end

%%%% Smoothing the EastSL out and calculate the change in shoreline

row_top = find((EastSL_4>0),1,'first');
row_bot = find((EastSL_4>0),1,'last');

EastSL_4(row_top-50:row_top-1) = EastSL_4(row_top);
EastSL_4(row_bot+50:row_top+1) = EastSL_4(row_bot);
row_top1 = find((EastSL_4>0),1,'first');
row_bot1 = find((EastSL_4>0),1,'last');

% use medfilt1 to smooth the ESL and WSL a 150th order one-dimensional median fit
EastSL_4_smooth(row_top1:row_bot1) = round(medfilt1(EastSL_4(row_top1:row_bot1),N));
EastSL_4_smooth(row_top1:row_top-1) = 0;
EastSL_4_smooth(row_bot+1:row_bot1) = 0;

for i=row_top:row_bot
if EastSL_4_smooth(i)~=0 && EastSL_4(i)~= 0 
if EastSL_4_smooth(i) > EastSL_4(i)
    H(i, EastSL_4(i):EastSL_4_smooth(i)) = 5; %fills in the island between old and new shoreline
end
if EastSL_4_smooth(i) < EastSL_4(i)
    H(i, EastSL_4_smooth(i):EastSL_4(i)) = 4;
end
end
end
%% Plot outline and smoothed outline
% figure(1)
% plot(EastSL,Color='blue')
% hold on
% plot(EastSL_smooth,Color='red')
% legend('ESL before smoothing','ESL after smoothing',fontsize=18)
% title('Smoothing ESL yr 27 with medfilt1',FontSize=22)
% 
% figure(2)
% plot(EastSL_1,Color='blue')
% hold on
% plot(EastSL_1_smooth,Color='red')
% legend('Elev 1 before','Elev 1 after',fontsize=18)
% title('Smoothing Elev 1 yr 27 with medfilt1',FontSize=22)
% 
% figure(3)
% plot(EastSL_2,Color='blue')
% hold on
% plot(EastSL_2_smooth,Color='red')
% legend('Elev 2 before','Elev 2 after',fontsize=18)
% title('Smoothing Elev 2 yr 27 with medfilt1',FontSize=22)
% 
% figure(4)
% plot(EastSL_3,Color='blue')
% hold on
% plot(EastSL_3_smooth,Color='red')
% legend('Elev 3 before','Elev 3 after',fontsize=18)
% title('Smoothing Elev 3 yr 27 with medfilt1',FontSize=22)
% 
% figure(5)
% plot(EastSL_4,Color='blue')
% hold on
% plot(EastSL_4_smooth,Color='red')
% legend('Elev 4 before','Elev 4 after',fontsize=18)
% title('Smoothing Elev 4 yr 27 with medfilt1',FontSize=22)

Hstar(:,:,1) = H;
end
