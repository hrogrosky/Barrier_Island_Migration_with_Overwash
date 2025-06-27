% Figure plotting script for previously saved data

%%Load data for elevation plots - change filenamestring (Line 4) and time of data (yearnum, Line 5) as needed
% filenamestring='St00te702BP2000SLRi0p0044872SLRa4p3077en05MMax1'; % Enter folder name here!!!!
% yearnum=0; 
% islandstring=filenamestring(1);
% cd /Users/beththomas/Documents/'Documents - Beth MacBook Air'/MATLAB/'Barrier Island'
% curdir=pwd;
% cd(filenamestring);
% filenamestringmat=strcat([filenamestring,'t',int2str(yearnum),'.mat']);
% load(filenamestringmat);
% cd ..

%% Elevation - (if ElevIm_during==1)

% figfontsize=15; % sometimes matlab chooses pretty small fonts...
% 
% MElevCont=0;
% SElevCont=0;
% 
% 
% Hfilt=round(imgaussfilt(H(:,:,1),4),0);
% climsH=[min(min(H(:,:,1))) max(max(H(:,:,1)))];
% climsP=[-1 1];
% figure   
% %
% colormap(jet())
% imagesc(H(:,:,1),climsH);
% grid on
% grid minor
% hold on
% if MElevCont==0
%    contour(H(:,:,1),[-0 -0],'color','m','linewidth',1.1) %try changing to [0 0]
%     hold on
%    %contour(H_initial(:,:,1),[-0 -0],'color','m','linewidth',1.1)
% end
% if SElevCont==1
%     contour(Hfilt(:,:,1),[-0 -0],'color','k','linewidth',1.1)
% end
% if MElevCont==1 && SElevCont==1
%     lgnd=legend('\color{white} marsh contour (-0.5m)','\color{white} sea level contour (0m)');
%     set(lgnd,'color','none','location','southeast');
% elseif MElevCont==1 && SElevCont==0
%     lgnd=legend('\color{white} marsh contour (-0.5m)');
%     set(lgnd,'color','none','location','southeast');
% elseif SElevCont==1 && MElevCont==0
%     lgnd=legend('\color{white} sea level contour (0m)');
%     set(lgnd,'color','none','location','southeast');
% end
% colorbar
% title(sprintf('Elevation at %d years',t/26))
% hold off
% %%%%
% figureHandle = gcf;
% set(gca,'fontsize',figfontsize);
% set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
% figfilenamestring=strcat(['FigElt',num2str(t/26),'_',filenamestring]);
% cd(filenamestring);
% saveas(gcf,figfilenamestring); % save as .fig
% saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
% cd ..
% %%%

%%%%%%%%%%
% % %% Load data for Elevation contours
filenamestring='St00te1300BP2000SLRi0p0044872SLRa4p3077en05MMax1'; % Enter folder name here!!!!
yearnumvec=[0 10 20]; 
islandstring=filenamestring(1);
cd /Users/beththomas/Documents/'Documents - Beth MacBook Air'/MATLAB/'Barrier Island'
curdir=pwd;

cd(filenamestring);
filenamestringmat=strcat([filenamestring,'t',int2str(yearnumvec(1)),'.mat']);
load(filenamestringmat);
cd ..
[Hstar,Centroid]=Qdata03312021(Hstar,P1,P2,P3,P4);%,ESL,WSL,MESL,MWSL,OppositeLocation);
Sv1f=round(imgaussfilt(H(:,:,1),4),0);
Sv1Cent=Centroid;

cd(filenamestring);
filenamestringmat=strcat([filenamestring,'t',int2str(yearnumvec(2)),'.mat']);
load(filenamestringmat);
cd ..
[Hstar,Centroid]=Qdata03312021(Hstar,P1,P2,P3,P4);%,ESL,WSL,MESL,MWSL,OppositeLocation);
Sv2f=round(imgaussfilt(H(:,:,1),4),0);
Sv2Cent=Centroid;

cd(filenamestring);
filenamestringmat=strcat([filenamestring,'t',int2str(yearnumvec(3)),'.mat']);
load(filenamestringmat);
cd ..
[Hstar,Centroid]=Qdata03312021(Hstar,P1,P2,P3,P4);%,ESL,WSL,MESL,MWSL,OppositeLocation);
Sv3f=round(imgaussfilt(H(:,:,1),4),0);
Sv3Cent=Centroid;

%% Contours of the shoreline (0m elevation)    

figfontsize=15; % sometimes matlab chooses pretty small fonts...

figure;                                                                                                                           %
colormap(jet());
hold on
contour(Sv1f(:,:,1),[1 1],'Color','blue','LineWidth',2);           %beginning contour (after initialization)
contour(Sv2f(:,:,1),[1 1],'Color','green','LineWidth',2);          %contour at 13 years (t=338)
contour(Sv3f(:,:,1),[1 1],'Color','red','LineWidth',2);            %contour at 27 years (t=702)
scatter(Sv1Cent(1),Sv1Cent(2),'o','MarkerEdgeColor','k','MarkerFaceColor','b')
scatter(Sv2Cent(1),Sv2Cent(2),'o','MarkerEdgeColor','k','MarkerFaceColor','g')
scatter(Sv3Cent(1),Sv3Cent(2),'o','MarkerEdgeColor','k','MarkerFaceColor','r')
% scatter(Sv1Cent(1),Sv1Cent(2),'o','MarkerEdgeColor','w',...
%                                    'MarkerFaceColor','b',...
%                                    'LineWidth',1.25)
% scatter(Sv2Cent(1),Sv2Cent(2),'o','MarkerEdgeColor','w',...
%                                    'MarkerFaceColor','g',...
%                                    'LineWidth',1.25)                               
% scatter(Sv3Cent(1),Sv3Cent(2),'o','MarkerEdgeColor','w',...
%                                    'MarkerFaceColor','r',...
%                                    'LineWidth',1.25)
set(gca, 'YDir','reverse')
legend(strcat([int2str(yearnumvec(1)),' years']),strcat([int2str(yearnumvec(2)),' years']),strcat([int2str(yearnumvec(3)),' years']),'Location','southeast')
title('Island Contours at 0m Elevation');
%%%%
set(gcf,'Position',[100,100,380,400]);
figureHandle = gcf;
set(gca,'fontsize',figfontsize);
set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
figfilenamestring=strcat(['FigCont',int2str(yearnumvec(1)),int2str(yearnumvec(2)),int2str(yearnumvec(3)),'_',filenamestring]);
cd(filenamestring);
saveas(gcf,figfilenamestring); % save as .fig
saveas(gcf,strcat([figfilenamestring,'.jpg'])); % save as .jpg
cd ..
%%


%mesh(H)