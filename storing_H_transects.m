function [top_hist1, mid_hist1, bot_hist1, top_hist2, mid_hist2, bot_hist2, H_cell_count] = storing_H_transects(Hstar,t,top_hist1, mid_hist1, bot_hist1, top_hist2, mid_hist2, bot_hist2,H_cell_count,top_transect1, mid_transect1, bot_transect1, top_transect2, mid_transect2, bot_transect2)

% filename='Parramore03312021.mat'; 
%  H = importdata(filename);

H=Hstar(:,:,1);

%stores the transects
top_hist1(t+1,:) = H(top_transect1,:);
mid_hist1(t+1,:) = H(mid_transect1,:);
bot_hist1(t+1,:) = H(bot_transect1,:);

top_hist2(t+1,:) = H(top_transect2,:);
mid_hist2(t+1,:) = H(mid_transect2,:);
bot_hist2(t+1,:) = H(bot_transect2,:);

%calculating cells above 20, 10-20, and 0-10
H_high = length(H(H>=20));
H_med = length(H(H>= 10)) - H_high;
H_low = length(H(H>=0)) - H_med;

H_cell_count(1,t+1) = H_high;
H_cell_count(2,t+1) = H_med;
H_cell_count(3,t+1) = H_low;
end


%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting the transects %%%%%%%%%%%%%%%%%%%%
% Pwin = [1600:2377];
% Pwin2 = [500:1277];
% Swin = [1250:2250];
% Swin2 = [300:1300];
% 
% t = tiledlayout(2,1,'TileSpacing','Compact');
% title(t,'Elevation Over Time on Smith Island','FontSize',18)
% 
% nexttile
% plot(top_hist1(1,Swin),'b',Linewidth=2)
% hold on
% plot(top_hist1(365,Swin),'g',Linewidth=2)
% hold on
% plot(top_hist1(703,Swin),'r',Linewidth=2)
% legend({'Initial', '14 yr', '27 yr'},'FontSize',14);
% xlabel('Column location (West to East)', 'FontSize',14);
% ylabel('Elevation (0.1 m)', 'FontSize',14);
% title('Top Transect','FontSize',18)
% 
% nexttile
% plot(bot_hist1(1,Swin2),'b',Linewidth=2)
% hold on
% plot(bot_hist1(365,Swin2),'g',Linewidth=2)
% hold on
% plot(bot_hist1(703,Swin2),'r',Linewidth=2)
% % legend({'Initial', '14 yr', '27 yr'},'FontSize',14);
% xlabel('Column location (West to East)', 'FontSize',14);
% ylabel('Elevation (0.1 m)', 'FontSize',14);
% title('Bottom Transect','FontSize',18)






%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting elevation categories over time %%%%%
% plot(log(H_cell_count(1,:)),'b.')
% hold on
% plot(log(H_cell_count(2,:)),'r.')
% hold on
% plot(log(H_cell_count(3,:)),'g.')
% legend({'Cells above 20', 'Cells between 10 and 20', 'Cells between 0 and 10'},'FontSize',14);
% xlabel('Time in two weeks', 'FontSize',14);
% xlim([0 1301]);
% ylim([5 15.5]);
% ylabel('log(number of cells)', 'FontSize',14);
% title('Number of Cells in Different Elevation Classes','FontSize',14)

%%%%%%%%%%%%%% Plotting error_vec %%%%%%%%%%%%%%
% figure;
% subplot(3,1,1)
% plot(log(error_vec(:,2)),'b')
% hold on
% xlabel('time (in 2 weeks)', 'FontSize',14);
% ylabel('Cells', 'FontSize',14);
% title('Number of cells < 0','FontSize',14)
% 
% subplot(3,1,2)
% plot(log(error_vec(:,3)),'b')
% hold on
% xlabel('time (in 2 weeks)', 'FontSize',14);
% ylabel('Cells', 'FontSize',14);
% title('Number of cells = 0','FontSize',14)
% 
% subplot(3,1,3)
% plot(log(error_vec(:,4)),'b')
% xlabel('time (in 2 weeks)', 'FontSize',14);
% ylabel('Cells', 'FontSize',14);
% title('Number of cells > 0','FontSize',14)

