function [ESL_hist] = storing_ESL(Hstar,t, time, ESL_hist)

%%%%%%% COMMENT THESE OUT WHEN NOT TESTING CODE %%%%%%%%
%  filename='Parramore03312021.mat'; 
%  H = importdata(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=Hstar(:,:,1);
[n1, n2]=size(H(:,:));
EastSL=zeros(n1,1);

%%%% To determine the eastern shore line (ESL) (from Overwash code)
for i=1:n1
    for j=n2:-1:1
        if H(i,j)>0
            EastSL(i)=j; 
            break
       end
    end
end

% ESL_hist = zeros(n1,time+1);
ESL_hist(:,t+1) = EastSL;
% ESL_diff = zeros(n1,3);
% ESL_diff(:,1) = abs(ESL_hist(:,365) - ESL_hist(:,1));
% ESL_diff(:,2) = abs(ESL_hist(:,703) - ESL_hist(:,365));
% ESL_diff(:,3) = abs(ESL_hist(:,703) - ESL_hist(:,1));


% figure;
% plot(ESL_hist(:,1),'b.');
% hold on 
% plot(ESL_hist(:,365),'r.');
% hold on 
% plot(ESL_hist(:,703),'k.');
% legend({'Initial', '14 years', '27 years'},'FontSize',14);
% xlabel('Row location (North to South)', 'FontSize',14);
% ylabel('Shore line column location', 'FontSize',14);
% title('Eastern shore line movement over time','FontSize',14)

% figure;
% plot(ESL_diff(:,1),'b');
% hold on 
% plot(ESL_diff(:,2),'r');
% hold on
% plot(ESL_diff(:,3),'k');
% xlim([500 2200]);
% [~, objh] = legend({'Difference between initial and year 14', 'Difference between year 14 and year 27',...
%     'Difference between initial and year 27'},'FontSize',14);
% objhl = findobj(objh, 'type', 'line');
% set(objhl, 'Linewidth', 3)
% xlabel('Row location (North to South)', 'FontSize',14);
% ylabel('Cell difference', 'FontSize',14);
% title('Eastern shore line difference after overwash','FontSize',14)