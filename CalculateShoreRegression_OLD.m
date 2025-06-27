%% Loads an image with three outlines, one for year 0, 14, and 27
I=imread('ParramoreOutlineReal_nolegend2.jpg');
[nr,nc]=size(I);
ParResize=imresize(I,[nr*3,nc]); %*6 for Smith
[nr,nc]=size(ParResize);
nrc=200; ncc=200;
ParCrop=imcrop(ParResize,[nrc,ncc+190,nr-nrc*8.8,nc-ncc*26.7]); %13.59
imshow(ParCrop);


figure;
contourf(flipud(ParCrop(:,:,2)),'Linestyle','none');


figure;
contourf(flipud(ParCrop(:,:,3)),'Linestyle','none');

%% 
[n1,n2,~]=size(ParCrop);
b1=4;
b2=10;
b3=20;
winr=[1 140];
wing=[1 160];
winb=[1 100];
min1=zeros(n1,1);
min2=zeros(n1,1);
min3=zeros(n1,1);
for ii=b1:n1-b1
    for jj=n2-b2:-1:1+b2
        if ParCrop(ii,jj,1)<=min([ParCrop(ii,jj-b2:jj-1,1) ParCrop(ii,jj+1:jj+b2,1)])
            if min1(ii)==0 && ParCrop(ii,jj,1)<winr(2)
                min1(ii)=jj;
            end
%             if min1(ii)>0 && min2(ii)==0 && min1(ii)>jj+b3 && ParCrop(ii,jj,1)<wing(2) && ParCrop(11,jj,1)>wing(1)
%                 min2(ii)=jj;
%             end
        end    
    end
end
for ii=b1:n1-b1
    for jj=n2-b2:-1:1+b2
        if ParCrop(ii,jj,2)<=min([ParCrop(ii,jj-b2:jj-1,2) ParCrop(ii,jj+1:jj+b2,2)])
            if min3(ii)==0 && ParCrop(ii,jj,2)<winb(2)
                min3(ii)=jj;
            end
%             if min1(ii)>0 && min2(ii)==0 && min1(ii)>jj+b3 && ParCrop(ii,jj,1)<wing(2) && ParCrop(11,jj,1)>wing(1)
%                 min2(ii)=jj;
%             end
        end    
    end
end
for ii=b1:n1-b1
    for jj=n2-b2:-1:1+b2
        if ParCrop(ii,jj,3)<=min([ParCrop(ii,jj-b2:jj-1,3) ParCrop(ii,jj+1:jj+b2,3)])
            if min2(ii)==0 && ParCrop(ii,jj,3)<wing(2)
                min2(ii)=jj;
            end
%             if min1(ii)>0 && min2(ii)==0 && min1(ii)>jj+b3 && ParCrop(ii,jj,1)<wing(2) && ParCrop(11,jj,1)>wing(1)
%                 min2(ii)=jj;
%             end
        end    
    end
end

markervec=[];
for ii=1000:n1-1000
    if abs(min2(ii)-min2(ii-1))>15
        markervec=[markervec ii];
    end
end

%%
markervec
min2test=min2;
for jj=1:2:length(markervec)
    min2test(markervec(jj):markervec(jj+1)-1)=min2test(markervec(jj)-1)-(1:markervec(jj+1)-markervec(jj))*(min2test(markervec(jj)-1)-min2test(markervec(jj+1)))/(markervec(jj+1)-markervec(jj)+1)
end
%% Plot of eastern shoreline
figure;
plot(min1,1:length(min1));
hold on;
plot(min3,1:length(min3));
plot(min2test,1:length(min2test));

%% Plot of shoreline regression - ignore ends of island!!
figure;
plot([1:length(min1)]*4,4*(min1-min2test),'k','Linewidth',1.5);
hold on;
plot([1:length(min1)]*4,4*(min1-min3),'b','Linewidth',1.5);
% xlim([235 2809]);
ylim([0 350]*4);

%% I think these are estimates of the change in latitude from the top to bottom of the island
% 37.482916-37.577392 % Parramore
37.181377-37.111183 % Smith

7792/3145

%% I think this plots the elevation of the modeled island
figure
contour(Sv3f);
[m1,m2]=size(Sv3f);
colorbar;
mim3=zeros(m1,1);
mim2=zeros(m1,1);
mim1=zeros(m1,1);
for ii=1:m1
    for jj=m2:-1:2
        if Sv3f(ii,jj)<0 && Sv3f(ii,jj-1)>=0 && mim3(ii)==0
            mim3(ii)=jj-1;
        end
    end
end
for ii=1:m1
    for jj=m2:-1:2
        if Sv2f(ii,jj)<0 && Sv2f(ii,jj-1)>=0 && mim2(ii)==0
            mim2(ii)=jj-1;
        end
    end
end
for ii=1:m1
    for jj=m2:-1:2
        if Sv1f(ii,jj)<0 && Sv1f(ii,jj-1)>=0 && mim1(ii)==0
            mim1(ii)=jj-1;
        end
    end
end

%% 
% figure;
plot([1:m1]*3.08,(mim1-mim3)*3.08,'--b');
hold on;
plot([1:m1]*3.08,(mim1-mim2)*3.08,'--k');
% xlim([449 3809]);
% ylim([0 1000]);

% figure;
% hold on;
% plot(mim1,1:m1);
% plot(mim2,1:m1);
% plot(mim3,1:m1);

3850-447
1400/3.08
1400/4
11000/3.08
11000/4
6200/3.08
6200/4

% 455:2013:3571;
% 350:1550:2750;

%%
disp(strcat(['Parr Mod, 1984-2011 Shoreline Reg., Upper ',num2str(mean(3.08*(mim1(455:2013)-mim3(455:2013)))/27),' m/yr']));
disp(strcat(['Parr Mod, 1984-1998 Shoreline Reg., Upper ',num2str(mean(3.08*(mim1(455:2013)-mim2(455:2013)))/14),' m/yr']));
disp(strcat(['Parr Mod, 1998-2011 Shoreline Reg., Upper ',num2str(mean(3.08*(mim2(455:2013)-mim3(455:2013)))/13),' m/yr']));
disp(strcat(['Parr Obs, 1984-2011 Shoreline Reg., Upper ',num2str(mean(4*(min1(350:1550)-min3(350:1550)))/27),' m/yr']));
disp(strcat(['Parr Obs, 1984-1998 Shoreline Reg., Upper ',num2str(mean(4*(min1(350:1550)-min2(350:1550)))/14),' m/yr']));
disp(strcat(['Parr Obs, 1998-2011 Shoreline Reg., Upper ',num2str(mean(4*(min2(350:1550)-min3(350:1550)))/13),' m/yr']));
disp(strcat(['Parr Mod, 1984-2011 Shoreline Reg., Lower ',num2str(mean(3.08*(mim1(2013:3571)-mim3(2013:3571)))/27),' m/yr']));
disp(strcat(['Parr Mod, 1984-1998 Shoreline Reg., Lower ',num2str(mean(3.08*(mim1(2013:3571)-mim2(2013:3571)))/14),' m/yr']));
disp(strcat(['Parr Mod, 1998-2011 Shoreline Reg., Lower ',num2str(mean(3.08*(mim2(2013:3571)-mim3(2013:3571)))/13),' m/yr']));
disp(strcat(['Parr Obs, 1984-2011 Shoreline Reg., Lower ',num2str(mean(4*(min1(1550:2750)-min3(1550:2750)))/27),' m/yr']));
disp(strcat(['Parr Obs, 1984-1998 Shoreline Reg., Lower ',num2str(mean(4*(min1(1550:2750)-min2(1550:2750)))/14),' m/yr']));
disp(strcat(['Parr Obs, 1998-2011 Shoreline Reg., Lower ',num2str(mean(4*(min2(1550:2750)-min3(1550:2750)))/13),' m/yr']));

%%
(mim1-mim2)*3.08


















%%
I=imread('SmithOutlineReal.jpg');
[nr,nc]=size(I);
SmithResize=imresize(I,[nr*5,nc*1.1]); %*6 for Smith
% [nr,nc]=size(SmithResize);
nrc=00; ncc=00;
SmithCrop=SmithResize;
figure;
% SmithCrop=imcrop(SmithResize,[nrc,ncc,nr-nrc*2,nc-ncc*13.59]);
% imshow(SmithCrop);
imshow(SmithResize);
% ParResize=imresize(I,[nr*3,nc]); %*6 for Smith
% [nr,nc]=size(ParResize);
% nrc=200; ncc=200;
% ParCrop=imcrop(ParResize,[nrc,ncc+190,nr-nrc*8.8,nc-ncc*26.7]); %13.59
% imshow(ParCrop);

%%
figure;
contourf(flipud(SmithCrop(:,:,1)),'Linestyle','none');
colorbar;
%
figure;
contourf(flipud(SmithCrop(:,:,2)),'Linestyle','none');
colorbar;
%
figure;
contourf(flipud(SmithCrop(:,:,3)),'Linestyle','none');
colorbar;

%%
[n1,n2,~]=size(SmithCrop);
b1=4;
b2=10;
b3=20;
winr=[1 140];
wing=[1 160];
winb=[1 70];
min1=zeros(n1,1);
min2=zeros(n1,1);
min3=zeros(n1,1);
for ii=b1:n1-b1
    for jj=n2-b2:-1:1+b2
        if SmithCrop(ii,jj,1)<=min([SmithCrop(ii,jj-b2:jj-1,1) SmithCrop(ii,jj+1:jj+b2,1)])
            if min1(ii)==0 && SmithCrop(ii,jj,1)<winr(2)
                min1(ii)=jj;
            end
%             if min1(ii)>0 && min2(ii)==0 && min1(ii)>jj+b3 && ParCrop(ii,jj,1)<wing(2) && ParCrop(11,jj,1)>wing(1)
%                 min2(ii)=jj;
%             end
        end    
    end
end
% figure;
% plot(min1,1:length(min1));
% hold on;
% % asdf

for ii=b1:n1-b1
    for jj=n2-b2:-1:1+b2
        if SmithCrop(ii,jj,2)<=min([SmithCrop(ii,jj-b2:jj-1,2) SmithCrop(ii,jj+1:jj+b2,2)])
            if min3(ii)==0 && SmithCrop(ii,jj,2)<winb(2)
                min3(ii)=jj;
            end
%             if min1(ii)>0 && min2(ii)==0 && min1(ii)>jj+b3 && ParCrop(ii,jj,1)<wing(2) && ParCrop(11,jj,1)>wing(1)
%                 min2(ii)=jj;
%             end
        end    
    end
end

% min3keep=min3;
min3temp=min3; %min3keep
for ii=228+5:3359-5 %n1-100
    if min3(ii)<mean([min3(ii-5:ii-1); min3(ii+1:ii+5)])-100
        min3temp(ii)=mean([min3(ii-5:ii-1); min3(ii+1:ii+5)])-100;
    end
end
% figure;
% plot(min3temp,1:length(min3temp));
% hold on;
min3=min3temp;
% plot(min3,1:length(min3));
% asdf
markervec=[];
for ii=100:n1-100
    if abs(min3(ii)-min3(ii-1))>40 % && ii>markervec(end)+
        markervec=[markervec ii];
    end
end

markervec
% asdf
markervec=[markervec(2:10) 459 markervec(11:37) markervec(39:end-3)];
markervec=[markervec(1:73) markervec(75:end)];
min3test=min3;
for jj=1:2:length(markervec)
    min3test(markervec(jj):markervec(jj+1)-1)=min3test(markervec(jj)-1)-(1:markervec(jj+1)-markervec(jj))*(min3test(markervec(jj)-1)-min3test(markervec(jj+1)))/(markervec(jj+1)-markervec(jj)+1);
end

%%
wing=[1 160];
min2=zeros(size(min2));
for ii=b1:n1-b1
    for jj=n2-b2:-1:1+b2
        if SmithCrop(ii,jj,3)<=min([SmithCrop(ii,jj-b2:jj-1,3) SmithCrop(ii,jj+1:jj+b2,3)])
            if min2(ii)==0 && SmithCrop(ii,jj,3)<wing(2) && jj<min1(ii)-20
                min2(ii)=jj;
            end
%             if min1(ii)>0 && min2(ii)==0 && min1(ii)>jj+b3 && ParCrop(ii,jj,1)<wing(2) && ParCrop(11,jj,1)>wing(1)
%                 min2(ii)=jj;
%             end
        end    
    end
end
markervec=[];
for ii=100:n1-100
    if abs(min2(ii)-min2(ii-1))>40 % && ii>markervec(end)+
        markervec=[markervec ii];
    end
end
% asdf

markervec
% asdf
markervec=[markervec(7:12) markervec(14:27) markervec(29) markervec(32:33) markervec(36:37) markervec(40) markervec(43:end-1)];
% markervec=[markervec(1:73) markervec(75:end)];
min2test=min2;
for jj=1:2:length(markervec)
    min2test(markervec(jj):markervec(jj+1)-1)=min2test(markervec(jj)-1)-(1:markervec(jj+1)-markervec(jj))*(min2test(markervec(jj)-1)-min2test(markervec(jj+1)))/(markervec(jj+1)-markervec(jj)+1);
end

figure;
plot(min1,1:length(min1));
hold on;
plot(min3test,1:length(min3test));
plot(min2test,1:length(min2test));


%asdf
% for ii=b1:n1-b1
%     for jj=n2-b2:-1:1+b2
%         if SmithCrop(ii,jj,3)<=min([SmithCrop(ii,jj-b2:jj-1,3) SmithCrop(ii,jj+1:jj+b2,3)])
%             if min2(ii)==0 && SmithCrop(ii,jj,3)<wing(2)
%                 min2(ii)=jj;
%             end
% %             if min1(ii)>0 && min2(ii)==0 && min1(ii)>jj+b3 && ParCrop(ii,jj,1)<wing(2) && ParCrop(11,jj,1)>wing(1)
% %                 min2(ii)=jj;
% %             end
%         end    
%     end
% end
% 
% markervec=[];
% for ii=1000:n1-1000
%     if abs(min2(ii)-min2(ii-1))>15
%         markervec=[markervec ii];
%     end
% end

%%
markervec
min2test=min2;
for jj=1:2:length(markervec)
    min2test(markervec(jj):markervec(jj+1)-1)=min2test(markervec(jj)-1)-(1:markervec(jj+1)-markervec(jj))*(min2test(markervec(jj)-1)-min2test(markervec(jj+1)))/(markervec(jj+1)-markervec(jj)+1)
end
%%
figure;
plot(min1,1:length(min1));
hold on;
plot(min3,1:length(min3));
plot(min2test,1:length(min2test));

%%
figure;
plot([1:length(min1)]*2.48,2.48*(min1-min2test),'k','Linewidth',1.5);
hold on;
plot([1:length(min1)]*2.48,2.48*(min1-min3test),'b','Linewidth',1.5);
% xlim([235 2809]);
% ylim([0 350]*4);

%%
37.482916-37.577392

%%
% load(
figure
contour(Sv3f);
[m1,m2]=size(Sv3f);
mim3=zeros(m1,1);
mim2=zeros(m1,1);
mim1=zeros(m1,1);
for ii=1:m1
    for jj=m2:-1:2
        if Sv3f(ii,jj)<0 && Sv3f(ii,jj-1)>=0 && mim3(ii)==0
            mim3(ii)=jj-1;
        end
    end
end
for ii=1:m1
    for jj=m2:-1:2
        if Sv2f(ii,jj)<0 && Sv2f(ii,jj-1)>=0 && mim2(ii)==0
            mim2(ii)=jj-1;
        end
    end
end
for ii=1:m1
    for jj=m2:-1:2
        if Sv1f(ii,jj)<0 && Sv1f(ii,jj-1)>=0 && mim1(ii)==0
            mim1(ii)=jj-1;
        end
    end
end

%%
figure;
plot([1:m1]*4.04,(mim1-mim3)*4.04,'--b');
hold on;
plot([1:m1]*4.04,(mim1-mim2)*4.04,'--k');
% xlim([449 3809]);
% ylim([0 1000]);

% figure;
% hold on;
% plot(mim1,1:m1);
% plot(mim2,1:m1);
% plot(mim3,1:m1);

% 3850-447
% 1400/3.08
% 1400/4
% 11000/3.08
% 11000/4
% 6200/3.08
% 6200/4

% 455:2013:3571;
% 350:1550:2750;

%%
disp(strcat(['Smith Mod, 1984-2011 Shoreline Reg., Upper ',num2str(mean(4.04*(mim1(423:1388)-mim3(423:1388)))/27),' m/yr']));
disp(strcat(['Smith Mod, 1984-1998 Shoreline Reg., Upper ',num2str(mean(4.04*(mim1(423:1388)-mim2(423:1388)))/14),' m/yr']));
disp(strcat(['Smith Mod, 1998-2011 Shoreline Reg., Upper ',num2str(mean(4.04*(mim2(423:1388)-mim3(423:1388)))/13),' m/yr']));
disp(strcat(['Smith Obs, 1984-2011 Shoreline Reg., Upper ',num2str(mean(2.48*(min1(224:1784)-min3(224:1784)))/27),' m/yr']));
disp(strcat(['Smith Obs, 1984-1998 Shoreline Reg., Upper ',num2str(mean(2.48*(min1(224:1784)-min2(224:1784)))/14),' m/yr']));
disp(strcat(['Smith Obs, 1998-2011 Shoreline Reg., Upper ',num2str(mean(2.48*(min2(224:1784)-min3(224:1784)))/13),' m/yr']));
disp(strcat(['Smith Mod, 1984-2011 Shoreline Reg., Lower ',num2str(mean(4.04*(mim1(1388:2354)-mim3(1388:2354)))/27),' m/yr']));
disp(strcat(['Smith Mod, 1984-1998 Shoreline Reg., Lower ',num2str(mean(4.04*(mim1(1388:2354)-mim2(1388:2354)))/14),' m/yr']));
disp(strcat(['Smith Mod, 1998-2011 Shoreline Reg., Lower ',num2str(mean(4.04*(mim2(1388:2354)-mim3(1388:2354)))/13),' m/yr']));
disp(strcat(['Smith Obs, 1984-2011 Shoreline Reg., Lower ',num2str(mean(2.48*(min1(1784:3345)-min3(1784:3345)))/27),' m/yr']));
disp(strcat(['Smith Obs, 1984-1998 Shoreline Reg., Lower ',num2str(mean(2.48*(min1(1784:3345)-min2(1784:3345)))/14),' m/yr']));
disp(strcat(['Smith Obs, 1998-2011 Shoreline Reg., Lower ',num2str(mean(2.48*(min2(1784:3345)-min3(1784:3345)))/13),' m/yr']));
disp(strcat(['Smith Mod, 1984-2011 Shoreline Reg., All ',num2str(mean(4.04*(mim1(423:2354)-mim3(423:2354)))/27),' m/yr']));
disp(strcat(['Smith Mod, 1984-1998 Shoreline Reg., All ',num2str(mean(4.04*(mim1(423:2354)-mim2(423:2354)))/14),' m/yr']));
disp(strcat(['Smith Mod, 1998-2011 Shoreline Reg., All ',num2str(mean(4.04*(mim2(423:2354)-mim3(423:2354)))/13),' m/yr']));
disp(strcat(['Smith Obs, 1984-2011 Shoreline Reg., All ',num2str(mean(2.48*(min1(224:3345)-min3(224:3345)))/27),' m/yr']));
disp(strcat(['Smith Obs, 1984-1998 Shoreline Reg., All ',num2str(mean(2.48*(min1(224:3345)-min2(224:3345)))/14),' m/yr']));
disp(strcat(['Smith Obs, 1998-2011 Shoreline Reg., All ',num2str(mean(2.48*(min2(224:3345)-min3(224:3345)))/13),' m/yr']));

%% Read in 3 years at a time

disp(strcat(['Smith Mod, 0 to 10 Shoreline Reg., Upper ',num2str(mean(4.04*(mim1(423:1388)-mim2(423:1388)))/10),' m/yr']));
disp(strcat(['Smith Mod, 10 to 20 Shoreline Reg., Upper ',num2str(mean(4.04*(mim2(423:1388)-mim3(423:1388)))/10),' m/yr']));

disp(strcat(['Smith Mod, 0 to 10 Shoreline Reg., Lower ',num2str(mean(4.04*(mim1(1388:2354)-mim2(1388:2354)))/10),' m/yr']));
disp(strcat(['Smith Mod, 10 to 20 Shoreline Reg., Lower ',num2str(mean(4.04*(mim2(1388:2354)-mim3(1388:2354)))/10),' m/yr']));

disp(strcat(['Smith Mod, 0 to 10 Shoreline Reg., All ',num2str(mean(4.04*(mim1(423:2354)-mim2(423:2354)))/10),' m/yr']));
disp(strcat(['Smith Mod, 10 to 20 Shoreline Reg., All ',num2str(mean(4.04*(mim2(423:2354)-mim3(423:2354)))/10),' m/yr']));

%%
(mim1-mim2)*3.08

%%
% load(
for SLRvalue=1:4
    if SLRvalue==1
        filenamestring='St00te780BP4000SLRi0p0044872SLRa0MMax1'; % Enter folder name here!!!!
        filenamestring2='St00te780BP4000SLRi0p0044872SLRa0MMax1';
        yearnumvec=[0 15 30];
        islandstring=filenamestring2(1);
        cd /Users/hrogrosky/Downloads/Parameterization_Runs_Code_6_11_2021
        curdir=pwd;

        cd(filenamestring);
        filenamestringmat=strcat([filenamestring2,'t',int2str(yearnumvec(3)),'.mat']);
        load(filenamestringmat);
        cd ..
    elseif SLRvalue==2
        filenamestring='drive_download_20210714T191339Z_Run2'; % Enter folder name here!!!!
        filenamestring2='St00te780BP4000SLRi0p0044872SLRa4p3077en05MMax1';
        yearnumvec=[0 15 30];
        islandstring=filenamestring2(1);
        cd /Users/hrogrosky/Downloads/Parameterization_Runs_Code_6_11_2021
        curdir=pwd;

        cd(filenamestring);
        filenamestringmat=strcat([filenamestring2,'t',int2str(yearnumvec(3)),'.mat']);
        load(filenamestringmat);
        cd ..
    elseif SLRvalue==3
        filenamestring='drive_download_20210714T191339Z_Run3'; % Enter folder name here!!!!
        filenamestring2='St00te780BP4000SLRi0p0044872SLRa0p00010462MMax1';
        yearnumvec=[0 15 30];
        islandstring=filenamestring2(1);
        cd /Users/hrogrosky/Downloads/Parameterization_Runs_Code_6_11_2021
        curdir=pwd;

        cd(filenamestring);
        filenamestringmat=strcat([filenamestring2,'t',int2str(yearnumvec(3)),'.mat']);
        load(filenamestringmat);
        cd ..
    elseif SLRvalue==4
        filenamestring='drive_download_20210714T191402Z_Run4'; % Enter folder name here!!!!
        filenamestring2='St00te780BP4000SLRi0p0044872SLRa0p00015692MMax1';
        yearnumvec=[0 15 30];
        islandstring=filenamestring2(1);
        cd /Users/hrogrosky/Downloads/Parameterization_Runs_Code_6_11_2021
        curdir=pwd;

        cd(filenamestring);
        filenamestringmat=strcat([filenamestring2,'t',int2str(yearnumvec(3)),'.mat']);
        load(filenamestringmat);
        cd ..
    end

%     figure
%     contour(Sv3f);
%     asdf
    [m1,m2]=size(Sv3f);
    mim3=zeros(m1,1);
    mim2=zeros(m1,1);
    mim1=zeros(m1,1);
    for ii=1:m1
        for jj=m2:-1:2
            if Sv3f(ii,jj)<0 && Sv3f(ii,jj-1)>=0 && mim3(ii)==0
                mim3(ii)=jj-1;
            end
        end
    end
    for ii=1:m1
        for jj=m2:-1:2
            if Sv2f(ii,jj)<0 && Sv2f(ii,jj-1)>=0 && mim2(ii)==0
                mim2(ii)=jj-1;
            end
        end
    end
    for ii=1:m1
        for jj=m2:-1:2
            if Sv1f(ii,jj)<0 && Sv1f(ii,jj-1)>=0 && mim1(ii)==0
                mim1(ii)=jj-1;
            end
        end
    end

    if SLRvalue==1
        disp(strcat(['Historic SLR']));
    elseif SLRvalue==2
        disp(strcat(['Low SLR']));
    elseif SLRvalue==3
        disp(strcat(['High SLR']));
    elseif SLRvalue==4
        disp(strcat(['Highest SLR']));
    end
    disp(strcat(['Smith Mod, 2020-2035 Shoreline Reg., Upper ',num2str(mean(4.04*(mim1(423:1388)-mim2(423:1388)))/15),' m/yr']));
    disp(strcat(['Smith Mod, 2035-2050 Shoreline Reg., Upper ',num2str(mean(4.04*(mim2(423:1388)-mim3(423:1388)))/15),' m/yr']));
    disp(strcat(['Smith Mod, 2020-2050 Shoreline Reg., Upper ',num2str(mean(4.04*(mim1(423:1388)-mim3(423:1388)))/30),' m/yr']));
    % disp(strcat(['Smith Obs, 1984-2011 Shoreline Reg., Upper ',num2str(mean(2.48*(min1(224:1784)-min3(224:1784)))/27),' m/yr']));
    % disp(strcat(['Smith Obs, 1984-1998 Shoreline Reg., Upper ',num2str(mean(2.48*(min1(224:1784)-min2(224:1784)))/14),' m/yr']));
    % disp(strcat(['Smith Obs, 1998-2011 Shoreline Reg., Upper ',num2str(mean(2.48*(min2(224:1784)-min3(224:1784)))/13),' m/yr']));
    disp(strcat(['Smith Mod, 2020-2035 Shoreline Reg., Lower ',num2str(mean(4.04*(mim1(1388:2354)-mim2(1388:2354)))/15),' m/yr']));
    disp(strcat(['Smith Mod, 2035-2050 Shoreline Reg., Lower ',num2str(mean(4.04*(mim2(1388:2354)-mim3(1388:2354)))/15),' m/yr']));
    disp(strcat(['Smith Mod, 2020-2050 Shoreline Reg., Lower ',num2str(mean(4.04*(mim1(1388:2354)-mim3(1388:2354)))/30),' m/yr']));
    % disp(strcat(['Smith Obs, 1984-2011 Shoreline Reg., Lower ',num2str(mean(2.48*(min1(1784:3345)-min3(1784:3345)))/27),' m/yr']));
    % disp(strcat(['Smith Obs, 1984-1998 Shoreline Reg., Lower ',num2str(mean(2.48*(min1(1784:3345)-min2(1784:3345)))/14),' m/yr']));
    % disp(strcat(['Smith Obs, 1998-2011 Shoreline Reg., Lower ',num2str(mean(2.48*(min2(1784:3345)-min3(1784:3345)))/13),' m/yr']));
    disp(strcat(['Smith Mod, 2020-2035 Shoreline Reg., All ',num2str(mean(4.04*(mim1(423:2354)-mim2(423:2354)))/15),' m/yr']));
    disp(strcat(['Smith Mod, 2035-2050 Shoreline Reg., All ',num2str(mean(4.04*(mim2(423:2354)-mim3(423:2354)))/15),' m/yr']));
    disp(strcat(['Smith Mod, 2020-2050 Shoreline Reg., All ',num2str(mean(4.04*(mim1(423:2354)-mim3(423:2354)))/30),' m/yr']));
    % disp(strcat(['Smith Obs, 1984-2011 Shoreline Reg., All ',num2str(mean(2.48*(min1(224:3345)-min3(224:3345)))/27),' m/yr']));
    % disp(strcat(['Smith Obs, 1984-1998 Shoreline Reg., All ',num2str(mean(2.48*(min1(224:3345)-min2(224:3345)))/14),' m/yr']));
    % disp(strcat(['Smith Obs, 1998-2011 Shoreline Reg., All ',num2str(mean(2.48*(min2(224:3345)-min3(224:3345)))/13),' m/yr']));
end
