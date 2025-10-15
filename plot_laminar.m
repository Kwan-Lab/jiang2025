%% Scripts to analyze psilocybin-rabies data
% by Alex Kwan, last revised 8/3/2024

clear all;
close all;

maindir = '/Users/alexkwan/Desktop/draft manuscripts/2025 rabies/MATLAB analysis';
lamdatadir = fullfile(maindir,'data_laminar');
lamfigdir = fullfile(maindir,'figs_laminar');

%% set up figure
set(groot, ...
    'DefaultFigureColor', 'w', ...
    'DefaultAxesLineWidth', 2, ...
    'DefaultAxesXColor', 'k', ...
    'DefaultAxesYColor', 'k', ...
    'DefaultAxesFontUnits', 'points', ...
    'DefaultAxesFontSize', 18, ... 
    'DefaultAxesFontName', 'Helvetica', ...
    'DefaultLineLineWidth', 1, ...
    'DefaultTextFontUnits', 'Points', ...
    'DefaultTextFontSize', 18, ... 
    'DefaultTextFontName', 'Helvetica', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.02 0.025],...
    'DefaultLineLineWidth',3);

% set the tickdirs to go out - need to run these lines in this specific order
set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesTickDirMode', 'manual');
set(0,'defaultfigureposition',[40 40 1200 800]);

% colormap
% add current path with underlying folders
addpath(genpath(pwd));

colors=cbrewer('seq','YlOrRd',256);
%colors=flipud(colors);

disp('-----------------------------');
disp('----- Starting analysis -----');
disp('-----------------------------');

%% load the afferent intensity data
txtdata = readmatrix(fullfile(lamdatadir,'23-11-9-Fezf2-increase.xlsx'),'Sheet','Sheet1','Range','1:1','OutputType','string');
ndata = readmatrix(fullfile(lamdatadir,'23-11-9-Fezf2-increase.xlsx'),'Sheet','Sheet1');

distance=ndata(:,1);       %distance from cortical surface

num_inputAll_inc=(size(ndata,2)-1)/2;    %number of regions extracted from Allen Connectivity database
for k=1:num_inputAll_inc
    sig_inputAll_ACAd(:,k)=ndata(:,2*k);  %intensity of fluorescent axons, from each input region examined that terminate in ACAd
    sig_inputAll_MOs(:,k)=ndata(:,2*k+1); %ditto for medial MOs
    name_inputAll{k}=txtdata{2*k}(1:end-6);      %region name
    rabies_dirAll(k)=1;                         %rabies showed increase after psilocybin
end

clear ndata txtdata;

txtdata = readmatrix(fullfile(lamdatadir,'23-11-9-Fezf2-decrease.xlsx'),'Sheet','Sheet1','Range','1:1','OutputType','string');
ndata = readmatrix(fullfile(lamdatadir,'23-11-9-Fezf2-decrease.xlsx'),'Sheet','Sheet1');

distance=ndata(:,1);       %distance from cortical surface

num_inputAll_dec=(size(ndata,2)-1)/2;    %number of regions extracted from Allen Connectivity database
for k=1:num_inputAll_dec
    sig_inputAll_ACAd(:,num_inputAll_inc+k)=ndata(:,2*k);  %intensity of fluorescent axons, from each input region examined that terminate in ACAd
    sig_inputAll_MOs(:,num_inputAll_inc+k)=ndata(:,2*k+1); %ditto for medial MOs
    name_inputAll{num_inputAll_inc+k}=txtdata{2*k}(1:end-6);      %region name
    rabies_dirAll(num_inputAll_inc+k)=-1;                        %rabies showed decrease after psilocybin
end

num_inputAll = num_inputAll_inc + num_inputAll_dec;

clear ndata txtdata;

%% Put in the % difference in inputs measured via rabies tracing
rabies_valAll = [69.7 43.3 43.0 24.4 23.3 17.4 16.4 16.0 16.0 13.5... %SSp-n-R to RSPv-R
                 10.5 10.1 9.50 5.52 5.45 5.30 4.98 4.95 4.86 4.56... %MOs-R to VISam-R
                 4.23 4.09 3.83 2.64 2.42 0.72... %ORBl-L to PF-R
                 -34.0 -31.9 -30.0 -30.0 -29.4 -28.7 -27.9 -27.8 -25.2 -21.2... %AON-R to SMT-R
                 -21.0 -20.5 -20.3 -20.1 -19.6 -17.2 -16.8 -16.7 -15.9 -14.7... %PL-R to ORBvl-L
                 -12.8 -11.2 -8.31 -8.24 -7.49 -5.01 -4.94 -4.11 -3.79 -3.75... %AM-R to CL-R
                 -3.37 -3.35 -2.86 -1.97]; %end TEa-R

%% combine the data for ACAd and MOs

% ACAd has thinner cortical depth, so need to interpolate to extend to
% match MOs
for k=1:num_inputAll
    sig_inputAll_ACAd_interp(:,k) = interp1(sig_inputAll_ACAd(1:701,k),linspace(1,701,801));
end

sig_inputAll = sig_inputAll_ACAd_interp + sig_inputAll_MOs;

%% Estimate of layer positions
% Capture reference atlas from Allen Institute, image 39 of 132
% in Photoshop, measured at boundary of ACAd and MOs
% length of L1: 0.398; L2/3: 0.795; L5: 1.019; L6a: 1.553; L6b: 0.267
% --> unfortunately seems to overestimate thickness of L1, when comparing
% to images of regions known to project to L1

%laminarPortion = [0.0987 0.2959 0.5486 0.9338 1];
%laminarDepth = laminarPortion * size(sig_inputAll,1);

% Estimate by eyeball
laminarPortion = [0.05 0.25 0.5 0.9 1];
laminarDepth = laminarPortion * size(sig_inputAll,1);

%% plot maximal intensity by input region

max_inputAll = max(sig_inputAll,[],1);
max_inputAll(isnan(max_inputAll)) = 0;

% input regions with maximal pixel intensity of
% Low-percentile regions provide very little inputs, so drop for further analysis
thresh = prctile(max_inputAll,20);

figure('Position',[40 40 1200 800]);
plot(1:num_inputAll,max_inputAll,'k.','MarkerSize',30);
hold on;
plot([1 num_inputAll],thresh*[1 1],'k--');
ylabel('Maximal pixel intensity value');
xlabel('Input region');
xticks([1:num_inputAll]);
xticklabels(name_inputAll);

saveas(gcf,fullfile(lamfigdir,'fig_max.png'));
savefig(gcf,fullfile(lamfigdir,'fig_max.fig'));

% only keep input regions that have decent amount of intensity
saveIdx = max_inputAll>thresh & ~isnan(max_inputAll);

sig_inputRegion = sig_inputAll(:,saveIdx);
name_inputRegion = name_inputAll(saveIdx);
rabies_dirRegion = rabies_dirAll(saveIdx);
rabies_valRegion = rabies_valAll(saveIdx);

num_inputRegion = sum(saveIdx);

%% plot afferent intensity data, sort by center of mass

%processing:
for k=1:num_inputRegion
    %1) smooth the intensity data along distance
    %sig = smooth(sig_inputRegion(:,k),50);
    
    %1) cubic smoothing spline
    x = [1:1:size(sig_inputRegion,1)];
    y = sig_inputRegion(:,k)';
    fit_inputRegion(:,k) = csaps(x,y,1e-6,x); %cubic smoothing spline

    %2) normalize for each region by its maximal intensity
    normfit_inputRegion(:,k)=fit_inputRegion(:,k)./max(fit_inputRegion(:,k));
end

% for k=1:num_inputRegion     % find center of mass
%     mass = normfit_inputRegion(:,k);
%     com(k) = nansum(distance.*mass)/nansum(mass);
% end
% [~,idxSort]=sort(com);

for k=1:num_inputRegion     % sort by peak location
    mass = sig_inputRegion(:,k);
    [~,peakloc(k)] = max(mass); 
end
[~,idxSort]=sort(peakloc);

figure('Position',[40 40 1200 800]); 
subplot(4,20,[1:16 21:36 41:56]);
image(1:num_inputRegion,distance./max(distance),normfit_inputRegion(:,idxSort),'CDataMapping','scaled');
colormap(colors);
caxis([0 1]);
hold on;
for kk=1:numel(laminarDepth)-1
    plot([0.5 num_inputRegion+0.5],laminarDepth(kk)*[1 1],'w--','LineWidth',2)
end
ylabel('Normalized cortical depth');
xlabel('Input region');
xticks([1:num_inputRegion]);
xticklabels(name_inputRegion(idxSort));

%make a color scale bar
subplot(3,15,30);
image(0,linspace(0,1,100),linspace(0,1,100)','CDataMapping','scaled');
colormap(colors);
caxis([0 1]);
set(gca,'YDir','normal');
set(gca,'XTick',[]);

set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',18,'LineWidth',2);

saveas(gcf,fullfile(lamfigdir,'fig_afferent.png'));
savefig(gcf,fullfile(lamfigdir,'fig_afferent.fig'));

%% plot example regions

depth = linspace(0,1,size(fit_inputRegion,1));

figure('Position',[40 40 1200 800]);
for kk=1:4
    if kk==1
        example=7; %RSPv-R
    elseif kk==2
        example=2; %VISrl-R
    elseif kk==3
        example=1; %SSp-bfd-R
    elseif kk==4
        example=5; %PO-R
    end

    subplot(1,4,kk); hold on;
    plot(depth,sig_inputRegion(:,example),'k.','MarkerSize',10);
    plot(depth,fit_inputRegion(:,example),'k','LineWidth',5);
    for jj=1:numel(laminarDepth)-1
        plot(laminarPortion(jj)*[1 1],[0 1000],'r--','LineWidth',2)
    end
    title(name_inputRegion{example});
    if kk==1
        xlabel('Normalized cortical depth');
    end
    xlim([0 1]);
    ylabel('Pixel intensity');
    ylim([0 100]);
    view([90 90]);
end
saveas(gcf,fullfile(lamfigdir,'fig_examples1.png'));
savefig(gcf,fullfile(lamfigdir,'fig_examples1.fig'));

%%
figure('Position',[40 40 1200 800]);
for kk=1:4
    if kk==1
        example=24;  %AIv-R
    elseif kk==2
        example=28;  %BLA-R
    elseif kk==3
        example=36;  %CLA-R
    elseif kk==4
        example=37;  %RE-R
    end

    subplot(1,4,kk); hold on;
    plot(depth,sig_inputRegion(:,example),'k.','MarkerSize',10);
    plot(depth,fit_inputRegion(:,example),'k','LineWidth',5);
    for jj=1:numel(laminarDepth)-1
        plot(laminarPortion(jj)*[1 1],[0 1000],'r--','LineWidth',2)
    end
    title(name_inputRegion{example});
    if kk==1
        xlabel('Normalized cortical depth');
    end
    xlim([0 1]);
    ylabel('Pixel intensity');
    ylim([0 100]);
    if kk==2 || kk==4
        ylim([0 300]);
    end
    view([90 90]);
end
saveas(gcf,fullfile(lamfigdir,'fig_examples2.png'));
savefig(gcf,fullfile(lamfigdir,'fig_examples2.fig'));

%% what predicts the psilocybin-evoked input difference

% does the summed axonal density correlate with drug-evoked difference in input?
sum_inputRegion = nansum(sig_inputRegion,1);

figure('Position',[40 40 1200 800]);
subplot(2,2,1);
plot(sum_inputRegion/max(sum_inputRegion),rabies_valRegion,'k.','MarkerSize',20);
xlim([-0.05 1.05]);
xlabel('Summed axonal density (a.u.)');
ylabel('Drug-evoked difference (%)');
axis square;
[r,p]=corrcoef(sum_inputRegion/max(sum_inputRegion),rabies_valRegion);
disp('Drug-evoked difference vs summed axonal density:');
disp(['   correlation coefficient = ' num2str(r(1,2))]);
disp(['   p-value for corr coef = ' num2str(p(1,2))]);

% does the location of peak axonal density correlate with drug-evoked difference in input?
subplot(2,2,2);
plot(peakloc/size(sig_inputRegion,1),rabies_valRegion,'k.','MarkerSize',20);
xlim([-0.05 1.05]);
xlabel('Location of peak axonal density');
axis square;
[r,p]=corrcoef(peakloc/size(sig_inputRegion,1),rabies_valRegion);
disp('Drug-evoked difference vs location of peak axonal density:');
disp(['   correlation coefficient = ' num2str(r(1,2))]);
disp(['   p-value for corr coef = ' num2str(p(1,2))]);

saveas(gcf,fullfile(lamfigdir,'fig_correlations.png'));
savefig(gcf,fullfile(lamfigdir,'fig_correlations.fig'));

%%

% does fraction of summed axonal density residing in layer 1
% correlate with drug-evoked difference in input?

L1_inputRegion = nansum(sig_inputRegion(1:round(laminarDepth(1)),:),1);
L23_inputRegion = nansum(sig_inputRegion(round(laminarDepth(1)):round(laminarDepth(2)),:),1);
L5_inputRegion = nansum(sig_inputRegion(round(laminarDepth(2)):round(laminarDepth(3)),:),1);
L6a_inputRegion = nansum(sig_inputRegion(round(laminarDepth(3)):round(laminarDepth(4)),:),1);

figure('Position',[40 40 1200 800]);
subplot(2,2,1);
plot(100*L1_inputRegion./sum_inputRegion,rabies_valRegion,'k.','MarkerSize',20);
xlim([0 65]);
xlabel('Fraction of summed axonal density in L1 (%)');
ylabel('Drug-evoked difference (%)');
axis square;
[r,p]=corrcoef(100*L1_inputRegion./sum_inputRegion,rabies_valRegion);
disp('Drug-evoked difference vs frac of summed axonal density in L1:');
disp(['   correlation coefficient = ' num2str(r(1,2))]);
disp(['   p-value for corr coef = ' num2str(p(1,2))]);

subplot(2,2,2);
plot(100*L23_inputRegion./sum_inputRegion,rabies_valRegion,'k.','MarkerSize',20);
xlim([0 65]);
xlabel('Fraction of summed axonal density in L2/3 (%)');
axis square;
[r,p]=corrcoef(100*L23_inputRegion./sum_inputRegion,rabies_valRegion);
disp('Drug-evoked difference vs frac of summed axonal density in L2/3:');
disp(['   correlation coefficient = ' num2str(r(1,2))]);
disp(['   p-value for corr coef = ' num2str(p(1,2))]);

subplot(2,2,3);
plot(100*L5_inputRegion./sum_inputRegion,rabies_valRegion,'k.','MarkerSize',20);
xlim([0 65]);
xlabel('Fraction of summed axonal density in L5 (%)');
ylabel('Drug-evoked difference (%)');
axis square;
[r,p]=corrcoef(100*L5_inputRegion./sum_inputRegion,rabies_valRegion);
disp('Drug-evoked difference vs frac of summed axonal density in L5:');
disp(['   correlation coefficient = ' num2str(r(1,2))]);
disp(['   p-value for corr coef = ' num2str(p(1,2))]);

subplot(2,2,4);
plot(100*L6a_inputRegion./sum_inputRegion,rabies_valRegion,'k.','MarkerSize',20);
xlim([0 65]);
xlabel('Fraction of summed axonal density in L6a (%)');
axis square;
[r,p]=corrcoef(100*L6a_inputRegion./sum_inputRegion,rabies_valRegion);
disp('Drug-evoked difference vs frac of summed axonal density in L6a:');
disp(['   correlation coefficient = ' num2str(r(1,2))]);
disp(['   p-value for corr coef = ' num2str(p(1,2))]);
saveas(gcf,fullfile(lamfigdir,'fig_correlations2.png'));
savefig(gcf,fullfile(lamfigdir,'fig_correlations2.fig'));

%%

% does relative axonal density in layer 1
% correlate with drug-evoked difference in input?

L1_inputRegion = nansum(sig_inputRegion(1:round(laminarDepth(1)),:),1);
L23_inputRegion = nansum(sig_inputRegion(round(laminarDepth(1)):round(laminarDepth(2)),:),1);
L5_inputRegion = nansum(sig_inputRegion(round(laminarDepth(2)):round(laminarDepth(3)),:),1);
L6a_inputRegion = nansum(sig_inputRegion(round(laminarDepth(3)):round(laminarDepth(4)),:),1);

% mean value of axonal density in each layer
L1_norm = L1_inputRegion/(round(laminarDepth(1))-1);
L23_norm = L23_inputRegion/(round(laminarDepth(2))-round(laminarDepth(1)));
L5_norm = L5_inputRegion/(round(laminarDepth(3))-round(laminarDepth(2)));
L6a_norm = L6a_inputRegion/(round(laminarDepth(4))-round(laminarDepth(3)));

% fraction based on mean axonal density
L1_frac = 100*L1_norm./(L1_norm+L23_norm+L5_norm+L6a_norm);
L23_frac = 100*L23_norm./(L1_norm+L23_norm+L5_norm+L6a_norm);
L5_frac = 100*L5_norm./(L1_norm+L23_norm+L5_norm+L6a_norm);
L6a_frac = 100*L6a_norm./(L1_norm+L23_norm+L5_norm+L6a_norm);

figure('Position',[40 40 1200 800]);
subplot(2,2,1);
plot(L1_frac,rabies_valRegion,'k.','MarkerSize',20);
xlim([0 65]);
xlabel('Relative density of axons in L1 (%)');
ylabel('Drug-evoked difference (%)');
axis square;
[r,p]=corrcoef(L1_frac,rabies_valRegion);
disp('Drug-evoked difference vs relative axonal density in L1:');
disp(['   correlation coefficient = ' num2str(r(1,2))]);
disp(['   p-value for corr coef = ' num2str(p(1,2))]);

subplot(2,2,2);
plot(L23_frac,rabies_valRegion,'k.','MarkerSize',20);
xlim([0 65]);
xlabel('Relative density of axons in L2/3 (%)');
axis square;
[r,p]=corrcoef(L23_frac,rabies_valRegion);
disp('Drug-evoked difference vs relative axonal density in L2/3:');
disp(['   correlation coefficient = ' num2str(r(1,2))]);
disp(['   p-value for corr coef = ' num2str(p(1,2))]);

subplot(2,2,3);
plot(L5_frac,rabies_valRegion,'k.','MarkerSize',20);
xlim([0 65]);
xlabel('Relative density of axons in L5 (%)');
ylabel('Drug-evoked difference (%)');
axis square;
[r,p]=corrcoef(L5_frac,rabies_valRegion);
disp('Drug-evoked difference vs relative axonal density in L5:');
disp(['   correlation coefficient = ' num2str(r(1,2))]);
disp(['   p-value for corr coef = ' num2str(p(1,2))]);

subplot(2,2,4);
plot(L6a_frac,rabies_valRegion,'k.','MarkerSize',20);
xlim([0 65]);
xlabel('Relative density of axons in L6a (%)');
axis square;
[r,p]=corrcoef(L6a_frac,rabies_valRegion);
disp('Drug-evoked difference vs relative axonal density in L6a:');
disp(['   correlation coefficient = ' num2str(r(1,2))]);
disp(['   p-value for corr coef = ' num2str(p(1,2))]);

saveas(gcf,fullfile(lamfigdir,'fig_correlations3.png'));
savefig(gcf,fullfile(lamfigdir,'fig_correlations3.fig'));
