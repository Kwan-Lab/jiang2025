%% Scripts to analyze psilocybin-rabies data
% by Alex Kwan, last revised 11/12/2024
% with additions from Neil Savalia

clear all;
close all;

maindir = '/Users/alexkwan/Desktop/draft manuscripts/2024 rabies/MATLAB analysis';
rabiesdir = fullfile(maindir,'data_rabies');
datadir = fullfile(maindir,'data_other');
figdir = fullfile(maindir,'figs');

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

% set plot color
PTcolor = [204 116 51]/255;
ITcolor = [86 143 146]/255;
PTcolor_psi = [244 196 71]/255; %following color scheme of Shao, 2024
ITcolor_psi = [126 83 236]/255; %following color scheme of Shao, 2024

%for Fezf2-CreER
% samples _01-_04 male saline, _05-_08 male psilocybin
% samples _09-_012 female saline, _13-_16 female psilocybin
% note: order in spreadsheet in order in which LifeCanvas processed the samples, but not the order by which we named the samples
% plus 1 pilot sample, an earlier test run before the cohort that was a female psilocybin
drug=[0 0 0 1 1 0 1 1 0 0 0 0 1 1 1 1 1]; %0 = saline, 1 = psilocybin
sex=[0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1]; %0 = male, 1 = female
celltype=[1*ones(1,17)]; %1 = PT, 2 = IT

%for PlexinD1-CreER
% samples _01-_04 male saline, _05-_08 male psilocybin
% samples _09-_012 female saline, _13-_16 female psilocybin
drug=[drug 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1]; %0 = saline, 1 = psilocybin
sex=[sex 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1]; %0 = male, 1 = female
celltype=[celltype 2*ones(1,16)]; %1 = PT, 2 = IT

disp('-----------------------------');
disp('----- Starting analysis -----');
disp('-----------------------------');

%% load the starter cell data for Fezf2-CreER mice
ndata = readmatrix(fullfile(rabiesdir,'fezf2_coexpression channel 488_561.xlsx'),'Sheet','All Samples');
region_id_fezf2=ndata(:,1);       %region ID
num_starter=ndata(:,7:22);  %number of starter cells

ndata = readmatrix(fullfile(rabiesdir,'fezf2_pilot-psilocybin-female-Points_in_region_coexpression_488_561.csv'));
num_starter(:,end+1)=ndata(:,6);

txtdata = readtable(fullfile(rabiesdir,'fezf2_coexpression channel 488_561.xlsx'),'Sheet','All Samples');
region_name=table2array(txtdata(:,3));   %region name

clear ndata txtdata;

% Remove item 1638 because Fezf2-CreER raw data contains a region named
% 'tspd-R' that was later absent in PlexinD1-CreER raw data. The region
% contains zero labeled cells anyways
region_id_fezf2=[region_id_fezf2(1:1637); region_id_fezf2(1639:end)];
num_starter=[num_starter(1:1637,:); num_starter(1639:end,:)];
region_name=[region_name(1:1637); region_name(1639:end)];

%% load the starter cell data for PlexinD1-CreER mice
ndata = readmatrix(fullfile(rabiesdir,'plexind1_coexpression channel 488_561.xlsx'),'Sheet','All Samples');
region_id_plexind1=ndata(:,1);       %region ID
num_starter(:,end+1:end+16)=ndata(:,7:22);  %number of starter cells

txtdata = readtable(fullfile(rabiesdir,'plexind1_coexpression channel 488_561.xlsx'),'Sheet','All Samples');
region_name=table2array(txtdata(:,3));   %region name

clear ndata txtdata;

% check that the .xls files have the same region list
if isequal(region_id_fezf2, region_id_plexind1)
    region_id = region_id_fezf2;
else
    error('region_id does not match!');
end

clear region_id_fezf2 region_id_plexind1;

%% plot results about starter cells

%analyze specific regions where starter cells may be
%MOp-R 10018, ACAv-R 10232, MOs-R 10024, ACAd-R 10226
%PL-R 10238, ILA-R 10245, ORBm-R 10264, ORBvl-R 10272
idx = (region_id==10018) + (region_id==10024) + (region_id==10226) + (region_id==10232) + ...
        (region_id==10238) + (region_id==10245) + (region_id==10258) + (region_id==10264) + (region_id==10272);
num_region = sum(idx);

regionList=[]; %names of these regions extracted directly from spreadsheet
temp = find(idx==1);
for j=1:numel(temp)
    regionList{j} = region_name{temp(j)};
end

figure('Position',[40 40 1200 1000]);
for kk=1:2
    % plot total number of starter cells in each mouse
    total_starter=sum(num_starter(idx==1,:),1);

    if kk==1
        color1 = PTcolor;
        color2 = PTcolor_psi;
    elseif kk==2
        color1 = ITcolor;
        color2 = ITcolor_psi;
    end

    subplot(4,2,kk); hold on;
    dat1 = total_starter(drug==0 & celltype==kk);
    p1=bar(-0.2,mean(dat1),0.3,'FaceColor',color1);
    for k=1:numel(dat1)
        plot(-0.2+(0.1*(rand-0.5)),dat1(k),'-o','MarkerEdgeColor',color1,'MarkerFaceColor','w','MarkerSize',10','LineWidth',2);
    end
    dat2 = total_starter(drug==1 & celltype==kk);
    p2=bar(0.2,mean(dat2,2),0.3,'FaceColor',color2);
    for k=1:numel(dat2)
        plot(0.2+(0.1*(rand-0.5)),dat2(k),'-o','MarkerEdgeColor',color2,'MarkerFaceColor','w','MarkerSize',10','LineWidth',2);
    end
    xticks([]);
    ylabel("Starter cells");
    view([90,90]);
    if kk==1
        title('PT (Fezf2-CreER)');
    elseif kk==2
        title('IT (PlexinD1-CreER)');
    end

    % plot proportion of starter cells in each region
    prop_starter = 100*num_starter./total_starter;

    subplot(4,2,[3 5 7]+kk-1); hold on;

    dat1 = prop_starter(idx==1,drug==0 & celltype==kk);
    p1=bar((1:num_region)-0.2,mean(dat1,2),0.3,'FaceColor',color1);
    for j=1:num_region
        for k=1:size(dat1,2)
            plot(j-0.2+(0.1*(rand-0.5)),dat1(j,k),'-o','MarkerEdgeColor',color1,'MarkerFaceColor','w','MarkerSize',10','LineWidth',2);
        end
    end
    dat2 = prop_starter(idx==1,drug==1 & celltype==kk);
    p2=bar((1:num_region)+0.2,mean(dat2,2),0.3,'FaceColor',color2);
    for j=1:num_region
        for k=1:size(dat2,2)
            plot(j+0.2+(0.1*(rand-0.5)),dat2(j,k),'-o','MarkerEdgeColor',color2,'MarkerFaceColor','w','MarkerSize',10','LineWidth',2);
        end
    end
    xticks(1:num_region);
    xticklabels(regionList);
    ylabel("Proportion of starter cells (%)");
    view([90,90]);
    legend([p1 p2],{'Saline','Psilocybin'},'location','northeast');

    disp('Fraction of starter cells, saline group, lies in ACAd, MOs, or PL:')
    if kk==1
        disp(['PT, saline:' num2str(mean(dat1(2,:)+dat1(3,:)+dat1(5,:))) '+-' num2str(std(dat1(2,:)+dat1(3,:)+dat1(5,:))./sqrt(sum(drug==0 & celltype==kk)))]);
        disp(['PT, psilocybin:' num2str(mean(dat2(2,:)+dat2(3,:)+dat2(5,:))) '+-' num2str(std(dat2(2,:)+dat2(3,:)+dat2(5,:))./sqrt(sum(drug==1 & celltype==kk)))]);
    else
        disp(['IT, saline:' num2str(mean(dat1(2,:)+dat1(3,:)+dat1(5,:))) '+-' num2str(std(dat1(2,:)+dat1(3,:)+dat1(5,:))./sqrt(sum(drug==0 & celltype==kk)))]);
        disp(['IT, psilocybin:' num2str(mean(dat2(2,:)+dat2(3,:)+dat2(5,:))) '+-' num2str(std(dat2(2,:)+dat2(3,:)+dat2(5,:))./sqrt(sum(drug==1 & celltype==kk)))]);
    end
end

disp('Number of starter cells');
disp(['PT, saline:' num2str(mean(total_starter(drug==0 & celltype==1))) '+-' num2str(std(total_starter(drug==0 & celltype==1))/sqrt(sum(drug==0 & celltype==1)))]);
disp(['PT, psi:' num2str(mean(total_starter(drug==1 & celltype==1))) '+-' num2str(std(total_starter(drug==1 & celltype==1))/sqrt(sum(drug==1 & celltype==1)))]);
disp(['p = ' num2str(ranksum(total_starter(drug==0 & celltype==1),total_starter(drug==1 & celltype==1)))]);
disp(['IT, saline:' num2str(mean(total_starter(drug==0 & celltype==2))) '+-' num2str(std(total_starter(drug==0 & celltype==2))/sqrt(sum(drug==0 & celltype==2)))]);
disp(['IT, psi:' num2str(mean(total_starter(drug==1 & celltype==2))) '+-' num2str(std(total_starter(drug==1 & celltype==2))/sqrt(sum(drug==1 & celltype==2)))]);
disp(['p = ' num2str(ranksum(total_starter(drug==0 & celltype==2),total_starter(drug==1 & celltype==2)))]);

saveas(gcf,fullfile(figdir,'fig_starter.png'));
savefig(gcf,fullfile(figdir,'fig_starter.fig'));

%% load the input cell data

%Fezf2-CreER
ndata = readmatrix(fullfile(rabiesdir,'fezf2_density channel 488.xlsx'),'Sheet','All Samples');
region_id_fezf2=ndata(:,1);  %region ID
%volume=ndata(:,6);     %volume of region
num_green=ndata(:,7:22);    %number of green cells
ndata = readmatrix(fullfile(rabiesdir,'fezf2_pilot-psilocybin-female-Points_in_region_488_endogenous_model_05_20_21_sensitive_auto_alignment_25_647_Full_trans_FullScale.csv'));
num_green(:,end+1)=ndata(:,6);

% Remove item 1638 because Fezf2-CreER raw data contains a region named
% 'tspd-R' that was later absent in PlexinD1-CreER raw data. The region
% contains zero labeled cells anyways
region_id_fezf2=[region_id_fezf2(1:1637); region_id_fezf2(1639:end)];
num_green=[num_green(1:1637,:); num_green(1639:end,:)];

% check that the .xls files have the same region list
if ~isequal(region_id, region_id_fezf2)
    error('region_id does not match!');
end
clear ndata region_id_fezf2;

%PlexinD1-CreER
ndata = readmatrix(fullfile(rabiesdir,'plexind1_density channel 488.xlsx'),'Sheet','All Samples');
region_id_plexind1=ndata(:,1);  %region ID
num_green(:,end+1:end+16)=ndata(:,7:22);    %number of green cells

% check that the .xls files have the same region list
if ~isequal(region_id, region_id_plexind1)
    error('region_id does not match!');
end
clear ndata region_id_plexind1;

%number of input cells = number of green cells minus number of green+red (starter) cells
num_input = num_green - num_starter;  

clear ndata txtdata;

%% find subset of data belonging to summary structure
% following definition in Wang Q et al., Cell, 2020

ndata = readmatrix(fullfile(datadir,'allen_ccf_csv.csv'));
region_idSummary=ndata(:,1);              %region ID - according to summary list

txtdata = readtable(fullfile(datadir,'allen_ccf_csv.csv'));
region_nameSummary=table2array(txtdata(:,4));   %region name - according to summary list
yesSummary=table2array(txtdata(:,11));           %Y if this is a summary region  - according to summary list

clear ndata txtdata;

% find 316 summary structures and only keep those
yesList = find(strcmp(yesSummary,'Y')); 
region_idSummary = region_idSummary(yesList);
region_nameSummary = region_nameSummary(yesList);

%analyze regions if they are one of the 316 summary structures - for list of inputs: frontal or long-range?
idx_L_frontal = zeros(size(region_id)); %frontal regions only
idx_R_frontal = zeros(size(region_id));
idx_L_global = zeros(size(region_id));  %long-range input regions only
idx_R_global = zeros(size(region_id));
for j=1:numel(region_id)  
    if region_id(j) == 18 || region_id(j) == 232 || region_id(j) == 24 || region_id(j) == 226 ...
             || region_id(j) == 238 || region_id(j) == 245 || region_id(j) == 258 || region_id(j) == 264 || region_id(j) == 272
        idx_L_frontal(j)=1;
    elseif region_id(j) == 10018 || region_id(j) == 10232 || region_id(j) == 10024 || region_id(j) == 10226 ...
             || region_id(j) == 10238 || region_id(j) == 10245 || region_id(j) == 10258 || region_id(j) == 10264 || region_id(j) == 10272
        idx_R_frontal(j)=1;
    elseif region_id(j) == 1101 || region_id(j) == 11101
        disp('Excluding fiber tract in input cell analysis.');
    elseif ~isempty(find(region_idSummary==region_id(j),1)) %if raw data is one of the summary structure
        idx_L_global(j)=1;
    elseif ~isempty(find(region_idSummary==region_id(j)-10000,1)) %if raw data is one of the summary structure, id offset by 10000 if it is on right hemisphere
        idx_R_global(j)=1;
    end
end
idx_L = idx_L_frontal + idx_L_global;   %all regions
idx_R = idx_R_frontal + idx_R_global;   

%somehow this matching process went from 316 summary --> 302 identified from xlsx, what's missing?
%this check shows: LING CENT CUL DEC FOTU PYR UVU NOD SIM AN PRM COPY PFL FL
% xlsnames=[]; %names of these regions extracted directly from spreadsheet
% idx_xls = find(idx_L==1);
% for j=1:numel(idx_xls)
%     xlsnames{j} = region_name{idx_xls(j)}(1:end-2);  %take out the '-L' from left hemisphere list
% end
% for j=1:numel(region_nameSummary)
%     if isempty(find(contains(xlsnames,region_nameSummary{j})))
%         disp(['The xls file is missing ' region_nameSummary{j}]);
%     end
% end

%% plot results about input cells

% plot total number of input cells in each mouse
total_input_R=nansum(num_input(idx_R==1,:),1);
total_input_L=nansum(num_input(idx_L==1,:),1);

figure('Position',[40 40 1400 1000]); 

for kk=1:2

    if kk==1
        color1 = PTcolor;
        color2 = PTcolor_psi;
    elseif kk==2
        color1 = ITcolor;
        color2 = ITcolor_psi;
    end

    subplot(4,2,kk); hold on;
    dat1 = total_input_R(drug==0 & celltype==kk);
    p1=bar(-0.2,mean(dat1),0.3,'FaceColor',color1);
    for k=1:numel(dat1)
        plot(-0.2+(0.1*(rand-0.5)),dat1(k),'-o','MarkerEdgeColor',color1,'MarkerFaceColor','w','MarkerSize',10','LineWidth',2);
    end
    dat2 = total_input_R(drug==1 & celltype==kk);
    p2=bar(0.2,mean(dat2,2),0.3,'FaceColor',color2);
    for k=1:numel(dat2)
        plot(0.2+(0.1*(rand-0.5)),dat2(k),'-o','MarkerEdgeColor',color2,'MarkerFaceColor','w','MarkerSize',10','LineWidth',2);
    end
    xticks([]);
    ylabel("Input cells in right hemisphere");
    ylim([0 6e5]);
    view([90,90]);
    if kk==1
        title('PT (Fezf2-CreER)');
    elseif kk==2
        title('IT (PlexinD1-CreER)');
    end

    subplot(4,2,3+kk-1); hold on;
    dat1 = total_input_L(drug==0 & celltype==kk);
    p1=bar(-0.2,mean(dat1),0.3,'FaceColor',color1);
    for k=1:numel(dat1)
        plot(-0.2+(0.1*(rand-0.5)),dat1(k),'-o','MarkerEdgeColor',color1,'MarkerFaceColor','w','MarkerSize',10','LineWidth',2);
    end
    dat2 = total_input_L(drug==1 & celltype==kk);
    p2=bar(0.2,mean(dat2,2),0.3,'FaceColor',color2);
    for k=1:numel(dat2)
        plot(0.2+(0.1*(rand-0.5)),dat2(k),'-o','MarkerEdgeColor',color2,'MarkerFaceColor','w','MarkerSize',10','LineWidth',2);
    end
    xticks([]);
    ylabel("Input cells in contralateral hemisphere");
    ylim([0 6e5]);
    view([90,90]);

    subplot(2,4,5+2*(kk-1)); hold on;
    plot(total_starter(drug==0 & celltype==kk),total_input_L(drug==0 & celltype==kk)+total_input_R(drug==0 & celltype==kk),'o','MarkerEdgeColor',color1,'MarkerSize',10','LineWidth',2);
    plot(total_starter(drug==1 & celltype==kk),total_input_L(drug==1 & celltype==kk)+total_input_R(drug==1 & celltype==kk),'o','MarkerEdgeColor',color2,'MarkerSize',10','LineWidth',2);
    xlabel("Starter cells");
    ylabel("Input cells");
    xlim([0 10000]);
    ylim([0 8e5]);
    axis square;

    xdata1 = total_starter(drug==0 & celltype==kk);
    ydata1 = (total_input_L(drug==0 & celltype==kk)+total_input_R(drug==0 & celltype==kk))./total_starter(drug==0 & celltype==kk);
    xdata2 = total_starter(drug==1 & celltype==kk);
    ydata2 = (total_input_L(drug==1 & celltype==kk)+total_input_R(drug==1 & celltype==kk))./total_starter(drug==1 & celltype==kk);

    fun = @(x,xdata) x(1)*exp(x(2)*xdata);
    x0 = [500,-0.001];
    x1 = lsqcurvefit(fun,x0,xdata1,ydata1)

    x_range=[0:1:10000];
    subplot(2,4,6+2*(kk-1)); hold on;
    plot(x_range,fun(x1,x_range),'k-','LineWidth',3);
    plot(xdata1,ydata1,'o','MarkerEdgeColor',color1,'MarkerSize',10','LineWidth',2);
    plot(xdata2,ydata2,'o','MarkerEdgeColor',color2,'MarkerSize',10','LineWidth',2);
    xlabel("Starter cells");
    ylabel("Input cells per starter cell");
    xlim([0 10000]);
    ylim([0 400]);
    axis square;

    crit=(xdata2<5000); %exclude far-right data point that is out of fit range
    ydiff=(ydata2(crit)-fun(x1,xdata2(crit)))./fun(x1,xdata2(crit));

    if kk==1
        disp('----- Fezf2-CreER -----');
    elseif kk==2
        disp('----- PlexinD1-CreER -----');
    end
    disp('Input cells per starter cell (excluding one far-right data point)');
    disp('Measured from rabies data, minus expected via fit from saline: ');
    disp([num2str(100*mean(ydiff)) '+-' num2str(100*std(ydiff)/sqrt(numel(ydiff))) '%']);
end


total_input = total_input_L + total_input_R;
disp('Number of input cells');
disp(['PT, saline:' num2str(mean(total_input(drug==0 & celltype==1))) '+-' num2str(std(total_input(drug==0 & celltype==1))/sqrt(sum(drug==0 & celltype==1)))]);
disp(['PT, psi:' num2str(mean(total_input(drug==1 & celltype==1))) '+-' num2str(std(total_input(drug==1 & celltype==1))/sqrt(sum(drug==1 & celltype==1)))]);
disp(['IT, saline:' num2str(mean(total_input(drug==0 & celltype==2))) '+-' num2str(std(total_input(drug==0 & celltype==2))/sqrt(sum(drug==0 & celltype==2)))]);
disp(['IT, psi:' num2str(mean(total_input(drug==1 & celltype==2))) '+-' num2str(std(total_input(drug==1 & celltype==2))/sqrt(sum(drug==1 & celltype==2)))]);

disp('Wilcoxon rank-sum test');
disp(['PT, saline vs psi: P = ' num2str(ranksum(total_input(drug==0 & celltype==1),total_input(drug==1 & celltype==1)))]);
disp(['IT, saline vs psi: P = ' num2str(ranksum(total_input(drug==0 & celltype==2),total_input(drug==1 & celltype==2)))]);

saveas(gcf,fullfile(figdir,'fig_input.png'));
savefig(gcf,fullfile(figdir,'fig_input.fig'));

%% identify regions that exceed a threshold number of red/green cells

% recall starter is defined as red/green cell + located in right frontal cortex
% this is an analysis of red/green cells in the entire brain

num_redgreen = num_starter;

data_redgreen = 100*num_redgreen./(sum(num_redgreen(idx_L==1,:),1)+sum(num_redgreen(idx_R==1,:),1)); % #red+green cells / total red+green cell

y_txt = 'Cells with colocalized red and green fluorescence (%)';
y_range_frontal = [-50 50];
y_range = [-10 10];
thresh = 0.05;    %plot only regions with value > threshold

%local regions (in frontal cortex)
tempL_frontal = find(idx_L_frontal==1);
tempR_frontal = find(idx_R_frontal==1);
idx_L_mod_frontal=[];
idx_R_mod_frontal=[];
for k=1:numel(tempL_frontal)
    %if this is a reasonably strong input in baseline on either hemisphere
    %for either cell type
    if nanmean(data_redgreen(tempL_frontal(k),drug==0 & celltype==1)) >= thresh || nanmean(data_redgreen(tempR_frontal(k),drug==0 & celltype==1)) >= thresh || ...
            nanmean(data_redgreen(tempL_frontal(k),drug==0 & celltype==2)) >= thresh || nanmean(data_redgreen(tempR_frontal(k),drug==0 & celltype==2)) >= thresh
        idx_L_mod_frontal=[idx_L_mod_frontal tempL_frontal(k)];       %keep to plot
        idx_R_mod_frontal=[idx_R_mod_frontal tempR_frontal(k)];       %keep to plot
    end
end
num_region_frontal = numel(idx_L_mod_frontal);
regionList_frontal=[]; %names of these regions extracted directly from spreadsheet
for j=1:num_region_frontal
    regionList_frontal{j} = region_name{idx_L_mod_frontal(j)}(1:end-2);  %take out the '-L' from left hemisphere list
end

%global regions (all other regions outside frontal cortex)
tempL = find((idx_L_global)==1);
tempR = find((idx_R_global)==1);
idx_L_mod_global=[];
idx_R_mod_global=[];
for k=1:numel(tempL)
    %if this is a reasonably strong input in baseline on either hemisphere
    %for either cell type
    if nanmean(data_redgreen(tempL(k),drug==0 & celltype==1)) >= thresh || nanmean(data_redgreen(tempR(k),drug==0 & celltype==1)) >= thresh || ...
            nanmean(data_redgreen(tempL(k),drug==0 & celltype==2)) >= thresh || nanmean(data_redgreen(tempR(k),drug==0 & celltype==2)) >= thresh
        idx_L_mod_global=[idx_L_mod_global tempL(k)];       %keep to plot
        idx_R_mod_global=[idx_R_mod_global tempR(k)];       %keep to plot
    end
end
num_region = numel(idx_L_mod_global);
regionList=[]; %names of these regions extracted directly from spreadsheet
for j=1:num_region
    regionList{j} = region_name{idx_L_mod_global(j)}(1:end-2);  %take out the '-L' from left hemisphere list
end

%% plot results -- red/green cell list

num_redgreen = num_starter;

data_redgreen = 100*num_redgreen./(sum(num_redgreen(idx_L==1,:),1)+sum(num_redgreen(idx_R==1,:),1)); % #red+green cells / total red+green cell

idx1=(drug==0 | drug==1) & (celltype==1 | celltype==2);  %all animals
color1=[0.2 0.2 0.2];

figure('Position',[40 40 800 1000]);
subplot(20,1,[1:6]);
hold on;
dat1 = data_redgreen(idx_L_mod_frontal,idx1);
p1=bar((1:num_region_frontal),-1*mean(dat1,2),0.5,'FaceColor',color1);
for k=1:size(dat1,1)
    plot(k*[1 1],-1*mean(dat1(k,:))+std(dat1(k,:))/sqrt(size(dat1,2))*[-1 1],'-','Color',color1,'LineWidth',3);
end
dat3 = data_redgreen(idx_R_mod_frontal,idx1);
bar((1:num_region_frontal),mean(dat3,2),0.5,'FaceColor',color1);
for k=1:size(dat3,1)
    plot(k*[1 1],mean(dat3(k,:))+std(dat3(k,:))/sqrt(size(dat3,2))*[-1 1],'-','Color',color1,'LineWidth',3);
end
xticks([1:num_region_frontal]);
xticklabels(regionList_frontal);
title('Contralateral                            Ipsilateral');
ylabel(y_txt);
ylim(y_range_frontal);
yticks([-50:10:50]);
yticklabels({'50','40','30','20','10','0','10','20','30','40','50'});
view([90,90]);
saveas(gcf,fullfile(figdir,'fig_starter-comp-list1.png'));
savefig(gcf,fullfile(figdir,'fig_starter-comp-list1.fig'));

%global regions
figure('Position',[40 40 800 1000]);
subplot(20,1,[1:20]);
hold on;
dat1 = data_redgreen(idx_L_mod_global,idx1);
p1=bar((1:num_region),-1*mean(dat1,2),0.5,'FaceColor',color1);
for k=1:size(dat1,1)
    plot(k*[1 1],-1*mean(dat1(k,:))+std(dat1(k,:))/sqrt(size(dat1,2))*[-1 1],'-','Color',color1,'LineWidth',3);
end
dat3 = data_redgreen(idx_R_mod_global,idx1);
bar((1:num_region),mean(dat3,2),0.5,'FaceColor',color1);
for k=1:size(dat3,1)
    plot(k*[1 1],mean(dat3(k,:))+std(dat3(k,:))/sqrt(size(dat3,2))*[-1 1],'-','Color',color1,'LineWidth',3);
end
% gray lines to guide eye
for k=1:floor(num_region/5)
    plot((k*5+0.5)*[1 1],y_range,'Color',[0.5 0.5 0.5],'LineWidth',1);
end
xticks([1:num_region]);
xticklabels(regionList);
title('Contralateral                            Ipsilateral');
ylabel(y_txt);
ylim(y_range);
yticks([-10:5:10]);
yticklabels({'10','5','0','5','10'});
view([90,90]);
saveas(gcf,fullfile(figdir,'fig_starter-comp-list2.png'));
savefig(gcf,fullfile(figdir,'fig_starter-comp-list2.fig'));

%% identify regions that exceed a threshold number of input cells
data = 100*num_input./(total_input_L+total_input_R); % #input cells / total input cell

y_txt = 'Input cells (%)';
y_range_frontal = [-14 14];
y_range = [-4 4];
thresh = 0.3;    %plot only regions with value > threshold

%local regions (in frontal cortex)
tempL_frontal = find(idx_L_frontal==1);
tempR_frontal = find(idx_R_frontal==1);
idx_L_mod_frontal=[];
idx_R_mod_frontal=[];
for k=1:numel(tempL_frontal)
    %if this is a reasonably strong input in baseline on either hemisphere
    %for either cell type
    if nanmean(data(tempL_frontal(k),drug==0 & celltype==1)) >= thresh || nanmean(data(tempR_frontal(k),drug==0 & celltype==1)) >= thresh || ...
            nanmean(data(tempL_frontal(k),drug==0 & celltype==2)) >= thresh || nanmean(data(tempR_frontal(k),drug==0 & celltype==2)) >= thresh
        idx_L_mod_frontal=[idx_L_mod_frontal tempL_frontal(k)];       %keep to plot
        idx_R_mod_frontal=[idx_R_mod_frontal tempR_frontal(k)];       %keep to plot
    end
end
num_region_frontal = numel(idx_L_mod_frontal);
regionList_frontal=[]; %names of these regions extracted directly from spreadsheet
for j=1:num_region_frontal
    regionList_frontal{j} = region_name{idx_L_mod_frontal(j)}(1:end-2);  %take out the '-L' from left hemisphere list
end

%global regions (all other regions outside frontal cortex)
tempL = find((idx_L_global)==1);
tempR = find((idx_R_global)==1);
idx_L_mod_global=[];
idx_R_mod_global=[];
for k=1:numel(tempL)
    %if this is a reasonably strong input in baseline on either hemisphere
    %for either cell type
    if nanmean(data(tempL(k),drug==0 & celltype==1)) >= thresh || nanmean(data(tempR(k),drug==0 & celltype==1)) >= thresh || ...
            nanmean(data(tempL(k),drug==0 & celltype==2)) >= thresh || nanmean(data(tempR(k),drug==0 & celltype==2)) >= thresh
        idx_L_mod_global=[idx_L_mod_global tempL(k)];       %keep to plot
        idx_R_mod_global=[idx_R_mod_global tempR(k)];       %keep to plot
    end
end
num_region = numel(idx_L_mod_global);
regionList=[]; %names of these regions extracted directly from spreadsheet
for j=1:num_region
    regionList{j} = region_name{idx_L_mod_global(j)}(1:end-2);  %take out the '-L' from left hemisphere list
end

%% plot results -- input cell list, various comparisons

for jj=1:3

    if jj==1 %PT vs IT
        idx1=(drug==0 & celltype==1);
        idx2=(drug==0 & celltype==2);
        color1=PTcolor; color2=ITcolor;
        label1='PT, saline'; label2='IT, saline';
    elseif jj==2 %PT, saline vs psilocybin
        idx1=(drug==0 & celltype==1);
        idx2=(drug==1 & celltype==1);
        color1=PTcolor; color2=PTcolor_psi;
        label1='PT, saline'; label2='PT, psilocybin';
    elseif jj==3 %IT, saline vs psilocybin
        idx1=(drug==0 & celltype==2);
        idx2=(drug==1 & celltype==2);
        color1=ITcolor; color2=ITcolor_psi;
        label1='IT, saline'; label2='IT, psilocybin';
    end

    %local regions
    figure('Position',[40 40 800 1000]);
    subplot(20,1,[1:6]);
    hold on;
    dat1 = data(idx_L_mod_frontal,idx1);
    p1=bar((1:num_region_frontal)-0.2,-1*mean(dat1,2),0.3,'FaceColor',color1);
    for k=1:size(dat1,1)
        plot((k-0.2)*[1 1],-1*mean(dat1(k,:))+std(dat1(k,:))/sqrt(size(dat1,2))*[-1 1],'-','Color',color1,'LineWidth',2);
    end
    dat2 = data(idx_L_mod_frontal,idx2);
    p2=bar((1:num_region_frontal)+0.2,-1*mean(dat2,2),0.3,'FaceColor',color2);
    for k=1:size(dat2,1)
        plot((k+0.2)*[1 1],-1*mean(dat2(k,:))+std(dat2(k,:))/sqrt(size(dat2,2))*[-1 1],'-','Color',color2,'LineWidth',2);
    end
    dat3 = data(idx_R_mod_frontal,idx1);
    bar((1:num_region_frontal)-0.2,mean(dat3,2),0.3,'FaceColor',color1);
    for k=1:size(dat3,1)
        plot((k-0.2)*[1 1],mean(dat3(k,:))+std(dat3(k,:))/sqrt(size(dat3,2))*[-1 1],'-','Color',color1,'LineWidth',2);
    end
    dat4 = data(idx_R_mod_frontal,idx2);
    bar((1:num_region_frontal)+0.2,mean(dat4,2),0.3,'FaceColor',color2);
    for k=1:size(dat4,1)
        plot((k+0.2)*[1 1],mean(dat4(k,:))+std(dat4(k,:))/sqrt(size(dat4,2))*[-1 1],'-','Color',color2,'LineWidth',2);
    end
    xticks([1:num_region_frontal]);
    xticklabels(regionList_frontal);
    title('Contralateral                            Ipsilateral');
    ylabel(y_txt);
    ylim(y_range_frontal);
    yticks([-12:3:12]);
    yticklabels({'12','9','6','3','0','3','6','9','12'});
    view([90,90]);
    legend([p1 p2],{label1,label2},'location','southeast');
    saveas(gcf,fullfile(figdir,['fig_input-comp' num2str(jj) '-list1.png']));
    savefig(gcf,fullfile(figdir,['fig_input-comp' num2str(jj) '-list1.fig']));

    %global regions
    figure('Position',[40 40 800 1000]);
    subplot(20,1,[1:20]);
    hold on;
    dat1 = data(idx_L_mod_global,idx1);
    p1=bar((1:num_region)-0.2,-1*mean(dat1,2),0.3,'FaceColor',color1);
    for k=1:size(dat1,1)
        plot((k-0.2)*[1 1],-1*mean(dat1(k,:))+std(dat1(k,:))/sqrt(size(dat1,2))*[-1 1],'-','Color',color1,'LineWidth',2);
    end
    dat2 = data(idx_L_mod_global,idx2);
    p2=bar((1:num_region)+0.2,-1*mean(dat2,2),0.3,'FaceColor',color2);
    for k=1:size(dat2,1)
        plot((k+0.2)*[1 1],-1*mean(dat2(k,:))+std(dat2(k,:))/sqrt(size(dat2,2))*[-1 1],'-','Color',color2,'LineWidth',2);
    end
    dat3 = data(idx_R_mod_global,idx1);
    bar((1:num_region)-0.2,mean(dat3,2),0.3,'FaceColor',color1);
    for k=1:size(dat3,1)
        plot((k-0.2)*[1 1],mean(dat3(k,:))+std(dat3(k,:))/sqrt(size(dat3,2))*[-1 1],'-','Color',color1,'LineWidth',2);
    end
    dat4 = data(idx_R_mod_global,idx2);
    bar((1:num_region)+0.2,mean(dat4,2),0.3,'FaceColor',color2);
    for k=1:size(dat4,1)
        plot((k+0.2)*[1 1],mean(dat4(k,:))+std(dat4(k,:))/sqrt(size(dat4,2))*[-1 1],'-','Color',color2,'LineWidth',2);
    end

    % gray lines to guide eye
    for k=1:floor(num_region/5)
        plot((k*5+0.5)*[1 1],y_range,'Color',[0.5 0.5 0.5],'LineWidth',1);
    end
    xticks([1:num_region]);
    xticklabels(regionList);
    title('Contralateral                            Ipsilateral');
    ylabel(y_txt);
    ylim(y_range);
    yticks([-3:1:3]);
    yticklabels({'3','2','1','0','1','2','3'});
    view([90,90]);
    legend([p1 p2],{label1,label2},'location','southeast');
    saveas(gcf,fullfile(figdir,['fig_input-comp' num2str(jj) '-list2.png']));
    savefig(gcf,fullfile(figdir,['fig_input-comp' num2str(jj) '-list2.fig']));

end

%% plot results -- drug-evoked difference

figure('Position',[40 40 1200 800]);

for kk=1:2

    if kk==1
        col=PTcolor_psi;
    elseif kk==2
        col=ITcolor_psi;
    end

    %fractional difference, drug-evoked
    dat_sal = data(:,drug==0 & celltype==kk);
    dat_psi = data(:,drug==1 & celltype==kk);
    frac_diff = 100*(nanmean(dat_psi,2) - nanmean(dat_sal,2)) ./ nanmean(dat_sal,2);
    %sort based on difference for PT neurons
    if kk==1
        [~,sort_idx] = sort(frac_diff,'descend');
    end

    %only plotting regions above threshold
    sort_idx_mod=[];
    for k=1:numel(sort_idx)
        %if this is a reasonably strong input in baseline for either cell type
        if nanmean(data(sort_idx(k),drug==0 & celltype==1)) >= thresh || nanmean(data(sort_idx(k),drug==0 & celltype==2)) >= thresh        
            if idx_L(sort_idx(k))==1 || idx_R(sort_idx(k))==1  %if it is one of the summary structure
                sort_idx_mod=[sort_idx_mod sort_idx(k)];       %keep to plot
            end
        end
    end
    num_region = numel(sort_idx_mod);

    subplot(2,1,kk)
    hold on;
    bar(1:numel(sort_idx_mod),frac_diff(sort_idx_mod),'FaceColor',col);

    % gray lines to guide eye
    for temp=-75:25:75
        plot([0 num_region+1],temp*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',1);
    end

    xticks([1:num_region]);
    xticklabels(region_name(sort_idx_mod));
    xtickangle(90);     % [NS 20240416] : had to code in the angle re-orientation on my comp, might be different handling for others?
    ylabel(['Drug-evoked difference (%)']);
    yticks([-100:25:100]);
    ylim([-100 100]);

    if kk==1
        title('PT (Fezf2-CreER)');
    elseif kk==2
        title('IT (PlexinD1-CreER)');
    end

    fracdiff_name = region_name(sort_idx_mod);
    fracdiff_val(:,kk) = frac_diff(sort_idx_mod);

    % statistical significance for the drug-evoked difference

    if kk==1
        disp('----- Fezf2-CreER ------');
    elseif kk==2
        disp('----- PlexinD1-CreER ------');
    end
    disp('CI for drug-evoked difference, shuffling saline and psi data:');

    %for each of these region, bootstrap to see if CI of diff is significantly different than 0
    for k=1:numel(sort_idx_mod)
        tempdata_sal = data(sort_idx_mod(k),drug==0 & celltype==kk);
        tempdata_psi = data(sort_idx_mod(k),drug==1 & celltype==kk);
        shuffle_diff=[];
        for j=1:10000   %split - take 3/4 of data set each time
            shuffle_data = tempdata_sal(randperm(numel(tempdata_sal)));
            split_sal = shuffle_data(1:round(3*end/4));
            shuffle_data = tempdata_psi(randperm(numel(tempdata_psi)));
            split_psi = shuffle_data(1:round(3*end/4));
            shuffle_diff(j) = 100*(nanmean(split_psi) - nanmean(split_sal)) ./ nanmean(split_sal);
        end

        %plot error bars on the figure - 90% CI
        ubound = prctile(shuffle_diff,95);
        lbound = prctile(shuffle_diff,5);
        if frac_diff(sort_idx_mod(k)) > 0
            plot([k k],[frac_diff(sort_idx_mod(k)) ubound],'-','Color',col);
        else
            plot([k k],[lbound frac_diff(sort_idx_mod(k))],'-','Color',col);
        end

        ubound(k) = prctile(shuffle_diff,97.5);
        lbound(k) = prctile(shuffle_diff,2.5);
        if 0 >= ubound(k) || 0 <= lbound(k)
            disp([region_name{sort_idx_mod(k)} ' is outside the 95% CI']);
        else
            ubound(k) = prctile(shuffle_diff,95);
            lbound(k) = prctile(shuffle_diff,5);
            if 0 >= ubound(k) || 0 <= lbound(k)
                disp([region_name{sort_idx_mod(k)} ' is outside the 90% CI']);
            end
        end
    end
end

for kk=1:2
    if kk==1
        disp('----- Fezf2-CreER -----');
    elseif kk==2
        disp('----- PlexinD1-CreER -----');
    end
    disp('After thresholding, input cells (%) in saline condition captured by our list of regions:');
    disp([num2str(numel(sort_idx_mod)) ' presynaptic regions, saline group: ' num2str(nansum(mean(data(sort_idx_mod,drug==0 & celltype==kk),2))) '+-' num2str(nanstd(mean(data(sort_idx_mod,drug==0 & celltype==kk),2))/sqrt(sum(drug==0 & celltype==kk))) '%']);
    disp([num2str(numel(sort_idx_mod)) ' presynaptic regions, psilocybin group: ' num2str(nansum(mean(data(sort_idx_mod,drug==1 & celltype==kk),2))) '+-' num2str(nanstd(mean(data(sort_idx_mod,drug==1 & celltype==kk),2))/sqrt(sum(drug==1 & celltype==kk))) '%']);
    %nansum(mean(data(idx_L == 1 | idx_R==1,drug==0 & celltype==kk),2))    %check that total input sum to 100%
    %nansum(mean(data(idx_L == 1 | idx_R==1,drug==1 & celltype==kk),2))    %check that total input sum to 100%
end

saveas(gcf,fullfile(figdir,'fig_input-difference.png'));
savefig(gcf,fullfile(figdir,'fig_input-difference.fig'));

%% plot results -- drug-evoked difference in scatterplot
% scatterplot including both cell types
figure('Position',[40 40 600 600]);
plot(fracdiff_val(:,1),fracdiff_val(:,2),'ko','MarkerSize',10,'LineWidth',2);
hold on;
plot([-80 80],[0 0],'Color',[0.5 0.5 0.5],'LineWidth',1);
plot([0 0],[-80 80],'Color',[0.5 0.5 0.5],'LineWidth',1);
xlabel('Drug-evoked difference, PT (%)');
ylabel('Drug-evoked difference, IT (%)');
xlim([-80 80]);
ylim([-80 80]);
axis square;
saveas(gcf,fullfile(figdir,'fig_input-difference-scatterplot.png'));
savefig(gcf,fullfile(figdir,'fig_input-difference-scatterplot.fig'));

[r,p]=corrcoef(fracdiff_val(:,1),fracdiff_val(:,2));
disp(['Drug-evoked difference, PT vs IT, correlation coefficient = ' num2str(r(1,2))]);
disp(['Drug-evoked difference, PT vs IT, p-value for corr coef = ' num2str(p(1,2))]);

%% input PT vs IT, separated based on cortical network membership

%analyze regions if they are one of the 316 summary structures - for cortex: which network? 
%membership based on Zingg et al., Cell, 2014, Fig. 7C (and to lesser extent 3A, 4A, 5A)

idx_medial = zeros(size(region_id));
idx_medial_contra = zeros(size(region_id));
idx_visaud = zeros(size(region_id));
idx_visaud_contra = zeros(size(region_id));
idx_sensorimotor = zeros(size(region_id));
idx_sensorimotor_contra = zeros(size(region_id));
idx_lateral = zeros(size(region_id));  
idx_lateral_contra = zeros(size(region_id));  
idx_vmpfc = zeros(size(region_id));  
idx_vmpfc_contra = zeros(size(region_id));  

%regions belonging to one of the subnetworks
for j=1:numel(region_id)  
    %Medial: ACAd, ACAv, RSPagl, RSPd, RSPv, VISa, VISrl
    if region_id(j) == 10226 || region_id(j) == 10232 || region_id(j) == 10298 || region_id(j) == 10325 || ...
            region_id(j) == 10332 || region_id(j) == 10346 || region_id(j) == 10353    
        idx_medial(j)=1;
    elseif region_id(j) == 226 || region_id(j) == 232 || region_id(j) == 298 || region_id(j) == 325 || ...
            region_id(j) == 332 || region_id(j) == 346 || region_id(j) == 353    
        idx_medial_contra(j)=1;
    %Visual-auditory: AUDd, AUDp, AUDpo, AUDv, VISal, VISam, VISl, VISp, VISpl, VISpm, VISli, VISpor
    elseif region_id(j) == 10122 || region_id(j) == 10136 || region_id(j) == 10143 || region_id(j) == 10150 || ...
            region_id(j) == 10164 || region_id(j) == 10171 || region_id(j) == 10178 || region_id(j) == 10185 || ...     
            region_id(j) == 10192 || region_id(j) == 10199|| region_id(j) == 10206 || region_id(j) == 10213   
        idx_visaud(j)=1; 
    elseif region_id(j) == 122 || region_id(j) == 136 || region_id(j) == 143 || region_id(j) == 150 || ...
            region_id(j) == 164 || region_id(j) == 171 || region_id(j) == 178 || region_id(j) == 185 || ...     
            region_id(j) == 192 || region_id(j) == 199|| region_id(j) == 206 || region_id(j) == 213   
        idx_visaud_contra(j)=1; 
    %Sensorimotor: MOp, MOs, SSp-n, SSp-bfd, SSp-ll, SSp-m, SSp-ul, SSp-tr, SSp-un, SSs
    elseif region_id(j) == 10018 || region_id(j) == 10024 || region_id(j) == 10044 || region_id(j) == 10051 || ...
            region_id(j) == 10065 || region_id(j) == 10072 || region_id(j) == 10079 || region_id(j) == 10086 || ...
            region_id(j) == 10093 || region_id(j) == 10100     
        idx_sensorimotor(j)=1;
    elseif region_id(j) == 18 || region_id(j) == 24 || region_id(j) == 44 || region_id(j) == 51 || ...
            region_id(j) == 65 || region_id(j) == 72 || region_id(j) == 79 || region_id(j) == 86 || ...
            region_id(j) == 93 || region_id(j) == 100     
        idx_sensorimotor_contra(j)=1;
    %Lateral: GU, VISC, AId, AIp, AIv, TEa, PERI, ECT, PIR, ENTI
    elseif region_id(j) == 10107 || region_id(j) == 10114 || region_id(j) == 10279 || region_id(j) == 10285 || ...
            region_id(j) == 10291 || region_id(j) == 10360 || region_id(j) == 10367 || region_id(j) == 10373 || ...     
            region_id(j) == 10416 || region_id(j) == 10494            
        idx_lateral(j)=1;
    elseif region_id(j) == 107 || region_id(j) == 114 || region_id(j) == 279 || region_id(j) == 285 || ...
            region_id(j) == 291 || region_id(j) == 360 || region_id(j) == 367 || region_id(j) == 373 || ...     
            region_id(j) == 416 || region_id(j) == 494     
        idx_lateral_contra(j)=1;
    %vmPFC: PL, ILA, ORBl, ORBm, ORBvl, AON, TT, DP
    elseif region_id(j) == 10238 || region_id(j) == 10245 || region_id(j) == 10258 || region_id(j) == 10264 || ...
            region_id(j) == 10272 || region_id(j) == 10390 || region_id(j) == 10398 || region_id(j) == 10410
        idx_vmpfc(j)=1;
    elseif region_id(j) == 238 || region_id(j) == 245 || region_id(j) == 258 || region_id(j) == 264 || ...
            region_id(j) == 272 || region_id(j) == 390 || region_id(j) == 398 || region_id(j) == 410
        idx_vmpfc_contra(j)=1;
    end
end

%regions that were summary structures 
idx_right = idx_R_frontal+idx_R_global;
idx_left = idx_L_frontal+idx_L_global;

%plot input % into PT and IT from each network
figure('Position',[40 40 800 500]);
hold on;
for jj=1:6
    if jj==1
        idx=idx_medial;
        label{jj}='Medial'; 
    elseif jj==2
        idx=idx_sensorimotor;
        label{jj}='Sensorimotor';
    elseif jj==3
        idx=idx_visaud;
        label{jj}='Visual-auditory';
    elseif jj==4
        idx=idx_lateral;
        label{jj}='Lateral';
    elseif jj==5
        idx=idx_vmpfc;
        label{jj}='Ventromedial PFC';
    elseif jj==6
        idx=idx_right & ~idx_medial & ~idx_visaud & ~idx_lateral & ~idx_vmpfc & ~idx_sensorimotor;
        label{jj}='Outside neocortex';
    end

    dat1 = nansum(data(idx==1,drug==0 & celltype==1),1);
    dat2 = nansum(data(idx==1,drug==0 & celltype==2),1);
    p1=bar(jj-0.2,mean(dat1),0.3,'FaceColor',PTcolor);
    for k=1:numel(dat1)
        plot(jj-0.2+(0.1*(rand-0.5)),dat1(k),'-o','MarkerEdgeColor',PTcolor,'MarkerFaceColor','w','MarkerSize',10','LineWidth',2);
    end
    p2=bar(jj+0.2,mean(dat2),0.3,'FaceColor',ITcolor);
    for k=1:numel(dat2)
        plot(jj+0.2+(0.1*(rand-0.5)),dat2(k),'-o','MarkerEdgeColor',ITcolor,'MarkerFaceColor','w','MarkerSize',10','LineWidth',2);
    end
    if jj==1
        disp('----- Wilcoxon rank sum test for PT and IT input cells %:');
    end
    [p,~]=ranksum(dat1,dat2);
    disp(['PT:' label{jj} ': input cells % ' num2str(mean(dat1)) '+-' num2str(std(dat1)/sqrt(numel(dat1)))]);
    disp(['IT:' label{jj} ': input cells % ' num2str(mean(dat2)) '+-' num2str(std(dat2)/sqrt(numel(dat2)))]);
    disp(['PT vs IT:' label{jj} ': input cells %, p = ' num2str(p)]);
end
xticks(1:jj);
xticklabels(label);
ylabel("Input cells (%)");
yticks([0:10:35]);
ylim([0 35]);
view([90,90]);
legend([p1 p2],{'PT, saline','IT, saline'},'location','southeast');

saveas(gcf,fullfile(figdir,'fig_networks.png'));
savefig(gcf,fullfile(figdir,'fig_networks.fig'));

%% input saline vs drug-evoked differences, separated based on cortical network membership

disp('-------------------------------');
figure('Position',[40 40 800 500]);
hold on;

for jj=1:6
    if jj==1
        idx=idx_medial;
        label{jj}='Medial'; 
    elseif jj==2
        idx=idx_sensorimotor;
        label{jj}='Sensorimotor';
    elseif jj==3
        idx=idx_visaud;
        label{jj}='Visual-auditory';
    elseif jj==4
        idx=idx_lateral;
        label{jj}='Lateral';
    elseif jj==5
        idx=idx_vmpfc;
        label{jj}='Ventromedial PFC';
    elseif jj==6
        idx=idx_right & ~idx_medial & ~idx_lateral & ~idx_vmpfc & ~idx_sensorimotor & ~idx_visaud;
        label{jj}='Outside neocortex';
    end

    dat_sal_PT = nansum(data(idx==1,drug==0 & celltype==1),1);
    dat_psi_PT = nansum(data(idx==1,drug==1 & celltype==1),1);
    dat_PT = 100*(mean(dat_psi_PT) - mean(dat_sal_PT)) ./ mean(dat_sal_PT);
    p1=bar(jj-0.2,dat_PT,0.3,'FaceColor',PTcolor_psi);
    dat_sal_IT = nansum(data(idx==1,drug==0 & celltype==2),1);
    dat_psi_IT = nansum(data(idx==1,drug==1 & celltype==2),1);
    dat_IT = 100*(mean(dat_psi_IT) - mean(dat_sal_IT)) ./ mean(dat_sal_IT);
    p2=bar(jj+0.2,dat_IT,0.3,'FaceColor',ITcolor_psi);

    %for each of these network, bootstrap to see if CI of PT vs IT diff is significantly different than 0
    if jj==1
        disp('----- CI for difference in PT and IT drug-evoked differences, shuffling saline and psi data:');
    end
    disp(['PT:' label{jj} ': input cells % ' num2str(dat_PT)]);
    disp(['IT:' label{jj} ': input cells % ' num2str(dat_IT)]);

    shuffle_diff=[];
    for j=1:10000  %split - take 3/4 of data set each time
        shuffle_data = dat_sal_PT(randperm(numel(dat_sal_PT)));
        split_sal = shuffle_data(1:round(3*end/4));
        shuffle_data = dat_psi_PT(randperm(numel(dat_psi_PT)));
        split_psi = shuffle_data(1:round(3*end/4));
        shuffle_diff_PT(j) = 100*(nanmean(split_psi) - nanmean(split_sal)) ./ nanmean(split_sal);

        shuffle_data = dat_sal_IT(randperm(numel(dat_sal_IT)));
        split_sal = shuffle_data(1:round(3*end/4));
        shuffle_data = dat_psi_IT(randperm(numel(dat_psi_IT)));
        split_psi = shuffle_data(1:round(3*end/4));
        shuffle_diff_IT(j) = 100*(nanmean(split_psi) - nanmean(split_sal)) ./ nanmean(split_sal);

        shuffle_diff(j) = shuffle_diff_PT(j)-shuffle_diff_IT(j);
    end

    %plot error bars on the figure - 90% CI
    ubound = prctile(shuffle_diff_PT,95);
    lbound = prctile(shuffle_diff_PT,5);
    if dat_PT > 0
        plot([jj-0.2 jj-0.2],[dat_PT ubound],'-','Color',PTcolor_psi);
    else
        plot([jj-0.2 jj-0.2],[lbound dat_PT],'-','Color',PTcolor_psi);
    end
    ubound = prctile(shuffle_diff_IT,95);
    lbound = prctile(shuffle_diff_IT,5);
    if dat_IT > 0
        plot([jj+0.2 jj+0.2],[dat_IT ubound],'-','Color',ITcolor_psi);
    else
        plot([jj+0.2 jj+0.2],[lbound dat_IT],'-','Color',ITcolor_psi);
    end

    %output the results of the statistcal tests
    ubound = prctile(shuffle_diff,97.5);
    lbound = prctile(shuffle_diff,2.5);
    if 0 >= ubound || 0 <= lbound
        disp(['PT vs IT:' label{jj} ': drug-evoked difference is outside the 95% CI']);
    else
        ubound = prctile(shuffle_diff,95);
        lbound = prctile(shuffle_diff,5);
        if 0 >= ubound || 0 <= lbound
            disp(['PT vs IT:' label{jj} ': drug-evoked difference is outside the 90% CI']);
        end
    end
end

xticks(1:jj);
xticklabels(label);
ylabel("Drug-evoked difference (%)");
yticks([-30:10:30]);
ylim([-35 35]);
view([90,90]);
legend([p1 p2],{'PT','IT'},'location','southeast');

saveas(gcf,fullfile(figdir,'fig_networks2.png'));
savefig(gcf,fullfile(figdir,'fig_networks2.fig'));

%% run a chi-squared test
% the different networks as categories
% test whether increased inputs and decreased inputs are distributed randomly or
% with bias among the different categories

for l = 1:2
    inputdiff=nan(1,numel(sort_idx_mod));   %1 if increase input, -1 if decrease input
    networkmem=nan(1,numel(sort_idx_mod)); %1 if medial, 2 if visaud, 3 if sensorimotor, 4 if lateral, 5 if vmpfc

    for j = 1:numel(sort_idx_mod)   % for the presynaptic regions

        if fracdiff_val(j,l)>0  % l=1 for PT cell, 2 for IT cells
            inputdiff(j) = 1;
        else
            inputdiff(j) = -1;
        end

        if idx_medial(sort_idx_mod(j))==1 || idx_medial_contra(sort_idx_mod(j))==1
            networkmem(j) = 1;
        elseif idx_visaud(sort_idx_mod(j))==1 || idx_visaud_contra(sort_idx_mod(j))==1
            networkmem(j) = 2;
        elseif idx_sensorimotor(sort_idx_mod(j))==1 || idx_sensorimotor_contra(sort_idx_mod(j))==1
            networkmem(j) = 3;
        elseif idx_lateral(sort_idx_mod(j))==1 || idx_lateral_contra(sort_idx_mod(j))==1
            networkmem(j) = 4;
        elseif idx_vmpfc(sort_idx_mod(j))==1 || idx_vmpfc_contra(sort_idx_mod(j))==1
            networkmem(j) = 5;
        end
    end

    % calculation of chi-square
    numInc = sum(inputdiff(~isnan(networkmem))==1);
    numDec = sum(inputdiff(~isnan(networkmem))==-1);
    for j = 1:5  %for each network
        E1(j) = numInc/(numInc+numDec)*sum(networkmem==j); % expected numbers for input increase in each network
        E2(j) = numDec/(numInc+numDec)*sum(networkmem==j);
        OE1(j) = sum(inputdiff==1 & networkmem==j) - E1(j); % observed minus expected
        OE2(j) = sum(inputdiff==-1 & networkmem==j) - E2(j);
        OEsqE1(j) = OE1(j)^2/E1(j); % (O-E)^2/E
        OEsqE2(j) = OE2(j)^2/E2(j); % (O-E)^2/E
    end
    % check: sum(E1)+sum(E2) == sum(~isnan(networkmem))
    % check, should be all close to zero: OE1 + OE2
    % check: sum(OE1) == 0
    % check: sum(OE2) == 0

    chisq = sum(OEsqE1)+sum(OEsqE2); 
    df = (2-1)*(5-1);
    %24.5777, df = 4, p-value of 0.000061 for PT neurons, using https://www.socscistatistics.com/pvalues/chidistribution.aspx
    %15.0044, df = 4, p-value of 0.004693 for IT neurons
    
    % bootstrap to generate null distribution
    for k = 1:1000000
        draw_boot = randsample(numel(sort_idx_mod),numel(sort_idx_mod));
        networkmem_boot = networkmem(draw_boot);

        % calculation of chi-square
        numInc = sum(inputdiff(~isnan(networkmem_boot))==1);
        numDec = sum(inputdiff(~isnan(networkmem_boot))==-1);
        for j = 1:5  %for each network
            E1(j) = numInc/(numInc+numDec)*sum(networkmem_boot==j); % expected numbers for input increase in each network
            E2(j) = numDec/(numInc+numDec)*sum(networkmem_boot==j);
            OE1(j) = sum(inputdiff==1 & networkmem_boot==j) - E1(j); % observed minus expected
            OE2(j) = sum(inputdiff==-1 & networkmem_boot==j) - E2(j);
            OEsqE1(j) = OE1(j)^2/E1(j); % (O-E)^2/E
            OEsqE2(j) = OE2(j)^2/E2(j); % (O-E)^2/E
        end

        chisq_boot(k) = sum(OEsqE1)+sum(OEsqE2);
    end

    if l == 1
        figure;
    end
    subplot(2,2,l); hold on;
    h1 = histogram(chisq_boot,[0:1:30],'FaceColor','k','EdgeColor','k');
    h1.Normalization = 'probability';
    plot(chisq*[1 1],[0 0.2],'k','LineWidth',3);
    xlabel('Network selectivity (\chi^2)');
    ylabel('Probability');
    if l == 1
        title('PT neurons');
    elseif l == 2
        title('IT neurons');
    end
end

saveas(gcf,fullfile(figdir,'fig_chisqtest.png'));
savefig(gcf,fullfile(figdir,'fig_chisqtest.fig'));

%% Drug-evoked difference vs. transcriptional in situ hybridization (from Fulcher et al. 2019 PNAS)
% this section of code written and added by NKS

disp('--- paused, before going onto analyses involving transcripts ---')
pause;

% add current path with underlying folders
addpath(genpath(pwd));

% cortical parcellation colors
reg_cmap = [255, 165, 0;         % PFC: orange 
            131, 57, 146;        % Anterolateral: purple
            102,220,169;         % Medial: aqua green
            200, 144, 212;       % Temporal: pink
            22,106,223;          % Somatomotor: medium blue
            196,206,80] / 255;   % Visual: yellow
reg_corder = {'prefrontal', 'anterolateral', 'medial', 'temporal', 'somatomotor', 'visual'};    % rgb order; needed to get/set rgb per dot in downstratm scatterplot

 % load regional energy values (transcript in situ hybridization from Fulcher et al 2019 PNAS)
load(fullfile(datadir,'mouseGradients','Data','AllenGeneDataset_40_19419.mat'));   % cortical parcels, lamina combined
ged = GeneExpData;  % using z-scored energy values below (ged.combZ.energy), but can look at raw (ged.comb.energy)
names_subreg = structInfo.acronym;
names_genes = geneInfo.acronym;
clear GeneExpData structInfo geneInfo

% Assign appropriate colors to the appropriate order/positions specified by Region Names
[names_region_sorted, ~, ~] = LabelCorticalAreas(names_subreg);     % script @ ./MouseGradients/DataProcessing/LabelCorticalAreas.m
cmap_sorted = [];
for i = 1:size(names_region_sorted) %goes 1-40
    cmap_sorted(i, :) = reg_cmap(logical(ismember(reg_corder, names_region_sorted(i))), :);
end
clear names_region_sorted reg_corder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% list of transcripts to plot %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   - can try GOI_list = any group of 3 from the following:   
%       - any of the following 5-HT receptors: 
%           {'Htr1a', 'Htr2a', 'Htr2c'}
%       - main GABAergic subtypes: 
%           {'Pvalb', 'Sst', 'Vip'}

GOI_list = {'Htr1a', 'Htr2a', 'Htr2c'}    % main 5-HT receptors of interest 
% GOI_list = {'Pvalb', 'Sst', 'Vip'}      % main GABAergic interneuron subtypes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% celltype x transcript(s) plots %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop through plotting for celltype x transcript
for gg = 1:numel(GOI_list)

    % setup the summary plot
    if gg == 1
        ff = figure('color', 'w', 'position', [100, 1000, 400 * numel(GOI_list), 720]);
    end
   
    % Select gene of interest
    GOI = GOI_list{gg};
    switch GOI          % some unecessary break-out/hard-codings below... just leaving for now... if one wanted to do combinations, ratios, etc this is where they would do it.
        % any available 5-ht receptor
        case {'Htr1a', 'Htr1b', 'Htr1d', 'Htr1f', 'Htr2a', 'Htr2b', 'Htr2c', 'Htr3a', 'Htr4', 'Htr5a', 'Htr5b', 'Htr6', 'Htr7'}
            values_GOI = (ged.combZ.energy(:, find(ismember(names_genes, GOI))));
        % main GABAergic subtypes
        case {'Pvalb', 'Sst', 'Vip'}
            values_GOI = (ged.combZ.energy(:, find(ismember(names_genes, GOI))));
        % catch all for exceptions
        otherwise
            values_GOI = (ged.combZ.energy(:, find(ismember(names_genes, GOI))));
    end
   
    % loop through cell type plotting 
    for kk = 1:2    % kk = 1 is pt, kk = 2 is it

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sub-select region data and names %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % calculate input cells data
        data = 100 * num_input ./ (total_input_L + total_input_R);  % #input cells / total input cell

        % set threshold (maybe don't threshold so we can just see all the regions available???)
        thresh = 0.3;    %plot only regions with value > threshold

        % get rabies regions
        %   - following ACK strategy above to subselect regions based on thresh
        %   - additionally need to combine L & R hemispheres to match ISH data
        tempL_all = find(idx_L == 1);
        tempR_all = find(idx_R == 1);
        idx_L_mod_all = [];
        idx_R_mod_all = [];
        for k = 1:numel(tempL_all)
            %if this is a reasonably strong input in baseline on either hemisphere for either cell type
            if nanmean(data(tempL_all(k), drug == 0 & celltype == 1)) >= thresh || nanmean(data(tempR_all(k), drug == 0 & celltype == 1)) >= thresh || ...
                    nanmean(data(tempL_all(k), drug == 0 & celltype == 2)) >= thresh || nanmean(data(tempR_all(k), drug == 0 & celltype == 2)) >= thresh
                idx_L_mod_all = [idx_L_mod_all tempL_all(k)];       %keep to plot
                idx_R_mod_all = [idx_R_mod_all tempR_all(k)];       %keep to plot
            end
        end
        num_region_LR = numel(idx_L_mod_all);  % now get a region list ignoring hemispheres
        regionList_LR = []; %names of these regions extracted directly from spreadsheet
        for j = 1:num_region_LR
            regionList_LR{j} = region_name{idx_L_mod_all(j)}(1:end-2);  %take out the '-L' from left hemisphere list
        end
        regionList_LR = regionList_LR';

        % calculate fractional difference, drug-evoked (sum across hemispheres)
        data_LR = data(idx_L_mod_all, :) + data(idx_R_mod_all, :);  % use sum
    %     data_LR = (data(idx_L_mod_all, :) + data(idx_R_mod_all, :)) ./ 2;  % use avg
        dat_sal = data_LR(:, drug == 0 & celltype == kk);
        dat_psi = data_LR(:, drug == 1 & celltype == kk);
        fracdiff = 100 * (nanmean(dat_psi, 2) - nanmean(dat_sal, 2)) ./ nanmean(dat_sal, 2);

        % arrange rabies input regions based on ISH region list/organization
        fracdiff_select = nan(size(values_GOI, 1), 1);
        for i = 1:size(fracdiff_select, 1)

            % current region
            if sum(strcmp(regionList_LR, names_subreg(i))) > 0
                cr = find(strcmp(regionList_LR, names_subreg(i)));
                fracdiff_select(i, :) = fracdiff(cr, :);
            else    % could just skip since init matrix with nan, but leaving for now...
                fracdiff_select(i, :) = NaN; 
            end

        end
        clear cr

        %%%%%%%%
        % plot %
        %%%%%%%%
        pnum = (numel(GOI_list) * (kk - 1) + gg);
        subplot(2, numel(GOI_list), pnum);

        % scatter
        xd = values_GOI;
        yd = nanmean(fracdiff_select, 2);
        splot = scatter(xd, yd, 225, cmap_sorted, 'filled', 'MarkerEdgeColor', 'k');

        % get & set dot labels
        for getset_dot_labels = 1

            tpos_x = xd + (.03 * range(xd));
            tpos_y = yd + (.05 * range(yd));

            % apply labels
            for i = 1:size(xd) %1-40
                DesiredColors = cmap_sorted(i,:);
                text(tpos_x(i), tpos_y(i), names_subreg(i), 'color', DesiredColors, 'FontSize', 12);
            end

        end
        clear getset_dot_labels tpos_* DesiredColors

        % axis labels and titles
        xlabel(['\it' GOI '\rm ISH energy (z)']);
        ylabel(['Drug-evoked difference (%)']);

        % ticks
        yticks([-80:40:80]);
        yy = [-80, 80];
        ylim(yy);

        if kk==1
            title('PT (Fezf2-CreER)');
        elseif kk==2
            title('IT (PlexinD1-CreER)');
        end

        % Set axes and x/y labels
        for getset_xy_labels = 1

            % axis positions and labeling are done together
            switch GOI
                case {'Htr1a'}
                    xx = [-2, 4];
                case {'Htr2a', 'Htr2c'}
                    xx = [-2, 3];
                case {'Pvalb'}
                    xx = [-2, 3];
                case {'Sst'}
                    xx = [-4, 2];
                case {'Vip'}
                    xx = [-2, 3];
                otherwise
                    xx = [(round(min(xd) * 2) / 2) - 1, (round(max(xd) * 2) / 2) + 1];
            end
            xlim(xx);
            xticks(-10:2:10);

        end
        clear getset_xy_labels

        % Stats
        [r, p] = corrcoef(xd, yd, 'rows', 'complete');
        if kk == 1
            cn = 'PT';
        else
            cn = 'IT';
        end
        disp(['Drug-evoked input difference vs. ' GOI ' expression, ' cn ' neurons, r = ' num2str(r(1,2)) ', P-value = ' num2str(p(1,2))]);
        
        % print stats to plot
        if p >= .001
            text((xx(1) + 0.03 * range(xx)), (yy(1) + 0.05 * range(yy)), ['{\itr} = ' sprintf('%.3f', r(1,2)) ', {\itP} = ' sprintf('%.3f', p(1,2))], 'FontSize', 16);
        else
            text((xx(1) + 0.03 * range(xx)), (yy(1) + 0.05 * range(yy)), ['{\itr} = ' sprintf('%.3f', r(1,2)) ', {\itP} < .001'], 'FontSize', 16);
        end

        %%%%%%%%%%%%
        % clean up %
        %%%%%%%%%%%%

        clear i j k
        clear thresh tempL_all tempR_all idx_L_mod_all idx_R_mod_all
        clear num_region_LR regionList_LR 
        clear data data_LR fracdiff fracdiff_select
        clear splot xd yd tpos_x tpos_
        clear r p cn pnum xx yy

    end
    clear kk ff

end
clear gg

%%%%%%%%%%%%%
% save plot %
%%%%%%%%%%%%%

% compile name of transcripts as '_' separated string
tstring = [];
for ss = 1:size(GOI_list, 2)
    if ss == 1
        tstring = sprintf('%s', GOI_list{ss});
    else
        tstring = sprintf('%s_%s', tstring, GOI_list{ss});
    end
end

% % save the file
saveas(gcf,fullfile(figdir,['fig_input-difference-by-ish-scatterplot_' tstring '.png']));
savefig(gcf,fullfile(figdir,['fig_input-difference-by-ish-scatterplot_' tstring '.fig']));

% % clean up
% clear GOI GOI_list values_GOI cmap_sorted ged
% clear dat_psi dat_sal
% clear names_subreg names_genes

%% Drug-evoked rabies input difference vs. drug-evoked cfos differences (from Davoudian et al. ACS Chem Neuro)
% WORKING 20240417

% add current path with underlying folders
addpath(genpath(pwd));

% cortical parcellation colors
reg_cmap = [255, 165, 0;         % PFC: orange 
            131, 57, 146;        % Anterolateral: purple
            102,220,169;         % Medial: aqua green
            200, 144, 212;       % Temporal: pink
            22,106,223;          % Somatomotor: medium blue
            196,206,80] / 255;   % Visual: yellow
reg_corder = {'prefrontal', 'anterolateral', 'medial', 'temporal', 'somatomotor', 'visual'};    % rgb order; needed to get/set rgb per dot in downstratm scatterplot

% load fulcher ish data (used to filter region names for consistent cortical parcellation color scheme)
load(fullfile(datadir,'mouseGradients','Data','AllenGeneDataset_40_19419.mat'));   % cortical parcels, lamina combined
names_subreg_fulcher = structInfo.acronym;
clear GeneExpData structInfo geneInfo

 % load regional cfos light sheet counts (from Davoudian et al ACS Chem Neuro)
cfostab = readtable(fullfile(datadir,'cfos_davoudian','suppTable2.csv'));   % parcels fos x animal
names_subreg = cfostab.abbreviation;

% normalize raw values within-mouse by brain-wise sum (within-column)
for cc = find(contains(cfostab.Properties.VariableNames, 'Count'))  % for data columns... variable names has word "Count"
    evalc(['cfostab.' cfostab.Properties.VariableNames{cc} ' = table2array(cfostab(:, ' num2str(cc) ')) ./ sum(table2array(cfostab(:, ' num2str(cc) ')), 1);']);
end
clear cc

% filter regions so consistent with fulcher labeling code (consistency needed for LabelCorticalAreas.m function call below)
names_filter = zeros(size(names_subreg));
for k = 1:numel(names_subreg_fulcher)
    names_filter(strcmp(names_subreg, names_subreg_fulcher{k})) = 1;
end
clear k

% cfos table & names, after normalizing and filtering down the regions.
cfostab = cfostab(logical(names_filter), :);
names_subreg = names_subreg(logical(names_filter));
clear names_filter

% Assign appropriate colors based on order/positions specified by names_subreg
[names_region_sorted, ~, ~] = LabelCorticalAreas(names_subreg);     % script @ ./MouseGradients/DataProcessing/LabelCorticalAreas.m
cmap_sorted = [];
for i = 1:size(names_region_sorted)     %full fulcher list goes 1-40
    cmap_sorted(i, :) = reg_cmap(logical(ismember(reg_corder, names_region_sorted(i))), :);
end
clear names_region_sorted reg_corder

%%%%%%%%%%%%%%%%%%%%%%%%%
% rabies vs. cfos plots %
%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate fractional difference between treatment (psilo) and baseline (saline)
cfos_tid = {'PSI'};     % "treatment" ID for finding psilocybin-treated mice
cfos_bid = {'SAL'};     % "baseline" ID for finding saline-treated mice
cfos_t = mean(table2array(cfostab(:, contains(cfostab.Properties.VariableNames, cfos_tid))), 2);
cfos_b = mean(table2array(cfostab(:, contains(cfostab.Properties.VariableNames, cfos_bid))), 2);
fracdiff_cfos = 100 * (nanmean(cfos_t, 2) - nanmean(cfos_b, 2)) ./ nanmean(cfos_b, 2);
    
% loop through cell type plotting 
for kk = 1:2    % kk = 1 is pt, kk = 2 is it
    
    % setup the summary plot
    if kk == 1
        ff = figure('color', 'w', 'position', [100, 1000, 400, 720]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sub-select region data and names %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % calculate input cells data
    data = 100 * num_input ./ (total_input_L + total_input_R);  % #input cells / total input cell

    % set threshold (maybe don't threshold so we can just see all the regions available???)
    thresh = 0.3;    %plot only regions with value > threshold

    % get rabies regions
    %   - following ACK strategy above to subselect regions based on thresh
    %   - additionally need to combine L & R hemispheres to match ISH data
    tempL_all = find(idx_L == 1);
    tempR_all = find(idx_R == 1);
    idx_L_mod_all = [];
    idx_R_mod_all = [];
    for k = 1:numel(tempL_all)
        %if this is a reasonably strong input in baseline on either hemisphere for either cell type
        if nanmean(data(tempL_all(k), drug == 0 & celltype == 1)) >= thresh || nanmean(data(tempR_all(k), drug == 0 & celltype == 1)) >= thresh || ...
                nanmean(data(tempL_all(k), drug == 0 & celltype == 2)) >= thresh || nanmean(data(tempR_all(k), drug == 0 & celltype == 2)) >= thresh
            idx_L_mod_all = [idx_L_mod_all tempL_all(k)];       %keep to plot
            idx_R_mod_all = [idx_R_mod_all tempR_all(k)];       %keep to plot
        end
    end
    num_region_LR = numel(idx_L_mod_all);  % now get a region list ignoring hemispheres
    regionList_LR = []; %names of these regions extracted directly from spreadsheet
    for j = 1:num_region_LR
        regionList_LR{j} = region_name{idx_L_mod_all(j)}(1:end-2);  %take out the '-L' from left hemisphere list
    end
    regionList_LR = regionList_LR';

    % calculate fractional difference, drug-evoked (sum across hemispheres)
    data_LR = data(idx_L_mod_all, :) + data(idx_R_mod_all, :);  % use sum
%     data_LR = (data(idx_L_mod_all, :) + data(idx_R_mod_all, :)) ./ 2;  % use avg
    dat_sal = data_LR(:, drug == 0 & celltype == kk);
    dat_psi = data_LR(:, drug == 1 & celltype == kk);
    fracdiff = 100 * (nanmean(dat_psi,2) - nanmean(dat_sal,2)) ./ nanmean(dat_sal,2);

    % arrange rabies input regions based on CFOS region list/organization
    fracdiff_select = nan(size(fracdiff_cfos, 1), 1);
    for i = 1:size(fracdiff_select, 1)

        % current region
        if sum(strcmp(regionList_LR, names_subreg(i))) > 0
            cr = find(strcmp(regionList_LR, names_subreg(i)));
            fracdiff_select(i, :) = fracdiff(cr, :);
        else    % could just skip since init matrix with nan, but leaving for now...
            fracdiff_select(i, :) = NaN; 
        end

    end

    %%%%%%%%
    % plot %
    %%%%%%%%
    
    subplot(2, 1, kk);

    % scatter
    xd = fracdiff_cfos;
    yd = nanmean(fracdiff_select, 2);
    splot = scatter(xd, yd, 225, cmap_sorted, 'filled', 'MarkerEdgeColor', 'k');

    % get & set dot labels
    for getset_dot_labels = 1

        tpos_x = xd + (.04 * range(xd));
        tpos_y = yd + (.05 * range(yd));

        % apply labels
        for i = 1:size(xd) %full list is 1-40
            DesiredColors = cmap_sorted(i,:);
            text(tpos_x(i), tpos_y(i), names_subreg(i), 'color', DesiredColors, 'FontSize', 12);
        end

    end
    clear getset_dot_labels tpos_* DesiredColors

    % axis labels and titles
    xlabel(['Drug-evoked difference, c-Fos (%)']);
    ylabel(['Drug-evoked difference, rabies (%)']);

    % set x
    xx = [-100, ((round(max(xd) / 100) * 100) + 100)];
    xticks(xx(1):100:xx(2)); 
    xlim(xx);
    
    % set y
    yy = [-80, 80];
    yticks([yy(1):40:yy(2)]);
    ylim(yy);

    if kk==1
        title('PT (Fezf2-CreER)');
    elseif kk==2
        title('IT (PlexinD1-CreER)');
    end

    % Stats
    [r, p] = corrcoef(xd, yd, 'rows', 'complete');
    if kk == 1
        cn = 'PT';
    else
        cn = 'IT';
    end
    disp(['Drug-evoked input difference vs. drug-evoked c-Fos difference, ' cn ' neurons, r = ' num2str(r(1,2)) ', P-value = ' num2str(p(1,2))]);

    % print stats to plot
    if p >= .001
        text((xx(1) + 0.03 * range(xx)), (yy(1) + 0.05 * range(yy)), ['{\itr} = ' sprintf('%.3f', r(1,2)) ', {\itP} = ' sprintf('%.3f', p(1,2))], 'FontSize', 16);
    else
        text((xx(1) + 0.03 * range(xx)), (yy(1) + 0.05 * range(yy)), ['{\itr} = ' sprintf('%.3f', r(1,2)) ', {\itP} < .001'], 'FontSize', 16);
    end

    %%%%%%%%%%%%
    % clean up %
    %%%%%%%%%%%%

    clear i j k
    clear thresh tempL_all tempR_all idx_L_mod_all idx_R_mod_all
    clear num_region_LR regionList_LR 
    clear data data_LR fracdiff fracdiff_select
    clear splot xd yd tpos_x tpos_
    clear r p cn pnum xx yy

end
clear kk 

%%%%%%%%%%%%%%
% save plots %
%%%%%%%%%%%%%%
saveas(gcf,fullfile(figdir,['fig_input-difference-by-cfos-difference-scatterplot.png']))
savefig(gcf,fullfile(figdir,['fig_input-difference-by-cfos-difference-scatterplot.fig']))

% clean up
clear cfostab cmaps_sorted fracdiff_cfos
clear cfos_tid cfos_bid cfos_t cfos_b
clear dat_psi dat_sal
clear names_subreg names_subreg_fulcher names_genes








    