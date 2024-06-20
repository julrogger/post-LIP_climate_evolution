%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Part 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%    Parse the Data   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script created: May 2023 (E. Judd)
% Last updated: June 2024 (E. Judd)

% Script to parse the (1) Permo/Triassic, (2) Triassic/Jurassic, and 
% (3) PETM sea surface temperature (SST) data from the PhanSST database.

% To run this script, you'll need to:
%   + run the script Part1_ParseData
%   + run the script Part2_CalculateDeltaT


%% PART 1: Load the data and define the plotting variables
load("GTS2020_PETM.mat")
load('./MatFiles/PTdeltaT.mat')
load('./MatFiles/TJdeltaT.mat')
load('./MatFiles/PETMdeltaT.mat')
load('./MatFiles/PTdata.mat')
load('./MatFiles/TJdata.mat')
load('./MatFiles/PETMdata.mat')
% Define the colormap for the DeltaT figures
P = [5:5:95];
cm = customcolormap(linspace(0,1,2),{'#6A6A6A','#F9F9F9'},floor(numel(P)/2));
% Define the colors for the site map figure
color.d18c = hex2rgb('#0D2847',1);
color.d18p = hex2rgb('#898F54',1);
color.mg = hex2rgb('#BDA437',1);
color.tex = hex2rgb('#8E554A',1);
% Define the marker sizes for the site map figure
ms = [12:-3:6];

%% PART 2: Plot the P/T
figure('Units','centimeters','Position',[5,5,7.3,6],'color','w'); ax = gca; hold on
x1 = GTS.UpperBoundary(GTS.Stage == "Changhsingian");
w1 = GTS.LowerBoundary(GTS.Stage == "Changhsingian")-x1;
x2 = GTS.UpperBoundary(GTS.Stage == "Induan");
w2 = GTS.LowerBoundary(GTS.Stage == "Induan")-x2;
x3 = GTS.UpperBoundary(GTS.Stage == "Olenekian");
w3 = GTS.LowerBoundary(GTS.Stage == "Olenekian")-x3;
x4 = GTS.UpperBoundary(GTS.Stage == "Anisian");
w4 = GTS.LowerBoundary(GTS.Stage == "Anisian")-x4;
DeltaChanghsingian = cell2mat(PTdeltaT.Changhsingian);
DeltaInduan = cell2mat(PTdeltaT.Induan);
DeltaOlenekian = cell2mat(PTdeltaT.Olenekian);
DeltaAnisian = cell2mat(PTdeltaT.Anisian);
y1 = prctile(DeltaChanghsingian,P);
y2 = prctile(DeltaInduan,P);
y3 = prctile(DeltaOlenekian,P);
y4 = prctile(DeltaAnisian,P);
for ii = 1:size(cm,1)
    rectangle('Position',[x1,y1(ii),w1,y1(end-ii+1)-y1(ii)],...
        'FaceColor',cm(ii,:),'EdgeColor','none')
    rectangle('Position',[x2,y2(ii),w2,y2(end-ii+1)-y2(ii)],...
        'FaceColor',cm(ii,:),'EdgeColor','none')
    rectangle('Position',[x3,y3(ii),w3,y3(end-ii+1)-y3(ii)],...
        'FaceColor',cm(ii,:),'EdgeColor','none')
    rectangle('Position',[x4,y4(ii),w4,y4(end-ii+1)-y4(ii)],...
        'FaceColor',cm(ii,:),'EdgeColor','none')
end
plot([x1,x1+w1],[y1(P==50),y1(P==50)],'k-','LineWidth',1.5)
plot([x2,x2+w2],[y2(P==50),y2(P==50)],'k-','LineWidth',1.5)
plot([x3,x3+w3],[y3(P==50),y3(P==50)],'k-','LineWidth',1.5)
plot([x4,x4+w4],[y4(P==50),y4(P==50)],'k-','LineWidth',1.5)
ylim([-9 18])
geologictimescale(x4,x1+w1,...
    'normal','reverse',gca,'standard','stages','on',8,1,'Arial',6,[.5,.5,.5],[.2,.2,.2])
ax.YAxis(1).Color = [.5 .5 .5];
ax.YAxis(2).Color = [.5 .5 .5];
ax.XAxis(1).Color = [.5 .5 .5];
xlabel('Age (Ma)','FontName','Arial','FontSize',8,'color','k')
ylabel(['\DeltaSST (',char(176),'C)'],'FontName','Arial','FontSize',8,'color','k')
export_fig(gcf,'./Figures/PNGs/PTdeltaT.png','-p0.01','-m10')
export_fig(gcf,'./Figures/PDFs/PTdeltaT.pdf','-painters')
copygraphics(gcf,'ContentType','image','BackgroundColor','none') 

%% PART 3: Plot the T/J
figure('Units','centimeters','Position',[5,5,7.3,6],'color','w')
ax = gca; hold on
x1 = GTS.UpperBoundary(GTS.Stage == "Rhaetian");
w1 = GTS.LowerBoundary(GTS.Stage == "Rhaetian")-x1;
x2 = GTS.UpperBoundary(GTS.Stage == "Hettangian");
w2 = GTS.LowerBoundary(GTS.Stage == "Hettangian")-x2;
x3 = GTS.UpperBoundary(GTS.Stage == "Sinemurian");
w3 = GTS.LowerBoundary(GTS.Stage == "Sinemurian")-x3;
x4 = GTS.UpperBoundary(GTS.Stage == "Pliensbachian");
w4 = GTS.LowerBoundary(GTS.Stage == "Pliensbachian")-x4;
DeltaRhaetian = cell2mat(cellfun(@(x) nanmedian(x,2), TJdeltaT.Rhaetian,'UniformOutput', false));
DeltaHettangian = cell2mat(cellfun(@(x) nanmedian(x,2), TJdeltaT.Hettangian,'UniformOutput', false));
DeltaSinemurian = cell2mat(cellfun(@(x) nanmedian(x,2), TJdeltaT.Sinemurian,'UniformOutput', false));
DeltaPliensbachian = cell2mat(cellfun(@(x) nanmedian(x,2), TJdeltaT.Pliensbachian,'UniformOutput', false));
y1 = prctile(DeltaRhaetian,P);
y2 = prctile(DeltaHettangian,P);
y3 = prctile(DeltaSinemurian,P);
y4 = prctile(DeltaPliensbachian,P);
for ii = 1:size(cm,1)
    rectangle('Position',[x1,y1(ii),w1,y1(end-ii+1)-y1(ii)],...
        'FaceColor',cm(ii,:),'EdgeColor','none')
    rectangle('Position',[x2,y2(ii),w2,y2(end-ii+1)-y2(ii)],...
        'FaceColor',cm(ii,:),'EdgeColor','none')
    rectangle('Position',[x3,y3(ii),w3,y3(end-ii+1)-y3(ii)],...
        'FaceColor',cm(ii,:),'EdgeColor','none')
    rectangle('Position',[x4,y4(ii),w4,y4(end-ii+1)-y4(ii)],...
        'FaceColor',cm(ii,:),'EdgeColor','none')
end
plot([x1,x1+w1],[y1(P==50),y1(P==50)],'k-','LineWidth',1.5)
plot([x2,x2+w2],[y2(P==50),y2(P==50)],'k-','LineWidth',1.5)
plot([x3,x3+w3],[y3(P==50),y3(P==50)],'k-','LineWidth',1.5)
plot([x4,x4+w4],[y4(P==50),y4(P==50)],'k-','LineWidth',1.5)
ylim([-9 18])
colormap(cm), ylim([-9 18]),caxis([5 95])
geologictimescale(x4,x1+w1,...
    'normal','reverse',gca,'standard','stages','on',8,1,'Arial',6,[.5,.5,.5],[.2,.2,.2])
ax.YAxis(1).Color = [.5 .5 .5];
ax.YAxis(2).Color = [.5 .5 .5];
ax.XAxis(1).Color = [.5 .5 .5];
xlabel('Age (Ma)','FontName','Arial','FontSize',8,'color','k')
ylabel(['\DeltaSST (',char(176),'C)'],'FontName','Arial','FontSize',8,'color','k')
export_fig(gcf,'./Figures/PNGs/TJdeltaT.png','-p0.01','-m10')
export_fig(gcf,'./Figures/PDFs/TJdeltaT.pdf','-painters')

%% PART 4: Plot the PETM
figure('Units','centimeters','Position',[5,5,7.3,6],'color','w')
ax = gca; hold on
Range = GTS.LowerBoundary(GTS.Stage=="Ypresian")-GTS.UpperBoundary(GTS.Stage=="Ypresian");
x1 = GTS.UpperBoundary(GTS.Stage == "Thanetian");
w1 = GTS.LowerBoundary(GTS.Stage == "Thanetian")-x1;
x2 = GTS.UpperBoundary(GTS.Stage == "PETM");
w2 = GTS.LowerBoundary(GTS.Stage == "PETM")-x2;
x3 = GTS.UpperBoundary(GTS.Stage == "Ypresian")+2*Range/3;
w3 = Range/3;
x4 = GTS.UpperBoundary(GTS.Stage == "Ypresian")+Range/3;
w4 = Range/3;
x5 = GTS.UpperBoundary(GTS.Stage == "Ypresian");
w5 = Range/3;
x6 = GTS.UpperBoundary(GTS.Stage == "Lutetian");
w6 = GTS.LowerBoundary(GTS.Stage == "Lutetian")-x6;
DeltaThanetian = cell2mat(cellfun(@(x) nanmedian(x,2), PETMdeltaT.Thanetian,'UniformOutput', false));
DeltaPETM = cell2mat(cellfun(@(x) nanmedian(x,2), PETMdeltaT.PETM,'UniformOutput', false));
DeltaYearly = cell2mat(cellfun(@(x) nanmedian(x,2), PETMdeltaT.YpresianEarly,'UniformOutput', false));
DeltaYmiddle = cell2mat(cellfun(@(x) nanmedian(x,2), PETMdeltaT.YpresianMiddle,'UniformOutput', false));
DeltaYlate = cell2mat(cellfun(@(x) nanmedian(x,2), PETMdeltaT.YpresianLate,'UniformOutput', false));
DeltaLutetian = cell2mat(cellfun(@(x) nanmedian(x,2), PETMdeltaT.Lutetian,'UniformOutput', false));
y1 = prctile(DeltaThanetian,P);
y2 = prctile(DeltaPETM,P);
y3 = prctile(DeltaYearly,P);
y4 = prctile(DeltaYmiddle,P);
y5 = prctile(DeltaYlate,P);
y6 = prctile(DeltaLutetian,P);
for ii = 1:size(cm,1)
    rectangle('Position',[x1,y1(ii),w1,y1(end-ii+1)-y1(ii)],...
        'FaceColor',cm(ii,:),'EdgeColor','none')
    rectangle('Position',[x2,y2(ii),w2,y2(end-ii+1)-y2(ii)],...
        'FaceColor',cm(ii,:),'EdgeColor','none')
    rectangle('Position',[x3,y3(ii),w3,y3(end-ii+1)-y3(ii)],...
        'FaceColor',cm(ii,:),'EdgeColor','none')
    rectangle('Position',[x4,y4(ii),w4,y4(end-ii+1)-y4(ii)],...
        'FaceColor',cm(ii,:),'EdgeColor','none')
    rectangle('Position',[x5,y5(ii),w5,y5(end-ii+1)-y5(ii)],...
        'FaceColor',cm(ii,:),'EdgeColor','none')
    rectangle('Position',[x6,y6(ii),w6,y6(end-ii+1)-y6(ii)],...
        'FaceColor',cm(ii,:),'EdgeColor','none')
end
plot([x1,x1+w1],[y1(P==50),y1(P==50)],'k-','LineWidth',1.5)
plot([x2,x2+w2],[y2(P==50),y2(P==50)],'k-','LineWidth',1.5)
plot([x3,x3+w3],[y3(P==50),y3(P==50)],'k-','LineWidth',1.5)
plot([x4,x4+w4],[y4(P==50),y4(P==50)],'k-','LineWidth',1.5)
plot([x5,x5+w5],[y5(P==50),y5(P==50)],'k-','LineWidth',1.5)
plot([x6,x6+w6],[y6(P==50),y6(P==50)],'k-','LineWidth',1.5)
colormap(cm), ylim([-9 18]),caxis([5 95])
geologictimescale(x6,x1+w1,...
    'normal','reverse',gca,'standard','stages','on',8,1,'Arial',6,[.5,.5,.5],[.2,.2,.2])
ax.YAxis(1).Color = [.5 .5 .5];
ax.YAxis(2).Color = [.5 .5 .5];
ax.XAxis(1).Color = [.5 .5 .5];
xlabel('Age (Ma)','FontName','Arial','FontSize',8,'color','k')
ylabel(['\DeltaSST (',char(176),'C)'],'FontName','Arial','FontSize',8,'color','k')
export_fig(gcf,'./Figures/PNGs/PETMdeltaT.png','-p0.01','-m10')
export_fig(gcf,'./Figures/PDFs/PETMdeltaT.PDF','-p0.01','-m10')

%% PART 5: Plot the colorbar
figure('Units','centimeters','Position',[5,5,7.3,6],'color','w')
ax = gca; hold on, yyaxis right
for ii = 1:size(cm,1)
    rectangle('Position',[1,P(ii),1,unique(diff(P))],...
        'FaceColor',cm(ii,:),'EdgeColor','none')
    rectangle('Position',[1,P(end-ii),1,unique(diff(P))],...
        'FaceColor',cm(ii,:),'EdgeColor','none')
end
rectangle('Position',[1,0,1,100],...
        'FaceColor','none','EdgeColor',[.5,.5,.5],'LineWidth',1)
plot([1,2],[50,50],'k','LineWidth',1.5)
ax.XAxis.Visible = 'off';
ax.YAxis(1).Visible = 'off';
ax.YAxis(2).Color = [.5 .5 .5];
ax.FontName = 'Arial'; ax.FontSize = 6;
ax.YTick = [0:25:100];
ax.YAxis(2).LineWidth = .25;
ylabel('Percentile','FontName','Arial','FontSize',8,'color',[.2,.2,.2])
ax.Position = [.13,.11,.045,.815];
export_fig(gcf,'./Figures/PNGs/ColorbarVertical.png','-p0.01','-m10')
export_fig(gcf,'./Figures/PDFs/ColorbarVertical.pdf','-painters')

%% PART 6: Location maps
PETMsites.ProxyType(contains(PETMsites.ProxyType,"d18")) = "d18c";
figure('Position',[3 540 1129 282],'Color','w')
% P/T
subplot('Position',[.01,.1,.32,.98])
ax = worldmap('World');mlabel off, plabel off, gridm off
geoshow(shaperead('landareas', 'UseGeoCoords', true),'EdgeColor','none','FaceColor',[.75 .75 .75]);
framem('FEdgeColor',[.75 .75 .75])
idx = unique(PTsites.ProxyType);
for ii = 1:numel(idx)
    plotm(PTsites.ModLat(PTsites.ProxyType == idx(ii)),...
        PTsites.ModLon(PTsites.ProxyType == idx(ii)),'o',...
        'MarkerEdgeColor',color.(idx{ii}),'MarkerFaceColor',color.(idx{ii}),...
        'MarkerSize',ms(ii))
end
title('Permian/Triassic','FontName','Arial','FontSize',15,'HorizontalAlignment','center','Color',[.5 .5 .5])
% T/J
subplot('Position',[.34,.1,.32,.98])
ax = worldmap('World');mlabel off, plabel off, gridm off
geoshow(shaperead('landareas', 'UseGeoCoords', true),'EdgeColor','none','FaceColor',[.75 .75 .75]);
framem('FEdgeColor',[.75 .75 .75])
idx = unique(TJsites.ProxyType);
for ii = 1:numel(idx)
    plotm(TJsites.ModLat(TJsites.ProxyType == idx(ii)),...
        TJsites.ModLon(TJsites.ProxyType == idx(ii)),'o',...
        'MarkerEdgeColor',color.(idx{ii}),'MarkerFaceColor',color.(idx{ii}),...
        'MarkerSize',ms(ii))
end
title('Triassic/Jurassic','FontName','Arial','FontSize',15,'HorizontalAlignment','center','Color',[.5 .5 .5])
% PETM
subplot('Position',[.67,.1,.32,.98])
ax = worldmap('World');mlabel off, plabel off, gridm off
geoshow(shaperead('landareas', 'UseGeoCoords', true),'EdgeColor','none','FaceColor',[.75 .75 .75]);
framem('FEdgeColor',[.75 .75 .75])
idx = unique(PETMsites.ProxyType);
for ii = 1:numel(idx)
    plotm(PETMsites.ModLat(PETMsites.ProxyType == idx(ii)),...
        PETMsites.ModLon(PETMsites.ProxyType == idx(ii)),'o',...
        'MarkerEdgeColor',color.(idx{ii}),'MarkerFaceColor',color.(idx{ii}),...
        'MarkerSize',ms(ii))
end
title('PETM','FontName','Arial','FontSize',15,'HorizontalAlignment','center','Color',[.5 .5 .5])
% Legend
labx = linspace(.25,.75,6);
dim = [labx(2) .2 0 0];
str = "\fontsize{15}δ^{18}O\fontsize{11}carbonate";
annotation('textbox',dim,'String',str,'FontName','Arial','Color',color.d18c,...
    'FontWeight','bold','VerticalAlignment','middle','HorizontalAlignment','center')
dim = [labx(3) .2 0 0];
str = "\fontsize{15}δ^{18}O\fontsize{11}phosphate";
annotation('textbox',dim,'String',str,'FontName','Arial','Color',color.d18p,...
    'FontWeight','bold','VerticalAlignment','middle','HorizontalAlignment','center')
dim = [labx(4) .2 0 0];
str = "\fontsize{15}Mg/Ca";
annotation('textbox',dim,'String',str,'FontName','Arial','Color',color.mg,...
    'FontWeight','bold','VerticalAlignment','middle','HorizontalAlignment','center')
dim = [labx(5) .2 0 0];
str = "\fontsize{15}TEX\fontsize{11}86";
annotation('textbox',dim,'String',str,'FontName','Arial','Color',color.tex,...
    'FontWeight','bold','VerticalAlignment','middle','HorizontalAlignment','center')
export_fig(gcf,'./Figures/PNGs/SupplementalMap.png','-p0.01','-m10')
export_fig(gcf,'./Figures/PDFs/SupplementalMap.pdf','-painters')

