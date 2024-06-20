%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Part 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%    Parse the Data   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script created: May 2023 (E. Judd)
% Last updated: June 2024 (E. Judd)

% Script to parse the (1) Permo/Triassic, (2) Triassic/Jurassic, and 
% (3) PETM sea surface temperature (SST) data from the PhanSST database.

% In order to run this script, you'll need to:
% + to download the PhanSST data from: 
%       Judd, E.J., Tierney, J.E., Huber, B.T. et al. The PhanSST global 
%           database of Phanerozoic sea surface temperature proxy data. 
%           Sci Data 9, 753 (2022). 
%           https://doi.org/10.1038/s41597-022-01826-0
%   + add the subfolder "GlobalFilesAndFunctions" to your filepath

%% PART 1: Load and pre-treat PhanSST data

% Step 1: Load data files 
%   (you will need to change you directory informaton based on where the
%   PhanSST file is stored)
iCloud = '/Users/emilyjudd/Library/Mobile Documents/com~apple~CloudDocs';
DataDir = [iCloud,'/PhanSST_RawData/OutputFiles'];
datafilename = [DataDir,'/PhanSST_v001_28Sep2023.csv'];
Data = loadphansst(datafilename);
load("GTS2020_PETM.mat")

% Step 2: Pre-treat data
% Round coordinates
Data.ModLat = round(Data.ModLat);
Data.ModLon = round(Data.ModLon);
% Add PETM
Preferences.PETM = "deepMIP";
Data.ProxyType(Data.ProxyType == "d18c" & Data.Taxon1 == "pf") = "d18cforam";
Data = addpetm(Data,Preferences);
% Corect standization of d18O phosphate data
Corrections.NBS120c = 21.7;
Corrections.Durango = 9.8;
idx = find(strcmpi(Data.ProxyType, 'd18p'));
Data.ProxyValue(idx) = Data.ProxyValue(idx) - ...
        (Data.NBS120c(idx) - Corrections.NBS120c);
idx = find(strcmpi(Data.ProxyType, 'd18p') & ...
    strcmpi(Data.AnalyticalTechnique, 'SIMS'));
Data.DiagenesisFlag(idx(isnan(Data.Durango(idx)))) = 1;
Data.ProxyValue(idx) = Data.ProxyValue(idx) - ...
        (Data.Durango(idx) - Corrections.Durango);
Data(Data.Period=="Jurassic"&Data.ProxyType=="d18p",:) = [];
% Remove diagenetically altered data
Data(Data.DiagenesisFlag == 1,:) = [];
Data(Data.ProxyType == 'tex' & ...
        (Data.dRI>0.5 | Data.MI>0.5 | Data.BIT>0.5),:) = [];
% Remove estuarine/freshwater data
Data(contains(lower(Data.Environment),["estuarine","freshwater"]),:) = [];
% Remove data only assigned to the stage level of the Ypresian
% (Because the Ypresian is so long, we're going to subdivide it into early,
% middle and late)
Data(Data.AgeFlag == 1 & Data.Age == GTS.Average(GTS.Stage == "Ypresian"),:) = [];


%% PART 2: Parse the Permo/Triassic data

% Step 1: Parse the P/T data
PT = Data(Data.Stage == "Changhsingian" | Data.Stage == "Induan" | ...
            Data.Stage == "Olenekian" | Data.Stage == "Anisian",:);

% Step 2: Find unique cooridates
coords = unique([PT.ModLat,PT.ModLon,PT.ProxyType],'rows');

% Step 3: Find sites that have the same proxy type & data from at least 2 
%         stages
% (Because no single site spans the full onset and recovery, here we have 
% opted to include sites who have data from at least 2 of the 4 possible 
% stages we'll be looking at.) 
PTsites = table(strings(size(coords,1),1), strings(size(coords,1),1), ...
    NaN(size(coords,1),1), NaN(size(coords,1),1), ...
    NaN(size(coords,1),1), NaN(size(coords,1),1), ...
    'VariableNames', ["SiteName","ProxyType","ModLat","ModLon","Nper","Ntri"]);
s=0;
keepidx = [];
for ii = 1:size(coords,1)
    idx = find(string(PT.ModLat) == coords(ii,1) & string(PT.ModLon) == coords(ii,2) & ...
        PT.ProxyType == coords(ii,3));
    s = s+numel(idx);
    if numel(unique(PT.Stage(idx))) > 1
        PTsites.SiteName(ii) = strjoin(unique(PT.SiteName(idx)),'/');
        PTsites.ProxyType(ii) = coords(ii,3);
        PTsites.ModLat(ii) = str2double(coords(ii,1));
        PTsites.ModLon(ii) = str2double(coords(ii,2));
        PTsites.Nchang(ii) = sum(PT.Stage(idx) == "Changhsingian");
        PTsites.Nind(ii) = sum(PT.Stage(idx) == "Induan");
        PTsites.Nolen(ii) = sum(PT.Stage(idx) == "Olenekian");
        PTsites.Nanis(ii) = sum(PT.Stage(idx) == "Anisian");
        keepidx = [keepidx;idx];
    end
end
PTsites(isnan(PTsites.ModLat),:) = [];
PT = PT(keepidx,:);

% Step 5: Save data
save("./MatFiles/PTdata.mat","PT","PTsites")

%% PART 3: Parse the Triassic/Jurassic data

% Step 1: Parse the T/J data
TJ = Data(Data.Stage == "Rhaetian" | Data.Stage == "Hettangian" | ...
            Data.Stage == "Sinemurian" | Data.Stage == "Pliensbachian",:);

% Step 2: Find unique cooridates
coords = unique([TJ.ModLat,TJ.ModLon,TJ.ProxyType,TJ.Taxon1,TJ.Taxon2],'rows');

% Step 3: Find sites that have the same proxy type & data from at least 2 
%         stages
TJsites = table(strings(size(coords,1),1), strings(size(coords,1),1), ...
    strings(size(coords,1),1), strings(size(coords,1),1), ...
    NaN(size(coords,1),1), NaN(size(coords,1),1), ...
    NaN(size(coords,1),1), NaN(size(coords,1),1), ...
    'VariableNames', ["SiteName","ProxyType","Taxon1","Taxon2","ModLat","ModLon","N1","N2"]);
s=0;
keepidx = [];
for ii = 1:size(coords,1)
    idx = find(string(TJ.ModLat) == coords(ii,1) & string(TJ.ModLon) == coords(ii,2) & ...
        TJ.ProxyType == coords(ii,3) & TJ.Taxon1 == coords(ii,4) & ...
        TJ.Taxon2 == coords(ii,5));
    s = s+numel(idx);
    if numel(unique(TJ.Stage(idx))) > 1
        TJsites.SiteName(ii) = strjoin(unique(TJ.SiteName(idx)),'/');
        TJsites.ProxyType(ii) = coords(ii,3);
        TJsites.Taxon1(ii) = coords(ii,4);
        TJsites.Taxon2(ii) = coords(ii,5);
        TJsites.ModLat(ii) = str2double(coords(ii,1));
        TJsites.ModLon(ii) = str2double(coords(ii,2));
        TJsites.N1(ii) = sum(TJ.Stage(idx) == "Rhaetian");
        TJsites.N2(ii) = sum(TJ.Stage(idx) == "Hettangian");
        TJsites.N3(ii) = sum(TJ.Stage(idx) == "Sinemurian");
        TJsites.N4(ii) = sum(TJ.Stage(idx) == "Pliensbachian");
        keepidx = [keepidx;idx];
    end
end
TJsites(isnan(TJsites.ModLat),:) = [];
TJ = TJ(keepidx,:);

save("./MatFiles/TJdata.mat","TJ","TJsites")

%% PART 4: Parse the PETM data

% Step 1: Parse the PETM data
PETM = Data(Data.Stage == "Thanetian" | Data.Stage == "PETM" | ...
            Data.Stage == "Ypresian" | Data.Stage == "Lutetian",:);

% Step 2: Add in stage position delineation for Ypresian data
% (Since we're subdivising the Ypresin into early, middle, and late)
Range = GTS.LowerBoundary(GTS.Stage=="Ypresian")-GTS.UpperBoundary(GTS.Stage=="Ypresian");
PETM.Stage(PETM.Stage == "Ypresian" & PETM.Age > GTS.LowerBoundary(GTS.Stage=="Ypresian") - Range/3) = ...
    "Ypresian early";
PETM.Stage(PETM.Stage == "Ypresian" & PETM.Age < GTS.LowerBoundary(GTS.Stage=="Ypresian") - Range/3 & ...
    PETM.Age > GTS.LowerBoundary(GTS.Stage=="Ypresian") - 2*Range/3) = ...
    "Ypresian middle";
PETM.Stage(PETM.Stage == "Ypresian" & PETM.Age < GTS.LowerBoundary(GTS.Stage=="Ypresian") - 2*Range/3) = ...
    "Ypresian late";

%  Step 3: Find unique cooridates
coords = unique([PETM.ModLat,PETM.ModLon,PETM.ProxyType],'rows');

% Step 4: Find sites that have the same proxy type & data from at least 2 
%         (sub)stages
PETMsites = table(strings(size(coords,1),1), strings(size(coords,1),1), ...
    NaN(size(coords,1),1), NaN(size(coords,1),1), ...
    NaN(size(coords,1),1), NaN(size(coords,1),1), ...
    'VariableNames', ["SiteName","ProxyType","ModLat","ModLon","N1","N2"]);
s=0;
keepidx = [];
for ii = 1:size(coords,1)
    idx = find(string(PETM.ModLat) == coords(ii,1) & string(PETM.ModLon) == coords(ii,2) & ...
        PETM.ProxyType == coords(ii,3));
    s = s+numel(idx);
    if numel(unique(PETM.Stage(idx))) > 1
        PETMsites.SiteName(ii) = strjoin(unique(PETM.SiteName(idx)),'/');
        PETMsites.ProxyType(ii) = coords(ii,3);
        PETMsites.ModLat(ii) = str2double(coords(ii,1));
        PETMsites.ModLon(ii) = str2double(coords(ii,2));
        PETMsites.N1(ii) = sum(PETM.Stage(idx) == "Thanetian");
        PETMsites.N2(ii) = sum(PETM.Stage(idx) == "PETM");
        PETMsites.N3(ii) = sum(PETM.Stage(idx) == "Ypresian early");
        PETMsites.N4(ii) = sum(PETM.Stage(idx) == "Ypresian middle");
        PETMsites.N5(ii) = sum(PETM.Stage(idx) == "Ypresian late");
        PETMsites.N6(ii) = sum(PETM.Stage(idx) == "Lutetian");
        keepidx = [keepidx;idx];
    end
end
PETMsites(isnan(PETMsites.ModLat),:) = [];
PETM = PETM(keepidx,:);

% Step 5: Save data
save("./MatFiles/PETMdata.mat","PETM","PETMsites")
