%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Part 2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    Calculate deltaT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script created: May 2023 (E. Judd)
% Last updated: June 2024 (E. Judd)

% Script to calculate the change in sea surface temperature across the 
% (1) Permo/Triassic, (2) Triassic/Jurassic, and (3) PETM 

% To run this script, you'll need to:
%   + run the script Part1_ParseData
%   + install stan (for instructions, see documentation for BAYMAG
%     at https://github.com/jesstierney/BAYMAG)
%   + add the subfolder "GlobalFilesAndFunctions" to your filepath

%% PART 1: The Permo/Triassic

% Step 1: load the data that we saved in the first script
clear
load("./MatFiles/PTdata.mat")

% Step 2: Preallocate the temperature table
PTtemp = table();
PTtemp.SiteName = PTsites.SiteName;
PTtemp.Changhsingian = cell(height(PTtemp),1);
PTtemp.Induan = cell(height(PTtemp),1);
PTtemp.Olenekian = cell(height(PTtemp),1);
PTtemp.Anisian = cell(height(PTtemp),1);

% Step 3: Calculate temperatures for each site
for ii = 1:height(PTsites)
    % Index permian and Triassic data
    sn = strsplit(PTsites.SiteName(ii),'/');
    Cidx = PT.ModLat == PTsites.ModLat(ii) & ...
        PT.ModLon == PTsites.ModLon(ii) & PT.Stage == "Changhsingian";
    Iidx = PT.ModLat == PTsites.ModLat(ii) & ...
        PT.ModLon == PTsites.ModLon(ii) & PT.Stage == "Induan";
    Oidx = PT.ModLat == PTsites.ModLat(ii) & ...
        PT.ModLon == PTsites.ModLon(ii) & PT.Stage == "Olenekian";
    Aidx = PT.ModLat == PTsites.ModLat(ii) & ...
        PT.ModLon == PTsites.ModLon(ii) & PT.Stage == "Anisian";
    % Calculate temperatures
    d18Osw = 0;
    if PTsites.ProxyType(ii) == "d18p"
        eq = "Puceat2010";
        mineralogy = "phosphate";
    elseif PTsites.ProxyType(ii) == "d18c"
        eq = "Kim1997";
        mineralogy = "calcite";
    end
    PTtemp.Changhsingian{ii} = d18O2temp(PT.ProxyValue(Cidx), eq, mineralogy, d18Osw);
    PTtemp.Induan{ii} = d18O2temp(PT.ProxyValue(Iidx), eq, mineralogy, d18Osw);
    PTtemp.Olenekian{ii} = d18O2temp(PT.ProxyValue(Oidx), eq, mineralogy, d18Osw);
    PTtemp.Anisian{ii} = d18O2temp(PT.ProxyValue(Aidx), eq, mineralogy, d18Osw);
end

% Step 4: Determine which stage each site needs to use as a reference frame
% (Not all sites have data over the same window, so we'll need to do some
% adjusting to get everything in the same reference frame.)
RefFrame = NaN(height(PTsites),1);
RefFrame(PTsites.Nchang ~= 0) = 1;
RefFrame(isnan(RefFrame) & PTsites.Nind ~= 0) = 2;
RefFrame(isnan(RefFrame) & PTsites.Nolen ~= 0) = 3;
% --> So 11 sites will have a Changhsingian reference frame
%         3 sites will have an Induan reference frame
%         1 site will have an Olenekian reference frame

% Step 5: Calculate the baseline value for each site, which we'll use to
%         determine the change in temperature
% (For sites with Changsingian data, we'll the median value from that stage
% and site for the baseline)
BaseLine = NaN(height(PTsites),1);
PTdeltaT = PTtemp;
for ii = 1:numel(RefFrame)
    if RefFrame(ii) == 1
        BaseLine(ii) = median(PTtemp.Changhsingian{ii});
        PTdeltaT.Changhsingian{ii} = PTtemp.Changhsingian{ii} - BaseLine(ii);
        PTdeltaT.Induan{ii} = PTtemp.Induan{ii} - BaseLine(ii);
        PTdeltaT.Olenekian{ii} = PTtemp.Olenekian{ii} - BaseLine(ii);
        PTdeltaT.Anisian{ii} = PTtemp.Anisian{ii} - BaseLine(ii);      
    end
end
% For the sites that don't have Changsingian data, we'll adjust their first
% stage's median value to the average for that stage to determine a
% baseline
IndIdx = RefFrame == 1;
dInduan = median(cell2mat(PTdeltaT.Induan(IndIdx)));
for ii = 1:numel(RefFrame)
    if RefFrame(ii) == 2
        BaseLine(ii) = median(PTtemp.Induan{ii}) - dInduan;
        PTdeltaT.Induan{ii} = PTtemp.Induan{ii} - BaseLine(ii);
        PTdeltaT.Olenekian{ii} = PTtemp.Olenekian{ii} - BaseLine(ii);
        PTdeltaT.Anisian{ii} = PTtemp.Anisian{ii} - BaseLine(ii);      
    end
end
OleIdx = RefFrame == 1 | RefFrame == 2;
dOlenekian = median(cell2mat(PTdeltaT.Olenekian(OleIdx)));
for ii = 1:numel(RefFrame)
    if RefFrame(ii) == 3
        BaseLine(ii) = median(PTtemp.Olenekian{ii}) - dOlenekian;
        PTdeltaT.Olenekian{ii} = PTtemp.Olenekian{ii} - BaseLine(ii);
        PTdeltaT.Anisian{ii} = PTtemp.Anisian{ii} - BaseLine(ii);      
    end
end

% Step 6: Save results
save("./MatFiles/PTdeltaT.mat","PTdeltaT")


%% PART 2: The Triassic/Jurassic

% Step 1: load the data that we saved in the first script
clear 
load("./MatFiles/TJdata.mat")

% Step 2: Preallocate the temperature table
TJtemps = table();
TJtemps.SiteName = TJsites.SiteName;
TJtemps.Rhaetian = cell(height(TJtemps),1);
TJtemps.Hettangian = cell(height(TJtemps),1);
TJtemps.Sinemurian = cell(height(TJtemps),1);
TJtemps.Pliensbachian = cell(height(TJtemps),1);

% Step 3: Calculate temperatures for each site
d18Osw = 0;
prior_mean=25;
prior_std=10;

for ii = 1:height(TJsites)
    % Index permian and Triassic data
    Nidx = TJ.ModLat == TJsites.ModLat(ii) & ...
        TJ.ModLon == TJsites.ModLon(ii) & ...
        TJ.ProxyType == TJsites.ProxyType(ii) & ...
        TJ.Taxon1 == TJsites.Taxon1(ii) & TJ.Taxon2 == TJsites.Taxon2(ii);
    N1idx = TJ.Stage(Nidx) == "Rhaetian";
    N2idx = TJ.Stage(Nidx) == "Hettangian";
    N3idx = TJ.Stage(Nidx) == "Sinemurian";
    N4idx = TJ.Stage(Nidx) == "Pliensbachian";
    
    % Calculate temperatures
    if TJsites.ProxyType(ii) == "d18c" & TJsites.Taxon1(ii) == "br"
        eq = "Kim1997";
        mineralogy = "calcite";
        temp = d18O2temp(TJ.ProxyValue(Nidx), eq, mineralogy, d18Osw);
    elseif TJsites.ProxyType(ii) == "d18c" & TJsites.Taxon2(ii) == "bi"
        eq = "Kim1997";
        mineralogy = "calcite";
        temp = d18O2temp(TJ.ProxyValue(Nidx), eq, mineralogy, d18Osw);
    elseif TJsites.ProxyType(ii) == "d18c" & TJsites.Taxon2(ii) == "ce"
        eq = "Daeron2019";
        mineralogy = "calcite";
        temp = d18O2temp(TJ.ProxyValue(Nidx), eq, mineralogy, d18Osw);
    elseif TJsites.ProxyType(ii) == "d18p"
        eq = "Puceat2010";
        mineralogy = "phosphate";
        temp = d18O2temp(TJ.ProxyValue(Nidx), eq, mineralogy, d18Osw);
    elseif TJsites.ProxyType(ii) == "tex"
        search_tol=std(TJ.ProxyValue(Nidx))*2;
        output = bayspar_tex_analog(TJ.ProxyValue(Nidx), ...
            prior_mean, prior_std, search_tol, 'SST', 1000, 1);
        temp = output.PredsEns;
    end
    TJtemps.Rhaetian{ii} = temp(N1idx,:);
    TJtemps.Hettangian{ii} = temp(N2idx,:);
    TJtemps.Sinemurian{ii} = temp(N3idx,:);
    TJtemps.Pliensbachian{ii} = temp(N4idx,:);
    fprintf('Completed site %.0f/%.0f\n',ii,height(TJsites))
end

% Step 4: Determine which stage each site needs to use as a reference frame
RefFrame = NaN(height(TJsites),1);
RefFrame(TJsites.N1 ~= 0) = 1;
RefFrame(isnan(RefFrame) & TJsites.N2 ~= 0) = 2;
RefFrame(isnan(RefFrame) & TJsites.N3 ~= 0) = 3;
RefFrame(isnan(RefFrame) & TJsites.N4 ~= 0) = 4;
% --> So 1 sites will have a Rhaetian reference frame
%        1 site will have an Hettangian reference frame
%        9 site will have an Sinemurian reference frame

% Step 5: Calculate the baseline value for each site, which we'll use to
%         determine the change in temperature
BaseLine = NaN(height(TJsites),1);
TJdeltaT = TJtemps;
for ii = 1:numel(RefFrame)
    if RefFrame(ii) == 1
        BaseLine(ii) = median(TJtemps.Rhaetian{ii},'all');
        TJdeltaT.Rhaetian{ii} = TJtemps.Rhaetian{ii} - BaseLine(ii);
        TJdeltaT.Hettangian{ii} = TJtemps.Hettangian{ii} - BaseLine(ii);
        TJdeltaT.Sinemurian{ii} = TJtemps.Sinemurian{ii} - BaseLine(ii);
        TJdeltaT.Pliensbachian{ii} = TJtemps.Pliensbachian{ii} - BaseLine(ii);
    end
end
% For the sites that don't have Rhaetian data, we'll adjust their first
% stage's median value to the average for that stage to determine a
% baseline
HetIdx = RefFrame == 1;
dHet = nanmedian(cell2mat(cellfun(@(x) nanmedian(x,'all'), ...
    TJdeltaT.Hettangian(HetIdx), 'UniformOutput', false)));
for ii = 1:numel(RefFrame)
    if RefFrame(ii) == 2
        BaseLine(ii) = median(TJtemps.Hettangian{ii},'all') - dHet;
        TJdeltaT.Hettangian{ii} = TJtemps.Hettangian{ii} - BaseLine(ii);
        TJdeltaT.Sinemurian{ii} = TJtemps.Sinemurian{ii} - BaseLine(ii);
        TJdeltaT.Pliensbachian{ii} = TJtemps.Pliensbachian{ii} - BaseLine(ii);
    end
end
SinIdx = RefFrame == 1 | RefFrame == 2;
dSin = nanmedian(cell2mat(cellfun(@(x) nanmedian(x,'all'), ...
    TJdeltaT.Sinemurian(SinIdx), 'UniformOutput', false)));
for ii = 1:numel(RefFrame)
    if RefFrame(ii) == 3
        BaseLine(ii) = median(TJtemps.Sinemurian{ii},'all') - dSin;
        TJdeltaT.Sinemurian{ii} = TJtemps.Sinemurian{ii} - BaseLine(ii);
        TJdeltaT.Pliensbachian{ii} = TJtemps.Pliensbachian{ii} - BaseLine(ii);
    end
end


% Step 6: Save results
save("./MatFiles/TJdeltaT.mat","TJdeltaT")

%% PART 3: The PETM

% Step 1: load the data that we saved in the first script
clear
load("./MatFiles/PETMdata.mat")

% Step 2: Preallocate the temperature table
PETMtemps = table();
PETMtemps.SiteName = PETMsites.SiteName;
PETMtemps.Thanetian = cell(height(PETMtemps),1);
PETMtemps.PETM = cell(height(PETMtemps),1);
PETMtemps.YpresianEarly = cell(height(PETMtemps),1);
PETMtemps.YpresianMiddle = cell(height(PETMtemps),1);
PETMtemps.YpresianLate = cell(height(PETMtemps),1);
PETMtemps.Lutetian = cell(height(PETMtemps),1);

% Step 4: Calculate temperatures for each site
d18Osw = 0;
prior_mean=25;
prior_std=10;
h = 0.75;

for ii = 1:height(PETMsites)
    % Index permian and Triassic data
    Nidx = PETM.ModLat == PETMsites.ModLat(ii) & ...
        PETM.ModLon == PETMsites.ModLon(ii) & ...
        PETM.ProxyType == PETMsites.ProxyType(ii);
    N1idx = PETM.Stage(Nidx) == "Thanetian";
    N2idx = PETM.Stage(Nidx) == "PETM";
    N3idx = PETM.Stage(Nidx) == "Ypresian early";
    N4idx = PETM.Stage(Nidx) == "Ypresian middle";
    N5idx = PETM.Stage(Nidx) == "Ypresian late";
    N6idx = PETM.Stage(Nidx) == "Lutetian";
    
    % Calculate temperatures
    if PETMsites.ProxyType(ii) == "d18a"
        eq = "Grossman1986";
        mineralogy = "aragonite";
        temp = d18O2temp(PETM.ProxyValue(Nidx), eq, mineralogy, d18Osw);
    elseif PETMsites.ProxyType(ii) == "d18cforam" 
        temp = median(predict_seatemp(PETM.ProxyValue(Nidx), d18Osw, prior_mean, ...
            prior_std, false, "none"),2);
    elseif PETMsites.ProxyType(ii) == "tex"
        search_tol=std(PETM.ProxyValue(Nidx))*2;
        output = bayspar_tex_analog(PETM.ProxyValue(Nidx), ...
            prior_mean, prior_std, search_tol, 'SST', 1000, 1);
        temp = output.PredsEns;
    elseif PETMsites.ProxyType(ii) == "mg"
        [~,ph] = omgph(PETMsites.ModLat(ii),PETMsites.ModLon(ii),0);
        output = baymag_predict(PETM.Age(Nidx),PETM.ProxyValue(Nidx),5,35,...
            ph,median(PETM.CleaningMethod(Nidx)),"all",prior_std,2,h);
        temp = median(output.ens,2);
    end
    PETMtemps.Thanetian{ii} = temp(N1idx,:);
    PETMtemps.PETM{ii} = temp(N2idx,:);
    PETMtemps.YpresianEarly{ii} = temp(N3idx,:);
    PETMtemps.YpresianMiddle{ii} = temp(N4idx,:);
    PETMtemps.YpresianLate{ii} = temp(N5idx,:);
    PETMtemps.Lutetian{ii} = temp(N6idx,:);
    fprintf('Completed site %.0f/%.0f\n',ii,height(PETMsites))
end

% Step 4: Determine which stage each site needs to use as a reference frame
RefFrame = NaN(height(PETMsites),1);
RefFrame(PETMsites.N1 ~= 0) = 1;
RefFrame(isnan(RefFrame) & PETMsites.N2 ~= 0) = 2;
RefFrame(isnan(RefFrame) & PETMsites.N3 ~= 0) = 3;
RefFrame(isnan(RefFrame) & PETMsites.N4 ~= 0) = 4;
RefFrame(isnan(RefFrame) & PETMsites.N5 ~= 0) = 5;
RefFrame(isnan(RefFrame) & PETMsites.N5 ~= 0) = 6;
% --> So 19 sites will have a Thanetian reference frame
%         3 sites will have a PETM reference frame
%         9 sites will have an early Ypresian reference frame
%         4 sites will have an middle Ypresian reference frame
%         1 site will have an late Ypresian reference frame


% Step 5: Calculate the baseline value for each site, which we'll use to
%         determine the change in temperature
BaseLine = NaN(height(PETMsites),1);
PETMdeltaT = PETMtemps;
for ii = 1:numel(RefFrame)
    if RefFrame(ii) == 1
        BaseLine(ii) = median(PETMtemps.Thanetian{ii},'all');
        PETMdeltaT.Thanetian{ii} = PETMtemps.Thanetian{ii} - BaseLine(ii);
        PETMdeltaT.PETM{ii} = PETMtemps.PETM{ii} - BaseLine(ii);
        PETMdeltaT.YpresianEarly{ii} = PETMtemps.YpresianEarly{ii} - BaseLine(ii);
        PETMdeltaT.YpresianMiddle{ii} = PETMtemps.YpresianMiddle{ii} - BaseLine(ii);
        PETMdeltaT.YpresianLate{ii} = PETMtemps.YpresianLate{ii} - BaseLine(ii);
        PETMdeltaT.Lutetian{ii} = PETMtemps.Lutetian{ii} - BaseLine(ii);
    end
end
% For the sites that don't have Thanetian data, we'll adjust their first
% stage's median value to the average for that stage to determine a
% baseline
PETMIdx = RefFrame == 1;
dPETM = nanmedian(cell2mat(cellfun(@(x) nanmedian(x,'all'), ...
    PETMdeltaT.PETM(PETMIdx), 'UniformOutput', false)));
for ii = 1:numel(RefFrame)
    if RefFrame(ii) == 2
        BaseLine(ii) = median(PETMtemps.PETM{ii},'all') - dPETM;
        PETMdeltaT.PETM{ii} = PETMtemps.PETM{ii} - BaseLine(ii);
        PETMdeltaT.YpresianEarly{ii} = PETMtemps.YpresianEarly{ii} - BaseLine(ii);
        PETMdeltaT.YpresianMiddle{ii} = PETMtemps.YpresianMiddle{ii} - BaseLine(ii);
        PETMdeltaT.YpresianLate{ii} = PETMtemps.YpresianLate{ii} - BaseLine(ii);
        PETMdeltaT.Lutetian{ii} = PETMtemps.Lutetian{ii} - BaseLine(ii);
    end
end
YeIdx = RefFrame == 1 | RefFrame == 2;
dYe = nanmedian(cell2mat(cellfun(@(x) nanmedian(x,'all'), ...
    PETMdeltaT.YpresianEarly(YeIdx), 'UniformOutput', false)));
for ii = 1:numel(RefFrame)
    if RefFrame(ii) == 3
        BaseLine(ii) = median(PETMtemps.YpresianEarly{ii},'all') - dYe;
        PETMdeltaT.YpresianEarly{ii} = PETMtemps.YpresianEarly{ii} - BaseLine(ii);
        PETMdeltaT.YpresianMiddle{ii} = PETMtemps.YpresianMiddle{ii} - BaseLine(ii);
        PETMdeltaT.YpresianLate{ii} = PETMtemps.YpresianLate{ii} - BaseLine(ii);
        PETMdeltaT.Lutetian{ii} = PETMtemps.Lutetian{ii} - BaseLine(ii);
    end
end
YmIdx = RefFrame == 1 | RefFrame == 2 | RefFrame == 3;
dYm = nanmedian(cell2mat(cellfun(@(x) nanmedian(x,'all'), ...
    PETMdeltaT.YpresianEarly(YmIdx), 'UniformOutput', false)));
for ii = 1:numel(RefFrame)
    if RefFrame(ii) == 4
        BaseLine(ii) = median(PETMtemps.YpresianEarly{ii},'all') - dYm;
        PETMdeltaT.YpresianMiddle{ii} = PETMtemps.YpresianMiddle{ii} - BaseLine(ii);
        PETMdeltaT.YpresianLate{ii} = PETMtemps.YpresianLate{ii} - BaseLine(ii);
        PETMdeltaT.Lutetian{ii} = PETMtemps.Lutetian{ii} - BaseLine(ii);
    end
end
YlIdx = RefFrame == 1 | RefFrame == 2 | RefFrame == 3 | RefFrame == 4;
dYl = nanmedian(cell2mat(cellfun(@(x) nanmedian(x,'all'), ...
    PETMdeltaT.YpresianEarly(YlIdx), 'UniformOutput', false)));
for ii = 1:numel(RefFrame)
    if RefFrame(ii) == 5
        BaseLine(ii) = median(PETMtemps.YpresianEarly{ii},'all') - dYl;
        PETMdeltaT.YpresianLate{ii} = PETMtemps.YpresianLate{ii} - BaseLine(ii);
        PETMdeltaT.Lutetian{ii} = PETMtemps.Lutetian{ii} - BaseLine(ii);
    end
end

% Step 6: Save results
save("./MatFiles/PETMdeltaT.mat","PETMdeltaT")

