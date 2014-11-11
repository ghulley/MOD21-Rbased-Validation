function returnvals = MODTES_v4th_tlr(varargin)

% MODIS WVS driver

% Process subsection over ASTER scenes (proc=1) or entire granule (proc=2)
proc = 1;

% Set satellite view angle limit (avoid emissivity lambertian effects)
valim = varargin{10};

% Surface type
%surface = varargin{11};

% Use MOD13A2 NDVI product (ndecide=0)
% Use own reflectances from L1B data (ndecide=1);
ndecide=0;


%% File locations

% loc = 'C:\Users\cghughes\Desktop\ESDR-Work\modtes\TES4MODIS_MODAPS';
% year = ' ';
% month = ' ';

modsens = varargin{1};

% Vegetation Index location
MOD13A2loc = varargin{2};

% Snow Index location
MOD10A2loc = varargin{7};

% site location
pointlat = varargin{8};
pointlon = varargin{9};

%% Use Tonooka (T) or Hulley (H) coefficients
coef = 'T';

% Version
%vwvs = 'v1';
vwvs = 'v2'; % Includes view angles till 65 degrees and satellite altitude fix (750 to 705 km)

wvsloc = 'C:\Users\ghulley\Documents\MATLAB\WVS_Coefficients\';
if strcmp(modsens,'MOD')
    wvslocs = [wvsloc,'MODIS_Terra\',vwvs,'\'];
elseif strcmp(modsens,'MYD')
    wvslocs = [wvsloc,'MODIS_Aqua\',vwvs,'\'];
end

% Open WVS coefficients
fid = fopen([wvslocs,'WVS_',modsens,'.txt'],'r');
frm = repmat('%g ',1,12);
c_wvs = fscanf(fid,frm,[12,3]);
fclose(fid);

% Open Band Model Paramater coefficients
fid = fopen([wvslocs,'BMP_',modsens,'.txt'],'r');
frm = repmat('%g ',1,3);
bmp = fscanf(fid,frm);
fclose(fid);

% Open Sky Irradiance coefficients
fid = fopen([wvslocs,'SKY_',modsens,'.txt'],'r');
frm = repmat('%g ',1,3);
c_sky = fscanf(fid,frm,[3,3]);
fclose(fid);

% MODTRAN module to use (MODTRAN5, MODTRAN5_2, MODTRAN5_3)
mmod = 'C:\Users\ghulley\Documents\MATLAB\MODTRAN5';

% MODIS Reprojection module
mrtloc = 'C:\cygwin\ModisTools\Modis\bin\';

% Response function location
srfloc = 'C:\Users\ghulley\Desktop\MODIS_SRF_GEN\output\';

% home location
cdhome = 'C:\Users\ghulley\Documents\MATLAB\work\mcode\';


for f = 1:1 %length(files_02);
    
    f
    
    [tok, rem]=strtok(varargin{3},'\');
    while 1
        if isempty(rem)
            break;
        else
            [tok, rem]=strtok(rem,'\');
        end
    end
    MOD02_name = tok;
    %
    %     [tok, rem]=strtok(varargin{6},'\');
    %     while 1
    %         if isempty(rem)
    %             break;
    %         else
    %             [tok, rem]=strtok(rem,'\');
    %         end
    %     end
    %     MOD07_name = tok;
    %
    %
    %     [tok, rem]=strtok(varargin{5},'\');
    %     while 1
    %         if isempty(rem)
    %             break;
    %         else
    %             [tok, rem]=strtok(rem,'\');
    %         end
    %     end
    %     MOD35_name = tok;
    %
    %     [tok, rem]=strtok(varargin{4},'\');
    %     while 1
    %         if isempty(rem)
    %             break;
    %         else
    %             [tok, rem]=strtok(rem,'\');
    %         end
    %     end
    %     MOD03_name = tok;
    
    obs = MOD02_name(10:26);
    
    file_MOD021km = varargin{3};
    file_MOD07 = varargin{6};
    file_MOD35 = varargin{5};
    file_MOD03 = varargin{4};
    %file_MOD05 = varargin{11};
    
    testopen(file_MOD03)
    
    SD_id = hdfsd( 'start' ,file_MOD03, 'read' );
    if SD_id==-1
        disp('Bad MOD03 file')
        status = hdfsd('end',SD_id);
        
        returnvals=[1 0 0 0 0 0 0 0 0 0 0 0 -999 -999 0 0];
        returnvals=returnvals'; %makes python happy - or rather, the translator between here and python
        
        continue
    end
    status = hdfsd('end',SD_id);
    
    try
        Lat = double(hdfread(file_MOD03,'Latitude'));
        Lon = double(hdfread(file_MOD03,'Longitude'));
    catch ME
        disp('MOD03 file cant read')
        pause(60)
        try
            Lat = double(hdfread(file_MOD03,'Latitude'));
            Lon = double(hdfread(file_MOD03,'Longitude'));
        catch ME
            
            returnvals=[0 1 0 0 0 0 0 0 0 0 0 0 -999 -999 0 0];
            returnvals=returnvals'; %makes python happy - or rather, the translator between here and python
            
            continue
        end
    end
    
    
    %% Destripe MODIS thermal bands (29, 31, 32)
    
    %Yd = MODIS_destripe_v2(file_MOD021km,obs);
    
    disp('Destriping data')
    
    testopen(file_MOD021km)
    
    
    SD_id = hdfsd( 'start' ,file_MOD021km, 'read' );
    if SD_id==-1
        disp('Bad MOD021KM file')
        status = hdfsd('end',SD_id);
        
        returnvals=[0 0 1 0 0 0 0 0 0 0 0 0 -999 -999 0 0];
        returnvals=returnvals'; %makes python happy - or rather, the translator between here and python
        
        continue
    end
    status = hdfsd('end',SD_id);
    
    try
        Rad_TIR = hdfread(file_MOD021km,'EV_1KM_Emissive');
    catch ME
        SD_id = hdfsd( 'start' ,file_MOD021km, 'read' );
        if SD_id==-1
            disp('Bad MOD021KM file')
            status = hdfsd('end',SD_id);
            continue
        end
        status = hdfsd('end',SD_id);
        
        try
            Rad_TIR = hdfread(file_MOD021km,'EV_1KM_Emissive');
        catch ME
            
            
            returnvals=[0 0 0 0 1 0 0 0 0 0 0 0 -999 -999 0 0];
            returnvals=returnvals'; %makes python happy - or rather, the translator between here and python
            
            continue
        end
    end
    
    Rad_TIR = shiftdim(Rad_TIR,1);
    
    disp('Destriping')
    Yd = modis_destripe(Rad_TIR,file_MOD021km,modsens);
    
    % clear RAD_TIR
    
    % Scale radiances
    S = hdfinfo(file_MOD021km);
    
    % Read in radiance scales for EV_1KM_Emissive
    rscales_EV_1KM_Emissive = double(S.Vgroup.Vgroup(2).SDS(3).Attributes(6).Value);
    roffsets_EV_1KM_Emissive = double(S.Vgroup.Vgroup(2).SDS(3).Attributes(7).Value);
    
    % Scale EV_1KM_Emissive values (destriped)
    Y = cell(1,3);
    Y{1} = rscales_EV_1KM_Emissive(9).*(double(Yd{1}) - roffsets_EV_1KM_Emissive(9));
    Y{2} = rscales_EV_1KM_Emissive(11).*(double(Yd{2}) - roffsets_EV_1KM_Emissive(11));
    Y{3} = rscales_EV_1KM_Emissive(12).*(double(Yd{3}) - roffsets_EV_1KM_Emissive(12));
    
    % Read MODIS IR water vapor (MOD07)
    testopen(file_MOD07)
    try
        PWV = double(hdfread(file_MOD07,'Water_Vapor'));
        PWV = PWV*0.001;
    catch ME
        SD_id = hdfsd( 'start' ,file_MOD07, 'read' );
        if SD_id==-1
            disp('Bad MOD07file')
            status = hdfsd('end',SD_id);
            
            returnvals=[0 0 0 0 0 1 0 0 0 0 0 0 -999 -999 0 0];
            returnvals=returnvals'; %makes python happy - or rather, the translator between here and python
            
            continue
        end
        status = hdfsd('end',SD_id);
        
        try
            PWV = double(hdfread(file_MOD07,'Water_Vapor'));
            PWV = PWV*0.001;
        catch ME
            
            returnvals=[0 0 0 0 0 0 1 0 0 0 0 0 -999 -999 0 0];
            returnvals=returnvals'; %makes python happy - or rather, the translator between here and python
            
            continue
        end
    end
    
%      % Read MODIS NIR water vapor (MOD05)
%     testopen(file_MOD05)
%     try
%         PWV_NIR = double(hdfread(file_MOD05,'Water_Vapor_Near_Infrared'));
%         PWV_NIR = PWV_NIR*0.001;
%         
%         PWneg = PWV_NIR<0;
%         PWV_NIR(PWneg) = nan;
%         
%     catch ME
%         SD_id = hdfsd( 'start' ,file_MOD05, 'read' );
%         if SD_id==-1
%             disp('Bad MOD05file')
%             status = hdfsd('end',SD_id);
%             
%               returnvals=[0 0 0 0 0 1 0 0 0 0 0 0 -999 -999];
%             returnvals=returnvals'; %makes python happy - or rather, the translator between here and python
%             
%             continue
%         end
%         status = hdfsd('end',SD_id);
%         
%        
%         try
%             PWV_NIR = double(hdfread(file_MOD05,'Water_Vapor_Near_Infrared'));
%             PWV_NIR = PWV_NIR*0.001;
%             
%             PWneg = PWV_NIR<0;
%             PWV_NIR(PWneg) = nan;
%             
%         catch ME
%             
%               returnvals=[0 0 0 0 0 1 0 0 0 0 0 0 -999 -999];
%             returnvals=returnvals'; %makes python happy - or rather, the translator between here and python
%             
%             disp('Could not read MOD05 PWV')
%             continue
%         end
%     end
    
    
    % Open Cloud mask
    testopen(file_MOD35)
    disp('Opened Cloudmask')
    try
        CM1 = double(hdfread(file_MOD35,'Cloud_Mask'));
        CM1 = squeeze(CM1(1,:,:));
    catch ME
        SD_id = hdfsd( 'start' ,file_MOD07, 'read' );
        if SD_id==-1
            disp('Bad MOD07file')
            status = hdfsd('end',SD_id);
            
            returnvals=[0 0 0 0 0 0 0 0 0 0 0 0 -999 -999 0 0];
            returnvals=returnvals'; %makes python happy - or rather, the translator between here and python
            
            continue
        end
        status = hdfsd('end',SD_id);
        
        try
            CM1 = double(hdfread(file_MOD35,'Cloud_Mask'));
            CM1 = squeeze(CM1(1,:,:));
        catch ME
            returnvals=[0 0 0 0 0 0 0 0 0 0 0 0 -999 -999 0 0];
            returnvals=returnvals'; %makes python happy - or rather, the translator between here and python
            
            continue
        end
    end
    
    
    % Eliminate geolocation fill values (-999)
    Lfill = (Lat==-999);
    Lat(Lfill) = nan;
    Lon(Lfill) = nan;
    
    lf = find(Lfill==1,1);
    
    % Set CM1 fill values to Lat,Lon and Y fill
    CMbad = (CM1==0);
    
    cs = size(CMbad);
    ls = size(Lat);
    
    
    if cs(1)>ls(1) && isempty(lf)
        disp('more fill values in cloud mask than lat/lon, radiances')
        
        for i = 1:cs(1)
            cs = find(CM1(i,:)==0);
            if length(cs)>3
                CM1(i,:) = nan;
            end
        end
        
        %CM1(CMbad) = nan;
        cf = find(isnan(CM1(:,1)));
        
        % Remove fill values from extra scan lines
        while ~isempty(cf)
            
            CM1(cf(1),:) = [];
            
            cf = find(isnan(CM1(:,1)));
        end
        
    else
        
        % Remove extra scan lines from L1B MODIS products (eg. 2040 lines instead of 2030)
        for i = 1:ls(1)
            cs = find(Lfill(i,:)==1);
            if length(cs)>3
                Lat(i,:) = nan;
            end
        end
        latf = find(isnan(Lat(:,1)));
        
        if ~isempty(latf)    % If there are nans
            
            % Remove nan values
            while ~isempty(latf)
                Lat(latf(1),:) = [];
                Lon(latf(1),:) = [];
                
                Y{1}(latf(1),:) = [];
                Y{2}(latf(1),:) = [];
                Y{3}(latf(1),:) = [];
                
                CM1(latf(1),:) = [];
                PWV_NIR(latf(1),:) = [];
                
                latf = find(isnan(Lat(:,1)));
            end
            
        end
        
    end
    
    % If sizes still not same, skip granule
    szt1 = size(CM1);
    szt2 = size(Y{1});
    
    if szt1(1)~=szt2(1)
        
        returnvals=[szt1(1) szt2(1) 0 0 0 0 0 0 0 0 0 0 -999 -999 0 0];
        returnvals=returnvals'; %makes python happy - or rather, the translator between here and python
        
        continue
    end
    
    %% Define granule sizes
    
    % Specify Loop indices for full area
    if proc==2
        % Full area
        sz1 = size(Y{1});
        sz5 = size(PWV);
        ar1_1 = 1:sz1(1); ar1_2 = 1:sz1(2);
        ar5_1 = 1:sz5(1); ar5_2 = 1:sz5(2);
        
    elseif proc==1   % check to see if MODIS granule covers area of interest
        
        %         ap = 100;
        %
        %         szmod = size(Lat);
        %         [dist,az] = distance(Lat,Lon,pointlat,pointlon);
        %         [C,IND] = min(dist(:));
        %         [I,J] = ind2sub(szmod,IND);
        %
        %         %         St = double(hdfread(file_MOD07,'Sensor_Zenith'));
        %         %         St = St.*0.01;
        %
        %         top = I-ap; bot = I+ap;
        %         if top<1
        %             top = 1;
        %         end
        %
        %         if bot>2030
        %             bot = 2030;
        %         end
        %
        %         % Specified area
        %         ar1_1 = top:bot;
        %         ar1_2 = J-ap:J+ap;
        %         ar5_1 = round(top/5):round(bot/5); ar5_2 = round((J-ap)/5):round((J+ap)/5);
        %         if ar5_1(1) == 0
        %             ar5_1=ar5_1+1;
        %         end
        %         sz1 = [length(ar1_1) length(ar1_2)];    % 1 km
        %         sz5 = [length(ar5_1) length(ar5_2)];    % 5 km
        %
        %         %         St = St(ar5_1,ar5_2);
        %         %
        %         %         if St(round(sz5(1)/2),round(sz5(1)/2))>45   % If center pixel viewing angle>45 deg, then skip
        %         %             disp('Center view angle too large')
        %         %             continue
        %         %         end
        
        ap = 100;
        
        szmod = size(Lat);
        [dist,az] = distance(Lat,Lon,pointlat,pointlon);
        [C,IND] = min(dist(:));
        [I,J] = ind2sub(szmod,IND);
        
        St = double(hdfread(file_MOD07,'Sensor_Zenith'));
        St = St.*0.01;
        
        top = I-ap; bot = I+ap;
        if top<1
            top = 1;
        end
        
        if bot>2030
            bot = 2030;
        end
        
        left = J-ap; right = J+ap;
        if left<1
            left = 1;
        end
        
        if right>1354
            right = 1354;
        end
        
        % Specified area
        ar1_1 = top:bot;
        ar1_2 = left:right;
        ar5_1 = ceil(top/5):floor(bot/5); ar5_2 = ceil(left/5):floor(right/5);
        
        sz1 = [length(ar1_1) length(ar1_2)];    % 1 km
        sz5 = [length(ar5_1) length(ar5_2)];    % 5 km
        
%        Sttest = St(ar5_1,ar5_2);
        
%         if Sttest(round(sz5(1)/2),round(sz5(2)/2))>valim   % If center pixel viewing angle>45 deg, then skip
%             disp('Center view angle too large')
%             myLST=0;myE29=0;myE31=0;myE32=0;myPWV=0;myView=0;th_row=0;th_col=0;
%             return
%         end
        
    end
    
    %% Get Viewing angle
    intype = 'bicubic';
    Senszen = double(hdfread(file_MOD07,'Sensor_Zenith'));
    Senszen = Senszen.*0.01;
    
    Senszen5 = Senszen;
    Senszen5 = Senszen5(ar5_1,ar5_2);
    
    Senszen1 = imresize(Senszen5,[sz1(1) sz1(2)],intype);
    
    % Slice area specified by viewing angle limit, valim
    vin1 = (Senszen1(1,:)>valim);
    vf1 = find(vin1==0);
    Senszen1 = Senszen1(:,vf1(1):vf1(length(vf1)));
    vfd1 = vf1(length(vf1)) - vf1(1)+1;
    
    vin5 = (Senszen5(1,:)>valim);
    vf5 = find(vin5==0);
    Senszen5 = Senszen5(:,vf5(1):vf5(length(vf5)));
    
    %% Get snow/water/ice map
    
    disp('Get Snow/water/ice Map')
    
    %**** For Greenland, assume all gray pixels
    
    % Check to see if snowmap already computed
    year = MOD02_name(11:14);
    fname = ['MOD10info.',MOD02_name(10:22),'.',num2str(pointlat),'_',num2str(pointlon),'.mat'];
    f10repo = ['C:\Users\ghulley\Documents\MATLAB\Uncertainty\Rbased_processing\MOD10_repository\',year,'\',fname];
    if exist(f10repo,'file')
        load(f10repo)
    else
        
        try
            [swi,snowmap] = get_MOD10A2_bitmasked(cdhome,Lat,Lon,obs,mrtloc,MOD10A2loc);
            
        catch ME
            returnvals=[21 0 0 0 0 0 0 0 0 0 0 0 -999 -999 0 0];
            returnvals=returnvals'; %makes python happy - or rather, the translator between here and python
            
            continue
        end
        
        
        % Resize to cloud mask size
        sizetemp = size(CM1);
        swi = imresize(swi,sizetemp,'nearest');
        snowmap = imresize(snowmap,sizetemp,'nearest');
        
        
        snow = snowmap==1;
        ice_lake = swi==100;
        ocean = swi==39;
        inland_water = swi==37;
        water_tot = snow + ice_lake + ocean + inland_water;
        %water_tot = ocean+inland_water;
        
        % If Terra night or Aqua day, flip back to original scene orientation
        if (strcmp(modsens,'MOD') && Lat(1,1)<Lat(2,1)) || (strcmp(modsens,'MYD') && Lat(1,1)<Lat(2,1))
            water_tot = rot90(water_tot,2);
        end
        
        water = water_tot(ar1_1,ar1_2);
        water = water(:,vf1(1):vf1(length(vf1)));
        
        % If scene consists mostly of water (ie. ocean), continue
        szwtest = size(water);
        fwater = find(water==1);
        wper = 100*(length(fwater)/(szwtest(1)*szwtest(2))); 
        
        save(f10repo,'water','wper');
       
    end
    
    if wper>99.95
        returnvals=[0 0 0 0 0 0 0 0 0 0 0 0 -999 -999 0 0];
        returnvals=returnvals'; %makes python happy - or rather, the translator between here and python
        continue
    end

    %cutsize = size(Senszen1);
    %water = ones(cutsize(1),cutsize(2));
    
    
    %% Get NDVI
    
    disp('Get Gray Pixels')
    
    % Check to see if NDVI already computed
    year = MOD02_name(11:14);
    fname = ['MOD13info.',MOD02_name(10:22),'.',num2str(pointlat),'_',num2str(pointlon),'.mat'];
    f13repo = ['C:\Users\ghulley\Documents\MATLAB\Uncertainty\Rbased_processing\MOD13_repository\',year,'\',fname];
    if exist(f13repo,'file')
        load(f13repo)
    else
        
        
        NDVI_mod = get_NDVI(cdhome,ndecide,Lfill,file_MOD021km,MOD02_name,Lat,Lon,obs,mrtloc,MOD13A2loc);
        
        % If Terra night or Aqua day, flip back
        if (strcmp(modsens,'MOD') && Lat(1,1)<Lat(2,1)) || (strcmp(modsens,'MYD') && Lat(1,1)<Lat(2,1))
            NDVI_mod = rot90(NDVI_mod,2);
        end
        
        % Get green vegetation pixels
        NDVI_sub = NDVI_mod(ar1_1,ar1_2);
        
        save(f13repo,'NDVI_sub');
        
    end
    
    %NDVI_sub = zeros(cutsize(1),cutsize(2));

    %% Get cloud mask at 1 km from PWV NIR product
    
    disp('Status: Computing cloud mask..')
    
    cmn = (CM1<0);
    CM1(cmn) = CM1(cmn)+256;
    cc_1 = bitget(CM1,1);     % Cloud mask flag (0-Not determined, 1-determined)
    c_status = cc_1==0;
    
    % 95% confident clear
    cbit2_1 = bitget(CM1,2);
    cloud = (cbit2_1==0);
    c95in = cloud | c_status;   % (** note: cloud=1, non-cloud=0)
    
    % 66% confident clear
    cbit3_1 = bitget(CM1,3);
    cloud_1 = (cbit2_1==0);
    cloud_2 = (cbit3_1==0);
    c66in = cloud_1.*cloud_2 | c_status;
    
    clear cbit3_1 cloud_1 cloud_2 cbit2_1 cloud
    
    % Use 66% confidence level (less conservative)
    
    extend = 1; % Expand clouds by two pixels
    c66 = Cloud_fill_v2_4MODIS(c66in,extend);
    %c95 = Cloud_fill_v2_4MODIS(c95in,2);
    
    % Edge cloud pixels
    cedge = c66-c66in;
    
    se = strel('disk',3);
    c66 = imclose(c66,se);  % Clouds plus edge clouds
    %c95 = imclose(c95,se);  % Clouds plus edge clouds
    
    c95 = c66 | c95in;
    
%     
%     figure()
%     imagesc(cloud)
%     colormap('gray')

    cf = logical(c66in);
    
    % Slice clouds to view angles<40
    c66 = c66(ar1_1,ar1_2);
    c66 = c66(:,vf1(1):vf1(length(vf1)));
    
    c66in = c66in(ar1_1,ar1_2);
    c66in = c66in(:,vf1(1):vf1(length(vf1)));
    
    cedge = cedge(ar1_1,ar1_2);
    cedge = cedge(:,vf1(1):vf1(length(vf1)));
    
    c95 = c95(ar1_1,ar1_2);
    c95 = c95(:,vf1(1):vf1(length(vf1)));
    
    cextra = c95-c66;
    
    % Total cloud percentage on cutout    
    szc = size(c66);
    ctotpercent = 100*(length(find(c66==1))./(szc(1)*szc(2)));
    
    clear cbit2_1 cbit2_3 cloud_1 cloud_2 c66temp1 c66temp2 c1 cloud cmn CM1 cc_1
    
    
    %% Compute surface brightness temperature using EMC/WVD (Tg)
    
    disp('Status: Computing Tg..')
    
    % Replace cloudy radiance values with nans
    for b = 1:3
        Y{b}(cf) = nan;
    end
    
    % Resize PWV to 1 km
    % TIR
    badpwv = PWV<0;
    PWV(badpwv) = nan;
    
    sz1 = size(Y{1});
    PWVin = imresize(PWV,[sz1(1),sz1(2)],'nearest');
    PWVins = PWVin;
    %PWVins = smooth2a(PWVin,25,25);
    
    % Compute surface brightness temperature
    Tg = MODIS_Tg(Y,PWVins,c_wvs,coef);
    
    for b = 1:3
        Tg{b} = Tg{b}(ar1_1,ar1_2);
    end
    
    % Input to MODTRAN
    % TIR
    %PWV1 = smooth2a(PWV,5,5);
    PWV1 = PWV;
    PWV1 = PWV1(ar5_1,ar5_2);
    
    %% Compute transmittance and path radiance for two gammas (1 and 0.7)
    g1 = 1;
    g2 = 0.7;
    
    disp('Atmospheric correction')
    
    atc = 1; % Spatially interpolate profiles (1) or leave in original resolution of 5km (0)
    rfrac = 5;  % (ie. rfrac*5km = 25 km);
    Hscale = 1; % Scale water profiles in MODTRAN with PWV1 estimate
    
    [t1_out,t2_out,path_out,sky_out] = modtran_4modis_v2_Rbased(atc,Hscale,PWV1,file_MOD07,g2,valim,ar5_1,ar5_2,rfrac,mmod,cdhome,srfloc,modsens,c_sky);
    
    
    
    t1 = t1_out; t2 = t2_out; path = path_out; sky = sky_out;
    
    %     ts = size(t1{1});
    %     cloud_refine = imresize(c66,[ts(1) ts(2)],'bicubic');
    %     tlow = cloud_refine>0.2;
    %
    %     for b = 1:3
    %         t1{b}(tlow) = nan;
    %         t2{b}(tlow) = nan;
    %         path{b}(tlow) = nan;
    %         sky{b}(tlow) = nan;
    %     end
    
    valpix = find(t1{1}>0,1);   % Find any valid values left
    if isempty(valpix)
        returnvals=[0 0 0 0 0 0 0 0 0 0 0 0 -999 -999 0 0];
        returnvals=returnvals'; %makes python happy - or rather, the translator between here and python
        
        %system('touch C:\aquadata\finished-modtes.txt');
        
        continue
    end
    
    
    
    % Fill in nan gaps
    r2 = 5;
    e = -2;
    t1i = interp_id(t1,r2,e);
    t2i = interp_id(t2,r2,e);
    skyi = interp_id(sky,r2,e);
    pathi = interp_id(path,r2,e);
    
    % Revert back to full size (5 km)
    t1t = cell(1,3); t2t = cell(1,3);
    patht = cell(1,3); skyt = cell(1,3);
    for b = 1:3
        t1t{b} = imresize(t1i{b},[sz5(1) sz5(2)],'bicubic');
        t2t{b} = imresize(t2i{b},[sz5(1) sz5(2)],'bicubic');
        skyt{b} = imresize(skyi{b},[sz5(1) sz5(2)],'bicubic');
        patht{b} = imresize(pathi{b},[sz5(1) sz5(2)],'bicubic');
    end
    
    % Slice atmos. parameters
    for b = 1:3
        
        t1t{b} = t1t{b}(:,vf5(1):vf5(length(vf5)));
        t2t{b} = t2t{b}(:,vf5(1):vf5(length(vf5)));
        patht{b} = patht{b}(:,vf5(1):vf5(length(vf5)));
        skyt{b} = skyt{b}(:,vf5(1):vf5(length(vf5)));
        
    end
    
    
    %% Start WVS processing
    
    % Slice clouds
    %c66 = c66(ar1_1,ar1_2);
    %c66 = c66(:,vf1(1):vf1(length(vf1)));
    
    %c95 = c95(ar1_1,ar1_2);
    %c95 = c95(:,vf1(1):vf1(length(vf1)));
    
    % Resize atmospheric parameters from 5km to 1km
    intype = 'bicubic';
    t1r = cell(1,3); t2r = cell(1,3); pathr = cell(1,3); skyr = cell(1,3);
    
    for b = 1:3
        %     t1t = imresize(t1ri{b},1/2,intype);
        %     t2t = imresize(t2ri{b},1/2,intype);
        %     patht = imresize(pathri{b},1/2,intype);
        %     skyt = imresize(skyri{b},1/2,intype);
        
        t1r{b} = imresize(t1t{b},[length(ar1_1) vfd1],intype);
        t2r{b} = imresize(t2t{b},[length(ar1_1) vfd1],intype);
        pathr{b} = imresize(patht{b},[length(ar1_1) vfd1],intype);
        skyr{b} = imresize(skyt{b},[length(ar1_1) vfd1],intype);
    end
    
    
    %clear t1t t2t patht skyt
    
    
    % Resize radiances
    for b = 1:3
        Y{b} = Y{b}(ar1_1,ar1_2);
        Y{b} = Y{b}(:,vf1(1):vf1(length(vf1)));
    end
    
    % Compute surface radiance
    surfrad = cell(1,3);
    for b = 1:3
        surfrad{b} = (Y{b} - pathr{b})./t1r{b};
    end
    
    for b = 1:3
        Tg{b} = Tg{b}(:,vf1(1):vf1(length(vf1)));
    end
    
    
    %% Get Gray pixels
    
    ndv = (NDVI_sub>0.3);
    
    sz1=size(ndv);
    gp = ones(sz1(1),sz1(2));  % Initialize gray pixel map
    gp(ndv) = 0;               % Set gray pixels to zero
    
    % Resize to exclude viewing angles larger than valim
    gp = gp(:,vf1(1):vf1(length(vf1)));
    sz2=size(gp);
    
    % Add in water pixels to gray pixel map
    water = logical(water);
    gp(water) = 0;
    
    gp_no_tlr = logical(gp);
    
    % Only use TLRs if there are bare pixels on scene
    gpcheck = find(gp==0);
    if length(gpcheck)~=sz2(1)*sz2(2)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Thermal Log Residuals to refine gray-pixel map
        
        t29 = t1r{1};
        [TLR_ng,TLR_gm] = MODTES_TLR(surfrad,gp_no_tlr,t29);
        
        
        % Correlation coefficient approach
        ccs = zeros(sz2(1),sz2(2));
        for i = 1:sz2(1)
            for j = 1:sz2(2)
                if isnan(TLR_ng{1}(i,j))
                    continue
                end
                TLR_test = [TLR_ng{1}(i,j) TLR_ng{2}(i,j) TLR_ng{3}(i,j)];
                %cmat(:,c) = TLR_test';
                
                cc = corrcoef([TLR_gm' TLR_test']);
                ccs(i,j) = cc(2,1);
                
                %c=c+1;
            end
        end
        
        %cc = corrcoef([TLR_gm' cmat]);
        
        tr4 = ccs>0.9;
        
        % Refine gray pixels by adding TLR gray spectra
        ttot = logical(tr4);
        gp(ttot) = 0;
        
    end
    
    gp = logical(gp);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sgp = size(gp_no_tlr); sgpt = sgp(1)*sgp(2);
    gpe = find(gp_no_tlr==0);
    gpg = find(gp_no_tlr==1);
    
    
    if isempty(gpe) || length(gpe)<10 % If few gray pixels or nearly all gray, use standard correction
        
        
        [emisf_STD,Ts_STD] = TES4MODIS_v2(surfrad,skyr,gp);
        Ts_WVS=0;
        Ts_WVS_no_tlr=0;
        for b = 1:3
            emisf_WVS{b}=0;
            emisf_WVS_no_tlr{b}=0;
        end
        
    else
        
        
        % clear ttot tr4 tr3 tr1 TLR_ng TLR_g TLR Imax Imin
        
        %% Compute WVS Surface radiances from WVS surface brightness temperatures
        
        disp('Apply WVS coefficients')
        
        % Compute Surface radiance, given surface brightness temperature (Tg)
        wave = [8.5288e-6 11.0183e-6 12.0323e-6];
        con1 = 3.741775e-22;
        con2 = 0.0143877;
        B = cell(1,3);
        for b = 1:3
            B{b} = con1./(wave(b)^5*pi.*(exp(con2./(wave(b).*Tg{b}))-1));
        end
        
        %% Compute gamma terms
        
        % Band model parameters
        if strcmp(coef,'T')
            bmp = [1.1869 1.6317 1.6967];
        end
        
        term_1 = cell(1,3);
        term_2t = cell(1,3);
        term_2 = cell(1,3);
        term_3 = cell(1,3);
        for b = 1:3
            
            gf = g2^bmp(b);
            
            term_1{b} = t2r{b}./(t1r{b}.^gf);
            term_2t{b} = (B{b} - pathr{b}./(1 - t1r{b}))./(Y{b} - pathr{b}./(1 - t1r{b}));
            term_3{b} = t2r{b}./t1r{b};
            
        end
        
        
        % Find where term2 is negative in all bands (computing g will lead to nan
        % if negative values not found)
        tn = 0;
        for b = 1:3
            tn = tn + (term_2t{b}<0);
        end
        
        tn = (tn>0);
        
        % Compute gammas for each band
        g = cell(1,3);
        for b = 1:3
            gf = g2^bmp(b);
            
            term_1{b}(tn) = nan;
            
            term_2t{b}(tn) = nan;
            term_2{b} = term_2t{b}.^(g1 - gf);
            
            term_3{b}(tn) = nan;
            
            g{b} = log(term_1{b}.*term_2{b})./log(term_3{b});
            
        end
        
        %clear term_1 term_2t term_2 term_3
        
        r2 = 10;
        scale = 1;
        smooth_scale = 10;
        
        %         if strcmp(surface,'graybody')
        %
        %             % Open geolocation at 1km
        %             Latt = Lat(ar1_1,ar1_2);
        %             Lont = Lon(ar1_1,ar1_2);
        %             Latt = Latt(:,vf1(1):vf1(length(vf1)));
        %             Lont = Lont(:,vf1(1):vf1(length(vf1)));
        %
        %             %find my closest point
        %             th_lat=Latt-pointlat;
        %             th_lat=th_lat.*th_lat;
        %
        %             th_lon=Lont-pointlon;
        %             th_lon=th_lon.*th_lon;
        %
        %             th_dist=th_lat+th_lon;
        %
        %             [sr,sc]=find(th_dist==min(min(th_dist)));
        %
        %             [gtemp_no_tlr gorig_no_tlr] = modis_ginterp_v2_graybody(g,gp_no_tlr,scale,c66,r2,sr,sc);
        %             [gtemp gorig] = modis_ginterp_v2_graybody(g,gp,scale,c66,r2,sr,sc);  % checks for negative gamma values over graybody pixels (for R-based validation)
        %
        %         else
        
        %[gtemp,gorig] = modis_ginterp_v2_rval(g,gp,scale,c66,cextra,r2);  % Fill in all bare pixels as normal
        
        %[gtemp, gorig] = modis_ginterp_v2_nested3_vcc_nnan_select(g, gp, scale, c66, cextra, r2);
        
        [gtemp,gorig] = modis_ginterp_v4(g,gp,tr4,scale,c66,cextra,r2);
         
        %[gtemp_no_tlr,gorig_no_tlr] = modis_ginterp_v2_rval(g,gp_no_tlr,scale,c66,r2);
        %[gtemp_no_tlr,gorig_no_tlr] = modis_ginterp_v2_nested3_vcc_nnan_select(g, gp_no_tlr, scale, c66, cextra, r2);
        
        %end
        
        %[gtemp gorig] = modis_ginterp_v2(g,gp,scale,c66,r2);
        %[gtemp_no_tlr gorig_no_tlr] = modis_ginterp_v2(g,gp_no_tlr,scale,c66,r2);
        
        gsize = size(gp);
        %gsize_no_tlr = size(gp_no_tlr);
        
        % Smooth
        gi = cell(1,3);
        %gi_no_tlr = cell(1,3);
        for b = 1:3
            
            % resize back to 1 km resolution
            gi{b} = imresize(gtemp{b},[gsize(1),gsize(2)],'nearest');
            %gi_no_tlr{b} = imresize(gtemp_no_tlr{b},[gsize_no_tlr(1),gsize_no_tlr(2)],'nearest');
            
            gi{b} = smooth2a(gi{b},smooth_scale,smooth_scale);
            %gi_no_tlr{b} = smooth2a(gi_no_tlr{b},smooth_scale,smooth_scale);
            
            
        end
        
        %% Apply WVS to atmospheric parameters
        
        sc = cell(1,3);
        
        if strcmp(coef,'T')
            % Tonooka
            sc{1} = [-0.017 1.760 -0.04503];
            sc{2} = [-0.003 1.740 -0.06227];
            sc{3} = [-0.005 1.725 -0.06658];
        else
            for b=1:3
                sc{b} = c_sky(:,b);
            end
        end
        
        
        % Compute improved atmospheric effect parameters
        ti = cell(1,3); pathi = cell(1,3); skyi = cell(1,3);
        %ti_no_tlr = cell(1,3); pathi_no_tlr = cell(1,3); skyi_no_tlr = cell(1,3);
        surfradi = cell(1,3);
        %surfradi_no_tlr = cell(1,3);
        
        for b = 1:3
            
            gf = g2^bmp(b);
            ti{b} = t1r{b}.^((gi{b} - gf)./(1 - gf)).*t2r{b}.^((1 - gi{b})./(1 - gf));
            %ti_no_tlr{b} = t1r{b}.^((gi_no_tlr{b} - gf)./(1 - gf)).*t2r{b}.^((1 - gi_no_tlr{b})./(1 - gf));
            
            pathi{b} = pathr{b}.*((1 - ti{b})./(1 - t1r{b}));
            %pathi_no_tlr{b} = pathr{b}.*((1 - ti_no_tlr{b})./(1 - t1r{b}));
            
            path_z = pathi{b}.*((1 - ti{b}.^cos(Senszen1.*pi/180))./(1 - ti{b}));
            %path_z_no_tlr = pathi_no_tlr{b}.*((1 - ti_no_tlr{b}.^cos(Senszen1.*pi/180))./(1 - ti_no_tlr{b}));
            
            skyi_c = sc{b}(1) + sc{b}(2).*path_z + sc{b}(3).*path_z.^2;
            %skyi_c_no_tlr = sc{b}(1) + sc{b}(2).*path_z_no_tlr + sc{b}(3).*path_z_no_tlr.^2;
            
            surfradi_c = (Y{b} - pathi{b})./ti{b};
            %surfradi_c_no_tlr = (Y{b} - pathi_no_tlr{b})./ti_no_tlr{b};
            
            % Remove negative surface radiances
            sneg = surfradi_c<=0;
            %sneg_no_tlr=surfradi_c_no_tlr<=0;
            
            surfradi_c(sneg) = nan;
            %surfradi_c_no_tlr(sneg_no_tlr) = nan;
            
            skyi_c(sneg) = nan;
            %skyi_c_no_tlr(sneg_no_tlr) = nan;
            
            surfradi{b} = surfradi_c;
            %surfradi_no_tlr{b}=surfradi_c_no_tlr;
            
            skyi{b} = skyi_c;
            %skyi_no_tlr{b}=skyi_c_no_tlr;
        end
        
        
        %% TES
        [emisf_WVS,Ts_WVS] = TES4MODIS_v2(surfradi,skyi,gp);
        
        %[emisf_WVS_no_tlr,Ts_WVS_no_tlr] = TES4MODIS_v2(surfradi_no_tlr, skyi_no_tlr,gp_no_tlr);
        
        [emisf_STD,Ts_STD] = TES4MODIS_v2(surfrad,skyr,gp);
        
    end
    
    %% Remove detector striping
    
    % Set clouds to nan
    cloud = emisf_WVS{1}==0;
    
    for b =1:3
        emisf_WVS{b}(cloud) = nan;
        %emisf_WVS_no_tlr{b}(cloud) = nan;
        emisf_STD{b}(cloud) = nan;
    end
    
    Ts_WVS(cloud) = nan;
    Ts_STD(cloud) = nan;
    
    % Correct for bad detector using emissivity data
    %%%commented out as we're doing a subset and don't know where lines are
    %[emisout_WVS Tsout_WVS] = MODIS_baddet(emisf_WVS,Ts_WVS);
    %[emisout_STD Tsout_STD] = MODIS_baddet(emisf_STD,Ts_STD);
    emisout_WVS=emisf_WVS;
    %emisout_WVS_no_tlr=emisf_WVS_no_tlr;
    
    Tsout_WVS=Ts_WVS;
    %Tsout_WVS_no_tlr=Ts_WVS_no_tlr;
    
    emisout_STD=emisf_STD;
    Tsout_STD=Ts_STD;
    
    % Open geolocation at 1km
    Lat = Lat(ar1_1,ar1_2);
    Lon = Lon(ar1_1,ar1_2);
    Lat = Lat(:,vf1(1):vf1(length(vf1)));
    Lon = Lon(:,vf1(1):vf1(length(vf1)));
    
    %find my closest point
    th_lat=Lat-pointlat;
    th_lat=th_lat.*th_lat;
    
    th_lon=Lon-pointlon;
    th_lon=th_lon.*th_lon;
    
    th_dist=th_lat+th_lon;
    
    [th_row,th_col]=find(th_dist==min(min(th_dist)));
    
    %get my LST
    myLST = Tsout_WVS(th_row,th_col);
    %myLST_no_tlr = Tsout_WVS_no_tlr(th_row,th_col);
    myLST_no_tlr = 0;
    myLST_std = Tsout_STD(th_row,th_col);
    
    myE29 = emisout_WVS{1}(th_row,th_col);
    myE31 = emisout_WVS{2}(th_row,th_col);
    myE32 = emisout_WVS{3}(th_row,th_col);
    
    %myE29_no_tlr = emisout_WVS_no_tlr{1}(th_row,th_col);
    %myE31_no_tlr = emisout_WVS_no_tlr{2}(th_row,th_col);
    %myE32_no_tlr = emisout_WVS_no_tlr{3}(th_row,th_col);
    
    myE29_no_tlr = 0;
    myE31_no_tlr = 0;
    myE32_no_tlr = 0;
    
    myE29_std = emisout_STD{1}(th_row,th_col);
    myE31_std = emisout_STD{2}(th_row,th_col);
    myE32_std = emisout_STD{3}(th_row,th_col);
    
    
    %for testing only
    myLat=Lat(th_row, th_col);
    myLon=Lon(th_row, th_col);
    
    % Cloud percent on 10x10 pixel around validation site  
    rowpix = th_row-10:th_row+10;
    colpix = th_col-10:th_col+10;
    sg = size(c66);
    % Remove edge rows
    t1 = rowpix<=0 | rowpix>sg(1);
    rowpix(t1) = [];
    % Remove edge columns
    t2 = colpix<=0 | colpix>sg(2);
    colpix(t2) = [];
    
    c66site = c66(rowpix,colpix);
    
    szc = size(c66site);
    csitepercent = 100*(length(find(c66site==1))./(szc(1)*szc(2)));
    
    
    %        NDVIout = NDVI_mod(ar1_1,ar1_2);
    %        NDVIout = NDVIout(:,vf1(1):vf1(length(vf1)));
    
    %        PWV = PWV(ar5_1,ar5_2);
    %        PWV = PWV(:,vf5(1):vf5(length(vf5)));
    
    %     % Save data
    %     filesave = 'C:\Users\ghulley\Documents\MATLAB\TES4MODIS_MODAPS\MODTES_v3\';
    %     filesavename = ['emis_WVS_.',obs,'.mat'];
    %
    %     if ~isdir(filesave)
    %         mkdir(filesave);
    %     end
    %
    %     file2s = [filesave,filesavename];
    %
    %     try save(file2s,'surfradi','skyi','emisout_WVS','Tsout_WVS','Lat','Lon','gi','t1r','pathr','PWV','Senszen1','vf1','ocean','gp');
    %
    %     catch ME
    %         pause(300)
    %         save(file2s,'surfradi','skyi','emisout_WVS','Tsout_WVS','Lat','Lon','gi','t1r','pathr','PWV','Senszen1','vf1','ocean','gp');
    %     end
    
    returnvals=[myLST_std myLST_no_tlr myLST myE29_std myE31_std myE32_std myE29_no_tlr myE31_no_tlr myE32_no_tlr myE29 myE31 myE32 th_row th_col ctotpercent csitepercent];
    returnvals=returnvals'; %makes python happy - or rather, the translator between here and python
    
    %system('touch C:\aquadata\finished-modtes.txt');
    
end
