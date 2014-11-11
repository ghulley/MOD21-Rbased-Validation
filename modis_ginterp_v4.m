function [gtemp,gorig] = modis_ginterp_v4(g,gp,tr4,scale,c66,cextra,r2)

% Testing inverse distance algorithm for MOD21

% For interpolation:
% Values of 1 = pixels to interpolate
% Values of 0 = graybody pixels (used for interpolation)
% Values of nan = don't interpolate pixel (e.g. cloud or bad data)
%
% figure()
% imagesc(g{1})
% caxis([0 3])
%
% figure()
% imagesc(tr4)

c66 = logical(c66);
cextra = logical(cextra);

ming = 0;
maxg1 = 2;
maxg2 = 3;

% Set negative g values to 1 in all bands for stability
tn1 = g{1}<ming | g{2}<ming | g{3}<ming;

for b = 1:3
    g{b}(tn1) = 1;
end

% Set large g values to 1 in all bands for stability
tn2 = g{1}>maxg1 | g{2}>maxg1 | g{3}>maxg1;

for b = 1:3
    g{b}(tn2) = 1;
end

% figure()
% imagesc(g{1})
% caxis([0 3])

% % Set g values above 2 to 1 for tlr pixels only
% gtemp1 = g{1}; gtemp2 = g{2}; gtemp3 = g{3};
% gtemp1(~tr4) = 1; gtemp2(~tr4) = 1; gtemp3(~tr4) = 1;
% tn3 = gtemp1>maxg1 | gtemp2>maxg1 | gtemp3>maxg1;
% for b = 1:3
%     g{b}(tn3) = 1;
% end


% gdiff = abs(g{1}-g{2});
% tn3 = gdiff>0.5;
% for b = 1:3
%     g{b}(tn3) = 1;
% end

gn = isnan(g{1}) | isnan(g{2}) | isnan(g{3});

for b = 1:3
    
    % Set any nans to 1
    g{b}(gn) = 1;
    
    % Set non-gray, clear pixels to 1 using gray pixel thresholder, gp
    g{b}(gp) = 1;
    
    % Set extra cloudy pixels using c95 confidence to 1
    %g{b}(cextra) = 1;
    
    % Set 66% cloudy values to nan
    g{b}(c66) = nan;
    
    % Set ocean pixels to nan
    %g{b}(oceanpix) = nan;
    
    
end

%
% figure()
% imagesc(g{1})
% caxis([0 3])


% Check if any gray pixels on scene
cgp = find(g{1}>1);
if isempty(cgp) || length(cgp)<50
    gs = size(g{1});
    gtemp{1} = ones(gs(1),gs(2));
    gorig{1} = ones(gs(1),gs(2));
    
    gtemp{2} = ones(gs(1),gs(2));
    gorig{2} = ones(gs(1),gs(2));
    
    gtemp{3} = ones(gs(1),gs(2));
    gorig{3} = ones(gs(1),gs(2));
    return
end

% return

%% Resize gammas using nearest neighbor interpolation
gtemp = cell(1,3);
for b = 1:3
    
    % Resize gamma by one quarter
    gtemp{b} = imresize(g{b},1/scale,'nearest');
    
end

gs = size(gtemp{1});

%% Remove pixels near clouds using cloud displacement vector

% Find all graybody pixels
[c1,c2] = find(gtemp{1}~=1);

sg = size(gtemp{1});

% 5x5 pixel displacement vector
d1 = [ 1 0; -1 0; 1 1; 0 1; -1 1; 1 -1; 0 -1; -1 -1]; % displacement vector
d2 = 2.*[ 1 0; -1 0; 1 1; 0 1; -1 1; 1 -1; 0 -1; -1 -1];
d3 = [ 1 2; -1 2; 2 1;-1 2; 2 1;2 1 ;2 -1; 2 -1];
d = [d1; d2; d3];

% Loop over all valid gamma values, skip cloud pixels (nan)
for j = 1:length(c1)
    
    if isnan(gtemp{1}(c1(j),c2(j)))
        continue
    end
    
    loc = [c1(j) c2(j)];
    %loc = [425,96];
    neighbors = d+repmat(loc,[24 1]);
    %neighbors = d+repmat(loc,[8 1]);
    
    % Remove edge rows
    [t1,t2] = find(neighbors(:,1)<=0 | neighbors(:,1)>sg(1));
    
    neighbors(t1,:) = [];
    
    % Remove edge columns
    [t1,t2] = find(neighbors(:,2)<=0 | neighbors(:,2)>sg(2));
    
    neighbors(t1,:) = [];
    
    % get neighbor values
    nbor = zeros(1,length(neighbors));
    for kn = 1:length(neighbors);
        nbor(kn) = gtemp{1}(neighbors(kn,1),neighbors(kn,2));
    end
    
    %         % Remove pixels with 2 or less neighbors
    %         if length(find(nbor>0))<=3
    %             for b=1:3
    %                 gtemp{b}(c1(j),c2(j)) = 1;
    %             end
    %         end
    
    % Set gray pixels with clouds in 5x5 vicinity to 1
    cn = isnan(nbor);
    if length(find(cn>0))>=1
        for b=1:3
            gtemp{b}(c1(j),c2(j)) = 1;
        end
    end
    
    
end
%
%     figure()
%     imagesc(gtemp{1})
%     caxis([0 3])


% Check if any gray pixels on scene
cgp = find(gtemp{1}>1);
if isempty(cgp) || length(cgp)<50
    gs = size(g{1});
    gtemp{1} = ones(gs(1),gs(2));
    gorig{1} = ones(gs(1),gs(2));
    
    gtemp{2} = ones(gs(1),gs(2));
    gorig{2} = ones(gs(1),gs(2));
    
    gtemp{3} = ones(gs(1),gs(2));
    gorig{3} = ones(gs(1),gs(2));
    return
end

%% Save copy

gorig = gtemp;


%% Method 1, average pixels in given radius around gray pixel


%r2 = 10;    % Radius length in pixels
e = -2;     % p value determines amount of weight given to nearest pixels

testlen1 = 0;

ntest = find(gtemp{1}==1,1);
while isempty(ntest)==0
    
    vc = cell(1,3);
    for b = 1:3
        
        vc{b} = reshape(gtemp{b},1,gs(1)*gs(2));
        
        % Remove bare pixels
        gz = (vc{b}==1);
        vc{b}(gz) = [];
        
        %     gn = isnan(vc{b});
        %     vc{b}(gn) = [];
        %
        %     gg = (vc{b}==100);
        %     vc{b}(gg) = [];
    end
    
    % Set vcc cells as (Nx1) arrays = transposed vc (1xN) arrays
    vcc = cell(1,3);
    for b = 1:3
        vcc{b} = (vc{b})';
    end
    
    % Trim NaN pixels from vcc cells
    vcc_nnan_select = ~isnan(vcc{1});       % non-NaNs logical selection list
    for b = 1:3
        vcc{b} = vcc{b}(vcc_nnan_select);   % select non-NaN elements of vc{b}
    end
    
    %[xg,yg] = find(gtemp{1}>0 & gtemp{1}~=1);                % gray pixels
    [xg,yg] = find(gtemp{1}~=1 & ~isnan(gtemp{1}));          % trim bare and NaN pixels
    [x,y] = find(gtemp{1}==1);                               % bare pixels
    
    testlen = length(y);
    
    if testlen==testlen1   % Remaining bare pixels are not within r2 of gray pixels
        
        for j = 1:length(y)
            
            D=[]; V1=[];V2=[];V3=[];wV =[];vcc1=[];vcc2=[];vcc3=[];
            D2_select = []; % logical selection list for pixels within r distance
            
            D = sqrt((x(j)-xg).^2 +(y(j)-yg).^2);
            
            D2_select = D < r2;
            D2 = D(D2_select);
            
            if isempty(D2)       % Increase distance if D2 is empty
                r3 = r2;
                while isempty(D2)
                    r3 = r3+(r2*2);
                    D2_select = D < r3;
                    D2 = D(D2_select);
                end
                
                
                wV = D2.^e;
                vcc1 = vcc{1}(D2_select); % set vcc's as (Nx1) arrays
                vcc2 = vcc{2}(D2_select);
                vcc3 = vcc{3}(D2_select);
                
                % Commented out trimming NaN pixels at this step,
                % accounted for in setting vcc's, xg, yg arrays
                % vnan = isnan(vcc1);
                % vcc1(vnan) = []; vcc2(vnan) = []; vcc3(vnan) = []; wV(vnan) = [];
                
                if isempty(vcc1)
                    continue
                end
                
                % Not transposing vcc's, they are (Nx1) vectors already
                % V1 = vcc1'.*(wV); V2 = vcc2'.*(wV); V3 = vcc3'.*(wV);
                V1 = vcc1.*(wV); V2 = vcc2.*(wV); V3 = vcc3.*(wV);
                
                %wV = D2.^e;
                V1 = sum(V1)/sum(wV); V2 = sum(V2)/sum(wV); V3 = sum(V3)/sum(wV);
                
                %V = vc{b}'.*(wV);
                %V = sum(V)/sum(wV);
                
                gtemp{1}(x(j),y(j)) = V1;
                gtemp{2}(x(j),y(j)) = V2;
                gtemp{3}(x(j),y(j)) = V3;
                
            end
        end
        
        break
    end
    
    %tic;
    for j = 1:length(y)
        
        %j;
        
        
        %D=[]; V1=[];V2=[];V3=[];wV =[];vcc1=[];vcc2=[];vcc3=[];
        
        D = sqrt((x(j)-xg).^2 +(y(j)-yg).^2);
        
        D2_select = D < r2;
        D2 = D(D2_select);
        
        if isempty(D2)       % Continue if no values found
            continue
        else
            
            wV = D2.^e;
            vcc1 = vcc{1}(D2_select);
            vcc2 = vcc{2}(D2_select);
            vcc3 = vcc{3}(D2_select);
            
            % Commented out % Remove nans
            % vnan = isnan(vcc1);
            % vcc1(vnan) = []; vcc2(vnan) = []; vcc3(vnan) = []; wV(vnan) = [];
            
            if isempty(vcc1)
                continue
            end
            
            % Not transposing vcc's
            % V1 = vcc1'.*(wV); V2 = vcc2'.*(wV); V3 = vcc3'.*(wV);
            V1 = vcc1.*(wV); V2 = vcc2.*(wV); V3 = vcc3.*(wV);
            
            %wV = D2.^e;
            wVs = sum(wV);
            V1 = sum(V1)/wVs; V2 = sum(V2)/wVs; V3 = sum(V3)/wVs;
            
            %V = vc{b}'.*(wV);
            %V = sum(V)/sum(wV);
            
            gtemp{1}(x(j),y(j)) = V1;
            gtemp{2}(x(j),y(j)) = V2;
            gtemp{3}(x(j),y(j)) = V3;
            
            %             subaxis(10,10,s)
            %             imagesc(gtemp{1})
            %             caxis([0.5 2.5])
            %             s=s+1;
            
            
            %  break
            
        end
        
    end
    
    % t=toc
    
    
    %             figure()
    %             %subaxis(4,4,s)
    %             imagesc(gtemp{1})
    %             caxis([0.5 2.5])
    %             s=s+1;
    
    ntest = find(gtemp{1}==1);
    
    testlen1 = testlen;
    
end


%% Method 2, average all pixels
%
%
%     %r2 = 100;    % Radius length in pixels
%     %r3 = 100;  % 'Ping' radius
%     e = -4;     % p value determines amount of weight given to nearest pixels
%
%
%     vc = cell(1,3);
%     for b = 1:3
%
%         vc{b} = reshape(gtemp{b},1,gs(1)*gs(2));
%
%         % Extract non-gray pixels
%         gz = (vc{b}==1);
%         vc{b}(gz) = [];
%
%         %     gn = isnan(vc{b});
%         %     vc{b}(gn) = [];
%         %
%         %     gg = (vc{b}==100);
%         %     vc{b}(gg) = [];
%     end
%
%     [xg,yg] = find(gtemp{1}~=1);                % gray pixels
%     [x,y] = find(gtemp{1}==1);                  % non-gray pixels
%
%
%     for j = 1:length(y)
%
%         D = sqrt((x(j)-xg).^2 +(y(j)-yg).^2);
%
%
%         for b = 1:3
%
%             wV = D.^e;
%             vcc = vc{b};
%             vnan = isnan(vcc);
%             vcc(vnan) = []; wV(vnan) = [];
%
%             if isempty(vcc)
%                 continue
%             end
%
%             V = vcc'.*(wV);
%             V = sum(V)/sum(wV);
%
%             gtemp{b}(x(j),y(j)) = V;
%         end
%
%     end


%% Smooth

% figure()
% imagesc(g31)
% caxis([0.5 2])
% title('interpolated')
%
% figure()
% imagesc(g29_si.^(1/bmp(1)))
% caxis([0.8 2])
% title('g29 smoothed')
%
% figure()
% imagesc(g31_si.^(1/bmp(2)))
% caxis([0.8 2])
% title('g31 smoothed')
%
% figure()
% imagesc(g32_si.^(1/bmp(3)))
% caxis([0.8 2])
% title('g32 smoothed')













