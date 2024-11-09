%% Settings
hSI.hFastZ.enable();
% Q: what's enableFieldCurveCorr? 
% 
%%
sample_pos_um = [37500, 37500, 50000];
hSI.hMotors.moveSample(sample_pos_um);
% hSI.hMotors.setRelativeZero();
% Can I avoid this step and use the absolute axis position? 
%% Script tile workflow
% 1) Make a scanfield for your tiles
% - If the tiles are identical you can make 1 scanfield but do not feed
% each tile this same object as it is a handle, instead make a copy via
% uniqueDuplicateSf = sf.copy(); - same parameters, but a unique object. 
%
% if tiles are different geometry then make different scanfields.
sf0 = scanimage.mroi.scanfield.fields.RotatedRectangle();
sf0.centerXY = [0 0]; 
sf0.sizeXY = [20 20];
sf0.pixelResolutionXY = [512, 512];

numTiles = 9;
% center_pos = [37500, 37500, 0];
tile_size = [10, 10];
tile_overlap = 0;
tile_spacing = 500; % um
num_z_sec = 250;
z_spacing_um = 1;
z_list_um = (0 : num_z_sec-1) * z_spacing_um;

% Create rest of tiles
clear hTiles;
for i = numTiles : -1 : 1
    % Get the sample point and the corner point of the scan field in the
    % sample coordinate (um)
    tmp_sf = sf0.copy();
    [samplePoint, cornerPoints] = scanimage.components.tileTools.tileGeneratorFcns.makeTilePoints(hSI.hCoordinateSystems, tmp_sf, 0);
    
    % Adust sample point and corner points in XY
    % For you
    % Try one dimension first
    samplePoint(1) = samplePoint(1) + (i - 1) * tile_spacing;
    cornerPoints(:, 1) = cornerPoints(:, 1) + (i - 1) * tile_spacing;
    
    tileParams = {false, samplePoint, cornerPoints, tmp_sf.affine, hSI.hChannels.channelsAvailable, [], tmp_sf.pixelResolutionXY, 1};
    hTile = scanimage.components.tiles.tile(tileParams{:});
%     hTile.zPos = z_list_um;
    
    
    % Add to tile array
    hTiles(i) = hTile;
end
%%
hSI.hTileManager.addScanTile(hTiles);
%% Set stack
hSI.hTileManager.isFastZ = true;
hSI.hStackManager.enable = true;
hSI.hStackManager.stackMode = 'fast';
hSI.hStackManager.stackFastWaveformType = 'step';
hSI.hStackManager.stackZStepSize = z_spacing_um;
hSI.hStackManager.numSlices = num_z_sec;

%% To Scan Tiles
% Start
hSI.hTileManager.initTileScanning(); % This starts automated scanning of ALL scan tiles defined.
%%
% Stop
hSI.hTileManager.scanAbortFlag = true; % Important because abort event happens automatically at end of acquisitions
hSI.hTileManager.abortTileScan(); % Calls hSI.abort eventually, clears up flags and workflow.
%%
% hSI.hTileManager.saveTiles('Scan');