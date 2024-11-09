%% Create ScanField
sf = scanimage.mroi.scanfield.fields.RotatedRectangle();

% Set CenterXY. These are Scan Angles!
sf.centerXY = [xVal yVal]; % Generally make this 0,0 this will center the tile at the desired sample location when you scan.

% Set SizeXY. These are Scan Angles! Bidirectional, i.e. a size of 20 means
% +/- 10 degrees around center. There are no checks on this value so make
% sure you do not create a size that is too big for your scanners to
% actually image.
sf.SizeXY = [xRange yRange];

% Set Resolution XY
sf.pixelResolutionXY = [pxPerLine, linesPerFrame];

% This gives you the scanning geometry and the necessary affine

% You can visualize this in the ROI Editor
hRoiGroup = scanimage.mroi.RoiGroup();
hRoi = scanimage.mroi.Roi;
hRoi.add(zPos, sf)
hRoiGroup.add(hRoi);
hSI.hRoiManager.roiGroupMroi.copyobj(hRoiGroup);

%% Make sure size is valid
hSI.hScan2D.fovCornerPoints; % Corner points of FOV in Scan Angles. Size of SF must fall within these params (bi-direction)

%% Generate Sample Points
% Uses the geometry of the scan and our coordinate system - specifically
% the motor alignment transform - to give you a real world size and center
% position of the soon to be tile from the scan geometry

[samplePoint, cornerPoints] = scanimage.components.tileTools.tileGeneratorFcns.makeTilePoints(hSI.hCoordinateSystems,sf,ZPosInMicrons);

% The Z position is where the Tile will be located in Z and is independent
% of the device used for traversing Z - make sure the Z value is in range
% of device! 

% hCoordinate systems gives this tool access to the transforms needed to
% get the real world coordinates

%% Create Tile Params
% tile(tfInMemory, samplePoint, cornerPoints, affine, channels, imData, resolutionXY, displayAvgFactor
%
% tfInMemory: Keep the image data in Memory(RAM) OR Write to Disc (C:\Users\dklab\AppData\Local\Temp\ScanImage_ImagePyramid) and Read
% only when in view. Uses scanimage.util.ImagePyramid to cache and recall. 
%
% samplePoint: CenterXY of the Tile Object in real world space i.e.
% microns, use
% scanimage.components.tileTools.tileGeneratorFcns.makeTilePoints to
% determine this
%
% cornerPoints: Effectively the Size of the Tile Object in real world
% coordinates. Use above function. 
%
% affine: This is the transform used to create a scanfield to image this
% tile. This is why its important to start with a scanfield or SI image.
%
% channels: Number of channels available for your tile to display data on.
% Should match available channels in the system.
%
% imData: Empty for scan tile, or cell array of image data {[] [] [] []}
%
% resolutionXY: Information about the resolution
%
% displayAvgFactor: Information about what the average factor was for the
% image data. 
tileParams = {false, samplePoint, cornerPoints, sf.affine, hSI.hChannels.channelsAvailable, [], sf.pixelResolutionXY, 1};

%% Create Tile
hTile = scanimage.components.tiles.tile(tileParams{:});

%% Create Tile Alternative - Using ROI Data Object as input (tile from SI Image)
% Images acquired in ScanImage are wrapped in stripeData objects that then
% contain roiData objects. This tool can automaticaly conver that to a
% tile. 
tileParams = roiData;
hCoordinateSystems = hSI.hCoordinateSystems;
extraParams = {hSI.hDisplay.displayRollingAverageFactor, ohSI.hChannels.channelsAvailable};
tiles = scanimage.components.tileTools.tileGeneratorFcns.defaultTileGenerator(hCoordinateSystems, tfInMemory, tileParams, extraParams);

% This is also wrapped in an easier function
hSI.hTileManager.makeTiles(tfInMemory, tileParams);

% Access last stripe acquired 
stripeData = hSI.hDisplay.rollingStripeDataBuffer{zIdx}{stripe}; % stripe == frame and the number will always be the rolling avg factor + 1. First index is the averaged frame. 
roiData = stripData.roiData{roiIdx}; % This is sufficient to creat a tile. 

% Deeper info. 
imageData = roiData.imageData{chan};
hRoi = roiData.hRoi; % This has a scanfield!
sf = hSI.hDisplay.rollingStripeDataBuffer{zIdx}{stripe}.roiData{roiIdx}.hRoi.scanfields;

% ^ You can use that sf to create a tile as well. 

%% Add Tile to system
% If Scan Tile - empty image data, desiginates an area to be scanned
hSI.hTileManager.addScanTile(hTile);
% Add as many as you want. hTile can be an array. 

% If Overview Tile - has image data, used to pin images to overview.
hSI.hTileManager.addOverviewTile(hTile);

%% To Scan Tiles
% Start
hSI.hTileManager.initTileScanning(); % This starts automated scanning of ALL scan tiles defined.

% Stop
hSI.hTileManager.scanAbortFlag = true; % Important because abort event happens automatically at end of acquisitions
hSI.hTileManager.abortTileScan(); % Calls hSI.abort eventually, clears up flags and workflow.

%% Save/Load Tile objects
% To Save
hSI.hTileManager.saveTiles(type);
% To Load
hSI.hTileManager.loadTiles(type);
% Type is either 'Scan' or 'Overview'





