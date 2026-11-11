%% PV MINI PROJECT: ALL FOUR UNITS USING DCT
% Run this script from a saved file. MATLAB + Image Processing Toolbox.
% Reference: pv_all_four_units.m. This is a separate DCT-based adaptation.
% Exactly four single-page result windows, one per unit; no tabs.
% All dataset JPEGs are processed. Whole dates stay in one split.
% DCT was selected for this project from earlier reconstruction experiments;
% this does not establish superiority for every enhancement/restoration task.
% Intensities are not calibrated temperatures. Candidate regions are not
% verified faults. Pseudocolor is an artificial intensity palette.
% Cosine-domain filters use reflecting boundaries; no other transform is used.

clear;
close all;
clc;

%% SETTINGS
projectFolder = fileparts(mfilename('fullpath'));
addpath(projectFolder);
datasetFolder = fullfile(projectFolder, 'PV_images_failures');
outputFolder = fullfile(projectFolder, 'PV_DCT_results');
saveResults = false;
showDetailedOutput = false;
representativeImage = 1;
N = 128;
budgets = [0.05 0.10 0.25];
noiseNames = ["Gaussian","SaltPepper","Speckle","Poisson"];
noiseSeed = 2026;
rng(noiseSeed,'twister');
assert(license('test','image_toolbox'), 'Image Processing Toolbox is required.');

manifest = pv_dataset_manifest(datasetFolder);
imageCount = height(manifest);
developmentCount = nnz(manifest.Split == "development");
evaluationCount = nnz(manifest.Split == "evaluation");
assert(ismember(representativeImage,1:imageCount), 'Invalid representativeImage.');

%% ORTHONORMAL DCT-II: ANALYSIS C = T*A*T', SYNTHESIS A = T'*C*T
[k,j] = ndgrid(0:N-1,0:N-1);
T = sqrt(2/N)*cos(pi*k.*(j+0.5)/N);
T(1,:) = T(1,:)/sqrt(2);
assert(norm(T*T'-eye(N),'fro')/N < 1e-10, 'DCT basis validation failed.');
[u,v] = meshgrid(0:N-1);
D = hypot(u,v);
frequencyResponses = {double(D<=12), 1./(1+(D/12).^4), ...
    exp(-D.^2/(2*12^2))};
% A fixed low-frequency model, rather than largest-magnitude truncation,
% defines the slowly varying background for candidate-region exploration.
backgroundResponse = double(D<=8);

%% RESULT CONTAINERS: ONLY THE REPRESENTATIVE IMAGE IS CACHED
reconstructionRows = cell(imageCount*numel(budgets),5);
restorationRows = cell(imageCount*numel(noiseNames)*numel(budgets),6);
baselineRows = cell(imageCount*numel(noiseNames)*2,5);
observationRows = cell(imageCount*4,4);
colorRows = cell(imageCount*4,4);
candidateRows = cell(imageCount,5);
regionRows = cell(0,9);
validationData = zeros(imageCount,2);
r = 0; s = 0; b = 0; o = 0; c = 0;
experimentTimer = tic;

for imageID = 1:imageCount
    raw = imread(manifest.FullPath(imageID));
    if size(raw,3) == 3
        grayImage = rgb2gray(raw);
    else
        grayImage = raw;
    end
    I = im2double(grayImage);
    A = imresize(I,[N N],'bilinear');
    C = T*A*T';
    full = T'*C*T;
    relativeError = norm(A-full,'fro')/max(norm(A,'fro'),eps);
    energyError = abs(sum(C(:).^2)-sum(A(:).^2))/max(sum(A(:).^2),eps);
    assert(relativeError < 1e-10 && energyError < 1e-10, ...
        'DCT reconstruction/energy validation failed for image %d.',imageID);
    validationData(imageID,:) = [relativeError energyError];

    %% UNIT 2: CLEAN RECONSTRUCTION AT MATCHED COEFFICIENT BUDGETS
    for fraction = budgets
        [B,retained] = reconstruct(C,T,fraction);
        [mseValue,psnrValue] = quality(A,B);
        r = r+1;
        reconstructionRows(r,:) = {imageID,fraction,retained,mseValue,psnrValue};
    end

    %% UNIT 3: SAME NOISY OBSERVATION FOR DCT AND SPATIAL BASELINE
    observations = makeNoise(A);
    for n = 1:numel(noiseNames)
        Y = observations{n};
        noisyC = T*Y*T';
        for fraction = budgets
            [B,retained] = reconstruct(noisyC,T,fraction);
            [mseValue,psnrValue] = quality(A,B);
            s = s+1;
            restorationRows(s,:) = ...
                {imageID,noiseNames(n),fraction,retained,mseValue,psnrValue};
        end
        [baseline,baselineName] = spatialRestore(Y,noiseNames(n));
        comparisons = {Y,baseline};
        methods = ["Observed",baselineName];
        for q = 1:2
            [mseValue,psnrValue] = quality(A,comparisons{q});
            b = b+1;
            baselineRows(b,:) = {imageID,noiseNames(n),methods(q),mseValue,psnrValue};
        end
    end

    % Controlled observation models, with a known reference and known filter.
    [blurred,blurRestored,blurResponse] = blurExperiment(A,T,D);
    [periodic,periodicRestored,notch] = periodicExperiment(A,T,u,v);
    observationViews = {blurred,blurRestored,periodic,periodicRestored};
    observationMethods = ["Cosine blur + noise","Known-model Wiener", ...
        "Injected cosine interference","Known-frequency notch"];
    for q = 1:4
        [mseValue,psnrValue] = quality(A,observationViews{q});
        o = o+1;
        observationRows(o,:) = {imageID,observationMethods(q),mseValue,psnrValue};
    end

    %% UNIT 4: CANDIDATE RESIDUALS, MORPHOLOGY AND COLOR NOISE
    background = T'*(C.*backgroundResponse)*T;
    [rawMask,cleanedMask,residual,threshold] = candidateMask(A,background);
    objects = regionprops(cleanedMask,residual,'Area','Centroid', ...
        'BoundingBox','MeanIntensity');
    candidateRows(imageID,:) = ...
        {imageID,threshold,nnz(rawMask),nnz(cleanedMask),numel(objects)};
    for q = 1:numel(objects)
        xy = objects(q).Centroid;
        originalX = (xy(1)-0.5)*size(I,2)/N+0.5;
        originalY = (xy(2)-0.5)*size(I,1)/N+0.5;
        regionRows(end+1,:) = {imageID,q,objects(q).Area,xy(1),xy(2), ...
            originalX,originalY,objects(q).MeanIntensity,objects(q).Area/(N*N)};
    end
    P = pseudocolor(A);
    [colorViews,colorMethods] = colorExperiment(P);
    for q = 1:4
        [mseValue,psnrValue] = quality(P,colorViews{q});
        c = c+1;
        colorRows(c,:) = {imageID,colorMethods(q),mseValue,psnrValue};
    end

    if imageID == representativeImage
        demo = struct('Raw',raw,'I',I,'A',A,'C',C,'Full',full, ...
            'OriginalType',class(grayImage),'Background',background, ...
            'RawMask',rawMask,'Mask',cleanedMask,'Residual',residual, ...
            'Objects',objects,'P',P,'Blurred',blurred,'BlurRestored',blurRestored, ...
            'BlurResponse',blurResponse,'Periodic',periodic, ...
            'PeriodicRestored',periodicRestored,'Notch',notch);
        demo.Noisy = observations;
        demo.ColorViews = colorViews;
    end
end
batchSeconds = toc(experimentTimer);

%% TRACEABLE TABLES AND DATE-BALANCED VALIDATION
reconstructionResults = metadata(cell2table(reconstructionRows,'VariableNames', ...
    {'ImageID','Fraction','StoredRealScalars','MSE','PSNR'}),manifest);
restorationResults = metadata(cell2table(restorationRows,'VariableNames', ...
    {'ImageID','Noise','Fraction','StoredRealScalars','MSE','PSNR'}),manifest);
baselineResults = metadata(cell2table(baselineRows,'VariableNames', ...
    {'ImageID','Noise','Method','MSE','PSNR'}),manifest);
observationResults = metadata(cell2table(observationRows,'VariableNames', ...
    {'ImageID','Method','MSE','PSNR'}),manifest);
colorResults = metadata(cell2table(colorRows,'VariableNames', ...
    {'ImageID','Method','MSE','PSNR'}),manifest);
candidateResults = metadata(cell2table(candidateRows,'VariableNames', ...
    {'ImageID','Threshold','RawPixels','CleanedPixels','RegionCount'}),manifest);
regionMeasurements = cell2table(regionRows,'VariableNames', ...
    {'ImageID','RegionID','AreaPixels128','CentroidX128','CentroidY128', ...
     'OriginalCentroidX','OriginalCentroidY','MeanResidual','ImageAreaFraction'});
if ~isempty(regionMeasurements)
    regionMeasurements = metadata(regionMeasurements,manifest);
else
    regionMeasurements.Split = strings(0,1);
    regionMeasurements.DateGroup = strings(0,1);
    regionMeasurements.RelativePath = strings(0,1);
end
validationResults = table(manifest.ImageID,validationData(:,1),validationData(:,2), ...
    'VariableNames',{'ImageID','FullRelativeError','EnergyError'});
[reconstructionSummary,reconstructionByDate,worstReconstructionCases] = ...
    pv_validation_summary(reconstructionResults,{'Fraction'});

% Choose a restoration budget using development images only, once per noise.
chosenRows = cell(numel(noiseNames),3);
evaluationDCT = restorationResults([], :);
for n = 1:numel(noiseNames)
    development = restorationResults(restorationResults.Split == "development" & ...
        restorationResults.Noise == noiseNames(n),:);
    S = pv_validation_summary(development,{'Fraction'});
    S = sortrows(S,{'DateMeanPSNR','Fraction'},{'descend','ascend'});
    chosenRows(n,:) = {noiseNames(n),S.Fraction(1),S.DateMeanPSNR(1)};
    keep = restorationResults.Split == "evaluation" & ...
        restorationResults.Noise == noiseNames(n) & ...
        abs(restorationResults.Fraction-S.Fraction(1)) < 1e-12;
    evaluationDCT = [evaluationDCT; restorationResults(keep,:)]; %#ok<AGROW>
end
chosenSettings = cell2table(chosenRows,'VariableNames', ...
    {'Noise','Fraction','DevelopmentDateMeanPSNR'});
evaluationDCT.Method = repmat("DCT",height(evaluationDCT),1);
comparisonColumns = {'ImageID','Noise','Method','MSE','PSNR', ...
    'Split','DateGroup','RelativePath'};
evaluationComparison = [evaluationDCT(:,comparisonColumns); ...
    baselineResults(baselineResults.Split == "evaluation",comparisonColumns)];
[restorationSummary,restorationByDate,worstRestorationCases] = ...
    pv_validation_summary(evaluationComparison,{'Noise','Method'});
observationResults.Experiment = repmat( ...
    ["Cosine blur";"Cosine blur";"Periodic";"Periodic"],imageCount,1);
colorResults.Noise = repmat( ...
    ["Gaussian";"Gaussian";"SaltPepper";"SaltPepper"],imageCount,1);
[observationSummary,observationByDate] = pv_validation_summary( ...
    observationResults(observationResults.Split == "evaluation",:),{'Experiment','Method'});
[colorSummary,colorByDate] = pv_validation_summary( ...
    colorResults(colorResults.Split == "evaluation",:),{'Noise','Method'});
syntheticControls = detectorControls(T,backgroundResponse);

%% CONCISE COMMAND WINDOW RESULTS
terminalRestoration = chosenSettings(:,{'Noise','Fraction'});
terminalRestoration.EvaluationMeanPSNR = zeros(4,1);
terminalRestoration.EvaluationDateMeanPSNR = zeros(4,1);
terminalRestoration.MeanGainOverObserved = zeros(4,1);
for n = 1:numel(noiseNames)
    R = restorationSummary(restorationSummary.Noise == noiseNames(n),:);
    dctRow = R(R.Method == "DCT",:);
    observedRow = R(R.Method == "Observed",:);
    terminalRestoration.EvaluationMeanPSNR(n) = dctRow.MeanPSNR;
    terminalRestoration.EvaluationDateMeanPSNR(n) = dctRow.DateMeanPSNR;
    terminalRestoration.MeanGainOverObserved(n) = dctRow.MeanPSNR-observedRow.MeanPSNR;
end
E = candidateResults(candidateResults.Split == "evaluation",:);
eval10 = reconstructionResults(reconstructionResults.Split == "evaluation" & ...
    abs(reconstructionResults.Fraction-0.10)<1e-12,:);
fprintf('Unit 1: %d images (%d development/%d evaluation); source %d x %d, %s.\n', ...
    imageCount,developmentCount,evaluationCount,size(demo.I,1),size(demo.I,2),demo.OriginalType);
fprintf('Unit 2: evaluation PSNR at 10%% = %.2f dB; inverse error %.2e; energy error %.2e.\n', ...
    mean(eval10.PSNR),max(validationData(:,1)),max(validationData(:,2)));
fprintf('Unit 3: mean date-balanced evaluation PSNR across four noise models = %.2f dB.\n', ...
    mean(terminalRestoration.EvaluationDateMeanPSNR));
fprintf('Unit 4: %d candidate regions in %d evaluation images (unverified intensity deviations).\n', ...
    sum(E.RegionCount),height(E));
if showDetailedOutput
    disp(reconstructionSummary); disp(terminalRestoration); disp(syntheticControls);
    disp(chosenSettings); disp(restorationSummary); disp(restorationByDate);
    disp(observationSummary); disp(colorSummary); disp(regionMeasurements);
    fprintf('Batch calculations: %.1f seconds.\n',batchSeconds);
end

%% WINDOW 1: UNIT 1 — DIGITAL IMAGE FUNDAMENTALS
unitFigure('Unit 1: Digital image fundamentals');
layout = tiledlayout(3,4,'TileSpacing','compact');
tile(layout,demo.Raw,'Original image'); tile(layout,round(demo.I*255)/255,'Grayscale: 8-bit quantization');
nexttile(layout); histogram(demo.I(:),256,'BinLimits',[0 1]);
xlabel('Intensity'); ylabel('Pixels'); title('Histogram');
perception = [0.2*ones(64),0.8*ones(64)];
perception(22:42,22:42) = 0.5; perception(22:42,86:106) = 0.5;
tile(layout,perception,'Both centre squares = 0.5');
tile(layout,imresize(demo.I,0.5,'nearest'),'Half-resolution sampling');
tile(layout,imresize(demo.I,0.25,'nearest'),'Quarter-resolution sampling');
tile(layout,round(demo.I*15)/15,'4-bit quantization');
tile(layout,round(demo.I*3)/3,'2-bit quantization');
tile(layout,imfilter(demo.I,ones(3)/9,'symmetric'),'Linear: averaging');
tile(layout,medfilt2(demo.I,[3 3],'symmetric'),'Nonlinear: median');
N4 = zeros(5); N4([2 4],3) = 1; N4(3,[2 4]) = 1; N4(3,3) = 2;
N8 = zeros(5); N8(2:4,2:4) = 1; N8(3,3) = 2;
nexttile(layout); imagesc(N4,[0 2]); axis image off; title('4-neighbourhood');
nexttile(layout); imagesc(N8,[0 2]); axis image off; title('8-neighbourhood');
colormap(gray(256));
sgtitle('Acquisition -> numerical image -> sampling/quantization -> operations');

unitResults.Unit1.Original = demo.Raw;
unitResults.Unit1.Grayscale = demo.I;
unitResults.Unit1.Sampling = {demo.I,imresize(demo.I,0.5,'nearest'), ...
    imresize(demo.I,0.25,'nearest')};
unitResults.Unit1.Quantization = {round(demo.I*255)/255, ...
    round(demo.I*15)/15,round(demo.I*3)/3};
unitResults.Unit1.Perception = perception;
unitResults.Unit1.Neighbours4 = N4;
unitResults.Unit1.Neighbours8 = N8;

%% WINDOW 2: UNIT 2 — DCT REPRESENTATION AND RECONSTRUCTION
unitFigure('Unit 2: DCT reconstruction');
layout = tiledlayout(2,3,'TileSpacing','compact');
tile(layout,demo.A,'Common 128 x 128 reference');
nexttile(layout); imagesc(log1p(abs(demo.C))); axis image off; colorbar;
title('DCT-II coefficients: log magnitude');
tile(layout,demo.Full,'Full inverse: validated below 1e-10');
unitResults.Unit2.Reference = demo.A;
unitResults.Unit2.Coefficients = demo.C;
unitResults.Unit2.FullReconstruction = demo.Full;
unitResults.Unit2.Summary = reconstructionSummary;
unitResults.Unit2.Approximations = cell(1,numel(budgets));
for budgetIndex = 1:numel(budgets)
    fraction = budgets(budgetIndex);
    B = reconstruct(demo.C,T,fraction);
    unitResults.Unit2.Approximations{budgetIndex} = B;
    [~,p] = quality(demo.A,B);
    tile(layout,B,sprintf('%.0f%% coefficients | %.2f dB',100*fraction,p));
end
sgtitle('One fixed orthonormal basis; coefficient values adapt to each image');

%% UNIT 3 COMPUTATIONS: FULL TOPIC RESULTS RETAINED WITHOUT EXTRA VIEWS
unitResults.Unit3.GrayTransforms = {demo.A,1-demo.A,log1p(9*demo.A)/log(10), ...
    demo.A.^0.5,demo.A.^2,imadjust(demo.A),histeq(demo.A,256),adapthisteq(demo.A), ...
    min(1,demo.A+0.15),max(0,demo.A-0.15),min(1,1.5*demo.A),min(1,demo.A/0.75)};
unitResults.Unit3.GrayTransformNames = {'Original','Negative','Log','Gamma 0.5', ...
    'Gamma 2','Contrast stretch','Histogram equalization','Adaptive equalization', ...
    'Add','Subtract','Multiply','Divide'};
mask1 = demo.A>mean(demo.A(:)); mask2 = demo.A>0.5;
unitResults.Unit3.Spatial = {imfilter(demo.A,ones(3)/9,'symmetric'), ...
    imfilter(demo.A,fspecial('gaussian',7,1.5),'symmetric'), ...
    medfilt2(demo.A,[3 3],'symmetric')};
unitResults.Unit3.SpatialNames = {'Mean','Gaussian','Median'};
unitResults.Unit3.Logic = {demo.A.*mask1,mask1 & mask2,mask1 | mask2,xor(mask1,mask2)};
unitResults.Unit3.LogicNames = {'Arithmetic masking','AND','OR','XOR'};
unitResults.Unit3.FrequencyResponses = frequencyResponses;
unitResults.Unit3.FrequencyOutputs = cell(1,3);
for q = 1:3
    unitResults.Unit3.FrequencyOutputs{q} = T'*(demo.C.*frequencyResponses{q})*T;
end
unitResults.Unit3.Noisy = demo.Noisy;
unitResults.Unit3.DCTRestored = cell(1,4);
unitResults.Unit3.BaselineRestored = cell(1,4);
for n = 1:4
    setting = chosenSettings(chosenSettings.Noise == noiseNames(n),:);
    unitResults.Unit3.DCTRestored{n} = reconstruct(T*demo.Noisy{n}*T',T,setting.Fraction);
    unitResults.Unit3.BaselineRestored{n} = spatialRestore(demo.Noisy{n},noiseNames(n));
end
unitResults.Unit3.BlurObserved = demo.Blurred;
unitResults.Unit3.BlurRestored = demo.BlurRestored;
unitResults.Unit3.PeriodicObserved = demo.Periodic;
unitResults.Unit3.PeriodicRestored = demo.PeriodicRestored;
unitResults.Unit3.RestorationSummary = restorationSummary;

%% WINDOW 3: UNIT 3 — ONE COMBINED RESULT PAGE
unitFigure('Unit 3: Enhancement, filtering and restoration results');
layout = tiledlayout(3,4,'TileSpacing','compact');
tile(layout,demo.A,'Reference intensity');
tile(layout,unitResults.Unit3.GrayTransforms{6},'Contrast enhancement');
tile(layout,unitResults.Unit3.GrayTransforms{7},'Histogram equalization');
tile(layout,unitResults.Unit3.Spatial{1},'Spatial mean filtering');
filterNames = {'Ideal cosine LPF','Butterworth cosine LPF','Gaussian cosine LPF'};
for q = 1:3
    tile(layout,unitResults.Unit3.FrequencyOutputs{q},filterNames{q});
end
tile(layout,demo.Noisy{1},'Gaussian noisy observation');
[~,p] = quality(demo.A,unitResults.Unit3.DCTRestored{1});
tile(layout,unitResults.Unit3.DCTRestored{1},sprintf('DCT restored | %.2f dB',p));
tile(layout,demo.BlurRestored,'Known-blur restoration');
tile(layout,demo.PeriodicRestored,'Periodic-noise reduction');
nexttile(layout);
comparisonPSNR = zeros(4,3);
for n = 1:4
    R = restorationSummary(restorationSummary.Noise == noiseNames(n),:);
    comparisonPSNR(n,1) = R.DateMeanPSNR(R.Method == "Observed");
    comparisonPSNR(n,2) = R.DateMeanPSNR(R.Method == "DCT");
    comparisonPSNR(n,3) = R.DateMeanPSNR(R.Method ~= "Observed" & R.Method ~= "DCT");
end
bar(comparisonPSNR); xticks(1:4); xticklabels(noiseNames); xtickangle(20);
ylabel('Date-balanced evaluation PSNR (dB)');
legend('Observed','DCT','Spatial baseline','Location','best');
title('All evaluation images'); grid on;
sgtitle('Unit 3: enhancement -> filtering -> restoration -> evaluation');

%% UNIT 4 COMPUTATIONS: ALL MORPHOLOGY/COLOR RESULTS REMAIN AVAILABLE
binary = false(N); binary(20:65,20:65) = true; binary(32:43,32:43) = false;
binary(80:100,75:105) = true; binary(110,110) = true;
se = strel('disk',2,0);
hit = bwhitmiss(binary,[0 0 0;0 1 0;0 0 0],[1 1 1;1 0 1;1 1 1]);
unitResults.Unit4.Morphology = {binary,imdilate(binary,se),imerode(binary,se), ...
    imopen(binary,se),imclose(binary,se),hit,imfill(binary,'holes'), ...
    bwareaopen(binary,8,8),bwmorph(binary,'skel',Inf),bwperim(binary,8)};
unitResults.Unit4.MorphologyNames = {'Synthetic controls','Dilation','Erosion', ...
    'Opening','Closing','Hit-or-miss','Hole filling','Small-object removal', ...
    'Skeleton','Boundary'};
unitResults.Unit4.Pseudocolor = demo.P;
unitResults.Unit4.HSV = rgb2hsv(demo.P);
unitResults.Unit4.YCbCr = rgb2ycbcr(demo.P);
[x,y] = meshgrid(linspace(0,1,N));
unitResults.Unit4.SyntheticFullColor = cat(3,x,y,0.25+0.5*(x>0.5));
changed = unitResults.Unit4.HSV; changed(:,:,2) = 0.5*changed(:,:,2);
unitResults.Unit4.SaturationChanged = hsv2rgb(changed);
unitResults.Unit4.Complement = 1-demo.P;
unitResults.Unit4.ColorSmoothed = zeros(size(demo.P));
unitResults.Unit4.ColorSharpened = zeros(size(demo.P));
for q = 1:3
    low = imfilter(demo.P(:,:,q),fspecial('gaussian',5,1),'symmetric');
    unitResults.Unit4.ColorSmoothed(:,:,q) = low;
    unitResults.Unit4.ColorSharpened(:,:,q) = demo.P(:,:,q)+0.7*(demo.P(:,:,q)-low);
end
unitResults.Unit4.ColorNoiseViews = demo.ColorViews;
unitResults.Unit4.ColorNoiseMethods = colorMethods;
unitResults.Unit4.RawMask = demo.RawMask;
unitResults.Unit4.CleanedMask = demo.Mask;
unitResults.Unit4.SignedResidual = demo.Residual;
unitResults.Unit4.RegionMeasurements = regionMeasurements;

%% WINDOW 4: UNIT 4 — ONE COMBINED RESULT PAGE
unitFigure('Unit 4: Morphology and color processing results');
layout = tiledlayout(3,4,'TileSpacing','compact');
tile(layout,demo.Raw,'Original source');
tile(layout,demo.RawMask,'Candidate intensity mask');
tile(layout,imdilate(demo.RawMask,se),'Dilation');
tile(layout,imerode(demo.RawMask,se),'Erosion');
tile(layout,imopen(demo.RawMask,se),'Opening');
tile(layout,demo.Mask,'Closing + small-region removal');
tile(layout,hit,'Hit-or-miss on synthetic controls');
nexttile(layout); imshow(demo.Raw); hold on;
for q = 1:numel(demo.Objects)
    box = demo.Objects(q).BoundingBox;
    box = [(box(1)-0.5)*size(demo.I,2)/N+0.5, ...
        (box(2)-0.5)*size(demo.I,1)/N+0.5, ...
        box(3)*size(demo.I,2)/N,box(4)*size(demo.I,1)/N];
    rectangle('Position',box,'EdgeColor','r','LineWidth',1.2);
end
title('Final candidates: not verified faults');
tile(layout,demo.P,'Intensity pseudocolor');
tile(layout,unitResults.Unit4.SaturationChanged,'HSV saturation transformation');
tile(layout,unitResults.Unit4.ColorSmoothed,'Color smoothing');
tile(layout,unitResults.Unit4.ColorSharpened,'Color sharpening');
sgtitle('Unit 4: intensity mask -> morphology -> candidate overlay; color processing');

%% OPTIONAL EXPORT: DEFAULT DOES NOT WRITE RESULT FILES
if saveResults
    if ~isfolder(outputFolder)
        mkdir(outputFolder);
    end
    tables = {manifest,reconstructionResults,reconstructionSummary,reconstructionByDate, ...
        validationResults,chosenSettings,restorationSummary,restorationByDate, ...
        worstRestorationCases,observationSummary,colorSummary,candidateResults, ...
        regionMeasurements,syntheticControls};
    files = {'manifest','reconstruction','reconstruction_summary','reconstruction_by_date', ...
        'validation','chosen_settings','restoration_summary','restoration_by_date', ...
        'worst_restoration_cases','observation_summary','color_summary','candidates', ...
        'candidate_regions','synthetic_controls'};
    for q = 1:numel(tables)
        writetable(tables{q},fullfile(outputFolder,[files{q} '.csv']));
    end
    save(fullfile(outputFolder,'results.mat'),'manifest','reconstructionResults', ...
        'restorationResults','baselineResults','observationResults','colorResults', ...
        'chosenSettings','validationResults','noiseSeed','N','budgets');
end

%% LOCAL FUNCTIONS
function [B,retained] = reconstruct(C,T,fraction)
% Largest-magnitude coefficients; storage counts exclude indices/encoding.
    retained = min(numel(C),max(1,round(fraction*numel(C))));
    [~,order] = sort(abs(C(:)),'descend');
    kept = zeros(size(C)); kept(order(1:retained)) = C(order(1:retained));
    B = T'*kept*T;
end

function R = metadata(R,manifest)
    R.Split = manifest.Split(R.ImageID);
    R.DateGroup = manifest.DateGroup(R.ImageID);
    R.RelativePath = manifest.RelativePath(R.ImageID);
end

function [mseValue,p] = quality(A,B)
% Unclipped arrays, reference peak 1. Display clipping never affects metrics.
    mseValue = mean((A(:)-B(:)).^2);
    if mseValue == 0
        p = Inf;
    else
        p = 10*log10(1/mseValue);
    end
end

function observations = makeNoise(A)
    gaussian = A+0.02*randn(size(A));
    impulses = imnoise(A,'salt & pepper',0.03);
    speckle = A.*(1+0.2*randn(size(A)));
    poisson = imnoise(A*1e-10,'poisson')*1e10; % Same convention as reference script.
    observations = {gaussian,impulses,speckle,poisson};
end

function [B,name] = spatialRestore(Y,noise)
    switch noise
        case "SaltPepper"
            B = medfilt2(Y,[3 3],'symmetric'); name = "Median";
        case "Gaussian"
            B = wiener2(Y,[5 5],0.02^2); name = "Local Wiener";
        case "Speckle"
            B = wiener2(Y,[5 5]); name = "Local Wiener";
        case "Poisson"
            B = wiener2(Y,[3 3]); name = "Local Wiener";
    end
end

function [Y,B,H] = blurExperiment(A,T,D)
% Define a synthetic reflecting-boundary blur by its known cosine response.
% This is a controlled observation model, not an estimated camera PSF.
    H = exp(-D.^2/(2*12^2));
    Y = T'*((T*A*T').*H)*T+0.005*randn(size(A));
    K = 0.005^2/max(var(A(:)),eps); % Known-reference synthetic calibration.
    B = T'*((T*Y*T').*H./(H.^2+K))*T;
end

function [Y,B,notch] = periodicExperiment(A,T,u,v)
% Known cosine modes make a controlled DCT-only periodic-noise experiment.
    N = size(A,1);
    [x,y] = meshgrid(0:N-1);
    Y = A+0.08*cos(pi*18*(x+0.5)/N)+0.06*cos(pi*26*(y+0.5)/N);
    notch = ones(N);
    notch(u==18 & v==0) = 0; notch(u==0 & v==26) = 0;
    B = T'*((T*Y*T').*notch)*T;
end

function [raw,cleaned,residual,threshold] = candidateMask(A,background)
    residual = A-background;
    valid = true(size(A));
    valid([1:5 end-4:end],:) = false;
    valid(:,[1:5 end-4:end]) = false;
    values = residual(valid);
    threshold = max(1e-8,mean(values)+3*std(values));
    raw = residual>threshold & valid;
    se = strel('disk',1,0);
    cleaned = bwareaopen(imclose(imopen(raw,se),se),4,8);
    cleaned = cleaned & valid;
end

function controls = detectorControls(T,H)
% Grid response is reported openly; no fault-detection accuracy is implied.
    N = size(T,1); [x,y] = meshgrid(0:N-1);
    uniform = 0.3*ones(N); bright = uniform;
    spot = (floor(N/2)-3):(floor(N/2)+4);
    bright(spot,spot) = 0.8;
    grid = uniform+0.12*double(mod(x,16)<2 | mod(y,16)<2);
    images = {uniform,bright,grid};
    pixels = zeros(3,1); regions = pixels; overlap = NaN(3,1);
    for q = 1:3
        A = images{q}; background = T'*((T*A*T').*H)*T;
        [~,mask] = candidateMask(A,background);
        pixels(q) = nnz(mask); cc = bwconncomp(mask,8); regions(q) = cc.NumObjects;
        if q == 2
            overlap(q) = nnz(mask(spot,spot))/numel(spot)^2;
        end
    end
    controls = table(["Uniform";"Bright spot";"Cell grid"],pixels,regions,overlap, ...
        'VariableNames',{'Control','CandidatePixels','CandidateRegions','KnownSpotCoverage'});
end

function P = pseudocolor(A)
    P = ind2rgb(uint8(round(255*min(1,max(0,A)))),turbo(256));
end

function [views,methods] = colorExperiment(P)
    gaussian = P+0.03*randn(size(P));
    impulse = P; z = rand(size(P)); impulse(z<0.015) = 0; impulse(z>0.985) = 1;
    restoredGaussian = zeros(size(P)); restoredImpulse = restoredGaussian;
    for q = 1:3
        restoredGaussian(:,:,q) = wiener2(gaussian(:,:,q),[5 5],0.03^2);
        restoredImpulse(:,:,q) = medfilt2(impulse(:,:,q),[3 3],'symmetric');
    end
    views = {gaussian,restoredGaussian,impulse,restoredImpulse};
    methods = ["Gaussian RGB observed","Channel Wiener", ...
        "Impulse RGB observed","Channel median"];
end

function f = unitFigure(name)
    f = figure('Name',name,'NumberTitle','off','WindowStyle','normal', ...
        'Position',[40 60 1400 850]);
end

function tile(layout,A,label)
    nexttile(layout);
    if ndims(A) == 3
        imshow(min(1,max(0,A)));
    else
        imshow(A,[0 1]);
    end
    title(label,'Interpreter','none');
end
