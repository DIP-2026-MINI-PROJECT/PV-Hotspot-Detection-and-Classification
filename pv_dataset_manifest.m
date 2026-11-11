function manifest = pv_dataset_manifest(datasetFolder)
%PV_DATASET_MANIFEST Build a sorted manifest for PV dataset JPEG images.
% Expected structure:
% datasetFolder/<date-folder>/Images/<subfolder>/*.jpg

if ~isfolder(datasetFolder)
    error('Dataset folder not found: %s', datasetFolder);
end

files = dir(fullfile(datasetFolder, '**', '*.jpg'));
files = files(~[files.isdir]);

% Keep images located under a folder named "Images".
keep = false(numel(files),1);
for k = 1:numel(files)
    folderParts = split(string(strrep(files(k).folder, '\', '/')), '/');
    keep(k) = any(strcmpi(folderParts, "Images"));
end
files = files(keep);

if isempty(files)
    error(['No .jpg images found under an Images folder in: %s\n' ...
           'Check the dataset folder and its extracted contents.'], datasetFolder);
end

fullPaths = strings(numel(files),1);
relativePaths = strings(numel(files),1);
dateGroups = strings(numel(files),1);

root = char(datasetFolder);
root = regexprep(root, [regexptranslate('escape', filesep) '+$'], '');
prefix = [root filesep];

for k = 1:numel(files)
    fullPaths(k) = string(fullfile(files(k).folder, files(k).name));

    folder = char(files(k).folder);
    if startsWith(folder, prefix)
        relFolder = folder(numel(prefix)+1:end);
    else
        relFolder = folder;
    end
    relativePaths(k) = string(fullfile(relFolder, files(k).name));

    % First directory below datasetFolder identifies the date.
    parts = split(string(strrep(relFolder, '\', '/')), '/');
    dateGroups(k) = parts(1);
end

% Stable ordering for reproducible ImageID values.
[relativePaths, order] = sort(relativePaths);
fullPaths = fullPaths(order);
dateGroups = dateGroups(order);

developmentDates = ["2018_April_27", "2018_April_28", ...
                    "2018_May_4", "2018_December_20"];
evaluationDates  = ["2018_December_21", ...
                    "2019_January_16", "2019_January_19"];

split = strings(numel(dateGroups),1);
for k = 1:numel(dateGroups)
    if any(dateGroups(k) == developmentDates)
        split(k) = "development";
    elseif any(dateGroups(k) == evaluationDates)
        split(k) = "evaluation";
    else
        error('Unassigned date folder "%s". Check the date split in pv_dataset_manifest.m.', ...
              dateGroups(k));
    end
end

manifest = table((1:numel(files))', fullPaths, relativePaths, ...
                 dateGroups, split, ...
    'VariableNames', {'ImageID','FullPath','RelativePath','DateGroup','Split'});
end
