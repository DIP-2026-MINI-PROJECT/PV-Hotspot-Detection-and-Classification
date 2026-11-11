function [summaryTable, byDateTable, worstCases] = pv_validation_summary(inputTable, groupVariables)
%PV_VALIDATION_SUMMARY Summarize PSNR by method/settings and collection date.
% Required columns: PSNR, DateGroup, ImageID, RelativePath, plus groupVariables.

if ~istable(inputTable)
    error('First input must be a table.');
end
groupVariables = cellstr(groupVariables);
required = [groupVariables, {'PSNR','DateGroup','ImageID','RelativePath'}];
missing = setdiff(required, inputTable.Properties.VariableNames);
if ~isempty(missing)
    error('Missing required table variable(s): %s', strjoin(missing, ', '));
end

% Overall summaries for each requested group.
[G, groupKeys] = findgroups(inputTable(:, groupVariables));
nGroups = max(G);
summaryRows = cell(nGroups, numel(groupVariables) + 6);

for g = 1:nGroups
    idx = (G == g);
    values = inputTable.PSNR(idx);
    dates = unique(inputTable.DateGroup(idx));
    dateMeans = zeros(numel(dates),1);
    for d = 1:numel(dates)
        dateMeans(d) = mean(values(inputTable.DateGroup(idx) == dates(d)), 'omitnan');
    end

    summaryRows(g,1:numel(groupVariables)) = table2cell(groupKeys(g,:));
    summaryRows(g,numel(groupVariables)+1:end) = ...
        {nnz(idx), mean(values,'omitnan'), std(values,'omitnan'), ...
         min(values,[],'omitnan'), max(values,[],'omitnan'), mean(dateMeans,'omitnan')};
end

summaryNames = [groupVariables, ...
    {'ImageCount','MeanPSNR','StdPSNR','MinPSNR','MaxPSNR','DateMeanPSNR'}];
summaryTable = cell2table(summaryRows, 'VariableNames', summaryNames);

% Per-date means for each requested group.
[GD, dateKeys] = findgroups(inputTable(:, [groupVariables, {'DateGroup'}]));
nDateGroups = max(GD);
dateRows = cell(nDateGroups, numel(groupVariables)+2);

for g = 1:nDateGroups
    idx = (GD == g);
    dateRows(g,1:numel(groupVariables)+1) = table2cell(dateKeys(g,:));
    dateRows(g,numel(groupVariables)+2) = {mean(inputTable.PSNR(idx),'omitnan')};
end

byDateNames = [groupVariables, {'DateGroup','DateMeanPSNR'}];
byDateTable = cell2table(dateRows, 'VariableNames', byDateNames);

% Lowest-PSNR image in each requested group.
worstRows = cell(nGroups, numel(groupVariables)+3);
for g = 1:nGroups
    idx = find(G == g);
    [~, localIndex] = min(inputTable.PSNR(idx));
    rowIndex = idx(localIndex);
    worstRows(g,1:numel(groupVariables)) = table2cell(groupKeys(g,:));
    worstRows(g,numel(groupVariables)+1:end) = ...
        {inputTable.PSNR(rowIndex), inputTable.ImageID(rowIndex), inputTable.RelativePath(rowIndex)};
end
worstNames = [groupVariables, {'WorstPSNR','ImageID','RelativePath'}];
worstCases = cell2table(worstRows, 'VariableNames', worstNames);
end
