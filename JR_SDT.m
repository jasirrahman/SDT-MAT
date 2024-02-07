%%TD-SDT

% Load data
project_dir = '/Users/jasirrahman/Desktop/McGinley/D'' + C + Survivial Function Update';
file_path = fullfile(project_dir, 'data', 'OCT_behavior.csv');
df = readtable(file_path);

% Fix noise durations
df = df(~isinf(df.noise_dur), :); % Remove rows where noise_dur is infinity
df.noise_dur(df.outcome == 3) = 14; % Set noise_dur to 14 where outcome is 3

% Create groups based on 'subject_id' and 'reward'
[groups, subject_id, reward] = findgroups(df.subject_id, df.reward);

% Initialize a cell array or table to store results
groupResults = cell(max(groups), 1);  % Adjust based on how you want to store results

% Loop over each group
for i = 1:max(groups)
    % Subset the data for the current group
    groupData = df(groups == i, :);

    % Apply the tr_sdt function to the subset
    groupResult = tr_sdtv2(groupData);

    % Store the result along with the group identifiers
    groupResults{i} = [table(subject_id(i), reward(i), 'VariableNames', {'subject_id', 'reward'}), groupResult];
end

% Optionally, combine the results into a single table
combinedResults = vertcat(groupResults{:});
