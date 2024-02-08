function TD_SDTv2()
    % Main function TD-SDT
    file = path;
    %check general file info
    info = tdmsinfo(file);
    %check channels
    tdmsreadprop(file);
    %load channels
    channelsToLoad = ["Trial Number", "NC Time (chord)", "React Time (chord)", "Trial Result"];
    %read necessary parts of tdms file into MATLAB
    data = tdmsread(file, ...
    ChannelGroupName='Untitled', ...
    ChannelNames = channelsToLoad, ...
    TimeStep = seconds(1)); 
    %alter data to transform from 1x1 cell to double
    timetable = data{1};
    df = timetable2table(timetable);

    % Fix noise durations
    df = df(~isinf(df.noise_dur), :); % Remove rows where noise_dur is infinity
    df = df(df.trial > 60, :); % Remove rows for first 60 trials
    df.noise_dur(df.outcome == 3) = 14; % Set noise_dur to 14 where outcome is 3

    % Create groups based on 'subject_id' and 'reward'
    [groups, subject_id, reward] = findgroups(df.subject_id, df.reward);

    % Initialize a cell array or table to store results
    groupResults = cell(max(groups), 1);

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

    % Define the filename for the Excel file
    outputFileName = 'combinedResults.xlsx';
    % Create the full file path for the Excel file
    outputPath = fullfile(project_dir, outputFileName);

    % Save combinedResults to an Excel file
    writetable(combinedResults, outputPath);
end

function res2 = tr_sdtv2(df)
    % tr_sdt v2: Signal detection theory analysis with survival analysis computation

    % Set signal duration and interpolation resolution
    signal_dur = 3;
    dt = 0.01;

    % Load event data
    T = [df.rt(df.outcome == 2);
         df.noise_dur(df.outcome == 0);
         df.noise_dur(df.outcome == 1);
         df.noise_dur(df.outcome == 3)];

    % Grab KM table from lifelines
    km_far = lifelines(df);

    % Scale the Kaplan-Meier estimates by the number of observations
    numObservations = length(T);
    km_far.KM_t = km_far.KM * numObservations;

    % Convert to cumulative probability
    km_far.CD = 1 - km_far.KM;
    km_far.CD_t = numObservations - km_far.KM_t;

    % Compute conditional probabilities
    numPlacesToShift = round(signal_dur / dt);
    km_far.CD_shift = [km_far.CD(numPlacesToShift+1:end); nan(numPlacesToShift, 1)];
    km_far.CD_conditional = (km_far.CD_shift - km_far.CD) ./ (1 - km_far.CD);
    km_far.CD_shift_t = [km_far.CD_t(numPlacesToShift+1:end); nan(numPlacesToShift, 1)];
    km_far.CD_conditional_t = (km_far.CD_shift_t - km_far.CD_t + 0.5) ./ (numObservations - km_far.CD_t + 1);

    % Select a subset of data based on conditions
    d = df((df.outcome == 0) | (df.outcome == 1), :);

    % Calculate hit rate (hr) and false alarm rate (far)
    hr = (sum(d.outcome == 0) + 0.5) / (sum(d.outcome == 0) + sum(d.outcome == 1) + 1);

    % Calculate false alarm rate (far)
    [~, closestIdx] = arrayfun(@(x) min(abs(km_far.time - x)), d.noise_dur);
    validIdx = ~isnan(closestIdx);
    far = mean(km_far.CD_conditional_t(closestIdx(validIdx)), 'omitnan');

    % Calculate SDT measures (d' and c)
    tr_d = norminv(hr) - norminv(far);
    tr_c = -0.5 * (norminv(hr) + norminv(far));

    % Return results
    res2 = table(hr, far, tr_d, tr_c, 'VariableNames', {'tr_hr', 'tr_far', 'tr_d', 'tr_c'});
end

function kmf_far = lifelines(df)
    % survival analysis computation: Computes the survival function of licking during noise

    % Resolution for interpolation
    dt = 0.01;

    T = [df.rt(df.outcome == 2);
         df.noise_dur(df.outcome == 0);
         df.noise_dur(df.outcome == 1);
         df.noise_dur(df.outcome == 3)];

    E = [ones(sum(df.outcome == 2), 1);
         zeros(sum(df.outcome == 0), 1);
         zeros(sum(df.outcome == 1), 1);
         zeros(sum(df.outcome == 3), 1)];

    % Sort the data
    [sorted_T, sort_index] = sort(T);
    sorted_E = E(sort_index);

    % Calculate the estimator
    unique_times = unique(sorted_T);
    n = length(sorted_T);
    S = zeros(length(unique_times), 1);
    ni = n;
    % Calculation of survival function
    for i = 1:length(unique_times)
        ti = unique_times(i);
        di = sum(sorted_T == ti & sorted_E == 1);
        S(i) = (1 - di / ni);
        ni = ni - sum(sorted_T == ti);
    end

    % Cumulative product to get the survival probability
    S = cumprod(S);

    % KM table: Interpolate at resolution set above
    time_points = linspace(0, 14, ((14/dt) + 1));
    KM = interp1(unique_times, S, time_points, 'previous');

    % Set survival function to = 1 at t = 0
    KM(1) = 1;

    % Output table
    kmf_far = table(time_points', KM', 'VariableNames', {'time', 'KM'});
end
