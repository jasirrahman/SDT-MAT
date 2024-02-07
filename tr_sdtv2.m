%tr_sdt v2

function res2 = tr_sdtv2(df)
    %set signal duration and interpolation resolution
    signal_dur = 3;
    dt = 0.01;

    %load event data
    T = [df.rt(df.outcome == 2);
        df.noise_dur(df.outcome == 0);
        df.noise_dur(df.outcome == 1);
        df.noise_dur(df.outcome == 3)];

    %grab KM table from lifelines
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

    % Calculate hit rate (hr)
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
