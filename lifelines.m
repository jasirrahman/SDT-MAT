%survival analysis computation
%computes the survival function of licking during noise
function kmf_far = lifelines(df)

 %resolution for interpolation
 dt = 0.01;

 T = [df.rt(df.outcome == 2);
        df.noise_dur(df.outcome == 0);
        df.noise_dur(df.outcome == 1);
        df.noise_dur(df.outcome == 3)];

 E = [ones(sum(df.outcome == 2), 1);
     zeros(sum(df.outcome == 0), 1);
     zeros(sum(df.outcome == 1), 1);
     zeros(sum(df.outcome == 3), 1)];

%sort the data
[sorted_T, sort_index] = sort(T);
sorted_E = E(sort_index);

%calculate the estimator
unique_times = unique(sorted_T);
n = length(sorted_T);
S = zeros(length(unique_times), 1);
ni = n;
%calculation of survival function
for i = 1:length(unique_times)
    ti = unique_times(i);
    di = sum(sorted_T == ti & sorted_E == 1);
    S(i) = (1 - di / ni);
    ni = ni - sum(sorted_T == ti);
end

% Cumulative product to get the survival probability
S = cumprod(S);

%km table
%interpolate at resolution set above
time_points = linspace(0, 14, ((14/dt) + 1));
KM = interp1(unique_times, S, time_points, 'previous');
%set survival function to = 1 at t = 0
KM(1) = 1;
%output table
kmf_far = table(time_points', KM', 'VariableNames', {'time', 'KM'});

end
