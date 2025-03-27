clear all
close all
clc

%Read in data
directory = '/Users/6008399/Documents/PhD/CESM_collapse/netcdf/';

%Restoring rate data from Snellius stored in netcdfs
lon = ncread('/Users/6008399/Documents/PhD/CESM_collapse/netcdf/lambda_TEMP_transient_1500_endtimes_1700_1702_spacing_2_window_70_SAMBA_all_gridpoints.nc','lon');
end_times = ncread('/Users/6008399/Documents/PhD/CESM_collapse/netcdf/lambda_TEMP_transient_1500_endtimes_1700_1702_spacing_2_window_70_SAMBA_all_gridpoints.nc','end_time');
data = ncread('/Users/6008399/Documents/PhD/CESM_collapse/netcdf/lambda_TEMP_transient_1500_endtimes_1700_1702_spacing_2_window_70_SAMBA_all_gridpoints.nc','lambda_values_1700');
depth = ncread('/Users/6008399/Documents/PhD/CESM_collapse/netcdf/lambda_TEMP_transient_1500_endtimes_1700_1702_spacing_2_window_70_SAMBA_all_gridpoints.nc','depth');
time = ncread('/Users/6008399/Documents/PhD/CESM_collapse/netcdf/lambda_TEMP_transient_1500_endtimes_1700_1702_spacing_2_window_70_SAMBA_all_gridpoints.nc','time_1700');

%lon = ncread('/Users/6008399/Documents/PhD/CESM_collapse/netcdf/lambda_SALT_transient_1500_endtimes_1670_1730_spacing_10_window_70_SAMBA_all_gridpoints.nc','lon');
%end_times = ncread('/Users/6008399/Documents/PhD/CESM_collapse/netcdf/lambda_SALT_transient_1500_endtimes_1670_1730_spacing_10_window_70_SAMBA_all_gridpoints.nc','end_time');
%data = ncread('/Users/6008399/Documents/PhD/CESM_collapse/netcdf/lambda_SALT_transient_1500_endtimes_1670_1730_spacing_10_window_70_SAMBA_all_gridpoints.nc','lambda_values_1700');
%depth = ncread('/Users/6008399/Documents/PhD/CESM_collapse/netcdf/lambda_SALT_transient_1500_endtimes_1670_1730_spacing_10_window_70_SAMBA_all_gridpoints.nc','depth');
%time = ncread('/Users/6008399/Documents/PhD/CESM_collapse/netcdf/lambda_SALT_transient_1500_endtimes_1670_1730_spacing_10_window_70_SAMBA_all_gridpoints.nc','time_1700');

%Lambda is negative following Boers et al., but should be positive
data = -data;

%Replace the very large numbers with nans (are masked in python but change
%to large number in matlab)
threshold = 9e10;  
data(abs(data) > threshold) = NaN;

figure
contourf(lon, depth, squeeze(data(:,:,36))')
colorbar

%Change shape to (time, depth, lon)
data = permute(data, [3, 2, 1]);
window = 70; %Sliding window which is used for calculating the restoring rate on Snellius
end_year = 1700; %CHANGE ACCORDING TO ENDTIME
min_theshold = 1;

fprintf('%f\n', time(1))
fprintf('%f\n', time(end))

figure
plot(time, data(:,10,10))

%% Loop over different start_year_change to test sensitivity to minimum time length for fit, mulitple CPs allowed per time series

CPend = 1580:end_year-window/2 - 20;

change_point_all_start = nan(10, length(CPend), length(depth), length(lon));
change_point_time_all_start = nan(10, length(CPend), length(depth), length(lon));
time_at_zero_all_start = nan(10, length(CPend), length(depth), length(lon));
trend_all_start = nan(10, length(CPend), length(depth), length(lon));
base_all_start = nan(10, length(CPend), length(depth), length(lon));
r_value_all_start = nan(10, length(CPend), length(depth), length(lon));

for year_i = 1:length(CPend)
    fprintf('Index: %d, Year: %d\n', year_i, CPend(year_i))
    for depth_i = 1:length(depth)
        for lon_i = 1:length(lon)

            % Extract the data for the current grid point and remove nans at begin
            data_point = squeeze(data(window/2+1 :end-window/2, depth_i, lon_i));
            time_point = time(window/2+1:end-window/2);

            if isnan(data_point(1))
                continue
            end
            
            % Find change points
            x = findchangepts(data_point, 'Statistic', 'linear', 'MinThreshold', 1);

            % Skip if no changepoint
            if isempty(x)
                continue
            end

            % Iterate over each change point
            for cp_i = 1:length(x)
                % Select change point
                change_point = x(cp_i);

                %Changepoint must be within change point range (1500-CPend)
                if time_point(change_point) > CPend(year_i)
                    continue
                end

                % Linear fit
                data_plot = data_point(change_point:end);
                time_plot = time_point(change_point:end);
                p = polyfit(time_plot, data_plot, 1);

                % If positive trend -> infinite tipping time (i.e. 10e10)
                if p(1) > 0
                    time_at_zero_all_start(cp_i, year_i, depth_i, lon_i) = 10e10;
                    change_point_time_all_start(cp_i, year_i, depth_i, lon_i) = time_point(change_point);
                    trend_all_start(cp_i, year_i, depth_i, lon_i) = p(1);
                    base_all_start(cp_i, year_i, depth_i, lon_i) = p(2);
                    r_matrix = corrcoef(time_plot, data_plot);
                    r_value_all_start(cp_i, year_i, depth_i, lon_i) = r_matrix(1, 2);
                    continue
                end

                % Calculate the r-value
                r_matrix = corrcoef(time_plot, data_plot);
                r_value = r_matrix(1, 2);

                time_zero = -p(2) / p(1); % Solving y = p(1)*x + p(2) = 0
                %time_zero cannot be smaller than starting year of input
                %data
                if time_zero < 1500
                    continue
                end

                % Save results
                change_point_time_all_start(cp_i, year_i, depth_i, lon_i) = time_point(change_point);
                time_at_zero_all_start(cp_i, year_i, depth_i, lon_i) = time_zero;
                trend_all_start(cp_i, year_i, depth_i, lon_i) = p(1);
                base_all_start(cp_i, year_i, depth_i, lon_i) = p(2);
                r_value_all_start(cp_i, year_i, depth_i, lon_i) = r_value;
            end
        end
    end
    fprintf('%.2f\n', time_point(1));
    fprintf('%.2f\n', time_point(end));

end

%%

%indices = [20, 21, 22, 23, 24, 25, 26, 27, 28];
indices = [1, 10, 20, 30, 40, 50, 60, 64, 66];
%indices = [1, 2, 4, 6, 8, 10, 12, 14, 16];

figure
for i = 1:length(indices)
    index = indices(i);

    % Flatten the 2D matrix into a 1D array
    time_at_zero_all_flat = reshape(time_at_zero_all_start(:,index,:,:), 1, []);

    % Create a subplot
    subplot(3, 3, i)
    histogram(time_at_zero_all_flat,100)
    title(['Estimated tipping times (' num2str(length(find(~isnan(time_at_zero_all_start(:,index,:,:))))) ' points, window ' num2str(window) 'y, change point > y' num2str(CPend(index)) ')'])
   

    % Calculate the 95% confidence interval
    ci_lower = prctile(time_at_zero_all_flat, 2.5);
    ci_upper = prctile(time_at_zero_all_flat, 97.5);

    % Add vertical dashed red lines to the histogram
    line([ci_lower ci_lower], ylim, 'Color', 'r', 'LineStyle', '--');
    line([ci_upper ci_upper], ylim, 'Color', 'r', 'LineStyle', '--');
end
