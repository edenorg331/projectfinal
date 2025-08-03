%% Parameters
params.minPeakHeight_accY = 1.4;
params.minPeakHeight_gyroY = 1.3;
params.minPeakHeight_accX = 1;
params.ForwardStepSize = 0.75;
params.BackwardStepSize = 0.6;
params.SideStepSize = 0.75/2;
params.StairHeight = 0.158;
params.StairStepSize = 0.30;
params.fs = 55;
params.max_allowed_offset = 0.3;
params.sgolay_order = 7;
params.sgolay_window = 31;
plot_params.gravity = 9.8;
plot_params.line_width = 1.5;
plot_params.marker_size = 8;
plot_params.color_map = containers.Map( ...
{'forward','backward','up stairs','down stairs','turn left','turn right','side left','side right'}, ...
{[0 0 1], [1 0 0], [0 0.6 0], [1 0 1], [0 1 1], [1 0.5 0], [0 1 1], [0 0 0]} ...
);


threshold_back_ratio = 0.8;
initial_freq_window_size = 5;

%% Setup FTP parameters
ftp_host = '10.100.102.2';
ftp_port = '2368';
ftp_user = 'pc';
ftp_pass = '123456';

remote_path = '/device/Download/';

local_folder = tempdir; % תיקיה זמנית במחשב
image_path = local_folder;

%% Loop - Real-time reading
buffer_window_sec = 0.5;
last_processed_time = -Inf;

prev_num_images = 0;
current_csv_file = '';
last_csv_date = datetime(0, 'ConvertFrom','posixtime');
  
%% Parameters for Real-time Processing
buffer_window_sec = 0.5; % how much time to process backwards
last_processed_time = -Inf; % initialize to inf before processing any data

%% Loop - Real-time reading
all_detected_steps_times = [];
all_trajectory = [0 0 0];
all_motions = {};
plot_every_n_cycles = 5; %  every 5 cycles update FixTrajectory + PlotMotionData
cycle_counter = 0;

while true
    %% Step 1: If needed — find latest SensorData_*.csv
    if isempty(current_csv_file)
        % ftp files
        list_cmd = sprintf('curl -s --user %s:%s ftp://%s:%s%s', ...
            ftp_user, ftp_pass, ftp_host, ftp_port, remote_path);
        [~, list_output] = system(list_cmd);

        % extract file name
        file_lines = split(string(list_output), newline);
        file_list = [];
        for i = 1:length(file_lines)
            line = strtrim(file_lines(i));
            if contains(line, 'SensorData_') && contains(line, '.csv')
                parts = split(line);
                filename_candidate = parts(end);
                file_list = [file_list; filename_candidate];
            end
        end

        if isempty(file_list)
            disp('No SensorData CSV files found on FTP.');
            pause(3);
            continue;
        else
            % take latest file
            latest_csv_name = sort(file_list);
            latest_csv_name = latest_csv_name(end);
            current_csv_file = latest_csv_name;
            disp(['[CSV] Initial CSV file: ' current_csv_file]);
        end
    end

    %% Step 2: download + read current_csv_file
    try
        download_csv_cmd = sprintf('curl -s --user %s:%s ftp://%s:%s%s%s -o "%s"', ...
            ftp_user, ftp_pass, ftp_host, ftp_port, remote_path, current_csv_file, ...
            fullfile(local_folder, current_csv_file));
        system(download_csv_cmd);

        % Read table
        local_csv_path = fullfile(local_folder, current_csv_file);
        data = readtable(local_csv_path);

        % Extract columns
        time = data{:,1};
        acc = data{:,2:4};
        gyro = data{:,5:7};
        baro = data{:,8};
        mag = data{:,9:11};

        if max(time) > 1e5
            time = time / 1000; % ms to s
        end

        % Select only new data
        if last_processed_time == -Inf
            idx_to_process = true(size(time));
        else
            idx_to_process = time >= (last_processed_time - buffer_window_sec);
        end

        time_proc = time(idx_to_process);
        acc_proc = acc(idx_to_process,:);
        gyro_proc = gyro(idx_to_process,:);
        baro_proc = baro(idx_to_process);
        mag_proc = mag(idx_to_process,:);

        last_processed_time = max(time_proc);

        disp(['[CSV] Processing ' num2str(length(time_proc)) ' samples (up to time ' num2str(last_processed_time) ' sec)']);

        % Smooth signals
        acc_smooth = sgolayfilt(acc_proc, params.sgolay_order, params.sgolay_window);
        gyro_smooth = sgolayfilt(gyro_proc, params.sgolay_order, params.sgolay_window);
        mag_smooth = sgolayfilt(mag_proc, params.sgolay_order, params.sgolay_window);
        baro_smooth = sgolayfilt(baro_proc, params.sgolay_order, params.sgolay_window);

        % Heading estimation
        heading_deg = rad2deg(atan2(mag_smooth(:,3),mag_smooth(:,1)));
        heading_deg(heading_deg < 0) = heading_deg(heading_deg < 0) + 360;

        % Detect steps
        [step_peaks_accY, steps_locs_accY, step_times_accY, ...
         step_peaks_accX, steps_locs_accX, peaks_accX, peaks_times_accX, acc_y] = ...
         detect_steps(acc_smooth, time_proc, params);

        % Detect turns
        [gyro_y, peaks_gyroY, peaks_locs_gyroY, peaks_times_gyroY, ...
         sorted_peaks_time, sorted_peaks_locs, sorted_peaks_height, ...
         sorted_source_labels, start_loc, end_loc] = ...
         detect_turns(gyro_smooth, time_proc, step_times_accY, steps_locs_accY, step_peaks_accY, params);

        % Image timestamps and assignment %%eden
        [image_files, image_timestamps] = extract_image_timestamps(image_path);
        [image_files, image_timestamps] = prepare_image_metadata(image_files, image_timestamps);
        step_images = assign_images_to_steps(image_files, image_timestamps, step_times_accY, folder, params.max_allowed_offset);

        % Detect stairs
        [stair_movement, height, delta_alt, StairsPeakHeightAccY] = ...
            StairsMovementDetection(mean(step_peaks_accY), baro_proc, sorted_peaks_locs, sorted_peaks_height, params);

        % Detect Forward Average Frequency and Step Height
        [initial_frequency, AverageStepForward, step_indices_accY, step_times_accY_sorted] = ...
            SetThresholds(sorted_source_labels, sorted_peaks_time, step_peaks_accY, initial_freq_window_size);

        % Detect directions
        directions = {}; %%eden
        directions_images = repmat({'unknown'}, length(directions), 1); 
        directions = DetectStepDirection(sorted_source_labels, sorted_peaks_height, ...
            step_indices_accY, step_times_accY_sorted, threshold_back_ratio, ...
            initial_frequency, stair_movement, StairsPeakHeightAccY, step_images, directions_images);

        % Build partial trajectory
        [trajectory, motions, directions_images] = BuildTrajectory( ...
            directions, step_images, delta_alt, ...
            heading_deg, start_loc, end_loc, time_proc, gyro_y, sorted_source_labels, params);


        % Add only new steps
        new_steps_idx = ~ismember(step_times_accY, all_detected_steps_times);
        new_step_times = step_times_accY(new_steps_idx);

        if ~isempty(new_step_times)
            disp(['[NEW STEPS] Adding ' num2str(sum(new_steps_idx)) ' new steps']);
            all_detected_steps_times = [all_detected_steps_times; new_step_times];
            all_trajectory = [all_trajectory; trajectory(2:end,:)];
            all_motions = [all_motions, motions];
        else
            disp('[NEW STEPS] No new steps detected in this batch.');
        end

        % Plot live trajectory
        figure(100); clf;
        plot3DTrajectoryWithMotionColors(all_trajectory, all_motions, plot_params);
        title(sprintf('Trajectory (Total Steps: %d)', length(all_motions)));
        drawnow;

        % Every N cycles — FixTrajectory and PlotMotionData
        cycle_counter = cycle_counter + 1;
        if mod(cycle_counter, plot_every_n_cycles) == 0
            disp('[FixTrajectory] Running FixTrajectory & PlotMotionData');
            [all_trajectory, all_motions] = FixTrajectory(all_trajectory, all_motions, sorted_source_labels, sorted_peaks_time);
            PlotMotionData(time_proc, acc_proc, acc_smooth, acc_y, step_times_accY, step_peaks_accY, ...
                height, sorted_peaks_time, sorted_peaks_locs, gyro_proc, gyro_y, ...
                peaks_times_gyroY, peaks_locs_gyroY, acc_proc(:,1), acc_smooth(:,1), peaks_times_accX, steps_locs_accX, ...
                heading_deg, all_trajectory, all_motions, plot_params);
        end

    catch ME
        disp('Error reading CSV or processing:');
        disp(ME.message);
    end

    %% Step 3: Download latest images — רק הורדה, לא להציג
    try
        image_list_cmd = sprintf('curl -s --user %s:%s ftp://%s:%s%s', ...
            ftp_user, ftp_pass, ftp_host, ftp_port, remote_path);
        [~, image_output] = system(image_list_cmd);

        image_lines = split(string(image_output), newline);
        image_list = [];
        for i = 1:length(image_lines)
            line = strtrim(image_lines(i));
            if contains(line, 'step_photo_') && contains(line, '.jpg')
                parts = split(line);
                img_candidate = parts(end);
                image_list = [image_list; img_candidate];
            end
        end

        num_images = numel(image_list);

        if num_images > prev_num_images
            disp(['[Images] Found ' num2str(num_images) ' images']);
            for i = 1:num_images
                img_name = strtrim(image_list(i));
                download_img_cmd = sprintf('curl -s --user %s:%s ftp://%s:%s%s%s -o "%s"', ...
                    ftp_user, ftp_pass, ftp_host, ftp_port, remote_path, img_name, ...
                    fullfile(local_folder, img_name));
                system(download_img_cmd);
            end
            prev_num_images = num_images;
        else
            disp(['[Images] No new images (Total: ' num2str(num_images) ')']);
        end
    catch ME
        disp('Error reading images:');
        disp(ME.message);
    end

    %% Pause
    pause(3);
end


%% Functions
function [step_peaks_accY, steps_locs_accY, step_times_accY, ...
          step_peaks_accX, steps_locs_accX, peaks_accX, peaks_times_accX, acc_y] = detect_steps(acc_smooth, time, params)
    acc_wo_gravity = acc_smooth;
    acc_wo_gravity(:,2) = acc_smooth(:,2) - 9.8;
    acc_y = acc_wo_gravity(:,2);
    [step_peaks_accY, steps_locs_accY] = findpeaks(acc_y, 'MinPeakHeight', params.minPeakHeight_accY, 'MinPeakDistance', params.fs);
    step_peaks_accY = acc_y(steps_locs_accY);
    step_times_accY = time(steps_locs_accY);
    [step_peaks_accX, steps_locs_accX] = findpeaks(abs(acc_smooth(:,1)), 'MinPeakHeight', params.minPeakHeight_accX, 'MinPeakDistance', params.fs);
    peaks_accX = acc_smooth(steps_locs_accX,1);
    peaks_times_accX = time(steps_locs_accX);
end

function [gyro_y, peaks_gyroY, peaks_locs_gyroY, peaks_times_gyroY, ...
          sorted_peaks_time, sorted_peaks_locs, sorted_peaks_height, sorted_source_labels, ...
          start_loc, end_loc] = detect_turns(gyro_smooth, time, step_times_accY, steps_locs_accY, step_peaks_accY, params)
    gyro_y = gyro_smooth(:,2);
    [peaks_gyroY, peaks_locs_gyroY, ~] = findpeaks(abs(gyro_y), 'MinPeakHeight', params.minPeakHeight_gyroY, 'MinPeakDistance', 1.5*params.fs, 'MinPeakWidth', params.fs);
    peaks_gyroY = gyro_y(peaks_locs_gyroY);
    peaks_times_gyroY = time(peaks_locs_gyroY);
    source_labels = [ones(1, length(step_times_accY)), 2 * ones(1, length(peaks_times_gyroY))];
    merged_peaks_locs = [steps_locs_accY; peaks_locs_gyroY];
    merged_peaks_time = [step_times_accY; peaks_times_gyroY];
    merged_peaks_height = [step_peaks_accY; peaks_gyroY];
    [sorted_peaks_time, sort_peaks_idx] = sort(merged_peaks_time);
    sorted_peaks_locs = merged_peaks_locs(sort_peaks_idx);
    sorted_peaks_height = merged_peaks_height(sort_peaks_idx);
    sorted_source_labels = source_labels(sort_peaks_idx);
    zero_locs_gyro_y = time(abs(gyro_y) < 0.05);
    start_loc = zeros(1,length(peaks_locs_gyroY));
    end_loc = zeros(1,length(peaks_locs_gyroY));
    for peak_i = 1:length(peaks_locs_gyroY)
        peak_time = peaks_times_gyroY(peak_i);
        zero_left = zero_locs_gyro_y(zero_locs_gyro_y < peak_time);
        if ~isempty(zero_left)
            start_loc(peak_i) = zero_left(end);
        end
        zero_right = zero_locs_gyro_y(zero_locs_gyro_y > peak_time);
        if ~isempty(zero_right)
            end_loc(peak_i) = zero_right(1);
        end
    end
end

function [image_files, image_timestamps] = extract_image_timestamps(image_path)
    image_files = image_path;
    image_timestamps = zeros(length(image_files), 1);
    for k = 1:length(image_files)
        tokens = regexp(char(image_files(k).name), 'step_photo_(\d+)', 'tokens');
        if ~isempty(tokens)
            image_timestamps(k) = str2double(tokens{1}{1});
        else
            warning('Could not extract timestamp from file name: %s', image_files(k).name);
            image_timestamps(k) = NaN;
        end
    end
end

function [image_files, image_timestamps] = prepare_image_metadata(image_files, image_timestamps)
    valid_idx = ~isnan(image_timestamps);
    image_files = image_files(valid_idx);
    image_timestamps = image_timestamps(valid_idx);
    [image_timestamps, sort_idx] = sort(image_timestamps);
    image_files = image_files(sort_idx);
    image_timestamps = image_timestamps / 1000;
end

function step_images = assign_images_to_steps(image_files, image_timestamps, step_times_accY, folder, max_allowed_offset)
    step_images = cell(length(step_times_accY), 1);
    % Match each step to the nearest image timestamp within a threshold
    for i = 1:length(step_times_accY)
        [min_diff, idx] = min(abs(image_timestamps - step_times_accY(i)));
        if min_diff < max_allowed_offset
            step_images{i} = fullfile(folder, image_files(idx).name);
        else
            step_images{i} = '';
        end
    end
    % Fix corrupted images
    for i = 1:length(step_images)
        if ~isempty(step_images{i})
            original = step_images{i};
            repaired = fullfile(folder, ['fixed_' extractAfter(original, length(folder)+2)]);
            fixed = fix_corrupted_jpeg(original, repaired);
            if fixed
                step_images{i} = repaired;
            end
        end
    end
end

function [initial_frequency, AverageStepForward, step_indices_accY, step_times_accY_sorted] = ...
    SetThresholds(sorted_source_labels, sorted_peaks_time, step_peaks_accY, initial_freq_window_size)

    % Find indices of steps from accelerometer Y
    step_indices_accY = find(sorted_source_labels == 1);
    step_times_accY_sorted = sorted_peaks_time(step_indices_accY);

    % Initialize defaults
    step_indices_accY = find(sorted_source_labels == 1);
    step_times_accY_sorted = sorted_peaks_time(step_indices_accY);
    initial_frequency = NaN;
    AverageStepForward = NaN;

    % Check if there are enough steps
    if length(step_times_accY_sorted) >= initial_freq_window_size
        % Compute average interval and frequency
        initial_times = step_times_accY_sorted(1:5);
        initial_intervals = diff(initial_times);
        initial_avg_interval = mean(initial_intervals);
        initial_frequency = 1 / initial_avg_interval;

        % Estimate average step amplitude
        AverageStepForward = mean(step_peaks_accY(1:5));
        %fprintf('average step forward %f\n',AverageStepForward);
    end
end

function [stair_movement, height, delta_alt, StairsPeakHeightAccY] = ...
    StairsMovementDetection(AverageStepForward, baro, sorted_peaks_locs, sorted_peaks_height, params)

    % Compute stair height threshold from average step
    StairsPeakHeightAccY = 2.5 * AverageStepForward;

    % Estimate altitude from barometer
    height = 44330 .* (1 - (baro / 1013.25).^0.1903);  % P0 = 1013.25 hPa

    % Parameters
    window_size = 2 * params.fs;  % samples before/after the step to analyze (adjust if needed)
    num_steps = length(sorted_peaks_locs);
    delta_alt = zeros(num_steps, 1);

    for i = 1:num_steps
        idx = sorted_peaks_locs(i);
        % Ensure bounds
        before_idx = max(1, idx - window_size);
        after_idx  = min(length(height), idx + window_size);
        
        alt_before = median(height(before_idx:idx));
        alt_after  = median(height(idx:after_idx));
        
        delta_alt(i) = alt_after - alt_before;
    end

    % Determine stair movement based on AccY and Altitude
    stair_movement = (abs(sorted_peaks_height) > StairsPeakHeightAccY) | (abs(delta_alt) >= 0.2);
    %disp(delta_alt);
end


function directions = DetectStepDirection( sorted_source_labels, sorted_peaks_height, ...
    step_indices_accY, step_times_accY_sorted, threshold_back_ratio, initial_frequency, ...
    stair_movement, StairsPeakHeightAccY, step_images, directions_images )

    directions = cell(length(sorted_source_labels), 1);
    fprintf('initial frequecncy is - %.5f\n', initial_frequency);
    for i = 1:length(sorted_source_labels)
        if sorted_source_labels(i) == 2  % Gyroscope-based detection
            if sorted_peaks_height(i) > 1.3
                directions{i} = 'turn left';
            elseif sorted_peaks_height(i) < -1.3
                directions{i} = 'turn right';
            else
                directions{i} = 'gyro event';
            end

        elseif sorted_source_labels(i) == 1  % Accelerometer-based detection
            accY_idx = find(step_indices_accY == i);

            if accY_idx > 5
                window_times = step_times_accY_sorted(accY_idx-1:accY_idx);
                window_intervals = diff(window_times);
                window_avg_interval = mean(window_intervals);
                current_frequency = 1.0 / window_avg_interval;
                fprintf('step %d: frequecncy = %.5f, StepHeight = %.5f\n', accY_idx, current_frequency, sorted_peaks_height(accY_idx));
                if (current_frequency < threshold_back_ratio * initial_frequency) && ...
                        (abs(sorted_peaks_height(i)) < 3)
                    directions{i} = 'backward';
                elseif ((current_frequency >= 0.85 * initial_frequency && ~stair_movement(i))) 
                    directions{i} = 'forward';
                elseif (stair_movement(i))
                    directions{i} = 'stairs';
                else
                    directions{i} = 'unknown';
                end
            else
                directions{i} = 'forward';  % Default for early steps
            end
    
            if i < length(step_images) && ~isempty(step_images{i}) && ~isempty(step_images{i+1})
                I1 = rgb2gray(imread(step_images{i}));
            
                if (sorted_source_labels(i+1) == 2)
                    I2 = rgb2gray(imread(step_images{i+2}));
                    j = j + 1;
                else
                    I2 = rgb2gray(imread(step_images{i+1}));
                end

                % Optical Flow Estimation between I1 and I2
                opticFlow = opticalFlowFarneback('NumPyramidLevels', 4, ...
                                     'PyramidScale', 0.7, ...
                                     'NumIterations', 7, ...
                                     'NeighborhoodSize', 7, ...
                                     'FilterSize', 9);
                estimateFlow(opticFlow, I1);     % initialize with first image
                flow = estimateFlow(opticFlow, I2);  % compute flow from I1 to I2
            
                % Compute average flow vector
                Vx = flow.Vx;
                Vy = flow.Vy;
            
                % Mask low-magnitude flow to reduce noise
                mag = sqrt(Vx.^2 + Vy.^2);
                threshold = 0.3; % adjust as needed
                mask = mag > threshold;
            
                % Weighted mean flow in x and y
                dx = mean(Vx(mask), 'omitnan');
                dy = mean(Vy(mask), 'omitnan');
            
                if isnan(dx) || isnan(dy)
                    fprintf("Optical flow unreliable at step %d\n", i);
                    continue;
                end
            
                % Interpret direction
                fprintf('dx- %.5f, dy - %.5f\n', dx,dy);
                
                if (abs(dx) > 1.2 * abs(dy) && abs(dx) > 0.5)
                    if dx > 0
                        directions{i} = 'side right'; %override sensors if side motion detected
                        directions_images{i} = 'side right';
                        disp('--------------------->Moved Right');
                    else
                        directions{i} = 'side left'; %override sensors if side motion detected
                        directions_images{i} = 'side left';
                        disp('<--------------------Moved Left');
                    end
                % else
                %     if dy > 0
                %         directions_images{i} = 'Backwards/down';
                %         disp('Moved Down (likely Backward)');
                %     else
                %         directions_images{i} = 'Forwards/up';
                %         disp('Moved Up (likely Forward)');
                %     end
                end
            end
        else
            directions{i} = 'unknown';
        end
     
    end
end

function [trajectory, motions, directions_images] = BuildTrajectory( ...
    directions, step_images, delta_alt, ...
    heading_deg, start_loc, end_loc, time, gyro_y, sorted_source_labels , params)

    fprintf("Trail summary:\n");

    pos = [0 0 0];
    trajectory = pos;
    motions = {};
    heading_angle = 0;
    sum_height = 0;
    dist = 0;
    curr_motion = '';
    pre_motion = '';
    curr_gyroY_peak_idx = 0;
    directions_images = repmat({'unknown'}, length(directions), 1);

    for i = 1:length(directions)
        pre_motion = curr_motion;
        
        if strcmp(directions{i}, 'stairs')
            fprintf("Step %d: Stairs", i);
            if delta_alt(i) > 0.15
                dx = params.StairStepSize * sind(heading_angle);
                dz = params.StairStepSize * cosd(heading_angle);
                pos(1) = pos(1) + dx;
                pos(2) = pos(2) + params.StairHeight;  % vertical (y)
                pos(3) = pos(3) + dz;

                sum_height = sum_height + params.StairHeight;
                dist = dist + params.StairStepSize;
                curr_motion = 'up stairs';
                fprintf(" + up\n");
            elseif delta_alt(i) < -0.15
                sum_height = sum_height - params.StairHeight;
                dist = dist + params.StairStepSize;
                curr_motion = 'down stairs';
                dx = params.StairStepSize * sind(heading_angle);
                dz = params.StairStepSize * cosd(heading_angle);
                pos(1) = pos(1) + dx;
                pos(2) = pos(2) - params.StairHeight;
                pos(3) = pos(3) + dz;

                fprintf(" + down\n");
            else
                fprintf("\n");
            end

        else
            curr_motion = directions{i};

            if strcmp(curr_motion, 'turn left') || strcmp(curr_motion, 'turn right')
                curr_gyroY_peak_idx = curr_gyroY_peak_idx + 1;
                start_time = start_loc(curr_gyroY_peak_idx);
                end_time   = end_loc(curr_gyroY_peak_idx);
                gyro_segment_time = time(time >= start_time & time <= end_time);
                gyro_segment_y = gyro_y(time >= start_time & time <= end_time);

                if length(gyro_segment_time) >= 2
                    delta_angle_rad = trapz(gyro_segment_time, gyro_segment_y);
                    delta_angle_deg = rad2deg(delta_angle_rad);
                    if delta_angle_deg > 180
                        delta_angle_deg = delta_angle_deg - 360;
                    elseif delta_angle_deg < -180
                        delta_angle_deg = delta_angle_deg + 360;
                    end
                    heading_angle = mod(heading_angle + delta_angle_deg, 360);
                    fprintf('Step %d: %s by %.1f degrees\n', i, curr_motion, delta_angle_deg);
                else
                    fprintf(' → Not enough points for integration\n');
                end
            end

            angle = heading_angle;
            if strcmp(curr_motion, 'forward')
                dx = params.ForwardStepSize * sind(angle);
                dz = params.ForwardStepSize * cosd(angle);
                pos(1) = pos(1) + dx;
                pos(3) = pos(3) + dz;
                dist = dist + params.ForwardStepSize;
            elseif strcmp(curr_motion, 'backward')
                dx = params.BackwardStepSize * sind(angle);
                dz = params.BackwardStepSize * cosd(angle);
                pos(1) = pos(1) - dx;
                pos(3) = pos(3) - dz;
                dist = dist - params.BackwardStepSize;
            elseif strcmp(curr_motion, 'side left')
                dx = params.SideStepSize * sind(angle - 90);
                dz = params.SideStepSize * cosd(angle - 90);
                pos(1) = pos(1) + dx;
                pos(3) = pos(3) + dz;
            elseif strcmp(curr_motion, 'side right')
                dx = params.SideStepSize * sind(angle + 90);
                dz = params.SideStepSize * cosd(angle + 90);
                pos(1) = pos(1) + dx;
                pos(3) = pos(3) + dz;
            end

            fprintf("Step %d: %s\n", i, curr_motion);
        end

        trajectory(end+1,:) = pos;
        motions{end+1} = curr_motion;

        if ~strcmp(pre_motion, curr_motion) && ~strcmp(pre_motion, '')
            fprintf('change in height %.4f, change in distance %.4f\n', sum_height, dist);
            sum_height = 0;
            dist = 0;
        end

        %fprintf('directions images - %s\n', directions_images{i});
    end
end

function [trajectory, motions] = FixTrajectory(trajectory, motions, sorted_source_labels, sorted_peaks_time)
    for i = 2:(length(sorted_source_labels)-1)
        if (sorted_source_labels(i) == 1)
            fprintf('dist is %.6f\n', (sorted_peaks_time(i+1) - sorted_peaks_time(i)));
     %       if ((sorted_peaks_time(i+1) - sorted_peaks_time(i)) > 4) && (sorted_source_labels(i+1) == 1) %%should rest more :)))
     %           motions(i) = '';
            if (~(strcmp(motions(i), motions(i-1))) && strcmp(motions(i-1), motions(i+1))) && (sorted_source_labels(i-1) == 1) && (sorted_source_labels(i+1) == 1) % if motion is the same before and after, fix the error
                motions(i) = motions(i-1);
            end
        end
    end
    %motions = motions(~cellfun(@isempty, motions));
    %trajectory = trajectory(~cellfun(@isempty, motions));
end

function fixed = fix_corrupted_jpeg(input_path, output_path)
    try
        I = imread(input_path);          % Try reading original image
        imwrite(I, output_path, 'jpg');  % Re-encode and save as valid JPEG
        fixed = true;
    catch
        warning('Could not fix image: %s', input_path);
        fixed = false;
    end
end

function PlotMotionData(time, acc, acc_smooth, acc_y, step_times_accY, step_peaks_accY, ...
    height, sorted_peaks_time, sorted_peaks_locs, gyro, gyro_y, ...
    peaks_times_gyroY, peaks_locs_gyroY, accX, acc_smoothX, peaks_times_accX, steps_locs_accX, ...
    heading_deg, trajectory, motions, plot_params)

    % Plot Accelerometer Y
    figure;
    plot(time, acc(:,2)-plot_params.gravity, 'r-', 'DisplayName', 'Raw AccY');
    hold on;
    plot(time, acc_y, 'b-', 'LineWidth', plot_params.line_width, 'DisplayName', 'Filtered AccY');
    plot(step_times_accY, step_peaks_accY, 'ro', 'MarkerSize', plot_params.marker_size, 'LineWidth', 2);
    legend;
    xlabel('Time'); ylabel('AccY'); grid on;

    % Plot Barometer height
    figure;
    plot(time, height, 'r-', 'DisplayName', 'Height');
    hold on;
    plot(sorted_peaks_time, height(sorted_peaks_locs), 'ro', 'MarkerSize', plot_params.marker_size, 'LineWidth', 2);
    legend('Baro', 'Detected Steps');
    xlabel('Time'); ylabel('Height (m)'); grid on;

    % Plot Accelerometer Z
    figure;
    plot(time, acc(:,3), 'r-', 'DisplayName', 'Raw AccZ');
    hold on;
    plot(time, acc_smooth(:,3), 'b-', 'LineWidth', plot_params.line_width, 'DisplayName', 'Filtered AccZ');
    legend; xlabel('Time'); ylabel('AccZ'); grid on;

    % Plot Accelerometer X
    figure;
    plot(time, acc(:,1), 'r-', 'DisplayName', 'Raw AccX');
    hold on;
    plot(time, acc_smoothX, 'b-', 'LineWidth', plot_params.line_width, 'DisplayName', 'Filtered AccX');
    plot(peaks_times_accX, acc_smoothX(steps_locs_accX), 'ro', 'MarkerSize', plot_params.marker_size, 'LineWidth', 2);
    legend; xlabel('Time'); ylabel('AccX'); grid on;

    % Plot Gyroscope Y
    figure;
    plot(time, gyro(:,2), 'r-', 'DisplayName', 'Raw GyroY');
    hold on;
    plot(time, gyro_y, 'b-', 'LineWidth', plot_params.line_width, 'DisplayName', 'Filtered GyroY');
    plot(peaks_times_gyroY, gyro_y(peaks_locs_gyroY), 'ro', 'MarkerSize', plot_params.marker_size, 'LineWidth', 2);
    legend; xlabel('Time'); ylabel('GyroY'); grid on;

    % Plot Heading
    figure;
    plot(time, heading_deg, 'b-', 'LineWidth', plot_params.line_width);
    hold on;
    plot(peaks_times_gyroY, heading_deg(peaks_locs_gyroY), 'ro', 'MarkerSize', plot_params.marker_size, 'LineWidth', 2);
    xlabel('Time'); ylabel('Heading (deg)'); legend('Heading', 'Peaks'); grid on;
end


function plot3DTrajectoryWithMotionColors(trajectory, motions, plot_params)

    if ~isfield(plot_params, 'line_width')
        plot_params.line_width = 2;
    end
    if ~isfield(plot_params, 'color_map')
        plot_params.color_map = containers.Map( ...
            {'forward','backward','up stairs','down stairs','turn left','turn right'}, ...
            {[0 0 1], [1 0 0], [0 0.6 0], [1 0 1], [0 1 1], [1 0.5 0]});
    end

    figure;
    hold on;
    legend_handles = [];
    legend_labels = {};

    for i = 1:(length(trajectory)-1)
        p1 = trajectory(i, :);
        p2 = trajectory(i+1, :);
        motion = motions{i};

        if isKey(plot_params.color_map, motion)
            c = plot_params.color_map(motion);
        else
            c = [0 0 0];
        end

        if any(strcmp(motion, {'turn left','turn right'}))
            h = plot3(p1(3), p1(1), p1(2), 'o', ...
                'MarkerFaceColor', c, 'MarkerEdgeColor', c, 'MarkerSize', 6);
        else
            h = plot3([p1(3) p2(3)], [p1(1) p2(1)], [p1(2) p2(2)], '--', ...
                'Color', c, 'LineWidth', plot_params.line_width);
        end

        if ~ismember(motion, legend_labels)
            legend_handles(end+1) = h;
            legend_labels{end+1} = motion;
        end
    end

    % color map
    num_points = size(trajectory, 1);
    cmap = jet(num_points);  % או parula, hot, hsv וכו'
    for i = 1:num_points
        pos = trajectory(i, :);  % [X Y Z]
        plot3(pos(3), pos(1), pos(2), '*', ...
            'MarkerFaceColor', cmap(i,:), ...
            'MarkerEdgeColor', cmap(i,:), ...
            'MarkerSize', 4);
    end

    % plot
    view(3);
    xlabel('Z (forward)');
    ylabel('X (sideways)');
    zlabel('Y (altitude)');
    title('3D Movement Trajectory with Time-Based Colors');
    axis equal;
    grid on;

    colormap(jet(num_points));
    colorbar('Ticks', [0 1], 'TickLabels', {'Start', 'End'}, 'Location', 'eastoutside');

    legend(legend_handles, legend_labels);
end
