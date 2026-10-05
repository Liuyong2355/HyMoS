function status = hymos_monitor(task, refresh_seconds)
% Foreground display; does not own or stop the asynchronous calculation.
if nargin < 2
    refresh_seconds = 0.5;
end
if ~isnumeric(refresh_seconds) || ~isreal(refresh_seconds) || ...
        numel(refresh_seconds) ~= 1 || isnan(refresh_seconds) || ...
        isinf(refresh_seconds) || refresh_seconds < 0.05 || refresh_seconds > 60
    error('refresh_seconds must be a finite real scalar in [0.05, 60] seconds');
end

previous_length = 0;
status = hymos_status(task);
while true
    state = char(status.state);
    residual_text = hymos_progress_number(status.residual);
    tolerance_text = 'n/a';
    % Preserve helper compatibility with status structs from older packages.
    if isfield(status, 'tolerance')
        tolerance_text = hymos_progress_number(status.tolerance);
    end
    line = sprintf( ...
        'State: %s | Iteration: %d | Residual: %s | Target: %s | Elapsed: %.2f s', ...
        state, status.iteration, residual_text, tolerance_text, status.elapsed_seconds);
    padding = repmat(' ', 1, max(0, previous_length - length(line)));
    fprintf('\r%s%s', line, padding);
    previous_length = length(line);

    if ~strcmp(state, 'pending') && ~strcmp(state, 'running')
        fprintf('\n');
        if strcmp(state, 'converged')
            disp('HyMoS-Baltamatica converged: target tolerance reached.');
        elseif strcmp(state, 'iteration_limit')
            disp('HyMoS-Baltamatica iteration limit reached: convergence was not achieved.');
        elseif strcmp(state, 'stopped')
            disp('HyMoS-Baltamatica stopped at a safe iteration boundary; result may be partial.');
        elseif strcmp(state, 'busy')
            disp('HyMoS-Baltamatica busy: another Couette solve is active; this task did not run.');
        elseif strcmp(state, 'invalid_config')
            disp('HyMoS-Baltamatica invalid configuration: calculation was not started.');
        else
            disp('HyMoS-Baltamatica failed: calculation did not complete successfully.');
        end
        disp(['Details: ' char(status.message)]);
        return;
    end
    status = hymos_wait(task, refresh_seconds);
end
end

function text = hymos_progress_number(value)
if isempty(value) || ~isnumeric(value) || ~isreal(value) || ...
        numel(value) ~= 1 || isnan(value) || isinf(value)
    text = 'n/a';
else
    text = sprintf('%.6e', value);
end
end
