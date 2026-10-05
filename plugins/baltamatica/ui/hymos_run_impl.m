function hymos_run_impl()
% Host-side convenience workflow; the Couette worker remains asynchronous.
% Keep all task references out of GUI UserData.
% Baltamatica 2025 evalin(base, exist(...)) checks the wrong scope inside
% functions. Read the value directly and catch a missing/invalid variable.
previous_status = [];
try
    previous_task = evalin('base', 'hymos_task');
    previous_status = hymos_status(previous_task);
catch
    % No existing HyMoS-Baltamatica task in the base workspace.
end
if ~isempty(previous_status)
    previous_state = char(previous_status.state);
    if strcmp(previous_state, 'pending') || strcmp(previous_state, 'running')
        disp('HyMoS-Baltamatica busy: hymos_task is still active; no new task submitted.');
        disp(['Workspace: ' char(previous_status.workspace)]);
        disp('Continue: hymos_monitor(hymos_task);  Stop: hymos_stop(hymos_task);');
        return;
    end
end

try
    text = evalin('base', 'hymos_config_text');
catch
    disp('HyMoS-Baltamatica submission failed: no applied parameters.');
    disp('Open hymos_setup(''1D'') and apply parameters first.');
    return;
end
try
    check = hymos_validate_text('couette', text);
    if ~check.valid
        disp(['HyMoS-Baltamatica submission failed: ' char(check.message)]);
        return;
    end
    task = hymos_submit_text('couette', text);
catch err
    disp('HyMoS-Baltamatica submission failed. Existing hymos_task was preserved.');
    disp(err);
    return;
end

try
    % This must precede all foreground waiting/printing. If Ctrl+C unwinds
    % the script, the base workspace still owns the task's CStruct reference.
    assignin('base', 'hymos_task', task);
catch err
    hymos_stop(task);
    disp('HyMoS-Baltamatica could not save hymos_task; safe stop requested for the new task.');
    disp(err);
    return;
end
status = hymos_status(task);
if strcmp(char(status.state), 'busy')
    status = hymos_monitor(task);
    return;
end
disp('HyMoS-Baltamatica submission successful. Task saved as hymos_task.');
disp(['Workspace: ' char(status.workspace)]);
disp('Foreground monitoring; the solver runs asynchronously in the background.');
disp('Ctrl+C interrupts monitoring. Query: s = hymos_status(hymos_task);');
disp('Stop calculation: hymos_stop(hymos_task);');
status = hymos_monitor(task);
end
