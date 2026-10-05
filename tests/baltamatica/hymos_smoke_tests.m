function hymos_smoke_tests(example_dir)
% Host integration checks for the public HyMoS-Baltamatica plugin release.
% Run after load_plugin('HyMoS-Baltamatica'). Tests read templates without modifying them.
if nargin == 0
    example_dir = 'examples/ParallelPlate1D';
end
names = {'Couette', 'Fourier', 'Couette-Fourier'};
configs = {'couette.txt', 'fourier.txt', 'couette_fourier.txt'};
for k = 1:3
    path = fullfile(example_dir, configs{k});
    v = hymos_validate('couette', path);
    if ~v.valid
        error([names{k} ': ' char(v.message)]);
    end
    task = hymos_submit('couette', path);
    s = hymos_wait(task, 60);
    if ~strcmp(char(s.state), 'converged')
        hymos_stop(task);
        error([names{k} ': ' char(s.state) ' - ' char(s.message)]);
    end
    r = hymos_result(task);
    if numel(r.x) ~= 128 || size(r.velocity, 2) ~= 3 || size(r.heat_flux, 2) ~= 2
        error([names{k} ': unexpected result dimensions']);
    end
    values = [r.x(:); r.density(:); r.temperature(:); r.velocity(:); r.stress_moments(:); r.heat_flux(:); r.residual(:)];
    if any(isnan(values)) || any(isinf(values)) || any(r.temperature <= 0) || any(r.density <= 0)
        error([names{k} ': invalid physical fields']);
    end
    if isempty(r.residual) || r.residual(end) > r.metadata.tolerance
        error([names{k} ': residual exceeds tolerance']);
    end
    if k == 2 && max(abs(r.velocity(:, 2))) > 1e-10
        error('Fourier: nonzero tangential velocity');
    end
    if k ~= 2 && ~(r.velocity(1, 2) < 0 && r.velocity(end, 2) > 0)
        error([names{k} ': wall-driven velocity direction']);
    end
    files = hymos_export(task);
    if exist(char(files.fields), 'file') ~= 2 || exist(char(files.residual), 'file') ~= 2 || exist(char(files.metadata), 'file') ~= 2 || exist(char(files.distribution), 'file') ~= 2
        error([names{k} ': missing export file']);
    end
    fprintf('%s: PASS | %d iterations | residual %.3e\n', names{k}, s.iteration, s.residual);
end
text = fileread(fullfile(example_dir, 'couette.txt'));
small = hymos_ui_config(text, 'Mesh/Nx', '8');
v = hymos_validate_text('couette', small);
if ~v.valid
    error(['Nx=8 validation: ' char(v.message)]);
end
task = hymos_submit_text('couette', small);
s = hymos_wait(task, 60);
if ~strcmp(char(s.state), 'converged')
    hymos_stop(task);
    error(['Nx=8: ' char(s.state)]);
end
r = hymos_result(task);
if numel(r.x) ~= 8
    error('Nx=8: wrong result dimensions');
end
invalid = hymos_ui_config(text, 'Mesh/Nx', '7');
v = hymos_validate_text('couette', invalid);
if v.valid || isempty(char(v.message))
    error('Nx=7: expected explanatory validation failure');
end
v = hymos_validate('couette', fullfile(example_dir, '__missing_input__.txt'));
if v.valid || isempty(char(v.message))
    error('Missing file: expected explanatory validation failure');
end
disp('Nx=8, rejected Nx=7, missing-input diagnostics: PASS');
disp('HYMOS_HOST_SMOKE_PASS');
end
