function [fourier, coupled] = demo_cases(example_dir)
% Run the two bundled HyMoS-Baltamatica examples without changing GUI-applied parameters.
% These checks cover convergence and field consistency, not reference accuracy.
if nargin == 0
    example_dir = 'examples/ParallelPlate1D';
end
configs = {'fourier.txt', 'couette_fourier.txt'};
names = {'Fourier', 'Couette-Fourier'};
results = cell(1, 2);
disp('Bundled-example functional checks');
for k = 1:2
    path = fullfile(example_dir, configs{k});
    check = hymos_validate('couette', path);
    if ~check.valid
        error([names{k} ': ' char(check.message)]);
    end
    task = hymos_submit('couette', path);
    status = hymos_wait(task, 60);
    if ~strcmp(char(status.state), 'converged')
        error([names{k} ': ' char(status.state) ' - ' char(status.message)]);
    end
    r = hymos_result(task);
    if numel(r.x) ~= 128 || size(r.velocity, 2) ~= 3
        error([names{k} ': unexpected result dimensions.']);
    end
    values = [r.x(:); r.density(:); r.temperature(:); r.velocity(:); r.heat_flux(:); r.residual(:)];
    if any(isnan(values)) || any(isinf(values))
        error([names{k} ': non-finite result values.']);
    end
    if r.residual(end) > r.metadata.tolerance || any(r.temperature <= 0) || any(r.density <= 0)
        error([names{k} ': residual or positivity check failed.']);
    end
    if k == 1 && max(abs(r.velocity(:, 2))) > 1e-10
        error('Fourier: tangential-velocity consistency check failed.');
    end
    if k == 2 && ~(r.velocity(1, 2) < 0 && r.velocity(end, 2) > 0)
        error('Couette-Fourier: wall-driven velocity direction check failed.');
    end
    files = hymos_export(task);
    if exist(char(files.fields), 'file') ~= 2 || exist(char(files.residual), 'file') ~= 2
        error([names{k} ': export-file check failed.']);
    end
    fprintf('%s: PASS | %d iterations | residual %.3e | %.2f s\n', ...
        names{k}, status.iteration, status.residual, status.elapsed_seconds);
    results{k} = r;
end
fourier = results{1};
coupled = results{2};
disp('Both examples: converged, finite fields, CSV export verified.');
demo_compare(fourier, coupled);
end

