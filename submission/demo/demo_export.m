function files = demo_export(task)
% Export through the actual HyMoS-Baltamatica interface and verify the generated files.
files = hymos_export(task);
paths = {char(files.fields), char(files.residual), char(files.metadata), char(files.distribution)};
for k = 1:numel(paths)
    if exist(paths{k}, 'file') ~= 2
        error(['Export file is missing: ' paths{k}]);
    end
end
disp('Export completed:');
disp('  couette_fields.csv  - macroscopic fields');
disp('  residual.csv        - iteration history');
disp('  result_metadata.txt - solver settings');
disp('  DisSol1th.dat       - distribution coefficients');
disp(['Output folder: ' char(files.workspace)]);
end
