function demo_plot(task, field)
% Plot a completed HyMoS-Baltamatica task. Call once for each field during recording.
r = hymos_result(task);
if ~strcmp(char(r.metadata.state), 'converged')
    error('The demonstration task has not converged.');
end
figure(1);
clf;
if strcmp(field, 'temperature')
    plot(r.x, r.temperature, 'b-', 'LineWidth', 2);
    xlabel('Position x1', 'Interpreter', 'none'); ylabel('Temperature', 'Interpreter', 'none');
    title('Couette flow: temperature', 'Interpreter', 'none');
elseif strcmp(field, 'velocity')
    plot(r.x, r.velocity(:, 2), 'b-', 'LineWidth', 2);
    xlabel('Position x1', 'Interpreter', 'none'); ylabel('Tangential velocity', 'Interpreter', 'none');
    title('Couette flow: tangential velocity', 'Interpreter', 'none');
elseif strcmp(field, 'residual')
    k = (0:numel(r.residual)-1)';
    semilogy(k, r.residual, 'b-', 'LineWidth', 2);
    hold on;
    semilogy(k, r.metadata.tolerance * ones(size(k)), 'r--', 'LineWidth', 1.5);
    hold off;
    xlabel('Iteration', 'Interpreter', 'none'); ylabel('Residual', 'Interpreter', 'none');
    title('Couette flow: convergence', 'Interpreter', 'none');
    legend('Residual', 'Tolerance', 'Interpreter', 'none');
else
    error('Choose temperature, velocity, or residual.');
end
grid on;
set(gca, 'FontSize', 14);

end

