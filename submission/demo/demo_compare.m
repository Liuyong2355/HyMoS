function demo_compare(fourier, coupled)
figure(1); clf;
plot(fourier.x, fourier.temperature, 'b-', 'LineWidth', 2); hold on;
plot(coupled.x, coupled.temperature, 'r-', 'LineWidth', 2); hold off;
xlabel('Position x1', 'Interpreter', 'none');
ylabel('Temperature', 'Interpreter', 'none');
title('Fourier / Couette-Fourier: temperature', 'Interpreter', 'none');
legend('Fourier', 'Couette-Fourier', 'Interpreter', 'none');
grid on; set(gca, 'FontSize', 14);
end
