N = 100;
beta_values = 0.2 + 0.02 * (1:20);
num_iterations = 1000;
num_plots = 4;
selected_beta_indices = [6, 11, 13, 20];

X_sequences = zeros(N, N, num_iterations, length(beta_values));
v_hat = zeros(1, length(beta_values));

for b = 1:length(beta_values)
    beta = beta_values(b);
    lattice = randi([0, 1], N, N) * 2 - 1;
    for i = 1:num_iterations
        for x = 1:N
            for y = 1:N
                neighbors = lattice(mod(x-2, N)+1, y) + lattice(mod(x, N)+1, y) + ...
                            lattice(x, mod(y-2, N)+1) + lattice(x, mod(y, N)+1);
                p = exp(beta * neighbors) ./ (exp(beta * neighbors) + exp(-beta * neighbors));
                lattice(x, y) = 2 * (rand < p) - 1;
            end
        end
        X_sequences(:, :, i, b) = lattice;
    end
    v_hat(b) = mean(sum(sum(X_sequences(:, :, :, b), 1), 2) / (N^2), 'all');
end

disp('20 values of v(beta):');
disp(v_hat);

figure;
for i = 1:num_plots
    beta_index = selected_beta_indices(i);
    lattice_1000 = X_sequences(:, :, end, beta_index);
    subplot(2, 2, i);
    imagesc(lattice_1000);
    colormap gray;
    title(['\beta = ', num2str(beta_values(beta_index))]);
end

figure;
plot(beta_values, v_hat, '-o');
xlabel('\beta');
ylabel('v(\beta)');
title('v(\beta) vs \beta');
grid on;