function gauss_fit = singleGaussian(tuning_data, color)
% Define stimulus indices (x)
stimulus = linspace(-floor(length(tuning_data)/2), floor(length(tuning_data)/2), length(tuning_data));

% Initialize parameters
A = max(tuning_data);  % Amplitude of the Gaussian
sigma = 1;             % Width of the Gaussian
offset = mean(tuning_data); % Offset (baseline) value

% Learning rate and number of iterations
learning_rate = 0.01;
num_iterations = 1000000;

% Gradient descent
for iter = 1:num_iterations
    % Calculate the current model response
    R_model = A * exp(-stimulus.^2 / (2 * sigma^2)) + offset;

    % Compute the loss (Mean Squared Error)
    loss = mean((R_model - tuning_data).^2);

    % Compute gradients w.r.t parameters
    grad_A = mean(2 * (R_model - tuning_data) .* exp(-stimulus.^2 / (2 * sigma^2)));
    grad_sigma = mean(2 * (R_model - tuning_data) .* ...
                      A .* exp(-stimulus.^2 / (2 * sigma^2)) .* (stimulus.^2 / sigma^3));
    grad_offset = mean(2 * (R_model - tuning_data));

    % Update parameters
    A = A - learning_rate * grad_A;
    sigma = sigma - learning_rate * grad_sigma;
    offset = offset - learning_rate * grad_offset;

    % Print loss every 1000 iterations for monitoring
%     if mod(iter, 1000) == 0
%         fprintf('Iteration %d, Loss: %.4f\n', iter, loss);
%     end
end
% Plot the data and the fitted model
% figure;
% plot(stimulus, tuning_data, 'o', 'Color', color);
hold on;
stimulus_n = linspace(stimulus(1), stimulus(end), 1000);
gauss_fit = A * exp(-stimulus_n.^2 / (2 * sigma^2)) + offset;
% gauss_fit = gauss_fit + max(tuning_data)-max(gauss_fit);
plot(stimulus_n, gauss_fit, '-', ...
     'Color', color);
% legend;
% xlabel('Stimulus');
% ylabel('Response');
% title('Gaussian Fit with Offset');
end