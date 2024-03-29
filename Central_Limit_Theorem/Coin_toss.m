central_limit_theorem_coins(100, 10000);

function central_limit_theorem_coins(N, num_tosses)
    % Simulate tossing N fair coins num_tosses times
    tosses = randi([0, 1], num_tosses, N); % 0 represents tails, 1 represents heads

    % Calculate the sum of each set of N coin tosses
    sums = sum(tosses, 2);

    % Plot the histogram
    figure;
    histogram(sums, 'Normalization', 'probability', 'BinWidth', 1, 'EdgeColor', 'w');

    % Set plot labels and title
    title(sprintf('Central Limit Theorem for Tossing %d Fair Coins %d times', N, num_tosses));
    xlabel('Heads');
    ylabel('Probability');
end
