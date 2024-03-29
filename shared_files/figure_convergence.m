function figure_convergence(X1, Y1, X2, Y2, X3, Y3)
    %CREATEFIGURE(X1, Y1, X2, Y2, X3, Y3)
    %  X1:  vector of semilogy x data
    %  Y1:  vector of semilogy y data
    %  X2:  vector of semilogy x data
    %  Y2:  vector of semilogy y data
    %  X3:  vector of semilogy x data
    %  Y3:  vector of semilogy y data

    %  Auto-generated by MATLAB on 26-May-2023 14:03:56

    % Create figure
    figure1 = figure;

    % Create axes
    axes1 = axes('Parent', figure1);
    hold(axes1, 'on');

    % Create semilogy
    semilogy(X1, Y1, 'DisplayName', 'Newton', 'Marker', 'o', 'LineWidth', 1, ...
        'Color', [0 0 0]);

    % Create semilogy
    semilogy(X2, Y2, 'DisplayName', 'q-Newton 1', 'Marker', 'square', 'LineWidth', 1, ...
        'Color', [0 0 1]);

    % Create semilogy
    semilogy(X3, Y3, 'DisplayName', 'q-Newton 2', 'Marker', 'diamond', 'LineWidth', 1, ...
        'Color', [1 0 0]);

    % Create ylabel
    ylabel('||F(u_n)|| / ||F(0)||');

    % Create xlabel
    xlabel('iteration');

    % Create title
    %title('Model 1');

    box(axes1, 'on');
    hold(axes1, 'off');
    % Set the remaining axes properties
    set(axes1, 'FontSize', 24, 'XGrid', 'on', 'YGrid', 'on', 'YMinorTick', 'on', ...
        'YScale', 'log');
    % Create legend
    legend1 = legend(axes1, 'show');
    set(legend1, ...
        'Position', [0.788194446074258 0.73314044951611 0.103645831036071 0.156380748649022]);
