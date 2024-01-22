function assembly_comparison(X1, Y1, X2, Y2)
    %CREATEFIGURE(X1, Y1, X2, Y2)
    %  X1:  vector of plot x data
    %  Y1:  vector of plot y data
    %  X2:  vector of plot x data
    %  Y2:  vector of plot y data

    %  Auto-generated by MATLAB on 04-May-2023 08:48:38

    % Create figure
    figure1 = figure;

    % Create axes
    axes1 = axes('Parent', figure1);
    hold(axes1, 'on');

    % Create plot
    plot(X1, Y1, 'DisplayName', 'Newton', 'Marker', 'o', 'LineWidth', 1, ...
        'Color', [0 0 0]);

    % Create plot
    plot(X2, Y2, 'DisplayName', 'q-Newton 1', 'Marker', 'square', 'LineWidth', 1, ...
        'Color', [0 0 1]);

    % Create ylabel
    ylabel('assembly time');

    % Create xlabel
    xlabel('iteration');

    % Create title
    title('mesh with 129600 elements');

    box(axes1, 'on');
    hold(axes1, 'off');
    % Set the remaining axes properties
    set(axes1, 'FontSize', 24);
    % Create legend
    legend(axes1, 'show');
