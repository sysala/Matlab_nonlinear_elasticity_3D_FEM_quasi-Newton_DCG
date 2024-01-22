function draw_mesh(coord, surf, Q_weak)

    % =========================================================================
    %
    %  This function draws a mesh reflecting the prescribed heterogeneity.
    %  Black color is used for the coal while the yellow one for the PUR.
    %
    %  input data:
    %    coord - coordinates of the nodes, size(coord)=(3,n_n) where n_n is a
    %            number of nodes
    %    surf  - array containing numbers of nodes defining each surface element,
    %            size(surf)=(n_p,n_s), n_s = number of surface elements
    %    elem_type - the type of finite elements; available choices:
    %                'P1', 'P2', 'Q1', 'Q2'
    %
    % ======================================================================
    %

    figure
    hold on
    patch('Faces', surf(1:3, Q_weak)', 'Vertices', coord', 'FaceVertexCData', ...
        0 * ones(size(coord, 2), 1), 'FaceColor', 'yellow', 'EdgeColor', 'none');
    patch('Faces', surf(1:3, ~Q_weak)', 'Vertices', coord', 'FaceVertexCData', ...
        0 * ones(size(coord, 2), 1), 'FaceColor', 'black', 'EdgeColor', 'none');

    %   ind=unique(surf(:));
    %   plot3( coord(1,ind),coord(2,ind),coord(3,ind), 'b.', 'MarkerSize',10);
    axis equal; % real ratios
    view(3); % standard view ve 3D
    hold off;
    axis off;
end
