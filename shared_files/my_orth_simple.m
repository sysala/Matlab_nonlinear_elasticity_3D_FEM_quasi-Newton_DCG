function [W_orth] = my_orth_simple(V, W, tol, maxsize)

    %--------------------------------------------------------------------------
    % This function orthonormalizes the vectors in W with respect to the V and
    % also to itself.
    %
    %  input data:
    %       V      - matrix of vectors to be orthogonalized against, size(V)=(n_p,n)
    %       W      - matrix of vectors to be orthogonalized, size(W)=(n_p,n)
    %       tol    - tolerance for the orthogonality, default value is 1e-6
    %       maxsize- maximum size of the matrix V, default value is 1000
    %
    %  output data:
    %       W_orth - matrix of orthogonalized vectors, size(W_orth)=(n_p,n)
    %--------------------------------------------------------------------------

    n = size(W, 2);

    W_orth = [];

    if ~exist('maxsize', 'var')
        % third parameter does not exist, so default it to something
        maxsize = 1000;
    end

    if ~exist('tol', 'var')
        % third parameter does not exist, so default it to something
        tol = 1e-6;
    end

    if size(V, 2) >= maxsize
        return
    end

    for i = 1:n
        temp = W(:, i) / norm(W(:, i));

        if ~isempty(V)
            temp = temp - V * (temp' * V)';
            temp = temp - V * (temp' * V)';
        end

        if ~isempty(W_orth)

            for j = 1:size(W_orth, 2)
                temp = temp - dot(temp, W_orth(:, j)) * W_orth(:, j);
            end

            for j = 1:size(W_orth, 2)
                temp = temp - dot(temp, W_orth(:, j)) * W_orth(:, j);
            end

        end

        n_temp = norm(temp);

        if n_temp > tol
            W_orth = [W_orth temp / n_temp];
        end

    end

end
