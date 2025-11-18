function q = fdr_bh(p)
    % Benjamini-Hochberg

    [p_sorted, sort_idx] = sort(p);
    n = length(p);
    q_sorted = p_sorted .* n ./ (1:n)';  
    q_sorted = min(1, q_sorted);  
    for i = n-1:-1:1
        q_sorted(i) = min(q_sorted(i), q_sorted(i+1));
    end

    q = zeros(size(p));
    q(sort_idx) = q_sorted;
end