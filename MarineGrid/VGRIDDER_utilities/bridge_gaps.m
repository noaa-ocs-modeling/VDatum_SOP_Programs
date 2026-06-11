function bridges_found=bridge_gaps(grid)

[rows, cols] = size(grid);
% 1. Number groups of connected 1s
groups = grlabel(grid);

% 2. Find all internal '0' tiles (skipping the 1-pixel outer border)
[R, C] = find(grid == 0 & ...
              (1:rows)' > 1 & (1:rows)' < rows & ...
              (1:cols) > 1 & (1:cols) < cols);

nR = length(R);
distinct_group_counts = zeros(nR, 1);

% 3. Calculate how many UNIQUE groups surround each candidate '0'
for i = 1:nR
    r = R(i);
    c = C(i);

    % Extract the 3x3 neighborhood of group IDs
    neighborhood = groups(r-1:r+1, c-1:c+1);

    % Find unique IDs and remove 0 (which represents empty black tiles)
    unique_groups = unique(neighborhood);
    unique_groups(unique_groups == 0) = []; 

    % Store the count of distinct groups touching this tile
    distinct_group_counts(i) = length(unique_groups);
end

% 4. Filter out tiles touching fewer than 2 distinct groups 
% (If it touches 0 or 1 group, it cannot possibly be a bridge)
valid_idx = find(distinct_group_counts >= 2);
R = R(valid_idx);
C = C(valid_idx);
distinct_group_counts = distinct_group_counts(valid_idx);

% 5. Sort candidates in descending order of distinct groups touched
[~, sort_order] = sort(distinct_group_counts, 'descend');
sorted_R = R(sort_order);
sorted_C = C(sort_order);
nR = length(R);

bridges_found = zeros(nR, 2); 
k = 0; % Counter for found bridges

for i=1:nR
    r=sorted_R(i);
    c=sorted_C(i);
    
    % Get the 4 orthogonal neighbor group IDs
    i_up = groups(r-1, c); 
    i_dn = groups(r+1, c); 
    i_rt = groups(r, c+1); 
    i_lf = groups(r, c-1); 
    % Get the 4 diagogonal neighbor group IDs
    i_uprt = groups(r-1, c+1); 
    i_uplf = groups(r-1, c-1); 
    i_dnrt = groups(r+1, c+1); 
    i_dnlf = groups(r+1, c-1); 
   
    % Store all 6 possible pairs of sides
    gr = [i_up, i_dn;  % Pair 1  orthogonal to orthogonal
          i_up, i_rt;  % Pair 2
          i_up, i_lf;  % Pair 3
          i_dn, i_rt;  % Pair 4
          i_dn, i_lf;  % Pair 5
          i_rt, i_lf;  % Pair 6
          i_up, i_dnrt; % Pair 7  orthogonal to diagonal
          i_up, i_dnlf; % Pair 8
          i_dn, i_uprt; % Pair 9
          i_dn, i_uplf; % Pair 10
          i_rt, i_uplf; % Pair 11
          i_rt, i_dnlf; % Pair 12
          i_lf, i_uprt; % Pair 13
          i_lf, i_dnrt]; % Pair 14
    
    gr_rows = 1:14;
    
    % Find rows where BOTH elements are non-zero AND not equal
    valid_pairs = (gr(:,1) ~= 0) & (gr(:,2) ~= 0) & (gr(:,1) ~= gr(:,2));    
    gr = gr(valid_pairs, :);
    gr_rows = gr_rows(valid_pairs);
    ngr = length(gr_rows);
    
    % If no sides can be bridged, move to the next candidate 0
    if ngr == 0 
        continue;
    end
    
    % --- Valid Bridge Found! ---
    grid(r, c) = 1; % Flip the tile to 1
    k = k + 1;
    bridges_found(k, 1:2) = [r, c];
    % update the groups matrix
    ia=gr(1,1);
    groups(r,c)=ia;
    %leave only orthogonal-orthogonal connections for group merging
    gr_rows(gr_rows>6)=[];
    ngr = length(gr_rows);    
    if ngr == 0 
        continue;
    end
    for j=1:ngr
        switch gr_rows(j)
            case 1
                ib = i_dn;
            case {2, 4}
                ib = i_rt;
            case {3, 5, 6}
                ib = i_lf;
        end
        % Merge group ib into ia
        if ia ~= ib
            groups(groups == ib) = ia;
        end
    end
end

% 2. Clean up the pre-allocated array by chopping off the unused rows
bridges_found = bridges_found(1:k, :);

end

function labels = grlabel(grid)
    % BWLABEL labels connected components in a binary matrix (4-connectivity).
    
    [rows, cols] = size(grid);
    labels = zeros(rows, cols);
    current_id = 0;
    
    % Direction offsets for 4-connectivity: Up, Down, Left, Right
    dr = [-1, 1, 0, 0];
    dc = [ 0, 0, -1, 1];
    
    for r = 1:rows
        for c = 1:cols
            % If we find an unlabeled '1', start a new group
            if grid(r, c) == 1 && labels(r, c) == 0
                current_id = current_id + 1;
                
                % Initialize a manual Queue for BFS
                % Maximum possible queue size is total number of elements
                q_r = zeros(rows * cols, 1);
                q_c = zeros(rows * cols, 1);
                head = 1;
                tail = 1;
                
                % Push starting element onto queue
                q_r(tail) = r;
                q_c(tail) = c;
                labels(r, c) = current_id;
                
                % Loop until the queue is empty
                while head <= tail
                    % Pop current coordinates
                    curr_r = q_r(head);
                    curr_c = q_c(head);
                    head = head + 1;
                    
                    % Check all 4 orthogonal neighbors
                    for i = 1:4
                        nr = curr_r + dr(i);
                        nc = curr_c + dc(i);
                        
                        % Ensure neighbor is within matrix bounds
                        if nr >= 1 && nr <= rows && nc >= 1 && nc <= cols
                            % If it's a '1' and hasn't been labeled yet
                            if grid(nr, nc) == 1 && labels(nr, nc) == 0
                                tail = tail + 1;
                                q_r(tail) = nr;
                                q_c(tail) = nc;
                                labels(nr, nc) = current_id;
                            end
                        end
                    end
                end
                
            end
        end
    end
end
