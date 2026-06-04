clc; clear; close all;

% initial positionx x0,y0,z0 of 21 entries each. The 21st entry is the
% position of the gap entry point

x0 = [37.3042355570610	17.4414337604173	33.5563940064744	7.06139153417222	18.0838304844253	41.1989169323702	44.2500031849004	52.5818889566870	56.8728136701576	51.0869161091717	21.2464354638316	8.02298395028173	53.4296029344811	38.4036428783863	59.8517779515722	27.7527079835487	1.90926705790630	51.6764048506426	31.1022650893168	23.9453845395823	30];
y0 = [41.1549850771227	15.7507253783850	48.8358056819503	59.2732494763547	33.3517783240820	11.1271264111924	45.4084173365506	21.0265541847184	13.1847763483154	54.2873463832546	47.1377113273386	9.53084435527321	2.18563816019498	53.8993880816862	25.2626534905591	37.6340635848952	32.8388557731536	37.9959763121511	6.99101609967396	58.2724486378843	65];
z0 = [4.32751150642279	3.21334574964582	4.64860581491173	7.98243089236830	7.27074428922344	4.87360664131823	7.24851562569565	7.06344977681379	6.59453880849495	6.57257425402966	5.59752263070131	3.22767919205823	4.31500673639087	7.62459371327626	7.80334224142276	3.92350787915914	5.82412492815636	7.73729464962396	5.46084763419448	6.80821135668840	5];

% Cost parameters
Gamma_a = 100; Gamma_p = 1e4; sf = 50; gamma_crit = 20*pi/180; dmin = 1;

% greedy algorithm
greedy_parent = hbt_parent(x0, y0, z0);

total_cost_greedy = compute_path_cost(greedy_parent, x0, y0, z0, Gamma_a, Gamma_p, sf, gamma_crit, dmin);
[best_parent, best_cost_ga1, cost_history] = hbt_ga(x0, y0, z0, greedy_parent, Gamma_a, Gamma_p, sf, gamma_crit, dmin);


[~,root_node] = max(y0);


figure()


        hold on
        % axis equal
        ParentNode_value =best_parent;
        for i =1:length(x0)-1
                plot3([x0(i),x0(ParentNode_value(i))],[y0(i),y0(ParentNode_value(i))],[z0(i),z0(ParentNode_value(i))],'Color','b','LineWidth',2)
            
   
        end
        plot3(x0,y0,z0,'o','MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',8)
        plot3(x0(root_node),y0(root_node),z0(root_node),'o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',8)

drawTunnel([x0(root_node),y0(root_node)+1,z0(root_node)])


% =========================================================================

% GREEDY HBT


function parent = hbt_parent(x, y, z)
% HBT_PARENT  Deterministic greedy hierarchical binary tree construction.

N = length(x);
root = find(y == max(y), 1);
parent = zeros(1, N);
assigned = false(1, N);
assigned(root) = true;
candidates = cell(1, N);
candidates{root} = find(~assigned);
queue = root;

while ~isempty(queue)
    k = queue(1); queue(1) = [];
    Q = candidates{k};
    if isempty(Q), continue; end

    dx = x(Q) - x(k); dy = y(Q) - y(k); dz = z(Q) - z(k);

    if numel(Q) == 1 || all_collinear(dx, dy, dz)
        children = Q(1); Lset = Q; Rset = [];
    else
     
        M = [dx(:), dy(:), dz(:)];      
        [~, ~, V] = svd(M, 0);
        n = V(:, 3);                    

        [~, imax] = max(abs(n));
        if n(imax) < 0, n = -n; end

        proj = n(1)*dx + n(2)*dy + n(3)*dz;

        tol = 1e-9 * max(1, max(abs(M(:))));
        Lmask = proj < -tol | (abs(proj) <= tol);   % on-plane -> Left
        Lset = Q(Lmask); Rset = Q(~Lmask);

        if isempty(Lset) || isempty(Rset)
            [~, idx] = min(dx.^2 + dy.^2 + dz.^2);
            children = Q(idx); Lset = Q; Rset = [];
        else
            kL = nearest(k, Lset, x, y, z);
            kR = nearest(k, Rset, x, y, z);
            children = [kL, kR];
        end
    end

    for c = children
        parent(c) = k; assigned(c) = true;
        if ismember(c, Lset)
            cset = Lset(Lset ~= c);
        else
            cset = Rset(Rset ~= c);
        end
        candidates{c} = cset(~assigned(cset));
        queue(end+1) = c; 
    end
end
end

function idx = nearest(k, S, x, y, z)
    d2 = (x(S)-x(k)).^2 + (y(S)-y(k)).^2 + (z(S)-z(k)).^2;
    [~, i] = min(d2); idx = S(i);
end

function tf = all_collinear(dx, dy, dz)
    if numel(dx) < 2, tf = true; return; end
    v = [dx(1); dy(1); dz(1)];
    if norm(v) < 1e-10, tf = true; return; end
    v = v / norm(v);
    dots = abs(dx*v(1) + dy*v(2) + dz*v(3));
    norms = sqrt(dx.^2 + dy.^2 + dz.^2);
    tf = all(abs(dots - norms) < 1e-8);
end

% =========================================================================
% COST FUNCTION

function [total_cost, node_cost] = compute_path_cost(parent, x, y, z, Gamma_a, Gamma_p, sf, gamma_crit, dmin)
N    = length(parent);
root = find(parent == 0, 1);

% Validity check
order = bfs_order(parent, root);
if length(order) ~= N
    total_cost = inf; node_cost = inf(1,N); return;
end

% Subtree sizes (n_i = number of UAVs traversing branch (i, parent(i)))
n_i = ones(1, N);
for idx = length(order):-1:1
    k = order(idx);
    if parent(k) ~= 0
        n_i(parent(k)) = n_i(parent(k)) + n_i(k);
    end
end

% Children list
ch_list = cell(1, N);
for i = 1:N
    if parent(i) ~= 0
        ch_list{parent(i)}(end+1) = i;
    end
end

% Node costs in BFS order
node_cost = zeros(1, N);
for idx = 1:length(order)
    i = order(idx);
    if i == root, continue; end
    pi_i     = parent(i);
    edge_len = norm([x(i)-x(pi_i), y(i)-y(pi_i), z(i)-z(pi_i)]);
    gamma_i  = compute_gamma(i, pi_i, parent, ch_list, x, y, z);
    angle_pen = (Gamma_a * n_i(i)) / (1 + exp(sf * (gamma_i - gamma_crit)^3) );
    node_cost(i) = node_cost(pi_i) + edge_len + angle_pen;
end

sep_total  = compute_all_separation_penalties(parent, x, y, z, dmin, Gamma_p);
total_cost = sum(node_cost) + sep_total;
end

function gamma_i = compute_gamma(i, pi_i, parent, ch_list, x, y, z)
v_i     = [x(i)-x(pi_i); y(i)-y(pi_i); z(i)-z(pi_i)];
pi_pi_i = parent(pi_i);
siblings = ch_list{pi_i};
sibling  = siblings(siblings ~= i);

if pi_pi_i == 0
    if ~isempty(sibling)
        j = sibling(1);
        v_j = [x(j)-x(pi_i); y(j)-y(pi_i); z(j)-z(pi_i)];
        gamma_i = angle_between(v_i, v_j);
    else
        gamma_i = pi;
    end
    return;
end

v_gp     = [x(pi_pi_i)-x(pi_i); y(pi_pi_i)-y(pi_i); z(pi_pi_i)-z(pi_i)];
angle_gp = angle_between(v_i, v_gp);
if ~isempty(sibling)
    j = sibling(1);
    v_j      = [x(j)-x(pi_i); y(j)-y(pi_i); z(j)-z(pi_i)];
    angle_sib = angle_between(v_i, v_j);
    gamma_i  = min(angle_gp, angle_sib);
else
    gamma_i = angle_gp;
end
end

function a = angle_between(u, v)
    a = acos(max(-1, min(1, dot(u,v) / (norm(u)*norm(v)))));
end

% Separation  penalty

function total_sep_penalty = compute_all_separation_penalties(parent, x, y, z, dmin, Gamma_p)
edges = find(parent ~= 0);
M     = length(edges);
total_sep_penalty = 0;
P1 = [x(edges);         y(edges);         z(edges)        ];
P2 = [x(parent(edges)); y(parent(edges)); z(parent(edges)) ];

for a = 1:M
    i    = edges(a); pi_i = parent(i);
    p1   = P1(:,a);  p2   = P2(:,a);
    for b = a+1:M
        k = edges(b); pk_i = parent(k);
        if pk_i == pi_i, continue; end
        if pk_i == i,    continue; end
        if k    == pi_i, continue; end
        d = seg_seg_dist(p1, p2, P1(:,b), P2(:,b));
        if d <= dmin
            total_sep_penalty = total_sep_penalty + Gamma_p;
        end
    end
end
end

function d = seg_seg_dist(p1, p2, p3, p4)
d1 = p2-p1; d2 = p4-p3; r = p1-p3;
a = dot(d1,d1); e = dot(d2,d2); f = dot(d2,r);
if a<1e-10 && e<1e-10, d=norm(r); return; end
if a<1e-10
    s=0; t=clamp01(f/e);
else
    c=dot(d1,r);
    if e<1e-10, t=0; s=clamp01(-c/a);
    else
        b=dot(d1,d2); denom=a*e-b*b;
        if denom>1e-10, s=clamp01((b*f-c*e)/denom); else s=0; end
        t=(b*s+f)/e;
        if t<0, t=0; s=clamp01(-c/a);
        elseif t>1, t=1; s=clamp01((b-c)/a); end
    end
end
d = norm(p1+s*d1-p3-t*d2);
end

function v = clamp01(x), v = max(0, min(1,x)); end

function order = bfs_order(parent, root)
N     = length(parent);
order = zeros(1,N);
head  = 1; tail = 1;
order(1) = root;
while head <= tail
    k = order(head); head = head+1;
    ch = find(parent == k);
    for c = ch
        tail = tail+1; order(tail) = c;
    end
end
order = order(1:tail);
end

% =========================================================================
% GA


function [best_parent, best_cost, cost_history] = hbt_ga(x, y, z, greedy_parent, Gamma_a, Gamma_p, sf, gamma_crit, dmin)
N        = length(x);
root     = find(y == max(y), 1);
M_pop    = 50;
max_it    = 10000;
k_tourn  = 2;
i_max = 500;

% --- Initialization
pop = cell(1, M_pop);
pop{1} = greedy_parent;
for m = 2:M_pop
    pop{m} = random_hbt(N, root);
end

% --- Initial fitness
fitness = inf(1, M_pop);
for m = 1:M_pop
    fitness(m) = compute_path_cost(pop{m}, x, y, z, Gamma_a, Gamma_p, sf, gamma_crit, dmin);
end

[best_cost, bi] = min(fitness);
best_parent  = pop{bi};
cost_history = nan(1, max_it);
stag_count   = 0;

% --- Main loop
for gen = 1:max_it
    % Selection
    idx1 = tournament_select(fitness, k_tourn);
    idx2 = tournament_select(fitness, k_tourn);
    while idx2 == idx1
        idx2 = tournament_select(fitness, k_tourn);
    end

    % Crossover — produces two guaranteed-valid offspring
    [o1, o2] = crossover(pop{idx1}, pop{idx2}, root, N);

    % Mutation
    o1 = mutate(o1, root, N);
    o2 = mutate(o2, root, N);

    % Evaluate
    f1 = compute_path_cost(o1, x, y, z, Gamma_a, Gamma_p, sf, gamma_crit, dmin);
    f2 = compute_path_cost(o2, x, y, z, Gamma_a, Gamma_p, sf, gamma_crit, dmin);

    % Elitist replacement
    pop     = [pop,    {o1}, {o2}];
    fitness = [fitness, f1,   f2 ];
    [fitness, sidx] = sort(fitness);
    pop     = pop(sidx);
    pop     = pop(1:M_pop);
    fitness = fitness(1:M_pop);

    % Track best
    cost_history(gen) = fitness(1);
    if fitness(1) < best_cost
        best_cost   = fitness(1);
        best_parent = pop{1};
        stag_count  = 0;
    else
        stag_count = stag_count + 1;
    end

    if stag_count >= i_max
        fprintf('Converged at generation %d | Best cost: %.4f\n', gen, best_cost);
        cost_history = cost_history(1:gen);
        return;
    end

    if mod(gen, 100) == 0
        fprintf('Gen %d | Best cost: %.4f\n', gen, best_cost);
    end
end
end


function idx = tournament_select(fitness, k)
    cands  = randperm(length(fitness), k);
    [~, i] = min(fitness(cands));
    idx    = cands(i);
end

% =========================================================================
% crossover

function [o1, o2] = crossover(chrA, chrB, root, N)

cx      = pick_crossover_node(chrA, chrB, root, N);
subA    = get_subtree(chrA, cx, N);   
subB    = get_subtree(chrB, cx, N);   



floatA  = setdiff(subA, subB);       
o1      = build_offspring(chrA, chrB, cx, subB, floatA, root, N);


floatB  = setdiff(subB, subA);
o2      = build_offspring(chrB, chrA, cx, subA, floatB, root, N);
end

function cx = pick_crossover_node(chrA, chrB, root, N)

non_root = setdiff(1:N, root);

order = non_root(randperm(length(non_root)));
cx = order(1);  
for i = 1:length(order)
    cand   = order(i);
    subA_i = get_subtree(chrA, cand, N);
    subB_i = get_subtree(chrB, cand, N);
    if ~isequal(sort(subA_i), sort(subB_i))
        cx = cand;
        return;
    end
end
end

function offspring = build_offspring(base, donor, cx, donor_sub, float_nodes, root, N)

offspring = base;


for n = donor_sub
    offspring(n) = donor(n);
end


offspring(cx) = donor(cx);


for n = float_nodes
    offspring(n) = 0;  
end



offspring = reattach_nodes(offspring, float_nodes, root, N);



offspring = enforce_degree(offspring, root, N);
end

function parent = reattach_nodes(parent, nodes, root, N)


nodes = nodes(randperm(length(nodes)));

for d = nodes

    reachable = get_reachable(parent, root, N);

    deg = compute_degree(parent, N);


    valid = find(reachable & deg < 2 & (1:N) ~= d);

    if isempty(valid)
        

        parent(d) = root;
    else
        parent(d) = valid(randi(length(valid)));
    end
end
end

function reachable = get_reachable(parent, root, N)
% only following nodes with valid parent links
reachable = false(1, N);
reachable(root) = true;
queue = root;
while ~isempty(queue)
    k = queue(1); queue(1) = [];
    ch = find(parent == k);
    for c = ch
        if ~reachable(c)
            reachable(c) = true;
            queue(end+1) = c; 
        end
    end
end
end

function deg = compute_degree(parent, N)
% Number of children of each node
deg = zeros(1, N);
for i = 1:N
    if parent(i) ~= 0
        deg(parent(i)) = deg(parent(i)) + 1;
    end
end
end

function parent = enforce_degree(parent, root, N)
% If any node has > 2 children, detach excess and reattach
changed = true;
while changed
    changed = false;
    for k = 1:N
        ch = find(parent == k);
        if length(ch) > 2
            excess  = ch(3:end);
            parent(excess) = 0;
            parent  = reattach_nodes(parent, excess, root, N);
            changed = true;
        end
    end
end

% Check for any remaining detached non-root nodes
detached = find(parent == 0 & (1:N) ~= root);
if ~isempty(detached)
    parent = reattach_nodes(parent, detached, root, N);
end
end

% =========================================================================
% MUTATION: parent-swap as described in the paper]

function parent = mutate(parent, root, N)
non_root   = setdiff(1:N, root);
% Candidates: non-root nodes that have at least one child
has_child  = non_root(arrayfun(@(k) any(parent == k), non_root));
if isempty(has_child), return; end

% Shuffle candidates and find first safe mutation
order = has_child(randperm(length(has_child)));
for trial = 1:length(order)
    k  = order(trial);
    ch = find(parent == k);
    kL = ch(1);     % pick first child

    pi_k = parent(k);

    % Perform swap on a copy first, validate, then commit
    p_new      = parent;
    p_new(k)   = kL;      
    p_new(kL)  = pi_k;   

    % If k had second child kR, kL adopts it
    if length(ch) > 1
        kR         = ch(2);
        p_new(kR)  = kL;
    end

    % Validate: check kL doesn't now have > 2 children
    kL_children = find(p_new == kL);
    if length(kL_children) <= 2

        if pi_k ~= 0
            pik_children = find(p_new == pi_k);
            if length(pik_children) > 2
                continue;  
            end
        end
        
        parent = p_new;
        return;
    end
  
end
% No valid mutation found — return unchanged
end


function sub = get_subtree(parent, node, N)
sub   = false(1, N);
sub(node) = true;
queue = node;
while ~isempty(queue)
    k     = queue(1); queue(1) = [];
    ch    = find(parent == k);
    for c = ch
        if ~sub(c)          
            sub(c) = true;
            queue(end+1) = c;
        end
    end
end
sub = find(sub);
end

% =========================================================================
% RANDOM VALID HBT GENERATION

function parent = random_hbt(N, root)
parent   = zeros(1, N);
deg      = zeros(1, N);
assigned = false(1, N);
assigned(root) = true;

% Random insertion order: root first, then all others shuffled
others = setdiff(1:N, root);
order  = [root, others(randperm(length(others)))];

for i = 2:N
    node  = order(i);
    valid = find(assigned & deg < 2);
    par   = valid(randi(length(valid)));
    parent(node)   = par;
    deg(par)       = deg(par) + 1;
    assigned(node) = true;
end
end



function drawTunnel(center)
    
    
    % Tunnel dimensions
    height = 2; % Height in meters
    width = 2;  % Width in meters
    depth = 10; % Depth in meters
    
    % Extract center coordinates
    x0 = center(1);
    y0 = center(2);
    z0 = center(3);
    
    % Define the 8 vertices of the tunnel
    % Front face vertices (starting at bottom-left and moving clockwise)
    front =[x0,y0,z0] + (rotx(90)* [
        - width/2,  - height/2, 0; % Bottom-left
         + width/2,  - height/2, 0; % Bottom-right
         + width/2,  + height/2, 0; % Top-right
        - width/2, + height/2, 0; % Top-left
    ]')';
    
    % Back face vertices (shifted by depth along x-axis)
    back = [x0,y0,z0] + (rotx(90)*[
        - width/2, - height/2,  - depth; % Bottom-left
        + width/2, - height/2, - depth; % Bottom-right
        + width/2, + height/2,  - depth; % Top-right
       - width/2,  + height/2,  - depth; % Top-left
    ]')';
    
    % Combine all vertices
    vertices = [front; back];
    
    % Define the faces of the tunnel using vertex indices
    faces = [
        1 2 3 4; % Front face
        5 6 7 8; % Back face
        1 2 6 5; % Bottom face
        4 3 7 8; % Top face
        1 5 8 4; % Left face
        2 6 7 3; % Right face
    ];
    
    % Plot the tunnel
   
    hold on;
    patch('Vertices', vertices, 'Faces', faces, ...
          'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'k', 'FaceAlpha', 0.8);
    
    % Plot the center point for reference
    % plot3(x0, y0, z0, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    % text(x0, y0, z0, ' Center', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    
    % Enhance the visualization
    xlabel('$$x,$$ m','Interpreter','latex');
    ylabel('$$y,$$ m','Interpreter','latex');
    zlabel('$$z,$$ m','Interpreter','latex');

    axis equal;
    grid on;
    view(3); % 3D view
    hold off;
end


