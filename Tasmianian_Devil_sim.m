function [positions,population,infect_percentage] = Tasmianian_Devil_sim(N, init_pos)
%       The Random_walk_cell_sim function uses a random walk given the 
%       inputs to derive and output the position of M nutrient cells for
%       N iterations. This function also outputs a tracking array which
%       shows the amount of nutrient cells inside the cell Membrane. 
    
%   Initialise variables
    Simulation_vid = VideoWriter('Simulation_vid');
    Simulation_vid.FrameRate = 10;
    open(Simulation_vid);

    K = length(init_pos(:,1));
    M = length(init_pos(1,:));
    positions = zeros(K, M, N);
    positions(:, :,1) = init_pos;
    population = zeros(N,1);
    population(1) = sum(sum(init_pos >= 1));
    infect_percentage = zeros(N,1);
    infect_percentage(1) = sum(sum(init_pos > 1))/population(1);
%   Visualise Simulation    
    fig = figure('visible','on');
    handle = imshow(positions(:,:,1), 'InitialMagnification', 'Fit');
    colormap(parula)

for dt = 2:N
    weeks = mod(dt,52);
    mating_season = (weeks >= 8 && weeks <= 20);

    infect_status = positions(:,:,dt-1);
    new_positions = zeros(K,M);  % Create a fresh layer for this timestep

    for i = 1:K
        for j = 1:M
            if positions(i,j,dt-1) >= 1
                A = Boundary_Conditions(positions,i,j,K,M,dt-1);
                A2 = Boundary_Conditions(positions,i,j,K,M,dt);
                [weighting_x, weighting_y, movement_type_weighting, infect_weight, cub_status, breeding_weighting, backup_step] = Cellular_autonomy_direction(A, A2, mating_season);

                [dx, dy] = Random_walk_two_dim(positions(:,:,dt-1), j, i, weighting_x, weighting_y, movement_type_weighting);
                new_i = i + dy;
                new_j = j + dx;

                % Check bounds
                if new_i < 1 || new_i > K || new_j < 1 || new_j > M || new_positions(new_i, new_j) >= 1
                    % Try backup directions
                    moved = false;
                    for f = 1:8
                        dy_try = backup_step(f,1);
                        dx_try = backup_step(f,2);
                        ni = i + dy_try;
                        nj = j + dx_try;
                        if ni >= 1 && ni <= K && nj >= 1 && nj <= M && new_positions(ni, nj) == 0
                            new_i = ni;
                            new_j = nj;
                            moved = true;
                            break;
                        end
                    end
                    if ~moved
                        % Stay in place
                        new_i = i;
                        new_j = j;
                    end
                end

                new_positions(new_i, new_j) = update_infection_status(infect_status(i,j), infect_weight);

                 % Breeding logic
                   if cub_status == 1
                    [B,y,x] = Boundary_Conditions(positions,new_i,new_j,K,M,dt);
                   for r = 1:length(y)
                        for c = 1:length(x)
                            row = y(r);
                            col = x(c);
                            if new_positions(row, col) == 0
                                new_positions(row, col) = breeding_outcome(cub_status, breeding_weighting);
                           end
                        end
                    end
                end
            end
        end
    end

    % Assign updated layer
    positions(:,:,dt) = new_positions;
    handle.CData = new_positions;
    drawnow;
    writeVideo(Simulation_vid, getframe(fig));

    % Logging
    population(dt) = sum(new_positions(:) >= 1);
    infect_percentage(dt) = sum(new_positions(:) > 1) / max(1, population(dt));
    fprintf("Week %d - Living: %d, Infected: %.2f%%\n", dt, population(dt), 100*infect_percentage(dt));
end

end

% -------------------------------------------------------------------------

function [updated_infect_status] = update_infection_status(infect_status,infect_weight)
%   The update_infection_status function takes the current infection status 
%   of a Tasmanian Devil and the infection weighting and determines what 
%   happens to that Tasmanian Devil's infection status at the next time 
%   iteration.

    probability_infection = abs(randn(1));
    if round(probability_infection,0) > 1 
        probability_infection = probability_infection - round(probability_infection,0);
    end

    if  infect_status == 1 && infect_weight ~= 0
        if infect_status == 1 && probability_infection < infect_weight
            updated_infect_status = 2;  
        else
            updated_infect_status = infect_status;
        end
    elseif infect_status >= 2 % Infected, Greater infect_status closer to death until infect_stat == 7
        if infect_status == 27
            updated_infect_status = 0; % Dead tasmanian devil
        else
            updated_infect_status = infect_status + 1;  
        end
    elseif infect_status == 1 && infect_weight == 0   
        updated_infect_status = infect_status;
    end
end

% -------------------------------------------------------------------------

function [Newborn_Status] = breeding_outcome(Cub_status,Cub_weighting)
%   The breeding_outcome function takes the current infection status of a
%   Tasmanian Devil and the infection weighting and determines what happens
%   to that Tasmanian Devil's infection status at the next time iteration.

    probability_infection = abs(randn(1));
    if round(probability_infection,0) > 1 
        probability_infection = probability_infection - round(probability_infection,0);
    end
    
    conception_probability = abs(randn(1));
    if round(conception_probability,0) > 1 
        conception_probability = conception_probability - round(conception_probability,0);
    end

    if  Cub_status == 1 && Cub_weighting ~= 0
        if conception_probability > Cub_weighting % Baby not born
            Newborn_Status = 0;
        elseif Cub_status == 1 && conception_probability <= Cub_weighting % Baby born without infection
            Newborn_Status = 1;
        end
    else
            Newborn_Status = 0; % No New tasmanian devil
    end
end

% -------------------------------------------------------------------------

function [A,I,J] = Boundary_Conditions(positions,i,j,K,M,dt)
%   Boundary_Conditions Initialises a 3 by 3 matrix cut out of the
%   positions matrix at time dt, before applying the boundary conditions of
%   the simulation to make sure the movement of Tasmanian Devils adhere to
%   the simulation conditions.

    A = zeros(3,3);
    if i <= 1 && j <= 1
        I = 1:2;
        J = 1:2;
        A(2:3,2:3) = positions(I, J, dt);  
    elseif i <= 1 && j >= M
        I = 1:2;
        J = M-1:M;
        A(2:3,1:2) = positions(I, J, dt);
    elseif i >= K && j <= 1
        I = K-1:K;
        J = 1:2;        
        A(1:2,2:3) = positions(I, J, dt);
    elseif i <= 1
        I = 1:2;
        J = j-1:j+1;
        A(2:3,:) = positions(I, J, dt); 
    elseif j <= 1 
        I = i-1:i+1;
        J = 1:2;
        A(:,2:3) = positions(I, J, dt);
   elseif i >= K && j >= M
        I = K-1:K;
        J = M-1:M;
       A(1:2,1:2) = positions(I, J, dt);
    elseif i >= K
        I = K-1:K;
        J = j-1:j+1;  
        A(1:2,:) = positions(I, J, dt); 
    elseif j >= M 
        I = i-1:i+1;
        J = M-1:M;
        A(:,1:2) = positions(I, J, dt);
    else 
        I = i-1:i+1;
        J = j-1:j+1;
        A(:,:) = positions(I, J, dt);
    end
end