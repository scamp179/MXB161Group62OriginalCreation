function [positions,pop_change_percentage] = Tasmianian_Devil_sim(N, init_pos)
%       The Random_walk_cell_sim function uses a random walk given the 
%       inputs to derive and output the position of M nutrient cells for
%       N iterations. This function also outputs a tracking array which
%       shows the amount of nutrient cells inside the cell Membrane. 
    
%   Initialise variables
    Simulation_vid = VideoWriter('Simulation_vid');
    open(Simulation_vid);

    K = length(init_pos(:,1));
    M = length(init_pos(1,:));
    positions = zeros(K, M, N);
    positions(:, :,1) = init_pos;
    init_pop = sum(sum(init_pos));

%   Visualise Simulation    
    fig = figure('visible','on');
    handle = imshow(positions(:,:,1), 'InitialMagnification', 'Fit');
    colormap(parula)

    for dt = 2:N
%       Determine if simulation is in mating season. 
        weeks = mod(dt,52);
        if weeks >= 8 && weeks <= 20
            mating_season = 1;
        else
            mating_season = 0;
        end
        infect_status = positions(:, :, dt-1);
        for i = 1:K
            for j = 1:M
%               Check Boundary Conditions for input into cellular autonoma                
                if positions(i,j,dt-1) ~= 0
                    A = Boundary_Conditions(positions,i,j,K,M,dt-1);
%                   Use cellular autonoma to determine weightings for
%                   random walk.
                    [weighting_x, weighting_y, movement_type_weighting,infect_weight,Cub] = Cellular_autonomy_direction(A,mating_season);
%                   Perform Random  Walk for single Tasmanian Devil.                    
                    [dx,dy] = Random_walk_two_dim(positions(:, :, dt-1), j, i, weighting_x, weighting_y, movement_type_weighting);

%                   Check Boundary Conditions for next step 
                    if i == 1 && dy == -1
                        dy = 0;
                    elseif i == K && dy == 1
                        dy = 0;
                    end
                    if j == 1 && dx == -1
                        dx = 0;
                    elseif j == M && dx == 1
                        dx = 0;
                    end
%                   Cannot move to cell if already occupied
                    if positions(i + dy, j + dx, dt) == 1
                        dx = 0;
                        dy = 0;
                    end
%                   Update position and infection status of Tasmainian 
%                   Devil for not time iteration.
                    positions(i + dy, j + dx, dt) = update_infection_status(infect_status(i,j),infect_weight);
                    
                    if Cub > 0  
                        I = i + dy;
                        J = j + dx;
                        [B,y,x] = Boundary_Conditions(positions,I,J,K,M,dt);
                        R = length(y);
                        C = length(x);
                        for r = 1 : R
                            for c = 1 : C
                                if B(r,c) == 0
                                    row = y(r);
                                    col = x(c);
                                    positions(row, col , dt) = update_infection_status(Cub,infect_weight);

                                end
                            end
                        end
                    end 

                else
                    positions(i,j,dt) = positions(i,j,dt-1);
                end
            end
        end
        drawnow;
        Simulation_frames = getframe(fig);
        writeVideo(Simulation_vid, Simulation_frames);
        handle.CData = positions(:,:,dt);
     end
    
    final_pop = sum(sum(positions(:, :,end)>0));
    pop_change_percentage = (final_pop/init_pop)*100;
    drawnow;
    Simulation_frames = getframe(fig);
    writeVideo(Simulation_vid, Simulation_frames);
    close(Simulation_vid);
end


% -------------------------------------------------------------------------

function [updated_infect_status] = update_infection_status(infect_status,infect_weight)
%   The update_infection_status takes the current infection status of a
%   Tasmanian Devil and the infection weighting and determines what happens
%   to that Tasmanian Devil's infection status at the next time iteration.

    probability_infection = abs(randn(1));
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

function [A,I,J] = Boundary_Conditions(positions,i,j,K,M,dt)
%   Boundary_Conditions Initialises a 3 by 3 matrix cut out of the
%   positions matrix at time dt, before applying the boundary conditions of
%   the simulation to make sure the movement of Tasmanian Devils adhere to
%   the simulation conditions.

    A = zeros(3,3);
    if i == 1 && j == 1
        I = i:i+1;
        J = j:j+1;
        A(2:3,2:3) = positions(I, J, dt);  
    elseif i == 1 && j == M
        I = i:i+1;
        J = j-1:j;
        A(2:3,1:2) = positions(I, J, dt);
    elseif i == K && j == 1
        I = i-1:i;
        J = j:j+1;        
        A(1:2,2:3) = positions(I, J, dt);
    elseif i == 1
        I = i:i+1;
        J = j-1:j+1;
        A(2:3,:) = positions(I, J, dt); 
    elseif j == 1 
        I = i-1:i+1;
        J = j:j+1;
        A(:,2:3) = positions(I, J, dt);
   elseif i == K && j == M
        I = i-1:i;
        J = j-1:j;
       A(1:2,1:2) = positions(I, J, dt);
    elseif i == K
        I = i-1:i;
        J = j-1:j+1;  
        A(1:2,:) = positions(I, J, dt); 
    elseif j == M 
        I = i-1:i+1;
        J = j-1:j;
        A(:,1:2) = positions(I, J, dt);
    else 
        I = i-1:i+1;
        J = j-1:j+1;
        A(:,:) = positions(I, J, dt);
    end
end