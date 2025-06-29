function [weighting_x, weighting_y, movement_type_weighting,Infection_weighting,Cub_stat,Breeding_weighting,future_neighboursless_step] = Cellular_autonomy_direction(A,A2,mating_season)
% Cellular_autonomy_direction determines direction and infection based movement weightings for an individual cell based on its 3x3 neighbourhood. 
%   
% INPUT:
%     A - 3x3 matrix surrounding a cell of previous iteration
%     A2 - 3x3 matrix surrounding a cell of current iteration
%     mating_season = Boolean: True if in mating season, False if not
%
% OUTPUT:
%     weighting_x is positive if going right (cols increasing)
%     weighting_y is positive if going up (rows decreasing)
%     movement_type_weighting is positive if tasmanian devil is moving in a
%     straight line (up,down,left,right)
%     infection_weighting is infection pressure based on infected neighbours
%     Cub_stat is the indication if cell can breed
%     Breeding_weighting is the percentage of producing a new cell
%     future_neighboursless_step is an array of possible movements

    % Sets B, C and F to be empty 3x3 matrices
    B = zeros(3,3);
    C = zeros(3,3);
    B2 = zeros(3,3);
    % Populates B, C and F based on different conditions applied to A
    for i = 1 : 3
        for j = 1 : 3
            if A(i,j) > 1
                B (i,j) = 1;
            else
                B (i,j) = 0; 
            end
            if A(i,j) >= 1
                C(i,j) = 1;
            else
                C(i,j) = 0;
            end
        
            if A2(i,j) < 0
                B2(i,j) = 1;
            else
                B2(i,j) = 0;
            end
        end
    end
    Y = 3;
    X = 3;
    % Define neighbouring indices for cardinal directions
    north = [Y 1:Y-1];     % indices of north neighbour
    east  = [2:X 1];       % indices of east neighbour
    south = [2:Y 1];       % indices of south neighbour
    west  = [X 1:X-1];     % indices of west neighbour


    % Count how many infected neighbours each cell has in its neighbourhood
    infected_neighbours = B(north, :) + B(south, :) + B(:, east) + B(:, west) ...
                    + B(north, east) + B(north, west) + B(south, east) + B(south, west);

    num_infected_neighbours= sum(sum(infected_neighbours== 0 ));
    % Define infection Rules based on total number of infected neighbours:
    Infection_rule_1 = sum(sum(num_infected_neighbours)) == 1; % surrounded by 1 infected 
    Infection_rule_2 = sum(sum(num_infected_neighbours)) == 2; % surrounded by 2 infected 
    Infection_rule_3 = sum(sum(num_infected_neighbours)) == 3; % surrounded by 3 infected 
    Infection_rule_4 = sum(sum(num_infected_neighbours)) == 4; % surrounded by 4 infected 
    Infection_rule_5 = sum(sum(num_infected_neighbours)) == 5; % surrounded by 5 infected 
    Infection_rule_6 = sum(sum(num_infected_neighbours)) == 6; % surrounded by 6 infected 
    Infection_rule_7 = sum(sum(num_infected_neighbours)) == 7; % surrounded by 7 infected 
    Infection_rule_8 = sum(sum(num_infected_neighbours)) == 8; % surrounded by 8 infected 

   % Count how many regular neighbours each cell has in its neighbourhood
    breeding_neighbours = C(north, :) + C(south, :) + C(:, east) + C(:, west) ...
                    + C(north, east) + C(north, west) + C(south, east) + C(south, west);
    % Breeding rules
    num_breeding_neighbours = sum(sum(breeding_neighbours == 1));
    Breeding_rule = ((num_breeding_neighbours >= 1) && (mating_season == 1));

    % Future neighbours
    straight_dir = [north;south;east;west];
    upupdowndown = [north;north;south;south;];
    leftrightleftright = [west;east;west;east];

    future_neighboursless_step = zeros(8,2);
    for i = 1 : 8
        if i == 1 || i == 3 || i == 5 || i == 7
            if B2(straight_dir(i), :) == 0
                if isequal(straight_dir(i),north)
                    dx = 0; 
                    dy = -1;
                elseif isequal(straight_dir(i),east)
                    dx = 1;
                    dy = 0;
                elseif isequal(straight_dir(i),south)
                    dx = 0; 
                    dy = 1;
                else
                    dx = -1;
                    dy = 0;
                end
                future_neighboursless_step(i,:) = [dy,dx];
            else
                future_neighboursless_step(i,:) = [0,0];
            end
        elseif i == 2 || i == 4 || i == 6 || i == 8
            if B2(upupdowndown(i), leftrightleftright(i)) == 0
                if isequal(upupdowndown(i),north)
                    dy = -1; 
                else
                    dy = 1;
                end
                if isequal(leftrightleftright(i),west)
                    dx = -1; 
                else
                    dx = 1;
                end
                future_neighboursless_step(i,:) = [dy,dx];
            else
                future_neighboursless_step(i,:) = [0,0];
            end
        else
            message = ('fuckin error');
            disp(message);
            return
        end
    end

    % sets variable based on mating season
    % If it is mating season, var = 1 (cells come together)
    % If it is not mating season, var = -1 (cells avoid each other)
    if mating_season == 1
        var = 1;
    else
        var = -1;
    end

    % count number of neighbours in all cardinal directions
    north_num_neigh = sum(C(1,:));
    south_num_neigh = sum(C(3,:)); 
    east_num_neigh = sum(C(:,3));
    west_num_neigh = sum(C(:,1));
    
    % Set directional weighting based on variable (up is +ve y, right is +ve x)
    weighting_y = var*((north_num_neigh - south_num_neigh)*0.3);
    weighting_x = var*((east_num_neigh - west_num_neigh)*0.3);

    % Determines movement type weighting
    % Straight movement is preferred: one direction non zero
    % No movement: both zero
    % Diagonal movement: both non zero
    if weighting_y == 0 || weighting_x == 0
        movement_type_weighting = 0.5;
    elseif  weighting_y == 0 && weighting_x == 0
        movement_type_weighting = 0;
    else 
        movement_type_weighting = -0.5;
    end

    % Determine infection weighting
    % If the centre cell is already infected, weighting = 0
    if A(2,2) >= 2
        Infection_weighting = 0;
    else
        % identify how many infected neighbours are present using previous rules
        if Infection_rule_8
            Infection_status = 8;
        elseif Infection_rule_7
            Infection_status = 7;
        elseif Infection_rule_6
            Infection_status = 6;
        elseif Infection_rule_5
            Infection_status = 5;
        elseif Infection_rule_4
            Infection_status = 4;
        elseif Infection_rule_3
            Infection_status = 3;
        elseif Infection_rule_2
            Infection_status = 2;
        elseif Infection_rule_1
            Infection_status = 1;
        else 
            Infection_status = rand>0.95; % 5% chance to get infection status
        end
        Infection_weighting = Infection_status* 0.12;
    end
    
        if Breeding_rule
                Cub_stat = 1;
            % if num_breeding_neighbours
                %mod(num_breeding_neighbours,8) ~= 0 % 8 is max number of neighbours
                %Breeding_weighting = round(mod(num_breeding_neighbours,8),2)/2;
            if (num_breeding_neighbours == 2)
                Breeding_weighting = 0.05;
            %elseif num_breeding_neighbours == 4
                %Breeding_weighting = 0.32;
            %elseif num_breeding_neighbours == 8
                %Breeding_weighting = 0.1;
            else
                Breeding_weighting = 0;
            end
        else
            Breeding_weighting = 0 ;
            Cub_stat = 0;
        end

end
