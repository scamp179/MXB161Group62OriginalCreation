function [weighting_x, weighting_y, movement_type_weighting,Infection_status] = Cellular_autonomy_direction(A,mating_season)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% weighting_x is positive if going right (cols increasing)
% weighting_y is positive if going up (rows decreasing)
% movement_type_weighting is positive if tasmanian devil is moving in a
% straight line (up,down,left,right)

    B = A > 1;
    if A >= 1
        A = 1;
    else
        A = 0;
    end
    Y = 3;
    X = 3;

    north = [Y 1:Y-1];     % indices of north neighbour
    east  = [2:X 1];       % indices of east neighbour
    south = [2:Y 1];       % indices of south neighbour
    west  = [X 1:X-1];     % indices of west neighbour


    % Count how many live neighbours each cell has in its Moore neighbourhood
    live_neighbours = A(north, :) + A(south, :) + A(:, east) + A(:, west) ...
                    + A(north, east) + A(north, west) + A(south, east) + A(south, west);

    infected_neighbours = B(north, :) + B(south, :) + B(:, east) + B(:, west) ...
                    + B(north, east) + B(north, west) + B(south, east) + B(south, west);

    if mating_season == true
        var = 1;
    else
        var = -1;
    end
                    

    %Make infection Rules:
    Infection_rule_1; % surrounded by 1 infected 
    Infection_rule_2; % surrounded by 2 infected 
    Infection_rule_3; % surrounded by 3 infected 
    Infection_rule_4; % surrounded by 4 infected 
    Infection_rule_5; % surrounded by 5 infected 
    Infection_rule_6; % surrounded by 6 infected 
    Infection_rule_7; % surrounded by 7 infected 
    Infection_rule_8; % surrounded by 8 infected 

  
    % These two rules determine the new state of every element
    if A == Movement_rule_left
        weighting_x = -0.5; weighting_y = 0; movement_type_weighting = 0.5;

    elseif A == Movement_rule_right
        weighting_x = 0.5; weighting_y = 0; movement_type_weighting = 0.5;

    elseif A == Movement_rule_up
        weighting_x = 0; weighting_y = -0.5; movement_type_weighting = 0.5;

    elseif A == Movement_rule_down
        weighting_x = 0; weighting_y = 0.5; movement_type_weighting = 0.5;

    elseif A == Movement_rule_up_left
        weighting_x = -0.5; weighting_y = -0.5; movement_type_weighting = -0.5;

    elseif A == Movement_rule_up_right
        weighting_x = 0.5; weighting_y = -0.5; movement_type_weighting = -0.5;

    elseif A == Movement_rule_down_left
        weighting_x = -0.5; weighting_y = 0.5; movement_type_weighting = -0.5;

    elseif A == Movement_rule_down_right
        weighting_x = 0.5; weighting_y = 0.5; movement_type_weighting = -0.5;

    end
    
    if A(2,2) >= 2
        Infection_status = 0;
    elseif A == Infection_rule_8
        Infection_status = 8;
    elseif A == Infection_rule_7
        Infection_status = 7;
    elseif A == Infection_rule_6
        Infection_status = 6;
    elseif A == Infection_rule_5
        Infection_status = 5;
    elseif A == Infection_rule_4
        Infection_status = 4;
    elseif A == Infection_rule_3
        Infection_status = 3;
    elseif A == Infection_rule_2
        Infection_status = 2;
    elseif A == Infection_rule_1
        Infection_status = 1;
    else 
        Infection_status = 0;
    end
    
    
end

%-------------------------------------------------------------------------

function [weighting_x, weighting_y, movement_type_weighting,Infection_status] = apply_movement_rules(A,mating_season)
    
%Make Movement Rules:

    Movement_rule_up = [0,1,0;0,1,0;0,0,0] ;        % a cell lives if it has 3 live neighbours
    Movement_rule_down = [0,0,0;0,1,0;0,1,0];    % a cell lives if it's alive already, and has 2 live neighbours
    Movement_rule_left = [0,0,0;1,1,0;0,0,0];
    Movement_rule_right = [0,0,0;0,1,1;0,0,0];

    Movement_rule_up_left = [1,0,0;0,1,0;0,0,0]; 
    Movement_rule_up_right = [0,0,1;0,1,0;0,0,0];% a cell lives if it has 3 live neighbours
    Movement_rule_down_left = [0,0,0;0,1,0;1,0,0];    % a cell lives if it's alive already, and has 2 live neighbours
    Movement_rule_down_right = [0,0,0;0,1,0;0,0,1]; 


end