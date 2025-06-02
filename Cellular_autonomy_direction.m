function [weighting_x, weighting_y, movement_type_weighting,Infection_weighting] = Cellular_autonomy_direction(A,mating_season)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% weighting_x is positive if going right (cols increasing)
% weighting_y is positive if going up (rows decreasing)
% movement_type_weighting is positive if tasmanian devil is moving in a
% straight line (up,down,left,right)
    B = zeros(3,3);
    C = zeros(3,3);
    F = zeros(3,3);
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
        
            if A(i,j) < 0
                F(i,j) = 1;
            else
                F(i,j) = 0;
            end
        end
    end
    Y = 3;
    X = 3;

    north = [Y 1:Y-1];     % indices of north neighbour
    east  = [2:X 1];       % indices of east neighbour
    south = [2:Y 1];       % indices of south neighbour
    west  = [X 1:X-1];     % indices of west neighbour


    % Count how many infected neighbours each cell has in its Moore neighbourhood
    infected_neighbours = B(north, :) + B(south, :) + B(:, east) + B(:, west) ...
                    + B(north, east) + B(north, west) + B(south, east) + B(south, west);


    %Make infection Rules:
    Infection_rule_0 = sum(sum(infected_neighbours)) == 0; % surrounded by 0 infected 
    Infection_rule_1 = sum(sum(infected_neighbours)) == 1; % surrounded by 1 infected 
    Infection_rule_2 = sum(sum(infected_neighbours)) == 2; % surrounded by 2 infected 
    Infection_rule_3 = sum(sum(infected_neighbours)) == 3; % surrounded by 3 infected 
    Infection_rule_4 = sum(sum(infected_neighbours)) == 4; % surrounded by 4 infected 
    Infection_rule_5 = sum(sum(infected_neighbours)) == 5; % surrounded by 5 infected 
    Infection_rule_6 = sum(sum(infected_neighbours)) == 6; % surrounded by 6 infected 
    Infection_rule_7 = sum(sum(infected_neighbours)) == 7; % surrounded by 7 infected 
    Infection_rule_8 = sum(sum(infected_neighbours)) == 8; % surrounded by 8 infected 

    if mating_season == true
        var = 1;
    else
        var = -1;
    end

    north_num_neigh = sum(C(1,:));
    south_num_neigh = sum(C(3,:)); 
    east_num_neigh = sum(C(:,3));
    west_num_neigh = sum(C(:,1));
    

    weighting_y = var*((north_num_neigh - south_num_neigh)*0.3);
    weighting_x = var*((east_num_neigh - west_num_neigh)*0.3);
    
    if weighting_y == 0 || weighting_x == 0
        movement_type_weighting = 0.5;
    elseif  weighting_y == 0 && weighting_x == 0
        movement_type_weighting = 0;
    else 
        movement_type_weighting = -0.5;
    end
    
    if A(2,2) >= 2
        Infection_weighting = 0;
    else
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
            Infection_status = 0;
        end
        Infection_weighting = Infection_status* 0.12;
    end
     
end