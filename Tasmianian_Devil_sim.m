function [positions,pop_change_percentage] = Tasmianian_Devil_sim(N, init_pos)
%       The Random_walk_cell_sim function uses a random walk given the 
%       inputs to derive and output the position of M nutrient cells for
%       N iterations. This function also outputs a tracking array which
%       shows the amount of nutrient cells inside the cell Membrane. 

    K = length(init_pos(:,1));
    M = length(init_pos(1,:));
    positions = zeros(K, M, N);
    positions(:, :,1) = init_pos;
    init_pop = Sum(init_pos);
    A = zeros(3,3);
    for dt = 2:N
        weeks = mod(dt,52);
        if weeks >= 8 && weeks <= 20
            mating_season = 1;
        else
            mating_season = 0;
        end
        for i = 1:M
            for j = 1:M
                if i == 1 && j ==1
                    A(2:3,:) = positions(i:i+1, j:j+1, dt-1);          
                elseif i == 1
                    A(2:3,:) = positions(i:i+1, j-1:j+1, dt-1); 
                elseif j == 1 
                    A(:,2:3) = positions(i-1:i+1, j:j+1, dt-1);
                else 
                    A(:,:) = positions(i-1:i+1, j-1:j+1, dt-1);
                end
                [weighting_x, weighting_y, movement_type_weighting,infect_weight] = Cellular_autonomy_direction(A,mating_season);
                infect_status = positions(:, :, dt-1);
                [dx,dy] = Random_walk_two_dim(positions(:, :, dt-1), i, j, weighting_x, weighting_y, movement_type_weighting);

                if positions(i + dx,j + dy, dt) > 0
                    dx = 0;
                    dy = 0;
                end
                positions(i + dx,j + dy, dt) = update_infection_status(infect_status,infect_weight);
            end
        end

    end
    final_pop = Sum(positions(:, :,end)>0);
    pop_change_percentage = (final_pop - init_pop)/100;
end


% -------------------------------------------------------------------------

function [updated_infect_status] = update_infection_status(infect_status,infect_weight)
    probability_infection = randn(0,1);
    if  infect_status == 1 && infect_weight ~= 0
        if infect_status == 1 && probability_infection < (infect_weight*0.05 + 0.5)
            updated_infect_status = 2;  
        else
            updated_infect_status = infect_status;
        end
    elseif infect_stat >= 2 % Infected, Greater infect_status closer to death until infect_stat == 7
        if infect_status == 27
            updated_infect_status = 0; % Dead tasmanian devil
        else
            updated_infect_status = infect_status + 1;  
        end
    elseif infect_status == 1 && infect_weight == 0   
        updated_infect_status = infect_status;
    end
end