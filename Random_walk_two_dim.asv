function [delta_x,delta_y] = Random_walk_two_dim(position, x, y, weighting_x, weighting_y, movement_type_weighting)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
    movement_type = 2*rand(1)-1; % [0 to -1] is diagonal, [1 to 0] is straight (up,down,left,right)

    Y = length(position(:,1));
    X = length(position(1,:));
    if movement_type < 0 - movement_type_weighting
        [delta_x,delta_y] = rand_walk_diag(x,y,weighting_x,weighting_y);
    elseif movement_type > 0 - movement_type_weighting %  movement_type > 0 - movement_type_weighting moves in straight lines
        [delta_x,delta_y] = rand_walk_straight(x,y,weighting_x,weighting_y,X,Y);
    else % No movement
        delta_x = 0;
        delta_y = 0;
    end
    
end

% -------------------------------------------------------------------------

function [delta_x,delta_y] = rand_walk_straight(x, y, weighting_x, weighting_y, X, Y)

    x_step = 2*rand(1)-1;
    y_step = 2*rand(1)-1;

    if abs(weighting_x*x_step) > abs(weighting_y*y_step)
        if x_step > 0 - weighting_x % step right
            % if weighted in + x direction then more likely to step in + x direction
            % if weighted in - x direction then more likely to step in - x direction
            if x == X % Boundary Condition Cut-off
                delta_x = 0;
    %            delta_y = 0; 
            else
                delta_x = 1;
    %            delta_y = 0;
            end
    
        elseif x_step < 0 - weighting_x % step left
            % if weighted in + x direction then more likely to step in + x direction
            % if weighted in - x direction then more likely to step in - x direction
            if x == 1 % Boundary Condition Cut-off
                delta_x = 0;
                delta_y = 0; 
            else
                delta_x = -1;
                delta_y = 0;
            end
        end
    elseif abs(weighting_x*x_step) < abs(weighting_y*y_step)
        if y_step > 0 - weighting_y % step up
            % if weighted in + x direction then more likely to step in + x direction
            % if weighted in - x direction then more likely to step in - x direction
            if y == 1 % Boundary Condition Cut-off
                delta_x = 0;
                delta_y = 0; 
            else
                delta_x = 0;
                delta_y = 1;
            end
    
        elseif y_step < 0 - weighting_y % step down
            % if weighted in + x direction then more likely to step in + x direction
            % if weighted in - x direction then more likely to step in - x direction
            if y == Y % Boundary Condition Cut-off
                delta_x = 0;
                delta_y = 0; 
            else
                delta_x = 0;
                delta_y = -1;
            end
    
        end
    else % no movement
        delta_x = 0;
        delta_y = 0;
    end
end

% -------------------------------------------------------------------------

function [delta_x,delta_y] = rand_walk_diag(x, y, weighting_x, weighting_y, X, Y)
    x_step = 2*rand(1)-1;
    y_step = 2*rand(1)-1;

    if x_step > 0 - weighting_x  && y_step > 0 - weighting_y % step right % up
        dx = 1; dy = -1;
        [delta_x,delta_y] = Boundary_Conditions_diag(x,y,1,Y,dx,dy);
        
    elseif x_step < 0 - weighting_x  && y_step > 0 - weighting_y% step left % up
        dx = -1; dy = -1;
        [delta_x,delta_y] = Boundary_Conditions_diag(x,y,1,1,dx,dy);

    elseif x_step < 0 - weighting_x  && y_step < 0 - weighting_y% step left % down
        dx = -1; dy = 1;
        [delta_x,delta_y] = Boundary_Conditions_diag(x,y,X,1,dx,dy);

    elseif  x_step > 0 - weighting_x  && y_step < 0 - weighting_y % step right % down
        dx = 1; dy = 1;      
        [delta_x,delta_y] = Boundary_Conditions_diag(x,y,X,Y,dx,dy);

    else % no movement
        delta_x = 0;
        delta_y = 0;
    end
end

% -------------------------------------------------------------------------

function [delta_x,delta_y] = Boundary_Conditions_diag(x,y,bound_x, bound_y, dx, dy)
 
     if y == bound_y && x ~= bound_x % Boundary Condition Cut-off
         delta_x = dx;
         delta_y = 0; 
    elseif y ~= bound_y && x == bound_x % Boundary Condition Cut-off
        delta_x = 0;
        delta_y = dy;
    elseif y == bound_y && x == bound_x % Boundary Condition Cut-off
        delta_x = 0;
        delta_y = 0;
     else
        delta_x = dx;
        delta_y = dy;     
     end
end
