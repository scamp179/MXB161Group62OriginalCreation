steps = 255;
cells = 521;

A = zeros(steps+1, cells);
A(1,260) = 1;

for n = 1:steps
    
    P = A(n, :);
    L = [0 P(1:cells-1)];
    R = [P(2:cells) 0];

    C000 = L == 0 & P == 0 & R == 0;
    C001 = L == 0 & P == 0 & R == 1;
    C010 = L == 0 & P == 1 & R == 0;
    C011 = L == 0 & P == 1 & R == 1;
    C100 = L == 1 & P == 0 & R == 0;
    C101 = L == 1 & P == 0 & R == 1;
    C110 = L == 1 & P == 1 & R == 0;
    C111 = L == 1 & P == 1 & R == 1;
    
    live_mask = C011 | C110 | C100 | C001;
    A(n+1, live_mask) = 1;
end

    figure
imshow(~A, 'InitialMagnification', 'fit') % live cells are black