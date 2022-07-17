function xhat = kalman(z, A, C, R, Q)
% Initializtion =====================================================

% Number of sensors
m = size(C, 1);

% Number of state values
n = size(C, 2);

% Number of observations
numobs = size(z, 2);

% Use linear least squares to estimate initial state from initial
% observation
xhat = zeros(n, numobs);
xhat(:,1) = C \ z(:,1);
% Initialize P, I
P = ones(size(A));
I = eye(size(A));

% Kalman Filter =====================================================

for k = 2:numobs
    
    % Predict
    xhat(:,k) = A * xhat(:,k-1);
    P         = A * P * A' + Q;
    
    % Update
    G         = P  * C' / (C * P * C' + R);
    P         = (I - G * C) * P;
    xhat(:,k) = xhat(:,k) + G * (z(:,k) - C * xhat(:,k));
    
end