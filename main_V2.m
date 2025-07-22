clc; clear all;

% Gravitational parameter of the Sun
mu_sun = 1.32712440018e11; 

% Input Dates
initial_date = datetime(2032,12,31);
launch_date   = datetime(2037,4,1);

% Convert to Julian Dates
jd_initial = juliandate(initial_date);
jd_launch = juliandate(launch_date);

initial_id = 1;
launch_id = initial_id + (jd_launch - jd_initial); 
launch_window_days = 365 * 6;  % 6-year search window for launch
return_window_days = 365 * 6;  % 6-year search window for return


% Define time span for the propagation
% Initial Date for the Propagation
t_max = 10000 * 86400; % Propagate for one year (in seconds)
dt = 86400; % Time step of one day (in seconds)
T_earth_d = 365;
T_venus_d = 225;
T_mars_d = 687;
T_311p_d = 1180;


% Data Allocation
r_dataEarth = zeros(numel(0:dt:t_max),3);
r_dataVenus = zeros(numel(0:dt:t_max),3);
r_dataMars = zeros(numel(0:dt:t_max),3);
r_data311p = zeros(numel(0:dt:t_max),3);
v_dataEarth = zeros(numel(0:dt:t_max),3);
v_dataVenus = zeros(numel(0:dt:t_max),3);
v_dataMars = zeros(numel(0:dt:t_max),3);
v_data311p = zeros(numel(0:dt:t_max),3);



[r_earth, v_earth] = deal([-2.401326544689433E+07,  1.451318959619919E+08, -1.070974048417062E+04], ...
                           [-2.988748235147376E+01, -4.974336062964746E+00, -3.772304106908209E-05]);

[r_venus, v_venus] = deal([6.809146217201895E+07,  8.389728980863130E+07, -2.774481338997282E+06], ...
                           [-2.729977128261288E+01,  2.191785497588214E+01,  1.876489163063463E+00]);

[r_mars, v_mars] = deal([-2.446380464187839E+08, -2.465813948372017E+07,  5.478487412892363E+06], ...
                           [3.329974849896336E+00, -2.203658800742296E+01, -5.434652819462649E-01]);

[r_311p, v_311p] = deal([2.818788636618620E+08, -1.606433927476021E+08,  2.199377304407741E+07], ...
                           [7.875018129294897E+00,  1.865766925782450E+01,  9.320867968128654E-01]);


% Convert to Keplerian elements
[a_earth, e_earth, i_earth, RAAN_earth, omega_earth, M_earth] = stateToKepler(r_earth, v_earth, mu_sun);
[a_venus, e_venus, i_venus, RAAN_venus, omega_venus, M_venus] = stateToKepler(r_venus, v_venus, mu_sun);
[a_mars, e_mars, i_mars, RAAN_mars, omega_mars, M_mars] = stateToKepler(r_mars, v_mars, mu_sun);
[a_311p, e_311p, i_311p, RAAN_311p, omega_311p, M_311p] = stateToKepler(r_311p, v_311p, mu_sun);


% Loop through time and calculate the state vector
i=0;
for t = 0:dt:t_max
    i=i+1;
    [r_earth_new, v_earth_new] = ephemerisFunction(a_earth, e_earth, i_earth, RAAN_earth, omega_earth, M_earth, mu_sun, t);
    r_dataEarth(i,:)=r_earth_new;
    v_dataEarth(i,:)=v_earth_new;
    [r_venus_new, v_venus_new] = ephemerisFunction(a_venus, e_venus, i_venus, RAAN_venus, omega_venus, M_venus, mu_sun, t);
    r_dataVenus(i,:)=r_venus_new;
    v_dataVenus(i,:)=v_venus_new;
    [r_mars_new, v_mars_new] = ephemerisFunction(a_mars, e_mars, i_mars, RAAN_mars, omega_mars, M_mars, mu_sun, t);
    r_dataMars(i,:)=r_mars_new;
    v_dataMars(i,:)=v_mars_new;
    [r_311p_new, v_311p_new] = ephemerisFunction(a_311p, e_311p, i_311p, RAAN_311p, omega_311p, M_311p, mu_sun, t);
    r_data311p(i,:)=r_311p_new;
    v_data311p(i,:)=v_311p_new;
end


% Initialize arrays to store delta-v values for each launch-arrival pair
delta_v_matrix = NaN(launch_window_days, return_window_days);

% Loop through each possible combination of launch and arrival days
for launch_offset = 1:launch_window_days
    for return_offset = 1:return_window_days
        % Calculate the arrival day based on the launch day
        arrival_id = launch_id + launch_offset + return_offset;
        
        % Get the launch position and velocity
        r_launch = r_data311p(launch_id + launch_offset, :);
        v_launch = v_data311p(launch_id + launch_offset, :);
        
        % Get the arrival position and velocity
        r_arrival = r_dataEarth(arrival_id, :);
        v_arrival = v_dataEarth(arrival_id, :);

        % Calculate time of flight in seconds
        tof_return = (arrival_id - (launch_id + launch_offset)) * 86400; % days to seconds

        % Lambert transfer to find the initial velocity for return
        [v1_return, ~] = lambertBattin_V2(r_launch, r_arrival, tof_return, mu_sun, false);

        % Calculate Δv for the return departure
        dv_departure_return = norm(v1_return - v_launch);

        % Store the Δv in the matrix
        delta_v_matrix(launch_offset, return_offset) = dv_departure_return;
    end
end

% Plot the Pork Chop plot
figure; 
imagesc(delta_v_matrix);  % Display the Δv matrix as an image
colorbar;  % Show the color bar representing the Δv
xlabel('Arrival Day Offset from Launch');
ylabel('Launch Day Offset');
title('Pork Chop Plot for Δv (km/s)');
colormap('jet');  % Optional: Choose the colormap you prefer
axis equal;

% Optional: plot the best Δv found
[min_dv, idx] = min(delta_v_matrix(:));
[best_launch, best_return] = ind2sub(size(delta_v_matrix), idx);

fprintf("Best return found with Δv = %.3f km/s after %d days\n", min_dv, best_launch + best_return);



% ==FUNCTIONS==

function [a, e, i, RAAN, omega, M] = stateToKepler(r, v, mu)
    % Compute specific angular momentum
    h = cross(r, v);
    h_mag = norm(h);

    % Semi-major axis
    r_mag = norm(r);
    v_mag = norm(v);
    energy = 0.5 * v_mag^2 - mu / r_mag;
    a = -mu / (2 * energy);

    % Eccentricity
    e_vec = cross(v, h) / mu - r / r_mag;
    e = norm(e_vec);

    % Inclination
    i = acos(h(3) / h_mag);

    % RAAN (Right Ascension of Ascending Node)
    K = [0 0 1]; 
    n = cross(K, h);
    n_mag = norm(n);
    RAAN = acos(n(1) / n_mag);
    if n(2) < 0
        RAAN = 2*pi - RAAN;
    end
    
    % Argument of Periapsis (omega)
    omega = acos(dot(n, e_vec) / (n_mag * e));
    if e_vec(3) < 0
        omega = 2*pi - omega;
    end

    % Mean anomaly (M)
    n = sqrt(mu / a^3);
    
    cosE = dot(e_vec, r) / (e * r_mag);
    E = acos(max(min(cosE,1),-1));
    if dot(r, v) < 0
        E = 2*pi - E;
    end
    M = E - e * sin(E);
end

function [r_new, v_new] = ephemerisFunction(a, e, i, RAAN, omega, M0, mu, t)
    % Compute the new mean anomaly (M)
    n = sqrt(mu / a^3); % Mean motion
    M = M0 + n * t; % New mean anomaly
    
    % Solve Kepler's equation to get the eccentric anomaly (E)
    E = M; % Initial guess
    tol = 1e-8;
    for iter = 1:100
        delta_M = M - (E - e * sin(E)); % Difference between M and E
        if abs(delta_M) < tol
            break;
        end
        E = E + delta_M / (1 - e * cos(E));
    end
    
    % Calculate the true anomaly
    theta = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));
    
    % Calculate the radius
    r_mag = a * (1 - e^2) / (1 + e * cos(theta)); % Radial distance
    
    % Position and velocity in orbital plane
    r_orbit = [r_mag * cos(theta), r_mag * sin(theta), 0];
    v_orbit = sqrt(mu * a) / r_mag * [-sin(E), sqrt(1 - e^2) * cos(E), 0];
    
    % Rotate the orbit to the true orientation
    R = rotationMatrix(i, RAAN, omega); % Rotation matrix from orbital plane to 3D
    r_new = (R * r_orbit')';
    v_new = (R * v_orbit')';
end

function R = rotationMatrix(i, RAAN, omega)
    % Rotation matrix to rotate from orbital plane to 3D space
    % Using the classical orbital mechanics rotation matrix approach
    Rz_Omega = [cos(RAAN), sin(RAAN), 0; -sin(RAAN), cos(RAAN), 0; 0, 0, 1];
    Rx_i = [1, 0, 0; 0, cos(i), sin(i); 0, -sin(i), cos(i)];
    Rz_omega = [cos(omega), sin(omega), 0; -sin(omega), cos(omega), 0; 0, 0, 1];
    
    R = (Rz_omega * Rx_i * Rz_Omega)'; % Combined rotation matrix
end

function [times, positions] = propagateKepler(r0, v0, dt, mu, N)
    times = linspace(0, dt, N);
    positions = zeros(N, 3);
    for i = 1:N
        [r, ~] = keplerPropagate(r0, v0, times(i), mu);
        positions(i,:) = r';
    end
end

function [r_new, v_new] = keplerPropagate(r0, v0, dt, mu)
    r0 = r0(:); v0 = v0(:);
    r0_mag = norm(r0); v0_mag = norm(v0);
    vr0 = dot(r0, v0) / r0_mag;
    alpha = 2/r0_mag - v0_mag^2 / mu;
    chi = sqrt(mu) * abs(alpha) * dt;
    tol = 1e-8; max_iter = 100; ratio = 1; iter = 0;

    while abs(ratio) > tol && iter < max_iter
        iter = iter + 1;
        [S, C] = stumpff(alpha * chi^2);
        F = r0_mag * vr0 / sqrt(mu) * chi^2 * C + ...
            (1 - alpha * r0_mag) * chi^3 * S + r0_mag * chi - sqrt(mu)*dt;
        dF = r0_mag * vr0 / sqrt(mu) * chi * (1 - alpha * chi^2 * S) + ...
             (1 - alpha * r0_mag) * chi^2 * C + r0_mag;
        ratio = F / dF;
        chi = chi - ratio;
    end

    [S, C] = stumpff(alpha * chi^2);
    f = 1 - chi^2 / r0_mag * C;
    g = dt - 1/sqrt(mu) * chi^3 * S;
    r_new = f * r0 + g * v0;
    r_mag = norm(r_new);
    fdot = sqrt(mu) / (r_mag * r0_mag) * (alpha * chi^3 * S - chi);
    gdot = 1 - chi^2 / r_mag * C;
    v_new = fdot * r0 + gdot * v0;
end

function [S, C] = stumpff(z)
    if z > 0
        sqrtz = sqrt(z);
        S = (sqrtz - sin(sqrtz)) / sqrtz^3;
        C = (1 - cos(sqrtz)) / z;
    elseif z < 0
        sqrtz = sqrt(-z);
        S = (sinh(sqrtz) - sqrtz) / sqrtz^3;
        C = (cosh(sqrtz) - 1) / (-z);
    else
        S = 1/6; C = 1/2;
    end
end

function [v1, v2] = lambertBattin_V2(r1, r2, dt, mu, long_way)
    r1mag = norm(r1); r2mag = norm(r2);
    cos_theta = dot(r1, r2)/(r1mag * r2mag);
    theta = acos(max(min(cos_theta,1),-1));
    if long_way, theta = 2*pi - theta; end
    A = sin(theta) * sqrt(r1mag * r2mag / (1 - cos(theta)));
    z = 0; z_up = 4*pi^2; z_low = -4*pi^2;
    tol = 1e-8;

    for i = 1:1000
        [S, C] = stumpff(z);
        if C == 0, z = z + 1e-3; continue; end
        y = r1mag + r2mag + A * (z * S - 1)/sqrt(C);
        if y < 0, z = (z + z_up)/2; continue; end
        x = sqrt(y / C);
        dt_new = (x^3 * S + A * sqrt(y)) / sqrt(mu);
        if abs(dt_new - dt) < tol, break; end
        if dt_new <= dt, z_low = z; else, z_up = z; end
        z = (z_up + z_low) / 2;
    end

    f = 1 - y/r1mag; g = A * sqrt(y / mu); gdot = 1 - y/r2mag;
    v1 = (r2 - f*r1)/g; v2 = (gdot*r2 - r1)/g;
end
