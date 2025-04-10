clc; clear all;


% Gravitational parameter of the Sun
mu_sun = 1.32712440018e11; 

% Define time span for the propagation
t_max = 365.25 * 86400; % Propagate for one year (in seconds)
dt = 86400; % Time step of one day (in seconds)

r_dataEarth = zeros(numel(0:dt:t_max),3);

[r_earth, v_earth] = deal([-2.401326544689433E+07,  1.451318959619919E+08, -1.070974048417062E+04], ...
                           [-2.988748235147376E+01, -4.974336062964746E+00, -3.772304106908209E-05]);

% Convert to Keplerian elements
[a_earth, e_earth, i_earth, RAAN_earth, omega_earth, M_earth] = stateToKepler(r_earth, v_earth, mu_sun);



% Loop through time and calculate the state vector
i=0;
for t = 0:dt:t_max
    i=i+1;
    [r_earth_new, v_earth_new] = ephemerisFunction(a_earth, e_earth, i_earth, RAAN_earth, omega_earth, M_earth, mu_sun, t);
    r_dataEarth(i,:)=r_earth_new;
end

figure; hold on;
plot3(0,0,0,'yo','MarkerSize',8); % Sun
plot3(r_dataEarth(:,1), r_dataEarth(:,2), r_dataEarth(:,3), 'b--', 'LineWidth', 1.2);

legend('Sun', 'Earth Orbit');
xlabel('X [km]'); ylabel('Y [km]');
axis equal; grid on; view(3);

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
    RAAN = atan2(h(1), -h(2));
    
    % Argument of Periapsis (omega)
    omega = atan2(dot(cross(h, e_vec), r) / r_mag, dot(e_vec, r) / r_mag);

    % Mean anomaly (M)
    n = sqrt(mu / a^3);
    E = acos(dot(e_vec, r) / (e * r_mag)); % Eccentric anomaly
    M = E - e * sin(E); % Mean anomaly
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
    
    R = Rz_Omega * Rx_i * Rz_omega; % Combined rotation matrix
end