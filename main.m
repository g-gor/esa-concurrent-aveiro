clear; clc;

%% Constants
mu_sun = 1.32712440018e11; % [km^3/s^2] Gravitational parameter of the Sun

%% Input Dates
launch_date = datetime(2034,7,6);     % Earth departure
flyby_date  = datetime(2035,1,1);    % Venus flyby
arrival_date = datetime(2037,3,17);   % Comet arrival

% Convert to Julian Dates
jd_launch = juliandate(launch_date);
jd_flyby = juliandate(flyby_date);
jd_arrival = juliandate(arrival_date);

tof1 = (jd_flyby - jd_launch) * 86400;    % [s]
tof2 = (jd_arrival - jd_flyby) * 86400;   % [s]

%% State Vectors (input from SPICE/Horizons)
% Earth at departure 
r_earth = [3.565437082223013E+07, -1.478533028879504E+08, 1.107340669310093E+04];
v_earth = [28.46935752110131, 6.859753702888095, 2.695344230252417E-04];

% Venus at flyby 
r_venus = [-8.563747067717355E+07,  6.470781423921068E+07,  5.831463654246382E+06];
v_venus = [-2.124682800426775E+01, -2.811814052660865E+01,  8.387179139458070E-01];

% Comet at arrival
r_comet = [-3.390260066184935E+07, 2.925066886870129E+08, 1.087589031898230E+06];
v_comet = [-22.15387280661516, -1.279809995212222, -1.919831648586088];

%% Solve Lambert Transfers

[v1_EM, v2_EM, dv1_EM, ~, vinf_in] = solveLambertOptimal(r_earth, v_earth, r_venus, v_venus, tof1, mu_sun);
[v1_MC, v2_MC, ~, dv2_MC, vinf_out] = solveLambertOptimal(r_venus, v_venus, r_comet, v_comet, tof2, mu_sun);

% Flyby turning angle
theta_flyby = acos(dot(vinf_in, vinf_out) / (norm(vinf_in)*norm(vinf_out)));
fprintf('Turning angle at Venus: %.2f deg\n', rad2deg(theta_flyby));

% Total Δv (departure + arrival)
dv_total = dv1_EM + dv2_MC;
fprintf('Departure Δv: %.3f km/s\n', dv1_EM);
fprintf('Arrival Δv: %.3f km/s\n', dv2_MC);
fprintf('Total mission Δv: %.3f km/s\n', dv_total);

%% Propagate and Plot
N = 500;

[~, r_leg1] = propagateKepler(r_earth, v1_EM, tof1, mu_sun, N);
[~, r_leg2] = propagateKepler(r_venus, v1_MC, tof2, mu_sun, N);

T_earth = 365.25 * 86400;  % [s]
T_venus  = 225 * 86400;  % [s]

[~, r_earth_full] = propagateKepler(r_earth, v_earth, T_earth, mu_sun, N);
[~, r_venus_full]  = propagateKepler(r_venus, v_venus, T_venus, mu_sun, N);

figure; hold on;
plot3(0,0,0,'yo','MarkerSize',8); % Sun
plot3(r_earth(1), r_earth(2), r_earth(3), 'bo', 'MarkerSize', 8);
plot3(r_venus(1), r_venus(2), r_venus(3), 'go', 'MarkerSize', 8);
plot3(r_comet(1), r_comet(2), r_comet(3), 'ro', 'MarkerSize', 8);
plot3(r_leg1(:,1), r_leg1(:,2), r_leg1(:,3), 'c-', 'LineWidth', 1.5);
plot3(r_leg2(:,1), r_leg2(:,2), r_leg2(:,3), 'm-', 'LineWidth', 1.5);
plot3(r_earth_full(:,1), r_earth_full(:,2), r_earth_full(:,3), 'b--', 'LineWidth', 1.2);
plot3(r_venus_full(:,1), r_venus_full(:,2), r_venus_full(:,3), 'g--', 'LineWidth', 1.2);


legend('Sun', 'Earth', 'Venus (Flyby)', 'Comet', 'Leg 1: Earth→Venus', 'Leg 2: Venus→Comet');
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Lambert Transfer with Venus Gravity Assist');
axis equal; grid on; view(3);

%% === Functions ===

function [v1, v2, dv1, dv2, vinf] = solveLambertOptimal(r1, v1_ref, r2, v2_ref, tof, mu)
    % Try both short and long way and pick best
    [v1a, v2a] = lambertBattin_V2(r1, r2, tof, mu, false);
    [v1b, v2b] = lambertBattin_V2(r1, r2, tof, mu, true);

    dv1a = norm(v1a - v1_ref); dv2a = norm(v2a - v2_ref);
    dv1b = norm(v1b - v1_ref); dv2b = norm(v2b - v2_ref);

    total_a = dv1a + dv2a;
    total_b = dv1b + dv2b;

    if total_a < total_b
        v1 = v1a; v2 = v2a; dv1 = dv1a; dv2 = dv2a;
    else
        v1 = v1b; v2 = v2b; dv1 = dv1b; dv2 = dv2b;
    end

    vinf = v2 - v2_ref;
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
