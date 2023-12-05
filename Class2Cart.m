function [r, v] = Class2Cart(a, e, i, omega, Omega, nu, mu, opt)
    arguments
        a;
        e;
        i;
        omega;
        Omega;
        nu;
        mu;
        opt.unit(1,1) = "deg";
    end

    % Constants
    deg_to_rad = pi / 180;

    % Convert angles from degrees to radians
    if opt.unit == "deg"
        i = i * deg_to_rad;
        omega = omega * deg_to_rad;
        Omega = Omega * deg_to_rad;
        nu = nu * deg_to_rad;
    end

    theta = omega+nu;

    % Compute the direction cosine matrix (DCM)
    DCM = [cos(Omega) * cos(theta) - sin(Omega) * sin(theta) * cos(i), -cos(Omega) * sin(theta) - sin(Omega) * cos(theta) * cos(i), sin(Omega) * sin(i);
           sin(Omega) * cos(theta) + cos(Omega) * sin(theta) * cos(i), -sin(Omega) * sin(theta) + cos(Omega) * cos(theta) * cos(i), -cos(Omega) * sin(i);
           sin(theta) * sin(i), cos(theta) * sin(i), cos(i)];

    % Compute the position and velocity vectors in the Cartesian coordinate system
    r_mag = a * (1 - e^2) / (1 + e * cos(nu));
    v_mag = (2*(-mu/(2*a)+mu/r_mag))^(1/2);
    h = (mu*a*(1-e^2))^(1/2);
    gamma_1 = round([acos(h/(r_mag*v_mag)), -acos(h/(r_mag*v_mag))],8);
    gamma_2 = round([asin((mu*e*sin(nu)/h)/v_mag), pi-asin((mu*e*sin(nu)/h)/v_mag)],8);
    gamma = quadcheck([gamma_1 gamma_2]);
    r = DCM * [r_mag; 0; 0];
    v = DCM * [v_mag*sin(gamma); v_mag*cos(gamma); 0];
end