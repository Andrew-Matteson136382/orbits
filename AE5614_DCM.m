clear all
close all
clc

%% Introduction %%
%----------------%
%Programmer:     Clifford Matteson
%Date:           09/26/2023
%Class:          AE 5614: Aerospace Mechancics II

mu = 3.986e5; % Standard gravitational parameter of Earth (km^3/s^2)

[r1, v1] = classical_to_cartesian(1025415, .0111, 30, 64, 94, 4, mu);
[a, e, i, omega, Omega, nu] = cartesian_to_classical(r1, v1, mu, 'deg');


function [r, v] = classical_to_cartesian(a, e, i, omega, Omega, nu, mu)
    % Constants
    deg_to_rad = pi / 180;

    % Convert angles from degrees to radians
    i = i * deg_to_rad;
    omega = omega * deg_to_rad;
    Omega = Omega * deg_to_rad;
    nu = nu * deg_to_rad;
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
    gamma = gamma_1 + gamma_2;
    r = DCM * [r_mag; 0; 0];
    v = DCM * [v_mag*sin(gamma); v_mag*cos(gamma); 0];
end

function [a, e, i, omega, Omega, nu] = cartesian_to_classical(r, v, mu, a_unit)
    h = cross(r, v);
    r_mag = norm(r);
    v_mag = norm(v);

    a = -mu/(2*(v_mag^2/2-mu/r_mag));
    e_cap = 1/mu*((v_mag^2-mu/r_mag)*r-dot(r,v)*v);
    e = norm(e_cap);
    i = acos(h(3) / norm(h));
    a_cap = (cross([0;0;1],h))/(norm(cross([0;0;1],h)));
    if a_cap(2)<0
        Omega = 2*pi-acos(dot([1;0;0],a_cap));
    else
        Omega = acos(dot([1;0;0],a_cap));
    end
    if e_cap(3)<0
        omega = 2*pi-acos(dot(a_cap,e_cap)/e);
    else
        omega = acos(dot(a_cap,e_cap)/e);
    end
    if dot(r,v)<0
        nu = 2*pi-acos(dot(e_cap,r)/(e*r_mag));
    else
        nu = acos(dot(e_cap,r)/(e*r_mag));
    end
    if a_unit == 'deg'
        i = i*180/pi;
        Omega = Omega*180/pi;
        omega = omega*180/pi;
        nu = nu*180/pi;
    else
    end
end
