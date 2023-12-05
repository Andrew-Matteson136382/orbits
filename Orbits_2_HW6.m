clear
close all
clc

%% Introduction %%
%----------------%
%Programmer:     A. Clifford Matteson
%Date:           10/31/2023
%Class:          AE 5614: Spaceflight Mechanics II


%% Constants %%
%-------------%

const.r_earth = 6378.14;
const.mu_earth = 3.986*10^5;
pi2deg = 180/pi;
deg2pi = pi/180;
const.AU2km = 1.496*10^8;
const.G = 6.6738*10^-20;
const.mu_sun = 1.327*10^11;
const.r_saturn = 58232;
const.mu_saturn = 3.7931*10^7;
const.g = 9.8067;
const.a_saturn = 9.54327*const.AU2km;
const.r_sun = 696300;
const.geo_orb = 35785 + const.r_earth;
const.iss_orb = 409 + const.r_earth;
const.mu_moon = 4903;
const.r_moon = 1737.5;


%% Equations %%
%-------------%

% Energey Equation Derivations
Spe_Eng.E1 = @(v, mu, r) (v^2/2)-mu/r;
Spe_Eng.E2 = @(mu, a) -mu/(2*a);
Spe_Eng.v = @(mu, r, a) sqrt(2*((-mu/(2*a))+(mu./r))); %
Spe_Eng.r = @(mu, v, a) -mu./(-(v.^2)/2-mu./(2*a)); %
Spe_Eng.a = @(mu, v, r) -mu/(2*((v.^2/2)-mu./r)); %

% Orbit Equation Derivations
Orbit.r = @(p, ecc, nu) p/(1+ecc*cos(nu));
Orbit.p = @(r, ecc, nu) r*(1+ecc*cos(nu));
Orbit.ecc = @(r, p, nu) (p/r-1)/cos(nu);
Orbit.nu = @(r, p, ecc) acos((p./r-1)/ecc);

% Parameter Derivations
Para.p1 = @(a, ecc) a*(1-ecc^2);
Para.p2 = @(h, mu) (h^2)/mu;
Para.ecc1 = @(p, a) sqrt(1-p/a);
Para.ecc2 = @(h, a, mu) sqrt(-h.^2./(a*mu)+1); %
Para.a1 = @(p, ecc) p/(1-ecc^2);
Para.a2 = @(h, mu, ecc) h^2/(mu*(1-ecc^2));
Para.h1 = @(p, mu) sqrt(mu*p);
Para.h2 = @(a, ecc, mu) sqrt(a*mu*(1-ecc^2));

% Theta Velocity Derivations
T_Hat.v1 = @(v, gamma) v*cos(gamma);
T_Hat.v2 = @(mu, h, ecc, nu) (mu/h)*(1+ecc*cos(nu));
T_Hat.v3 = @(h,r) h/r;
T_Hat.gamma1 = @(the_v, v) acos(the_v/v);
T_Hat.gamma2 = @(h, r, v) acos(h./(r.*v)); %
T_Hat.nu = @(the_v, mu, h, ecc) acos(((the_v*h/mu)-1)/ecc);
T_Hat.h = @(r, v, gamma) r.*v.*cos(gamma);

% Radial Velocity Derivations
R_Hat.v1 = @(v, gamma) v*sin(gamma);
R_Hat.gamma1 = @(r_v, v) asin(r_v/v);
R_Hat.nu = @(r_v, mu, h, ecc) asin((r_v.*h)./(mu*ecc));

% Flight Relationships Derivations
Fli_Rel.ecc = @(r, v, mu, gamma) sqrt(((((r*v^2)/mu-1)^2)*cos(gamma)^2) ...
    +sin(gamma)^2);
Fli_Rel.nu = @(r, v, mu, gamma) atan((((r*v^2)/mu)*sin(gamma)*cos(gamma)) ...
    /(((r*v^2)/mu)*(cos(gamma)^2)-1));

% Eccentricity Derivation
Eccen.ecc = @(E, h, mu) sqrt(1+(2*E*h^2)/(mu^2));

% Elliptical Orbit Derivations
Ell_Orb.P = @(a, mu) 2*pi*sqrt((a^3)/mu);
Ell_Orb.r = @(a, ecc, E) a*(1-ecc*cos(E));
Ell_Orb.E = @(a, r, e) acos((-r./a+1)/e);
Ell_Orb.ra = @(a, e) a*(1+e);
Ell_Orb.rp = @(a, e) a*(1-e);
Ell_Orb.e1 = @(rp, a) 1-rp/a;
Ell_Orb.e2 = @(ra, a) -1+ra/a;
Ell_Orb.a1 = @(rp, e) rp/(1-e);
Ell_Orb.a2 = @(ra, e) ra/(1+e);
Ell_Orb.a3 = @(P, mu) ((P/(2*pi))^2*mu)^(1/3);

% Hyperbolic Orbits Derivations
Hype_Orb.rp = @(a, ecc) a*(1-ecc);
Hype_Orb.v = @(mu, a) sqrt(-mu/a);
Hype_Orb.nu = @(p, r, ecc) acos(((p/r)-1)/ecc);
Hype_Orb.a = @(mu, v) -mu/(v^2);
Hype_Orb.del = @(e) 2*asin(1/e);
Hype_Orb.e = @(del) 1/sin(de/2);
Hype_Orb.e1 = @(rp, a) 1-rp/a;

% Lambert's Theorem Derivations
LambertTime.c = @(r1, r2, phi) sqrt(r1^2+r2^2-2*r1*r2*cos(phi));
LambertTime.s = @(r1, r2, c) (r1+r2+c)/2;

% Elliptical Transfers Derivations
Ell_Trans.alpha = @(s, a) 2*asin(sqrt(s/(2*abs(a))));
Ell_Trans.beta = @(s, a, c) 2*asin(sqrt((s-c)/(2*abs(a))));
Ell_Trans.Tdel_1A = @(mu, a, alpha, beta) ((alpha-sin(alpha))-(beta-sin(beta)))/sqrt(mu/abs(a)^3);
Ell_Trans.Tdel_1B = @(mu, a, alpha, beta) ((alpha-sin(alpha))-(beta-sin(beta))+2*pi)/sqrt(mu/abs(a)^3);
Ell_Trans.Tdel_2A = @(mu, a, alpha, beta) ((alpha-sin(alpha))+(beta-sin(beta)))/sqrt(mu/abs(a)^3);
Ell_Trans.Tdel_2B = @(mu, a, alpha, beta) ((alpha-sin(alpha))+(beta-sin(beta))+2*pi)/sqrt(mu/abs(a)^3);
% 1 = (0 <= x < 180) , 2 =  (180 <= x < 360)
% A = Antifocus in region, B = Antifocus not in region

% Hyperbolic Transfers Derivations
Hype_Trans.a_prime = @(s, a) 2*asinh(sqrt(s/(2*abs(a))));
Hype_Trans.Tra_b_prime = @(s, a, c) 2*asinh(sqrt((s-c)/(2*abs(a))));
Hype_Trans.Tra_Tdel_2H = @(mu, a, a_prime, b_prime) ((sinh(a_prime)-a_prime)+ ...
    (sinh(b_prime)-b_prime))/sqrt(mu/abs(a)^3);
Hype_Trans.Tra_Tdel_1H = @(mu, a, a_prime, b_prime) ((sinh(a_prime)-a_prime)- ...
    (sinh(b_prime)-b_prime))/sqrt(mu/abs(a)^3);

% True & Ecc Anomaly Derivations
True_Ecc.E = @(nu, ecc) 2*atan(tan(nu./2)/sqrt((1+ecc)/(1-ecc)));
True_Ecc.nu = @(E, ecc) 2*atan(sqrt((1+ecc)/(1-ecc))*tan(E/2));

% Kepler % Rad
Kepler.n = @(mu, a) sqrt(mu/(a^3));
Kepler.Tdel = @(E1, E2, ecc, n) (E2-E1-ecc.*(sin(E2)-sin(E1)))/n;
Kepler.ti = @(M, n, tp) M/n+tp;


% Fly By Conics Derivations
Fly_By_Con.Phi = @(v_arr, gamma, v_pla) atan((v_arr*sin(gamma))/(v_arr*cos(gamma)-v_pla));

% Law of Cosine Derivatio
Law_Cos = @(v1, v2, theta) sqrt(v1^2+v2^2-2*v1*v2*cos(theta));

% Rocket Equation Derivations
Rocket.delV = @(isp, mo, mf) isp*0.00980665*log(mo/mf);
Rocket.m0 = @(mf, delV, isp) mf*exp(-delV/(isp*0.00980665));
Rocket.mf = @(m0, delV, isp) m0./exp(delV/(isp*0.00980665));
Rocket.isp = @(m0, mf, delV) delV/(0.00980665*log(m0/mf));
Rocket.delM = @(mf, delV, isp) mf.*(exp(delV./(isp*0.00980665))-1);

% F and G Functions
fgfunc.f = @(r,p,delNu) 1-r/p*(1-cos(delNu));
fgfunc.g = @(r0,r1,mu,p,delNu) r1*r0/((mu*p)^(1/2))*sin(delNu);

%% Loop Setup %%

% TOF (Time of Flight)
departure = 86400*1.5;
increment = 1800;
arrival = 86400*20;
TOF = (departure:increment:arrival);
lenTOF = length(TOF);

% TD (Time Delay)
start = 0;
inc = 30;
stop = 86400/4;
TD = (start:inc:stop);
lenTD = length(TD);

% Position Array 4D
PosArr = single(zeros(lenTD,lenTOF,4,3));
VelArr = single(zeros(lenTD,lenTOF,4,3));

% Orbit Array 3D
OrbArr = single(zeros(lenTD,lenTOF,2));

% VDel Array 3D
VDelArr = single(zeros(lenTD,lenTOF,3));

% Black Box velocity 3D
BB_VArr = single(zeros(lenTD,lenTOF,2,3));

% True Anomally Array 4D
NuArr = single(zeros(lenTD,lenTOF,4));

%% ISS Setup %%
%-------------%

% ISS inputs %
ISS_Epoch = Date2Julian(2023,291.53041186);
ISS_Rev = 15.50357024420925;
ISS_P = 1/ISS_Rev*86400;
ISS_a = Ell_Orb.a3(ISS_P, const.mu_earth);
ISS_i = 51.6404*deg2pi;
ISS_Omega = 77.5756*deg2pi;
ISS_ecc = 0.0004521;
ISS_omega = 114.7526*deg2pi;
ISS_M = 245.3933*deg2pi;

% Finds ISS's Mean Motion (n)
ISS_n = Kepler.n(const.mu_earth, ISS_a);

% Finds the initial change in time wrt tp (s)
ISS_fE = @(E) ISS_M-(E+ISS_ecc*sin(E));
ISS_fdE = @(E) ISS_ecc*cos(E)-1;
ISS_E = newtons_method(ISS_fE, ISS_fdE, ISS_M);
ISS_tDel = Kepler.Tdel(0, ISS_E, ISS_ecc, ISS_n);


%% Luna Setup %%
%-------------%

% Luna inputs %
Luna_Epoch = 2460235.5;
Luna_a = 3.843233396654867*10^5;
Luna_ecc = 3.913088756831110*10^-2;
Luna_i = 2.832505621905977*10^1*deg2pi;
Luna_Omega = 4.628106711973097*deg2pi;
Luna_omega = 3.380320838890249*10^2*deg2pi;
Luna_M = 2.643966828517082*10^2*deg2pi;
Luna_tp = 2460242.743662629277;

% Finds Luna's Mean Motion (n)
Luna_n = Kepler.n(const.mu_earth, Luna_a);

% Finds Luna's initial time (s)
Luna_ti = Kepler.ti(Luna_M, Luna_n, Luna_tp);

% Luna is 12 hours behide so propagate forward ~13 hours
% Finds Luna's ISS synchronous Mean Anomaly (M) (Rad)
Luna_M = Luna_n*(Luna_ti+12*3600+43*60+47.5847-Luna_tp);

% Finds Luna's ISS synchronous initial time (s)
Luna_ti = Kepler.ti(Luna_M, Luna_n, Luna_tp);

%% Black Box Setup %%

JJ = 1;
n = 10;
tol = 1*10^-6;
kmax = 10;

%% Loop %%

for i = 1:lenTD
    for j = 1:lenTOF

    %% Luna TD %%

    % Finds Luna's New Mean Anomaly wrt a time delay (M)
    Luna_DelM_TD = Luna_n*(Luna_ti+TD(i)-Luna_tp);

    % Finds Eccentric Anomaly for Luna wrt a time delay (M) 
    Luna_fE_TD = @(E) Luna_DelM_TD-(E+Luna_ecc*sin(E));
    Luna_fdE_TD = @(E) Luna_ecc*cos(E)-1;
    Luna_E_TD = newtons_method(Luna_fE_TD, Luna_fdE_TD, Luna_DelM_TD);

    % Finds r magnitude for Luna wrt a time delay (km)
    Luna_r_TD = Ell_Orb.r(Luna_a, Luna_ecc, Luna_E_TD);

    % Finds True Anomaly for Luna wrt a time delay (rad)
    Luna_nu_tan_TD = [True_Ecc.nu(Luna_E_TD, Luna_ecc), pi+True_Ecc.nu(Luna_E_TD, Luna_ecc)];
    Luna_nu_cos_TD = [Orbit.nu(Luna_r_TD, Para.p1(Luna_a, Luna_ecc), Luna_ecc), -Orbit.nu(Luna_r_TD, Para.p1(Luna_a, Luna_ecc), Luna_ecc)];
    Luna_nu_TD = quadcheck([Luna_nu_tan_TD Luna_nu_cos_TD]);
    NuArr(i,j,1) = Luna_nu_TD;

    % Finds the R and V vectors for Luna wrt a time delay 
    [PosArr(i,j,2,:), VelArr(i,j,2,:)] = Class2Cart(Luna_a, Luna_ecc, Luna_i, Luna_omega, Luna_Omega, Luna_nu_TD, const.mu_earth,"unit","Rad");

    %% ISS TD %%

    % Finds ISS's New Mean Anomaly wrt a time delay (M)
    ISS_DelM_TD = ISS_n*(ISS_tDel+TD(i));

    % Finds Eccentric Anomaly for ISS wrt a time delay (M) 
    ISS_fE_TD = @(E) ISS_DelM_TD-(E+ISS_ecc*sin(E));
    ISS_fdE_TD = @(E) ISS_ecc*cos(E)-1;
    ISS_E_TD = newtons_method(ISS_fE_TD, ISS_fdE_TD, ISS_DelM_TD);

    % Finds r magnitude for ISS wrt a time delay (M) 
    ISS_r_TD = Ell_Orb.r(ISS_a, ISS_ecc, ISS_E_TD);

    % Finds True Anomaly for ISS wrt a time delay (M) 
    ISS_nu_tan_TD = [True_Ecc.nu(ISS_E_TD, ISS_ecc), pi+True_Ecc.nu(ISS_E_TD, ISS_ecc)];
    ISS_nu_cos_TD = [Orbit.nu(ISS_r_TD, Para.p1(ISS_a, ISS_ecc), ISS_ecc), -Orbit.nu(ISS_r_TD, Para.p1(ISS_a, ISS_ecc), ISS_ecc)];
    ISS_nu_TD = quadcheck([ISS_nu_tan_TD ISS_nu_cos_TD]);
    NuArr(i,j,2) = ISS_nu_TD;

    % Finds the R and V vectors for the ISS wrt a time delay 
    % Care about ISS_V_TD
    [PosArr(i,j,1,:), ISS_V_TD] = Class2Cart(ISS_a, ISS_ecc, ISS_i, ISS_omega, ISS_Omega, ISS_nu_TD, const.mu_earth,"unit","Rad");
    VelArr(i,j,1,:) = ISS_V_TD;

    %% Luna TD + TOF %%

    % Finds Luna's New Mean Anomaly wrt a time delay & time of flight (M) 
    Luna_DelM_TOF = Luna_n*(Luna_ti+TD(i)+TOF(j)-Luna_tp);

    % Finds Eccentric Anomaly for Luna wrt a time delay & time of flight (M) 
    Luna_fE_TOF = @(E) Luna_DelM_TOF-(E+Luna_ecc*sin(E));
    Luna_fdE_TOF = @(E) Luna_ecc*cos(E)-1;
    Luna_E_TOF = newtons_method(Luna_fE_TOF, Luna_fdE_TOF, Luna_DelM_TOF);

    % Finds r magnitude for Luna wrt a time delay & time of flight (km)
    Luna_r_TOF = Ell_Orb.r(Luna_a, Luna_ecc, Luna_E_TOF);

    % Finds True Anomaly for Luna wrt a time delay & time of flight (rad)
    Luna_nu_tan_TOF = [True_Ecc.nu(Luna_E_TOF, Luna_ecc), pi+True_Ecc.nu(Luna_E_TOF, Luna_ecc)];
    Luna_nu_cos_TOF = [Orbit.nu(Luna_r_TOF, Para.p1(Luna_a, Luna_ecc), Luna_ecc), -Orbit.nu(Luna_r_TOF, Para.p1(Luna_a, Luna_ecc), Luna_ecc)];
    Luna_nu_TOF = quadcheck([Luna_nu_tan_TOF Luna_nu_cos_TOF]);
    NuArr(i,j,3) = Luna_nu_TOF;

    % Finds the R and V vectors for Luna wrt a time delay & time of flight 
    [PosArr(i,j,4,:), Luna_V_TOF] = Class2Cart(Luna_a, Luna_ecc, Luna_i, Luna_omega, Luna_Omega, Luna_nu_TOF, const.mu_earth,"unit","Rad");
    VelArr(i,j,4,:) = Luna_V_TOF;

    %% ISS TD + TOF %%

    % Finds ISS's New Mean Anomaly wrt a time delay & time of flight (M) 
    ISS_DelM_TOF = ISS_n*(ISS_tDel+TD(i)+TOF(j));

    % Finds Eccentric Anomaly for ISS wrt a time delay & time of flight (M) 
    ISS_fE_TOF = @(E) ISS_DelM_TOF-(E+ISS_ecc*sin(E));
    ISS_fdE_TOF = @(E) ISS_ecc*cos(E)-1;
    ISS_E_TOF = newtons_method(ISS_fE_TOF, ISS_fdE_TOF, ISS_DelM_TOF);

    % Finds r magnitude for ISS wrt a time delay & time of flight (M) 
    ISS_r_TOF = Ell_Orb.r(ISS_a, ISS_ecc, ISS_E_TOF);

    % Finds True Anomaly for ISS wrt a time delay & time of flight (M) 
    ISS_nu_tan_TOF = [True_Ecc.nu(ISS_E_TOF, ISS_ecc), pi+True_Ecc.nu(ISS_E_TOF, ISS_ecc)];
    ISS_nu_cos_TOF = [Orbit.nu(ISS_r_TOF, Para.p1(ISS_a, ISS_ecc), ISS_ecc), 2*pi-Orbit.nu(ISS_r_TOF, Para.p1(ISS_a, ISS_ecc), ISS_ecc)];
    ISS_nu_TOF = quadcheck([ISS_nu_tan_TOF ISS_nu_cos_TOF]);
    NuArr(i,j,4) = ISS_nu_TOF;

    % Finds the R and V vectors for the ISS wrt a time delay & time of flight
    [PosArr(i,j,3,:), ISS_V_TOF] = Class2Cart(ISS_a, ISS_ecc, ISS_i, ISS_omega, ISS_Omega, ISS_nu_TOF, const.mu_earth,"unit","Rad");
    VelArr(i,j,3,:) = ISS_V_TOF;


    % Black box, finds transfers a, e, and velocities at points
    % Had to modify to compute norms of N-D arrays
    [OrbArr(i,j,1), OrbArr(i,j,2), BB_V1, BB_V2, BB_conv] = Lambert(PosArr(i,j,1,:), PosArr(i,j,4,:), TOF(j), const.mu_earth, JJ, n, tol, kmax);
    BB_VArr(i,j,1,:) = BB_V1;
    BB_VArr(i,j,2,:) = BB_V2;

    VDelArr(i,j,1) = norm(BB_V1 - ISS_V_TD,"fro");
    VDelArr(i,j,2) = norm(BB_V2 - Luna_V_TOF,"fro");

    end
end

figure(1)
subplot(1,3,1)
contourf(TD/3600,TOF/86400,VDelArr(:,:,1)')
ylabel('Time of Flight (Days)')
xlabel('Time of Departure (Hours)')
title("TLI at " + ISS_Epoch + " (km/s)")

subplot(1,3,2)
contourf(TD/3600,TOF/86400,VDelArr(:,:,2)')
ylabel('Time of Flight (Days)')
xlabel('Time of Departure (Hours)')
title("Vinf (-) to moon at " + ISS_Epoch + " (km/s)")


VDelArr(:,:,3) = VDelArr(:,:,2) + VDelArr(:,:,1);
subplot(1,3,3)
contourf(TD/3600,TOF/86400,VDelArr(:,:,3)')
ylabel('Time of Flight (Days)')
xlabel('Time of Departure (Hours)')
title("Total VDel at " + ISS_Epoch + " (km/s)")

% Finds the value and instant with the lowest total delta V
vtli = min(VDelArr(:,:,1),[],"all");
[tlilTD,tliTOF] = find(VDelArr(:,:,1)==vtli);
vinf = min(VDelArr(:,:,2),[],"all");
[vinfTD,vinflTOF] = find(VDelArr(:,:,2)==vinf);
tdel = min(VDelArr(:,:,3),[],"all");
[tdelTD,tdelTOF] = find(VDelArr(:,:,3)==tdel);

tmp = squeeze(BB_VArr(vinfTD,vinflTOF,2,:))' - squeeze(VelArr(vinfTD,vinflTOF,4,:))';
save("output.mat","tmp")

figure(10)
[X,Y,Z] = sphere;
surf(X*const.r_earth,Y*const.r_earth,Z*const.r_earth)
axis('equal')
hold on
surf(const.r_moon*X+PosArr(tdelTD,tdelTOF,2,1), Y*const.r_moon + PosArr(tdelTD,tdelTOF,2,2), Z*const.r_moon + PosArr(tdelTD,tdelTOF,2,3))
surf(const.r_moon*X+PosArr(tdelTD,tdelTOF,4,1), Y*const.r_moon + PosArr(tdelTD,tdelTOF,4,2), Z*const.r_moon + PosArr(tdelTD,tdelTOF,4,3))

plot3(PosArr(tdelTD,tdelTOF,1,1), PosArr(tdelTD,tdelTOF,1,2), PosArr(tdelTD,tdelTOF,1,3),'o','MarkerSize',10,"LineWidth",2)
plot3(PosArr(tdelTD,tdelTOF,3,1), PosArr(tdelTD,tdelTOF,3,2), PosArr(tdelTD,tdelTOF,3,3),'o','MarkerSize',10,'LineWidth',2)

%% Part 2: BB analysis

% Inital BB velocity vec and magnitude
BB_V1 = squeeze(BB_VArr(tdelTD,tdelTOF,1,:))';
Bb_V1_mag = vecnorm(BB_V1);

% BB orbit semimajor axis and parameter
BB_a = OrbArr(tdelTD,tdelTOF,1);
BB_p = OrbArr(tdelTD,tdelTOF,2);
BB_ecc = Para.ecc1(BB_p, BB_a);

% BB Inital and final Position vectors and mag
BB_Pos1 = squeeze(PosArr(tdelTD,tdelTOF,1,:))';
BB_Pos1Mag = vecnorm(BB_Pos1);
BB_Pos2 = squeeze(PosArr(tdelTD,tdelTOF,4,:))';
BB_Pos2Mag = vecnorm(BB_Pos2);

% Inital true anomally 
BB_Pos1_Nu1 = [Orbit.nu(BB_Pos1Mag,BB_p, BB_ecc), -Orbit.nu(BB_Pos1Mag,BB_p, BB_ecc)];
BB_Pos1_Nu2 = [asin((1-((BB_p/BB_Pos1Mag-1)/BB_ecc)^2)^(1/2)), pi-asin((1-((BB_p/BB_Pos1Mag-1)/BB_ecc)^2)^(1/2))];
BB_Pos1_Nu = quadcheck([BB_Pos1_Nu1 BB_Pos1_Nu2]);

% Final true anomally 
BB_Pos2_Nu1 = [Orbit.nu(BB_Pos2Mag,BB_p, BB_ecc), 2*pi-Orbit.nu(BB_Pos2Mag,BB_p, BB_ecc)];
BB_Pos2_Nu2 = [asin((1-((BB_p/BB_Pos2Mag-1)/BB_ecc)^2)^(1/2)), pi-asin((1-((BB_p/BB_Pos2Mag-1)/BB_ecc)^2)^(1/2))];
BB_Pos2_Nu = quadcheck([BB_Pos2_Nu1 BB_Pos2_Nu2]);

BB_TranInc = (0:pi/180:2*pi-BB_Pos2_Nu-BB_Pos1_Nu); % Increments for BB during transfer
BBTran = zeros(length(BB_TranInc), 3);

for i = 1:length(BB_TranInc)
    BB_r = Orbit.r(BB_p,BB_ecc,BB_TranInc(i)+BB_Pos1_Nu);
    BB_f = fgfunc.f(BB_r, BB_p, BB_TranInc(i));
    BB_g = fgfunc.g(BB_Pos1Mag, BB_r, const.mu_earth, BB_p, BB_TranInc(i));
    BBTran(i,:) = BB_Pos1*BB_f+BB_V1*BB_g;
end

plot3(BBTran(:,1),BBTran(:,2),BBTran(:,3),'g')

%% Moon Ploting

moon_nu0 = NuArr(tdelTD,tdelTOF,1); % Inital angle
moon_nu1 = NuArr(tdelTD,tdelTOF,3); % Final angle
moon_p = Para.p1(Luna_a,Luna_ecc); % Parameter of moon
moonPos0 = PosArr(tdelTD,tdelTOF,2,:); % Initial Moon pos and vel
moonVel = VelArr(tdelTD,tdelTOF,2,:);
moon_r0 = vecnorm(moonPos0); % Inital radius of moon

moon_TranInc = (0:pi/180:moon_nu1-moon_nu0); % Increments for moon during transfer
moon_OrbInc = (0:pi/180:2*pi);

% Array to save output values to
moonTran = zeros(length(moon_TranInc), 3);
moonOrb = zeros(length(moon_OrbInc), 3);

for i = 1:length(moon_TranInc)
    moon_r = Orbit.r(moon_p,Luna_ecc,moon_nu0+moon_TranInc(i));
    moon_f = fgfunc.f(moon_r, moon_p, moon_TranInc(i));
    moon_g = fgfunc.g(moon_r0, moon_r, const.mu_earth, moon_p, moon_TranInc(i));
    moonTran(i,:) = moonPos0*moon_f+moonVel*moon_g;
end

for i = 1:length(moon_OrbInc)
    moon_r = Orbit.r(moon_p,Luna_ecc,moon_nu0+moon_OrbInc(i));
    moon_f = fgfunc.f(moon_r, moon_p, moon_OrbInc(i));
    moon_g = fgfunc.g(moon_r0, moon_r, const.mu_earth, moon_p, moon_OrbInc(i));
    moonOrb(i,:) = moonPos0*moon_f+moonVel*moon_g;
end

plot3(moonOrb(:,1),moonOrb(:,2),moonOrb(:,3),'b',LineWidth=1)
plot3(moonTran(:,1),moonTran(:,2),moonTran(:,3),'r',LineWidth=2)

%% ISS Plotting

ISS_nu0 = NuArr(tdelTD,tdelTOF,2); % Inital angle
ISS_nu1 = NuArr(tdelTD,tdelTOF,4); % Final angle
ISSPos0 = PosArr(tdelTD,tdelTOF,1,:); % Initial Moon pos and vel
ISS_p = Para.p1(ISS_a,ISS_ecc); % Parameter of moon
ISSVel = VelArr(tdelTD,tdelTOF,1,:);
ISS_r0 = vecnorm(ISSPos0); % Inital radius of moon

ISS_TranInc = (0:pi/180:ISS_nu1-ISS_nu0); % Increments for moon during transfer
ISS_OrbInc = (0:pi/180:2*pi);

% Array to save output values to
ISSTran = zeros(length(ISS_TranInc), 3);
ISSOrb = zeros(length(ISS_OrbInc), 3);

for i = 1:length(ISS_TranInc)
    ISS_r = Orbit.r(ISS_p,ISS_ecc,ISS_TranInc(i)+ISS_nu0);
    ISS_f = fgfunc.f(ISS_r, ISS_p, ISS_TranInc(i));
    ISS_g = fgfunc.g(ISS_r0, ISS_r, const.mu_earth, ISS_p, ISS_TranInc(i));
    ISSTran(i,:) = ISSPos0*ISS_f+ISSVel*ISS_g;
end

for i = 1:length(ISS_OrbInc)
    ISS_r = Orbit.r(ISS_p,ISS_ecc,ISS_OrbInc(i)+ISS_nu0);
    ISS_f = fgfunc.f(ISS_r, ISS_p, ISS_OrbInc(i));
    ISS_g = fgfunc.g(ISS_r0, ISS_r, const.mu_earth, ISS_p, ISS_OrbInc(i));
    ISSOrb(i,:) = ISSPos0*ISS_f+ISSVel*ISS_g;
end

plot3(ISSOrb(:,1),ISSOrb(:,2),ISSOrb(:,3),'b',LineWidth=1)
plot3(ISSTran(:,1),ISSTran(:,2),ISSTran(:,3),'r',LineWidth=4)
