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
Orbit.nu2 = @(r, p, ecc) asin(sqrt(1-((-1+p/r)/(ecc))^2));

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
R_Hat.v2 = @(mu, ecc, nu, h) mu*ecc*sin(nu)/h;

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

%% NRHO

Luna_Epoch = 2460235.5;
Luna_a = 3.843233396654867*10^5;
Luna_ecc = 3.913088756831110*10^-2;
Luna_i = 2.832505621905977*10^1*deg2pi;
Luna_Omega = 4.628106711973097*deg2pi;
Luna_omega = 3.380320838890249*10^2*deg2pi;
Luna_M = 2.643966828517082*10^2*deg2pi;
Luna_tp = 2460242.743662629277;

load('output.mat')
ECI_vinf = tmp;

Omega = 0;
Theta = 0;

DCM_E_M = [cos(Omega)*cos(Theta)-sin(Omega)*cos(Luna_i)*sin(Theta), -sin(Theta)*cos(Omega)-sin(Omega)*cos(Luna_i)*cos(Theta), sin(Omega)*sin(Luna_i);
    sin(Omega)*cos(Theta)+cos(Omega)*cos(Luna_i)*sin(Theta), -sin(Omega)*sin(Theta)+cos(Omega)*cos(Luna_i)*cos(Theta), -cos(Omega)*sin(Luna_i);
    sin(Luna_i)*sin(Theta), sin(Luna_i)*cos(Theta), cos(Luna_i)];
MCI_vinf = ECI_vinf*DCM_E_M';
MCI_vinf_Mag = vecnorm(MCI_vinf);

MCI_hyp_a = Hype_Orb.a(const.mu_moon,MCI_vinf_Mag);

r_guess = 20000;
gamma_guess = -pi/4;
Omega_guess = pi/4;

MCI_hyp_v = Spe_Eng.v(const.mu_moon, r_guess, MCI_hyp_a);

MCI_hyp_e = Fli_Rel.ecc(r_guess, MCI_hyp_v, const.mu_moon, gamma_guess);
MCI_hyp_p = Para.p1(MCI_hyp_a, MCI_hyp_e);
MCI_hyp_h = Para.h1(MCI_hyp_p, const.mu_moon);

MCI_hyp_vr = R_Hat.v1(MCI_hyp_v, gamma_guess);
MCI_hyp_vt = T_Hat.v1(MCI_hyp_v, gamma_guess);
MCI_hyp_vh = 0;
MCI_hyp_vel = [MCI_hyp_vr; MCI_hyp_vt; MCI_hyp_vh];

NRHO_P = 6.562*60*60*24;
NRHO_rp = 3366;
NRHO_a = Ell_Orb.a3(NRHO_P, const.mu_moon);
NRHO_e = Ell_Orb.e1(NRHO_rp, NRHO_a);
NRHO_ra = Ell_Orb.ra(NRHO_a, NRHO_e);
NRHO_p = Para.p1(NRHO_a, NRHO_e);
NRHO_nu1 = [Orbit.nu(r_guess, NRHO_p, NRHO_e), -Orbit.nu(r_guess, NRHO_p, NRHO_e)];
NRHO_nu2 = [Orbit.nu2(r_guess, NRHO_p, NRHO_e), pi-Orbit.nu2(r_guess, NRHO_p, NRHO_e)];
NRHO_nu = quadcheck([NRHO_nu1 NRHO_nu2]);
NRHO_h = Para.h1(NRHO_p, const.mu_moon);
NRHO_vel = [R_Hat.v2(const.mu_moon, NRHO_e, NRHO_nu, NRHO_h); T_Hat.v2(const.mu_moon, NRHO_h, NRHO_e, NRHO_nu);0];
NRHO_theta = NRHO_nu+pi/2;
DCM_HRHO2cart = DCMArr(pi/2, Omega_guess, NRHO_theta);
NRHO_cartvel = DCM_HRHO2cart*NRHO_vel

MCI_hyp_nhat = [cos(Omega_guess), sin(Omega_guess), 0];


MCI_hyp_hhat = (cross(MCI_hyp_nhat,MCI_vinf))/(vecnorm(cross(MCI_hyp_nhat,MCI_vinf)));

MCI_hyp_i = acos(dot(MCI_hyp_hhat,[0,0,1]));

MCI_hyp_nu1 = [Fli_Rel.nu(r_guess, MCI_hyp_v, const.mu_moon, gamma_guess), pi+Fli_Rel.nu(r_guess, MCI_hyp_v, const.mu_moon, gamma_guess)];
MCI_hyp_nu2 = [acos(1/sqrt(1+((r_guess*MCI_hyp_v^2/const.mu_moon*sin(gamma_guess)*cos(gamma_guess))/(-1+r_guess*MCI_hyp_v^2/const.mu_moon*cos(gamma_guess)^2))^2)), -acos(1/sqrt(1+((r_guess*MCI_hyp_v^2/const.mu_moon*sin(gamma_guess)*cos(gamma_guess))/(-1+r_guess*MCI_hyp_v^2/const.mu_moon*cos(gamma_guess)^2))^2))];
MCI_hyp_nu = quadcheck([MCI_hyp_nu1 MCI_hyp_nu2]);

MCI_hyp_TolDel = Hype_Orb.del(MCI_hyp_e);

MCI_hyp_phi = asin((vecnorm(cross(MCI_hyp_nhat,MCI_vinf)))/(MCI_vinf_Mag));
if dot(MCI_hyp_nhat,MCI_vinf)<0
    MCI_hyp_phi = pi - MCI_hyp_phi;
end

MCI_hyp_beta = pi/2-(MCI_hyp_TolDel)/2;

MCI_hyp_alpha = pi - MCI_hyp_phi;

MCI_hyp_omega = pi-MCI_hyp_beta-MCI_hyp_alpha;

DCM_hyp2cart = DCMArr(MCI_hyp_i, Omega_guess, MCI_hyp_nu+MCI_hyp_omega);

MCI_hyp_cartvel = DCM_hyp2cart*MCI_hyp_vel

Vdel = NRHO_cartvel-MCI_hyp_cartvel








