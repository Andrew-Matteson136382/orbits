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

%% NRHO

NRHO_P = 6.562*60*60*24;
NRHO_rp = 3366;
NRHO_a = Ell_Orb.a3(NRHO_P, const.mu_earth)
NRHO_e = Ell_Orb.e1(NRHO_rp, NRHO_a)
