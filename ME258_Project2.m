%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Randy Qafaiti
% ME 258
% Heat Transfer and Phase Change Thermophysics

% Project #2
% Due on 4/15/2022 at 11:59 PM

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear % Clears variables
close all % Closes opened figures
clc % Clears command line window

%% Figure Initialization

% Strings for saving figures
root = 'C:\Users\Randy Qafaiti\Documents\MATLAB\UC Berkeley\ME 258_Convection and Phase Change Thermophysics\Project2\';
root = pwd; % It stores current folder directory, and I guess it might work for both of us.
ext = '.png';

% Graph variables
fs = 16; % font size
lw = 2; % line width
f = 0; % figure counter

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User-defined parameters

% Pipe parameters
D_p = 15.0e-3; % pipe inside diameter [m]
L_evap = 25.0e-2; % evaproator section length [m]
L_adiab = 25.0e-2; % adiabatic section length [m]
L_cond = 25.0e-2; % condenser section length [m]

% Wick parameters
t_w = 0.457e-3; % wick thickness [m]
omega = 0.68; % wick porosity
r_p = 0.063e-3; % wick pore radius [m]
kappa = 1.94e-10; % wick permeability [m^2]
theta = 10; % contact angle [deg]
k_m = 386; % thermal conductivity of metal used in wick [W/(m*K)]

% Select working fluid (ONLY SET ONE FLUID EQUAL TO 1 AT A TIME)
water = 1;
R245fa = 0;

% Air temperature
T_a_in = 41.0; % [C]
T_a_out = 25.5; % [C]

% Ratio of fin area to section area
A_a_in_A_e = 5.0; % inlet fin area to evaporator area
A_a_out_A_c = 5.0; % outlet fin area to condenser area

% Air convection coefficients
h_a_in = 93.0; % [W/(m^2*C)]
h_a_out = 93.0; % [W/(m^2*C)]

%% Saturation properties

% Develop curve-fits to saturation properties k_l(T), rho_l(T), rho_v(T),
% mu_l(T), mu_v(T), T_sat(T), and P_sat(T) of water and R245fa over the
% temperature range 0C - 80C.

% Samples data within temperature range 0C - 80C
T_water_sample = [0.01, 20, 40, 60, 80] + 273.15; % temperature [C]
P_water_sample = [0.00061, 0.00234, 0.00738, 0.0195, 0.04741]*1e6; % Pressure [Pa]
k_l_water_sample = [561.0, 598.4, 630.6, 654.3, 670.0]*1e-3; % thermal conductivity of liquid [W/(m*K)]
rho_l_water_sample = [999.8, 998.2, 992.2, 983.2, 971.8]; % density of liquid [kg/m^3]
rho_v_water_sample = 1./[205.99, 57.757, 19.515, 7.6672, 3.4052]; % density of vapor [kg/m^3]
mu_l_water_sample = [1791.2, 1001.6, 653.0, 466.4, 354.3]*1e-6; % viscosity of liquid [Pa*s]
mu_v_water_sample = [9.22, 9.73, 10.31, 10.93, 11.59]*1e-6; % viscosity of vapor [Pa*s]
sigma_water_sample = [75.65, 72.74, 69.60, 66.24, 62.67]*1e-3; % surface tension of liquid [N/m]
h_lv_water = [0.00, 83.91, 167.53, 251.18, 335.01]*1e3; % latent heat of vaporization [J/kg]

T_R245fa_sample = [0, 20, 40, 60, 80] + 273.15; % temperature [C]
P_R245fa_sample = [0.05359, 0.12380, 0.25179, 0.46353, 0.78882]*1e6; % Pressure [Pa]
k_l_R245fa_sample = [88.8, 82.3, 76.1, 70.2, 64.5]*1e-3; % thermal conductivity of liquid [W/(m*K)]
rho_l_R245fa_sample = [1404.0, 1352.2, 1297.0, 1237.0, 1170.3]; % density of liquid [kg/m^3]
rho_v_R245fa_sample = 1./[0.30756, 0.13972, 0.07103, 0.03922, 0.02295]; % density of vapor [kg/m^3]
mu_l_R245fa_sample = [581.6, 431.5, 330.1, 256.7, 200.3]*1e-6; % viscosity of liquid [Pa*s]
mu_v_R245fa_sample = [9.47, 10.16, 10.87, 11.61, 12.43]*1e-6; % viscosity of vapor [Pa*s]
sigma_R245fa_sample = [17.25, 14.69, 12.13, 9.59, 7.13]*1e-3; % surface tension of liquid [N/m]
h_lv_R245fa = [200.00, 226.20, 253.24, 281.26, 310.50]*1e3; % latent heat of vaporization [J/kg]

% Average temperature
T_avg = (T_a_in + T_a_out)/2; % [C]

% Compute lagrange interpolating polynomials and compute saturation
% properties at average temperature T_avg

if water == 1 && R245fa == 0
    k_l = lagrangepoly(T_water_sample, k_l_water_sample, T_avg); % thermal conductivity of liquid [W/(m*K)]
    rho_l = lagrangepoly(T_water_sample, rho_l_water_sample, T_avg); % density of liquid [kg/m^3]
    rho_v = lagrangepoly(T_water_sample, rho_v_water_sample, T_avg); % density of vapor [kg/m^3]
    mu_l = lagrangepoly(T_water_sample, mu_l_water_sample, T_avg); % viscosity of liquid [Pa*s]
    mu_v = lagrangepoly(T_water_sample,  mu_v_water_sample, T_avg); % viscosity of vapor [Pa*s]
    %{
    Psatfit = coeffvalues(fit(1./T_water_sample', P_water_sample', 'exp1')); % exponential fit form: P = a*e^(b/T)
    C0 = log(Psatfit(1));
    C1 = -1*Psatfit(2);
    %}
    C1 = log(P_water_sample(1)/P_water_sample(end))*(T_water_sample(1)*T_water_sample(end))/(T_water_sample(1) - T_water_sample(end));
    C0 = log(P_water_sample(end))+(C1/T_water_sample(end));
    P_l = exp(C0 - C1/T_avg);
    sigma_poly = lagrangepoly(T_water_sample,  sigma_water_sample);
    sigma = polyval(sigma_poly, T_avg); % surface tension of liquid [N/m]
    h_lv = lagrangepoly(T_water_sample, h_lv_water, T_avg); % latent heat of vaporization [J/kg]
    MW = 18.01528e-3; % molecular weight [kg/mol]
elseif R245fa == 1 && water == 0
    k_l = lagrangepoly(T_R245fa_sample, k_l_R245fa_sample, T_avg); % thermal conductivity of liquid [W/(m*K)]
    rho_l = lagrangepoly(T_R245fa_sample, rho_l_R245fa_sample, T_avg); % density of liquid [kg/m^3]
    rho_v = lagrangepoly(T_R245fa_sample, rho_v_R245fa_sample, T_avg); % density of vapor [kg/m^3]
    mu_l = lagrangepoly(T_R245fa_sample, mu_l_R245fa_sample, T_avg); % viscosity of liquid [Pa*s]
    mu_v = lagrangepoly(T_R245fa_sample,  mu_v_R245fa_sample, T_avg); % viscosity of vapor [Pa*s]
    %{
    Psatfit = coeffvalues(fit(1./T_R245fa_sample', P_R245fa_sample', 'exp1')); % exponential fit form: P = a*e^(b/T)
    C0 = log(Psatfit(1));
    C1 = -1*Psatfit(2);
    %}
    sigma_poly = lagrangepoly(T_R245fa_sample,  sigma_R245fa_sample);
    sigma = polyval(sigma_poly, T_avg); % surface tension of liquid [N/m]
    h_lv = lagrangepoly(T_R245fa_sample, h_lv_R245fa, T_avg); % latent heat of vaporization [J/kg]
    MW = 100.21e-3; % molecular weight [kg/mol]
else
    error('Select one and only one fluid')
end

% P_test = Psatfunc(C0,C1,[0.01, 20:20:80]) % for troubleshooting Psat fit

%% Task [2] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (i) Guess qdot
qdot_guess = 1;
TOL = 1e-3; % tolerance to end iteration
maxit = 100; % max iteration limit

for i = 1:maxit
    % (ii) Compute mdot, D_v, A_v, L_eff, V_v, Re_v, f_v, Ustar_e, Ustar_c,
    % A_c, A_e, and A_w. Compute kappa if necessary.
    mdot = qdot_guess/h_lv; % mean vapor mass flow rate
    D_v = D_p - 2*t_w; % hydraulic diameter of vapor flow passage
    A_v = 0.25*pi*D_v^2; % cross sectional area of vapor flow passage
    L_eff = 0.5*L_evap + L_adiab + 0.5*L_cond; % vapor flow path effective length
    V_v = mdot/(rho_v*A_v); % mean vapor velocity
    Re_v = rho_v*V_v*D_v/mu_v; % Reynolds number of vapor
    if Re_v < 2300
        f_v = 16/Re_v; % friction factor
    else
        f_v = 0.023*Re_v^-0.2; % friction factor
    end
    k_w = omega*k_l + (1 - omega)*k_m; % wick effectice thermal conductivity
    Ustar_e = 1/((t_w/k_w) + (1/(h_a_in*A_a_in_A_e))); % effective heat transfer coefficient in evaporator
    Ustar_c = 1/((t_w/k_w) + (1/(h_a_out*A_a_out_A_c))); % effective heat transfer coefficient in condenser
    A_e = pi*D_v*L_evap; % surface area of evaporator vapor flow passage
    A_c = pi*D_v*L_cond; % surface area of condenser vapor flow passage
    A_w = 0.25*pi*(D_p^2 - D_v^2); % cross sectional area of wick
    
    % (iii) Solve for evaporator interface radius of curvature r_e
    r_e = 2*sigma/(4*f_v*(L_eff/D_v)*0.5*rho_v*V_v^2 + mu_l*L_eff*mdot/(kappa*rho_l*A_w));
    if r_e < r_p/cosd(theta)
        error('No solution possible for the specificied heat input qdot')
    end
    
    % (iv) Solve for T2 and Compute P2
    T2 = (qdot_guess + Ustar_c*A_c*T_a_out)/(Ustar_c*A_c);
    P2 = Psatfunc(C0,C1,T2); % saturation pressure at T2
    
    % (v) Compute P1
    P1 = P2 + 4*f_v*(L_eff/D_v)*0.5*rho_v*V_v^2;
    
    % (vi) Compute T1
    T1_range = 0:0.0001:80; % [C] Select smaller step size for accuracy
    R_universal = 8.31446261815324; % [J/(K*mol)]
    R = R_universal/MW; % [J/(K*kg)]
    P1_range = Psatfunc(C0,C1,T1_range).*exp(-2*sigma./(r_e*rho_l*R*(T1_range + 273.15)));
    [ClosestValueDifference, index_ClosestValueDifference] = min(abs(P1_range - P1));
    T1 = T1_range(index_ClosestValueDifference);
    
    % (vii) Solve for new value of qdot and compute error
    qdot_new = Ustar_e*A_e*(T_a_in - T1);
    E = abs(qdot_new - qdot_guess);
    
    if E < TOL
        break
    end
    
    if Re_v < 2300
        der_r_e = -2*sigma/((32*mu_v*L_eff/(D_v^2*h_lv*rho_v*A_v) ...
            + mu_l*L_eff/(kappa*rho_l*A_w*h_lv))*qdot_guess^2);
    else
        der_r_e = -2*sigma*(1.8*0.0046*mu_v^0.2*L_eff*qdot_guess^0.8/(D_v^1.2*h_lv^1.8*rho_v*A_v^1.8) ...
            + mu_l*L_eff/(kappa*rho_l*A_w*h_lv))/(0.0046*mu_v^0.2*L_eff*qdot_guess^1.8/(D_v^1.2*h_lv^1.8*rho_v*A_v^1.8) ...
            + mu_l*L_eff*qdot_guess/(kappa*rho_l*A_w*h_lv))^2;
    end
    
    qdot_guess = qdot_guess - r_e/der_r_e; % Newton-Raphson scheme

end

qdot = qdot_guess;

%% Functions
function [P] = Psatfunc(C0,C1,T)
P = exp(C0 - C1./T); % saturation pressure [Pa]
end

function [T] = Tsatfunc(C0,C1,P)
T = C1./(C0 - log(P)); % saturation temprature [C]
end