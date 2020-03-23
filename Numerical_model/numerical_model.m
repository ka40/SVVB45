function [t,symmetric_response, Ya] = numerical_model(start_inputs,model_matching_inputs, step_input_sym, step_input_asym_da, step_input_asym_dr, stop_time)
%NUMERICAL_MODEL Return the output of 
%   Detailed explanation goes here
% Citation 550 - Linear simulation

% xcg = 0.25*c

% Stationary flight condition

hp_0    = 2000;      	  % pressure altitude in the stationary flight condition [m]
V_0     = 86;            % true airspeed in the stationary flight condition [m/sec]
alpha_0 = 0;       	  % angle of attack in the stationary flight condition [rad]
th_0    = 0;       	  % pitch angle in the stationary flight condition [rad]

% Aircraft mass
m      = 8200;         	  % mass [kg]

% aerodynamic properties
e           = 0.8;            % Oswald factor [ ]
CD0         = 0.04;            % Zero lift drag coefficient [ ]
CL_alpha    = 5.084;            % Slope of CL-alpha curve [ ]

% Longitudinal stability

% Flight data:
% C_m_alpha= - 0,0299
% C_m_delta= -0,0719
% 
% Ref data:
% C_m_alpha = -0,0258
% C_m_delta= -0,0565

Cm_alpha     = -0.5626;            % longitudinal stabilty [ ] -0.5626
Cm_delta_e   = -1.1642;            % elevator effectiveness [ ] -1.1642
[hp_0, V_0,alpha_0,th_0,m, phi_0, p_0, r_0] = start_inputs{:};
% Aircraft geometry

S      = 30.00;	          % wing area [m^2]
Sh     = 0.2*S;           % stabiliser area [m^2]
Sh_S   = Sh/S;	          % [ ]
lh     = 0.71*5.968;      % tail length [m]
c      = 2.0569;	      % mean aerodynamic cord [m]
lh_c   = lh/c;	          % [ ]
b      = 15.911;          % wing span [m]
bh     = 5.791;	          % stabilser span [m]
A      = b^2/S;           % wing aspect ratio [ ]
Ah     = bh^2/Sh;         % stabilser aspect ratio [ ]
Vh_V   = 1;               % [ ]
ih     = -2*pi/180;       % stabiliser angle of incidence [rad]

% Constant values concerning atmosphere and gravity

rho_0   = 1.2250;          % air density at sea level [kg/m^3] 
lambda = -0.0065;         % temperature gradient in ISA [K/m]
Temp_0  = 288.15;          % temperature at sea level in ISA [K]
R      = 287.05;          % specific gas constant [m^2/sec^2K]
g      = 9.81;            % [m/sec^2] (gravity constant)

rho    = rho_0*((1+(lambda*hp_0/Temp_0)))^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
W      = m*g;				                        % [N]       (aircraft weight)

% Constant values concerning aircraft inertia

mu_c_0    = m/(rho*S*c);
mu_b_0    = m/(rho*S*b);
Kxx2      = 0.019;
Kyy2      = 1.25*1.114;
Kzz2      = 0.042;
Kxz       = 0.002;


% Aerodynamic constants

Cmac   = 0;                     % Moment coefficient about the aerodynamic centre [ ]
CNwa   = CL_alpha;   		        % Wing normal force slope [ ]
CNha   = 2*pi*Ah/(Ah+2);        % Stabiliser normal force slope [ ]
depsda = 4/(A+2);               % Downwash gradient [ ]

% Lift and drag coefficient

CL = 2*W/(rho*V_0^2*S);               % Lift coefficient [ ]
CD = CD0 + (CL_alpha*alpha_0)^2/(pi*A*e);  % Drag coefficient [ ]

% Stabiblity derivatives

Cx_0        = W*sin(th_0)/(0.5*rho*V_0^2*S);
Cx_u        = -0.02792;
Cx_alpha    = -0.47966;
Cx_alphadot = +0.08330;
Cx_q        = -0.28170;
Cx_delta_e  = -0.03728;

Cz_0        = -W*cos(th_0)/(0.5*rho*V_0^2*S);
Cz_u        = -0.37616;
Cz_alpha    = -5.74340;
Cz_alphadot = -0.00350;
Cz_q        = -5.66290;
Cz_delta_e  = -0.69612;

Cm_u        = +0.06990;
Cm_alphadot = +0.17800;
Cm_q        = -8.79415;

Cy_beta      = -0.7500;
Cy_betadot   =  0     ;
Cy_p         = -0.0304;
Cy_r         = +0.8495;
Cy_delta_a   = -0.0400;
Cy_delta_r   = +0.2300;

Cl_beta      = -0.10260;
Cl_p         = -0.71085;
Cl_r         = +0.23760;
Cl_delta_a   = -0.23088;
Cl_delta_r   = +0.03440;

Cn_beta      =  +0.1348;
Cn_betadot   =   0     ;
Cn_p         =  -0.0602;
Cn_r         =  -0.2061;
Cn_delta_a   =  -0.0120;
Cn_delta_r   =  -0.0939;

% Change input parameters

[Cm_u,Cm_q,Cl_beta,Cl_p,Cl_r] = model_matching_inputs{:};

% MATRIX COEFFICIENT CALCULATION
% Symmetric
Us = [-2*mu_c_0*c/(V_0^2), 0                             , 0 , 0                        ;
       0                 , (Cz_alphadot - 2*mu_c_0)*c/V_0, 0 , 0                        ;
       0                 , 0                             , -1, 0                        ;
       0                 , Cm_alphadot*c/V_0             , 0 , -2*mu_c_0*Kyy2*(c/V_0)^2 ]  ;

Vs = [-Cx_u/V_0, -Cx_alpha, -Cz_0, -Cx_q*c/V_0            ;
      -Cz_u/V_0, -Cz_alpha, Cx_0 , -(Cz_q+2*mu_c_0)*c/V_0 ;
      0        , 0        , 0    , -1                     ;
      -Cm_u/V_0, -Cm_alpha, 0    , -Cm_q*c/V_0            ]   ;

Ws = [-Cx_delta_e ;
      -Cz_delta_e ;
      0           ;
      -Cm_delta_e ]   ;
  
% Asymmetric
Ua = [(Cy_betadot - 2*mu_b_0)*b/V_0, 0                           , 0                       , 0                           ;
       0                           , -1                          , 0                       , 0                           ;
       0                           , 0                           , -2*mu_b_0*Kxx2*(b/V_0)^2, 2*mu_b_0*Kxz*(b/V_0)^2      ;
       Cn_betadot*b/V_0            , 0                           , 2*mu_b_0*Kxz*(b/V_0)^2          , -2*mu_b_0*Kzz2*(b/V_0)^2 ]   ;

Va = [-Cy_beta, -CL, -Cy_p*b/(2*V_0), -(Cy_r - 4*mu_b_0)*b/(2*V_0) ;
      0       , 0  , -1             , 0                              ;
      -Cl_beta, 0  , -b*Cl_p/(2*V_0), -b*Cl_r/(2*V_0)                ;
      -Cn_beta, 0  , -b*Cn_p/(2*V_0), -b*Cn_r/(2*V_0)                ]   ;

Wa = [-Cy_delta_a, -Cy_delta_r ;
      -0         , 0           ;
      -Cl_delta_a, -Cl_delta_r ;
      -Cn_delta_a, -Cn_delta_r ]    ; 

  
%MATRIX CALCULATION
% Symmetric matrices for x = [u ; alpha ; theta; q] and u = [delta_e]
As = inv(Us)*Vs   ;
  
Bs = inv(Us)*Ws   ;
  
Cs = eye(4)       ;
  
Ds = zeros(4,1)   ;

% Asymmetric matrices for x = [beta ; phi ; p ; r]
Aa = inv(Ua)*Va   ;

Ba = inv(Ua)*Wa   ;

Ca = eye(4)       ;

Da = zeros(4,2)   ;



%SIMULATION
% Set state space model
sys_sym  = ss(As,Bs,Cs,Ds);       %Create symmetric state-space model
sys_asym = ss(Aa,Ba,Ca,Da);       %Create asymmetric state space object

% Simulate step input
t = [0:0.1:stop_time];
% symmetric_response  = step(sys_sym,t,stepDataOptions('StepAmplitude',step_input_sym));
symmetric_response  = lsim(sys_sym,step_input_sym,t);
symmetric_response(:,1)=symmetric_response(:,1)+V_0;
symmetric_response(:,2)=symmetric_response(:,2)+alpha_0;
symmetric_response(:,3)=symmetric_response(:,3)+th_0;
symmetric_response(:,4)=symmetric_response(:,4);

% Simulate asymmetric step input

delta_a = step_input_asym_da;
delta_r = step_input_asym_dr;
step_input_asym = [delta_a, delta_r];
Ya = lsim(sys_asym,step_input_asym,t);
Ya(:,1)=Ya(:,1);
Ya(:,2)=Ya(:,2)+phi_0;
Ya(:,3)=Ya(:,3)+p_0;
Ya(:,4)=Ya(:,4)+r_0;
% [beta ; phi ; p ; r]

end