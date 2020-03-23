clear; clc;

%% Inputs 
load FTISxprt-20200311_flight1.mat
addpath('../Custom_functions')
addpath('../Numerical_model')

T0 = 288.15;
h = 10400;
p0 = 101325;
M = 0.0289644;     % Molar mass air
R = 287;           % Universal gass constant
a = -6.5/1000;     % deg/m 
g = 9.80665;

%% Convert flight data into correct units

sz = size(flightdata.time.data);
for index = 1:sz(2)
    flightdata.Dadc1_tas.data(index) = toms(flightdata.Dadc1_tas.data(index));                  % True air speed [kts] to [m/s]
    flightdata.vane_AOA.data(index) = deg2rad(flightdata.vane_AOA.data(index));                 % angle of attack [deg] to [rad]
    flightdata.Ahrs1_Pitch.data(index) = deg2rad(flightdata.Ahrs1_Pitch.data(index));           % theta [deg] to [rad]
    flightdata.Ahrs1_Roll.data(index) = deg2rad(flightdata.Ahrs1_Roll.data(index));             % phi [deg] to [rad]
    flightdata.Ahrs1_bRollRate.data(index) = deg2rad(flightdata.Ahrs1_bRollRate.data(index));   % p [deg/s] to [rad/s]
    flightdata.Ahrs1_bPitchRate.data(index) = deg2rad(flightdata.Ahrs1_bPitchRate.data(index)); % q [deg/s] to [rad/s]
    flightdata.Ahrs1_bYawRate.data(index) = deg2rad(flightdata.Ahrs1_bYawRate.data(index));     % r [deg/s] to [rad/s]
    flightdata.delta_e.data(index) = deg2rad(flightdata.delta_e.data(index));     % Elevator deflection [deg] to [rad]
    flightdata.delta_a.data(index) = deg2rad(flightdata.delta_a.data(index));     % Elevator deflection [deg] to [rad]
    flightdata.delta_r.data(index) = deg2rad(flightdata.delta_r.data(index));     % Elevator deflection [deg] to [rad]
    flightdata.Dadc1_alt.data(index) = tom(flightdata.Dadc1_alt.data(index));     % Altitude [ft] to [m]
end

%% Eigenmotion plots and model matching
sym_plots = false;
unchanged_plots = true;
changed_plots = false;
unsym_plots = true;

m = 6065.3; % at index 31430
fuel_flow = 0.132; % Combined average fuel flow [kg/s]
unchanged_model_inputs = {0.06290,-8.79415,-0.10260,-0.71085,+0.23760};  % [Cm_u,Cm_q,Cl_beta,Cl_p,Cl_r]
Cm_u    = +0.07090; %+0.07090
Cm_q    = -4.39415; %-4.39415
Cl_beta = -0.1460;
Cl_p    = -0.93085;
Cl_r    = +0.1260;
model_matching_inputs = {Cm_u,Cm_q,Cl_beta,Cl_p,Cl_r}; % [Cm_u,Cm_q,Cl_beta,Cl_p,Cl_r]

if sym_plots
    %% Data Analysis - Short Period
    t_dr = (60*55+44.8)/0.1;        % start time SP index 33350
    step = 100;
    phi_0 = flightdata.Ahrs1_Roll.data(t_dr);
    p_0 = flightdata.Ahrs1_bRollRate.data(t_dr);
    r_0 = flightdata.Ahrs1_bYawRate.data(t_dr);
    start_inputs = {flightdata.Dadc1_alt.data(t_dr), flightdata.Dadc1_tas.data(t_dr), flightdata.vane_AOA.data(t_dr), flightdata.Ahrs1_Pitch.data(t_dr), m - fuel_flow*(t_dr-31430)/10, phi_0, p_0, r_0}; 
    % [hp_0, V_0,alpha_0,th_0,m, phi_0, p_0, r_0]
    step_input_sym = flightdata.delta_e.data(t_dr:t_dr+step);
    step_input_asym_da = flightdata.delta_a.data(t_dr:t_dr+step);
    step_input_asym_dr = flightdata.delta_r.data(t_dr:t_dr+step);
    [t, symmetric_response, Ya0] = numerical_model(start_inputs, unchanged_model_inputs, step_input_sym, step_input_asym_da, step_input_asym_dr, step/10);
    
    if unchanged_plots
        figure('Renderer', 'painters', 'Position', [10 10 600 400])
        subplot(2,2,1)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_Pitch.data(t_dr:t_dr + step), t, symmetric_response(:,3),'k','linewidth',1.3)
        title('Pitch','fontsize',14)
        xlabel('time (s)','fontsize',14)
        ylabel('angle (rad)','fontsize',14)

        subplot(2,2,2)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_bPitchRate.data(t_dr:t_dr + step),t, symmetric_response(:,4),'k','linewidth',1.3)
        title('Pitch Rate','fontsize',14)
        xlabel('time (s)','fontsize',14)
        ylabel('angle per sec (rad/s)','fontsize',14)

        subplot(2,2,3)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.vane_AOA.data(t_dr:t_dr + step),t, symmetric_response(:,2),'k','linewidth',1.3)
        title('Angle of Attack','fontsize',14)
        xlabel('time (s)','fontsize',14)
        ylabel('angle (rad)','fontsize',14)

        subplot(2,2,4)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Dadc1_tas.data(t_dr:t_dr + step),t, symmetric_response(:,1),'k','linewidth',1.3)
        title('True airspeed','fontsize',14)
        xlabel('time (s)','fontsize',14)
        ylabel('velocity (m/s)','fontsize',14)
        saveas(gcf,'sym_short_period.png')
        sgtitle('Short Period')
    end
    
    %% Data Analysis - Phugoid
    t_dr = (60*52+23)/0.1;      % start time Phugoid index 31430
    step = 1200;
    start_inputs = {flightdata.Dadc1_alt.data(t_dr), flightdata.Dadc1_tas.data(t_dr), flightdata.vane_AOA.data(t_dr), flightdata.Ahrs1_Pitch.data(t_dr), m - fuel_flow*(t_dr-31430)/10, phi_0, p_0, r_0};
    step_input_sym = flightdata.delta_e.data(t_dr:t_dr+step);
    step_input_asym_da = flightdata.delta_a.data(t_dr:t_dr+step);
    step_input_asym_dr = flightdata.delta_r.data(t_dr:t_dr+step);
    [t, symmetric_response, Ya] = numerical_model(start_inputs,unchanged_model_inputs, step_input_sym, step_input_asym_da, step_input_asym_dr, step/10);
    [t, symmetric_response2, Ya] = numerical_model(start_inputs,model_matching_inputs, step_input_sym, step_input_asym_da, step_input_asym_dr, step/10);
    
    if unchanged_plots
        %Unchanged model coefficients
        figure(2)
        sgtitle('Phugoid')
        subplot(2,2,1)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_Pitch.data(t_dr:t_dr + step), t, symmetric_response(:,3),'k','linewidth',1.3)
        title('Pitch angle','fontsize',14)
        xlabel('time (s)','fontsize',14)
        ylabel('angle (rad)','fontsize',14)

        subplot(2,2,2)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_bPitchRate.data(t_dr:t_dr + step),t, symmetric_response(:,4),'k','linewidth',1.3)
        title('Pitch Rate','fontsize',14)
        xlabel('time (s)','fontsize',14)
        ylabel('angle per sec (rad/s)','fontsize',14)

        subplot(2,2,3)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.vane_AOA.data(t_dr:t_dr + step),t, symmetric_response(:,2),'k','linewidth',1.3)
        title('Angle of Attack','fontsize',14)
        xlabel('time (s)','fontsize',14)
        ylabel('angle (rad)','fontsize',14)

        subplot(2,2,4)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Dadc1_tas.data(t_dr:t_dr + step),t, symmetric_response(:,1),'k','linewidth',1.3)
        title('True airspeed','fontsize',14)
        xlabel('time (s)','fontsize',14)
        ylabel('velocity (m/s)','fontsize',14)
    end
    
    if changed_plots
        %Changed model coefficients
        figure(3)
        sgtitle('Phugoid with changed coefficients')
        subplot(2,2,1)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_Pitch.data(t_dr:t_dr + step), t, symmetric_response2(:,3),'k','linewidth',1.3)
        title('Pitch angle','fontsize',14)
        xlabel('time (s)','fontsize',14)
        ylabel('angle (rad)','fontsize',14)

        subplot(2,2,2)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_bPitchRate.data(t_dr:t_dr + step),t, symmetric_response2(:,4),'k','linewidth',1.3)
        title('Pitch Rate','fontsize',14)
        xlabel('time (s)','fontsize',14)
        ylabel('angle per sec (rad/s)','fontsize',14)

        subplot(2,2,3)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.vane_AOA.data(t_dr:t_dr + step),t, symmetric_response2(:,2),'k','linewidth',1.3)
        title('Angle of Attack','fontsize',14)
        xlabel('time (s)','fontsize',14)
        ylabel('angle (rad)','fontsize',14)

        subplot(2,2,4)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Dadc1_tas.data(t_dr:t_dr + step),t, symmetric_response2(:,1),'k','linewidth',1.3)
        title('True airspeed','fontsize',14)
        xlabel('time (s)','fontsize',14)
        ylabel('velocity (m/s)','fontsize',14)
        
    end
end


if unsym_plots
    %% Data Analysis - Dutch Roll
    t_dr = (60*58+30)/0.1;        % start time Dutch roll index 35100
    step = 350;
    phi_0 = flightdata.Ahrs1_Roll.data(t_dr);
    p_0 = flightdata.Ahrs1_bRollRate.data(t_dr);
    r_0 = flightdata.Ahrs1_bYawRate.data(t_dr);
    start_inputs = {flightdata.Dadc1_alt.data(t_dr), flightdata.Dadc1_tas.data(t_dr), flightdata.vane_AOA.data(t_dr), flightdata.Ahrs1_Pitch.data(t_dr), m - fuel_flow*(t_dr-31430)/10, phi_0, p_0, r_0}; 
    step_input_sym = flightdata.delta_e.data(t_dr:t_dr+step);
    step_input_asym_da = flightdata.delta_a.data(t_dr:t_dr+step);
    step_input_asym_dr = flightdata.delta_r.data(t_dr:t_dr+step);
    [t, symmetric_response, Ya0] = numerical_model(start_inputs, unchanged_model_inputs, step_input_sym, step_input_asym_da, step_input_asym_dr, step/10);
    [t, symmetric_response, Ya0_2] = numerical_model(start_inputs, model_matching_inputs, step_input_sym, step_input_asym_da, step_input_asym_dr, step/10);

    if unchanged_plots
 
        figure(4)
        sgtitle('Dutch Roll')

        subplot(2,2,1)
        plot(t, Ya0(:,1),'k','linewidth',1.3)
        title('Sideslip angle')
        xlabel('time (s)','fontsize',14)
        ylabel('velocity (m/s)','fontsize',14)

        subplot(2,2,2)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_Roll.data(t_dr:t_dr + step), t, Ya0(:,2),'k','linewidth',1.3)
        title('Roll')
        xlabel('time (s)','fontsize',14)
        ylabel('angle (rad)','fontsize',14)

        subplot(2,2,3)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_bRollRate.data(t_dr:t_dr + step), t, Ya0(:,3),'k','linewidth',1.3)
        title('Roll Rate')
        xlabel('time (s)','fontsize',14)
        ylabel('angle per sec (rad/s)','fontsize',14)

        subplot(2,2,4)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_bYawRate.data(t_dr:t_dr + step),t, Ya0(:,4),'k','linewidth',1.3)
        title('Yaw Rate')
        xlabel('time (s)','fontsize',14)
        ylabel('angle (rad)','fontsize',14)
    end
    
    if changed_plots
        figure(5)
        sgtitle('Dutch Roll with changed coefficients')

        subplot(2,2,1)
        plot(t, Ya0_2(:,1),'k','linewidth',1.3)
        title('Sideslip angle')
        xlabel('time (s)','fontsize',14)
        ylabel('velocity (m/s)','fontsize',14)

        subplot(2,2,2)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_Roll.data(t_dr:t_dr + step), t, Ya0_2(:,2),'k','linewidth',1.3)
        title('Roll')
        xlabel('time (s)','fontsize',14)
        ylabel('angle (rad)','fontsize',14)

        subplot(2,2,3)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_bRollRate.data(t_dr:t_dr + step), t, Ya0_2(:,3),'k','linewidth',1.3)
        title('Roll Rate')
        xlabel('time (s)','fontsize',14)
        ylabel('angle per sec (rad/s)','fontsize',14)

        subplot(2,2,4)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_bYawRate.data(t_dr:t_dr + step),t, Ya0_2(:,4),'k','linewidth',1.3)
        title('Yaw Rate')
        xlabel('time (s)','fontsize',14)
        ylabel('angle (rad)','fontsize',14)
    end
    
    %% Data Analysis - Aperiodic Roll
    t_dr = (60*56+30)/0.1;        % start time Aperiodic Roll start index 33900
    step = 1200;
    phi_0 = flightdata.Ahrs1_Roll.data(t_dr);
    p_0 = flightdata.Ahrs1_bRollRate.data(t_dr);
    r_0 = flightdata.Ahrs1_bYawRate.data(t_dr);
    start_inputs = {flightdata.Dadc1_alt.data(t_dr), flightdata.Dadc1_tas.data(t_dr), flightdata.vane_AOA.data(t_dr), flightdata.Ahrs1_Pitch.data(t_dr), m - fuel_flow*(t_dr-31430)/10, phi_0, p_0, r_0}; 
    step_input_sym = flightdata.delta_e.data(t_dr:t_dr+step);
    step_input_asym_da = flightdata.delta_a.data(t_dr:t_dr+step);
    step_input_asym_dr = flightdata.delta_r.data(t_dr:t_dr+step);
    [t, symmetric_response, Ya0] = numerical_model(start_inputs, unchanged_model_inputs, step_input_sym, step_input_asym_da, step_input_asym_dr, step/10);
    [t, symmetric_response, Ya0_2] = numerical_model(start_inputs, model_matching_inputs, step_input_sym, step_input_asym_da, step_input_asym_dr, step/10);
    
    if unchanged_plots
        figure(6)
        sgtitle('Aperiodic Roll')
        subplot(2,2,1)
        plot(t, Ya0(:,1),'k','linewidth',1.3)
        title('Sideslip angle')
        xlabel('time (s)','fontsize',14)
        ylabel('velocity (m/s)','fontsize',14)

        subplot(2,2,2)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_Roll.data(t_dr:t_dr + step), t, Ya0(:,2),'k','linewidth',1.3)
        title('Roll')
        xlabel('time (s)','fontsize',14)
        ylabel('angle (rad)','fontsize',14)

        subplot(2,2,3)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_bRollRate.data(t_dr:t_dr + step), t, Ya0(:,3),'k','linewidth',1.3)
        title('Roll Rate')
        xlabel('time (s)','fontsize',14)
        ylabel('angle per sec (rad/s)','fontsize',14)

        subplot(2,2,4)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_bYawRate.data(t_dr:t_dr + step),t, Ya0(:,4),'k','linewidth',1.3)
        title('Yaw Rate')
        xlabel('time (s)','fontsize',14)
        ylabel('angle (rad)','fontsize',14)
    
    end
    
    if changed_plots
        figure(7)
        sgtitle('Aperiodic Roll with changed coefficients')
        subplot(2,2,1)
        plot(t, Ya0_2(:,1),'k','linewidth',1.3)
        title('Sideslip angle')
        xlabel('time (s)','fontsize',14)
        ylabel('velocity (m/s)','fontsize',14)

        subplot(2,2,2)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_Roll.data(t_dr:t_dr + step), t, Ya0_2(:,2),'k','linewidth',1.3)
        title('Roll')
        xlabel('time (s)','fontsize',14)
        ylabel('angle (rad)','fontsize',14)

        subplot(2,2,3)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_bRollRate.data(t_dr:t_dr + step), t, Ya0_2(:,3),'k','linewidth',1.3)
        title('Roll Rate')
        xlabel('time (s)','fontsize',14)
        ylabel('angle per sec (rad/s)','fontsize',14)

        subplot(2,2,4)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_bYawRate.data(t_dr:t_dr + step),t, Ya0_2(:,4),'k','linewidth',1.3)
        title('Yaw Rate')
        xlabel('time (s)','fontsize',14)
        ylabel('angle (rad)','fontsize',14)
        
    end
      
    %% Data Analysis - Spiral start index 

    t_dr = (60*63)/0.1;        % start time Spiral start index 37800
    step = 1200;
    phi_0 = flightdata.Ahrs1_Roll.data(t_dr);
    p_0 = flightdata.Ahrs1_bRollRate.data(t_dr);
    r_0 = flightdata.Ahrs1_bYawRate.data(t_dr);
    start_inputs = {flightdata.Dadc1_alt.data(t_dr), flightdata.Dadc1_tas.data(t_dr), flightdata.vane_AOA.data(t_dr), flightdata.Ahrs1_Pitch.data(t_dr), m - fuel_flow*(t_dr-31430)/10, phi_0, p_0, r_0}; 
    step_input_sym = flightdata.delta_e.data(t_dr:t_dr+step);
    step_input_asym_da = flightdata.delta_a.data(t_dr:t_dr+step);
    step_input_asym_dr = flightdata.delta_r.data(t_dr:t_dr+step);
    [t, symmetric_response, Ya0] = numerical_model(start_inputs, unchanged_model_inputs, step_input_sym, step_input_asym_da, step_input_asym_dr, step/10);
    [t, symmetric_response, Ya0_2] = numerical_model(start_inputs, model_matching_inputs, step_input_sym, step_input_asym_da, step_input_asym_dr, step/10);
    
    if unchanged_plots
        figure(8)
        sgtitle('Spiral')
        subplot(2,2,1)
        plot(t, Ya0(:,1),'k','linewidth',1.3)
        title('Sideslip angle')
        xlabel('time (s)','fontsize',14)
        ylabel('velocity (m/s)','fontsize',14)

        subplot(2,2,2)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_Roll.data(t_dr:t_dr + step), t, Ya0(:,2),'k','linewidth',1.3)
        title('Roll')
        xlabel('time (s)','fontsize',14)
        ylabel('angle (rad)','fontsize',14)

        subplot(2,2,3)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_bRollRate.data(t_dr:t_dr + step), t, Ya0(:,3),'k','linewidth',1.3)
        title('Roll Rate')
        xlabel('time (s)','fontsize',14)
        ylabel('angle per sec (rad/s)','fontsize',14)

        subplot(2,2,4)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_bYawRate.data(t_dr:t_dr + step),t, Ya0(:,4),'k','linewidth',1.3)
        title('Yaw Rate')
        xlabel('time (s)','fontsize',14)
        ylabel('angle (rad)','fontsize',14)
    
    end
    
    if changed_plots
        figure(9)
        sgtitle('Spiral with changed coefficients')
        subplot(2,2,1)
        plot(t, Ya0_2(:,1),'k','linewidth',1.3)
        title('Sideslip angle')
        xlabel('time (s)','fontsize',14)
        ylabel('velocity (m/s)','fontsize',14)

        subplot(2,2,2)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_Roll.data(t_dr:t_dr + step), t, Ya0_2(:,2),'k','linewidth',1.3)
        title('Roll')
        xlabel('time (s)','fontsize',14)
        ylabel('angle (rad)','fontsize',14)

        subplot(2,2,3)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_bRollRate.data(t_dr:t_dr + step), t, Ya0_2(:,3),'k','linewidth',1.3)
        title('Roll Rate')
        xlabel('time (s)','fontsize',14)
        ylabel('angle per sec (rad/s)','fontsize',14)

        subplot(2,2,4)
        plot(flightdata.time.data(t_dr:t_dr+step)-t_dr*0.1-8.9, flightdata.Ahrs1_bYawRate.data(t_dr:t_dr + step),t, Ya0_2(:,4),'k','linewidth',1.3)
        title('Yaw Rate')
        xlabel('time (s)','fontsize',14)
        ylabel('angle (rad)','fontsize',14)
    end
end