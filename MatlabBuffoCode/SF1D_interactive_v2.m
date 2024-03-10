%% Vectors for each time chunk of iteration for plotting
Tplot=[];
Splot=[];
phiplot=[];
wplot=[];
Stotplot=[];
vel_hold=[];
salt_hold=[];
t_hold=[];

%% Time loops
t_i=0;                  %% Initial Time (sec) - most commonly=0
t_f=1175000;           %% Final Time (sec) - length of run
time=[t_i:t_f:t_f];
for k=1:length(time)-1;

%% Input Parameters 
%% Change these as desired - these govern the majority of the results
Ttop=270;%250;                %% Ice top temp (K)
Tbottom=273.15; %272.44;          %% Ocean Temp (K) [MgSO4 12.3ppt-273.0K, 100ppt-271.31K, 282ppt-267.5K NaCl 34ppt-271.2K]
interface_in=0.50;       %% Initial interface depth (m)
Tm=273.15;               %% Melt temperature of pure ice (K)
k_i=2;                   %% Thermal conductivity (ice) W/m*K
k_br=.6;                 %% '                            ' (brine)
k_s=2*10^-9;             %% Diffusivity for Salt
dt=47;                   %% Time step (s)
dz=.01;                  %% Spatial discretization (m)
H=1.00;                  %% Domain size (m)
c_br=3985;               %% Specific heat of seawater (J/K)
c_i=2000;                %% '              ' ice
L=334774;                %% Latent heat of fusion ice<->water (J/Kg)
S_bottom=0;%5;           %% Ocean Salinity (ppt)
rho_i=917;               %% Density of Ice (Kg/m^3)
rho_br=1004;             %% Density of Ocean (used in volume averaging - 1D grav. drainage uses delta S) 34ppt NaCl-1027
                         %  12.3ppt MgSO4-1012, 100pppt-1103, 282ppt-1323
t_start=time(k);                %dont change (change t_i and t_f above)
t_end=time(k+1);                %dont change (change t_i and t_f above)
tf=t_end-t_start;               %dont change (change t_i and t_f above)
STol=0.01;               %% Error Tolerance for Salinity
PhiTol=0.001;            %% Error Tolerance for Porosity
TTol=0.01;               %% Error Tolerance for Temperature
m=2;                     %% Cementation exponent
g=0;%9.8;                   %% Earth Gravity
%g=1.32;                 %% Europa Gravity (m/s^2)
mu=1.88*10^-3;           %% Viscosity of Ocean
phi_c=0.06823;           %% Critical Porosity


%% BC=1 for Dirilecht Boundary Condition, BC=2 Neumann BC
%% For icy worlds with very low surface temps Neumann BC is typically preferable for 
%% stability reasons. Also the model will likely be innacurate for very shallow ice
%% forming near a turbulent surface in contact with vacuum. Easiest to specify a starting
%% shell interface ~10m, and use Neumann condition. Can change thermal gradient 'T_grad' below
%% to accomodate various thermal profiles in overlying ice

%BC=1;
BC=1;

Depth=[0:-1:-(H/dz)];

%% Initial Condition Vectors for Phi and T and w and S

if t_start==0
    T_initial=0*[0:dz:H]+Tbottom;
    S_initial=0*[0:dz:H]+S_bottom;
    phi_initial=0*[0:dz:H]+1;
    w_initial=0*[0:dz:H];
    history=[];
    Interface_depth=interface_in;
    T_grad=(Tbottom-Ttop)/Interface_depth;
else
    T_initial=IC(:,1)';
    S_initial=IC(:,2)';
    phi_initial=IC(:,3)';
    %w_initial=IC(:,4)';
end

%% Vectors for plotting
Temperature=[];

Liquid_Fraction=[];

Salinity=[];

Darcy_Velocity=[];



%% Looping Over Time
for n=1:tf/dt
    if n==1;

        [NewTemp,NewPhi,NewS,Neww]=one_D_adv_ARD_interactive_v2(...
        k_i,k_br,k_s,dt,dz,H,c_i,c_br,L,STol,TTol,PhiTol,phi_initial,...
        T_initial,S_initial,S_bottom,rho_i,...
        rho_br,Tm,m,w_initial,n,Ttop,Tbottom,T_grad,BC,g,rho_br);

        Temperature=NewTemp;

        Liquid_Fraction=NewPhi;

        Salinity=NewS;

        Darcy_Velocity=Neww;

    else
        [NewTemp,NewPhi,NewS,Neww]=one_D_adv_ARD_interactive_v2(...
        k_i,k_br,k_s,dt,dz,H,c_i,c_br,L,STol,TTol,PhiTol,Liquid_Fraction',...
        Temperature',Salinity',S_bottom,rho_i,...
        rho_br,Tm,m,Darcy_Velocity',n,Ttop,Tbottom,T_grad,BC,g,rho_br);
        
        %% catalogueing ice that has dropped below the critical porosity
        crit_interface=0;
        for i=1:H/dz+1
            if NewPhi(i)<=phi_c;
                crit_interface=i;%crit_interface+1;
            else
            %break
            end
        end
        
        if crit_interface==0;
            Temperature=NewTemp;
            Liquid_Fraction=NewPhi;
            Salinity=NewS;
            Darcy_Velocity=Neww;
        else
            % Temperature=[NewTemp(crit_interface+1:end); Tbottom*(1+0*[1:crit_interface])'];
            % Liquid_Fraction=[NewPhi(crit_interface+1:end); 1*(1+0*[1:crit_interface])'];
            % Salinity=[NewS(crit_interface+1:end); S_bottom*(1+0*[1:crit_interface])'];
            % Darcy_Velocity=[Neww(crit_interface+1:end); Neww(end)*(1+0*[1:crit_interface])'];
            Temperature=NewTemp;
            Liquid_Fraction=NewPhi;
            Salinity=NewS;
            Darcy_Velocity=Neww;
            history=[history; NewTemp(1:crit_interface) NewS(1:crit_interface) ...
                NewPhi(1:crit_interface)];
            
            Interface_depth=Interface_depth+dz*crit_interface;
            T_grad=(Tbottom-Ttop)/Interface_depth;
            
            vel_hold=[vel_hold; Interface_depth, n*dt];
            salt_hold=[salt_hold, NewS(crit_interface)*NewPhi(crit_interface)];
            t_hold=[t_hold, n*dt];
        end
        
        Tplot=[Tplot Temperature];
    end
end

%% Catalogues properties of the mushy layer at times determined by 'block' in original 'time' array (line 16)
Tplot=[Tplot Temperature];
Splot=[Splot Salinity];
phiplot=[phiplot Liquid_Fraction];
wplot=[wplot Darcy_Velocity];
Stotplot=[Stotplot Liquid_Fraction.*Salinity];


% save(strcat('HECCtest_Europa_',num2str(interface_in),'m_salinity_',num2str(S_bottom)...
%     ,'_Tgrad_',num2str(T_grad),'_time_',num2str(t_start),'_to_',num2str(t_end),...
%     '.mat'));

%% Initial conditions for next block of iteration
IC=[Temperature Salinity Liquid_Fraction]

%% Catalogued ice properties [Temperature at time of freeze out: Salinity at time of freeze out: Porosity at time of freeze out]
[history]
end

%% Bulk salinity vs depth
% figure
% plot(history(:,2).*history(:,3),[-interface_in-dz:-dz:-interface_in-dz*length(history(:,2).*history(:,3))])
% title('Bulk Salinity vs. Depth')
% xlabel('Bulk Salinity (ppt)')
% ylabel('Depth (m)')
% 
% figure
% plot(history(:,1),[-interface_in-dz:-dz:-interface_in-dz*length(history(:,1))])
% title('Temperature vs. Depth')
% xlabel('Temperature in K')
% ylabel('Depth (m)')

% figure
% plot(Tplot(:,1),[0:-dz:-H])
% title('Temperature vs. Depth')
% xlabel('Temperature in K')
% ylabel('Depth (m)')

figure
plot([0:1:24999], Tplot(10,:))
title('Temperature vs. Depth')
xlabel('Temperature in K')
ylabel('Depth (m)')

% figure
% plot(Splot,[0:-dz:-1])
% title('Salinity vs. Depth')
% xlabel('Salinity in ppt')
% ylabel('Depth (m)')

% figure
% plot(phiplot,[0:-dz:-1])
% title('Phi vs. Depth')
% xlabel('Phi')
% ylabel('Depth (m)')
% 
% figure
% plot(wplot,[0:-dz:-1])
% title('Velocity vs. Depth')
% xlabel('Velocity in K')
% ylabel('Depth (m)')

% velocity=[];
% for i=2:length(vel_hold)
%     velocity(i-1)=(vel_hold(i,1)-vel_hold(i-1,1))/(vel_hold(i,2)-vel_hold(i-1,2));
% end

% figure
% plot(vel_hold(2:end,2),velocity)
% title('Growth Velocity vs. Time')
% xlabel('Growth Velocity (m/s)')
% ylabel('Bulk Salinity (ppt)')

% k_eff=(1/S_bottom)*salt_hold(2:end);
% log_k_eff=log((1./k_eff)-1);
% figure
% plot(100*velocity,log_k_eff)
% title('Growth Velocity vs. Log(Keff)')
% xlabel('Growth Velocity (cm/s)')
% ylabel('ln((1/Keff)-1)')