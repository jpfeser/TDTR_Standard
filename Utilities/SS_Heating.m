function dT_SS = SS_Heating(lambda,C,h,eta,r_pump,r_probe,absorbance,A_tot_powermeter)

if nargin==0
    absorbance =0.5; %Aluminum
    %abslayer =40;
    lambda=[20 0.1 1.35 142]; %W/m-K
    C=[2.4 0.1 1.6 1.6]*1e6; %J/m^3-K
    h=[100 1 500 1e6]*1e-9; %m ("Al" thickness estimated at 84nm=81nm Al+3nm oxide, from picosecond acoustic)
    eta=ones(1,numel(lambda)); %isotropic layers, eta=kx/ky;
    %eta(3)=sqrt(20*0.1)/0.1;
    
    r=10e-6;
    r_pump=r; %pump 1/e^2 radius, m
    r_probe=r; %probe 1/e^2 radius, m
    A_tot_powermeter = 2*(10)*1e-3; %laser power (Watts) . . . (assumes the chopper isn't on!)
    t_rep = 12.5e-9;
    suggested_powers = [A_tot_powermeter/2,A_tot_powermeter/4]
    p_pulse_heating = (A_tot_powermeter/2)*t_rep/(C(1)*h(1)*2*pi*r^2)
end

%-----------DO NOT MODIFY BELOW HERE------------------------------
f=0; %laser Modulation frequency, Hz
A_abs = absorbance*A_tot_powermeter;
kmin=1/(10000*max(r_pump,r_probe));
kmax=1/sqrt(r_pump^2+r_probe^2)*1.5;
dT_SS  =rombint_VV3(@(kvect) TDTR_TEMP_VV3(kvect,f,lambda,C,h,eta,r_pump,r_probe,A_abs),kmin,kmax,1)
%-----------------------------------------------------------------