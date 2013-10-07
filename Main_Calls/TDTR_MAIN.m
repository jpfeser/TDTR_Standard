%% These first 11 lines of code are not important
keepdir=0;
if keepdir==1
    save('bla.mat','filename','pathname','keepdir')
else
    save('bla.mat','keepdir')
end
clear all
pause(0.1)
load('bla.mat','filename','pathname','keepdir');
tic %Start the timer

%% -------------ENTER  THERMAL SYSTEM PARAMTERS HERE--------------
%definitions of THERMAL SYSTEM PARAMTERS:
%lambda: through-plane thermal conductivity (w/m-K)
%C: heat capacity (J/m3-K)
%h: layer thickness (m)
%eta: anistropic ratio, lambda_{radial} / lambda_{through-plane}
%f:  pump modulation frequency
%r_pump: pump 1/e^2 radius, m
%r_probe: pump 1/e^2 radius, m
%A_pump: pump laser power (Watts)
%TCR:  coefficient of thermal reflectance
%tau_rep: %laser repetition period, s (1/repetition rate)

abslayer =10;
lambda=[170*abslayer 170 0.1 140]; %W/m-K
C=[2.42*abslayer 2.42 0.1 1.7]*1e6; %J/m^3-K
h=[1 (80-abslayer) 1 1e6]*1e-9; %m ("Al" thickness estimated at 84nm=81nm Al+3nm oxide, from picosecond acoustic)
eta=ones(1,numel(lambda)); %isotropic layers, eta=kx/ky;

f=9.8e6; %laser Modulation frequency, Hz
r = 10e-6; %spot sizes, 1/e^2 radius, m
%usually assume they are the same size, otherwise alter below
r_pump=r; %pump 1/e^2 radius, m
r_probe=r; %probe 1/e^2 radius, m

A_pump=30e-3; %laser power (Watts) . . . only used for amplitude est.
TCR=1e-4; %coefficient of thermal reflectance . . . only used for amplitude est.
tau_rep=1/80e6; %laser repetition period, s
%------------- CALCULATE STEADY STATE HEATING--------------
absorbance = 0.1; %fraction of incident light absorbed at surface
A_tot_powermeter = 20e-3; %incident light intensity, Watts
dT_SS = SS_Heating(lambda,C,h,eta,r_pump,r_probe,absorbance,A_tot_powermeter) %avg. steady state temperature rise

%------------- SELECT TIME DELAYS TO SIMULATE OR FIT OVER--------------

%choose time delays for Sensitivity plots
tdelay_sense_min = 100e-12;%earliest time delay for sensitivity plots
tdelay_sense_max = 4000e-12;%latest time delay for sensitivity plots
tdelay = logspace(log10(tdelay_sense_min),log10(tdelay_sense_max),20)'; %vector of time delays (used to generate sensitivity plots)

%Choose range of time delays to fit data, (sec)
tdelay_min=100e-12; %earliest time delay for data fitting
tdelay_max=4000e-12; %latest time delay for data fitting

%----------------------------------PROGRAM OPTIONS BEGIN--------
%Generate Sensitivity Plots?
senseplot=1;
%Import Data? 0 for no, 1 for yes
importdata=1;
%If so, which variable(s) are you fitting for?
X_Fit_Variables=[lambda(4) lambda(3)];%, lambda(3)]; %initial guess for solution, could be for simulated data (if nothing imported) or real data
%Calculate Errorbars? 0 for no, 1 for yes (takes longer)
ebar=1;

%% ERRORBAR SETTINGS
if ebar==1
    %Initialize values
    C_consider=ones(length(C),1);
    L_consider=ones(length(C),1);
    h_consider=ones(length(C),1);
    r_probe_consider=1;
    r_pump_consider=1;
    phase_consider=1;
    CErr=zeros(length(C),length(X_Fit_Variables));
    lambdaErr=zeros(length(lambda),length(X_Fit_Variables));
    h_Err=zeros(length(lambda),length(X_Fit_Variables));
    r_probeErr=zeros(1,length(X_Fit_Variables));
    r_pumpErr=zeros(1,length(X_Fit_Variables));
    
    %define percent uncertainty in each layer/parameter
    Cperc_err=0.05*ones(size(C)); %percent uncertainty in specific heat %i.e. 0.02 -> 2% uncertainty
    lambdaperc_err=0.1*ones(size(lambda));% percent uncertainty in thermal conductivy
    h_perc_err=0.05*ones(size(h));  % percent uncertainty in layer thickness
    r_err=0.1;  % percent uncertainty in beam radius
    degphase=0.2;  %phase error in degree
    
    %VARIABLES CONSIDERED TO CONTRIBUTE ZERO ERRORBAR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Substrate (semi-infinite thickness)
    h_consider(length(h))=0; 
    %
    %Interfaces (no C or h)
    C_consider(3)=0; %thermal interface layer has no capacitance
    h_consider(3)=0; %thermal interface conductance at fixed h/vary lambda
    %
    %Fitting Variables (varying the initial guess, shouldn't change the
    %fit)
    L_consider(4)=0; %solving for this
    L_consider(3)=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%------------------SHOULDN'T NEED TO ALTER BELOW HERE---------------------

[deltaR_data,ratio_data]=TDTR_REFL(tdelay,TCR,tau_rep,f,lambda,C,h,eta,r_pump,r_probe,A_pump);

%% Make Sensitivity Plots
if senseplot==1
    MakeSensitivityPlots
end

%% IMPORT DATA
%--------------Import Data---------------------------
if importdata==1
    if keepdir==1
        filename=input('enter file name string:\n')
    else
        [filename,pathname]=uigetfile('*.*','Choose data file');
    end
    DATAMATRIX=dlmread(strcat(pathname,filename),'\h');
    tdelay_raw=DATAMATRIX(:,2)*1e-12; %imported in picoseconds.  Change to seconds.
    Vin_raw=DATAMATRIX(:,3); %Units (uV ?)
    Vout_raw=DATAMATRIX(:,4);
    ratio_raw=DATAMATRIX(:,5);
    [tdelay_waste,Vin_data] = extract_interior(tdelay_raw,Vin_raw,tdelay_min,tdelay_max);
    [tdelay_waste,Vout_data] = extract_interior(tdelay_raw,Vout_raw,tdelay_min,tdelay_max);
    [tdelay_data,ratio_data] = extract_interior(tdelay_raw,ratio_raw,tdelay_min,tdelay_max);
    
    %--------------Perform Fit (skips if no data import)--------------------------
    Xsol=fminsearch(@(X) TDTR_FIT(X,ratio_data,tdelay_data,TCR,tau_rep,f,lambda,C,h,eta,r_pump,r_probe,A_pump),X_Fit_Variables)
    X_Fit_Variables=Xsol;
    tdelay=tdelay_data;
    fprintf('Data fit completed\n')
else
    [tdelay_waste,Vin_data] = extract_interior(tdelay,real(deltaR_data),tdelay_min,tdelay_max);
    [tdelay_waste,Vout_data] = extract_interior(tdelay,imag(deltaR_data),tdelay_min,tdelay_max);
    [tdelay,ratio_data] = extract_interior(tdelay,ratio_data,tdelay_min,tdelay_max);
    Xsol = X_Fit_Variables;
end

%% COMPUTE ERRORBARS
%--------------Compute Errorbars---------------------
if ebar==1
    fprintf('Calculating Errobar\n')
    fprintf('YErr(n,:) = uncertainty (absolute) in X due to uncertainty in parameter Y(n)\n')
    Xsol=fminsearch(@(X) TDTR_FIT(X,ratio_data,tdelay,TCR,tau_rep,f,lambda,C,h,eta,r_pump,r_probe,A_pump),X_Fit_Variables)
    %Initial Guess, Variable to fit experiment to ... must change "TDTR_FIT"program also
    for ii=1:length(lambda)
        %-------Specific Heat--------------
        if C_consider(ii)==1
            Ctemp=C;
            Ctemp(ii)=C(ii)*(1+Cperc_err(ii));
            Xsoltemp=fminsearch(@(X) TDTR_FIT(X,ratio_data,tdelay,TCR,tau_rep,f,lambda,Ctemp,h,eta,r_pump,r_probe,A_pump),X_Fit_Variables);
            for n=1:length(X_Fit_Variables)
                CErr(ii,n)=abs(Xsoltemp(n)-X_Fit_Variables(n)); %Error in X(n) due to variable C(ii)
            end
            CErr
        end
        
        %-------Thermal Conductivity--------------
        if L_consider(ii)==1
            lambdatemp=lambda;
            lambdatemp(ii)=lambda(ii)*(1+lambdaperc_err(ii));;
            Xsoltemp=fminsearch(@(X) TDTR_FIT(X,ratio_data,tdelay,TCR,tau_rep,f,lambdatemp,C,h,eta,r_pump,r_probe,A_pump),X_Fit_Variables);
            for n=1:length(X_Fit_Variables)
                lambdaErr(ii,n)=abs(Xsoltemp(n)-X_Fit_Variables(n)); %Error in X(n) due to variable lambda(ii)
            end
            lambdaErr
        end
        
        %-------Layer Thickness--------------
        if h_consider(ii)==1
            h_temp=h;
            h_temp(ii)=h(ii)*(1+h_perc_err(ii));
            Xsoltemp=fminsearch(@(X) TDTR_FIT(X,ratio_data,tdelay,TCR,tau_rep,f,lambda,C,h_temp,eta,r_pump,r_probe,A_pump),X_Fit_Variables);
            for n=1:length(X_Fit_Variables)
                h_Err(ii,n)=abs(Xsoltemp(n)-X_Fit_Variables(n)); %Error in X(n) due to variable h(ii)
            end
            h_Err
        end
    end
    %-------Probe Radius-------------- %UPDATE Consider as simulatenous error
    %with pump since they are not independent.
    if r_probe_consider==1
        r_probetemp=r_probe*(1+r_err);
        r_pumptemp=r_pump*(1+r_err);
        Xsoltemp=fminsearch(@(X) TDTR_FIT(X,ratio_data,tdelay,TCR,tau_rep,f,lambda,C,h,eta,r_pumptemp,r_probetemp,A_pump),X_Fit_Variables);
        for n=1:length(X_Fit_Variables)
            r_probeErr(1,n)=abs(Xsoltemp(n)-X_Fit_Variables(n)); %Error in X(n) due to variable r_probe(ii)
        end
        r_probeErr
    end
    %-------Pump Radius-------------- %UPDATE Consider as same.
    % if r_pump_consider==1
    % r_pumptemp=r_pump*(1+r_err);
    % Xsoltemp=fminsearch(@(X) TDTR_FIT(X,ratio_data,tdelay,TCR,tau_rep,f,lambda,C,h,eta,r_pumptemp,r_probe,A_pump),X_Fit_Variables);
    % for n=1:length(X_Fit_Variables)
    %     r_pumpErr(1,n)=abs(Xsoltemp(n)-X_Fit_Variables(n)); %Error in X(n) due to variable r_pump(ii)
    % end
    % r_pumpErr
    % end
    %--------Phase Error-------------
    if phase_consider==1
        radphase=pi/180*degphase;
        Vtemp=(Vin_data+sqrt(-1)*Vout_data)*exp(sqrt(-1)*radphase);
        Vin_phaseshifted=real(Vtemp);
        Vout_phaseshifted=imag(Vtemp);
        ratio_phaseshifted=-Vin_phaseshifted./Vout_phaseshifted;
        Xsoltemp=fminsearch(@(X) TDTR_FIT(X,ratio_phaseshifted,tdelay,TCR,tau_rep,f,lambda,C,h,eta,r_pump,r_probe,A_pump),X_Fit_Variables);
        for n=1:length(X_Fit_Variables)
            phaseErr(1,n)=abs(Xsoltemp(n)-X_Fit_Variables(n));
        end
        phaseErr
    end
    
    ErrSummary_perc=[CErr;lambdaErr;h_Err;r_probeErr;phaseErr]./(ones(3*length(lambda)+2,1)*Xsol); %percent error broken by variable
    
    fprintf('Errorbar calculated\n')
    fprintf('Errorbar breakdown:\n')
    fprintf('Percent Err from C:\n')
    CErr
    fprintf('Abs Err from lambda:\n')
    lambdaErr
    fprintf('Abs Err from h:\n')
    h_Err
    fprintf('Abs Err from spot size:\n')
    r_probeErr
    fprintf('Abs Err from phase:\n')
    phaseErr
    %----------------------------------------------------
    fprintf('Total Percent Error:\n')
    kErr_perc=sqrt(sum(ErrSummary_perc.^2,1)) %total percent error in each fitted parameter
    fprintf('Total Absolute Error:\n')
    kErr_abs=kErr_perc.*Xsol %total absolute error in each fitted parameter
end
fprintf('Fitting Solution:\n')
Xsol
toc

%% EXPORT RESULTS
if (importdata|ebar)
    [Z,ratio_model]=TDTR_FIT(Xsol,ratio_data,tdelay,TCR,tau_rep,f,lambda,C,h,eta,r_pump,r_probe,A_pump);
    dlmwrite('last_fit.txt',[tdelay*1e12,ratio_model],'delimiter','\t')
    
    fignum=203;
    figure(fignum)
    semilogx(tdelay,ratio_data,'ko',tdelay,ratio_model,'k-')
    xlabel('time delay (s)','FontSize',18)
    ylabel('-Vin/Vout','FontSize',18)
    title('Data Fit')
    legend('experiment','model')
    set(gca,'FontSize',18)
    axis([min(tdelay) max(tdelay) 0 1.2*max([ratio_data;ratio_model])])
end

%----------------------------------------------------
saveint=input('Want to save results?\n(0=no, 1=yes)\n');
if saveint==1
    if (importdata&ebar)
    save(strcat(pathname,filename(1:end-4),'Results.mat'))
    print(fignum,'-depsc',strcat(pathname,filename(1:end-4),'FIT.eps'))
    else
        save(strcat('Results',datestr(now,'yyyymmddTHHMMSS'),'.mat'))
    end
end

%----------------------------------------------------
fprintf('Program Completed\n')
beep
beep
beep







