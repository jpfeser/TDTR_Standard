%These first 10 lines of code are not important
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

%-------------TYPE THERMAL SYSTEM PARAMTERS HERE--------------
%Anticipated system properties (initial guess for fitting, if you do fitting/errorbar estimation)
abslayer =10;
lambda=[170*abslayer 170 0.7 32 55]; %W/m-K
C=[2.42*abslayer 2.42 1.4 1.4 1.7]*1e6; %J/m^3-K
t=[1 (80-abslayer) 20 1e3 1e6]*1e-9; %m ("Al" thickness estimated at 84nm=81nm Al+3nm oxide, from picosecond acoustic)
eta=ones(1,numel(lambda)); %isotropic layers, eta=kx/ky;
eta(3)=1;

f=9.8e6; %laser Modulation frequency, Hz
r = 10e-6;

absorbance = 0.1;
A_tot_powermeter = 20e-3;

tau_rep=1/80e6; %laser repetition period, s
r_pump=r; %pump 1/e^2 radius, m
r_probe=r; %probe 1/e^2 radius, m
A_pump=30e-3; %laser power (Watts) . . . only used for amplitude est.
TCR=1e-4; %coefficient of thermal reflectance . . . only used for amplitude est.

dT_SS = SS_Heating(lambda,C,t,eta,r_pump,r_probe,absorbance,A_tot_powermeter)
%------------- THERMAL SYSTEM PARAMTERS END HERE--------------

%choose time delays for Sensitivity plots
tdelay=logspace(log10(100e-12),log10(4e-9),20)'; %vector of time delays (used to generate sensitivity plots)

%Choose range of time delays to fit, sec
tdelay_min=100e-12;
tdelay_max=4000e-12;

%----------------------------------PROGRAM OPTIONS BEGIN--------
%Generate Sensitivity Plots?
senseplot=1;
%Import Data? 0 for no, 1 for yes
importdata=1;
%If so, which variable(s) are you fitting for?
Xguess=[lambda(4) lambda(3)];%, lambda(3)]; %initial guess for solution, could be for simulated data (if nothing imported) or real data
%Calculate Errorbars? 0 for no, 1 for yes (takes longer)
ebar=0;
%----------------------------------PROGRAM OPTIONS END--------

%---------------------ERRORBAR OPTIONS----------------------
if ebar==1
    %Initialize values
    C_consider=ones(length(C),1);
    L_consider=ones(length(C),1);
    t_consider=ones(length(C),1);
    r_probe_consider=1;
    r_pump_consider=1;
    phase_consider=1;
    CErr=zeros(length(C),length(Xguess));
    lambdaErr=zeros(length(lambda),length(Xguess));
    tErr=zeros(length(lambda),length(Xguess));
    r_probeErr=zeros(1,length(Xguess));
    r_pumpErr=zeros(1,length(Xguess));
    
    %define parameters NOT to consider in error analysis (saves time)
    t_consider(length(t))=0; %last layer is semi-inf
    C_consider(3)=0; %thermal interface layer has no capacitance
    t_consider(3)=0; %thermal interface conductance at fixed t/vary lambda
    L_consider(3)=0; %solving for this 
    %L_consider(4)=0;
    

    %define percent uncertainty in each layer/parameter
    Cperc_err=0.05*ones(size(C)); %percent uncertainty in specific heat %i.e. 0.02 -> 2% uncertainty
    lambdaperc_err=0.1*ones(size(lambda));% percent uncertainty in thermal conductivy
    tperc_err=0.05*ones(size(t));  % percent uncertainty in layer thickness
    r_err=0.1;  % percent uncertainty in beam radius
    degphase=0.2;  %phase error in degree
    
end

%-----------------Make sensitivity plots--------------
[deltaR_data,ratio_data]=TDTR_REFL_VV3(tdelay,TCR,tau_rep,f,lambda,C,t,eta,r_pump,r_probe,A_pump);
if senseplot==1
    MakeSensitivityPlots
end

%--------------Import Data---------------------------
if importdata==1
    if keepdir==1
        filename=input('enter file name string:\n')
    else
        [filename,pathname]=uigetfile('*.*','Choose data file');
    end
    DATAMATRIX=dlmread(strcat(pathname,filename),'\t');
    tdelay_raw=DATAMATRIX(:,2)*1e-12; %imported in picoseconds.  Change to seconds.
    Vin_raw=DATAMATRIX(:,3); %Units (uV ?)
    Vout_raw=DATAMATRIX(:,4);
    ratio_raw=DATAMATRIX(:,5);
    [tdelay_data,Vin_data] = extract_interior(tdelay_raw,Vin_raw,tdelay_min,tdelay_max);
    [tdelay_data,Vout_data] = extract_interior(tdelay_raw,Vout_raw,tdelay_min,tdelay_max);
    [tdelay_data,ratio_data] = extract_interior(tdelay_raw,ratio_raw,tdelay_min,tdelay_max);
    
%--------------Perform Fit (skips if no data import)--------------------------
    Xsol=fminsearch(@(X) TDTR_FIT_VV3(X,ratio_data,tdelay_data,TCR,tau_rep,f,lambda,C,t,eta,r_pump,r_probe,A_pump),Xguess)
    Xguess=Xsol;
    tdelay=tdelay_data;
    fprintf('Data fit completed\n')
else
    [tdelay_waste,Vin_data] = extract_interior(tdelay,real(deltaR_data),tdelay_min,tdelay_max);
    [tdelay_waste,Vout_data] = extract_interior(tdelay,imag(deltaR_data),tdelay_min,tdelay_max);
    [tdelay,ratio_data] = extract_interior(tdelay,ratio_data,tdelay_min,tdelay_max);
    Xsol = Xguess;
end

%--------------Compute Errorbars---------------------
if ebar==1
fprintf('Calculating Errobar\n')    
fprintf('YErr(n,:) = uncertainty (absolute) in X due to uncertainty in parameter Y(n)\n')
Xsol=fminsearch(@(X) TDTR_FIT_VV3(X,ratio_data,tdelay,TCR,tau_rep,f,lambda,C,t,eta,r_pump,r_probe,A_pump),Xguess)
%Initial Guess, Variable to fit experiment to ... must change "TDTR_FIT_VV3"program also
    for ii=1:length(lambda)
        %-------Specific Heat--------------
        if C_consider(ii)==1
            Ctemp=C;
            Ctemp(ii)=C(ii)*(1+Cperc_err(ii));
            Xsoltemp=fminsearch(@(X) TDTR_FIT_VV3(X,ratio_data,tdelay,TCR,tau_rep,f,lambda,Ctemp,t,eta,r_pump,r_probe,A_pump),Xguess);
            for n=1:length(Xguess)
                CErr(ii,n)=abs(Xsoltemp(n)-Xguess(n)); %Error in X(n) due to variable C(ii)
            end
            CErr
        end
        
        %-------Thermal Conductivity--------------
        if L_consider(ii)==1
            lambdatemp=lambda;
            lambdatemp(ii)=lambda(ii)*(1+lambdaperc_err(ii));;
            Xsoltemp=fminsearch(@(X) TDTR_FIT_VV3(X,ratio_data,tdelay,TCR,tau_rep,f,lambdatemp,C,t,eta,r_pump,r_probe,A_pump),Xguess);
            for n=1:length(Xguess)
                lambdaErr(ii,n)=abs(Xsoltemp(n)-Xguess(n)); %Error in X(n) due to variable lambda(ii)
            end
            lambdaErr
        end
        
        %-------Layer Thickness--------------
        if t_consider(ii)==1
            ttemp=t;
            ttemp(ii)=t(ii)*(1+tperc_err(ii));
            Xsoltemp=fminsearch(@(X) TDTR_FIT_VV3(X,ratio_data,tdelay,TCR,tau_rep,f,lambda,C,ttemp,eta,r_pump,r_probe,A_pump),Xguess);
            for n=1:length(Xguess)
                tErr(ii,n)=abs(Xsoltemp(n)-Xguess(n)); %Error in X(n) due to variable t(ii)
            end
            tErr
        end
    end
%-------Probe Radius-------------- %UPDATE Consider as simulatenous error
%with pump since they are not independent.
if r_probe_consider==1
r_probetemp=r_probe*(1+r_err);
r_pumptemp=r_pump*(1+r_err);
Xsoltemp=fminsearch(@(X) TDTR_FIT_VV3(X,ratio_data,tdelay,TCR,tau_rep,f,lambda,C,t,eta,r_pumptemp,r_probetemp,A_pump),Xguess);
for n=1:length(Xguess)
    r_probeErr(1,n)=abs(Xsoltemp(n)-Xguess(n)); %Error in X(n) due to variable r_probe(ii)
end
r_probeErr
end
%-------Pump Radius-------------- %UPDATE Consider as same.
% if r_pump_consider==1
% r_pumptemp=r_pump*(1+r_err);
% Xsoltemp=fminsearch(@(X) TDTR_FIT_VV3(X,ratio_data,tdelay,TCR,tau_rep,f,lambda,C,t,eta,r_pumptemp,r_probe,A_pump),Xguess);
% for n=1:length(Xguess)
%     r_pumpErr(1,n)=abs(Xsoltemp(n)-Xguess(n)); %Error in X(n) due to variable r_pump(ii)
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
    Xsoltemp=fminsearch(@(X) TDTR_FIT_VV3(X,ratio_phaseshifted,tdelay,TCR,tau_rep,f,lambda,C,t,eta,r_pump,r_probe,A_pump),Xguess);
for n=1:length(Xguess)
    phaseErr(1,n)=abs(Xsoltemp(n)-Xguess(n));
end
phaseErr
end

ErrSummary_perc=[CErr;lambdaErr;tErr;r_probeErr;phaseErr]./(ones(3*length(lambda)+2,1)*Xsol); %percent error broken by variable

fprintf('Errorbar calculated\n')
fprintf('Errorbar breakdown:\n')
fprintf('Percent Err from C:\n')
CErr
fprintf('Abs Err from lambda:\n')
lambdaErr
fprintf('Abs Err from h:\n')
tErr
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


[Z,ratio_model]=TDTR_FIT_VV3(Xsol,ratio_data,tdelay_data,TCR,tau_rep,f,lambda,C,t,eta,r_pump,r_probe,A_pump);
dlmwrite('last_fit.txt',[tdelay_data*1e12,ratio_model],'delimiter','\t')
%----------------------------------------------------
saveint=input('Want to save results?\n(0=no, 1=yes)\n');
if saveint==1
    save(strcat(pathname,filename(1:end-4),'Results.mat'))
    fignum=203;
    figure(fignum)
    semilogx(tdelay_data,ratio_data,'ko',tdelay_data,ratio_model,'k-')
    xlabel('time delay (s)','FontSize',18)
    ylabel('-Vin/Vout','FontSize',18)
    title('Data Fit')
    legend('experiment','model')
    set(gca,'FontSize',18)
    axis([min(tdelay_data) max(tdelay_data) 0 1.2*max([ratio_data;ratio_model])])
    
    print(fignum,'-depsc',strcat(pathname,filename(1:end-4),'FIT.eps'))
end
%----------------------------------------------------
fprintf('Program Completed\n')
beep
beep
beep






