%Calculates the Reflectivity Signal and Ratio
%In order to speed up the code, it is parallelized...the convention is...
%tdelay is COLUMN vector of desired delay times
%Mvect (the fourier components) are ROW vectors
%Matrices have size, length(tdelay) x length(Mvect)

%PARAMETERS (scalars unless told otherwise):  
%TCR: temperature coefficient of reflectivity  (1/K)...doesn't effect ratio
%tau_rep:  laser repetition time (sec) = 1/f_rep
%f:  modulation frequency 
%lambda:  VECTOR of thermal conductivities (W/mK) (layer 1 is the topmost layer)
%C:  VECTOR of specific heats (J/m3-K)
%t:  VECTOR of layer thicknesses (last layer is alway treated semi-inf, but
%       you still have to enter something).  (m)
%eta:  VECTOR of anisotropic ratio (kx/ky), use ones(length(lambda)) for
%       isotropic
%r_pump:  pump 1/e2 radius (m)
%r_probe: probe 1/e2 radius (m)
%A_pump:  pump intensity (W)...doesn't effect ratio

function [deltaT,ratio]=TDTR_REFL(tdelay,TCR,tau_rep,f,lambda,C,h,eta,r_pump,r_probe,A_pump)

fmax=10/min(abs(tdelay)); %maximum frequency considered (see RSI paper)
ii=sqrt(-1);

M=10*ceil(tau_rep/min(abs(tdelay))); %Highest Fourier component considered
mvect=-M:M; %Range of Fourier components to consider (Vector)
fudge1=exp(-pi*((mvect/tau_rep+f)/fmax).^2);%artificial decay (see RSI paper)
fudge2=exp(-pi*((mvect/tau_rep-f)/fmax).^2);

dT1=zeros(1,length(mvect))';
dT2=zeros(1,length(mvect))';
kmax=1/sqrt(r_pump^2+r_probe^2)*1.5;

    dT1=rombint_multi(@(kvect) TDTR_TEMP(kvect,mvect/tau_rep+f,lambda,C,h,eta,r_pump,r_probe,A_pump),0,kmax,length(mvect));
    dT2=rombint_multi(@(kvect) TDTR_TEMP(kvect,mvect/tau_rep-f,lambda,C,h,eta,r_pump,r_probe,A_pump),0,kmax,length(mvect));
    
expterm=exp(ii*2*pi/tau_rep*(tdelay*mvect));
Retemp=(ones(length(tdelay),1)*(dT1'.*fudge1+dT2'.*fudge2)).*expterm;
Imtemp=-ii*(ones(length(tdelay),1)*(dT1-dT2)').*expterm;

Resum=sum(Retemp,2); %Sum over all Fourier series components
Imsum=sum(Imtemp,2);

%TCR =1;
deltaRm=TCR*(Resum+ii*Imsum); %
deltaR=deltaRm.*exp(ii*2*pi*f*tdelay); %Reflectance Fluxation (Complex)
deltaT = deltaR;

ratio=-real(deltaR)./imag(deltaR);