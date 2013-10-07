%The main program tries to minimize "Z" by optimizing the variable(s) X
%This program Lets you:
%   1) Define the vector X: example, if lambda(3) is what you want to solve for,
%   then set X=lambda(3)....if you whish to simulatenous solve for more than one
%   variable, you can just define multiple variables (eg. X(1)=lambda(3), X(2)=lambda(4))
%
%   2) Define the fit.  Typically, this is the sum of the squares of the
%   residuals, but you might want to weight these by the sensitivity,
%   particularly if you don't intend to calculate the errorbars!

function [Z,ratio_model]=TDTR_FIT(X,ratio_data,tdelay,TCR,tau_rep,f,lambda,C,h,eta,r_pump,r_probe,A_pump)

%Define the variables to be fit
%lambda(1)=X(1);
%lambda(5) = X(1);
lambda(4) = X(1);
lambda(3) = X(2);
%lambda(3) = X(2);
%lambda(3)=X(3);
%eta(2)=X(2);
%eta(1)=eta(2);

% don't touch this -------------------
[deltaR_model,ratio_model]=TDTR_REFL(tdelay,TCR,tau_rep,f,lambda,C,h,eta,r_pump,r_probe,A_pump);
%Uncomment the next three lines to see the non-linear optimization in
%action!
figure(10)
semilogx(tdelay,ratio_data,'ob',tdelay,ratio_model,'g'); 
axis([100e-12 10e-9 0 max([ratio_data;ratio_model])])
pause(0.1)
X
res=(ratio_model-ratio_data).^2;
Z=sum(res);
%--------------------------------------