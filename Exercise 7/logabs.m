function y=logabs(x,gain,dyn);
%LOGABS   Logaritmisk kompresjon av absoluttverdi

    y=abs(x); %Ta absoluttverdi 
    y=y+1e-30; %Legg til eit lite positivt tal for å unngå logaritme av 0
    y=20*log10(y); %Gjer amplitude om til dB
    y=y+gain; %Legg til gain (i dB)
    y=y/dyn;  %Skaler dynamisk område
    y=max(0,y); %Sett negative verdiar lik 0
    y=min(1,y); %sett verdiar >1  lik 1