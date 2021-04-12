clear all;

WE=0:0.01:3;

for s=1:999
    fileName = sprintf('WG_%03d.mat',s);
    load(fileName);
    Iwe=s;
    FCfitt_all(Iwe)=FCfitt;
    FCDfitt_all(Iwe)=FCDfitt;
    GBCfitt_all(Iwe)=GBCfitt;
    ICAfitt_all(Iwe)=ICAfitt;
    KSclu_all(Iwe)=KSclu;
    KLPstates_all(Iwe)=KLPstates;
    Ignition_all(Iwe)=Ignition;
    Ignitionthreshold_all(Iwe)=Ignitionthreshold;
    Peakrate_all(Iwe)=Maxrate;
    Decay_all(Iwe)=Decay;
end

Ignition_all=Ignition_all/max(Ignition_all);
Ignitionthreshold_all=Ignitionthreshold_all/max(Ignitionthreshold_all);
Peakrate_all=Peakrate_all/max(Peakrate_all);

save ratioEI_G.mat WE GBCfitt_all KSclu_all ICAfitt_all KLPstates_all FCDfitt_all FCfitt_all Ignition_all Ignitionthreshold_all Peakrate_all Decay_all;

figure
plot(WE,FCfitt_all,'k--');
hold on;
plot(WE,ICAfitt_all,'k');
plot(WE,GBCfitt_all,'c--');
plot(WE,FCDfitt_all,'r');
plot(WE,KSclu_all,'r--');
plot(WE,KLPstates_all,'c');
plot(WE,Ignition_all,'--b');
plot(WE,Ignitionthreshold_all,'b');
plot(WE,Peakrate_all,'--g');
plot(WE,Decay_all,'g');

figure
plot(WE,Ignition_all.*Ignitionthreshold_all.*Peakrate_all);
