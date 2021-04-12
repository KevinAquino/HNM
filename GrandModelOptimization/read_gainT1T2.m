clear all;

ALPHA=-0.9:0.05:0;
BETA=0:0.05:2;

for s=1:999
    fileName = sprintf('WgainT1T2_%03d.mat',s);
    if exist(fileName, 'file') == 2
        load(fileName);
        [Ialpha Ibeta]=ind2sub([length(ALPHA),length(BETA)],s);
        FCfitt_all(Ialpha,Ibeta)=FCfitt;
        GBCfitt_all(Ialpha,Ibeta)=GBCfitt;
        FCDfitt_all(Ialpha,Ibeta)=FCDfitt;
        Ignition_all(Ialpha,Ibeta)=Ignition;
        Decay_all(Ialpha,Ibeta)=Decay;
        ErrDecay_all(Ialpha,Ibeta)=ErrDecay;
        ErrDecay_all2(Ialpha,Ibeta)=ErrDecay2;
        Excitability_all(Ialpha,Ibeta)=Excitability;
        ErrDecayBold_all(Ialpha,Ibeta)=ErrDecayBold;
    else
        s
        [Ialpha Ibeta]=ind2sub([length(ALPHA),length(BETA)],s);
        FCfitt_all(Ialpha,Ibeta)=0;
        GBCfitt_all(Ialpha,Ibeta)=0;
        FCDfitt_all(Ialpha,Ibeta)=1;
        Ignition_all(Ialpha,Ibeta)=0;
        Decay_all(Ialpha,Ibeta)=0;
        ErrDecay_all(Ialpha,Ibeta)=1;
        ErrDecay_all2(Ialpha,Ibeta)=1;
        Excitability_all(Ialpha,Ibeta)=0;
        ErrDecayBold_all(Ialpha,Ibeta)=0;
    end
end

save ratio_gainT1T2.mat ALPHA BETA Excitability_all GBCfitt_all FCDfitt_all FCfitt_all Ignition_all ErrDecay_all ErrDecay_all2 Decay_all ErrDecayBold_all;

figure;
imagesc(BETA,ALPHA,FCfitt_all);
figure;
imagesc(BETA,ALPHA,GBCfitt_all);
figure;
imagesc(BETA,ALPHA,FCDfitt_all);
figure;
imagesc(BETA,ALPHA,Ignition_all);
figure;
imagesc(BETA,ALPHA,Decay_all);
figure;
imagesc(BETA,ALPHA,ErrDecay_all);
figure;
imagesc(BETA,ALPHA,ErrDecay_all2);
figure;
imagesc(BETA,ALPHA,ErrDecayBold_all);


aa=Ignition_all;
cc=Decay_all;
dd=ErrDecay_all;
dd2=ErrDecay_all2;
bb=FCDfitt_all;
aa1=FCfitt_all;
aa2=GBCfitt_all;
hh=ErrDecayBold_all;
for i=1:length(BETA)
    [FCD1(i) j]=min(bb(:,i));
    jmin=j;
    jmax=j;
    index(i)=j;
    FCD1(i)=mean(bb(jmin:jmax,i));
    Ig1(i)=mean(aa(jmin:jmax,i));
    GBC1(i)=mean(aa2(jmin:jmax,i));
    FC1(i)=mean(aa1(jmin:jmax,i));
    De1(i)=mean(cc(jmin:jmax,i));
    ErrDe1(i)=mean(dd(jmin:jmax,i));
    ErrDe12(i)=mean(dd2(jmin:jmax,i));
    ErrDeBOLD(i)=mean(hh(jmin:jmax,i));
end
% figure
% subplot(5,1,1);
% plot(BETA,Ig1);
% subplot(5,1,2);
% plot(BETA,FCD1);
% subplot(5,1,3);
% plot(BETA,De1);
% subplot(5,1,4);
% plot(BETA,ErrDe1);
% subplot(5,1,5);
% plot(BETA,ErrDe12);
% 
% figure
% plot(BETA,ErrDeBOLD);
% 
% figure
% subplot(2,1,1);
% plot(BETA,FC1);
% subplot(2,1,2);
% plot(BETA,GBC1);

%%%% Paper
figure;
subplot(6,2,1);
imagesc(BETA,ALPHA,FCfitt_all);
subplot(6,2,2);
imagesc(BETA,ALPHA,GBCfitt_all);

subplot(6,2,3);
plot(BETA,FC1);
subplot(6,2,4);
plot(BETA,GBC1);

subplot(6,2,5);
imagesc(BETA,ALPHA,FCDfitt_all);
subplot(6,2,6);
imagesc(BETA,ALPHA,ErrDecayBold_all);

subplot(6,2,7);
plot(BETA,FCD1);
subplot(6,2,8);
plot(BETA,ErrDe12);

subplot(6,2,9);
imagesc(BETA,ALPHA,Ignition_all);
subplot(6,2,10);
imagesc(BETA,ALPHA,Decay_all);

subplot(6,2,11);
plot(BETA,Ig1);
subplot(6,2,12);
plot(BETA,De1);