function [J,Cct,PLV,dwPLI] = calcConn(pn,onsetLabel)
%calculate the coherence between the significantly activated chns of one
%patient
%
load(env.get(pn,onsetLabel,'eeg_emg_trial'),'eeg_trial_task','S');
eeg = [eeg_trial_task];%ts x chn x trial
clear eeg_trial_task
load(env.get(pn,onsetLabel,'rtr2_time'),'fdr_dif_r2');
I = sum(fdr_dif_r2<0.001,1)>0;
eeg = eeg(:,I,:);
[tsLen,chnNo,trialNo] = size(eeg);
%parameter
win = 0.2;
stepwin = 0.1;
fs = 1000;
winLen = win*fs;
nfft = 2.^nextpow2(winLen);
stepLen = stepwin * fs;
Nwin = fix((tsLen-winLen)/stepLen)+1;
t = stepwin:stepwin:stepwin*Nwin;
f = linspace(0,fs,nfft);
findx = (f>0) & (f<150);
f =f(findx);
%multitaper
taper = [3 5];
multitapers = dpss(winLen,taper(1),taper(2));%ts x taper

%fft
J = zeros(Nwin,length(f),taper(2),trialNo,chnNo);
Cct = nan(Nwin,chnNo,chnNo,length(f));PLV = Cct;dwPLI = Cct;
for i=1:Nwin
    idx = ((i-1)*stepLen+1):((i-1)*stepLen+winLen);
    eegtmp = permute(...
        repmat(eeg(idx,:,:),1,1,1,taper(2)),...
        [1 4 3 2]);%ts x taper x trial x chn
    eegtmp = bsxfun(@times,eegtmp,multitapers);
    Jtmp = fft(eegtmp,nfft);%nfft x taper x trial x chn
    Jtmp = Jtmp(findx,:,:,:);
    J(i,:,:,:,:) = Jtmp;
    [Cct(i,:,:,:),PLV(i,:,:,:),dwPLI(i,:,:,:)]...
        = coherence.simpleCoh('ct',Jtmp);%chn x chn x f
end
[mC,mPLV,mPLI] = samfnmultvar(@(x) C2mC(x,f,Nwin,t,pn),Cct,PLV,dwPLI);
S = matdoc('win',win,'stepwin',stepwin,'Nwin',Nwin,'f',f,'t',t);
save(fullfile(env.get(pn),[onsetLabel,'_','coh.mat']),'J','Cct','PLV','dwPLI','S','-v7.3');
end

