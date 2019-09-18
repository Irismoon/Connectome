function visualConn(pn)
%load coherence calculated from calcConn
onsetLabel = 'emg_onset';
load(fullfile(env.get(pn),[onsetLabel,'_','coh.mat']),'Cct','PLV','dwPLI','S');%t x chn x chn x f
%
C = permute(Cct(:,:,:,3),[2 3 1]);%chn x chn x t
win = [10:2:36];
t = -2.9:0.1:6.9;
baseC = mean(abs(C(:,:,1:5)),3);
C = abs(C(:,:,win));
[chnNo,~,Nwin] = size(C);
% No = 10;
% indexPair = zeros(No,2,Nwin);% the index of one pair
% cohPair = zeros(No,Nwin);%the coherence strength of one pair
% for i = 1:Nwin
%     Ctmp = C(:,:,i);
%     Ctmp = abs(Ctmp -baseC);
%     [Ctmp,I] = sort(reshape(Ctmp,[],1),'descend');
%     idx = ~isnan(Ctmp);
%     [Ctmp,I] = samfnmultvar(@(x) x(idx),Ctmp,I);
%     %     thld = quantile(Ctmp,0.5);
%     %     idx = find(Ctmp>thld);
%     %     idx = idx(end-No:end);
%     %     [Ctmp,I] = samfnmultvar(@(x) x(idx),Ctmp,I);
%     [indexPair(:,1,i),indexPair(:,2,i)] = ind2sub([chnNo,chnNo],I(1:No));%row and col is a pair
%     cohPair(:,i) = Ctmp(1:No);
% end
% cohPair = reshape(zscore(cohPair(:)),No,Nwin);
% cohPair = cohPair+0.1+abs(min(cohPair(:)));
% figure,
% %parameter
% %load brain model
% load(fullfile(env.get(pn,'wholecortex')),'cortex');
% inflColor = repmat([0.5 0.5 0.5],length(cortex.vert),1);%vertex x 3
% alpha = 0.3;
% hemi = whichhemi(pn);
% [cmapstruct,viewstruct] = getdispstruct([],0,hemi,{'brain'});
% %load elec position
% pos = getElecInfo(pn,'pos');
% pos = pos(actchn(pn,onsetLabel),:);
% %at each time, plot its 10 strongest coherent pairs
% chnpos = zeros(No,2,3);
% median(expH.getEmgDelay(pn))
% for i=1:Nwin
%     clf
%     ha = trisurf(cortex.tri, cortex.vert(:, 1), cortex.vert(:, 2), cortex.vert(:, 3), 'FaceVertexCData', inflColor,...
%         'FaceColor', 'interp', 'CDataMapping', 'direct', 'linestyle', 'none','FaceAlpha',alpha);
%     generalview;
%     %
%     chnpos(:,1,:) = pos(indexPair(:,1,i),:);%positive, 10 x 3
%     chnpos(:,2,:) = pos(indexPair(:,2,i),:);%negative
%     uchn = unique(indexPair(:,:,i));
%     %     uchnpos = pos(uchn,:);
%     color = brewermap(2,'Set1');
%     color = repmat(color(1,:),length(pos),1);
%     color(uchn,:) = repmat(color(2,:),length(uchn),1);
%     plotBalls(pos,color,[]);hold on;
%     for j=1:No
%         line(chnpos(j,:,1),chnpos(j,:,2),chnpos(j,:,3),'Color','green','LineW',cohPair(j,i));
%     end
%     title(t(win(i)));
%     pause;
% end
%
dC = abs(C - baseC);%chn x chn x win
dC = trimSingleShaft(pn,dC,actchn(pn,onsetLabel));
dC = trimMultShaft(pn,dC,actchn(pn,onsetLabel));
dCtmp = dC(:);
idxtmp = isnan(dCtmp);
dCtmp(idxtmp) = [];
dCtmp = zscore(dCtmp);
iqr = quantile(dCtmp,[0.25 0.75]);
% mask = ~triu(true(chnNo,chnNo));
dC(~idxtmp) = dCtmp;
dC(dC<(iqr(2)+5*diff(iqr)))= nan;
cohPair = dC(~isnan(dC));
[indexPair(:,1),indexPair(:,2),indexPair(:,3)] = ind2sub([chnNo chnNo Nwin],find(~isnan(dC)));
assert(issorted(indexPair(:,3)),'error!');

load(fullfile(env.get(pn,'wholecortex')),'cortex');
inflColor = repmat([0.5 0.5 0.5],length(cortex.vert),1);%vertex x 3
alpha = 0.3;
hemi = whichhemi(pn);
[cmapstruct,viewstruct] = getdispstruct([],0,hemi,{'brain'});
%load elec position
pos = getElecInfo(pn,'pos');
pos = pos(actchn(pn,onsetLabel),:);
%at each time, plot its 10 strongest coherent pairs
median(expH.getEmgDelay(pn))
figure,
ha = trisurf(cortex.tri, cortex.vert(:, 1), cortex.vert(:, 2), cortex.vert(:, 3), 'FaceVertexCData', inflColor,...
    'FaceColor', 'interp', 'CDataMapping', 'direct', 'linestyle', 'none','FaceAlpha',alpha);
generalview;
color = repmat(rgb('red'),chnNo,1);
plotBalls(pos,color,[]);hold on;
for i=1:Nwin
    delete(findobj(gca,'Type','line'));
    idx = find(indexPair(:,3) == i );
    if ~isempty(idx)
        pair = indexPair(idx,1:2);
        chnpos = zeros(length(idx),2,3);
        chnpos(:,1,:) = pos(pair(:,1),:);
        chnpos(:,2,:) = pos(pair(:,2),:);
        for j=1:length(idx)
            line(chnpos(j,:,1),chnpos(j,:,2),chnpos(j,:,3),...
                'Color','green','LineW',cohPair(idx(j))*0.3);
        end
        title(t(win(i)));
        pause;
    end
end



