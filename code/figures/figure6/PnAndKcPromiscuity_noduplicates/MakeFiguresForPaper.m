function MakeFiguresForPaper(whichPanels)
% function MakeFiguresForPaper(whichPanels)
%
% whichPanels=1 => Figure S6C: PN response promsicuity
% whichPanels=2 => Figure S6D: KC response promsicuity
% whichPanels=3 => Figure S6A: PN responsivity
% whichPanels=4 => Figure S6B: Figure S6B: KC responsivity

Rk = @(X,d) ComputeResponsivityFromSpikeTimesLikeKai(X,d,1.5,2.1,3.1,0.2,1.5,4);

g = OdorGroupingFunction1(GetOdorsList);
whichOdors = find(g~=-1); % Everything but Paraffin oil.
odorBnds = arrayfun(@(i) find(g==i,1)-0.5, [1:5]);

if (any(whichPanels==1) | any(whichPanels==3))
  pnSpt = LoadTocSpikeTimes('rawpn');
  Rpn = Rk(pnSpt, [7 44 174]);
end

if (any(whichPanels==2) | any(whichPanels==4))
  kcSpt = LoadTocSpikeTimes('rawkc');
  Rkc = Rk(kcSpt, [7 44 209]);
end

if (any(whichPanels==1))
  sfigure(FindFigureCreate('Figure S6C: PN response promsicuity')); clf; 
  set(gcf,'Color',[1 1 1],'NumberTitle','off');
  ResizeFigure(gcf,8,8,'inches');
  hist(sum(Rpn(:,whichOdors),2),[0:43]);
  ylim([0 15]);
  xlim([-1 44]);
  set(gca,'ytick',0:3:15,'xtick',[0 43],'ytick',0:3:15,'ticklength',[0 0],'FontSize',12);
  set(findobj(gca,'Type','patch'), 'FaceColor','r');
  xlabel('# odors responded to', 'FontSize', 14);
  ylabel('frequency observed', 'FontSize', 14);
  box on;
  TightenAxesToFigure;
end

if (any(whichPanels==2))
  sfigure(FindFigureCreate('Figure S6D: KC response promsicuity')); clf; 
  set(gcf,'Color',[1 1 1],'NumberTitle','off');
  ResizeFigure(gcf,8,8,'inches');
  hist(sum(Rkc(:,whichOdors),2),[0:43]);
  ylim([0 100]);
  xlim([-1 44]);
  set(gca,'ytick',0:20:100,'xtick',[0 43],'tickLength',[0 0],'FontSize',12);
  set(findobj(gca,'Type','patch'), 'FaceColor','b');
  xlabel('# odors responded to', 'FontSize', 14);
  ylabel('frequency observed', 'FontSize', 14);
  box on;
  TightenAxesToFigure;
end

if (any(whichPanels==3))
  sfigure(FindFigureCreate('Figure S6A: PN responsivity')); clf; 
  set(gcf,'Color',[1 1 1],'NumberTitle','off');
  ResizeFigure(gcf, 12, 3,'inches');
  imagesc(Rpn(:, whichOdors)');
  colormap(1-gray(2));
  ytick = [0 odorBnds 44];
  ytick = (ytick(1:end-1)+ytick(2:end))/2;
  set(gca,'xtick',[1 174],'ytick',ytick,'yticklabel',{'4x 1-','1-','2-','3-','4-','5-, 8-'});
  arrayfun(@(u) set(line([0 175],[u u]),'Color','r'), odorBnds);  
  set(gca,'FontSize',12,'ticklength',[0 0]);
  xlabel('PN index','FontSize', 14);
  ylabel('mixture level', 'FontSize', 14);
  axis equal;
  axis tight;
end

if (any(whichPanels==4))
  sfigure(FindFigureCreate('Figure S6B: KC responsivity')); clf; 
  set(gcf,'Color',[1 1 1],'NumberTitle','off');
  ResizeFigure(gcf,12,3,'inches');
  imagesc(Rkc(:, whichOdors)');
  colormap(1-gray(2));
  ytick = [0 odorBnds 44];
  ytick = (ytick(1:end-1)+ytick(2:end))/2;
  set(gca,'xtick',[1 209],'ytick',ytick,'yticklabel',{'4x 1-','1-','2-','3-','4-','5-, 8-'});
  arrayfun(@(u) set(line([0 209],[u u]),'Color','r'), odorBnds);  
  set(gca,'FontSize',12,'ticklength',[0 0]);
  xlabel('KC index','FontSize', 14);
  ylabel('mixture level', 'FontSize', 14);
  axis equal;
  axis tight;
end
