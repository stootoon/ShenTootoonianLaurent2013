function MakeFiguresForPaper(whichPanels)
% MakeFiguresForPaper(whichPanels)
%
% whichPanels = 1 => Figure 7A: Example KC rasters .

kcSpt = ConvertSpikeTimesFromSparseToFull(LoadTocSpikeTimes('rawkc'));
whichOdors = {'W','WYZ', 'ALL'};
odorIds = cellfun(@(x) odor_name_to_index(['odor' x]), whichOdors);
t0 = 2;
t1 = 4;

if (any(whichPanels==1))
  whichKcs = [207 183 75];

  sfigure(FindFigureCreate('Figure 7A')); clf; set(gcf,'NumberTitle','off','Color',[1 1 1]);
  ResizeFigure(gcf,12,2,'inches');
  Q = ComputeSubplotPositions(3,3,[],0.0,0.0,0.05,0.175,0.01,0.02);
  
  spt = {};
  ip = 1;
  spikeWidth = 0.015;
  for i = 1:numel(whichKcs)
    for j = 1:numel(whichOdors)
      spt = GetSubsetOfTocMatrix(kcSpt, [7 44 209], {1:7, odorIds(j),  whichKcs(i)});
      subplotp(Q, ip);
      ispt = find(spt>=t0 & spt<=t1);
      tt = spt(ispt);
      [ord, trial] = find(spt>=t0 & spt<=t1);
      
      Xspt = bsxfun(@plus, tt(:), [0 spikeWidth spikeWidth 0])';
      Yspt = bsxfun(@plus, trial(:)-1, [0 0 1 1])';
      patch(Xspt, Yspt,'k');
      xlim([t0 t1]);
      ylim([0 7]);
      axis ij;
      set(gca,'xtick',[],'ytick',[]);
      set(line([t0 t1],[7 7]),'LineWidth',0.25,'Color',[0 0 0]);
      axis off;
      ip = ip + 1;
      if (i==1)
        title(['Odor ' whichOdors{j}],'FontSize',18);
      end
    end
  end
end
