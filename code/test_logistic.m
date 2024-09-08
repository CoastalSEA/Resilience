  x = (0:0.1:4);
  params = struct('A',0,'B',3,'C',1,'K',1,'nu',2,'M',0.5,'Q',2);
  n = [.2,.5,1,2,5];
  testparam = 'nu';
  ispos = false;
  figure('Tag','PlotFig'); 
  axes;
  hold on
  if ispos
    plot([0,0.8,1,2],[0,0.5,0.8,1],'DisplayName','Linear segments')
  else
    plot([0,0.8,1,4],[1,0.8,0.2,0],'DisplayName','Linear segments')
  end
  for i=1:5
     params.(testparam) = n(i);
      y = general_logistic(x,params,~ispos);      
      plot(x,y,'DisplayName',num2str(n(i)))
  end
  legend
  title(testparam)
  hold off