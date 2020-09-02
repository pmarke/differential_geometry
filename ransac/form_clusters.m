function clusters = form_clusters(measurments,t,vel_max,t_max)
siz = size(measurments);
meas = [];

clusters = {};
clus = [];

for ii = 1:length(t)
   for jj = 1:siz(2)
       S.pt = measurments(:,jj,ii);
       S.t = t(ii);
       meas = [meas,S];
   end
end

done = false;
while (~done)
    
  S1 = meas(1);
  meas(1) = [];
  clus = S1;
  index = [];
  clus_index = 1;
  done2 = false;
  
  while (~done2)
  
      S1 = clus(clus_index);
      for ii = 1:length(meas)

          if ( norm(S1.pt - meas(ii).pt) < vel_max*abs(S1.t-meas(ii).t) && abs(S1.t-meas(ii).t) < t_max)
              clus = [clus, meas(ii)];
              index = [index, ii];
          end

      end
      meas(index) = [];
      index = [];
      clus_index = clus_index +1;
      
      if (clus_index > length(clus))
          done2 = true;
      end
  end
  
  if (length(clus) > 5)
      clusters = {clusters, clus};
  end
  
  if(length(meas) < 5)
      done = true;
  end
    
end

end