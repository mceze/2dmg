function plot_leanelem( leanmesh_matlab, elemID)
% plots a matlab output of a xf_LeanMesh

run(leanmesh_matlab)

figure(1)
i = elemID+1
xplt=[coord(elem(i,1),1);coord(elem(i,2),1);coord(elem(i,3),1);coord(elem(i,1),1)];
yplt=[coord(elem(i,1),2);coord(elem(i,2),2);coord(elem(i,3),2);coord(elem(i,1),2)];
plot(xplt,yplt,'g','linewidth',2)
  
hold on


end

