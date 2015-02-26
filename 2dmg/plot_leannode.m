function plot_leannode( leanmesh_matlab, nodeID)
% plots a matlab output of a xf_LeanMesh

run(leanmesh_matlab)

figure(1)
i = nodeID+1
plot(coord(i,1),coord(i,2),'r*','linewidth',2)
  
hold on


end

