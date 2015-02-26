function plot_leanface( leanmesh_matlab, iface)
% plots a matlab output of a xf_LeanMesh

run(leanmesh_matlab)

figure(1)
i = iface+1
face(i,:)
xplt=[coord(face(i,1),1)];
yplt=[coord(face(i,1),2)];
uplt=[coord(face(i,2),1)-coord(face(i,1),1)];
vplt=[coord(face(i,2),2)-coord(face(i,1),2)];
quiver(xplt,yplt,uplt,vplt,'y-','linewidth',2)
hold on


end

