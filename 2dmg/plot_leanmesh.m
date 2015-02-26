function plot_leanmesh( leanmesh_matlab, pltfront_orient, plt_circumcircle,fignum)
% plots a matlab output of a xf_LeanMesh

run(leanmesh_matlab)

h=figure(fignum)
close(h)
figure(fignum)
for i=1:size(face,1)
    xplt=[coord(face(i,1),1);coord(face(i,2),1)];
    yplt=[coord(face(i,1),2);coord(face(i,2),2)];
    plot(xplt,yplt,'k-')
    hold on
end

%for i=1:size(elem,1)
%    xplt=[coord(elem(i,1),1);coord(elem(i,2),1);coord(elem(i,3),1);coord(elem(i,1),1)];
%    yplt=[coord(elem(i,1),2);coord(elem(i,2),2);coord(elem(i,3),2);coord(elem(i,1),2)];
%    plot(xplt,yplt,'b-')
%    hold on    
%end

if (plt_circumcircle == 1)
    for i=1:size(elem,1)
        xplt=[coord(elem(i,1),1) coord(elem(i,2),1) coord(elem(i,3),1)];
        yplt=[coord(elem(i,1),2) coord(elem(i,2),2) coord(elem(i,3),2)];
        [r,cn]=circumcircle([xplt;yplt],1);
        hold on
    end
end

%for i=1:size(front,1)
%    xplt=[coord(front(i,1),1);coord(front(i,2),1)];
%    yplt=[coord(front(i,1),2);coord(front(i,2),2)];
%    plot(xplt,yplt,'r--','linewidth',2)
%    hold on
%end

if (pltfront_orient == 1)
    for i=1:size(front,1)
        xplt=[coord(front(i,1),1)];
        yplt=[coord(front(i,1),2)];
        uplt=[coord(front(i,2),1)-coord(front(i,1),1)];
        vplt=[coord(front(i,2),2)-coord(front(i,1),2)];
        quiver(xplt,yplt,uplt,vplt,'m-','linewidth',2)
        hold on
    end
end

end

