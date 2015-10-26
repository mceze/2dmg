function plot_leanmesh( leanmesh_matlab, pltfront_orient, plt_circumcircle,fignum,plt_elemnum)
% plots a matlab output of a xf_LeanMesh

run(leanmesh_matlab)

h=figure(fignum)
hold off
figure(fignum)
for i=1:size(face,1)
    xplt=[coord(face(i,1),1);coord(face(i,2),1)];
    yplt=[coord(face(i,1),2);coord(face(i,2),2)];
    xc(i) = mean([coord(face(i,1),1);coord(face(i,2),1)]);
    yc(i) = mean([coord(face(i,1),2);coord(face(i,2),2)]);
    facenum(i,:)=cellstr(num2str(i-1)); 
    plot(xplt,yplt,'b-')
    hold on
end
if (plt_elemnum == 1)
    for i=1:size(elem,1)
        xplt(i)=mean([coord(elem(i,1),1);coord(elem(i,2),1);coord(elem(i,3),1)]);
        yplt(i)=mean([coord(elem(i,1),2);coord(elem(i,2),2);coord(elem(i,3),2)]);
        elemnum(i,:)=cellstr(num2str(i-1)); 
    end
    text(xplt',yplt',elemnum(:,1),'color','r')
    text(xc',yc',facenum(:,1),'color','k')
    hold on 
end

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

