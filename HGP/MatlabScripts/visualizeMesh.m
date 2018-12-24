function [ ] = visualizeMesh( mesh )

%mesh figure
    figure
    fig1 = tsurf(mesh.F, mesh.V);

    hold on
    axis equal
    grid off
    axis off
    axis vis3d
    light
    if ~isempty(mesh.conesIndices)

        colors = jet(max(mesh.Visualize.seam3d(:,3)));
        colors = colors(randperm(max(mesh.Visualize.seam3d(:,3))),:);
        for i=1:length(mesh.Visualize.seam3d)
            P1 = mesh.V(mesh.Visualize.seam3d(i,1),:);
            P2 = mesh.V(mesh.Visualize.seam3d(i,2),:);
            pts = [P1; P2];
            plot3(pts(:,1), pts(:,2), pts(:,3),'Color',colors(mesh.Visualize.seam3d(i,3),:),'LineWidth',3)
        end

        conesPlot = mesh.V(mesh.conesIndices,:);
        colors2 = jet(length(mesh.conesIndices));
        colors2 = colors2(randperm(length(mesh.conesIndices)),:);

        scatter3(conesPlot(:,1),conesPlot(:,2),conesPlot(:,3),90,colors2,'filled')
        light('Position',[-1 0 0],'Style','local')
        set(fig1,'FaceColor',[0.8 0.8 0.8])

        light('Position',[0 0 0],'Style','local')
    end

%UV figure
    figure
    hold on
    for i=1:length(mesh.halfEdges)
        P1 = mesh.UV(mesh.halfEdges(i,1),:);
        P2 = mesh.UV(mesh.halfEdges(i,2),:);
        P3 = mesh.UV(mesh.halfEdges(i,3),:);
        pts = [P1; P2];
        plot(pts(:,1), pts(:,2),'k')
        pts = [P2; P3];
        plot(pts(:,1), pts(:,2),'k')
        pts = [P3; P1];
        plot(pts(:,1), pts(:,2),'k')  
    end
    
    if ~isempty(mesh.conesIndices)
        for i=1:length(mesh.Visualize.seam2d)
            P1 = mesh.UV(mesh.Visualize.seam2d(i,1),:);
            P2 = mesh.UV(mesh.Visualize.seam2d(i,2),:);
            pts = [P1; P2];
            plot(pts(:,1), pts(:,2),'Color',colors(mesh.Visualize.seam2d(i,3),:),'LineWidth',2)
        end
        axis equal

        for i=1:length(conesPlot)
           conesToPlot = unique(mesh.Visualize.conesMapMatrix(:,i)); 
           if conesToPlot(1)==0
               conesToPlot(1) = [];
           end
           scatter(mesh.UV(conesToPlot,1),mesh.UV(conesToPlot,2),40,colors2(i,:),'filled')
        end    
    end

end

