function Surface = CutP1(T, phi, level)

    nodes = T.Connectivity;

    sizePhi = size(phi); % nnod x nsurf
    nsurf = sizePhi(2);
    
    
    Surface(nsurf).surfh = [];
    Surface(nsurf).SurfacePoints = [];
    Surface(nsurf).CutElements = [];
    Surface(nsurf).nCutEle = [];

    for isurf = 1:nsurf
        iPhi = phi(:,isurf);
        
        % Find Cut elements
        CutEle = (1:size(nodes,1))';
        rind1 = all(iPhi(nodes)>0,2);
        rind2 = all(iPhi(nodes)<0,2);
        rind = any([rind1,rind2],2);
        CutEle(rind) = [];
        nCutEle = length(CutEle);
        if nCutEle == 0
            error('No elements are cut by the surface phi!')
        end
        
        Surface(isurf).CutElements = CutEle;
        Surface(isurf).nCutEle = nCutEle;
        % preallocate space for coords
        % Maximum size is two triangles per Tet
        % Set all to nan rather than zero to be able to find empty elements later
        Surface(isurf).SurfacePoints = ones(nCutEle*2,3)*nan;
        Surface(isurf).surfh(nCutEle*2).iel = [];
        
        
        
        % Loop over all cut elements and extract zero level points
        Surface = ExtractSurface(isurf,iPhi,Surface,T,level);
    end
end

function Surface = ExtractSurface(iSurf,phi,Surface,T,level)
    xnod = T.XC; ynod = T.YC; znod = T.ZC;
    CutEle = Surface(iSurf).CutElements;
    
    %Initialize triangle counter. Number of triangles
    nTriEle = 1;
    
    %% Loop over all cut elements and extract zero level points
    faultyElements = NaN(length(CutEle),1);
    c = 0;
    for iel = CutEle'
        c = c+1;
        edges = T.Element(iel).edges;
        
        %% Loop over all edges
        Pt = NaN(6,3);
        normal = Pt;
        Q = phi(edges);
        rind1 = all(Q>0,2);
        rind2 = all(Q<0,2);
        rind3 = all(Q==0,2);
        rind = any([rind1,rind2,rind3],2);
        ie = (1:6)';
        cutEdge = ie;
        cutEdge(rind) = [];
        
        E0 = edges(cutEdge,:);
        Q1 = phi(E0);
        ind0 = Q1==0;
        E0ind = E0(ind0);
        if ~isempty(E0ind) && (length(unique(E0ind)) == 1 || length(unique(E0ind)) == 2)
            faultyElements(c) = c;
            continue
        end
        for j = cutEdge'
            ivv = edges(j,:);
            xc = xnod(ivv); yc = ynod(ivv); zc = znod(ivv);
            XC =[xc,yc,zc];
            q  = phi(ivv);
            if sum(sign(q)) == 0 % sign(q) is a 1-by-2 with elements 0,1 or -1
                [u,k] = sort(q);
                Xi = XC(k,:);
                n = Xi(2,:)-Xi(1,:);
                normal(j,:) = n/norm(n);
                xq = Xi(1,:)+(Xi(2,:)-Xi(1,:)).*((level-u(1))/(u(2)-u(1)));
                %                 plot3(xq(1),xq(2),xq(3),'b*')
                Pt(j,:) = xq;
            elseif any(q == 0)
                ind = q==0;
                Pt(j,:) = XC(ind,:);
                [~,k] = sort(q);
                Xi = XC(k,:);
                n = Xi(2,:)-Xi(1,:);
                normal(j,:) = n/norm(n);
            end
        end
        Pt1 = Pt(~any(isnan(Pt),2),:);
        P = unique(Pt1,'rows');
        
        if size(P,1) < 3
            continue
        end
        
        normal = mean(normal(all(~isnan(normal),2),:),1);
        normal = normal/norm(normal);
        
         %% Extract triangles
        if size(P,1) == 3 % Triangle
            [Surface, nTriEle] = assembleElement(P,Surface,iSurf,nTriEle,normal,iel,1);

        elseif size(P,1) == 4 % Quadliteral
            [Surface, nTriEle] = assembleElement(P,Surface,iSurf,nTriEle,normal,iel,2);

        else
            patch(P(:,1),P(:,2),P(:,3),'k')
            error('Unknown element')
        end
        
    end
    faultyElements = faultyElements(~isnan(faultyElements));
    faultyElements = unique(faultyElements);
    
    CutEle(faultyElements) = [];
    
    
    Surface(iSurf).SurfacePoints = Surface(iSurf).SurfacePoints(~(sum(isnan(Surface(iSurf).SurfacePoints),2) == 3),:);
    nele = length([Surface(iSurf).surfh.iel]);
    Surface(iSurf).surfh(nele+1:end) = [];
    Surface(iSurf).CutElements = CutEle;
    Surface(iSurf).nCutEle = length(CutEle);
    
    %% Viz surface
    
%     x = Surface(iSurf).SurfacePoints(:,1);
%     y = Surface(iSurf).SurfacePoints(:,2);
%     z = Surface(iSurf).SurfacePoints(:,3);
%     xnod = reshape(x,3,[]);
%     ynod = reshape(y,3,[]);
%     znod = reshape(z,3,[]);
%     xfigure(1); hold on;
%     patch(xnod,ynod,znod,'c')
%     axis equal;
%     view(3)
%     hold off
    
    
end

function [Surface, nTriEle] = assembleElement(P,Surface,iSurf,nTriEle,normal,iel,ntri)

    if size(P,1) == 3
        %% 
        Xp = P;
        tt = [1,2,3];
        v1 = Xp(2,:)-Xp(1,:);
        v2 = Xp(3,:)-Xp(1,:);
        n = cross(v1,v2);
        n = n/norm(n);
        if normal*n' < 0
            tt = [tt(1),tt(3),tt(2)];
            n = -n;
        end

        % Add to surfh
        Surface(iSurf).surfh(nTriEle).iel  = iel;
        Surface(iSurf).surfh(nTriEle).Xe = Xp(tt,:);
        Surface(iSurf).surfh(nTriEle).xp = Xp(tt,1);
        Surface(iSurf).surfh(nTriEle).yp = Xp(tt,2);
        Surface(iSurf).surfh(nTriEle).zp = Xp(tt,3);
        Surface(iSurf).surfh(nTriEle).faceNormal = n;
        Surface(iSurf).surfh(nTriEle).ElementNormal = normal;

        % Assembly
        lo = (nTriEle)*3-2;
        up = (nTriEle)*3;
        Surface(iSurf).SurfacePoints(lo:up,1:3) = Xp(tt,:);
    else
        %%
        % Rotate the polygon down to the x-y plane
        k1= RotatePointsToPlane(P);

        % New order of polygon points
        Xe = P(k1,:);

        % Split into triangles
        tt = [ones(ntri,1),(2:ntri+1)',(3:ntri+2)'];

        for itri = 1:ntri
            % triangle itri
            ti = tt(itri,:);
            Xp = Xe(ti,:);
            v1 = Xp(2,:)-Xp(1,:);
            v2 = Xp(3,:)-Xp(1,:);
            n = cross(v1,v2);
            n = n/norm(n);
            if normal*n' < 0
                ti = [ti(1),ti(3),ti(2)];
                n = -n;
            end
            % Add to surfh

            Surface(iSurf).surfh(nTriEle+itri-1).iel  = iel;
            Surface(iSurf).surfh(nTriEle+itri-1).Xe = Xe(ti,:);
            Surface(iSurf).surfh(nTriEle+itri-1).xp = Xe(ti,1);
            Surface(iSurf).surfh(nTriEle+itri-1).yp = Xe(ti,2);
            Surface(iSurf).surfh(nTriEle+itri-1).zp = Xe(ti,3);
            Surface(iSurf).surfh(nTriEle+itri-1).faceNormal = n;
            Surface(iSurf).surfh(nTriEle+itri-1).ElementNormal = normal;

            % Assembly
            lo = (nTriEle+itri-1)*3-2;
            up = (nTriEle+itri-1)*3;
             Surface(iSurf).SurfacePoints(lo:up,1:3) = Xe(ti,:);
        end
    end



    nTriEle = nTriEle + ntri;

end

function [k1] = RotatePointsToPlane(Xe)
origin = Xe(1,:);
localz = cross(Xe(2,:)-origin, Xe(3,:)-origin);
unitz = localz/norm(localz,2);
%calculate local x vector in plane
localx = Xe(2,:)-origin;
unitx = localx/norm(localx,2);
%calculate local y
localy = cross(localz, localx);
unity = localy/norm(localy,2);
TM = [unitx(:), unity(:), unitz(:), origin(:); 0 0 0 1];
CM = [Xe, ones(size(Xe,1),1)];
Xe2D = TM \ CM';
Xe2D = Xe2D(1:3,:)';
% Measure the curvature
sx = max(Xe2D(:,1))-min(Xe2D(:,1));
sy = max(Xe2D(:,2))-min(Xe2D(:,2));
sz = max(Xe2D(:,3))-min(Xe2D(:,3));
curvature = (sz / mean([sx,sy]))*100;


% New order of polygon points
k = convhull(Xe2D(:,1),Xe2D(:,2));
k1 = k(1:end-1);

if curvature > 10
    error(['excessive curvature (',num2str(curvature),') in cut elements! Increase number of elements!'])
    patch(Xe(:,1),Xe(:,2),Xe(:,3),'r','EdgeColor','none'),hold on;
    plot3(Xe(:,1),Xe(:,2),Xe(:,3),'-ok'); hold off;
end

end
