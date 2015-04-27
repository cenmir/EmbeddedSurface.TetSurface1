function [surfh, surfX, SurfEle, nCutEle ] = CutP1(T, phi, level)

    nodes = T.Connectivity;
    xnod = T.XC;
    ynod = T.YC;
    znod = T.ZC;



    SurfEle = find(sum(phi(nodes)>0,2) > 0 & sum(phi(nodes)>0,2)< 4);
%     SurfEle = [1:size(nodes,1)]';
%     SurfEle = find(any(phi(nodes) < 0,2) & any(phi(nodes) > 0,2));
    

    %% Viz exact surface
%     [xx,yy,zz] = meshgrid(linspace(0,1,40),linspace(0,1,40),linspace(0,1,40));
%     R = 0.45;
%     xc = 0.5; yc = 0.5; zc = 0.5;
%     surfaceFunction = @(x,y,z) ((x-xc).^2+(z - zc).^2+(y - yc).^2).^.5-R;
%     vv = surfaceFunction(xx,yy,zz);
%     fv = isosurface(xx,yy,zz,vv,0);
% 
%     h = T.vizMesh(SurfEle); hold on;
%     h.patch.EdgeAlpha = .3;
% 
%     p = patch(fv); isonormals(xx,yy,zz,vv,p);
%     p.FaceColor = 'red';
%     p.EdgeColor = 'none';
%     daspect([1,1,1])
%     view(3); axis tight
%     camlight
%     lighting gouraud
%     p.FaceAlpha = 0.3;


    %%
    nsurfEle = length(SurfEle);
    if nsurfEle == 0
        error('No elements are cut by the surface phi!')
    end
    
%     nprec = 1e-8;
%     isequalAbs = @(x,y) ( abs(x-y) <= nprec );
%     isequalPoint = @(a,b) (sum(isequalAbs(a, b))==3);
%     isodd = @(x) mod(x,2);
    
    % preallocate space for coords
    % Maximum size is two triangles per Tet
    % Set all to nan rather than zero to be able to find empty elements later
    surfX = ones(nsurfEle*2,3)*nan;
    surfh(nsurfEle*2).iel = [];

    nTriEle = 1;% number of triangles found
    
%     xfigure(1);clf;hold on; axis equal;
    % Loop over all cut elements and extract zero level points
    for iel = SurfEle'
        
        iv = nodes(iel,:);
        %% Viz element
%         ele = iel;
%         fele = [4*ele-3;4*ele-2;4*ele-1;4*ele-0];
%         patch(T.XC(T.Faces(fele(:),:)'),T.YC(T.Faces(fele(:),:)'),T.ZC(T.Faces(fele(:),:)'),'w','FaceColor','none');
        
        %% new
%         u = phi(iv);
%         [u,~] = sort(u);
%         %determine if the element is cut
%         if u(1) <= level && u(end) >= level
%             %% Compute cut points
%             % Sides of the tetrahedral
%             s = zeros(4,3);
%             s(1,:) = [iv(1),iv(2),iv(4)];
%             s(2,:) = iv(2:4);
%             s(3,:) = [iv(1),iv(3),iv(4)];
%             s(4,:) = iv(1:3);
%             
%             ncut = 0;
%             fp = 1;
%             tP = zeros(4*2,3);
%             normal = NaN(4*2,3);
%             for iside = 1:4
%                 is = s(iside,:);
%                 c = phi(is);
%                 [u,k] = sort(c);
%                 u3 = u(3); u2 = u(2);  u1 = u(1); %u3 is max, u2 is min
%                 k1 = k(1); k2 = k(2); k3 = k(3);
%                 p1 = T.Points(is(k1),:); p2 = T.Points(is(k2),:); p3 = T.Points(is(k3),:);
%                 
%                 if u3 == u2 && u2 == 1 && u3 == u1 % The whole side of the element lies on the interface
%                     % TODO
%                     disp('case 5')
%                     disp('u:')
%                     disp([u1;u2;u3])
%                     disp('The whole side of the element lies on the interface')
%                     error('This case has not been implemented!')
%                 elseif u3 < level && u2 < level && u1 < level
%                     %This side is not cut!
%                 elseif u3 >= u2 && u2 >= u1
%                     if u1 < level && u2 < level
%                         % u3 - u1
%                         xp1 = p1 + ( p3-p1 ).*((level-u1)/(u3-u1));
%                         % u3 - u2
%                         xp2 = p2 + ( p3-p2 ).*((level-u2)/(u3-u2));
%                         ncut = ncut +1;
%                         normal(fp:fp+1,:) = [p1-p3;p2-p3];
%                         tP(fp:fp+1,:) = [xp1;xp2];
%                         fp = fp + 2;
%                     elseif u1 < level && u2 > level
%                         % u1 - u3
%                         xp1 = p1 + ( p3-p1 ).*((level-u1)/(u3-u1));
%                         % u1 - u2
%                         xp2 = p1 + ( p2-p1 ).*((level-u1)/(u2-u1));
%                         ncut = ncut +1;
%                         normal(fp:fp+1,:) = [p3-p1;p2-p1];
%                         tP(fp:fp+1,:) = [xp1;xp2];
%                         fp = fp + 2;
%                     elseif u2 == level && u1 < 0 && u3 > 0
%                         disp('case 6')
%                         disp('u2 == level && u1 < 0 && u3 > 0')
%                         disp([u1;u2;u3])
%                         error('Not implemented')
%                     end
%                 elseif u1 == u2 && (u1 <= level && u3 > level)
%                     disp('case 3')
%                     disp('u1 == u2 is parallel')
%                     disp([u1;u2;u3])
%                     error('Not implemented')
%                 elseif u2 == u3 && (u2 >= level && u1 < level)
%                     disp('case 4')
%                     disp('u2 == u3 is parallel')
%                     disp([u1;u2;u3])
%                     error('Not implemented')
%                 else
%                     disp([u1;u2;u3])
%                     error('This case has not been identified')
%                 end
%             end
%             %% normal
%             normal = mean(normal(all(~isnan(normal),2),:),1);
%             normal = normal/norm(normal);
%             
%             %% Assemble cut points
%             switch ncut
%                 case 3 % Triangles
%                     %%
%                     Xe = zeros(3,3); %Pre-allocate the triangle coordinates
%                     Xe(1:2,:) = tP(1:2,:); % The first two points of one of the sides of the tet are the first two points of the triangle
%                     for i = 3:6 %We only need to look for one more point that is not already in the list
%                         xpi = tP(i,:);
%                         if ~(sum(isequalAbs(xpi, Xe(1,:)))==3) && ~(sum(isequalAbs(xpi, Xe(2,:)))==3) % if the point xpi is not in the list, add it and exit the loop.
%                             Xe(3,:) = xpi;
%                             break
%                         end
%                     end
%                     
%                     % Normal orientation
%                     tt = [1,2,3];
%                     v1 = Xe(2,:)-Xe(1,:);
%                     v2 = Xe(3,:)-Xe(1,:);
%                     n = cross(v1,v2);
%                     n = n/norm(n);
%                     if normal*n' < 0
%                         tt = [tt(1),tt(3),tt(2)];
%                         n = -n;
%                     end
%                     Xe = Xe(tt,:);
%                     patch(Xe(:,1),Xe(:,2),Xe(:,3),'b')
%                     quiver3(mean(Xe(:,1)),mean(Xe(:,2)),mean(Xe(:,3)),n(1),n(2),n(3),0.1,'Color','k')
%                     
%                     
%                     % Add to surfh
%                     surfh(nTriEle).iel  = iel;
%                     surfh(nTriEle).Xe = Xe;
%                     surfh(nTriEle).xp = Xe(:,1);
%                     surfh(nTriEle).yp = Xe(:,2);
%                     surfh(nTriEle).zp = Xe(:,3);
%                     surfh(nTriEle).faceNormal = n;
%                     surfh(nTriEle).ElementNormal = normal;
%                     
%                     %Assembly
%                     lo = nTriEle*3-2;
%                     up = nTriEle*3;
%                     
%                     surfX(lo:up,1:3) = Xe;
%                     nTriEle = nTriEle +1; % One triangle added
%                     
%                 case 4 % Quadliterals
%                     %%
%                 
%             end
%             
%         end
        
        %% OLD
        edges = T.Element(iel).edges;
        %% Viz edges
        %         ed = unique(edges);
        %         qq = phi(ed);
        %         for j = 1:4
        %             str = [num2str(ed(j)),' ',num2str(qq(j))];
        %             text(xnod(ed(j)),ynod(ed(j)),znod(ed(j)),str,'BackgroundColor','w')
        %         end

        %% Loop over all edges
        P = NaN(6,3);
        normal = P;
        for j = 1:size(edges,1)
            ivv = edges(j,:);
            xc = xnod(ivv); yc = ynod(ivv); zc = znod(ivv);
            XC =[xc,yc,zc];
            q  = phi(ivv);
            %             if sign(q(1)) ~= sign(q(2))
            if sum(sign(q)) == 0 % sign(q) is a 1-by-2 with elements 0,1 or -1
                [u,k] = sort(q);
                Xi = XC(k,:);
                n = Xi(2,:)-Xi(1,:);
                normal(j,:) = n/norm(n);
                P(j,:) = Xi(1,:)+(Xi(2,:)-Xi(1,:)).*((level-u(1))/(u(2)-u(1)));
            end
        end

        %         indP = find(~any(isnan(P),2));
        P = P(~any(isnan(P),2),:);



        normal = mean(normal(all(~isnan(normal),2),:),1);
        normal = normal/norm(normal);

%                 quiver3(mean(P(:,1)),mean(P(:,2)),mean(P(:,3)),normal(1),normal(2),normal(3),0.1,'Color','b')
        %         for i = 1:4
        %            text(xnod(iv(i)),ynod(iv(i)),znod(iv(i)),num2str(phi(iv(i))),'BackgroundColor','w')
        %         end
        %         plot3(P(:,1),P(:,2),P(:,3),'b*')
%                 patch(P(:,1),P(:,2),P(:,3),'b');
                

        %% Extract triangles
        if size(P,1) == 3
            %% Triangle
            [surfh, nTriEle, surfX] = assembleElement(P,surfX,surfh,nTriEle,normal,iel,1);

        elseif size(P,1) == 4
            %% Quadliteral
            [surfh, nTriEle, surfX] = assembleElement(P,surfX,surfh,nTriEle,normal,iel,2);

        else
            patch(P(:,1),P(:,2),P(:,3),'k')
            error('Unknown element')
        end



    end
    surfX= surfX(~(sum(isnan(surfX),2) == 3),:);
    nele = length([surfh.iel]);
    surfh(nele+1:end) = [];
    nCutEle = length(SurfEle);


end

function [surfh, nTriEle, surfX] = assembleElement(P,surfX,surfh,nTriEle,normal,iel,ntri)

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
        surfh(nTriEle).iel  = iel;
        surfh(nTriEle).Xe = Xp(tt,:);
        surfh(nTriEle).xp = Xp(tt,1);
        surfh(nTriEle).yp = Xp(tt,2);
        surfh(nTriEle).zp = Xp(tt,3);
        surfh(nTriEle).faceNormal = n;
        surfh(nTriEle).ElementNormal = normal;

        % Assembly
        lo = (nTriEle)*3-2;
        up = (nTriEle)*3;
        surfX(lo:up,1:3) = Xp(tt,:);

%         p = Xp(tt,:);
%         patch(p(:,1),p(:,2),p(:,3),'b')
%         quiver3(mean(Xp(:,1)),mean(Xp(:,2)),mean(Xp(:,3)),n(1),n(2),n(3),0.1,'Color','r')

    else
        %%
        % Rotate the polygon down to the x-y plane
        k1= RotatePointsToPlane(P);

        % New order of polygon points
        Xe = P(k1,:);
        
%         xfigure(423);clf;
%         np = size(Xe,1);
%         for i = 1:4
%             text(Xe(i,1),Xe(i,2),Xe(i,3),num2str(i),'BackgroundColor','w')
%         end
%         plot3(Xe([1:np,1],1),Xe([1:np,1],2),Xe([1:np,1],3),'-k')
%         pause

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
            
%             patch(Xe(ti,1),Xe(ti,2),Xe(ti,3),'r','EdgeColor','k')
%             quiver3(mean(Xe(ti,1)),mean(Xe(ti,2)),mean(Xe(ti,3)),n(1),n(2),n(3),0.1,'Color','b')

            % Add to surfh

            surfh(nTriEle+itri-1).iel  = iel;
            surfh(nTriEle+itri-1).Xe = Xe(ti,:);
            surfh(nTriEle+itri-1).xp = Xe(ti,1);
            surfh(nTriEle+itri-1).yp = Xe(ti,2);
            surfh(nTriEle+itri-1).zp = Xe(ti,3);
            surfh(nTriEle+itri-1).faceNormal = n;
            surfh(nTriEle+itri-1).ElementNormal = normal;

            % Assembly
            lo = (nTriEle+itri-1)*3-2;
            up = (nTriEle+itri-1)*3;
            surfX(lo:up,1:3) = Xe(ti,:);

%             p = Xe(tt,:);
%             patch(p(:,1),p(:,2),p(:,3),'r')

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
    warning(['excessive curvature (',num2str(curvature),') in cut elements! Increase number of elements!'])
    patch(Xe(k,1),Xe(k,2),Xe(k,3),'r','EdgeColor','none'),hold on;
    plot3(Xe(k(1:end,1),1),Xe(k(1:end,1),2),Xe(k(1:end,1),3),'-ok'); hold off;
end

end
