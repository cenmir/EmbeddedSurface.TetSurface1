function [tri,X,surfh] = TriangulateP1(o,H,showBBox)
%TRIANGULATEP1 Creates a triangulation from a set of P1 triangles
%   Triangulates the P1 Surface triangle set surfX
%   
%   [tri,X] = TriangulateP1(H,showBBox)
%   H is the Hex object
%   showBBox is used for debugging, set to 1 to show
%   tri is the triangulation matrix m-by-3
%   X is the set of triangle coordinates same as (surfX)
%   Usage:
%   T.TriangulateP1( ... )
%   or
%   TriangulateP1(T, ... )
%   T is the Hex1Mesh class
% Strategy
% The overal strategy is to loop over each set of triangle points
% and compare each point to the previous added list of points. If
% the current point exists in the list, then it was previously added
% and its position in the list is added to the list. Otherwise the
% number of nodes is incremented and that number is added to the
% triangulation list.
%   Example:
%   First set is initiated to [1,2,3]
%   In the next iteration a point in the set is already in the list,
%   the point number is determined to be 2, so two new nodes are
%   added:
%   [1,2,3;
%    4,5,3]
%   And so on...
%
% To improve performance a bounding box is sweeping over every triangle
% (every third point of the list). The bounding box ecapsulates a number of
% points which are the search space. So instead of for each point
% compare the distance to every other point, we compute the distance to
% a constant set of other points.



    surfX = o.Points;
    surfh = o.Surface;


    %% sort surfh
    % ntri = length(surfh);
    % XM = zeros(ntri,3);
    % for itri = 1:ntri
    %     XM(itri,:) = mean(surfh(itri).Xe,1);
    % end
    % [~,index] = sortrows(XM);
    % surfh = surfh(index);
    % X = [surfh(:).Xe];
    % x = X(1:3,1:3:end);
    % x = reshape(x,[],1);
    % y = X(1:3,2:3:end);
    % y = reshape(y,[],1);
    % z = X(1:3,3:3:end);
    % z = reshape(z,[],1);
    % surfX = [x,y,z];


    X = surfX;
    maxXC = max(X(:,1));
    minXC = min(X(:,1));
    maxYC = max(X(:,2));
    minYC = min(X(:,2));
    maxZC = max(X(:,3));
    minZC = min(X(:,3));

    DXC = (maxXC-minXC)/H.nxe;
    DYC = (maxYC-minYC)/H.nye;
    DZC = (maxZC-minZC)/H.nze;
    
    %% Determine size of bounding box
%     [DXC,DYC,DZC] = computeBoundingBox(X,DXC,DYC,DZC);
    
    disp('Triangulating surface...')
    tic
    [tri,surfh] = PrivateTriangulate(X,showBBox,surfh);
    toc
    %% Weld
    disp('Welding points...')
    tic
    [X,tri] = weldPoints(X,tri);
    toc
end

function [P,T] = weldPoints(X,tri)
    mt = 0;
    map = zeros(size(X,1),2);
    % Renumbering nodes by creating map
    % Map which nodes to replace with other nodes
    % Second column idices replaces first column
    c = 0;
    for iel = 1:size(tri)
        itri = tri(iel,:);
        if iel == 1
            map(1:3,:)=[(1:3)',(1:3)'];
            c = 3;
            mt = 4;
        else
            for indt=1:3
                it = itri(indt);
                if it >= mt && ~any(it==map(:,1))
                    map(c+1,:) = [it,mt];
                    mt=mt+1;
                    c=c+1;
                end
                
            end
            
        end
    end
    map = map(all(map~=0,2),:);
    
    for im = 1:size(map,1)
        tri(tri==map(im,1)) = map(im,2);
    end
    T = tri;
    P =  X(map(:,1),:);
end

function bandInds = computeBoundingBox(X,iP,ShowBox)
    ipnts = X(iP:iP+2,:);
    xc0 = min(ipnts(:,1));
    xc1 = max(ipnts(:,1));
    yc0 = min(ipnts(:,2));
    yc1 = max(ipnts(:,2));
    zc0 = min(ipnts(:,3));
    zc1 = max(ipnts(:,3));
    bandIndsx = X(:,1)>=xc0 & X(:,1)<=xc1;
    bandIndsy = X(:,2)>=yc0 & X(:,2)<=yc1;
    bandIndsz = X(:,3)>=zc0 & X(:,3)<=zc1;
    bandInds = all([bandIndsx,bandIndsy,bandIndsz],2);
    if ShowBox
        drawBBox(X,bandInds,xc1,xc0,yc1,yc0,zc1,zc0,iP)
    end
end

function [tri,surfh] = PrivateTriangulate(X,ShowBox,surfh)
    %% Begin triangulation
    nnod = size(X,1);
    nele = nnod/3;
    tri = zeros(nele, 3);
    itri = 1;
    for iP = 1:3:size(X,1)
        if itri == 1
            tri(1,1:3) = [1,2,3];
            xe = X(tri(itri,:),:);
            n = cross( xe(2,:)-xe(1,:) , xe(3,:)-xe(1,:) );
            n = n/norm(n);
            surfh(itri).facenormal = n;
            itri = itri +1;
        else
            %% Bounding Box
            bandInds = computeBoundingBox(X,iP,ShowBox);
            PrevPoints = [1:(iP+2-1)]';
            ninds = sum(bandInds);
            if ( length(PrevPoints) >= ninds )
                searchSpaceBand = X(bandInds,:);
                SearchIndsBand = find(bandInds); 
            end
            %% Loop over the three next points in the set of triangles
            % the local triangle node index is set to 1
            tind = 1;
            for iPnt = iP:iP+2
                % If the neighborhood point space is smaller than the space of
                % points we have found so far, then define the searchSpace as
                % the neighborhood points. Otherwise define the searchSpace as
                % the points we have found so far.
                % This is done because in the first few iterations the points
                % found so far are fewer than all the neighborhood points, but
                % later the neighborhgood is pretty mush constant where as the
                % points found so far would only increase.

                PrevPoints = [1:iPnt-1]';
                if ( length(PrevPoints) > ninds )
                    searchSpace = searchSpaceBand;
                    SearchInds = SearchIndsBand;
                else
                    searchSpace = X(PrevPoints,:);
                    SearchInds = PrevPoints;
                end

                % Chose a point in the current set of trinagle points to measure the distance from
                Xp = X(iPnt,:);

                % Choose a set of points to measure the distance to. The set of
                % points is defined by the searchSpace.
                P = ones(size(searchSpace,1),1)*Xp;

                % The distance vectors between the chosen point and the searchSpace
                % points.
                D = searchSpace-P;

                % Eucledian distance
                NDist = sqrt(D(:,1).^2+D(:,2).^2+D(:,3).^2);

                % Find the indices to the point that are identical (with room
                % for numerical error)
                indt = SearchInds(NDist<eps*100);


                % if duplicate points exist; set the local triangle index to
                % the index of the searchSpace and increment the triangle
                % index.
                % If no duplicate points are found; set the local triangle
                % index to the index of the current triangle set point
                if ~isempty(indt)
                    indt = indt(1);
                    tri(itri,tind) = indt;
                    tind = tind +1;
                else
                    tri(itri,tind) = iPnt;
                    tind = tind +1;
                end
            end
            xe = X(tri(itri,:),:);
            n = cross( xe(2,:)-xe(1,:) , xe(3,:)-xe(1,:) );
            n = n/norm(n);
            surfh(itri).facenormal = n;
            itri = itri +1;
        end
    end
end

function [tri,normals]=PrivateTriangulateWNormals(X,DXC, DYC, DZC)
    tri = zeros(nele, 3);
    normals = tri;
    itri = 1;
    for iP = 1:3:size(X,1)
        if itri == 1
            tri(1,1:3) = [1,2,3];
            xe = X(tri(itri,:),:);
            n = cross( xe(2,:)-xe(1,:) , xe(3,:)-xe(1,:) );
            n = n/norm(n);
            normals(itri,:)=n;
            %          quiver3(mean(xe(:,1)),mean(xe(:,2)),mean(xe(:,3)),n(1),n(2),n(3),0.1)
            itri = itri +1;
        else
            %% Bounding Box
            ix = X(iP,1);
            iy = X(iP,2);
            iz = X(iP,3);
            xc1 = ix+(DXC);
            xc0 = ix-(DXC);
            yc1 = iy+(DYC);
            yc0 = iy-(DYC);
            zc1 = iz+(DZC);
            zc0 = iz-(DZC);
            bandIndsx = X(:,1)>xc0 & X(:,1)<xc1;
            bandIndsy = X(:,2)>yc0 & X(:,2)<yc1;
            bandIndsz = X(:,3)>zc0 & X(:,3)<zc1;
            bandInds = all([bandIndsx,bandIndsy,bandIndsz],2);
            if ShowBox
                drawBBox(X,bandInds,iP,xc0,yc0,zc0,xc1,yc1,zc1)
            end

            if ( iP*3 >= sum(bandInds) )
                searchSpaceBand = X(bandInds,:);
                SearchIndsBand = find(bandInds);
            end




            %% Loop over the three next points in the set of triangles
            % the local triangle node index is set to 1
            tind = 1;
            for iPnt = iP:iP+2
                % If the neighborhood point space is smaller than the space of
                % points we have found so far, then define the searchSpace as
                % the neighborhood points. Otherwise define the searchSpace as
                % the points we have found so far.
                % This is done because in the first few iterations the points
                % found so far are fewer than all the neighborhood points, but
                % later the neighborhgood is pretty mush constant where as the
                % points found so far would only increase.

                PrevPoints = [1:iPnt-1]';


                ninds = sum(bandInds);
                if ( length(PrevPoints) > ninds )
                    searchSpace = searchSpaceBand;
                    SearchInds = SearchIndsBand;
                else
                    searchSpace = X(PrevPoints,:);
                    SearchInds = PrevPoints;
                end
                %

                % Chose a point in the current set of trinagle points to measure the distance from
                Xp = X(iPnt,:);

                % Choose a set of points to measure the distance to. The set of
                % points is defined by the searchSpace.
                P = ones(size(searchSpace,1),1)*Xp;

                % The distance vectors between the chosen point and the searchSpace
                % points.
                D = searchSpace-P;

                % Eucledian distance
                NDist = sqrt(D(:,1).^2+D(:,2).^2+D(:,3).^2);

                % Find the indices to the point that are identical (with room
                % for numerical error)
                indt = SearchInds(NDist<eps*100);


                % if duplicate points exist; set the local triangle index to
                % the index of the searchSpace and increment the triangle
                % index.
                % If no duplicate points are found; set the local triangle
                % index to the index of the current triangle set point
                if ~isempty(indt)
                    indt = indt(1);
                    tri(itri,tind) = indt;
                    tind = tind +1;
                else
                    tri(itri,tind) = iPnt;
                    tind = tind +1;
                end
            end

            xe = X(tri(itri,:),:);
            n = cross( xe(2,:)-xe(1,:) , xe(3,:)-xe(1,:) );
            n = n/norm(n);
            normals(itri,:)=n;
            %         quiver3(mean(xe(:,1)),mean(xe(:,2)),mean(xe(:,3)),n(1),n(2),n(3),0.1)
            itri = itri +1;


        end

    end
end

function drawBBox(X,bandInds,iP,xc0,yc0,zc0,xc1,yc1,zc1)
    %% Draw bounding box
    BX = [xc0,yc0,zc0;
        xc1,yc0,zc0;
        xc1,yc1,zc0;
        xc0,yc1,zc0;
        xc0,yc0,zc1;
        xc1,yc0,zc1;
        xc1,yc1,zc1;
        xc0,yc1,zc1];
    faces = [1,2,3,4;
        5,6,7,8;
        1,2,6,5;
        2,3,7,6;
        3,4,8,7;
        4,1,5,8];
    xfigure(23);
    hold off
    plot3(BX(:,1),BX(:,2),BX(:,3),'sk');hold on;
    patch('Faces',faces,'Vertices',BX,'FaceColor','none');
    axis equal;

    ssb = X(bandInds,:);
    plot3(ssb(:,1),ssb(:,2),ssb(:,3),'ok');
    
    PrevPoints = [1:(iP+2)];
    pnts = iP:iP+2;
    ssp = X(PrevPoints,:);
    plot3(ssp(:,1),ssp(:,2),ssp (:,3),'*b')
    plot3(X(pnts,1),X(pnts,2),X (pnts,3),'+r','MarkerSize',12)
    legend('','','Local Searchspace','Previous Points','Current Point')
    disp('Press any key or click on the figure to continue')
    w = waitforbuttonpress;
end

function rv = isenabled(mode, varargin)
%   ISENABLED  Checks if mode exists in the cell-array varargin.
    if nargin < 1
        error('No arguments')
    end
    varargin = varargin{:};
    ind = find(strcmpi(varargin,mode), 1);
    if ~isempty(ind)
        rv = 1;
    else
        rv = 0;
    end
end
