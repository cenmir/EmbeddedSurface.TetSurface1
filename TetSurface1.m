classdef TetSurface1 < matlab.mixin.Copyable
    %TetSurface1 Tetrahedral embedded P1 Surface
    %   o = TetSurface1(H, phi, level)
    %   Input:
    %   H               Hex1 mesh object
    %   phi             surface function m-by-1
    %   level           Iso surface level
    % 
    %   Created by Mirza Cenanovic 2015-04-23
    
    properties
        Surface
        Count
        CutElements
        nCutEle
        
    end
    
    
    methods
        
        function o = TetSurface1(varargin)
            %   o = TetSurface1(H, phi, level)
            %   Input:
            %   H               Hex1 mesh object
            %   phi             surface function m-by-1
            %   level           Iso surface level
            %% Check if proper data has been sent
            narginchk(3,3);
            
            T = varargin{1};
            
            if ~(isa(T,'Mesh.Tet1'))
                error('Object is not Mesh.Tet1!')
            end
            
            if isa(varargin{2},'function_handle')
                surfaceFunction = varargin{2};
                phi = surfaceFunction(T.XC, T.YC, T.ZC);
            else
                error('Wrong type of phi!')
            end
            
            level = varargin{3};
            
            %% Cut element
%             disp(['Extracting ',num2str(size(phi,2)),' surface(s)...'])
%             tic
            surf = CutP1(T, phi, level);
%             toc
            o.Count = length(surf);
            o.Surface = EmbeddedSurface.ImplicitSurface.empty(0);
            AllCutElements = {surf.CutElements};
            for ic = 1:(o.Count)
                o.Surface(ic) = EmbeddedSurface.ImplicitSurface(surf(ic).surfh, surf(ic).SurfacePoints, surf(ic).CutElements).Triangulate(T);
                o.Surface(ic).SurfaceFunction = surfaceFunction;
                o.Surface(ic).IsoLevel = level;
            end
            o.CutElements = unique(vertcat(AllCutElements{:}));
            o.nCutEle = length(o.CutElements);
            
        end
        
        function h = vizP1Surf(o, varargin)
            [h(1).fig, xf1] = xfigure('TetSurface1');
            for isurf = 1:o.Count
                Points = o.Surface(isurf).Points;
                Connectivity = o.Surface(isurf).Connectivity;
                ntri = length(o.Surface(isurf).surfh);
                nCutEle = o.Surface(isurf).nCutEle;
                h(isurf).fig = h(1).fig;
                
                
                FV.Vertices = Points;
                FV.Faces = Connectivity;
                h(isurf).patch = patch(FV,'FaceColor','c','FaceLighting','gouraud');
                axis equal; hold on;
                h(isurf).light = light;
                % shading interp
                [az,el]=view();
                h(isurf).az = az;
                h(isurf).el = el;
                
                
                isNotEmpty = @(x)~isempty(x);
                p = inputParser;
                addParameter(p,'showMesh',[],isNotEmpty);
                parse(p,varargin{:});
                T = p.Results.showMesh;
                
                param = 'showMesh';
                if isenabled(param,varargin)
                    h(isurf).mesh = showMesh(T);
                end
                
                title([num2str(o.Count),' surfaces '])
                
                
                helpText = xf1.uiTextHelp.String;
                helpText{end+1} = 's - Toggle shading';
                helpText{end+1} = 'l - Toggle light';
                helpText{end+1} = 'e - Toggle edges';
                helpText{end+1} = 'v - Toggle vertex normal';
                helpText{end+1} = 'f - Toggle face normal';
                xf1.uiTextHelp.String = helpText;
                xf1.uiTextHelp.Position = xf1.uiTextHelp.Position + [0,0,0,60];
                
                TR = triangulation(FV.Faces, FV.Vertices);
                FN = faceNormal(TR);
                VN = vertexNormal(TR);
                IC = incenter(TR);
                h(isurf).FQ = quiver3(IC(:,1),IC(:,2),IC(:,3),FN(:,1),FN(:,2),FN(:,3),1,'Color','b');
                h(isurf).FQ.Visible = 'off';
                
                h(isurf).VQ = quiver3(Points(:,1),Points(:,2),Points(:,3), VN(:,1),VN(:,2),VN(:,3),1,'Color','r');
                h(isurf).VQ.Visible = 'off';
                
                h(isurf).fig.KeyPressFcn = {@GKPFVizMesh,xf1,h};
                
            end
            
            
        end
        
        [Surface , Points, CutElements , nCutEle  ] = CutP1(T, phi, level)
        
        [tri,X,surfh] = TriangulateP1(o,T)
    end
    
    
    
     methods(Hidden)
      function lh = addlistener(varargin)
         lh = addlistener@handle(varargin{:});
      end
      function notify(varargin)
         notify@handle(varargin{:});
      end
      function delete(varargin)
         delete@handle(varargin{:});
      end
      function Hmatch = findobj(varargin)
         Hmatch = findobj@handle(varargin{:});
      end
      function p = findprop(varargin)
         p = findprop@handle(varargin{:});
      end
      function TF = eq(varargin)
         TF = eq@handle(varargin{:});
      end
      function TF = ne(varargin)
         TF = ne@handle(varargin{:});
      end
      function TF = lt(varargin)
         TF = lt@handle(varargin{:});
      end
      function TF = le(varargin)
         TF = le@handle(varargin{:});
      end
      function TF = gt(varargin)
         TF = gt@handle(varargin{:});
      end
      function TF = ge(varargin)
         TF = ge@handle(varargin{:});
      end
   end
    
end

function GKPFVizMesh(src,evnt,xfigure_This,h)
    %GKPF(src,evnt,xfigure_This,h)
    %   h is a struct that contains user defined graphical handles to be used
    %   in this function.
    %   xfigure_This is a handle to xfigure, returned by xfigure.

    %% Do not touch
    %Leave this be, we need the standard keypress function!
    xfigure_KPF(src, evnt, xfigure_This); %Do not touch
    
    
    %% Add own keys below
    if strcmpi(evnt.Character,'e')
        if strcmpi(get(h.patch,'EdgeColor'),'none')
            set(h.patch,'EdgeColor','k')
        else
            set(h.patch,'EdgeColor','none')
        end
    end

    if strcmpi(evnt.Character,'s')
        for i = 1:length(h)
            if strcmpi(get(h(i).patch,'FaceLighting'),'flat')
                set(h(i).patch,'FaceLighting','gouraud')
            else
                set(h(i).patch,'FaceLighting','flat')
            end
        end
    end
    
    if strcmpi(evnt.Character,'l')
        for i = 1:length(h)
            for ilight = h(i).light
                if strcmpi(ilight.Visible,'on')
                    ilight.Visible = 'off';
                else
                    ilight.Visible = 'on';
                end
            end
        end
    end
    
    if strcmpi(evnt.Character,'v')
        for i = 1:length(h)
            if strcmpi(h(i).VQ.Visible, 'off')
                h.VQ.Visible = 'on';
            else
                h.VQ.Visible = 'off';
            end
        end
    end
    
    if strcmpi(evnt.Character,'f')
        if strcmpi(h.FQ.Visible, 'off')
            h.FQ.Visible = 'on';
        else
            h.FQ.Visible = 'off';
        end
    end


end
function mesh = showMesh(T)
    ele = 1:size(T.Connectivity,1);
    ele = ele(:);
    fele = [4*ele-3;4*ele-2;4*ele-1;4*ele-0;];
    mesh = patch(T.XC(T.Faces(fele(:),:)'),T.YC(T.Faces(fele(:),:)'),T.ZC(T.Faces(fele(:),:)'),'w','FaceColor','none');
end
