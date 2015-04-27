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
        Points
        CutElements
        nCutEle
        nele
        nnod
        Connectivity
    end
    
    
    methods
        
        function o = TetSurface1(T, phi, level)
            %   o = TetSurface1(H, phi, level)
            %   Input:
            %   H               Hex1 mesh object
            %   phi             surface function m-by-1
            %   level           Iso surface level
            %% Check if proper data has been sent
            if ~(isa(T,'Mesh.Tet1'))
                error('Object is not Mesh.Tet1!')
            end
            
            if (length(phi) ~= T.nnod)
                error('Number of elements in phi must match number of vertices')
            end
            
            %% Cut element
            disp('Extracting surface points...')
            tic
            [o.Surface , o.Points, o.CutElements , o.nCutEle  ] = CutP1(T, phi, level);
            toc
            
            %% Triangulate
            [tri,X,surfh] = TriangulateP1(o,T);
            o.Connectivity = tri;
            o.Points = X;
            o.Surface = surfh;
            o.nele = size(tri,1);
            o.nnod = size(X,1);
            
        end
        
        h = vizP1Surf(o, varargin);
        
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