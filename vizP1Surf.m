function h = vizP1Surf(o, varargin)
%VIZP1SURF Visualize P1 surface
%   h = vizP1Surf(T)
%   h = vizP1Surf(T, parameters)
%   parameters:
%   showMesh


if nargin > 1
    [h.fig,xf1] = xfigure(varargin{1});
else
    [h.fig,xf1] = xfigure;
end

axis equal; hold on;
FV.Vertices = o.Points;
FV.Faces = o.Connectivity;
shading interp
h.light = light;
h.patch = patch(FV,'FaceColor','c','FaceLighting','gouraud');
[az,el]=view(136,16);
h.az = az;
h.el = el;

param = 'showMesh';
if isenabled(param,varargin)
    h.mesh = showMesh(o);
end



ntri = length(o.Surface);
nCutEle = o.nCutEle;

title([num2str(ntri),' triangles on ',num2str(nCutEle),' cut Tet1 elements'])


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
h.FQ = quiver3(IC(:,1),IC(:,2),IC(:,3),FN(:,1),FN(:,2),FN(:,3),1,'Color','b');
h.FQ.Visible = 'off';

h.VQ = quiver3(o.Points(:,1),o.Points(:,2),o.Points(:,3), VN(:,1),VN(:,2),VN(:,3),1,'Color','r');
h.VQ.Visible = 'off';

h.fig.KeyPressFcn = {@GKPFVizMesh,xf1,h};

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
        if strcmpi(get(h.patch,'FaceLighting'),'flat')
            set(h.patch,'FaceLighting','gouraud')
        else
            set(h.patch,'FaceLighting','flat')
        end
    end
    
    if strcmpi(evnt.Character,'l')
        for ilight = h.light
            if strcmpi(ilight.Visible,'on')
                ilight.Visible = 'off';
            else
                ilight.Visible = 'on';
            end
        end
    end
    
    if strcmpi(evnt.Character,'v')
        if strcmpi(h.VQ.Visible, 'off')
            h.VQ.Visible = 'on';
        else
            h.VQ.Visible = 'off';
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
    fele = [6*ele-5;6*ele-4;6*ele-3;6*ele-2;6*ele-1;6*ele-0;];
    mesh = patch(T.XC(T.Faces(fele(:),:)'),T.YC(T.Faces(fele(:),:)'),T.ZC(T.Faces(fele(:),:)'),'w','FaceColor','none');
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

function param = getoption(mode, varargin)
%   GETOPTION  Gets the parameter for the option 'mode'
%
%   getoption(mode,varargin) returns the parameter.
%   example:
%
%          varargin = {'contour', 'grid', 'view', [20,20]};
%          getoption('view',varargin)
%          ans =
%               20   20
%
%   NOTES: Input verification/ errorahandling should be added in the caller function to
%   keep things stable.
%
varargin = varargin{:};
ind1 = find(strcmpi(varargin,mode), 1);
if ~isempty(ind1)
    %Errorhandling    
    if ~iscell(varargin)
        varargin = {varargin};
    end
    if ind1+1 <= length(varargin)
        param = varargin{ind1+1};
    else
%         error(['No options are followed by the property ''', mode,''' '])
        param = [];
    end
else
    param = [];
end
end