%
% This script is messy, but really important. Don't lose it. It draws the
% AutoCAD which specifies the APCW structures.
%
function writeAutoCAD(varargin)
% Last modified: 2016-12-08
% writeAutoCAD() is a function based on AutoCAD_APCW(). It became apparent
% that a clean rewrite would make future comprehension and modification
% easier. At launch, the function asks for a .param file, which describes
% the structure to be written. A structure consists of the following
% elements:
% * coupler
% * nanobeam
% * node array
% * cooled nanobeam
% * photonic crystal
[s1,s2,s3,s4,s5,s6,s7] = importParameters(varargin);
fid = fopen(s1.filename,'w');
defineLISPfunctions(fid);
initializeDrawing(fid,s1);
s6 = compensateLength(s1,s2,s3,s4,s5,s6,s7);
for dev = 1:s1.num
    [~,~] = drawStructure(fid,dev,s1,s2,s3,s4,s5,s6,s7);
end
createAlignmentMarkers(fid,s1,round(mean(1:s1.num))*s1.spacing + 160/2);
if s1.pg
    createPerpendicularVeeGrooves(fid,s1,s2,s3); % draws alignment vees
    createFlipChip(fid,s1,s2,s3);
end
zoomExtents(fid);
hideLayer(fid,'grid');
% save file
fn = [pwd(),'/',strrep(s1.filename,'lsp','dxf')]; % save in current directory, as dxf.
saveDXF(fid,fn);
fclose(fid);
% alert Andrew that it's done
beep(); pause(0.2); beep();
end

function addDevMarkers(fid,xc,yc,vgw,dev)
% this function adds markers near the coupler identifying which device
% you're on. xc,yc is the location of the coupler
setCurrentLayer(fid,'window');
rad = 5; % 5 um radius
x0 = xc - vgw/2;
y0 = yc + 4 * 2.5*rad;
ctr = 0;
while ctr < dev
    if mod(ctr,4)==0 % typewriter carriage return
        y0 = y0 - 4 * 2.5 * rad;
        x0 = x0 - 2.5 * rad;
    end
    writeCircle(fid,x0,y0,rad);
    y0 = y0 + 2.5*rad;
    ctr = ctr + 1;
end
setCurrentLayer(fid,'nitride');
end

function addPHCMarkers(fid,xc,yc,wf,lf)
% this function adds markers along the PHC for SEM imaging; they will be
% circles representing in binary the distance
dsp = .800; % displace marker 800 nm from edge
rad = .150; % 100 nm radius;
x0 = xc - wf/2 - dsp - rad;
y0 = yc - lf/2;
mkrs = 10:2:lf;
setCurrentLayer(fid,'text');
for ix1 = 1:numel(mkrs)
    % find coordinates of current marker
    x_curr = x0;
    y_curr = y0 + mkrs(ix1);
    % convert number to binary
    str = dec2bin(mkrs(ix1));
    for ix2 = length(str):-1:1
        if str(ix2) == '1'
            writeCircle(fid,x_curr,y_curr,rad);
        else
            writeCircle(fid,x_curr,y_curr,rad/2);
        end
        x_curr = x_curr - 2.5*rad;
    end
end
setCurrentLayer(fid,'nitride');
end

function addPHCText(fid,xc,yc,wf,lf,s)
% function writes PHC text and adds it near photonic crystal. Only works in
% autocad for windows.
dspX = 10; % displace 10 um from waveguide
dspY = 0;  % no Y displacement
scl = 2;   % 2 um size letters
rot = -90; % rotate CW 90 degrees
x0 = xc - wf/2 - dspX - scl; % start position x
y0 = yc + lf/2 + dspY;       % start position y
str = ['a = ',num2str(s.a*1000),' nm'];
writeText(fid,x0,y0,scl,rot,str);
x0 = x0 - 1.5*scl;
str = ['gap = ',num2str(s.gap*1000),' nm'];
writeText(fid,x0,y0,scl,rot,str);
x0 = x0 - 1.5*scl;
str = ['amplitude = ',num2str(s.A*1000),' nm'];
writeText(fid,x0,y0,scl,rot,str);
x0 = x0 - 1.5*scl;
str = ['start width = ',num2str(s.wf*1000),' nm'];
writeText(fid,x0,y0,scl,rot,str);
x0 = x0 - 1.5*scl;
str = ['width = ',num2str(s.w*1000),' nm'];
writeText(fid,x0,y0,scl,rot,str);
x0 = x0 - 1.5*scl;
str = ['taper: ',num2str(s.nt)];
writeText(fid,x0,y0,scl,rot,str);
x0 = x0 - 1.5*scl;
str = ['nominal: ',num2str(s.n)];
writeText(fid,x0,y0,scl,rot,str);
end

function addScaleGrid(fid,x,y,n)
% this function adds markers with LRHS location (x, y) to act as a scale bar
rad = 0.075; % this radius seems to still be visible
ptch = 0.300; % fix pitch
dx = 0; dy = 0;
setCurrentLayer(fid,'text')
for ix = 1:n
    for iy = 1:n
        xc = x - dx;
        yc = y + dy;
        writeCircle(fid,xc,yc,rad);
        dy = dy + ptch;
    end
    dx = dx + ptch;
    dy = 0;
end
setCurrentLayer(fid,'nitride')
end

function coords = addSpring(s,coords)
% SET AMPLITUDE TO ZERO IF NOT EXPLICITLY DEFINED ALREADY
if ~isfield(s,'samp')
    s.samp = 0;
end
xc = s.veeGrooveWidth/4;
if s.samp == 0
    N =0;
else
    N = 100;
end
% ADDS COSINE MODULATION TO TETHERS
locs = diff(coords(:,1) < s.veeGrooveWidth/2);
p = [find(locs==-1),find(locs==1)]; clear locs;
% BOTTOM
x0 = coords(p(1),1); y0 = coords(p(1),2);
xf = coords(p(1)+1,1); yf = coords(p(1)+1,2);
x_bottom = xc + linspace(-5,5,N);
y_bottom = interp1([x0,xf],[y0,yf],x_bottom);
y_bottom = y_bottom + s.samp * (-1 + ...
    cos(4*pi*(x_bottom-xc)./(max(x_bottom)-min(x_bottom))));
% TOP
x0 = coords(p(2),1); y0 = coords(p(2),2);
xf = coords(p(2)+1,1); yf = coords(p(2)+1,2);
x_top = xc + linspace(5,-5,N);
y_top = interp1([x0,xf],[y0,yf],x_top);
y_top = y_top + s.samp * (-1 + ...
    cos(4*pi*(x_top-xc)./(max(x_top)-min(x_top))));
% CREATE OUTPUT
coords_out(:,1) = [coords(1:p(1),1);...
    x_bottom';...
    coords(p(1)+1:p(2),1);...
    x_top';...
    coords(p(2)+1:end,1)];
coords_out(:,2) = [coords(1:p(1),2);...
    y_bottom';...
    coords(p(1)+1:p(2),2);...
    y_top';...
    coords(p(2)+1:end,2)];
coords = coords_out;
end

function addString(fid,str)
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
writeOutput(fid,str);
end

function anonfun = amplitudeModulation(s)
% Takes struct input and returns anonymous function
n_seg = numel(s); % number of segments
edges = [0,cumsum([s.n])]; % bounds for piecewise function
anonfun = '@(y) ';
for j = 1:n_seg
    if j > 1
        anonfun = [anonfun,' + '];
    end
    if j < n_seg
        bln = ['((y>=',num2str(edges(j)),...
            ') & (y<',num2str(edges(j+1)),'))'];
    else
        bln = ['(y>=',num2str(edges(j)),')'];
    end
    x = s(j).ix + edges(j);
    m = repmat(x,[numel(x),1])'.^(repmat((1:numel(x))-1,[numel(x),1]));
    sol = s(j).amp;
    vec = m\sol';
    fn = '(';
    for k = 0:(numel(vec)-1)
        if k > 0
            fn = [fn,' + '];
        end
        fn = [fn,num2str(vec(k+1)),'.*y.^',num2str(k)];
    end
    fn = [fn,')'];
    anonfun = [anonfun,bln,'.*',fn];
end
anonfun = eval(anonfun);
end

function clearSelectionSets(fid,namedObjects)
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
for jj = 1:numel(namedObjects)
    currObj = namedObjects{jj};
    writeOutput(fid,['(setq ',currObj,' nil)'])
end
end

function str = colorLayer(layerColor,layerName)
str = ['(command "_layer" "color" "',layerColor,'" "',layerName,'" "")'];
end

function s = coerceStartDimensions(s,w0,bw0)
% check that beams stitch correctly to previous part
if abs(s.A(1).amp(1) - s.w(1).amp(1) - bw0) > 0.0005
    % Oops, done fucked something up. Need to correct w0.
    s.w(1).amp(1) = bw0;
end
w_match = 2*(s.A(1).amp + s.w(1).amp) + s.gap; % doesn't work if gap is changing
w_match = w_match(1);
if abs(w_match - w0) > 0.0005
    % Oops, done fucked something up. Need to correct amplitude.
    s.A(1).amp(1) = (w_match - s.gap)/2 - s.w(1).amp;
end
% force end values to stitch as well
s.A(end).amp(end) = s.A(1).amp(1);
s.w(end).amp(end) = s.w(1).amp(1);
end

function s6 = compensateLength(s1,s2,s3,s4,s5,s6,s7)
for dev = 1:s1.num
    [~,y0] = setOffsets(dev,s1.spacing,s2(dev).overhang);
    [~,yc] = drawStructure(0,dev,s1,s2,s3,s4,s5,s6,s7);
    % calculate length, add length to cooled nanobeam
    len = 2*(yc-y0);
    target_len = (s1.wl-2*s2(dev).overhang);
    dlen = target_len-len;
    s6(dev).length = s6(dev).length +dlen/2;
end
end

function s = convertToModulatingFunction(s)
% This function updates files with old formatting
old_s = s;
s.A = struct('n',[],'order',[],'ix',[],'amp',[]);
s.w = s.A;
% create amplitude mod function
s.A(1).n = old_s.nt; % taper
s.A(2).n = old_s.n; % nominal crystal
s.A(1).order = 1; % taper
s.A(2).order = 0; %nominal crystal
s.A(1).ix = [0,old_s.nt];
s.A(3) = s.A(1); % duplicate for taper down
s.A(1).amp = [0,old_s.A];
s.A(3).amp = fliplr(s.A(1).amp);
s.A(2).ix = 0;
s.A(2).amp = old_s.A;
% duplicate fields for w
s.w = s.A;
s.w(1).amp = ((old_s.w-old_s.wf)/old_s.A)*(s.w(1).amp) + old_s.wf;
s.w(2).amp = ((old_s.w-old_s.wf)/old_s.A)*(s.w(2).amp) + old_s.wf;
s.w(3).amp = ((old_s.w-old_s.wf)/old_s.A)*(s.w(3).amp) + old_s.wf;
end

function convertNamedObjectToBlock(fid,namedObject,ctr)
basePoint = [num2str(ctr(1)),',',num2str(ctr(2))];
addString(fid,['(command "_block" "blk',... % create a block
    namedObject,... % called blk(namedObject)
    '" "',...
    basePoint,... % with a coordinate system centered at this point
    '" ',...
    namedObject,... % from this selection set
    ' "")']);
end

function createAlignmentMarkers(fid,s1,x0)
mrsize = 20;
y0 = 3020;
extentX = s1.ww;
extentY = 5500;
% extent markers
setCurrentLayer(fid,'extent');
writeMarker(fid, x0 - extentX + mrsize/2, y0 - extentY + mrsize/2, mrsize, 1);
writeMarker(fid, x0 + extentX - mrsize/2, y0 - extentY + mrsize/2, mrsize, 1);
writeMarker(fid, x0 - extentX + mrsize/2, y0 + extentY - mrsize/2, mrsize, 1);
writeMarker(fid, x0 + extentX - mrsize/2, y0 + extentY - mrsize/2, mrsize, 1);
% edge markers
setCurrentLayer(fid,'marks');
wid = (s1.num+2)*s1.spacing;
height = 3520;
% yoff = -1304.9750;
yoff = -800;

writeMarker(fid, x0-wid/2, yoff, mrsize, 0);
writeMarker(fid, x0+wid/2, yoff, mrsize, 0);
writeMarker(fid, x0-wid/2, height+yoff, mrsize, 0);
writeMarker(fid, x0+wid/2, height+yoff, mrsize, 0);

setCurrentLayer(fid,'marks');
writeMarker(fid, x0-wid/2, yoff+480, mrsize, 1);
writeMarker(fid, x0+wid/2, yoff+480, mrsize, 1);
writeMarker(fid, x0-wid/2, height+yoff-480, mrsize, 1);
writeMarker(fid, x0+wid/2, height+yoff-480, mrsize, 1);

setCurrentLayer(fid,'nitride');
end

function [xf,yf,wf,tcp] = createCooledSingleNanobeam(fid,s6,x0,y0,w0)
% set defaults if not specified
s6 = setCooledSingleNanobeamDefaults(s6);
% calculate intermediate quantities
lrest = s6.length - s6.lstart - s6.lstop;
numNodes = round((lrest - s6.lnode)/s6.spacing) + 1;
% draw
if numNodes > 1
    % lstart and lstop fixed; otherwise spaced evenly
    s6.spacing = (lrest - s6.lnode)/(numNodes - 1); % finds spacing for int. num of nodes
    nanobeam.length = s6.lstart;
    nanobeam.w0 = s6.w0;
    nanobeam.num = 2; % start and end points
    [xf,yf,wf] = createSingleNanobeam(fid,nanobeam,x0,y0,w0); % draw first snb
    nanobeam.length = s6.spacing - s6.lnode;
    for m = 1:(numNodes - 1) % draw all but one of the nodes
        [xf,yf,wf,tcp(m)] = createNodeArray(fid,s6,xf,yf,wf);
        [xf,yf,wf] = createSingleNanobeam(fid,nanobeam,xf,yf,wf);
    end
    [xf,yf,wf,tcp(m+1)] = createNodeArray(fid,s6,xf,yf,wf); % draw last node
    nanobeam.length = s6.lstop;
    [xf,yf,wf] = createSingleNanobeam(fid,nanobeam,xf,yf,wf); % draw last snb
elseif numNodes == 1
    % one node array, in center
    nanobeam.length = (s6.length - s6.lnode)/2;
    nanobeam.w0 = s6.w0;
    nanobeam.num = 2;
    [xf,yf,wf] = createSingleNanobeam(fid,nanobeam,x0,y0,w0); % draw first snb
    [xf,yf,wf,tcp] = createNodeArray(fid,s6,xf,yf,wf); % draw lonely node
    [xf,yf,wf] = createSingleNanobeam(fid,nanobeam,xf,yf,wf); % draw second snb
elseif numNodes < 1
    % no cooling
    nanobeam.length = s6.length;
    nanobeam.w0 = s6.w0;
    nanobeam.num = 2;
    [xf,yf,wf] = createSingleNanobeam(fid,1,nanobeam,x0,y0,w0); % draw first snb
end
end

function [xf,yf,wf] = createCoupler(fid,s,x0,y0)
coords = drawTether(s);
coords = mirrorCoords(coords);
coords = displaceCoords(coords,x0,y0,0);
[xf,yf,wf] = writePLine(fid,coords);
addString(fid,'(setq currobj (ssadd (entlast) currobj))'); % Add this last
% object to list
if s.tccb % if traffic circle coupler
    filletLastObject(fid,s.tccr); % fillets the coupler
    setCurrentLayer(fid,'holes');
    writeCircle(fid,x0,y0+s.prelength+s.tetherCenterWidth/2,s.hole); % draws a hole
    addString(fid,'(setq currobj (ssadd (entlast) currobj))'); % Add this last
    setCurrentLayer(fid,'nitride');
end
xd = x0 - 0.8; % displace 800 nm from start point
yd = y0;
addScaleGrid(fid,xd,yd,3);
% DRAW PRE-BIAS RECTANGLES
prebiasCoupler(fid,coords);
end

function [xf,yf,wf,bwf] = createDoubleNanobeam(fid,s,x0,y0,w0,bw0)
% if (fid == 0)
%     writeOutput = @(fid,s) writeToStandardOut(fid,s);
% else
%     writeOutput = @(fid,s) writeToFile(fid,s);
% end
if s.newt % if new taper % kludge kludge kludge 2015/09/21
    % calculate initial gap
    g_c = (s.gap + s.w - s.A);
    s.gap = g_c - s.wf;
end
xf = x0;
yf = y0+s.length;
bwf = bw0;
wf = w0;
% writeOutput(fid,'; DRAWING DOUBLE NANOBEAM');
writeRectangle(fid,x0,y0,xf-bw0,yf);
addString(fid,'(setq currobj (ssadd (entlast) currobj))');
writeRectangle(fid,x0-w0,y0,xf-w0+bw0,yf);
addString(fid,'(setq currobj (ssadd (entlast) currobj))');
end

function [xc,yc,namedObjects] = createDoubleNanobeamCrystal(fid,dev,s,xf,yf,wf,namedObjects)
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
addString(fid,'(setq currobj nil)');        % Free currobj
addString(fid,'(setq currobj (ssadd))');    % I FORGET WHY THIS IS HERE
writeOutput(fid,'; DRAWING PHOTONIC CRYSTAL');
%% CREATE Y JUNCTION
if s.yneg% if inverted tone y
    setCurrentLayer(fid,'negYjunc');
    [xf,yf,wf,bwf] = createNegYJunction(fid,s,xf,yf,wf);
    setCurrentLayer(fid,'nitride');
else
    [xf,yf,wf,bwf] = createYJunction(fid,s,xf,yf,wf);
    addString(fid,'(setq currobj (ssadd (entlast) currobj))');
end
addString(fid,['(setq yjunction',num2str(dev),' currobj)']);
namedObjects{end+1} = ['yjunction',num2str(dev)];
%% CREATE BLANK DOUBLE NANOBEAM -- Removed 2016-05-20
if false
    addString(fid,'(setq currobj nil)');
    addString(fid,'(setq currobj (ssadd))');
    [xf,yf,wf,bwf] = createDoubleNanobeam(fid,s,xf,yf,wf,bwf);
    addString(fid,['(setq doublenanobeam',num2str(dev),' currobj)']);
    namedObjects{end+1} = ['doublenanobeam',num2str(dev)];
else
    s.length = 0;
end
%% CREATE DOUBLE NANOBEAM CRYSTAL
addString(fid,'(setq currobj nil)');
addString(fid,'(setq currobj (ssadd))');
[xc,yc,wf,lf] = createDPHC(fid,s,xf,yf,wf,bwf,dev);
addString(fid,['(setq phc',num2str(dev),' currobj)']); % not a mirrored object.
%% BOX OFF PHC REGION (fracturing cover)
setCurrentLayer(fid,'cover');
scl = 1; % how much bigger does the box need to be than minimally?
% writeRectangle(fid,xc-scl*wf/2,yc-scl*lf/2,xc+scl*wf/2,yc+scl*lf/2);
writeRectangle(fid,xc-scl*wf/2,yc-scl*lf/2,xc+scl*wf/2,yc+scl*lf/2);
setCurrentLayer(fid,'nitride');
%% ADD PERIOD MARKERS
addPHCMarkers(fid,xc,yc,wf,lf);
%% ADD PHC TEXTUAL INFO
addPHCText(fid,xc,yc,wf,lf,s);
end

function [xc,yc,wf,lf] = createDPHC(fid,s,x0,y0,w0,bw0,dev)
% ACM 2015-02-10
% This function is completely rewritten. The idea is to express the
% photonic crystal as a total number of lattice constants, and envelope
% modulating functions. Currently, the segments are parameterized as
% a struct with the following fields:
% n = number of cells in segment
% order = polynomial order
% ix = a list of cells...
% amp = ...and their corresponding amplitudes

% check to see if input has the old format (ramp, nominal, ramp)
if ~isstruct(s.A)
    s = convertToModulatingFunction(s);
end

% make sure important dimensions match
if ~s.newt % doesn't fix things if using the new taper (newt)
    s = coerceStartDimensions(s,w0,bw0);
end

% get anonymous functions
amp_fun = amplitudeModulation(s.A);
wid_fun = amplitudeModulation(s.w);

% if new taper
if s.newt
    % g_nominal + w_nominal - a_nominal is fixed
    ctr_cell = (s.n + 2 * s.nt)/2; % center of structure; this may not work if arb. modulation
    g_c = s.gap + wid_fun(ctr_cell) - amp_fun(ctr_cell);
    gap_fun = @(x) g_c - wid_fun(x) + amp_fun(x);
else
    gap_fun = @(x) s.gap + 0.*x;
end
% CREATE RETURN VALUES FOR FUNCTION
ctr_coords = displaceCoords([0,0],x0,y0,w0) + [0,sum([s.A(:).n])*s.a/2];
xc = ctr_coords(1); yc = ctr_coords(2);
% wf = max(coordsRHS(:,1) - coordsLHS(:,1)); % maximum width
ext_fun = @(x) gap_fun(x)/2 + wid_fun(x) + amp_fun(x).*cos(2*pi*x);
max_excur = fminbnd(@(x) -ext_fun(x),0,2*s.nt + s.n);
wf = 2 * ext_fun(max_excur);
lf = 2*(s.length + s.ltaper) + s.a*sum([s.A(:).n]); % length of Y + dnb + phc

% draw, mirror, displace, PLine, name
if s.neg % fork for cell-based fracturing
    %%% THIS IS EXPERIMENTAL; THE PLAN IS TO USE THE _BLOCK AND _INSERT
    %%% COMMANDS TO TAKE ADVANTAGE OF THE HIERARCHICAL PROCESSING IN
    %%% Layout BEAMER. IDEALLY THIS WAY EACH CELL OF THE PHOTONIC CRYSTAL
    %%% WILL BE FRACTURED IDENTICALLY.
    % sym = checkTaperSymmetry(s); % true if x-symmetric
    % s.w and s.A structs reveal where the nominal crystal is.
    % nominal = find([s.A(:).order]==0); % the segment which is specified
    % to be zero-order is flat. There should just be one of these, because
    % the default behavior is to use first-order polynomials to stitch if
    % the input function is complicated.
    % PHCBlockCtr = 0; % counter for the blocks in the photonic crystal
    %%% DRAW TAPER IN %%%
    rng = s.A(1).ix;
    if s.segTaper % MODIFIED 2016-05-19 to encourage better fracturing
        addString(fid,['(setq tap',num2str(dev),' nil)']); % clear variable
        addString(fid,['(setq tap',num2str(dev),' (ssadd))']); % Not entirely clear why this is necessary
        setCurrentLayer(fid,'phc');
        % THIS FUNCTION DOES THE OPERATION BELOW, BUT SEGMENTS BY LATTICE CONSTANT
        wf = drawSegNegTaper(fid,dev,s,rng,amp_fun,wid_fun,gap_fun,x0,y0,w0,xc,yc);
    else
        [coords_int,coords_ext] = drawNegTaper(s,rng,amp_fun,wid_fun,gap_fun);
        coords_int = displaceCoords(coords_int,x0,y0,w0);
        coords_ext = displaceCoords(coords_ext,x0,y0,w0);
        wf = 2 * (max(coords_ext(:,1))-xc); % final width
        setCurrentLayer(fid,'phc');
        [~,~,~] = writePLine(fid,coords_int); % DRAW INTERIOR BOX
        addString(fid,['(setq tap',num2str(dev),' nil)']); % clear variable
        addString(fid,['(setq tap',num2str(dev),' (ssadd))']); % Not entirely clear why this is necessary
        addString(fid,['(setq tap',num2str(dev),' (ssadd (entlast) tap',num2str(dev),'))']); % create tap1 selection set
        mirrorNamedObject(fid,'(entlast)',[xc,yc],[xc,yc+1]); % reflect
        addString(fid,['(setq tap',num2str(dev),' (ssadd (entlast) tap',num2str(dev),'))']); % add entity to selection set
        [~,~,~] = writePLine(fid,coords_ext);
        addString(fid,['(setq tap',num2str(dev),' (ssadd (entlast) tap',num2str(dev),'))']); % add entity to selection set
        mirrorNamedObject(fid,'(entlast)',[xc,yc],[xc,yc+1]); % reflect
        addString(fid,['(setq tap',num2str(dev),' (ssadd (entlast) tap',num2str(dev),'))']); % add entity to selection set
    end
    convertNamedObjectToBlock(fid,['tap',num2str(dev)],[xc,yc]); % reference point is center of crystal
    insertBlock(fid,['blktap',num2str(dev)],[xc,yc],1,1,0);
    %%% DRAW TAPER OUT %%%
    % ASSUMES SYMMETRY ABOUT [xc,yc]
    insertBlock(fid,['blktap',num2str(dev)],[xc,yc],1,1,180);
    %%% DRAW NOMINAL CRYSTAL %%%
    [coords_int,coords_ext] = drawNegUnitCell(s,amp_fun,wid_fun,gap_fun);
    addString(fid,['(setq uc',num2str(dev),' nil)']); % clear variable
    addString(fid,['(setq uc',num2str(dev),' (ssadd))']); % Not entirely clear why this is necessary
    [~,~,~] = writePLine(fid,coords_int);
    addString(fid,['(setq uc',num2str(dev),' (ssadd (entlast) uc',num2str(dev),'))']); % add entity to selection set
    [~,~,~] = writePLine(fid,coords_ext);
    addString(fid,['(setq uc',num2str(dev),' (ssadd (entlast) uc',num2str(dev),'))']); % add entity to selection set
    convertNamedObjectToBlock(fid,['uc',num2str(dev)],[0,s.a/2]); % reference point is center of cell
    for jj = 1:s.n
        % for whatever reason, the insertion point is at the corner of
        % the unit cell. I don't understand why this happens.
        insertBlock(fid,['blkuc',num2str(dev)],displaceCoords([0,(s.nt+jj)*s.a],x0,y0,w0),1,1,0);
        insertBlock(fid,['blkuc',num2str(dev)],displaceCoords([0,(s.nt+jj-1)*s.a],x0,y0,w0),1,1,180);
    end
    %%% DRAW PHC COVER %%%
    setCurrentLayer(fid,'phccover');
    xtall = s.a * (2*s.nt + s.n);
    writeRectangle(fid,xc-wf/2,yc+xtall/2,xc+wf/2,yc-xtall/2)
    setCurrentLayer(fid,'nitride')
else % traditional fracture
    coordsRHS = drawDPHC(s,amp_fun,wid_fun,gap_fun);
    coordsLHS = [-coordsRHS(:,1),coordsRHS(:,2)];
    coordsRHS = displaceCoords(coordsRHS,x0,y0,w0);
    coordsLHS = displaceCoords(coordsLHS,x0,y0,w0);
    [~,~,~] = writePLine(fid,coordsRHS);
    addString(fid,'(setq currobj (ssadd (entlast) currobj))');
    % [~,~,~] = writePLine(fid,coordsLHS);
    % addString(fid,'(setq currobj (ssadd (entlast) currobj))');
    % TEST! mirror to create LHS
    mirrorNamedObject(fid,'(entlast)',[xc,yc],[xc,yc+1]);
    addString(fid,'(setq currobj (ssadd (entlast) currobj))');
end
end

function createFlipChip(fid,s1,s2,s3)
edge = 25; gap = s1.tg; % tiling gap
for jj = 1:numel(s2)
    % BOTTOM SIDE
    grooveLength = 2e3;
    couplerDisp = 5;
    holeLatticeConstant = 5;
    [x0,y0] = setOffsets(jj,s1.spacing,s2(jj).overhang); % get dev1 coupler position
    x0 = x0 + s3(1).nominalWidth/2 - s1.fcgw/2; % move to edge of flip chip groove
    xf = x0 + s1.fcgw;
    yf = y0 - grooveLength; % 2 mm length
    y0 = y0 - couplerDisp; % 5 micron recess from coupler
    makeGrooveAtPosition(x0,y0,xf,yf)
    % TOP SIDE
    y0 = y0 + couplerDisp - 2*s2(1).overhang + s1.wl;
    yf = y0 + grooveLength;
    y0 = y0 + couplerDisp;
    makeGrooveAtPosition(x0,y0,xf,yf)
end
makeWindow();
    function makeGrooveAtPosition(x0,y0,xf,yf)
        setCurrentLayer(fid,'flipChip'); % we're going to be writing on the flipChip layer
        % MAKE GROOVE
        writeRectangle(fid,x0,y0,xf,yf);
        setCurrentLayer(fid,'flipChipSubtract');
        createSliceGrid(fid,x0,y0,xf,yf,edge,gap);
        % MAKE NET OF HOLES
        setCurrentLayer(fid,'flipChip');
        rc = 0.3 * holeLatticeConstant;
        if y0 < 0
            yc = y0 + holeLatticeConstant/2;
            while yc < y0 + s1.fcgw/2; % this only works for one side, but for testing
                xc = x0 + holeLatticeConstant/2;
                while xc < xf
                    writeCircle(fid,xc,yc,rc);
                    xc = xc + holeLatticeConstant;
                end
                yc = yc + holeLatticeConstant;
            end
        else
            yc = y0 - holeLatticeConstant/2;
            while yc > y0 - s1.fcgw/2 % this only works for one side, but for testing
                xc = x0 + holeLatticeConstant/2;
                while xc < xf
                    writeCircle(fid,xc,yc,rc);
                    xc = xc + holeLatticeConstant; % 5 um spacing
                end
                yc = yc - 5;
            end
        end
    end
    function makeWindow()
        [x0,y0] = setOffsets(1,s1.spacing,s2(1).overhang);
        outside = s1.ww - (s2(1).vgw + 2 * s2(1).srw) * s1.num;
        outside = outside/2;
        radius = s1.wfr;
        xa = x0 - s2(1).vgw/2 - s2(1).srw - s2(1).vgb - outside;
        xb = xa + 2*outside + s2(1).vgw + s1.spacing * (s1.num-1);
        setCurrentLayer(fid,'flipChipWindow'); % we're going to be writing on the flipChip layer; should probably change to some generic back window layer
        writeRectangle(fid,xa,0,xb,s1.wl);
        filletLastObject(fid,radius); % fillets the coupler
        setCurrentLayer(fid,'flipChipSubtract');
        createSliceGrid(fid,xa,0,xb,s1.wl,edge,gap);
    end
end

function createIndicators(fid,s,xc,yc)
% draws markers which show where the crystal is
setCurrentLayer(fid,'indicators');
pitch = s.ha;
rad = s.hr;
d = 20;
t = s.srw*sqrt(2);
x1 = xc-s.vgw/2-s.srw;
x2 = xc+s.vgw/2+s.srw;
coords = [0,2*d;...
    0,d;...
    -d,0;...
    0,-d;...
    0,-2*d;...
    0,0;...
    0,-d + t;...
    -d+t,0;...
    0,d-t;...
    0,0;...
    ];
coords(:,1) = coords(:,1) + x1;
coords(:,2) = coords(:,2) + yc;
writePLine(fid,coords);
filletLastObject(fid,1);
coords(:,1) = coords(:,1) - x1;
coords(:,1) = -coords(:,1);
coords(:,1) = coords(:,1) + x2;
writePLine(fid,coords);
filletLastObject(fid,1);
% put holes in it so it doesn't collect Si
setCurrentLayer(fid,'holes');
hole_path = [0,d-t/2;-d+t/2,0;0,-d+t/2];
len = sqrt(sum((hole_path(1,:)-hole_path(2,:)).^2));
hole_coords = [linspace(hole_path(1,1),hole_path(2,1),round(len/pitch));...
    linspace(hole_path(1,2),hole_path(2,2),round(len/pitch))]';
hole_coords(1,:) = [];
hole_coords(end,:) = [];
hole_coords = vertcat(hole_coords,...
    [linspace(hole_path(2,1),hole_path(3,1),round(len/pitch));...
    linspace(hole_path(2,2),hole_path(3,2),round(len/pitch))]');
hole_coords(end,:) = [];
for jj = 1:numel(hole_coords(:,1))
    x0 = x1+hole_coords(jj,1);
    y0 = yc+hole_coords(jj,2);
    writeCircle(fid,x0,y0,rad);
end
hole_coords(:,1) = -hole_coords(:,1);
for jj = 1:numel(hole_coords(:,1))
    x0 = x2+hole_coords(jj,1);
    y0 = yc+hole_coords(jj,2);
    writeCircle(fid,x0,y0,rad);
end
setCurrentLayer(fid,'nitride');
end

function str = createLayer(layerName)
str = ['(command "_layer" "new" "',layerName,'" "")'];
str = sprintf(str);
end

function [xf,yf,wf,tcp] = createNodeArray(fid,s,x0,y0,w0)
% EDITED 2015-07-26: IF IN WINDOW, POSSIBLY BULGED, AND NEED TO ADD
s.w0 = w0; % assert start width
if y0 > 0 % if in window
    s.veeGrooveWidth = s.veeGrooveWidth + s.vgb * 2;
end
s.tetherYDisp = s.tetherYDisp -...
    (s.numTethers - 1) * (s.aout - s.a)/2;
nodeArrayCover = nodeArrayToCover(s);
% draw cover
coords = drawNodeArray(nodeArrayCover);
% EDITED 2015-07-15 TO ADD SPRING
coords = addSpring(s,coords);
coords = mirrorCoords(coords);
coords = displaceCoords(coords,x0,y0,w0);
tcp = mean(coords(coords(:,1)==max(coords(:,1)),2)); % tether center position
[xf,yf,wf] = writePLine(fid,coords);
addString(fid,'(setq currobj (ssadd (entlast) currobj))');
prebiasNodeArray(fid,coords);
% draw node array
setCurrentLayer(fid,'fda1');
coords = drawNegativeTethers(s);
coords = mirrorCoords(coords);
coords = displaceCoords(coords,x0,y0+s.blankLength + s.taperLength,w0);
% [xf,yf,wf] = writePLine(fid,coords);
% j = 1;
% while j < numel(coords)/2 % need to do this to close boxes for subtractive
%     % layer; relies on the idea that each element
%     % has 4 points
%     [~,~,~] = writePLine(fid,coords(j:j+3,:));
%     addString(fid,'(setq currobj (ssadd (entlast) currobj))');
%     j  = j + 4;
% end
% NOW EACH HAS 204 POINTS...
[nr,~] = size(coords);
nHoles = s.numTethers - 1;
plineSeg = (nr/nHoles)/2;
ctr = 1;
while ctr <= nr
    start = ctr;
    stop = ctr + plineSeg - 1;
    [~,~,~] = writePLine(fid,coords(start:stop,:));
    addString(fid,'(setq currobj (ssadd (entlast) currobj))');
    ctr = ctr + plineSeg;
end
xd = x0 - 1.2; % a little farther than for
yd = y0;
setCurrentLayer(fid,'nitride');
addScaleGrid(fid,xd,yd,5);
end

function createPerpendicularVeeGrooves(fid,s1,s2,s3)
[x0,y0] = setOffsets(1,s1.spacing,s2(1).overhang); % get dev1 coupler position
x0 = x0 + s3(1).nominalWidth/2 - s2(1).vgw/2; % move to edge of vee groove
x0 = x0 - 2*s1.spacing;
xf = x0-2e3; % 2mm long grooves
yf = y0 - s2(1).vgw;
edge = 25; gap = s1.tg; % tiling gap
% FIRST GROOVE
makeGrooveAtPosition(x0,y0,xf,yf);
makeSmallerGrooveAtPosition(x0,y0-s2(1).vgw/2+s1.fcgw/2,xf,y0-s2(1).vgw/2-s1.fcgw/2);
% SECOND GROOVE
y0 = y0 - 2*s2(1).overhang + s1.wl;
yf = yf - 2*s2(1).overhang + s1.wl + 2*s2(1).vgw;
makeGrooveAtPosition(x0,y0,xf,yf);
makeSmallerGrooveAtPosition(x0,y0+s2(1).vgw/2-s1.fcgw/2,xf,y0+s2(1).vgw/2+s1.fcgw/2);
% THIRD GROOVE
[x0,y0] = setOffsets(1,s1.spacing,s2(1).overhang); % get dev1 coupler position
x0 = x0 + (s1.num-1)*s1.spacing + s3(1).nominalWidth/2 + s2(1).vgw/2; % go to dev16 coupler position, then to edge of groove
x0 = x0 + 2 * s1.spacing;
xf = x0 + 2e3;
yf = y0 - s2(1).vgw;
makeGrooveAtPosition(x0,y0,xf,yf);
makeSmallerGrooveAtPosition(x0,y0-s2(1).vgw/2+s1.fcgw/2,xf,y0-s2(1).vgw/2-s1.fcgw/2);
% FOURTH GROOVE
y0 = y0 - 2*s2(1).overhang + s1.wl;
yf = yf - 2*s2(1).overhang + s1.wl + 2*s2(1).vgw;
makeGrooveAtPosition(x0,y0,xf,yf);
makeSmallerGrooveAtPosition(x0,y0+s2(1).vgw/2-s1.fcgw/2,xf,y0+s2(1).vgw/2+s1.fcgw/2);makeSmallerGrooveAtPosition(x0,y0-s2(1).vgw/2+s1.fcgw/2,xf,y0-s2(1).vgw/2-s1.fcgw/2);
    function makeGrooveAtPosition(x0,y0,xf,yf)
        setCurrentLayer(fid,'window'); % we're going to be writing on the window layer
        % MAKE GROOVE
        writeRectangle(fid,x0,y0,xf,yf)
        % ADD TILES
        setCurrentLayer(fid,'grid');
        createSliceGrid(fid,x0,y0,xf,yf,edge,gap);
        % RETURN TO DEFAULT LAYER
        setCurrentLayer(fid,'nitride');
    end
    function makeSmallerGrooveAtPosition(x0,y0,xf,yf)
        setCurrentLayer(fid,'flipChip'); % we're going to be writing on the flipChip layer
        % MAKE GROOVE
        writeRectangle(fid,x0,y0,xf,yf);
        setCurrentLayer(fid,'flipChipSubtract');
        createSliceGrid(fid,x0,y0,xf,yf,edge,gap);
    end
end

function [xf,yf,wf] = createSingleNanobeam(fid,s,x0,y0,w0)
% if (fid == 0)
%     writeOutput = @(fid,s) writeToStandardOut(fid,s);
% else
%     writeOutput = @(fid,s) writeToFile(fid,s);
% end
s.w0 = w0; % assert start width

coords = drawSingleNanobeam(s);
coords = mirrorCoords(coords);
coords = displaceCoords(coords,x0,y0,w0);
[xf,yf,wf] = writePLine(fid,coords);
addString(fid,'(setq currobj (ssadd (entlast) currobj))'); % Add this last
prebiasWaveguide(fid,coords);
% object to list
end

function createRailHoles(fid,s1,s2,x0)
% define hole spacing, number, etc.
n = floor((s1.wl - 2 * s2.he)/s2.ha); % length of window minus twice the
% recess divided by lattice constant
s2.ha = (s1.wl - 2 * s2.he)/n; % correct lattice constant
n = n - 1;
% draw first hole
x0 = x0 + s2.vgw/2 + s2.vgb + s2.srw/(s2.nhr+1);
y0 = s2.ha/2 + s2.he;
dx = s2.srw/(s2.nhr+1);
dy = s2.ha;
nx = s2.nhr;
ny = n;
r = s2.hr;
setCurrentLayer(fid,'holes');
writeCircle(fid,x0,y0,r);
% this array function is used for tiling too; consider making subfunction?
addString(fid,['(command "_array" "last" "" "R" "',num2str(ny),...
    '" "',num2str(nx),'" "',num2str(dy),'" "',num2str(dx),'")']);
x0 = x0 - s2.vgw - 2 * s2.vgb - 2*s2.srw/(s2.nhr+1);
dx = -s2.srw/(s2.nhr+1);
writeCircle(fid,x0,y0,r);
addString(fid,['(command "_array" "last" "" "R" "',num2str(ny),...
    '" "',num2str(nx),'" "',num2str(dy),'" "',num2str(dx),'")']);
setCurrentLayer(fid,'nitride');
end

function createRailHolesWithSprings(fid,s1,s2,x0)
% ADDED 2015-07-26 to accomodate
% define hole spacing, number, etc.
n = floor(s1.wl/s2.ha);
s2.ha = s1.wl/n; % correct lattice constant
n = n - 1;
% draw first hole
x0 = x0 + s2.vgw/2 + s2.srw/(s2.nhr+1);
y0 = s2.ha/2;
dx = s2.srw/(s2.nhr+1);
dy = s2.ha;
% nx = s2.nhr;
ny = n;
r = s2.hr;
setCurrentLayer(fid,'holes');
y = y0:dy:(y0+ny*dy);
f = @(x) (((x>10)&(x<(10 + 2*s2.sp)))...
    |((x>s1.wl-10-2*s2.sp)&(x<(s1.wl-10)))) ...
    .* (-1 + cos(4*pi*(x-10)./40)) .* s2.sa;
x = f(y);
for kk = 1:numel(x)
    writeCircle(fid,x(kk)+x0,y(kk),r);
    writeCircle(fid,x(kk)+x0+dx,y(kk),r);
end
x0 = x0 - s2.vgw - 2*s2.srw/(s2.nhr+1);
dx = -s2.srw/(s2.nhr+1);zoom

for kk = 1:numel(x)
    writeCircle(fid,-x(kk)+x0,y(kk),r);
    writeCircle(fid,-x(kk)+x0+dx,y(kk),r);
end
setCurrentLayer(fid,'nitride');
end

function createSliceGrid(fid,x0,y0,xf,yf,edge,gap)
x = sort([x0,xf]); x0 = x(1); xf = x(2);
y = sort([y0,yf]); y0 = y(1); yf = y(2);
nx = ceil((xf - x0 - gap)/(edge + gap));
ny = ceil((yf - y0 - gap)/(edge + gap));
dx = (xf - x0 - gap)/nx;
dy = (yf - y0 - gap)/ny;
writeRectangle(fid,x0+gap,y0+gap,x0+dx,y0+dy);
addString(fid,['(command "_array" "last" "" "R" "',num2str(ny),...
    '" "',num2str(nx),'" "',num2str(dy),'" "',num2str(dx),'")']);
end

function createWindow(fid,s1,s2,tcp,dev,x0)
% this compensates the vee groove bulge
setCurrentLayer(fid,'window');
outside = s1.ww - (s2(dev).vgw + 2 * s2(dev).srw) * s1.num;
outside = outside/2;
edge = 25;
gap = s1.tg;
radius = s1.wfr;
if dev == 1 % draw void to left
    xa = x0 - s2(dev).vgw/2 - s2(dev).srw - s2(dev).vgb;
    xb = xa - outside;
    setCurrentLayer(fid,'window');
    writeRectangle(fid,xa,0,xb,s1.wl);
    filletLastObject(fid,radius);
    setCurrentLayer(fid,'grid');
    createSliceGrid(fid,xa,0,xb,s1.wl,edge,gap)
end
if dev == s1.num % draw void to right
    xa = x0 + s2(dev).vgw/2 + s2(dev).srw + s2(dev).vgb;
    xb = xa + outside;
    setCurrentLayer(fid,'window');
    writeRectangle(fid,xa,0,xb,s1.wl);
    filletLastObject(fid,radius);
    setCurrentLayer(fid,'grid');
    createSliceGrid(fid,xa,0,xb,s1.wl,edge,gap)
    return; % nothing left to draw
end
setCurrentLayer(fid,'window');
writeRectangle(fid,x0+s2(dev).vgw/2+s2(dev).srw +s2(dev).vgb,0,...
    x0+s1.spacing-s2(dev+1).vgw/2-s2(dev+1).srw-s2(dev+1).vgb,s1.wl);
filletLastObject(fid,radius);
%% TETHERS
% REMOVED 2016/04/07. NOT BEING USED
% % The tethers are new (2015/12/08). Hope to solve yield problems. Currently
% % implemented to draw in the window to the right, but should probably split
% % left/right (in case tethers are in different places in each).
% setCurrentLayer(fid,'tethers');
% for jj = 1:numel(tcp)
%     xx1 = x0+s2(dev).vgw/2+s2(dev).srw + s2(dev).vgb;
%     dx = s1.spacing... % device to device length minus
%         - ((s2(dev).vgw...
%         + s2(dev+1).vgw)/2 ... % half of each vee groove
%         + s2(dev).srw...
%         + s2(dev+1).srw... % each siderail
%         + s2(dev).vgb...
%         + s2(dev+1).vgb);   % each siderail bulge
%     yy = tcp(jj);
%     rad = 1; % one micron fillet radius
%     thk = 1;
%     coords = [xx1,yy + rad + thk/2;
%         xx1, yy + thk/2;
%         xx1 + dx, yy + thk/2;
%         xx1 + dx, yy + rad + thk/2;
%         xx1 + dx, yy - rad - thk/2;
%         xx1 + dx, yy - rad - thk/2;
%         xx1 + dx, yy - thk/2;
%         xx1, yy - thk/2;
%         xx1, yy - rad - thk/2
%         ];
%     writePLine(fid,coords);
%     filletLastObject(fid,rad);
%     % writeRectangle(fid,xx1,yy-1,xx1+dx,yy+1);
% end
%% TILES
setCurrentLayer(fid,'grid');
createSliceGrid(fid,x0+s2(dev).vgw/2+s2(dev).srw+s2(dev).vgb,0,...
    x0+s1.spacing-s2(dev+1).vgw/2-s2(dev+1).srw-s2(dev+1).vgb,s1.wl,edge,gap);
setCurrentLayer(fid,'nitride');
end

function createWindowWithStruts(fid,s1,s2,dev,x0)
setCurrentLayer(fid,'window');
outside = s1.ww - (s2.vgw + 2 * s2.srw) * s1.num;
outside = outside/2;
edge = 25;
gap = s1.tg;
sThk = 10; % support rail thickness
triBase = 65; % right isoceles triangle side length
bump = 5; % size of bump near photonic crystal
% radius = 1;
radius = s1.wfr;
if dev == 1 % draw void to left
    xa = x0 - s2.vgw/2 - s2.srw;
    xb = xa - outside;
    setCurrentLayer(fid,'window');
    % writeRectangle(fid,xa,0,xb,s1.wl);
    % REPLACED 2015/07/14 to add minimalistic spiderweb support
    % DRAW OUTER WINDOW
    clear coords;
    coords(:,1) = [xb,xa-sThk-triBase,xa,...
        xa,xa-bump,xa,...
        xa,xa-sThk-triBase,xb];
    coords(:,2) = [0,0,triBase+sThk,...
        s1.wl/2 + [-bump,0,bump],...
        s1.wl-triBase-sThk,s1.wl,s1.wl];
    writePLine(fid,coords);
    filletLastObject(fid,radius);
    % DRAW TRIANGLE IN CORNERS
    clear coords;
    coords(:,1) = [xa,xa-triBase,xa];
    coords(:,2) = [0,0,triBase];
    writePLine(fid,coords);
    filletLastObject(fid,radius);
    clear coords;
    coords(:,1) = [xa,xa-triBase,xa];
    coords(:,2) = s1.wl-[0,0,triBase];
    writePLine(fid,coords);
    filletLastObject(fid,radius);
    setCurrentLayer(fid,'grid');
    createSliceGrid(fid,xa,0,xb,s1.wl,edge,gap)
end
if dev == s1.num % draw void to right
    xa = x0 + s2.vgw/2 + s2.srw;
    xb = xa + outside;
    setCurrentLayer(fid,'window');
    % writeRectangle(fid,xa,0,xb,s1.wl);
    % REPLACED 2015/07/14 to add minimalistic spiderweb support
    % DRAW OUTER WINDOW
    clear coords;
    coords(:,1) = [xb,xa+sThk+triBase,xa,...
        xa,xa+bump,xa,...
        xa,xa+sThk+triBase,xb];
    coords(:,2) = [0,0,triBase+sThk,...
        s1.wl/2 + [-bump,0,bump],...
        s1.wl-triBase-sThk,s1.wl,s1.wl];
    writePLine(fid,coords);
    filletLastObject(fid,radius);
    % DRAW TRIANGLE IN CORNERS
    clear coords;
    coords(:,1) = [xa,xa+triBase,xa];
    coords(:,2) = [0,0,triBase];
    writePLine(fid,coords);
    filletLastObject(fid,radius);
    clear coords;
    coords(:,1) = [xa,xa+triBase,xa];
    coords(:,2) = s1.wl-[0,0,triBase];
    writePLine(fid,coords);
    filletLastObject(fid,radius);
    setCurrentLayer(fid,'grid');
    createSliceGrid(fid,xa,0,xb,s1.wl,edge,gap)
    return; % nothing left to draw
end

setCurrentLayer(fid,'window');
% writeRectangle(fid,x0+s2.vgw/2+s2.srw,0,x0+s1.spacing-s2.vgw/2-s2.srw,s1.wl);
% REPLACED 2015/07/14 to add minimalistic spiderweb support
% DRAW OUTER WINDOW
xa = x0+s2.vgw/2+s2.srw;
xb = x0+s1.spacing-s2.vgw/2-s2.srw;
clear coords;
coords(:,1) = [[xa,xa+sThk+triBase,xb-sThk-triBase,xb],...
    [xb,xb-bump,xb],...
    fliplr([xa,xa+sThk+triBase,xb-sThk-triBase,xb]),...
    [xa,xa+bump,xa]
    ];
coords(:,2) = [[sThk+triBase,0,0,sThk+triBase],...
    [-bump,0,bump]+s1.wl/2,...
    s1.wl-[sThk+triBase,0,0,sThk+triBase],...
    [bump,0,-bump]+s1.wl/2
    ];
writePLine(fid,coords);
filletLastObject(fid,radius);
% DRAW TRIANGLES
clear coords;
coords(:,1) = [xa,xa+triBase,xa]; coords(:,2) = [0,0,triBase];
writePLine(fid,coords); filletLastObject(fid,radius);
coords(:,2) = s1.wl - coords(:,2);
writePLine(fid,coords); filletLastObject(fid,radius);
coords(:,1) = [xb,xb-triBase,xb];
writePLine(fid,coords); filletLastObject(fid,radius);
coords(:,2) = s1.wl - coords(:,2);
writePLine(fid,coords); filletLastObject(fid,radius);
setCurrentLayer(fid,'grid');
createSliceGrid(fid,x0+s2.vgw/2+s2.srw,0,x0+s1.spacing-s2.vgw/2-s2.srw,s1.wl,edge,gap);
setCurrentLayer(fid,'nitride');
end

function createVeeGrooveWithSprings(fid,s1,s2,x0)
% EDIT 7/26/2015 TO ADD SPRINGS IN SIDERAILS
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
writeOutput(fid,'; DRAWING SPRINGY VEE GROOVE')
x0 = x0 - s2.vgw/2; xf = x0 + s2.vgw;
overlap = 5;
y0 = -(s2.fel(1) + overlap); yf = (s2.fel(2) + overlap) + s1.wl;
s2 = setVeeGrooveDefaults(s2);
setCurrentLayer(fid,'vgroove');
% THIS IS THE PART THAT NEEDS TO CHANGE; THE FDA, ETC CAN BE THE SAME AS
% THE ABOVE FUNCTION.
% AT THE EDGE OF THE WINDOW IS WHERE WE WANT THE SPRING.
yspring = linspace(0,2 * s2.sp); %
xspring = -s2.sa * (cos(4 * pi * yspring./max(yspring)) - 1);
x = [x0,... lower left
    x0 + xspring,... lower spring
    x0 + xspring,... upper spring
    x0,... upper left
    xf,... upper right
    s2.vgw + x0 - xspring,... upper spring
    s2.vgw + x0 - xspring,... lower spring
    xf]; % lower right
y = [y0,... lower left
    10 + yspring,... lower spring
    s1.wl - 10 - fliplr(yspring),... upper spring
    yf,...  upper left
    yf,... upper right
    fliplr(s1.wl - 10 - fliplr(yspring)),... upper spring
    fliplr(10 + yspring),... lower spring
    y0]; % lower right
writePLine(fid,[x;y]');
% ALSO NEED TO ADD SECTIONS WHICH CUT EXTERIOR; COULD SWITCH TO WINDOW
% LAYER HERE, BUT SHOULDN'T MATTER?
setCurrentLayer(fid,'window');
x = x0 + xspring - s2.srw;
y = 10 + yspring;
writePLine(fid,[x;y]');
x = x0 + s2.vgw - xspring + s2.srw;
y = 10 + yspring;
writePLine(fid,[x;y]');
x = x0 + xspring - s2.srw;
y = s1.wl - 10 - yspring;
writePLine(fid,[x;y]');
x = x0 + s2.vgw - xspring + s2.srw;
y = s1.wl - 10 - yspring;
writePLine(fid,[x;y]');
% FROM THIS POINT ON, ONE-TO-ONE-COPY OF ABOVE FUNCTION.
setCurrentLayer(fid,'window');
% fda is for edge of chip;
y0 = -(s2.vgl(1)-s2.lfda); yf = -s2.fel(1);
writeRectangle(fid,x0,y0,xf,yf);
y0 = s1.wl+s2.fel(2); yf = s1.wl + s2.vgl(2);
writeRectangle(fid,x0,y0,xf,yf);
setCurrentLayer(fid,'grid');
edge = 25; gap = s1.tg; % tiling gap
x0 = x0 + s2.vgw/2;
createSliceGrid(fid,x0, -(s2.vgl(1)-s2.lfda),x0+s2.vgw/2,s1.wl+s2.vgl(2),edge,gap);
createSliceGrid(fid,x0, -(s2.vgl(1)-s2.lfda),x0-s2.vgw/2,s1.wl+s2.vgl(2),edge,gap);
setCurrentLayer(fid,'vgrooveFDA');
FDAoverlap = 2;
writeRectangle(fid,x0+s2.vgw/2,FDAoverlap-s2.vgl(1),x0-s2.vgw/2,FDAoverlap+s2.lfda-s2.vgl(1));
setCurrentLayer(fid,'nitride');
end

function createVeeGrooveThatBulges(fid,s1,s2,x0)
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
writeOutput(fid,'; DRAWING BULGING VEE GROOVE');
x0 = x0 - s2.vgw/2; xf = x0 + s2.vgw;
radius = 10; % should match window!!
radius = 0; % should match window!!
overlap = 5 + radius;
y0 = -(s2.fel(1) + overlap); yf = (s2.fel(2) + overlap) + s1.wl;
s2 = setVeeGrooveDefaults(s2);
setCurrentLayer(fid,'vgroove');
bulge = s2.vgb;
x = [x0,x0,...
    x0 - bulge,x0 - bulge,...
    x0,x0,...
    xf,xf,...
    xf + bulge,xf + bulge,...
    xf,xf];
y = [y0,0,...
    0,s1.wl,...
    s1.wl,yf,...
    yf,s1.wl,...
    s1.wl,0,...
    0,y0];
writePLine(fid,[x;y]');
filletLastObject(fid,radius);
setCurrentLayer(fid,'window');
% fda is for edge of chip;
y0 = -(s2.vgl(1)-s2.lfda); yf = -s2.fel(1);
writeRectangle(fid,x0,y0,xf,yf);
y0 = s1.wl+s2.fel(2); yf = s1.wl + s2.vgl(2);
writeRectangle(fid,x0,y0,xf,yf);
% here we gotta augment the region which gets tiled to incorporate
setCurrentLayer(fid,'grid');
edge = 25; gap = s1.tg; % tiling gap
x0 = x0 + s2.vgw/2;
% bottom vee
createSliceGrid(fid,x0, -(s2.vgl(1)-s2.lfda),x0-s2.vgw/2,0,edge,gap);
createSliceGrid(fid,x0, -(s2.vgl(1)-s2.lfda),x0+s2.vgw/2,0,edge,gap);
% in window
createSliceGrid(fid,x0,0,x0-s2.vgb-s2.vgw/2,s1.wl,edge,gap);
createSliceGrid(fid,x0,0,x0+s2.vgb+s2.vgw/2,s1.wl,edge,gap);
% top vee
createSliceGrid(fid,x0,s1.wl,x0-s2.vgw/2,s1.wl+s2.vgl(2),edge,gap);
createSliceGrid(fid,x0,s1.wl,x0+s2.vgw/2,s1.wl+s2.vgl(2),edge,gap);
setCurrentLayer(fid,'vgrooveFDA');
FDAoverlap = 2;
writeRectangle(fid,x0+s2.vgw/2,FDAoverlap-s2.vgl(1),x0-s2.vgw/2,FDAoverlap+s2.lfda-s2.vgl(1));
setCurrentLayer(fid,'nitride');
end

function [xf,yf,wf,bwf] = createYJunction(fid,s,x0,y0,w0)
s.w0 = w0;
if s.newt % if new taper
    % calculate initial gap
    g_c = (s.gap + s.w - s.A);
    s.gap = g_c - s.wf;
    coords = drawNewYJunction(s);
else
    coords = drawYJunction(s);
end
% writeOutput(fid,'; DRAWING Y JUNCTION');

coords = mirrorCoords(coords);
coords = displaceCoords(coords,x0,y0,w0);
% return width is a little tricky here
[xf,yf,wf] = writePLine(fid,coords);
bwf = s.wf;
%
end

function [xf,yf,wf,bwf] = createNegYJunction(fid,s,x0,y0,w0)
s.w0 = w0;
if s.newt % if new taper--not supported for NegY yet
    %     % calculate initial gap
    %     g_c = (s.gap + s.w - s.A);
    %     s.gap = g_c - s.wf;
    %     coords = drawNewYJunction(s);
else
    coords = drawYJunction(s);
    for jj = 1:2
        [coordsInt,coordsExt] = toneInvertYJunction(coords);
        if jj == 2
            coordsInt(:,1) = -coordsInt(:,1);
            coordsExt(:,1) = -coordsExt(:,1);
        end
        coordsInt = displaceCoords(coordsInt,x0,y0,w0);
        coordsExt = displaceCoords(coordsExt,x0,y0,w0);
        writePLine(fid,coordsInt);
        addString(fid,'(setq currobj (ssadd (entlast) currobj))');
        writePLine(fid,coordsExt);
        addString(fid,'(setq currobj (ssadd (entlast) currobj))');
    end
    yf = y0 + s.ltaper;
    wf = s.gap + 2*s.wf;
    bwf = s.wf;
    xf = x0 + (-s.w0 + wf)/2;
end
% coords = mirrorCoords(coords);
% coords = displaceCoords(coords,x0,y0,w0);
    function [coordsInt,coordsExt] = toneInvertYJunction(coords)
        excur = (s.A + s.w + s.gap/2) * 2; % excursion needs to change if crystal changes
        coordsExt = [coords(1:2,:);
            [excur,coords(2,2)];
            [excur,coords(1,2)]]; % exterior
        coordsInt = [coords(3:4,:);coords(4,1),coords(3,2)]; % interior triangle
    end
end

function defineLISPfunctions(fid)
% These may or may not be necessary; I forget if they are used.
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
writeOutput(fid,['(defun x-coords (input-coords) (mapcar ''(lambda (x)',...
    ' (nth 0 x)) input-coords))']);
writeOutput(fid,['(defun y-coords (input-coords) (mapcar ''(lambda (x)',...
    ' (nth 1 x)) input-coords))']);
writeOutput(fid,['(defun remove-nil (in) (apply ''append (mapcar',...
    ' ''(lambda (x) (if (not (null x)) (list x) ) ) in) ))']);
writeOutput(fid,['(defun get-coords (input-entity) (remove-nil (mapcar',...
    ' ''(lambda (x) (if (= (car x) 10) (cdr x))) input-entity)))']);
% function to zoom to named objects
writeOutput(fid,[
    '(defun c:zoom-named (entname) ',...
    '(command "_zoom" "object" entname "")',...
    '(princ)',...
    ')'
    ]);
writeOutput(fid,'(if (not (numberp (vl-string-search "Mac" (getvar "platform")))) (load "TXTEXP.LSP"))');
end

function coords = displaceCoords(coords,x,y,w)
coords(:,1) = coords(:,1) + x - (w/2);
coords(:,2) = coords(:,2) + y;
end

function coords = drawDPHC(s,fa,fw,fg)
% npoints = 20; % number of points per cycle
npoints = 37; % number of points per cycle
ncycles = sum([s.A(:).n]); % number of cycles
% Draw RHS
x = linspace(0,ncycles,ncycles*npoints);
y = fg(x)/2 + fw(x) + fa(x).*cos(2*pi*x);
y_int = fg(x)/2;
% prepend and postpend % REPLACED 2015/09/21 TO ACCOMODATE NEW TAPER
% x = [x(1),x,x(end)];
% y = [0,y,0];
% % shift by half gap
% y = y + s.gap/2;
% scale by lattice constant
x = x*s.a;
coords(:,2) = [x,fliplr(x)];
coords(:,1) = [y,fliplr(y_int)];
end

function coords = drawFillet(s)
xi = s.p0(1); xf = s.pf(1);
yi = s.p0(2); yf = s.pf(2);
rx = xf - xi; ry = yf - yi;
switch s.orientation
    case 'vertical'
        xc = xf;
        yc = yi;
        theta_list = linspace(pi,pi/2,s.num);
    case 'horizontal'
        xc = xi;
        yc = yf;
        theta_list = linspace(3*pi/2,2*pi,s.num);
end
j = 1;
coords = zeros(numel(theta_list),2);
for theta = theta_list
    coords(j,:) = [rx*cos(theta)+xc,ry*sin(theta)+yc];
    j = j + 1;
end
end

function coords = drawNegativeTethers(s)
s = setNodeArrayDefaults(s);
% j = 1;
% coords = zeros((s.numTethers - 1) * 4,2);
coords = [];
for m = 1:(s.numTethers-1)
    coords_local = [s.width/2, s.tetherCenterWidth + (m - 1) * s.a;
        s.veeGrooveWidth/2, (m - 1) * s.aout + s.tetherYDisp + s.tetherEdgeWidth;
        s.veeGrooveWidth/2, m * s.aout + s.tetherYDisp;
        s.width/2, m*s.a];
    coords_local = addSpring(s,coords_local);
    % coords(j:j+3,:) = coords_local;
    coords = [coords;coords_local];
    % j = j + 4;
end
end

function wf = drawSegNegTaper(fid,dev,s,rng,amp_fun,wid_fun,gap_fun,x0,y0,w0,xc,yc)
%
% HAD TO ADD THIS KLUDGE TO COAX REASONABLE FRACTURE.
%
npts = 37;
f = @(x) gap_fun(x)/2 + wid_fun(x) + amp_fun(x).*cos(2*pi*x);
ext_max = 2*max(f(rng));
for jj = (rng(1)+1):rng(end) % loop over
    x_curr = linspace(jj-1,jj,npts);
    ext_int = max(f(x_curr));
    for kk = 1:3 % interior, exterior and
        if kk == 1
            cds = drawInteriorCoords(x_curr);
        elseif kk == 2
            cds = drawExteriorCoords(x_curr);
        else
            cds = drawRectangleFill(x_curr);
        end
        cds = displaceCoords(cds,x0,y0,w0);
        [~,~,~] = writePLine(fid,cds);
        addToSelectionSet;
        mirrorObject;
        addToSelectionSet;
    end
end

    function coords_int = drawInteriorCoords(xx)
        coords_int = [xx;gap_fun(xx)/2]; % interior
        coords_int = [[xx(1);0],coords_int,[xx(end);0]]; % bookend
        coords_int = coords_int'; coords_int = fliplr(coords_int);
        coords_int(:,2) = coords_int(:,2)*s.a; % convert to correct periodicity
    end

    function coords_ext = drawExteriorCoords(xx)
        coords_ext = [xx;f(xx)];
        coords_ext = [[xx(1);ext_int],coords_ext,[xx(end);ext_int]];
        coords_ext = coords_ext'; coords_ext = fliplr(coords_ext);
        coords_ext(:,2) = coords_ext(:,2)*s.a; % convert to the correct periodicity
    end

    function coords_rect = drawRectangleFill(xx)
        coords_rect = [xx;ones(size(xx)) * ext_int];
        coords_rect = [[xx(1);ext_max],coords_rect,[xx(end);ext_max]]; % bookend
        coords_rect = coords_rect'; coords_rect = fliplr(coords_rect);
        coords_rect(:,2) = coords_rect(:,2)*s.a; % convert to correct periodicity
    end

    function addToSelectionSet
        addString(fid,['(setq tap',num2str(dev),' (ssadd (entlast) tap',num2str(dev),'))']); % create tap1 selection set
    end

    function mirrorObject
        mirrorNamedObject(fid,'(entlast)',[xc,yc],[xc,yc+1]); % reflect
    end

wf = 2 * (max(cds(:,1))-xc); % final width
end

function [coords_int,coords_ext] = drawNegTaper(s,rng,amp_fun,wid_fun,gap_fun)
%
% This function is cryptic, but I promise what it does is draw a box on the
% interior and on the exterior of the photonic crystal.
%
npts = 37;
x_coords = linspace(rng(1),rng(2),npts*(rng(2)-rng(1)+1));
coords_int = [x_coords;gap_fun(x_coords)/2]; % interior
coords_int = [[x_coords(1);0],coords_int,[x_coords(end);0]];
f = @(x) gap_fun(x)/2 + wid_fun(x) + amp_fun(x).*cos(2*pi*x);
coords_ext = [x_coords;f(x_coords)];
ext_max = 2*max(coords_ext(2,:));
coords_ext = [[x_coords(1);ext_max],coords_ext,[x_coords(end);ext_max]];
coords_int = coords_int'; coords_int = fliplr(coords_int);
coords_ext = coords_ext'; coords_ext = fliplr(coords_ext);
% convert to the correct periodicity
coords_int(:,2) = coords_int(:,2)*s.a;
coords_ext(:,2) = coords_ext(:,2)*s.a;
end

function [coords_int,coords_ext] = drawNegUnitCell(s,amp_fun,wid_fun,gap_fun)
% draw a unit cell centered about the origin
npts = 37;
x = linspace(s.nt+1,s.nt+2,npts); % cell just after the taper
coords_int = [gap_fun(x)/2;x];
coords_int = coords_int';
coords_int = [0,coords_int(1,2);coords_int;0,coords_int(end,2)];
f = @(x) gap_fun(x)/2 + wid_fun(x) + amp_fun(x).*cos(2*pi*x);
coords_ext = [f(x);x];
coords_ext = coords_ext';
ext_max = 2*max(coords_ext(:,1));
% ext_max = max(coords_ext(:,1));
coords_ext = [ext_max,coords_ext(1,2);coords_ext;ext_max,coords_ext(end,2)];
coords_int(:,2) = coords_int(:,2) * s.a;
coords_ext(:,2) = coords_ext(:,2) * s.a;
yc = (max(coords_ext(:,2))+min(coords_ext(:,2)))/2;
coords_int(:,2) = coords_int(:,2) - yc;
coords_ext(:,2) = coords_ext(:,2) - yc;
end

function coords = drawNodeArray(s)
% set defaults if not specified
s = setNodeArrayDefaults(s);
% center inner and outer tethers
s.tetherYDisp = s.tetherYDisp - (s.numTethers-1)*(s.aout - s.a)/2;
% draw RHS
% single nanobeam, linear taper
s1.length = s.taperLength;
s1.w0 = s.w0;
s1.wf = s.width;
s1.fun = @(w0,wf,y) w0 + (wf - w0).* y.^1;
s1.num = 2;
coords = drawSingleNanobeam(s1);
% tether array
s2 = s;
s2.nominalWidth = s.width;
s2.prelength = s2.blankLength;
s2.postlength = s2.blankLength;
s2 = rmfield(s2,{'blankLength','w0','wf','taperLength','width'});
coords2 = drawTether(s2);
coords2(:,2) = coords2(:,2) + coords(end,2);
coords = [coords; coords2];
% single nanobeam, linear taper
s3.length = s.taperLength;
s3.w0 = s.width;
s3.wf = s.wf;
s3.fun = @(w0,wf,y) w0 + (wf - w0).* y.^1;
s3.num = 2;
coords2 = drawSingleNanobeam(s3);
coords2(:,2) = coords2(:,2) + coords(end,2);
coords = [coords; coords2];
end

function drawRegionLayer(fid,dev,s1,s2,x0)
% New Sept 22 2015.
mainfield_size = 160;
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
writeOutput(fid,'; DRAWING REGION LAYER')
x0 = x0 - s2.vgw/2 - s2.vgb; xf = x0 + s2.vgw + s2.vgb * 2;
y0 = s2.overhang; yf = s1.wl - s2.overhang;
% calculate center and extent
x_mean = mean([x0,xf]);
y_mean = mean([y0,yf]);
x_ext = abs(diff([x0,xf]));
y_ext = abs(diff([y0,yf]));
% round to odd number of mainfields
x_n = ceil(x_ext./mainfield_size);
y_n = ceil(y_ext./mainfield_size);
if ~mod(x_n,2) % if odd
    x_n = x_n + 1;
end
if ~mod(y_n,2) % if odd
    y_n = y_n + 1;
end
% correct x and y
x0 = x_mean - (x_n * mainfield_size)/2;
xf = x_mean + (x_n * mainfield_size)/2;
y0 = y_mean - (y_n * mainfield_size)/2;
yf = y_mean + (y_n * mainfield_size)/2;
% draw rectangle
setCurrentLayer(fid,['regionLayer',num2str(dev)]);
writeRectangle(fid,x0,y0,xf,yf);
setCurrentLayer(fid,'nitride');
end

function coords = drawSingleNanobeam(s)
% set defaults if not specified
s = setSingleNanobeamDefaults(s);
% draw RHS
coords(:,1) = s.fun(s.w0/2,s.wf/2,linspace(0,1,s.num));
coords(:,2) = linspace(0,s.length,s.num);
end

function [xc,yc] = drawStructure(fid,dev,s1,s2,s3,s4,s5,s6,s7)
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
[x0,y0] = setOffsets(dev,s1.spacing,s2(dev).overhang);
addString(fid,['; DEVICE #',num2str(dev)]); % First thing to be drawn, adds
% comment.
%% label which device
addDevMarkers(fid,x0,y0,s2(dev).vgw,dev);
%% draw coupler
addString(fid,'(setq currobj nil)');        % Free currobj
addString(fid,'(setq currobj (ssadd))');    % I FORGET WHY THIS IS HERE
writeOutput(fid,'; DRAWING COUPLER')
[xf,yf,wf] = createCoupler(fid,s3(dev),x0,y0);
addString(fid,['(setq coupler',num2str(dev),' currobj)']); % and call it
namedObjects = {['coupler',num2str(dev)]};
%% draw single nanobeam
addString(fid,'(setq currobj nil)');        % Free currobj
addString(fid,'(setq currobj (ssadd))');    % I FORGET WHY THIS IS HERE
writeOutput(fid,'; DRAWING SINGLE NANOBEAM');
[xf,yf,wf] = createSingleNanobeam(fid,s4(dev),xf,yf,wf);
addString(fid,['(setq singlenanobeam',num2str(dev),' currobj)']);
namedObjects{end+1} = ['singlenanobeam',num2str(dev)];
%% draw node array
addString(fid,'(setq currobj nil)');        % Free currobj
addString(fid,'(setq currobj (ssadd))');    % I FORGET WHY THIS IS HERE
writeOutput(fid,'; DRAWING NODE ARRAY');
s5(dev).vgb = s2(dev).vgb;
[xf,yf,wf,~] = createNodeArray(fid,s5(dev),xf,yf,wf);
addString(fid,['(setq anglednodearray',num2str(dev),' currobj)']);
% coupler
namedObjects{end+1} = ['anglednodearray',num2str(dev)];
%% draw cooled single nanobeam
addString(fid,'(setq currobj nil)');
addString(fid,'(setq currobj (ssadd))');
writeOutput(fid,'; DRAWING COOLED SINGLE NANOBEAM')
[xf,yf,wf,tcp] = createCooledSingleNanobeam(fid,s6(dev),xf,yf,wf);
addString(fid,['(setq cooledsinglenanobeam',num2str(dev),' currobj)']);
namedObjects{end+1} = ['cooledsinglenanobeam',num2str(dev)];
%% draw photonic crystal
if s7(dev).dnb
    [xc,yc,namedObjects] = createDoubleNanobeamCrystal(fid,dev,s7(dev),xf,yf,wf,namedObjects);
else
    [xc,yc,namedObjects] = createSingleNanobeamCrystal(fid,dev,s7(dev),xf,yf,wf,namedObjects);
end
%% draw indicators
createIndicators(fid,s2(dev),xc,yc);
%% draw reflections of named objects
reflectObjects(fid,xc,yc,namedObjects);
%% add positions of tethers in second half of window
tcp = sort([tcp,yc+(yc-tcp)]);
%% draw vee groove, window
if s1.window
    if s1.struts % if struts
        createVeeGrooveThatBulges(fid,s1,s2(dev),xc);
        % createVeeGroove(fid,s1,s2(dev),xc);
        createRailHoles(fid,s1,s2(dev),xc);
        createWindowWithStruts(fid,s1,s2(dev),dev,xc);
    elseif s2(dev).vgs % if spring
        createVeeGrooveWithSprings(fid,s1,s2(dev),xc);
        createRailHolesWithSprings(fid,s1,s2(dev),xc);
        createWindow(fid,s1,s2,tcp,dev,xc);
    else
        createVeeGrooveThatBulges(fid,s1,s2(dev),xc);
        % createVeeGroove(fid,s1,s2(dev),xc);
        createRailHoles(fid,s1,s2(dev),xc);
        createWindow(fid,s1,s2,tcp,dev,xc);
        % createWindowThatBulges(fid,s1,s2(dev),dev,xc);
    end
else
    createVeeGrooveThatBulges(fid,s1,s2(dev),xc);
end
%% draw region layer
drawRegionLayer(fid,dev,s1,s2(dev),xc);
%% clear selection sets
clearSelectionSets(fid,namedObjects);
end

function coords = drawTether(s)
% set defaults if not specified
s = setTetherDefaults(s);
% allows tethers to fan
s.tetherYDisp = s.tetherYDisp - (s.numTethers - 1)*(s.aout - s.a)/2;
% draw RHS
coords = [s.nominalWidth/2 0; s.nominalWidth/2 s.prelength];
% EDITED 2015-07-15
% add spring
if (s.fillet ~= 0)
    coords = [coords; tethersWithFillets(s)];
else
    coords = [coords; tethersWithoutFillets(s)];
end
coords = [coords; [s.nominalWidth/2 coords(end)] + [0,s.postlength]];
end

function coords = drawNewYJunction(s)
% Modified y junction 2015/09/21
s.offset = (s.gap + s.wf)/2;
coords = [s.w0/2, 0;
    2 * s.offset, s.lsplit; % this is an additional point
    s.offset + s.wf/2, s.ltaper;
    s.offset - s.wf/2, s.ltaper;
    0, s.lsplit
    ];
end

function coords = drawYJunction(s)
s.offset = (s.gap + s.wf)/2;
coords =[0.5*s.w0, 0; s.offset + 0.5 * s.wf, s.ltaper;
    s.offset - 0.5*s.wf, s.ltaper; 0, s.lsplit];
end

function eraseAll(fid)
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
writeOutput(fid,'; ERASE ALL');
writeOutput(fid,'(command "_erase" "all" "")');
end

function filletLastObject(fid,radius)
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
str = ['(command "._fillet" "radius" "',num2str(radius),...
    '" "._fillet" "P" (entlast))'];
writeOutput(fid,str);
end

function str = getStringVal(strcl,mtchstr)
% str = getStringVal(strcl,mtchstr) returns a string or numeric value
% associated with a parameter in the .param file.
try str = strcl{cellfun(@(x) strncmp(x,mtchstr,length(mtchstr)),...
        strcl)}(length(mtchstr)+1:end);
catch
    warning([mtchstr,' undefined!'])
    str = '';
end
if isempty(regexp(str,'[A-Za-z]','once'))
    str = str2num(str);
elseif strncmp(str,' ',1)
    str(1) = [];
end
end

function str = getStringEVal(strcl,mtchstr)
% str = getStringVal(strcl,mtchstr) returns a string or numeric value
% associated with a parameter in the .param file.
str = strcl{cellfun(@(x) strncmp(x,mtchstr,length(mtchstr)),...
    strcl)}(length(mtchstr)+1:end);
if isempty(regexp(str,'[A-Za-z]','once'))
    str = str2num(str);
elseif strncmp(str,' ',1)
    str = eval(str);
end
end

function hideLayer(fid,layerName)
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
writeOutput(fid,'; DRAW LAYERS, SET COLORS');
str = ['(command "_layer" "OFF" "',layerName,'" "")'];
writeOutput(fid,str);
end

function [s1,s2,s3,s4,s5,s6,s7] = importParameters(fn)
% s = importParameters() prompts the user for a .param file and returns a
% struct of structs which define the global geometry for the file. It is
% based on populateStructs(), but printed here to improve portability.

%% Get file, read lines into cell array
if isempty(fn)
    [fn,pn] = uigetfile('*.param','Select param file');
    fn = [pn,'/',fn];
end
fid = fopen(fn{1});
strcl = textscan(fid,'%s','delimiter','\n');
strcl = strcl{1};
fclose(fid);

%% Find sectioning markers
sectionIndices = cellfun(@(x) strncmp(x,'//',length('//')),strcl);
ix = 1:numel(sectionIndices);
sectionIndices = ix(sectionIndices);
sectionIndices = [sectionIndices,numel(strcl)];

%% Populate global struct
% Warning: cannot scan global parameters! By this logic, vee groove
% parameters shouldn't be included in the vee groove struct.
currIx = cellfun(@(x) strncmp(x,'// global',length('// global')),...
    strcl(sectionIndices));
subcl = strcl(sectionIndices(currIx):sectionIndices(circshift(currIx,1)));
s1.filename = getStringVal(subcl,'filename:');        % filename
s1.window = interpretTruth(getStringVal(subcl,...     % window?
    'window:'));
if paramExists(subcl,'struts:')
    s1.struts = interpretTruth(getStringVal(subcl,...
        'struts:'));
else
    s1.struts = false;
end
s1.num = getStringVal(subcl,'number of devices:');    % number of devices
s1.spacing = getStringVal(subcl,'device spacing:');   % structure spacing
s1.wl = getStringEVal(subcl,'window length:');         % length of window
s1.ww = getStringVal(subcl,'window width:');          % width of window
s1.tg = getStringVal(subcl,'tiling gap:');
if paramExists(subcl,'window fillet radius:')
    s1.wfr = getStringVal(subcl,'window fillet radius:');
else
    s1.wfr = 5; % defaults to 5 um unless script specifies
end

if ~paramExists(subcl,'perpendicular grooves:') % FOR FLIP CHIP
    s1.pg = false;
else
    s1.pg = interpretTruth(getStringVal(subcl,'perpendicular grooves:'));
end

if ~paramExists(subcl,'flip chip groove width:') % FOR FLIP CHIP
    s1.fcgw = 0;
else
    s1.fcgw = getStringVal(subcl,'flip chip groove width:');
end


%% Check for scanned parameter
currIx = cellfun(@(x) strncmp(x,'// scan',length('// scan')),...
    strcl(sectionIndices));
sbln = sum(currIx); % scan boolean
if sbln
    subcl = strcl(sectionIndices(currIx):sectionIndices(circshift(currIx,1)));
    scanStr = getStringVal(subcl,'scanned parameter:');
    scanStart = getStringVal(subcl,'start value:');
    scanEnd = getStringVal(subcl,'end value:');
    scanVal = arrayfun(@num2str,linspace(scanStart,scanEnd,s1.num),...
        'UniformOutput',0);
end

%% Loop over scanned parameter, make replacements
s2 = struct([]);
s3 = s2; s4 = s2; s5 = s2; s6 = s2; s7 = s2;
for dev = 1:s1.num
    if sbln
        devcl = replaceScanString(strcl,scanStr,scanVal{dev});
    else
        devcl = strcl;
    end
    %% Populate grooves and safety rails
    currIx = cellfun(@(x) strncmp(x,'// grooves',length('// grooves')),...
        devcl(sectionIndices));
    subcl = devcl(sectionIndices(currIx):sectionIndices(circshift(currIx,1)));
    if paramExists(subcl,'back vee groove:')
        s2(dev).backGroove = interpretTruth(getStringVal(subcl,... % back vee groove?
            'back vee groove:'));
    else
        s2(dev).backGroove = true; % COERCED 2017-05-08
    end
    s2(dev).vgw = getStringVal(subcl,'vee groove width:');     % width of vee groove
    if paramExists(subcl,'vee groove bulge:')
        s2(dev).vgb = getStringVal(subcl,'vee groove bulge:');
    else
        s2(dev).vgb = 0;
    end
    if paramExists(subcl,'vee groove spring:')
        s2(dev).vgs = interpretTruth(getStringVal(subcl,'vee groove spring:'));
        if s2(dev).vgs
            s2(dev).sa = getStringVal(subcl,'spring amplitude:');
            s2(dev).sp = getStringVal(subcl,'spring pitch:');
        end
    else
        s2(dev).vgs = false;
    end
    s2(dev).vgl = getStringEVal(subcl,'vee groove length:');    % length
    s2(dev).lfda = getStringVal(subcl,'fda length:');          % ???
    s2(dev).overhang = getStringVal(subcl,...                  % length of device
        'length of device past window:');                      % before window
    s2(dev).fel = getStringVal(subcl,'fine exposure length:');
    if ischar(s2(dev).fel)
        s2(dev).fel = str2num(s2(dev).fel);
    end
    s2(dev).srw = getStringVal(subcl,'safety rail width:');
    s2(dev).ha = getStringVal(subcl,'hole lattice constant:');
    s2(dev).nhr = getStringVal(subcl,'number of hole rows:');
    s2(dev).hr = getStringVal(subcl,'hole radius:');
    if paramExists(subcl,'hole recess:')
        s2(dev).he = getStringVal(subcl,'hole recess:');
    else
        s2(dev).he = 0;
    end
    %% Populate coupler
    currIx = cellfun(@(x) strncmp(x,'// coupler',length('// coupler')),...
        devcl(sectionIndices));
    subcl = devcl(sectionIndices(currIx):sectionIndices(circshift(currIx,1)));
    s3(dev).nominalWidth = getStringEVal(subcl,'nominal width:');
    s3(dev).veeGrooveWidth = s2(dev).vgw;
    s3(dev).tetherCenterWidth = getStringEVal(subcl,'tether center width:');
    s3(dev).prelength = getStringVal(subcl,'prelength:');
    s3(dev).postlength = getStringVal(subcl,'postlength:');
    s3(dev).fillet = getStringVal(subcl,'fillet radius:');
    s3(dev).tccb = interpretTruth(getStringVal(subcl,'tcc:'));
    if s3(dev).tccb
        s3(dev).hole = getStringEVal(subcl,'hole radius:');
        s3(dev).tccr = getStringVal(subcl,'tcc radius:');
    end
    
    %% Populate nanobeam
    currIx = cellfun(@(x) strncmp(x,'// single nanobeam',...
        length('// single nanobeam')),devcl(sectionIndices));
    subcl = devcl(sectionIndices(currIx):sectionIndices(circshift(currIx,1)));
    s4(dev).length = getStringVal(subcl,'length:');
    if paramExists(subcl,'start width:')
        s4(dev).w0 = getStringVal(subcl,'start width:');
    else
        s4(dev).w0 = s3(dev).nominalWidth;
    end
    s4(dev).wf = getStringVal(subcl,'end width:');
    s4(dev).fun = getStringVal(subcl,'function:');
    s4(dev).fun = eval(s4(dev).fun);
    s4(dev).num = getStringVal(subcl,'number of points:');
    
    %% Populate node array
    currIx = cellfun(@(x) strncmp(x,'// angled node array',...
        length('// angled node array')),devcl(sectionIndices));
    subcl = devcl(sectionIndices(currIx):sectionIndices(circshift(currIx,1)));
    if paramExists(subcl,'start width:')
        s5(dev).w0 = getStringVal(subcl,'start width:');
    else
        s5(dev).w0 = s4(dev).wf;
    end
    s5(dev).wf = getStringVal(subcl,'end width:');
    s5(dev).taperLength = getStringVal(subcl,'taper length:');
    s5(dev).numTethers = getStringVal(subcl,'number of tethers:');
    s5(dev).blankLength = getStringVal(subcl,'blank length:');
    s5(dev).width = getStringVal(subcl,'width at tethering:');
    s5(dev).tetherCenterWidth = getStringVal(subcl,'tether center width:');
    s5(dev).a = getStringVal(subcl,'lattice constant:');
    s5(dev).aout = getStringVal(subcl,'ext. lattice constant:');
    s5(dev).veeGrooveWidth = s2(dev).vgw;
    s5(dev).tetherYDisp = getStringVal(subcl,'tether displacement:');
    % Need to coerce to end width; see param file
    
    %% Populate single nanobeam
    currIx = cellfun(@(x) strncmp(x,'// cooled single nanobeam',...
        length('// cooled single nanobeam')),devcl(sectionIndices));
    subcl = devcl(sectionIndices(currIx):sectionIndices(circshift(currIx,1)));
    if paramExists(subcl,'start width:')
        s6(dev).w0 = getStringVal(subcl,'start width:');
    else
        s6(dev).w0 = s5(dev).wf;
    end
    s6(dev).wf = getStringVal(subcl,'end width:');
    s6(dev).taperLength = s5(dev).taperLength;
    s6(dev).numTethers = s5(dev).numTethers;
    s6(dev).blankLength = s5(dev).blankLength;
    s6(dev).width = s5(dev).width;
    s6(dev).tetherCenterWidth = s5(dev).tetherCenterWidth;
    s6(dev).a = s5(dev).a;
    s6(dev).aout = s5(dev).aout;
    s6(dev).veeGrooveWidth = s5(dev).veeGrooveWidth;
    s6(dev).lstart = getStringVal(subcl,'start length:');
    s6(dev).lstop =  getStringVal(subcl,'stop length:');
    s6(dev).tetherYDisp = getStringVal(subcl,'tether displacement:');
    s6(dev).spacing = getStringVal(subcl,'spacing:');
    if paramExists(subcl,'length:')
        s6(dev).length = getStringVal(subcl,'length:');
    else
        s6(dev).length = s1.wl; % gets replaced in length calculation anyway
    end
    s6(dev).lnode = lengthNodeArray(s5(dev));
    if paramExists(subcl,'spring amplitude:')
        s6(dev).samp = getStringVal(subcl,'spring amplitude:');
    else
        s6(dev).samp = 0;
    end
    s6(dev).vgb = s2(dev).vgb; % vee groove bulge
    
    %% photonic crystal
    currIx = cellfun(@(x) strncmp(x,'// photonic crystal',...
        length('// photonic crystal')),devcl(sectionIndices));
    subcl = devcl(sectionIndices(currIx):sectionIndices(circshift(currIx,1)));
    % Here we fork between single and double nanobeam crystals
    if paramExists(subcl,'type:')
        dbl = strncmp(getStringVal(subcl,'type:'),'double',length('double'));
    else
        dbl = true;
    end
    s7(dev).dnb = dbl;
    s7(dev).neg = false; % THIS IS TO DO CELL FRACTURING 2015-12-04
    if paramExists(subcl,'inverse crystal:')
        if interpretTruth(getStringVal(subcl,'inverse crystal:'));
            s7(dev).neg = true;
            s7(dev).segTaper = true; % SEGMENTS TAPER; SHOULD PROBABLY BE IN .param FILE, BUT I AM A SINNER
        end
    end
    % SHOULD WE DO Y JUNCTION TONE REVERSE AS WELL?
    s7(dev).yneg = false;
    if paramExists(subcl,'inverse y:')
        if interpretTruth(getStringVal(subcl,'inverse y:'))
            s7(dev).yneg = true;
            
        end
    end
    if dbl
        % y-junction
        if paramExists(subcl,'start width:')
            s7(dev).w0 = getStringVal(subcl,'start width:');
        else
            s7(dev).w0 = s6(dev).wf;
        end
        s7(dev).wf = getStringVal(subcl,'end width:');
        % s7(dev).gap = getStringVal(subcl,'gap:');
        s7(dev).ltaper = getStringVal(subcl,'taper length:');
        s7(dev).lsplit = getStringVal(subcl,'split length:');
        if paramExists(subcl,'new taper:')
            s7(dev).newt = interpretTruth(getStringVal(subcl,'new taper:'));
        else
            s7(dev).newt = 0;
        end
        % double nanobeam
        if paramExists(subcl,'length:')
            s7(dev).length = getStringVal(subcl,'length:');
        else
            s7(dev).length = 0;
        end
        % crystal
        s7(dev).a = getStringVal(subcl,'lattice constant:');
        s7(dev).gap = getStringVal(subcl,'gap:');
        s7(dev).A = getStringEVal(subcl,'amplitude:');
        if ischar(s7(dev).A)
            % load .mat file, select target device
            load(s7(dev).A);
            s7(dev).A = s(dev).A;
        end
        s7(dev).w = getStringEVal(subcl,'width:');
        if ischar(s7(dev).w)
            % load .mat file, select target device
            load(s7(dev).w);
            s7(dev).w = s(dev).w;
        end
        s7(dev).nt = getStringVal(subcl,'number taper:');
        s7(dev).n = getStringEVal(subcl,'number:');
    else % single beam
        % NEED TO REWORK THIS SECTION; WANT TO CREATE DOUBLE-BEAM CAVITIES
        % TOO. AMPLITUDE MODULATION, ETC. SHOULD BE NECESSARY.
    end
end
end

function initializeDrawing(fid,s1)
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
%% include metadata
writeOutput(fid,[';; NUMBER OF LAYERS: ',num2str(s1.num)]);
%% instantiate layers
layerNames = {'negYjunc','nitride','slices','holes','vgroove','vgrooveFDA',...
    'window','marks','extent','cover','grid','pads','text','fda1',...
    'prebias','phc','phccover','tethers','indicators','flipChip',...
    'flipChipSubtract','flipChipWindow'};
nlayers = numel(layerNames);
for kk = 1:s1.num
    layerNames{nlayers+kk} = ['regionLayer',num2str(kk)];
end
layerColors = {'green','white','red','blue','green','red','yellow','yellow',...
    'yellow','yellow','176','red','white','1','6','5'};
eraseAll(fid);                                  % erase all
makeLayers(fid,layerNames,layerColors);         % make layers
zoomToRectangle(fid,[-5e-2 -5e-2],[5e-2 5e-2]); % standard zoom
% PURGE BLOCKS
writeOutput(fid,'(command "_purge" "B" "*" "N")')
end

function insertBlock(fid,blockName,ctr,xscale,yscale,rotationAngle)
% should probably edit this command to specify specifically which numbers
% are where in the command. Ordinarily done by mnemonic.
insertionPoint = [num2str(ctr(1)),',',num2str(ctr(2))];
str = ['(command "_insert" "'...
    blockName,...
    '" "',...
    insertionPoint,...
    '" "',...
    num2str(xscale),...
    '" "',...
    num2str(yscale),...
    '" "',...
    num2str(rotationAngle),...
    '")'];
addString(fid,str)
end

function bln = interpretTruth(str)
% bln = interpretTruth(str) converts 'yes' and 'no' values to their boolean
% equivalents in the .param file.
if strncmp(str,'yes',length('yes'))
    bln = true;
elseif strncmp(str,'no',length('no'))
    bln = false;
end
end

function length = lengthNodeArray(s)
length = s.a * (s.numTethers-1) + s.tetherCenterWidth + ...
    2 * s.blankLength + 2 * s.taperLength;
end

function makeLayers(fid,layerNames,layerColors)
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
writeOutput(fid,'; DRAW LAYERS, SET COLORS');
cellfun(@(x) writeOutput(fid,createLayer(x)),layerNames);
layerNames(numel(layerColors)+1:end) = []; % only spec colors for first N layers
cellfun(@(x,y) writeOutput(fid,colorLayer(x,y)),layerColors,layerNames);
setCurrentLayer(fid,layerNames{1}); % make first layer current
end

function mirrorNamedObject(fid,name,x0,x1)
addString(fid,['(command "_mirror" ',name,' "" "',num2str(x0(1)),...
    ',',num2str(x0(2)),'" "',num2str(x1(1)),',',num2str(x1(2)),'" "No")']);
end

function coords = mirrorCoords(coords)
coords = [coords; flipud([-coords(:,1) coords(:,2)])];
end

function sOut = nodeArrayToCover(sIn)
sOut = sIn;
sOut = setNodeArrayDefaults(sOut);
sOut.tetherCenterWidth = sOut.tetherCenterWidth + ...
    (sOut.numTethers - 1) * sOut.a;
sOut.tetherEdgeWidth = sOut.tetherEdgeWidth + ...
    sOut.aout * (sOut.numTethers - 1);
sOut.numTethers = 1;
end

function bln = paramExists(cell,target)
bln = sum(cellfun(@(x) strncmp(x,target,length(target)),cell));
end

function prebiasNodeArray(fid,coords)
% This routine finds anchor points for the node array and draws a polyline
% around it. Kind of stupid, very by-hand. Gets the job done?
setCurrentLayer(fid,'prebias');
% LHS
min_x = min(coords(:,1));
bln_min = find(coords(:,1)==min_x);
A3 = [min_x,min(coords(bln_min,2))]; % anchor 3
ix = find(coords(:,1)==A3(:,1)&coords(:,2)==A3(:,2));
while coords(ix,1) == A3(1)
    ix = ix + 1; % decrement
end
A2 = coords(ix,:);
A4 = [min_x,max(coords(bln_min,2))];
ix = find(coords(:,1)==A4(:,1)&coords(:,2)==A4(:,2));
while coords(ix,1) == A4(1)
    ix = ix - 1; % increment
end
A5 = coords(ix,:);
% RHS
max_x = max(coords(:,1));
bln_min = find(coords(:,1)==max_x);
A10 = [max_x,min(coords(bln_min,2))];
ix = find(coords(:,1)==A10(:,1)&coords(:,2)==A10(:,2));
while coords(ix,1) == A10(1)
    ix = ix - 1; % decrement
end
A11 = coords(ix,:);
A9 = [max_x,max(coords(bln_min,2))];
ix = find(coords(:,1)==A9(:,1)&coords(:,2)==A9(:,2));
while coords(ix,1) == A9(1)
    ix = ix + 1; % decrement
end
A8 = coords(ix,:);
% OTHER POINTS
min_y = min(coords(:,2));
max_y = max(coords(:,2));
A1 = A2; A1(:,2) = min_y;
A6 = A5; A6(:,2) = max_y;
A7 = A8; A7(:,2) = max_y;
A12 = A11; A12(:,2) = min_y;
cds = [A1;A2;A3;A4;A5;A6;A7;A8;A9;A10;A11;A12;];
% POLYLINE
writePLine(fid,cds);
addString(fid,'(setq currobj (ssadd (entlast) currobj))');
setCurrentLayer(fid,'nitride');
end

function prebiasCoupler(fid,coords)
setCurrentLayer(fid,'prebias');
% Draw horizontal rectangle
min_x = min(coords(:,1));
max_x = max(coords(:,1));
bln_min = find(coords(:,1)==min_x);
bln_max = find(coords(:,1)==max_x);
min_y = min(coords([bln_min;bln_max],2));
max_y = max(coords([bln_min;bln_max],2));
writeRectangle(fid,min_x,min_y,max_x,max_y);
addString(fid,'(setq currobj (ssadd (entlast) currobj))');
% Draw vertical rectangle
min_y = min(coords(:,2));
max_y = max(coords(:,2));
bln_min = find(coords(:,2)==min_y);
bln_max = find(coords(:,2)==max_y);
min_x = min(coords([bln_min;bln_max],1));
max_x = max(coords([bln_min;bln_max],1));
writeRectangle(fid,min_x,min_y,max_x,max_y);
addString(fid,'(setq currobj (ssadd (entlast) currobj))');
setCurrentLayer(fid,'nitride');
end

function prebiasWaveguide(fid,coords)
setCurrentLayer(fid,'prebias');
% Draw vertical rectangle
min_y = min(coords(:,2));
max_y = max(coords(:,2));
bln_min = find(coords(:,2)==min_y);
bln_max = find(coords(:,2)==max_y);
min_x = min(coords([bln_min;bln_max],1));
max_x = max(coords([bln_min;bln_max],1));
writeRectangle(fid,min_x,min_y,max_x,max_y);
addString(fid,'(setq currobj (ssadd (entlast) currobj))');
setCurrentLayer(fid,'nitride');
end

function reflectObjects(fid,xc,yc,namedObjects)
for j = 1:numel(namedObjects)
    % mirrorNamedObject(fid,namedObjects{j},[xc,yc],[xc+1,yc]);
    convertNamedObjectToBlock(fid,namedObjects{j},[xc,yc]);
    insertBlock(fid,['blk',namedObjects{j}],[xc,yc],1,1,0);
    insertBlock(fid,['blk',namedObjects{j}],[xc,yc],1,1,180);
end
end

function cl = replaceScanString(cl,str,val)
% cl = replaceScanString(cl,str,val) takes the input cell and replaces the
% scan parameter anywhere it appears.
cl = cellfun(@(x) strrep(x,str,val),cl,'UniformOutput',0);
end

function saveDXF(fid,fn)
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
writeOutput(fid,['(command "_saveas" "DXF" "Version" "2000" "" "',...
    fn,...
    '")']);
end

function s = setCooledSingleNanobeamDefaults(s)
if ~isfield(s,'tetherEdgeWidth')
    s.tetherEdgeWidth = s.tetherCenterWidth;
end
end

function setCurrentLayer(fid,layerName)
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
writeOutput(fid,['(command "_layer" "set" "',layerName,'" "")']);
end

function s = setNodeArrayDefaults(s)
if ~isfield(s,'blankLength')
    s.blankLength = 1000;
end
if ~isfield(s,'aout')
    s.aout = s.a;
end
if ~isfield(s,'tetherEdgeWidth')
    s.tetherEdgeWidth = s.tetherCenterWidth;
end
if strncmp(s.wf,'inherit',length('inherit'))
    s.wf = s.w0;
end
end

function s = setSingleNanobeamDefaults(s)
if ~isnumeric(s.w0)
    s.w0 = 0;
end
if ~isfield(s,'num')
    s.num = 1;
end
if ~isfield(s,'fun')
    s.fun = @(w0,wf,y) w0 + (wf - w0) .*y.^1;
end
if ~isfield(s,'wf')
    s.wf = s.w0;
end
end

function s = setTetherDefaults(s)
if ~isfield(s,'tetherEdgeWidth');
    s.tetherEdgeWidth = s.tetherCenterWidth;
end
if ~isfield(s,'tetherYDisp')
    s.tetherYDisp = 0;
end
if ~isfield(s,'numTethers')
    s.numTethers = 1;
end
if ~isfield(s,'a')
    s.a = 0;
end
if ~isfield(s,'aout')
    s.aout = s.a;
end
if ~isfield(s,'fillet')
    s.fillet = false;
end
end

function [x0,y0] = setOffsets(dev,spacing,overhang)
x0 = dev*spacing;
y0 = overhang;
end

function s = setVeeGrooveDefaults(s)
% determine if there's a back groove
if ~(s.backGroove)||(length(s.vgl) < 2)
    s.vgl(2) = 0;
end
% correct fine exposure length if necessary
for j = 1:2
    if s.fel(j) >= s.vgl(j) - s.lfda*(j==1)
        s.fel(j) = s.vgl(j) - s.lfda*(j==1);
    end
end
end

function coords = tethersWithFillets(s)
j = 1;
npts_fillet = 5;
coords = zeros(s.numTethers * (4 + 2*npts_fillet) - 1,2);
for m = 1:(s.numTethers)
    coords(j,:) = [s.veeGrooveWidth/2-s.fillet, s.prelength + ...
        (m - 1) * s.aout + s.tetherYDisp];
    % specify first fillet
    fillet.p0 = coords(j,:);
    fillet.pf = coords(j,:) + s.fillet * [1, -1];
    fillet.orientation = 'horizontal';
    fillet.num = npts_fillet;
    j = j + 1;
    coords(j:j+4,:) = drawFillet(fillet);
    j = j + 5;
    coords(j,:) = [s.veeGrooveWidth/2, s.prelength + ...
        s.tetherEdgeWidth + s.tetherYDisp + (m - 1) * s.aout + s.fillet];
    % specify second fillet
    fillet.p0 = coords(j,:);
    fillet.pf = coords(j,:) + s.fillet * [-1, -1];
    fillet.orientation = 'vertical';
    fillet.num = npts_fillet;
    j = j + 1;
    coords(j:j+4,:) = drawFillet(fillet);
    j = j + 5;
    coords(j,:) = [s.nominalWidth/2, s.prelength + ...
        s.tetherCenterWidth + (m - 1)*s.a];
    j = j + 1;
    if (m < s.numTethers)
        coords(j,:) = coords(j - 1,:) + [0,s.a];
        j = j + 1;
    end
end
end

function coords = tethersWithoutFillets(s)
coords = zeros(s.numTethers*4 - 1,2);
for m = 1:(s.numTethers)
    j = (m - 1) * 4 + 1;
    coords(j:j+2,:) = [
        s.veeGrooveWidth/2, s.prelength + (m - 1) * s.aout + s.tetherYDisp;
        s.veeGrooveWidth/2, s.prelength + s.tetherEdgeWidth + ...
        (m - 1) * s.aout + s.tetherYDisp;
        s.nominalWidth/2, s.prelength + s.tetherCenterWidth + (m - 1) * s.a
        ];
    if m < s.numTethers,  %don't add 'a' length after last tether
        coords(j+3,:) = [s.nominalWidth/2 , s.prelength + m * s.a];
    end
end
end

function writeCircle(fid,x0,y0,r)
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
writeOutput(fid,['(command "_circle" "',num2str(x0),',',num2str(y0),...
    '" "',num2str(r),'")']);
end

function writeMarker(fid,x0,y0,size,polarity)
if polarity == 1 % "positive"
    writeRectangle(fid,x0-size/2,y0-size/2,x0+size/2,y0+size/2);
else % "negative"
    m = size/2;
    x = [-8*m, m, m, -m, -m, m, m, 8*m, 8*m, -8*m];
    x = x + x0;
    y = [-8*m, -8*m, -m, -m, m, m, -8*m, -8*m, 8*m, 8*m];
    y = y + y0;
    writePLine(fid,[x;y]');
end
end

function [xf,yf,wf] = writePLine(fid,coords)
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
coords = round(coords,6); % round to nearest pm
list = '(list';
[l,~] = size(coords);
for j = 1:l
    list = [list,' \n(list ',num2str(coords(j,1)),' ',...
        num2str(coords(j,2)),')'];
end
list = [list,')'];
writeOutput(fid,sprintf(['(setq coords ',list,')']));
writeOutput(fid,'(apply ''command (append ''("pline") coords ''("c")))');
indices = find(coords(:,2)==max(coords(:,2)));
xf = coords(indices(1),1);
yf = coords(indices(1),2);
wf = max(pdist(coords(indices,:)));
end

function writeRectangle(fid,x0,y0,xf,yf)
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
writeOutput(fid,['(command "_rectangle" "',num2str(x0),',',num2str(y0),...
    '" "',num2str(xf),',',num2str(yf),'")']);
end

function writeText(fid,x0,y0,scl,rot,str)
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
coords = ['"',num2str(x0),',',num2str(y0),'" '];
hgt = ['"',num2str(scl),'" '];
rot = ['"',num2str(rot),'" '];
str = ['"',str,'" '];
setCurrentLayer(fid,'text');
writeOutput(fid,['(command-s "_text" ',coords,hgt,rot,str,')']); % writes string
writeOutput(fid,'(sssetfirst nil (ssget "_L"))'); % selects last object
writeOutput(fid,'(if (not (numberp (vl-string-search "Mac" (getvar "platform")))) (c:txtexp))'); % explodes text into polyines
setCurrentLayer(fid,'nitride');
end

function writeToFile(fid,str)
fprintf(fid,'%s\n',str);
end

function writeToStandardOut(~,str)
if strncmp(str,';',1)
    fprintf('%s\n',str);
end
end

function zoomToRectangle(fid,a,b)
if (fid == 0)
    writeOutput = @(fid,s) writeToStandardOut(fid,s);
else
    writeOutput = @(fid,s) writeToFile(fid,s);
end
writeOutput(fid,['(command "_zoom" (list ',num2str(a(1)),' ',...
    num2str(a(2)),') (list ',num2str(b(1)),' ',num2str(b(2)),'))']);
end

function zoomExtents(fid)
addString(fid,'(command "_zoom" "extents")'); % zoom extents
end