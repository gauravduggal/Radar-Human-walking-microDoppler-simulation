function [segment,seglength,T] = HumanWalkingModel(showplots, formove,...
                              animove, genmovie, Height, rv, nt, numcyc,trajectory,radar_loc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Global Human Walking Model
%
% Based on "A global human walking model with real-time kinematic
% personification", by R. Boulic, N.M. Thalmann, and D. Thalmann
% The Visual Computer, vol.6, pp.344-358, 1990
% This model is based on biomechanical experimental data.
% 
% By V.C. Chen and Yang Hai
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% total number of pulses
np = nt*numcyc;
fprintf('Total Number of Samples = %g\n',np);
maxy = 10;% max(trajectory(2,1:end-1));
miny = 0;%min(trajectory(2,1:end-1));
maxx = 10;%max(trajectory(1,1:end-1));
minx = 0;%min(trajectory(1,1:end-1));

% spatial characteristics
rlc = 1.346*sqrt(rv); % relative length of a cycle
% temporal characteristics
dc = rlc/rv; % duration of a cycle
ds = 0.752*dc-0.143; % duration of support
dsmod = ds/dc; % relative duration of support
T = dc*numcyc; % total time duration 

% body segments' length (meter)
headlen = 0.130*Height;
shoulderlen = (0.259/2)*Height;
torsolen = 0.288*Height;
hiplen = (0.191/2)*Height;
upperleglen = 0.245*Height;
lowerleglen = 0.246*Height;
footlen = 0.143*Height;
upperarmlen = 0.188*Height;
lowerarmlen = 0.152*Height;
Ht = upperleglen + lowerleglen;
sprintf('The Walking Velocity is %g m/s',rv*Ht);

% time scaling
dt = 1/nt;
t = [0:dt:1-dt];
t3 = [-1:dt:2-dt];

% designate gait characteristic
if rv <= 0
    error('Velocity must be positive')
elseif rv < 0.5
    gait = 'a';
elseif rv < 1.3
    gait = 'b';
elseif rv <= 3
    gait = 'c';
else
    error('Relative velocity must be less than 3')
end

% Locations of body segments: 
%      3 translation trajectory coords. give the body segments location
%      relative to rv

% calculate vertical translation: offset from the current height (Hs) of 
% the origin of the spine Os
av = 0.015*rv;
verttrans = -av+av*sin(2*pi*(2*t-0.35));
maxvl = max(verttrans);
minvl = min(verttrans);
diffvl = maxvl-minvl;
if showplots == 'y'
    figure
    plot(t,verttrans)
    title('Vertical translation offset from the origin of the spine')
    xlabel('gait cycle')
    ylabel([sprintf('Position (m) [V_R = %g]',rv)])
    axis([0 1 minvl-0.2*diffvl maxvl+0.2*diffvl])
    drawnow
end

% calculate lateral translation: Os oscillates laterally to ensure the
% weight transfer from one leg to the other.
if gait == 'a'
    al = -0.128*rv^2+0.128*rv;
else
    al = -0.032;
end
lattrans = al*sin(2*pi*(t-0.1));
maxvl = max(lattrans);
minvl = min(lattrans);
diffvl = maxvl-minvl;
if showplots == 'y'
    figure
    plot(t,lattrans)
    title('Lateral translation offset from the origin of the spine')
    xlabel('gait cycle')
    ylabel([sprintf('Offset Position (m) [V_R = %g]',rv)])
    axis([0 1 minvl-0.2*diffvl maxvl+0.2*diffvl])
    drawnow
end

% calculate translation forward/backward: acceleration and deceleration
% phases. When rv grows this effect decreases.
if gait == 'a'
    aa = -0.084*rv^2+0.084*rv;
else
    aa = -0.021;
end
phia = 0.625-dsmod;
transforback = aa*sin(2*pi*(2*t+2*phia));
maxvl = max(transforback);
minvl = min(transforback);
diffvl = maxvl-minvl;
if showplots == 'y'
    figure
    plot(t,transforback)
    title('Translation forward/backward from the origin of the spine')
    xlabel('gait cycle')
    ylabel([sprintf('Offset Position (m) [V_R = %g]',rv)])
    axis([0 1 minvl-0.2*diffvl maxvl+0.2*diffvl])
    drawnow
end

% two rotations of the pelvis

% calculate rotation forward/backward: to make forward motion of the leg,
% the center of gravity of the body must move. To do this, flexing movement
% of the back relatively to the pelvis must be done.
if gait == 'a'
    a1 = -8*rv^2+8*rv;
else
    a1 = 2;
end
rotforback = -a1+a1*sin(2*pi*(2*t-0.1));
maxvl = max(rotforback);
minvl = min(rotforback);
diffvl = maxvl-minvl;
if showplots == 'y'
    figure
    plot(t,rotforback)
    title('Rotation Forward/Backward')
    xlabel('gait cycle')
    ylabel([sprintf('Rotation forward/backward (degree) [V_R = %g]',rv)])
    axis([0 1 minvl-0.2*diffvl maxvl+0.2*diffvl])
    drawnow
end

% calculate rotation left/right: the pelvis falls on th side of the
% swinging leg.
a2 = 1.66*rv;
temp1 = -a2+a2*cos(2*pi*(10*t(1:round(nt*0.15))/3));
temp2 = -a2-a2*cos(2*pi*(10*(t(round(nt*0.15)+1:round(nt*0.50))-0.15)/7));
temp3 =  a2-a2*cos(2*pi*(10*(t(round(nt*0.50)+1:round(nt*0.65))-0.5)/3));
temp4 =  a2+a2*cos(2*pi*(10*(t(round(nt*0.65)+1:nt)-0.65)/7));
rotleftright = [temp1,temp2,temp3,temp4];
maxvl = max(rotleftright);
minvl = min(rotleftright);
diffvl = maxvl-minvl;
if showplots == 'y'
    figure
    plot(t,rotleftright)
    title('Rotation Left/Right')
    xlabel('gait cycle')
    ylabel([sprintf('Rotation left/right (degree) [V_R = %g]',rv)])
    axis([0 1 minvl-0.2*diffvl maxvl+0.2*diffvl])
    drawnow
end

% calculate torsion rotation: pelvis rotates relatively to the spine to
% perform the step
a3 = 4*rv;
torrot = -a3*cos(2*pi*t);
maxvl = max(torrot);
minvl = min(torrot);
diffvl = maxvl-minvl;
if showplots == 'y'
    figure
    plot(t,torrot)
    title('Torsion Rotation')
    xlabel('gait cycle')
    ylabel([sprintf('Torsion angle (degree) [V_R = %g]',rv)])
    axis([0 1 minvl-0.2*diffvl maxvl+0.2*diffvl])
    drawnow
end

% leg flexing/extension: at the hip, at the knee, and at the ankle.

% calculate flexing at the hip
if gait == 'a'
    x1 = -0.1;
    x2 = 0.5;
    x3 = 0.9;
    y1 = 50*rv;
    y2 = -30*rv;
    y3 = 50*rv;
end
if gait == 'b'
    x1 = -0.1;
    x2 = 0.5;
    x3 = 0.9;
    y1 = 25;
    y2 = -15;
    y3 = 25;
end
if gait == 'c'
    x1 = 0.2*(rv-1.3)/1.7-0.1;
    x2 = 0.5;
    x3 = 0.9;
    y1 = 5*(rv-1.3)/1.7+25;
    y2 = -15;
    y3 = 6*(rv-1.3)/1.7+25;
end
if x1+1 == x3
    x = [x1-1,x2-1,x1,x2,x3,x2+1,x3+1];
    y = [y1,y2,y1,y2,y3,y2,y3];
else
    x = [x1-1,x2-1,x3-1,x1,x2,x3,x1+1,x2+1,x3+1];
    y = [y1,y2,y3,y1,y2,y3,y1,y2,y3];
end
temp = pchip(x,y,t3);
flexhip = temp(nt+1:2*nt);
maxvl = max(flexhip);
minvl = min(flexhip);
diffvl = maxvl-minvl;
if showplots == 'y'
    figure
    plot(t,flexhip)
    hold on
    plot(x1,y1,'ro',x2,y2,'ro',x3,y3,'ro')
    title('Flexing at the Hip')
    xlabel('gait cycle')
    ylabel([sprintf('Flexing angle (degree) [V_R = %g]',rv)])
    if x1 < 0
        axis([x1 1 minvl-0.2*diffvl maxvl+0.2*diffvl])
    else
        axis([0 1 minvl-0.2*diffvl maxvl+0.2*diffvl])
    end
    drawnow
end

% calculate flexing at the knee: there are 4 control points.
if gait == 'a'
    x1 = 0.17;
    x2 = 0.4;
    x3 = 0.75;
    x4 = 1;
    y1 = 3;
    y2 = 3;
    y3 = 140*rv;
    y4 = 3;
end
if gait == 'b'
    x1 = 0.17;
    x2 = 0.4;
    x3 = 0.75;
    x4 = 1;
    y1 = 3;
    y2 = 3;
    y3 = 70;
    y4 = 3;
end
if gait == 'c'
    x1 = -0.05*(rv-1.3)/1.7+0.17;
    x2 = 0.4;
    x3 = -0.05*(rv-1.3)/1.7+0.75;
    x4 = -0.03*(rv-1.3)/1.7+1;
    y1 = 22*(rv-1.3)/1.7+3;
    y2 = 3;
    y3 = -5*(rv-1.3)/1.7+70;
    y4 = 3*(rv-1.3)/1.7+3;
end
x = [x1-1,x2-1,x3-1,x4-1,x1,x2,x3,x4,x1+1,x2+1,x3+1,x4+1];
y = [y1,y2,y3,y4,y1,y2,y3,y4,y1,y2,y3,y4];

temp = pchip(x,y,t3);
flexknee = temp(nt+1:2*nt);
maxvl = max(flexknee);
minvl = min(flexknee);
diffvl = maxvl-minvl;
if showplots == 'y'
    figure
    plot(t,flexknee)
    hold on
    plot(x1,y1,'ro',x2,y2,'ro',x3,y3,'ro',x4,y4,'ro')    
    title('Flexing at the Knee')
    xlabel('gait cycle')
    ylabel([sprintf('Flexing angle (degree) [V_R = %g]',rv)])
    if x1 < 0
        axis([x1 1 minvl-0.2*diffvl maxvl+0.2*diffvl])
    else
        axis([0 1 minvl-0.2*diffvl maxvl+0.2*diffvl])
    end
    drawnow
end

% calculate flexing at the ankle: there are 5 control points
if gait == 'a'
    x1 = 0;
    x2 = 0.08;
    x3 = 0.5;
    x4 = dsmod;
    x5 = 0.85;
    y1 = -3;
    y2 = -30*rv-3;
    y3 = 22*rv-3;
    y4 = -34*rv-3;
    y5 = -3;
    
end
if gait == 'b'
    x1 = 0;
    x2 = 0.08;
    x3 = 0.5;
    x4 = dsmod;
    x5 = 0.85;
    y1 = -3;
    y2 = -18;
    y3 = 8;
    y4 = -20;
    y5 = -3;
end
if gait == 'c'
    x1 = 0;
    x2 = 0.08;
    x3 = -0.1*(rv-1.3)/1.7+0.5;
    x4 = dsmod;
    x5 = 0.85;
    y1 = 5*(rv-1.3)/1.7-3;
    y2 = 4*(rv-1.3)/1.7-18;
    y3 = -3*(rv-1.3)/1.7+8;
    y4 = -8*(rv-1.3)/1.7-20;
    y5 = 5*(rv-1.3)/1.7-3;
end
x = [x1-1,x2-1,x3-1,x4-1,x5-1,x1,x2,x3,x4,x5,x1+1,x2+1,x3+1,x4+1,x5+1];
y = [y1,y2,y3,y4,y5,y1,y2,y3,y4,y5,y1,y2,y3,y4,y5];
temp = pchip(x,y,t3);
flexankle = temp(nt+1:2*nt);
maxvl = max(flexankle);
minvl = min(flexankle);
diffvl = maxvl-minvl;
if showplots == 'y'
    figure
    plot(t,flexankle)
    hold on
    plot(x1,y1,'ro',x2,y2,'ro',x3,y3,'ro',x4,y4,'ro',x5,y5,'ro')       
    title('Flexing at the Ankle')
    xlabel('gait cycle')
    ylabel([sprintf('Flexing angle (degree) [V_R = %g]',rv)])
    axis([0 1 minvl-0.2*diffvl maxvl+0.2*diffvl])
    drawnow
end

% trajectory of upper body

% calculate motion (torsion) of the thorax: there are 4 control points
x1 = 0.1;
x2 = 0.4;
x3 = 0.6;
x4 = 0.9;
y1 = (4/3)*rv;
y2 = (-4.5/3)*rv;
y3 = (-4/3)*rv;
y4 = (4.5/3)*rv;
x = [x1-1,x2-1,x3-1,x4-1,x1,x2,x3,x4,x1+1,x2+1,x3+1,x4+1];
y = [y1,y2,y3,y4,y1,y2,y3,y4,y1,y2,y3,y4];
temp = pchip(x,y,t3);
motthor = temp(nt+1:2*nt);
maxvl = max(motthor);
minvl = min(motthor);
diffvl = maxvl-minvl;
if showplots == 'y'
    figure
    plot(t,motthor)
    hold on
    plot(x1,y1,'ro',x2,y2,'ro',x3,y3,'ro',x4,y4,'ro')    
    title('Torsion of the Thorax')
    xlabel('gait cycle')
    ylabel([sprintf('Torsion angle (degree) [V_R = %g]',rv)])
    axis([0 1 minvl-0.2*diffvl maxvl+0.2*diffvl])
    drawnow
end

% calculate flexing at the shoulder
as = 9.88*rv;
flexshoulder = 3-as/2-as*cos(2*pi*t);
maxvl = max(flexshoulder);
minvl = min(flexshoulder);
diffvl = maxvl-minvl;
if showplots == 'y'
    figure
    plot(t,flexshoulder)
    title('Flexing at the Shoulder')
    xlabel('gait cycle')
    ylabel([sprintf('Flexing angle (degree) [V_R = %g]',rv)])
    axis([0 1 minvl-0.2*diffvl maxvl+0.2*diffvl])
    drawnow
end

% calculate flexing at the elbow
if gait == 'a'
    x1 = 0.05;
    x2 = 0.5;
    x3 = 0.9;
    y1 = 6*rv+3;
    y2 = 34*rv+3;
    y3 = 10*rv+3;
end
if gait == 'b'
    x1 = 0.05;
    x2 = 0.01*(rv-0.5)/0.8+0.5;
    x3 = 0.9;
    y1 = 8*(rv-0.5)/0.8+6;
    y2 = 24*(rv-0.5)/0.8+20;
    y3 = 9*(rv-0.5)/0.8+8;
end
if gait == 'c'
    x1 = 0.05;
    x2 = 0.04*(rv-1.3)/1.7+0.51;
    x3 = -0.1*(rv-1.3)/1.7+0.9;
    y1 = -6*(rv-1.3)/1.7+14;
    y2 = 26*(rv-1.3)/1.7+44;
    y3 = -6*(rv-.13)/1.7+17;
end
x = [x1-1,x2-1,x3-1,x1,x2,x3,x1+1,x2+1,x3+1];
y = [y1,y2,y3,y1,y2,y3,y1,y2,y3];
temp = pchip(x,y,t3);
flexelbow = temp(nt+1:2*nt);
maxvl = max(flexelbow);
minvl = min(flexelbow);
diffvl = maxvl-minvl;
if showplots == 'y'
    figure
    plot(t,flexelbow)
    hold on
    plot(x1,y1,'ro',x2,y2,'ro',x3,y3,'ro')    
    title('Flexing at the Elbow')
    xlabel('gait cycle')
    ylabel([sprintf('Flexing angle (degree) [V_R = %g]',rv)])
    axis([0 1 minvl-0.2*diffvl maxvl+0.2*diffvl])
    drawnow
end

% Handling lower body flexing at ankles, knees, and hips
% handle flexing at the left ankle
psi = linspace(0,0,nt);
the(1:round(nt/2)) = flexankle(round(nt/2)+1:nt)*pi/180;
the(round(nt/2)+1:nt) = flexankle(1:round(nt/2))*pi/180;
phi = linspace(0,0,nt);
for i = 1:nt
    Rxyz = XYZConvention(psi(i),the(i),phi(i));
    temp = Rxyz*[footlen;0;0];
    ltoe(1,i) = temp(1);
    ltoe(2,i) = temp(2)+hiplen;
    ltoe(3,i) = temp(3)-(upperleglen+lowerleglen);
end

% handle flexing at the right ankle
psi = linspace(0,0,nt);
the = flexankle*pi/180;
phi = linspace(0,0,nt);
for i = 1:nt
    Rxyz = XYZConvention(psi(i),the(i),phi(i));
    temp = Rxyz*[footlen;0;0];
    rtoe(1,i) = temp(1);
    rtoe(2,i) = temp(2)-hiplen;
    rtoe(3,i) = temp(3)-(upperleglen+lowerleglen);
end

% handle flexing at the left knee
psi = linspace(0,0,nt);
the(1:round(nt/2)) = flexknee(round(nt/2)+1:nt)*pi/(-180);
the(round(nt/2)+1:nt) = flexknee(1:round(nt/2))*pi/(-180);
phi = linspace(0,0,nt);
for i = 1:nt
    Rxyz = XYZConvention(psi(i),the(i),phi(i));
    temp = Rxyz*[0;0;-lowerleglen];
    lankle(1,i) = temp(1);
    lankle(2,i) = temp(2)+hiplen;
    lankle(3,i) = temp(3)-upperleglen;
    temp = Rxyz*[ltoe(1,i);0;ltoe(3,i)+upperleglen];
    ltoe(1,i) = temp(1);
    ltoe(2,i) = temp(2)+hiplen;
    ltoe(3,i) = temp(3)-upperleglen;
end

% handle flexing at the right knee
psi = linspace(0,0,nt);
the = flexknee*pi/(-180);
phi = linspace(0,0,nt);
for i = 1:nt
    Rxyz = XYZConvention(psi(i),the(i),phi(i));
    temp = Rxyz*[0;0;-lowerleglen];
    rankle(1,i) = temp(1);
    rankle(2,i) = temp(2)-hiplen;
    rankle(3,i) = temp(3)-upperleglen;
    temp = Rxyz*[rtoe(1,i);0;rtoe(3,i)+upperleglen];
    rtoe(1,i) = temp(1);
    rtoe(2,i) = temp(2)-hiplen;
    rtoe(3,i) = temp(3)-upperleglen;
end

% handle flexing at the left hip
psi = linspace(0,0,nt);
the(1:round(nt/2)) = flexhip(round(nt/2)+1:nt)*pi/180;
the(round(nt/2)+1:nt) = flexhip(1:round(nt/2))*pi/180;
phi = linspace(0,0,nt);
for i = 1:nt
    Rxyz = XYZConvention(psi(i),the(i),phi(i));
    temp = Rxyz*[0;0;-upperleglen];
    lknee(1,i) = temp(1);
    lknee(2,i) = temp(2)+hiplen;
    lknee(3,i) = temp(3);
    temp = Rxyz*[lankle(1,i);0;lankle(3,i)];
    lankle(1,i) = temp(1);
    lankle(2,i) = temp(2)+hiplen;
    lankle(3,i) = temp(3);
    temp = Rxyz*[ltoe(1,i);0;ltoe(3,i)];
    ltoe(1,i) = temp(1);
    ltoe(2,i) = temp(2)+hiplen;
    ltoe(3,i) = temp(3);
end

% handle flexing at the right hip
psi = linspace(0,0,nt);
the = flexhip*pi/180;
phi = linspace(0,0,nt);
for i = 1:nt
    Rxyz = XYZConvention(psi(i),the(i),phi(i));
    temp = Rxyz*[0;0;-upperleglen];
    rknee(1,i) = temp(1);
    rknee(2,i) = temp(2)-hiplen;
    rknee(3,i) = temp(3);
    temp = Rxyz*[rankle(1,i);0;rankle(3,i)];
    rankle(1,i) = temp(1);
    rankle(2,i) = temp(2)-hiplen;
    rankle(3,i) = temp(3);
    temp = Rxyz*[rtoe(1,i);0;rtoe(3,i)];
    rtoe(1,i) = temp(1);
    rtoe(2,i) = temp(2)-hiplen;
    rtoe(3,i) = temp(3);
end

% Handling lower body rotation
psi = rotleftright*pi/(-180);
the = linspace(0,0,nt);
phi = torrot*pi/180;
for i = 1:nt
    Rxyz = XYZConvention(psi(i),the(i),phi(i));
    temp = Rxyz*[0;hiplen;0];
    lhip(1,i) = temp(1);
    lhip(2,i) = temp(2);
    lhip(3,i) = temp(3);
    temp = Rxyz*[0;-hiplen;0];
    rhip(1,i) = temp(1);
    rhip(2,i) = temp(2);
    rhip(3,i) = temp(3);
    temp = Rxyz*[lknee(1,i);lknee(2,i);lknee(3,i)];
    lknee(1,i) = temp(1);
    lknee(2,i) = temp(2);
    lknee(3,i) = temp(3);
    temp = Rxyz*[rknee(1,i);rknee(2,i);rknee(3,i)];
    rknee(1,i) = temp(1);
    rknee(2,i) = temp(2);
    rknee(3,i) = temp(3);
    temp = Rxyz*[lankle(1,i);lankle(2,i);lankle(3,i)];
    lankle(1,i) = temp(1);
    lankle(2,i) = temp(2);
    lankle(3,i) = temp(3);
    temp = Rxyz*[rankle(1,i);rankle(2,i);rankle(3,i)];
    rankle(1,i) = temp(1);
    rankle(2,i) = temp(2);
    rankle(3,i) = temp(3);
    temp = Rxyz*[ltoe(1,i);ltoe(2,i);ltoe(3,i)];
    ltoe(1,i) = temp(1);
    ltoe(2,i) = temp(2);
    ltoe(3,i) = temp(3);
    temp = Rxyz*[rtoe(1,i);rtoe(2,i);rtoe(3,i)];
    rtoe(1,i) = temp(1);
    rtoe(2,i) = temp(2);
    rtoe(3,i) = temp(3);
end

% Handling upper body flexing at elbows and shoulders
% handling flexing at the left elbow
psi = linspace(0,0,nt);
the(1:round(nt/2)) = flexelbow(round(nt/2)+1:nt)*pi/180;
the(round(nt/2)+1:nt) = flexelbow(1:round(nt/2))*pi/180;
phi = linspace(0,0,nt);
for i = 1:nt
    Rxyz =XYZConvention(psi(i),the(i),phi(i));
    temp = Rxyz*[0;0;-lowerarmlen];
    lhand(1,i) = temp(1);
    lhand(2,i) = temp(2)+shoulderlen;
    lhand(3,i) = temp(3)+(torsolen-upperarmlen);
end

% handle flexing at the right elbow
psi = linspace(0,0,nt);
the = flexelbow*pi/180;
phi = linspace(0,0,nt);
for i = 1:nt
    Rxyz =XYZConvention(psi(i),the(i),phi(i));
    temp = Rxyz*[0;0;-lowerarmlen];
    rhand(1,i) = temp(1);
    rhand(2,i) = temp(2)-shoulderlen;
    rhand(3,i) = temp(3)+(torsolen-upperarmlen);
end

% handle flexing at the left shoulder
psi = linspace(0,0,nt);
the(1:round(nt/2)) = flexshoulder(round(nt/2)+1:nt)*pi/180;
the(round(nt/2)+1:nt) = flexshoulder(1:round(nt/2))*pi/180;
phi = linspace(0,0,nt);
for i = 1:nt
    Rxyz =XYZConvention(psi(i),the(i),phi(i));
    temp = Rxyz*[0;0;-upperarmlen];
    lelbow(1,i) = temp(1);
    lelbow(2,i) = temp(2)+shoulderlen;
    lelbow(3,i) = temp(3)+torsolen;
    temp = Rxyz*[lhand(1,i);0;lhand(3,i)-torsolen];
    lhand(1,i) = temp(1);
    lhand(2,i) = temp(2)+shoulderlen;
    lhand(3,i) = temp(3)+torsolen;
end

% handle flexing at the right shoulder
psi = linspace(0,0,nt);
the = flexshoulder*pi/180;
phi = linspace(0,0,nt);
for i = 1:nt
    Rxyz =XYZConvention(psi(i),the(i),phi(i));
    temp = Rxyz*[0;0;-upperarmlen];
    relbow(1,i) = temp(1);
    relbow(2,i) = temp(2)-shoulderlen;
    relbow(3,i) = temp(3)+torsolen;
    temp = Rxyz*[rhand(1,i);0;rhand(3,i)-torsolen];
    rhand(1,i) = temp(1);
    rhand(2,i) = temp(2)-shoulderlen;
    rhand(3,i) = temp(3)+torsolen;
end

% Handling upper body rotation
psi = linspace(0,0,nt);
the = rotforback*pi/180;
phi = motthor*pi/180;
for i = 1:nt
    Rxyz =XYZConvention(psi(i),the(i),phi(i));
    temp = Rxyz*[0;0;torsolen+headlen];
    head(1,i) = temp(1);
    head(2,i) = temp(2);
    head(3,i) = temp(3);
    temp = Rxyz*[0;0;torsolen];
    neck(1,i) = temp(1);
    neck(2,i) = temp(2);
    neck(3,i) = temp(3);
    temp = Rxyz*[0;shoulderlen;torsolen];
    lshoulder(1,i) = temp(1);
    lshoulder(2,i) = temp(2);
    lshoulder(3,i) = temp(3);
    temp = Rxyz*[0;-shoulderlen;torsolen];
    rshoulder(1,i) = temp(1);
    rshoulder(2,i) = temp(2);
    rshoulder(3,i) = temp(3);
    temp = Rxyz*[lelbow(1,i);lelbow(2,i);lelbow(3,1)];
    lelbow(1,i) = temp(1);
    lelbow(2,i) = temp(2);
    lelbow(3,i) = temp(3);
    temp = Rxyz*[relbow(1,i);relbow(2,i);relbow(3,1)];
    relbow(1,i) = temp(1);
    relbow(2,i) = temp(2);
    relbow(3,i) = temp(3);
    temp = Rxyz*[lhand(1,i);lhand(2,i);lhand(3,i)];
    lhand(1,i) = temp(1);
    lhand(2,i) = temp(2);
    lhand(3,i) = temp(3);
    temp = Rxyz*[rhand(1,i);rhand(2,i);rhand(3,i)];
    rhand(1,i) = temp(1);
    rhand(2,i) = temp(2);
    rhand(3,i) = temp(3);
end

% The origin of the body coordinate system
base = [linspace(0,0,nt);linspace(0,0,nt);linspace(0,0,nt)];

% Handling translation
% base = base+[transforback;lattrans;verttrans];
% neck = neck+[transforback;lattrans;verttrans];
% head = head+[transforback;lattrans;verttrans];
% lshoulder = lshoulder+[transforback;lattrans;verttrans];
% rshoulder = rshoulder+[transforback;lattrans;verttrans];
% lelbow = lelbow+[transforback;lattrans;verttrans];
% relbow = relbow+[transforback;lattrans;verttrans];
% lhand = lhand+[transforback;lattrans;verttrans];
% rhand = rhand+[transforback;lattrans;verttrans];
% lhip = lhip+[transforback;lattrans;verttrans];
% rhip = rhip+[transforback;lattrans;verttrans];
% lknee = lknee+[transforback;lattrans;verttrans];
% rknee = rknee+[transforback;lattrans;verttrans];
% lankle = lankle+[transforback;lattrans;verttrans];
% rankle = rankle+[transforback;lattrans;verttrans];
% ltoe = ltoe+[transforback;lattrans;verttrans];
% rtoe = rtoe+[transforback;lattrans;verttrans];

% Animation of walking human
% if formove == 'y'
%     base(1,:) = base(1,:)+linspace(0,rlc-rlc/(nt+1),nt);
%     neck(1,:) = neck(1,:)+linspace(0,rlc-rlc/(nt+1),nt);
%     head(1,:) = head(1,:)+linspace(0,rlc-rlc/(nt+1),nt);
%     lshoulder(1,:) = lshoulder(1,:)+linspace(0,rlc-rlc/(nt+1),nt);
%     rshoulder(1,:) = rshoulder(1,:)+linspace(0,rlc-rlc/(nt+1),nt);
%     lelbow(1,:) = lelbow(1,:)+linspace(0,rlc-rlc/(nt+1),nt);
%     relbow(1,:) = relbow(1,:)+linspace(0,rlc-rlc/(nt+1),nt);
%     lhand(1,:) = lhand(1,:)+linspace(0,rlc-rlc/(nt+1),nt);
%     rhand(1,:) = rhand(1,:)+linspace(0,rlc-rlc/(nt+1),nt);
%     lhip(1,:) = lhip(1,:)+linspace(0,rlc-rlc/(nt+1),nt);
%     rhip(1,:) = rhip(1,:)+linspace(0,rlc-rlc/(nt+1),nt);
%     lknee(1,:) = lknee(1,:)+linspace(0,rlc-rlc/(nt+1),nt);
%     rknee(1,:) = rknee(1,:)+linspace(0,rlc-rlc/(nt+1),nt);
%     lankle(1,:) = lankle(1,:)+linspace(0,rlc-rlc/(nt+1),nt);
%     rankle(1,:) = rankle(1,:)+linspace(0,rlc-rlc/(nt+1),nt);
%     ltoe(1,:) = ltoe(1,:)+linspace(0,rlc-rlc/(nt+1),nt);
%     rtoe(1,:) = rtoe(1,:)+linspace(0,rlc-rlc/(nt+1),nt);
% end
neck_d = zeros(3,nt*numcyc);
head_d = zeros(3,nt*numcyc);
base_d = zeros(3,nt*numcyc);
lknee_d = zeros(3,nt*numcyc);
rknee_d = zeros(3,nt*numcyc);
lankle_d = zeros(3,nt*numcyc);
rankle_d = zeros(3,nt*numcyc);
lhip_d =zeros(3,nt*numcyc);
rhip_d =zeros(3,nt*numcyc);
relbow_d = zeros(3,nt*numcyc);
lelbow_d =zeros(3,nt*numcyc);
rhand_d = zeros(3,nt*numcyc);
lhand_d =zeros(3,nt*numcyc);
rshoulder_d =zeros(3,nt*numcyc);
lshoulder_d = zeros(3,nt*numcyc);
rtoe_d = zeros(3,nt*numcyc);
ltoe_d = zeros(3,nt*numcyc);
f = figure;
framenum = 1;
colormap([0.7 0.7 0.7])
f.Position = [0 0 1920 1080];
for i = 1:nt*numcyc-2
  if formove == 't'
    idx = mod(i,nt)+1;  
    ox = trajectory(1,i)+transforback(idx);
    oy = trajectory(2,i)+lattrans(idx);
    phi = -atan2(trajectory(2,i+1)-trajectory(2,i),trajectory(1,i+1)-trajectory(1,i));

    oz = verttrans(idx)+1;
    
%     size(base(:,idx))
    R = XYZConvention(0,0,phi);

    base_d(:,i) = RotateandTranslate(base(:,idx),R,ox,oy,oz);
    neck_d(:,i) = RotateandTranslate(neck(:,idx),R,ox,oy,oz);
    head_d(:,i) = RotateandTranslate(head(:,idx),R,ox,oy,oz);
    lknee_d(:,i) = RotateandTranslate(lknee(:,idx),R,ox,oy,oz);
    rknee_d(:,i) = RotateandTranslate(rknee(:,idx),R,ox,oy,oz);
    lankle_d(:,i) = RotateandTranslate(lankle(:,idx),R,ox,oy,oz);
    rankle_d(:,i) = RotateandTranslate(rankle(:,idx),R,ox,oy,oz);
    lhip_d(:,i) = RotateandTranslate(lhip(:,idx),R,ox,oy,oz);
    rhip_d(:,i) = RotateandTranslate(rhip(:,idx),R,ox,oy,oz);
    relbow_d(:,i) = RotateandTranslate(relbow(:,idx),R,ox,oy,oz);
    lelbow_d(:,i) = RotateandTranslate(lelbow(:,idx),R,ox,oy,oz);
    lhand_d(:,i) = RotateandTranslate(lhand(:,idx),R,ox,oy,oz);
    rhand_d(:,i) = RotateandTranslate(rhand(:,idx),R,ox,oy,oz);
    rshoulder_d(:,i) = RotateandTranslate(rshoulder(:,idx),R,ox,oy,oz);
    lshoulder_d(:,i) = RotateandTranslate(lshoulder(:,idx),R,ox,oy,oz);
    ltoe_d(:,i) = RotateandTranslate(ltoe(:,idx),R,ox,oy,oz);
    rtoe_d(:,i) = RotateandTranslate(rtoe(:,idx),R,ox,oy,oz);
  end
  if animove == 'y' && mod(i,50)==0
    clf
    hold on
    [x,y,z] = cylinder2P(neck_d(:,i)',head_d(:,i)',0.04,10);
    surf(x,y,z) % neck
    [x,y,z] = ...
        ellipsoid2P(base_d(:,i)',neck_d(:,i)',0.15,0.15,torsolen/2,20);
    surf(x,y,z) % torso
    [x,y,z] = ...
        ellipsoid2P(lankle_d(:,i)',lknee_d(:,i)',0.06,0.06,lowerleglen/2,10);
    surf(x,y,z) % left lower leg
    [x,y,z] = ...
        ellipsoid2P(rankle_d(:,i)',rknee_d(:,i)',0.06,0.06,lowerleglen/2,10);
    surf(x,y,z) % right lower leg
    [x,y,z] = ...
        ellipsoid2P(ltoe_d(:,i)',lankle_d(:,i)',0.05,0.05,footlen/2,10);
    surf(x,y,z) % left foot
    [x,y,z] = ...
        ellipsoid2P(rtoe_d(:,i)',rankle_d(:,i)',0.05,0.05,footlen/2,10);
    surf(x,y,z) % right foot
    [x,y,z] = ...
        ellipsoid2P(lknee_d(:,i)',lhip_d(:,i)',0.07,0.07,upperleglen/2,10);
    surf(x,y,z) % left upper leg
    [x,y,z] = ...
        ellipsoid2P(rknee_d(:,i)',rhip_d(:,i)',0.07,0.07,upperleglen/2,10);
    surf(x,y,z) % right upper leg
    [x,y,z] = ...
        ellipsoid2P(base_d(:,i)',lhip_d(:,i)',0.07,0.07,hiplen/2,10);
    surf(x,y,z) % left hip
    [x,y,z] = ...
        ellipsoid2P(base_d(:,i)',rhip_d(:,i)',0.07,0.07,hiplen/2,10);
    surf(x,y,z) % right hip
    [x,y,z] = ...
      ellipsoid2P(lelbow_d(:,i)',lshoulder_d(:,i)',0.06,0.06,upperarmlen/2,10);
    surf(x,y,z) % left upper arm
    [x,y,z] = ...
      ellipsoid2P(relbow_d(:,i)',rshoulder_d(:,i)',0.06,0.06,upperarmlen/2,10);
    surf(x,y,z) % right upper arm         
    [x,y,z] = ...
      ellipsoid2P(neck_d(:,i)',lshoulder_d(:,i)',0.06,0.06, shoulderlen/2,10);
    surf(x,y,z) % left shoulder
    [x,y,z] = ...
      ellipsoid2P(neck_d(:,i)',rshoulder_d(:,i)',0.06,0.06, shoulderlen/2,10);
    surf(x,y,z) % right shoulder          
    [x,y,z] = ...
      ellipsoid2P(lhand_d(:,i)',lelbow_d(:,i)',0.05,0.05,lowerarmlen/2,10);
    surf(x,y,z) % left lower arm
    [x,y,z] = ...
      ellipsoid2P(rhand_d(:,i)',relbow_d(:,i)',0.05,0.05,lowerarmlen/2,10);
    surf(x,y,z) % right lower arm
    [x,y,z] = sphere(10);
    surf(x*0.02+ltoe_d(1,i),y*0.02+ltoe_d(2,i),z*0.01+ltoe_d(3,i)) % left toe
    [x,y,z] = sphere(10);
    surf(x*0.02+rtoe_d(1,i),y*0.02+rtoe_d(2,i),z*0.01+rtoe_d(3,i)) % right toe
    [x,y,z] = sphere(10);
    surf(x*0.05+lankle_d(1,i),y*0.05+lankle_d(2,i),z*0.05+lankle_d(3,i)) % lankle
    [x,y,z] = sphere(10);
    surf(x*0.05+rankle_d(1,i),y*0.05+rankle_d(2,i),z*0.05+rankle_d(3,i)) % rankle
    [x,y,z] = sphere(10);
    surf(x*0.05+lknee_d(1,i),y*0.05+lknee_d(2,i),z*0.05+lknee_d(3,i)) % lknee
    [x,y,z] = sphere(10);
    surf(x*0.05+rknee_d(1,i),y*0.05+rknee_d(2,i),z*0.05+rknee_d(3,i)) % rknee
    [x,y,z] = sphere(10);
    surf(x*0.1+head_d(1,i),y*0.1+head_d(2,i),z*headlen/2+head_d(3,i)) % head
    [x,y,z] = sphere(10);
    surf(x*0.05+lhip_d(1,i),y*0.05+lhip_d(2,i),z*0.05+lhip_d(3,i)) % lhip
    [x,y,z] = sphere(10);
    surf(x*0.05+rhip_d(1,i),y*0.05+rhip_d(2,i),z*0.05+rhip_d(3,i)) % rhip
    [x,y,z] = sphere(10); % left shoulder
    surf(x*0.05+lshoulder_d(1,i),y*0.05+lshoulder_d(2,i),z*0.05+lshoulder_d(3,i))
    [x,y,z] = sphere(10); % right shoulder
    surf(x*0.05+rshoulder_d(1,i),y*0.05+rshoulder_d(2,i),z*0.05+rshoulder_d(3,i))
    [x,y,z] = sphere(10); 
    surf(x*0.05+lelbow_d(1,i),y*0.05+lelbow_d(2,i),z*0.05+lelbow_d(3,i)) % lelbow
    [x,y,z] = sphere(10);
    surf(x*0.05+relbow_d(1,i),y*0.05+relbow_d(2,i),z*0.05+relbow_d(3,i)) % relbow
    [x,y,z] = sphere(10);
    surf(x*0.05+lhand_d(1,i),y*0.05+lhand_d(2,i),z*0.05+lhand_d(3,i)) % lhand
    [x,y,z] = sphere(10);
    surf(x*0.05+rhand_d(1,i),y*0.05+rhand_d(2,i),z*0.05+rhand_d(3,i)) % rhand
    stem3(radar_loc(1),radar_loc(2),radar_loc(3))
    plot(trajectory(1,:),trajectory(2,:))
    light
    lighting gouraud
    shading interp
    axis equal
    if formove == 'y'
        axis([-2,2*numcyc,-2,2,-2,2])
    elseif formove =='t'
        xlim([minx-1 maxx+1])
        ylim([miny-1 maxy+1])
        zlim([0 2])
    else
        zlim([-2 2])
        axis([-2,2,-2,2,-2,2])
    end
%     axis off
%     grid off
    grid on
    axis on
    set(gcf,'Color',[1 1 1])
    view([45,20])
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
%     zoom(2)
    drawnow
    if genmovie == 'y'
        M(framenum) = getframe(gcf);
        framenum = framenum+1;
    end
  end

%   if (mod(i,nt) == 0 && numcyc>1)
%     if formove == 'y'
%             temp = base(:,i-nt+1:i);
%             temp(1,:) = temp(1,:)+linspace(rlc,rlc,nt);
%             base = [base,temp];
%             temp = neck(:,i-nt+1:i);
%             size(temp)
%             temp(1,:) = temp(1,:)+linspace(rlc,rlc,nt);
%             neck = [neck,temp];
%             temp = head(:,i-nt+1:i);
%             temp(1,:) = temp(1,:)+linspace(rlc,rlc,nt);
%             head = [head,temp];
%             temp = lshoulder(:,i-nt+1:i);
%             temp(1,:) = temp(1,:)+linspace(rlc,rlc,nt);
%             lshoulder = [lshoulder,temp];
%             temp = rshoulder(:,i-nt+1:i);
%             temp(1,:) = temp(1,:)+linspace(rlc,rlc,nt);
%             rshoulder = [rshoulder,temp];
%             temp = lelbow(:,i-nt+1:i);
%             temp(1,:) = temp(1,:)+linspace(rlc,rlc,nt);
%             lelbow = [lelbow,temp];
%             temp = relbow(:,i-nt+1:i);
%             temp(1,:) = temp(1,:)+linspace(rlc,rlc,nt);
%             relbow = [relbow,temp];
%             temp = lhand(:,i-nt+1:i);
%             temp(1,:) = temp(1,:)+linspace(rlc,rlc,nt);
%             lhand = [lhand,temp];
%             temp = rhand(:,i-nt+1:i);
%             temp(1,:) = temp(1,:)+linspace(rlc,rlc,nt);
%             rhand = [rhand,temp];
%             temp = lhip(:,i-nt+1:i);
%             temp(1,:) = temp(1,:)+linspace(rlc,rlc,nt);
%             lhip = [lhip,temp];
%             temp = rhip(:,i-nt+1:i);
%             temp(1,:) = temp(1,:)+linspace(rlc,rlc,nt);
%             rhip = [rhip,temp];
%             temp = lknee(:,i-nt+1:i);
%             temp(1,:) = temp(1,:)+linspace(rlc,rlc,nt);
%             lknee = [lknee,temp];
%             temp = rknee(:,i-nt+1:i);
%             temp(1,:) = temp(1,:)+linspace(rlc,rlc,nt);
%             rknee = [rknee,temp];
%             temp = lankle(:,i-nt+1:i);
%             temp(1,:) = temp(1,:)+linspace(rlc,rlc,nt);
%             lankle = [lankle,temp];
%             temp = rankle(:,i-nt+1:i);
%             temp(1,:) = temp(1,:)+linspace(rlc,rlc,nt);
%             rankle = [rankle,temp];
%             temp = ltoe(:,i-nt+1:i);
%             temp(1,:) = temp(1,:)+linspace(rlc,rlc,nt);
%             ltoe = [ltoe,temp];
%             temp = rtoe(:,i-nt+1:i);
%             temp(1,:) = temp(1,:)+linspace(rlc,rlc,nt);
%             rtoe = [rtoe,temp];
% 
%             
%     else
%             temp = base(:,i-nt+1:i);
%             base = [base,temp];
%             temp = neck(:,i-nt+1:i);
%             neck = [neck,temp];
%             temp = head(:,i-nt+1:i);
%             head = [head,temp];
%             temp = lshoulder(:,i-nt+1:i);
%             lshoulder = [lshoulder,temp];
%             temp = rshoulder(:,i-nt+1:i);
%             rshoulder = [rshoulder,temp];
%             temp = lelbow(:,i-nt+1:i);
%             lelbow = [lelbow,temp];
%             temp = relbow(:,i-nt+1:i);
%             relbow = [relbow,temp];
%             temp = lhand(:,i-nt+1:i);
%             lhand = [lhand,temp];
%             temp = rhand(:,i-nt+1:i);
%             rhand = [rhand,temp];
%             temp = lhip(:,i-nt+1:i);
%             lhip = [lhip,temp];
%             temp = rhip(:,i-nt+1:i);
%             rhip = [rhip,temp];
%             temp = lknee(:,i-nt+1:i);
%             lknee = [lknee,temp];
%             temp = rknee(:,i-nt+1:i);
%             rknee = [rknee,temp];
%             temp = lankle(:,i-nt+1:i);
%             lankle = [lankle,temp];
%             temp = rankle(:,i-nt+1:i);
%             rankle = [rankle,temp];
%             temp = ltoe(:,i-nt+1:i);
%             ltoe = [ltoe,temp];
%             temp = rtoe(:,i-nt+1:i);
%             rtoe = [rtoe,temp];
%     end
%   end
end


% generate movie file
if genmovie == 'y'
    v = VideoWriter('HumanWalkingModel.avi');
    open(v)
    writeVideo(v,M)
    close(v)
    %movie2avi(M,'HumanWalkingModel.avi','FPS',20,'compression','None')
end

% output data
seglength(1).name = 'Head Length';
seglength(1).length = headlen;
seglength(2).name = 'Shoulder Length';
seglength(2).length = shoulderlen;
seglength(3).name = 'Torso Length';
seglength(3).length = torsolen;
seglength(4).name = 'Hip Length';
seglength(4).length = hiplen;
seglength(5).name = 'Upper Leg Length';
seglength(5).length = upperleglen;
seglength(6).name = 'Lower Leg Length';
seglength(6).length = lowerleglen;
seglength(7).name = 'Foot Length';
seglength(7).length = footlen;
seglength(8).name = 'Upper Arm Length';
seglength(8).length = upperarmlen;
seglength(9).name = 'Lower Arm Length';
seglength(9).length = lowerarmlen;
segment(1).name = 'Base';
segment(1).PositionData = base_d;
segment(2).name = 'Neck';
segment(2).PositionData = neck_d;
segment(3).name = 'Head';
segment(3).PositionData = head_d;
segment(4).name = 'Left Shoulder';
segment(4).PositionData = lshoulder_d;
segment(5).name = 'Right Shoulder';
segment(5).PositionData = rshoulder_d;
segment(6).name = 'Left Elbow';
segment(6).PositionData = lelbow_d;
segment(7).name = 'Right Elbow';
segment(7).PositionData = relbow_d;
segment(8).name = 'Left Hand';
segment(8).PositionData = lhand_d;
segment(9).name = 'Right Hand';
segment(9).PositionData = rhand_d;
segment(10).name = 'Left Hip';
segment(10).PositionData = lhip_d;
segment(11).name = 'Right Hip';
segment(11).PositionData = rhip_d;
segment(12).name = 'Left Knee';
segment(12).PositionData = lknee_d;
segment(13).name = 'Right Knee';
segment(13).PositionData = rknee_d;
segment(14).name = 'Left Ankle';
segment(14).PositionData = lankle_d;
segment(15).name = 'Right Ankle';
segment(15).PositionData = rankle_d;
segment(16).name = 'Left Toe';
segment(16).PositionData = ltoe_d;
segment(17).name = 'Right Toe';
segment(17).PositionData = rtoe_d;

end
function [rt_point] = RotateandTranslate(point,R,x,y,z)
    rt_point = R*point+[x;y;z];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

