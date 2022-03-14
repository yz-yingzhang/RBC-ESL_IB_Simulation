%
% Driver for 2D Immersed Boundary Method (IB) using the coarse-grained model
%   implement the physical BCs
%

%---------- INITIALIZATION -------------%
% domain parameters
LD   = 2*pi;              % nondimensionalize length scale
Lsqx = 4*pi;              % side length of the box
Lsqy = 2*pi;
Lsqx = Lsqx/LD;
Lsqy = Lsqy/LD;

Ne = 128;                 % number of Eulerian grid points
hx = Lsqx/Ne;             % Eulerian mesh width
hy = Lsqy/Ne;

% Reynolds number
Re = 0.01;

% simulation temporal parameters
dt = 0.02*hx;            % time step size
Nt = round(10/dt);       % number of time steps

% Lagrangian parameters
K      = 20;            % RBC stiffness
Kb     = 5e-4;          % RBC bending constant
Ka     = 185;           % area conservation constant

% tether points -- ESL
Kp = 0;                 % porous slip parameter
KT = 3200;              % tether force
NT = 2*Ne;
dthetaT = Lsq/NT;
lIdx   = (0:(NT-1))';

% construct Eulerian mesh
X = (0:(Ne-1))' * hx;       % Coordinate Axis
Y = (0:(Ne-1))' * hy;
[X,Y] = meshgrid(X,Y);     % X,Y are coordinate matrices
X = X';
Y = Y';
x = reshape(X, Ne*Ne, 1);   % Eulerian x's for computation
y = reshape(Y, Ne*Ne, 1);   % Eulerian y's for computation

% initial condition for fluid flow and force:
u0  = zeros(Ne*Ne,2);
f   = zeros(Ne*Ne,2);
f(:,1) = sin(pi*y);
fB = f;

% ESL
a = 2*pi;
w = 0.5/6/2;
d = 0.5/6/1.3;
Tl1     = [(lIdx*dthetaT) w*sin(a*lIdx*dthetaT)+d];
Tl2     = [(lIdx*dthetaT) (1-(w*sin(a*lIdx*dthetaT)+d))];
% initial IB points (ESL) that are connected to the tether points
XTl10   = Tl1;                                
XTl20   = Tl2;

% biconcave shape of the RBC
Nl = Ne*2;
dtheta = 2*pi/Nl;
lIdx = (0:(Nl-1))';
r   = 2.8;
alpha = 1.38581894;
xi = lIdx*dtheta-pi/2;

Xl0 = [(r*0.5*alpha*(0.207+2.003*sin(xi).*sin(xi)-1.123*sin(xi).^4).*cos(xi)) (r*alpha*sin(xi))];
Xl0 = Xl0./LD;

% 30 degree initial orientation
[theta,r] = cart2pol(Xl0(:,1), Xl0(:,2));
theta = theta+atan(-tan(pi/3));
[xr,yr]=pol2cart(theta,r);
Xl0(:,1) = xr;
Xl0(:,2) = yr;

% shift the position of the RBC
Xl0(:,2) = Xl0(:,2) + abs(c-min(Xl0(:,2))) + 2.5*hy;
Xl0(:,1) = Xl0(:,1) + Tl1(idx1(1),1) + dic*(ic-1);

dXl0 = 0;
for(q = 1:Nl-1)
   dXl0 = dXl0 + sqrt((Xl0(q+1,1)-Xl0(q,1))^2 + (Xl0(q+1,2)-Xl0(q,2))^2); 
end
dXl0 = dXl0 + sqrt((Xl0(end,1)-Xl0(1,1))^2 + (Xl0(end,2)-Xl0(1,2))^2);
dtheta = dXl0/Nl;

% target area
TA = polyarea(Xl0(:,1), Xl0(:,2));

%---------- END OF INITIALIZATION -------------%

% declare fluid variables for simulation
uX = zeros(Ne*Ne,Nt+1);
uY = uX;
u  = u0;
uh = u0;
uX(:,1) = u0(:,1);
uY(:,1) = u0(:,2);
p  = zeros(Ne*Ne,Nt+1); 
ph = p(:,1);

% Lagrangian points:
Xlx(:,1) = Xl0(:,1);
Xly(:,1) = Xl0(:,2);
Xl       = Xl0;
FT1 = zeros(NT,2);
FT2 = zeros(NT,2);
Up1 = FT1;
Up2 = FT2;
Up12 = Up1;
Up22 = Up2;

%---- operators for the fluid solver ----%
[D0x, D0y]   = D02DDirPer(Ne, hx, hy);
L            = Lu2DDirPer(Ne, hx, hy);

[G0x, G0y] = D02DNeumann(Ne, hx, hy);
Lp         = L2DNeumannPer(Ne, hx, hy);

IL1 = (Re/dt)*speye(size(L)) - (1/2) * L;
IL2 = (Re/dt)*speye(size(L)) + (1/2) * L;
ILp = (Re/dt)*speye(size(Lp)) - (1/2) * Lp;


% plot initial position of the material
figure(1)
plot(mod(Xl(:,1), Lsqx), Xl(:,2),'ko');
hold on;
quiver(reshape(X,Ne,Ne), reshape(Y,Ne,Ne), reshape(u(:,1),Ne,Ne), reshape(u(:,2),Ne,Ne) );
hold off;
axis([0 Lsqx 0 Lsqy]);
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
set(gca,'fontsize',16)
pause

% start time marching
for( i = 1:Nt )
    
    % print time
    fprintf(1, ['i = ' num2str(i) '\n'] );
    
    %---------- Substep 1 ----------%  
    % compute U^n
    u1 = u(:,1);
    u2 = u(:,2);
    
    % ESL
    [Por_Mat,nX,nY] = getPorousSlipV(dthetaT, XTl1(:,1), XTl1(:,2), FT1, Kp);
    Up1(:,1) = -(Por_Mat(:,1)+Por_Mat(:,2)).*nX;
    Up1(:,2) = -(Por_Mat(:,1)+Por_Mat(:,2)).*nY;
    
    [Por_Mat,nX,nY] = getPorousSlipV(dthetaT, XTl2(:,1), XTl2(:,2), FT2, Kp);
    Up2(:,1) = -(Por_Mat(:,1)+Por_Mat(:,2)).*nX;
    Up2(:,2) = -(Por_Mat(:,1)+Por_Mat(:,2)).*nY;

    [idxs, delta] = evalDeltaPhysBCs( XTl1, Ne, hx, hy );
    UTn1 = (hx*hy) * [dot(u1(idxs), delta, 2) dot(u2(idxs), delta, 2)];
    
    [idxs, delta] = evalDeltaPhysBCs( XTl2, Ne, hx, hy );
    UTn2 = (hx*hy) * [dot(u1(idxs), delta, 2) dot(u2(idxs), delta, 2)];
    
    % RBC
    Xltmp = Xl;
    Xltmp(:,1) = mod(Xl(:,1), Lsq);
    [idxs, delta] = evalDeltaPhysBCs( Xltmp, Ne, hx, hy );
    Un = (hx*hy) * [dot(u1(idxs), delta, 2) dot(u2(idxs), delta, 2)];
    
    % march Lagrangian Markers by forward Euler
    XTl12 = XTl1 + dt*(UTn1 + Up1);
    XTl22 = XTl2 + dt*(UTn2 + Up2);
    Xl2 = Xl + dt*Un;
    
    % update Lagrangian force using Xl2
    FT1 = getLagForceTether( XTl12, Tl1, KT );
    FT2 = getLagForceTether( XTl22, Tl2, KT );
    DA = polyarea(Xl2(:,1), Xl2(:,2)) - TA;
    F = getLagForceRLA( Xl2, Nl, dtheta, K, Kb, Ka, DA );
    
    % spread Lagrangian force onto Eulerian grid
    fT1 = spreadLagForcePhysBCs( FT1, XTl12, XTl12, dthetaT, Ne, hx, hy );
    fT2 = spreadLagForcePhysBCs( FT2, XTl22, XTl22, dthetaT, Ne, hx, hy );
    fIB = spreadLagForcePhysBCs( F, Xl2, dtheta, Ne, hx, hy );
    ftot = fB + fIB + fT1 + fT2;
    
    %---------- Substep 2 ----------%
    % solve the non-dimensionalized unsteady stokes equation for t + dt
    [uT, pT] = solveNonDimenNSeqn( u, ph, f, ftot, D0x, D0y, G0x, G0y,...
            Lp, IL1, IL2, ILp );
    
    %---------- Substep 3 ----------%
    % compute UT using updated u and Xl2    
    % ESL
    [Por_Mat,nX,nY] = getPorousSlipV(dthetaT, XTl12(:,1), XTl12(:,2), FT1, Kp);
    Up12(:,1) = -(Por_Mat(:,1)+Por_Mat(:,2)).*nX;
    Up12(:,2) = -(Por_Mat(:,1)+Por_Mat(:,2)).*nY;
    
    [Por_Mat,nX,nY] = getPorousSlipV(dthetaT, XTl22(:,1), XTl22(:,2), FT2, Kp);
    Up22(:,1) = -(Por_Mat(:,1)+Por_Mat(:,2)).*nX;
    Up22(:,2) = -(Por_Mat(:,1)+Por_Mat(:,2)).*nY;
    
    [idxs, delta] = evalDeltaPhysBCs( XTl12, Ne, hx, hy );
    UTT1 = (hx*hy) * [dot(u1(idxs), delta, 2) dot(u2(idxs), delta, 2)];
    
    [idxs, delta] = evalDeltaPhysBCs( XTl22, Ne, hx, hy );
    UTT2 = (hx*hy) * [dot(u1(idxs), delta, 2) dot(u2(idxs), delta, 2)];
    
    % RBC
    Xltmp = Xl2;
    Xltmp(:,1) = mod(Xl2(:,1), Lsq);
    [idxs, delta] = evalDeltaPhysBCs( Xltmp, Ne, hx, hy );
    u1 = uT(:,1);
    u2 = uT(:,2);
    UT = (hx*hy) * [dot(u1(idxs), delta, 2) dot(u2(idxs), delta, 2)];
    
    % march Lagrandian markers using UT and U^n
    XTl1 = XTl1 + dt/2 * (UTn1 + UTT1 + Up1 + Up12);
    XTl2 = XTl2 + dt/2 * (UTn2 + UTT2 + Up2 + Up22);
    Xl = Xl + dt/2 * (Un + UT);
    
    % overwrite for next time step
    f = ftot;
    ph = pT;
    uh = u;
    u = uT;
  
    % plot fluid velocity field and IB points
    figure(1)
    plot(mod(Xl(:,1), Lsqx), Xl(:,2),'ko');
    hold on;
    quiver(reshape(X,Ne,Ne), reshape(Y,Ne,Ne), reshape(u(:,1),Ne,Ne), reshape(u(:,2),Ne,Ne) );
    hold off;
    axis([0 Lsqx 0 Lsqy]);
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$y$', 'Interpreter', 'latex')
    title(['t = ' num2str(dt * i)]);
    set(gca,'fontsize',16)
    
end
