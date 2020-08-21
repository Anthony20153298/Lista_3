//Motor DC
//Autor: Vasquez Rivera Anthony A.
//Sistema de Control Avanzado

ap= [-675,-26.25;1050,-0.093]// Matrix A
bp = [250,0;0,250]// Matrix B
cp = [1,0;0,1]// Matrix C
dp = 0*ones(2,2)// Matrix D
Gs= syslin("c",ap, bp, cp, dp)//Linear system definition
[tf]=ss2tf(Gs)//Conversion from state-space to transfer function

// Calculation of poles and zeros of the plant //
plzr(tf);// Pole-zero plot
scf(1);

// Definition of stability and performance barriers //
w=logspace(-3,3,400);
a=200;b=130;c=70;
x1=40*ones(1,a); x2=60*zeros(1,b); x3=-40*ones(1,c);
xt=[x1 x2 x3];
x4=[5*ones(1,100) 0*zeros(1,300)];
scf(2);
plot2d("ln", w, x4,3, rect=[10^-1 -60 10^3 60])
plot2d("ln", w, xt)
xgrid(12)
xtitle("Stability Barriers","Frequency w(rad/s)", "Amplitude (dB)");
s= svd(ap);

// Plot of maximum and minimum singular values ​​//
tri = trzeros(Gs)//Transmission zeros and normal rank
w = logspace(-3,3);
svi = svplot(Gs,w);//Singular-value sigma-plot
scf(3);
plot2d("ln", w, 20*log(svi')/log(10))
xgrid(12)
xtitle("Design Plant Singular Values","Frequency (rad/s)", "Amplitude (dB)");

// We add an integrator to the plant //
[ns,nc]=size(bp); // ns = number of entries; nc = number of controls
Ai=[ap,bp;0*ones(nc,ns),0*ones(nc,nc)]; //Matrix A with integrator
Bi=[0*ones(ns,nc); eye(nc,nc)];//Matrix B with integrator
Ci=[cp 0*ones(nc,nc)];//Matrix C with integrator
Di=0*ones(nc,nc);//Matrix D with integrator
sysi=syslin('c',Ai,Bi,Ci,Di);
I=eye(nc);//Identity matrix

// Calculation and plotting of singular values ​​with integrator //
tri = trzeros(sysi)
w = logspace(-3,3);
svi = svplot(sysi,w);
scf(4);
plot2d("ln", w, 20*log(svi')/log(10))
xgrid(12)
xtitle("Design Plant Singular Values","Frequency (rad/s)", "Amplitude (dB)");

// LQR controller //
C=0.7*Ci'*Ci;        // Weighting matrix Q
rho=1e-1;       // Value of rho
R = rho*eye(nc);    // Associated with the cost matrix C
B=Bi*inv(R)*Bi';    // Recalculating B
A=Ai;

// Riccati in the system //
X=riccati(A,B,C,'c','eigen'); 

// The gain of the controller //
G=inv(R)*Bi'*X;

// Kalman filter layout //
ll= inv(cp*inv(-ap)*bp+dp);      // Matrix L
lh = -inv(ap)*bp*ll;               
Lp=[lh,                        
   ll];      
pnint = eye(nc,nc)              // Values ​​for the duality of the Filter
mu = 0.1;                       // LQR controller, applying Riccati
THETA = mu*eye(nc,nc)  
Ah=Ai'; // Calculation of Ah
Bh=Ci'*inv(THETA)*Ci;          // Calculation of Bh
Ch=Lp*Lp'; //Calculation of Ch
Xh=riccati(Ah,Bh,Ch,'c','eigen'); // Riccati application to the system

// Calculation of the gain of H
H=(inv(THETA)*Ci*Xh)';
sysh = syslin('c',Ai,H,Ci,Di);

// Calculation of singular values ​​of the filter //
trh = trzeros(sysh)
w = logspace(-3,3);
svh = svplot(sysh,w);
scf(5);
plot2d("ln", w, 20*log(svh')/log(10))
xgrid(12)
xtitle("Singular Values – Kalman Filter","Frequency (rad/s)","Amplitude (dB)");

//Compensator//
Ak = [ Ai-Bi*G-H*Ci  0*ones(ns+nc,nc)
       G          0*ones(nc,nc) ]//Matrix A with compensator
Bk = [ H
       0*ones(nc,nc) ]//Matrix B with compensator
Ck = [0*ones(nc, ns+nc) eye(nc,nc) ]//Matrix C with compensator
Dk = 0*ones(nc,nc);//Matrix D with compensator
sysk=syslin('c',Ak,Bk,Ck,Dk);

// Calculation of singular values ​​of the compensator //
trk = trzeros(sysk)
w = logspace(-3,3);
svk = svplot(sysk,w);
scf(6);
plot2d("ln", w, 20*log(svk')/log(10))
xgrid(12)

// We analyze in an open loop //
Abo = [ ap                     bp*Ck
       0*ones(ns+nc+nc,ns)    Ak    ]//Matrix A in open loop
Bbo = [ 0*ones(ns,nc)
       Bk ]
    //Matrix B in open loop
Cbo = [ cp  0*ones(nc,ns+nc+nc) ]//Matrix C in open loop
Dbo = 0*ones(nc,nc);//Matrix D in open loop
sysbo = syslin('c',Abo,Bbo,Cbo,Dbo);

//Singular values ​​of open loop //
vsbo = svplot(sysbo,w);
scf(7)
plot2d("ln", w, 20*log(vsbo')/log(10))
xgrid(12)
xtitle("Singular values plot del bucle abierto","Frequency (rad/s)", "Amplitude (dB)");
xtitle("Compensator Singular Values","Frequency (rad/s)", "Amplitude (dB)");

// Sensitivity analysis of S //
SS= syslin("c",Abo-Bbo*Cbo, Bbo, -Cbo, eye(nc,nc))
ssi = svplot(SS,w);
scf(8)
plot2d("ln", w, 20*log(ssi')/log(10))
xgrid(12)
xtitle("Singular values plot sensibility S","Frequency (rad/s)", "Amplitude (dB)");

// Sensitivity analysis of T //
ST= syslin('c',Abo-Bbo*Cbo, Bbo, Cbo, Dbo)
sti = svplot(ST,w);
scf(9)
plot2d("ln", w, 20*log(sti')/log(10))
xgrid(12)
xtitle("Singular values plot Sensibility T","Frequency (rad/s)", "Amplitude (dB)");

// Analysis of overlapping S and T //
scf(10)
plot2d("ln", w, [20*log(ssi')/log(10)])
plot2d("ln", w, [20*log(sti')/log(10)])
xgrid(12)
xtitle("Singular values plot S y T","Frequency (rad/s)", "Amplitude (dB)");


