% Create input files for antarctic LES, start left

kmax = 40;
zsize = 100;
radmax = 55;

rowl = max(kmax,radmax);
dz = zsize/kmax;

A = zeros(rowl,14);

z = linspace(0.5*dz, (rowl-0.5)*dz, rowl);
th = zeros(rowl,1);
u = zeros(rowl,1);
v = zeros(rowl,1);
ug = zeros(rowl,1);
wls = zeros(rowl,1);
q_bg = zeros(rowl,1);

c0 = zeros(rowl,1);
c1 = zeros(rowl,1);
c2 = zeros(rowl,1);
pl0 = zeros(rowl,1);
pl1 = zeros(rowl,1);
fdn_above_0 = zeros(rowl,1);
fdn_above_1 = zeros(rowl,1);

for jj=1:rowl
    ug(jj) = 3.0;
    q_bg(jj) = 3.7*10^(-4);
    wls(jj) = -0.01* z(jj) / 100;
    
    if(z(jj) > 60)
       u(jj) = ug(jj);
       v(jj) = 0;
    end
end

B = textread('terms_g4_top100m_microhh');
B = B(:,1:7);

A(:,1) = z;
A(:,2) = th;
A(:,3) = u;
A(:,4) = v;
A(:,5) = ug;
A(:,6) = wls;
A(:,7) = q_bg;

A(:,8:14) = B;

% for jj=1:radmax
% c0(jj) = A(jj,1);
% c1(jj) = A(jj,2);
% c2(jj) = A(jj,3);
% pl0(jj) = A(jj,4);
% pl1(jj) = A(jj,5);
% fdn_above_0(jj) = A(jj,6);
% fdn_above_1(jj) = A(jj,7);
% end

fileName = 'antarctic_les_left_g4.prof';
fID = fopen(fileName,'w');
fprintf(fID,'%16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s\n','z','th','u','v','ug','wls','q_bg','c0','c1','c2','pl0','pl1','fdn_above_0','fdn_above_1');
fprintf(fID,'%1.10E %1.10E %1.10E %1.10E %1.10E %1.10E %1.10E %1.10E %1.10E %1.10E %1.10E %1.10E %1.10E %1.10E\n',A');
fclose(fID);


