clear;
clc;
tc=cputime;

%%%%%%Paramaters%%%%%%
%ep: small number;
%If: flood intermittency;
%timeyear: coefficient convert year to second;
%g: gravitaional coefficient;
%rho: water density;
%R: submerged specific gravity;
%L: channel length;
%SI: initial slope;
%Duration: duration of calculation (year);
%qwf: flow discharge per unit width at the entrance (m2/s);
%qtfT: total sediment supply per unit width (m2/s);
%alpr: coefficient in M-S resistance formula;
%nk: coefficient relating roughness height;
%Phi: phi value of grain size;
%D: characteristic grain size of each size range;
%ng: number of sediment group;
%qgfT: gravel supply rate per unit width (m2/s);
%qsfT: sand supply rate per unit width (m2/s);
%qtfT: total sediment supply rate per unit width (m2/s);
%qtf: sediment supply rate per unit width for each range (m2/s);
%fgravel: grain size distribution of gravel supply;
%fsand: grain size distribution of sand supply;
%fsupply: grain size distribution of total supply;
%fbase: grain size distribution of substrate base;
%na: coefficient relating active layer thickness;
%p: bed porosity;
%au: coefficient related to Exner discretization;
%Ls: thickness of storage layer;
%sigma: subsidence rate (m/year);
ep=1e-6;
If=1;
% time_a_minute
% timeyear=If*60;
timeyear=1;
g=9.8;
rho=1000;
R=1.65;
L=7.5;
SI=0.005;
Duration=300*60;
Nt=60;
qwf=22.8/1000/0.6;
alpr=8.1;
nk=1.9;
nt=0.5;
na=2;
p=0.35;
au=0.75;
Ls=0.02;

Phin=linspace(-2,7,10);
Phi=(Phin(1:9)+Phin(2:10))/2;
D=2.^Phi/1000;
ng=size(D,2);
% fgravel=[<0.5,<1,<2,<4,<8,<16,<32,<64,<128]/100
% % sediment supply
% ¦Ò=2.17
fgravel=[0,0,0,16,42,0,0,0,0]/100;
fsand=[0,0,42,0,0,0,0,0,0]/100;

% ¦Ò=4.23
% fgravel=[0,0,0,8,8,38,0,0,0]/100;
% fsand=[0,38,8,0,0,0,0,0,0]/100;

% ¦Ò=8.24
% fgravel=[0,0,0,4,3,10,35,0,0]/100;
% fsand=[35,10,3,0,0,0,0,0,0]/100;


qgfT=1;
qsfT=1;
qtf=qgfT*fgravel+qsfT*fsand;
qtf_unit=qtf;
qtfT=sum(qtf,2);
fsupply=qtf./qtfT;
% fbase=fsupply;

% bed material
fgravel1=[0,0,0,100,0,0,0,0,0]/100;
fsand1=[0,0,0,0,0,0,0,0,0]/100;
qgfT1=1;
qsfT1=1;
qtf1=qgfT1*fgravel1+qsfT1*fsand1;
qtf_unit1=qtf1;
qtfT1=sum(qtf1,2);
fsupply1=qtf1./qtfT1;
fbase=fsupply1;


%dt: time step (year);
%dth: time step for hydraulic calculation (year);
%nh: dt/dth;
%dx: cell size;
%n: cell number;
dt=0.1;
dth=0.01;
nh=dt/dth;
n=76;
dx=L/(n-1);
Courant=zeros(n,1);
interval=round(Duration/dt/Nt);
num_morpho_steps=round(Duration/dt);

% clear water
% qtfT=zeros(num_morpho_steps,1);
% qtf=zeros(num_morpho_steps,9);

% constant
% qsf_low=120;
% qsf_peak=120;
% peak_time=50; % min
% qsf_low=qsf_low/(1000*2650*(1-p));
% qsf_peak=qsf_peak/(1000*2650*(1-p));
% peak_time_index=round(num_morpho_steps*peak_time/300);
% sediment_feed_end_index=round(num_morpho_steps*100/300);
% qtfT=zeros(num_morpho_steps,1);
% qtfT(1:peak_time_index)=linspace(qsf_low,qsf_peak,peak_time_index);
% qtfT(peak_time_index:sediment_feed_end_index)=linspace(qsf_peak,qsf_low,(sediment_feed_end_index-peak_time_index+1));
% qtf=zeros(num_morpho_steps,9);
% for i=1:9
%     qtf(1:peak_time_index,i)=linspace(qsf_low,qsf_peak,peak_time_index)*qtf_unit(i);
%     qtf(peak_time_index:sediment_feed_end_index,i)=linspace(qsf_peak,qsf_low,(sediment_feed_end_index-peak_time_index+1))*qtf_unit(i);
% end

% r=0.25
qsf_low=45;
qsf_peak=195;
peak_time=25; % min
qsf_low=qsf_low/(1000*2650*(1-p));
qsf_peak=qsf_peak/(1000*2650*(1-p));
peak_time_index=round(num_morpho_steps*peak_time/300);
sediment_feed_end_index=round(num_morpho_steps*100/300);
qtfT=zeros(num_morpho_steps,1);
qtfT(1:peak_time_index)=linspace(qsf_low,qsf_peak,peak_time_index);
qtfT(peak_time_index:sediment_feed_end_index)=linspace(qsf_peak,qsf_low,(sediment_feed_end_index-peak_time_index+1));
qtf=zeros(num_morpho_steps,9);
for i=1:9
    qtf(1:peak_time_index,i)=linspace(qsf_low,qsf_peak,peak_time_index)*qtf_unit(i);
    qtf(peak_time_index:sediment_feed_end_index,i)=linspace(qsf_peak,qsf_low,(sediment_feed_end_index-peak_time_index+1))*qtf_unit(i);
end

% r=0.5
% qsf_low=45;
% qsf_peak=195;
% peak_time=50; % min
% qsf_low=qsf_low/(1000*2650*(1-p));
% qsf_peak=qsf_peak/(1000*2650*(1-p));
% peak_time_index=round(num_morpho_steps*peak_time/300);
% sediment_feed_end_index=round(num_morpho_steps*100/300);
% qtfT=zeros(num_morpho_steps,1);
% qtfT(1:peak_time_index)=linspace(qsf_low,qsf_peak,peak_time_index);
% qtfT(peak_time_index:sediment_feed_end_index)=linspace(qsf_peak,qsf_low,(sediment_feed_end_index-peak_time_index+1));
% qtf=zeros(num_morpho_steps,9);
% for i=1:9
%     qtf(1:peak_time_index,i)=linspace(qsf_low,qsf_peak,peak_time_index)*qtf_unit(i);
%     qtf(peak_time_index:sediment_feed_end_index,i)=linspace(qsf_peak,qsf_low,(sediment_feed_end_index-peak_time_index+1))*qtf_unit(i);
% end

% r=0.75
% qsf_low=45;
% qsf_peak=195;
% peak_time=75; % min
% qsf_low=qsf_low/(1000*2650*(1-p));
% qsf_peak=qsf_peak/(1000*2650*(1-p));
% peak_time_index=round(num_morpho_steps*peak_time/300);
% sediment_feed_end_index=round(num_morpho_steps*100/300);
% qtfT=zeros(num_morpho_steps,1);
% qtfT(1:peak_time_index)=linspace(qsf_low,qsf_peak,peak_time_index);
% qtfT(peak_time_index:sediment_feed_end_index)=linspace(qsf_peak,qsf_low,(sediment_feed_end_index-peak_time_index+1));
% qtf=zeros(num_morpho_steps,9);
% for i=1:9
%     qtf(1:peak_time_index,i)=linspace(qsf_low,qsf_peak,peak_time_index)*qtf_unit(i);
%     qtf(peak_time_index:sediment_feed_end_index,i)=linspace(qsf_peak,qsf_low,(sediment_feed_end_index-peak_time_index+1))*qtf_unit(i);
% end

%%%%%%Initial condition%%%%%%
%Fi: proportion of sediment on bed surface;
%FS: proportion of sand on bed surface;
%Dsg: geometrical mean grain size of bed surface;
%D90:grain size that 90 percent is finer;
%La: thickness of active layer;
Fi=ones(n,1)*fbase;
FS=sum(Fi(:,1:3),2);
Dsg=2.^(Fi*Phi')/1000;
D90=ones(n,1);
Cdg=[zeros(n,1),cumsum(Fi,2)];
for iD=1:n
    ind=find(Cdg(iD,:)>=0.9,1);
    D1=Phin(ind-1);
    D2=Phin(ind);
    P1=Cdg(iD,ind-1);
    P2=Cdg(iD,ind);
    D90(iD,1)=2^((0.9-P1)/(P2-P1)*(D2-D1)+D1)/1000;
end
La=na*D90;
%zb: bed elevation;
%zbini: initial bed elevation;
%s: bed slope;
x=linspace(0,L,n)';
zb=SI*(L-x);
zbini=zb;
s=[(zb(1:n-1)-zb(2:n))/dx;(zb(n-1)-zb(n))/dx];
t=0;
i=0;
%qw: flow discharge per unit width;
%h: water depth;
%Cf: resistance coefficient;
%u: flow velocity;
%FR: Froude number;

qw=qwf*ones(n+2,1);
h=qw.^(3/5).*(nk.*D90(1)).^(1/10)./(alpr)^(3/5)./(g*[s(1);s;s(n)]).^(3/10);


%qw=qwf*zeros(n+2,1);
%qw(1)=qwf
%h=zeros(n+2,1);
%h=qw.^(3/5).*(nk.*D90(1)).^(1/10)./(alpr)^(3/5)./(g*[s(1);s;s(n)]).^(3/10);



Cf=(nk.*D90./h(2:n+1)).^(1/3)./alpr.^2;
u=qw(2:n+1)./h(2:n+1);
FR=u./(g*h(2:n+1)).^0.5;

%Store: information of substrate stratigraphy;
%indup: number of sublayer at every node;
%Pup:proportion of sediment on uppermost sublayer;
%Lup: thickness of uppermost sublayer;
Store=cell(n,1);
indup=ceil((zb-La+4*Ls)./Ls);
for is=1:n
    Store{is,1}=zeros(300,ng);
    Store{is,1}(1:indup(is),:)=ones(indup(is),1)*fbase;
end
fup=ones(n,1)*fbase;
Lup=zb+4*Ls-La-Ls*(indup-1);

%%%%%%Variables%%%%%%
%left&right conditions for flux calculation
hl=zeros(n+1,1);
hr=zeros(n+1,1);
qwl=zeros(n+1,1);
qwr=zeros(n+1,1);
ul=zeros(n+1,1);
ur=zeros(n+1,1);
%numerical flux
Fl=zeros(n+1,2);
Fr=zeros(n+1,2);
Fhll=zeros(n+1,2);
F=zeros(n+1,2);
%source terms
ss=zeros(n,1);
sf=zeros(n,1);
%taub: bed shear stress;
%taur:reference shear stress;
%taurg:reference shear stress for surface geometric mean grain size;
%shirg: reference shields number for surface geometric mean grain size;
%ustar: shear velocity;
%b:component to compute taur;
%phi:taub on taur;
%Wstar: dimensionless bedload transport rate;
%qb: sediment transport rate;
%qbT: total volume sediment transport rate;
%qsT: total transport rate of sand;
%qgT: total transport rate of gravel;
%pbi: grain size distribution of bedload;
taub=zeros(n,1);
taur=zeros(n,ng);
taurg=zeros(n,1);
shirg=zeros(n,1);
ustar=zeros(n,1);
b=zeros(n,ng);
phi=zeros(n,ng);
Wstar=zeros(n,ng);
qb=zeros(n,ng);
qbT=zeros(n,1);
qsT=zeros(n,1);
qgT=zeros(n,1);
pbi=zeros(n,ng);
%fI: interfacial exchange fractions;
%dzb: change of bed evolution;
%dLa: change of activer layer thickness;
%delta: change of substrate elevation;
fI=zeros(n,ng);
dzb=zeros(n,1);
dLa=zeros(n,1);
dFi=zeros(n,ng);
delta=zeros(n,1);
%qwt: flow discharge per unit width at different time;
%zbt: bed elevation at different time;
%st: bed slope at different time;
%sct: (central scheme) bed slope at different time;
%Dsgt: Dsg at different time;
%FSt: FS at different time;
%Fit: Fi at different time;
%qbTt: total sediment transport rate per unit width at different time;
%Dg_loadt: Dg_load at different time;
qwt=zeros(n,Nt);
zbt=zeros(n,Nt);
st=zeros(n,Nt);
sct=zeros(n,Nt);
Dsgt=zeros(n,Nt);
FSt=zeros(n,Nt);
Fit=cell(Nt,1);
qbTt=zeros(n,Nt);
Dg_loadt=zeros(n,Nt);
ht=zeros(n,Nt);
zwt=zeros(n,Nt);
dzbt=zeros(n,Nt);
it=1;


%%%%%%Simulation%%%%%%
while t<Duration
    %%%%%%Hydraulics£ºShallow Water Equation%%%%%%
    for ih=1:nh
        %%%BC: upstream given discharge and mass balance, downstream bed resistance%%%
        qw(1)=qwf;
        hsub=h(1)+(qw(1)-F(1,1))*dth*timeyear/dx;
        hsuper=qw(1).^(3/5).*(nk.*D90(1)).^(1/10)./(alpr)^(3/5)./(g*s(1)).^(3/10);
        h(1)=(FR(1)>=1)*hsuper+(FR(1)<1)*hsub;
        qw(n+2)=qw(n+1);
        %h(n+2)=h(n+1);
        h(n+2)=qw(n+2).^(3/5).*(nk.*D90(n)).^(1/10)./(alpr)^(3/5)./(g*s(n)).^(3/10);
        %%%HLL scheme%%%
        %left&right conditions
        hl=h(1:n+1);
        hr=h(2:n+2);
        qwl=qw(1:n+1);
        qwr=qw(2:n+2);
        ul=qwl./hl;
        ur=qwr./hr;
        %wave speed
        cal=sqrt(g*hl);
        car=sqrt(g*hr);
        us=0.5*(ul+ur)+cal-car;
        cas=0.5*(cal+car)+0.25*(ul-ur);
        sl=min(ul-cal,us-cas)*ones(1,2);
        sr=max(ur+car,us+cas)*ones(1,2);
        %numerical flux
        Fl=[qwl,qwl.^2./hl+0.5*g*hl.^2];
        Fr=[qwr,qwr.^2./hr+0.5*g*hr.^2];
        Fhll=(sr.*Fl-sl.*Fr+sl.*sr.*[hr-hl,qwr-qwl])./(sr-sl);
        F=(sl<0).*(sr>0).*Fhll+(sl>=0).*Fl+(sr<=0).*Fr;
        %%%source term%%%
        ss=g.*h(2:n+1).*s;
        Cf=(nk.*D90./h(2:n+1)).^(1/3)./alpr.^2;
        sf=Cf.*qw(2:n+1)./h(2:n+1).^2;
        %%%time advance%%%
        h(2:n+1)=h(2:n+1)+(F(1:n,1)-F(2:n+1,1))/dx*dth*timeyear;
        qw(2:n+1)=(qw(2:n+1)+(F(1:n,2)-F(2:n+1,2))/dx*dth*timeyear+ss*dth*timeyear)./(1+sf*dth*timeyear);
        Courant=max(Courant,(qw(2:n+1)./h(2:n+1)+(g.*h(2:n+1)).^0.5)*dth*timeyear/dx);
    end
    u=qw(2:n+1)./h(2:n+1);
    FR=u./(g*h(2:n+1)).^0.5;
    
    %%%%%%Sediment transport: Wilcock and Crowe (2003)%%%%%%
    taub=rho.*Cf.*u.^2;
    ustar=(taub./rho).^0.5;
    shirg=nt*(0.021+0.015*exp(-20*FS));
    taurg=shirg.*R.*rho.*g.*Dsg;
    b=0.67./(1+exp(1.5-(ones(n,1)*D)./(Dsg*ones(1,ng))));
    taur=(taurg*ones(1,ng)).*((ones(n,1)*D)./(Dsg*ones(1,ng))).^b;
    phi=(taub*ones(1,ng))./taur;
    Wstar=(phi<1.35).*0.002.*(phi).^7.5+(phi>=1.35).*14.*(1-0.894./(phi).^0.5).^4.5;
    qb=Wstar.*Fi.*((ustar.^3)*ones(1,ng))./R./g;
    qbT=sum(qb,2);
    qsT=sum(qb(:,1:3),2);
    qgT=sum(qb(:,4:9),2);
    pbi=qb./(qbT*ones(1,ng));
    Dg_load=2.^(pbi*Phi')/1000;
    
    %%%%%%Exner equation%%%%%%
    %%%evolution of bed elevation
    % qbTback=[qtfT;qbT(1:n-1)];
    qbTback=[qtfT(i+1);qbT(1:n-1)];
    qbTit=qbT;
    qbTfrnt=[qbT(2:n);2*qbT(n)-qbT(n-1)];
    qbTdif=au*(qbTit-qbTback)+(1-au)*(qbTfrnt-qbTit);
    dzb=-qbTdif/dx/(1-p)*dt*timeyear;
    dzb(n)=0;
    zb=zb+dzb;
    %%%evolution of surface fraction
    delta=dzb-dLa;
    fI=((delta<=0)*ones(1,ng)).*fup+((delta>0)*ones(1,ng)).*(Fi+pbi)/2;
    qbback=[qtf(i+1,:);qb(1:n-1,:)];
    qbit=qb;
    qbfrnt=[qb(2:n,:);2*qb(n,:)-qb(n-1,:)];
    qbdif=au*(qbit-qbback)+(1-au)*(qbfrnt-qbit);
    dFi=((-qbdif+qbTdif*ones(1,ng).*fI)/dx*dt*timeyear/(1-p)-(Fi-fI).*(dLa*ones(1,ng)))./(La*ones(1,ng));
    Fi=Fi+dFi;
    Fi=(Fi>0).*Fi;
    Fi=Fi./(sum(Fi,2)*ones(1,ng));
    
    %%%%%%Stratigraphy storage%%%%%%
    indn=find((delta<=Ls-Lup)&(delta>=-Lup));
    indinc=find(delta>Ls-Lup);
    inddec=find(delta<-Lup);
    fup(indn,:)=(fup(indn,:).*(Lup(indn)*ones(1,ng))+fI(indn,:).*(delta(indn)*ones(1,ng)))./((Lup(indn)+delta(indn))*ones(1,ng));
    Lup(indn)=Lup(indn)+delta(indn);
    if size(indinc,1)>0
        for is=1:size(indinc,1)
            ii=indinc(is);
            fup(ii,:)=(fup(ii,:).*Lup(ii)+fI(ii,:).*(Ls-Lup(ii)))./Ls;
            Store{ii,1}(indup(ii),:)=fup(ii,:);
            inc=ceil((delta(ii)-(Ls-Lup(ii)))/Ls);
            Store{ii,1}(indup(ii)+1:indup(ii)+inc,:)=ones(inc,1)*fI(ii,:);
            indup(ii)=indup(ii)+inc;
            fup(ii,:)=fI(ii,:);
            Lup(ii)=delta(ii)-(Ls-Lup(ii))-(inc-1)*Ls;
        end
    end
    if size(inddec,1)>0
        for is=1:size(inddec,1)
            id=inddec(is);
            dec=ceil((-Lup(id)-delta(id))/Ls);
            Store{id,1}(indup(id)-dec+1:indup(id),:)=zeros(dec,ng);
            indup(id)=indup(id)-dec;
            fup(id,:)=Store{id,1}(indup(id),:);
            Lup(id)=dec*Ls+Lup(id)+delta(id);
        end
    end
    
    %%%%%%Parameter update%%%%%%
    t=t+dt;
    i=i+1;
    s=[(zb(1:n-1)-zb(2:n))/dx;(zb(n-1)-zb(n))/dx];
    sc=[(zb(1)-zb(2))/dx;(zb(1:n-2)-zb(3:n))/2/dx;(zb(n-1)-zb(n))/dx];
    FS=sum(Fi(:,1:3),2);
    Dsg=2.^(Fi*Phi')/1000;
    Cdg=[zeros(n,1),cumsum(Fi,2)];
    for iD=1:n
        ind=find(Cdg(iD,:)>=0.9,1);
        D1=Phin(ind-1);
        D2=Phin(ind);
        P1=Cdg(iD,ind-1);
        P2=Cdg(iD,ind);
        D90(iD,1)=2^((0.9-P1)/(P2-P1)*(D2-D1)+D1)/1000;
    end
    dLa=na*D90-La;
    La=na*D90;
    
    %%%record information every several years
    if (mod(i,interval)==0)
        qwt(:,it)=qw(2:n+1);
        zbt(:,it)=zb;
        st(:,it)=s;
        sct(:,it)=sc;
        Dsgt(:,it)=Dsg;
        FSt(:,it)=FS;
        Fit{it,1}=Fi;
        qbTt(:,it)=qbT;
        Dg_loadt(:,it)=Dg_load;
        ht(:,it)=h(2:n+1);
        zwt(:,it)=h(2:n+1)+zb;
        dzbt(:,it)=zb-zbini;
        disp(it);
        it=it+1;
    end
    
    %figure(1);plot(x,qw,'*r')
    %figure(1);plot(x,h,'*r')
    % figure(1);plot(x,zb,'*r')
    %figure(1);plot(x,sc,'*r')
    %figure(1);plot(x,Dsg,'*r')
    %figure(1);plot(x,Dg_load,'*r')
    % drawnow
end

for is=1:n
    Store{is,1}(indup(is),:)=fup(is,:);
end
% save results
save([mfilename,'.mat']);
tc=(cputime-tc)/60