function [sys,x0,str,ts] = trace_planning(t,x,u,flag)
switch flag
 case 0
  [sys,x0,str,ts] = mdlInitializeSizes; % Initialization
  
 case 2
  sys = mdlUpdates(t,x,u); % Update discrete states
  
 case 3
  sys = mdlOutputs(t,x,u); % Calculate outputs
 
%  case 4
%   sys = mdlGetTimeOfNextVarHit(t,x,u); % Get next sample time 

 case {1,4,9} % Unused flags
  sys = [];
  
 otherwise
  error(['unhandled flag = ',num2str(flag)]); % Error handling
end
% End of dsfunc.

%==============================================================
% Initialization
%==============================================================

function [sys,x0,str,ts] = mdlInitializeSizes

sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 6;
sizes.NumOutputs     = 5;
sizes.NumInputs      = 15;
sizes.DirFeedthrough = 1; % Matrix D is non-empty.
sizes.NumSampleTimes = 1;
sys = simsizes(sizes); 
x0 =[0.001;0.0001;0.0001;0.00001;0.0001;0.0001];    
global U;
U=0;%��������ʼ��,���������һ�������켣����������ȥ����UΪһά��
% global x;
% x = zeros(md.ne + md.pye + md.me + md.Hu*md.me,1);   
% Initialize the discrete states.
str = [];             % Set str to an empty matrix.
ts  = [0.02 0];       % sample time: [period, offset]
%End of mdlInitializeSizes
		      
%==============================================================
% Update the discrete states
%==============================================================
function sys = mdlUpdates(t,x,u)
  
sys = x;
%End of mdlUpdate.

%==============================================================
% Calculate outputs
%==============================================================
function sys = mdlOutputs(t,~,u)
    global a b; 
    %global u_piao;
    global U;
    %global kesi;
    tic
    Nx=6;%״̬���ĸ���
    Nu=1;%�������ĸ���
    Ny=2;%������ĸ���
    Np =15;%Ԥ�ⲽ��
    Nc=1;%���Ʋ���
    Ny_cepian = 4;
    Row=1000;%�ɳ�����Ȩ��
    c = 1;
    fprintf('Update start, t=%6.3f\n',t)
   
    %����ӿ�ת��,x_dot�����һ���ǳ�С�������Ƿ�ֹ���ַ�ĸΪ������
   % y_dot=u(1)/3.6-0.000000001*0.4786;%CarSim�������km/h��ת��Ϊm/s
    y_dot=u(1)/3.6;
    x_dot=u(2)/3.6+0.0001;%CarSim�������km/h��ת��Ϊm/s
    phi=u(3)*3.141592654/180; %CarSim�����Ϊ�Ƕȣ��Ƕ�ת��Ϊ����
    phi_dot=u(4)*3.141592654/180;
    Y=u(6);%��λΪm
    X=u(5);%��λΪ��
    Y_dot=u(7);
    X_dot=u(8);
%% ������������
%syms sf sr;%�ֱ�Ϊǰ���ֵĻ�����,��Ҫ�ṩ
    Sf=0.2; Sr=0.2;
%syms lf lr;%ǰ���־��복�����ĵľ��룬�������в���
    lf=1.232;lr=1.468;
%syms C_cf C_cr C_lf C_lr;%�ֱ�Ϊǰ���ֵ��ݺ����ƫ�նȣ��������в���
    Ccf=66900;Ccr=62700;Clf=66900;Clr=62700;
%syms m g I;%mΪ����������gΪ�������ٶȣ�IΪ������Z���ת���������������в���
    m=1723;g=9.8;I=4175;
    w=3.75;
    %sf=34;%ssss�Ǳ���ֱ����ʻ�ľ���
    %sf=55.8;
    %%
    %��Ҫ����
    sf=33.333;
    %
    d=78.329;
%  load('dlc.mat')
%% �ο��켣����
    shape=2.4;%�������ƣ����ڲο��켣����
    dx1=25;dx2=21.95;%û���κ�ʵ�����壬ֻ�ǲ�������
    dy1=4.05;dy2=5.7;%û���κ�ʵ�����壬ֻ�ǲ�������
    Xs1=27.19;Xs2=56.46;%��������
    X_predict=zeros(Np,1);%���ڱ���Ԥ��ʱ���ڵ�����λ����Ϣ�����Ǽ��������켣�Ļ���
    phi_ref=zeros(Np,1);%���ڱ���Ԥ��ʱ���ڵ������켣
    
    %  ���¼���kesi,��״̬�������������һ��   
    kesi=zeros(Nx+Nu,1);
    kesi(1)=y_dot;%u(1)==X(1)
    kesi(2)=x_dot;%u(2)==X(2)
    kesi(3)=phi; %u(3)==X(3)
    kesi(4)=phi_dot;
    kesi(5)=Y;
    kesi(6)=X;
    kesi(7)=U(1);
    delta_f=U(1);
    fprintf('Update start, u(1)=%4.2f\n',U(1));
    T=0.02;%���沽��
%     T_all=20;%�ܵķ���ʱ�䣬��Ҫ�����Ƿ�ֹ���������켣Խ��
         %Ȩ�ؾ������� 
    Q_cell=cell(Np,Np);
    for i=1:1:Np
        for j=1:1:Np
            if i==j 
                %Q_cell{i,j}=[200 0;0 100];
                 Q_cell{i,j}=[200 0;0 100];

            else 
                Q_cell{i,j}=zeros(Ny,Ny);      
  
            end
        end 
    end 

 M_cell = cell(Nc,Nc);
   for i=1:1:Nc
        for j=1:1:Nc
            if j<=i
               M_cell{i,j}=1;  
            else 
               M_cell{i,j}=zeros(Nu,Nu);              
            end
        end 
   end 
    MM = cell2mat(M_cell);
    
%     R=0.65*10^6*eye(Nu*Nc);
  R=10*eye(Nu*Nc);
    %�����Ҳ����Ҫ�ľ����ǿ������Ļ��������ö���ѧģ�ͣ��þ����복������������أ�ͨ���Զ���ѧ��������ſ˱Ⱦ���õ�
     a = [               1.0 - (150.44*T)/x_dot, -1.0*T*(phi_dot + (0.0011608*(92044.0*phi_dot - 62700.0*y_dot))/x_dot^2 - (77.655*(1.232*phi_dot + y_dot))/x_dot^2),                                        0,       -1.0*T*(x_dot - 11.17/x_dot),   0,   0
 T*(phi_dot + (77.655*delta_f)/x_dot),                                                            1.0 - (77.655*T*delta_f*(1.232*phi_dot + y_dot))/x_dot^2,                                        0, T*(y_dot + (95.671*delta_f)/x_dot),   0,   0
                                    0,                                                                                                                   0,                                      1.0,                                  T,   0,   0
                     (4.6097*T)/x_dot,             T*((39.483*(1.232*phi_dot + y_dot))/x_dot^2 + (0.00023952*(270244.0*phi_dot - 184099.0*y_dot))/x_dot^2),                                        0,             1.0 - (113.37*T)/x_dot,   0,   0
                           T*cos(phi),                                                                                                          T*sin(phi),  T*(x_dot*cos(phi) - 1.0*y_dot*sin(phi)),                                  0, 1.0,   0
                     -1.0*T*sin(phi),                                                                                                          T*cos(phi), -1.0*T*(y_dot*cos(phi) + x_dot*sin(phi)),                                  0,   0, 1.0 ];
 
     b=[                                              
                                                         77.655*T
 -1.0*T*(155.31*delta_f - (77.655*(1.232*phi_dot + y_dot))/x_dot)
                                                                0
                                                         39.483*T   
                                                                0 
                                                                0
 ];
    d_k=zeros(Nx,1);%����ƫ��
    state_k1=zeros(Nx,1);%Ԥ�����һʱ��״̬�������ڼ���ƫ��
    %���¼�Ϊ������ɢ������ģ��Ԥ����һʱ��״̬��
    %ע�⣬Ϊ����ǰ�����ı��ʽ��a,b�����������a,b�س�ͻ����ǰ�����ı��ʽ��Ϊlf��lr
    state_k1(1,1)=y_dot+T*(-x_dot*phi_dot+2*(Ccf*(delta_f-(y_dot+lf*phi_dot)/x_dot)+Ccr*(lr*phi_dot-y_dot)/x_dot)/m);
    state_k1(2,1)=x_dot+T*(y_dot*phi_dot+2*(Clf*Sf+Clr*Sr+Ccf*delta_f*(delta_f-(y_dot+phi_dot*lf)/x_dot))/m);
    state_k1(3,1)=phi+T*phi_dot;
    state_k1(4,1)=phi_dot+T*((2*lf*Ccf*(delta_f-(y_dot+lf*phi_dot)/x_dot)-2*lr*Ccr*(lr*phi_dot-y_dot)/x_dot)/I);
    state_k1(5,1)=Y+T*(x_dot*sin(phi)+y_dot*cos(phi));
    state_k1(6,1)=X+T*(x_dot*cos(phi)-y_dot*sin(phi));
    d_k=state_k1-a*kesi(1:6,1)-b*kesi(7,1);%����falcone��ʽ��2.11b�����d(k,t)
    d_piao_k=zeros(Nx+Nu,1);%d_k��������ʽ���ο�falcone(B,4c)
    d_piao_k(1:6,1)=d_k;
    d_piao_k(7,1)=0;
    %% ��ƫ��e(t)���� �����Ҫ��������Ĳ�ƫ����Լ�� ����ʵ�ֳ����ڼ��˹����µİ�ȫ��ʻ
    y_dot = state_k1(1,1);  %%  y_dot  x_dot  phi_dot��Ҫ�������� Alpha_fl r ..�Ȳ�ƫ��
    x_dot = state_k1(2,1);  %%  y_dot  x_dot  phi_dot��Ҫ�������� Alpha_fl r ..�Ȳ�ƫ��
    phi_dot = state_k1(4,1);%%  y_dot  x_dot  phi_dot��Ҫ�������� Alpha_fl r ..�Ȳ�ƫ��
    Alpha_fl = atan(( (y_dot+lf*phi_dot)*cos(U) - (x_dot-c*phi_dot)*sin(U) )/( (y_dot+lf*phi_dot)*sin(U) + (x_dot-c*phi_dot)*cos(U)) );
    Alpha_fr = atan(( (y_dot+lf*phi_dot)*cos(U) - (x_dot+c*phi_dot)*sin(U) )/( (y_dot+lf*phi_dot)*sin(U) + (x_dot+c*phi_dot)*cos(U)) );
    Alpha_rl = atan(( (y_dot-lr*phi_dot)*cos(U) - (x_dot-c*phi_dot)*sin(U) )/( (y_dot-lr*phi_dot)*sin(U) + (x_dot-c*phi_dot)*cos(U)) );
    Alpha_rr = atan(( (y_dot-lr*phi_dot)*cos(U) - (x_dot+c*phi_dot)*sin(U) )/( (y_dot-lr*phi_dot)*sin(U) + (x_dot+c*phi_dot)*cos(U)) );
    state2_k1(1,1) = Alpha_fl;
    state2_k1(2,1) = Alpha_fr; 
    state2_k1(3,1) = Alpha_rl;
    state2_k1(4,1) = Alpha_rr;
D_d = [  -1
 -1
 -1
 -1]; %% ��ƫ��������� D
C_d =  [ -(1.0*(cos(delta_f)/(cos(delta_f)*(phi_dot - 1.0*x_dot) - 1.0*sin(delta_f)*(1.232*phi_dot + y_dot)) + (sin(delta_f)*(cos(delta_f)*(1.232*phi_dot + y_dot) + sin(delta_f)*(phi_dot - 1.0*x_dot)))/(cos(delta_f)*(phi_dot - 1.0*x_dot) - 1.0*sin(delta_f)*(1.232*phi_dot + y_dot))^2))/((cos(delta_f)*(1.232*phi_dot + y_dot) + sin(delta_f)*(phi_dot - 1.0*x_dot))^2/(cos(delta_f)*(phi_dot - 1.0*x_dot) - 1.0*sin(delta_f)*(1.232*phi_dot + y_dot))^2 + 1.0),             (sin(delta_f)/(cos(delta_f)*(phi_dot - 1.0*x_dot) - 1.0*sin(delta_f)*(1.232*phi_dot + y_dot)) - (1.0*cos(delta_f)*(cos(delta_f)*(1.232*phi_dot + y_dot) + sin(delta_f)*(phi_dot - 1.0*x_dot)))/(cos(delta_f)*(phi_dot - 1.0*x_dot) - 1.0*sin(delta_f)*(1.232*phi_dot + y_dot))^2)/((cos(delta_f)*(1.232*phi_dot + y_dot) + sin(delta_f)*(phi_dot - 1.0*x_dot))^2/(cos(delta_f)*(phi_dot - 1.0*x_dot) - 1.0*sin(delta_f)*(1.232*phi_dot + y_dot))^2 + 1.0), 0,              -(1.0*((1.232*cos(delta_f) + sin(delta_f))/(cos(delta_f)*(phi_dot - 1.0*x_dot) - 1.0*sin(delta_f)*(1.232*phi_dot + y_dot)) - (1.0*(cos(delta_f)*(1.232*phi_dot + y_dot) + sin(delta_f)*(phi_dot - 1.0*x_dot))*(cos(delta_f) - 1.232*sin(delta_f)))/(cos(delta_f)*(phi_dot - 1.0*x_dot) - 1.0*sin(delta_f)*(1.232*phi_dot + y_dot))^2))/((cos(delta_f)*(1.232*phi_dot + y_dot) + sin(delta_f)*(phi_dot - 1.0*x_dot))^2/(cos(delta_f)*(phi_dot - 1.0*x_dot) - 1.0*sin(delta_f)*(1.232*phi_dot + y_dot))^2 + 1.0), 0, 0
        (cos(delta_f)/(sin(delta_f)*(1.232*phi_dot + y_dot) + cos(delta_f)*(phi_dot + x_dot)) - (sin(delta_f)*(1.0*cos(delta_f)*(1.232*phi_dot + y_dot) - sin(delta_f)*(phi_dot + x_dot)))/(sin(delta_f)*(1.232*phi_dot + y_dot) + cos(delta_f)*(phi_dot + x_dot))^2)/((1.0*cos(delta_f)*(1.232*phi_dot + y_dot) - sin(delta_f)*(phi_dot + x_dot))^2/(sin(delta_f)*(1.232*phi_dot + y_dot) + cos(delta_f)*(phi_dot + x_dot))^2 + 1.0),                              -(1.0*(sin(delta_f)/(sin(delta_f)*(1.232*phi_dot + y_dot) + cos(delta_f)*(phi_dot + x_dot)) + (1.0*cos(delta_f)*(1.0*cos(delta_f)*(1.232*phi_dot + y_dot) - sin(delta_f)*(phi_dot + x_dot)))/(sin(delta_f)*(1.232*phi_dot + y_dot) + cos(delta_f)*(phi_dot + x_dot))^2))/((1.0*cos(delta_f)*(1.232*phi_dot + y_dot) - sin(delta_f)*(phi_dot + x_dot))^2/(sin(delta_f)*(1.232*phi_dot + y_dot) + cos(delta_f)*(phi_dot + x_dot))^2 + 1.0), 0,                                             ((1.232*cos(delta_f) - 1.0*sin(delta_f))/(sin(delta_f)*(1.232*phi_dot + y_dot) + cos(delta_f)*(phi_dot + x_dot)) - ((1.0*cos(delta_f)*(1.232*phi_dot + y_dot) - sin(delta_f)*(phi_dot + x_dot))*(cos(delta_f) + 1.232*sin(delta_f)))/(sin(delta_f)*(1.232*phi_dot + y_dot) + cos(delta_f)*(phi_dot + x_dot))^2)/((1.0*cos(delta_f)*(1.232*phi_dot + y_dot) - sin(delta_f)*(phi_dot + x_dot))^2/(sin(delta_f)*(1.232*phi_dot + y_dot) + cos(delta_f)*(phi_dot + x_dot))^2 + 1.0), 0, 0
        -(1.0*(cos(delta_f)/(cos(delta_f)*(phi_dot - 1.0*x_dot) + sin(delta_f)*(1.468*phi_dot - 1.0*y_dot)) + (1.0*sin(delta_f)*(1.0*sin(delta_f)*(phi_dot - 1.0*x_dot) - cos(delta_f)*(1.468*phi_dot - 1.0*y_dot)))/(cos(delta_f)*(phi_dot - 1.0*x_dot) + sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2))/((1.0*sin(delta_f)*(phi_dot - 1.0*x_dot) - cos(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2/(cos(delta_f)*(phi_dot - 1.0*x_dot) + sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2 + 1.0), (sin(delta_f)/(cos(delta_f)*(phi_dot - 1.0*x_dot) + sin(delta_f)*(1.468*phi_dot - 1.0*y_dot)) - (cos(delta_f)*(1.0*sin(delta_f)*(phi_dot - 1.0*x_dot) - cos(delta_f)*(1.468*phi_dot - 1.0*y_dot)))/(cos(delta_f)*(phi_dot - 1.0*x_dot) + sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2)/((1.0*sin(delta_f)*(phi_dot - 1.0*x_dot) - cos(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2/(cos(delta_f)*(phi_dot - 1.0*x_dot) + sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2 + 1.0), 0, ((1.468*cos(delta_f) - 1.0*sin(delta_f))/(cos(delta_f)*(phi_dot - 1.0*x_dot) + sin(delta_f)*(1.468*phi_dot - 1.0*y_dot)) + (1.0*(1.0*sin(delta_f)*(phi_dot - 1.0*x_dot) - cos(delta_f)*(1.468*phi_dot - 1.0*y_dot))*(cos(delta_f) + 1.468*sin(delta_f)))/(cos(delta_f)*(phi_dot - 1.0*x_dot) + sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2)/((1.0*sin(delta_f)*(phi_dot - 1.0*x_dot) - cos(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2/(cos(delta_f)*(phi_dot - 1.0*x_dot) + sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2 + 1.0), 0, 0
        (1.0*(cos(delta_f)/(1.0*cos(delta_f)*(phi_dot + x_dot) - sin(delta_f)*(1.468*phi_dot - 1.0*y_dot)) + (1.0*sin(delta_f)*(sin(delta_f)*(phi_dot + x_dot) + cos(delta_f)*(1.468*phi_dot - 1.0*y_dot)))/(1.0*cos(delta_f)*(phi_dot + x_dot) - sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2))/((sin(delta_f)*(phi_dot + x_dot) + cos(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2/(1.0*cos(delta_f)*(phi_dot + x_dot) - sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2 + 1.0),                -(sin(delta_f)/(1.0*cos(delta_f)*(phi_dot + x_dot) - sin(delta_f)*(1.468*phi_dot - 1.0*y_dot)) - (cos(delta_f)*(sin(delta_f)*(phi_dot + x_dot) + cos(delta_f)*(1.468*phi_dot - 1.0*y_dot)))/(1.0*cos(delta_f)*(phi_dot + x_dot) - sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2)/((sin(delta_f)*(phi_dot + x_dot) + cos(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2/(1.0*cos(delta_f)*(phi_dot + x_dot) - sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2 + 1.0), 0,                        -((1.468*cos(delta_f) + sin(delta_f))/(1.0*cos(delta_f)*(phi_dot + x_dot) - sin(delta_f)*(1.468*phi_dot - 1.0*y_dot)) - ((sin(delta_f)*(phi_dot + x_dot) + cos(delta_f)*(1.468*phi_dot - 1.0*y_dot))*(cos(delta_f) - 1.468*sin(delta_f)))/(1.0*cos(delta_f)*(phi_dot + x_dot) - sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2)/((sin(delta_f)*(phi_dot + x_dot) + cos(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2/(1.0*cos(delta_f)*(phi_dot + x_dot) - sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2 + 1.0), 0, 0
     ]; %% ��ƫ����� = C_d ����״̬��kesi���� + D_d *delat_u + e_k
 
    CC_d_cell = cell(1,2);
    CC_d_cell{1, 1} = C_d;
    CC_d_cell{1, 2} = D_d;
    Cc_d = cell2mat(CC_d_cell);
    side_slip_error = state2_k1 - C_d*kesi(1:6,1) - D_d*kesi(7,1);%����������¹��������µ�ֵ ������Ҫ����δ��ÿ��Ԥ��ʱ���ڵĲ�ֵ����ò���ǿ�����Ϊ���ȣ��˴�ֻ��Ϊ�˱任��������ʾ���ɿ����Ļ���ģ��Ԥ������㷨������Ƶ���ʽ��
    sideslip_piao_k=zeros(Ny_cepian, 1);%d_k��������ʽ���ο�falcone(B,4c)
    sideslip_piao_k(1:4,1)=side_slip_error;
  
    
    
    %% 
    A_cell=cell(2,2);
    B_cell=cell(2,1);
    A_cell{1,1}=a;
    A_cell{1,2}=b;
    A_cell{2,1}=zeros(Nu,Nx);
    A_cell{2,2}=eye(Nu);
    B_cell{1,1}=b;
    B_cell{2,1}=eye(Nu);
    %A=zeros(Nu+Nx,Nu+Nx);
    A=cell2mat(A_cell);
    B=cell2mat(B_cell);
    C=[0 0 1 0 0 0 0;0 0 0 0 1 0 0;];
    PSI_cell=cell(Np,1);
    PSI_sideslip_cell = cell(Np,1);
    THETA_cell=cell(Np,Nc);
    THETA_sideslip_cell=cell(Np,Nc);   
    THETA_sideslip_cell2=cell(Np,Nc);  
    GAMMA_cell=cell(Np,Np);
    GAMMA_sideslip_cell=cell(Np,Np);
    PHI_cell=cell(Np,1);
    PHI_sideslip_cell=cell(Np,1);
    for p=1:1:Np
        PHI_cell{p,1} = d_piao_k;%��������˵�������Ҫʵʱ���µģ�����Ϊ�˼�㣬������һ�ν���
        PHI_sideslip_cell{p,1} = sideslip_piao_k;
        for q=1:1:Np
            if q<=p
                GAMMA_cell{p,q}=C*A^(p-q);
                GAMMA_sideslip_cell{p,q}=Cc_d*A^(p-q);
            else 
                GAMMA_cell{p,q} = zeros(Ny,Nx+Nu);
                GAMMA_sideslip_cell{p,q} = zeros(Ny_cepian,Nx+Nu);
            end 
        end
    end
    for j=1:1:Np
        PSI_cell{j,1}=C*A^j;
        PSI_sideslip_cell{j,1} = Cc_d*A^j;
        for k=1:1:Nc
            if k<=j
                THETA_cell{j,k}=C*A^(j-k)*B;
                THETA_sideslip_cell{j,k} = Cc_d*A^(j-k)*B;
            else 
                THETA_cell{j,k}=zeros(Ny,Nu);
                THETA_sideslip_cell{j,k} = zeros(Ny_cepian,Nu);
            end
        end
    end
    
     for j=1:1:Np
        for k=1:1:Nc
            if k == j+1
                THETA_sideslip_cell2{j,k} = D_d;
            else 
                THETA_sideslip_cell2{j,k} = zeros(Ny_cepian,Nu);
            end
        end
     end
    
    PSI=cell2mat(PSI_cell);%size(PSI)=[Ny*Np Nx+Nu] 
    THETA=cell2mat(THETA_cell);%size(THETA)=[Ny*Np Nu*Nc]
    GAMMA=cell2mat(GAMMA_cell);%��д��GAMMA
    PHI=cell2mat(PHI_cell);
    
    PHI_sideslip = cell2mat(PHI_sideslip_cell);
    GAMMA_sideslip = cell2mat(GAMMA_sideslip_cell);
    PSI_sideslip = cell2mat(PSI_sideslip_cell);
    THETA_sideslip  = cell2mat(THETA_sideslip_cell) + cell2mat(THETA_sideslip_cell2);
    
    Q = cell2mat(Q_cell);
    S = 100;
    H_cell=cell(2,2);
    H_cell{1,1}=2*(THETA'*Q*THETA+R+MM'*S*MM);
    H_cell{1,2}=zeros(Nu*Nc,1);
    H_cell{2,1}=zeros(1,Nu*Nc);
    H_cell{2,2}=Row;
    H=cell2mat(H_cell);
    Yita_ref_cell = cell(Np,1);

     for p=1:1:Np

               X_DOT=x_dot*cos(phi)-y_dot*sin(phi);%��������ϵ�������ٶ�,���Ǵ������ϵ�µ�
               X_predict(Np,1)=X+X_DOT*p*T;%���ȼ����δ��X��λ�ã�          
               x1=X+X_DOT*p*T-sf;
               x2=x1-d-10;
               Y_ref(p,1)=0.*(x1<0)+(10*w/(d^3)* x1^3-15*w/(d^4)* x1^4+6*w/(d^5)*x1^5).*((0<x1)&&(x1<=d))+3.75.*(x1>d);%��������·��
               phi_ref(p,1)=0.*(x1<0)+(atan(30*w/d^3* x1^2-60*w/d^4* x1^3+30*w/d^5*x1^4)).*((0<x1)&&(x1<=d))+0.*(x1>d);%������ڽ�  
               %Y_ref(p,1)=0.*(x1<0)+(10*w/(d^3)* x1^3-15*w/(d^4)* x1^4+6*w/(d^5)*x1^5).*((0<x1)&&(x1<=d))+3.75.*((d<x1)&&(x1<=d+10))+(3.75-((10*w/(d^3)* x2^3-15*w/(d^4)* x2^4+6*w/(d^5)*x2^5))).*((d+10<=x1)&&(x1<=2*d+10))+0.*(x1>=2*d+10);%��������·��
%                phi_ref(p,1)=0.*(x1<0)+(atan(30*w/d^3* x1^2-60*w/d^4* x1^3+30*w/d^5*x1^4)).*((0<x1)&&(x1<=d))+0.*(x1>d);%������ڽ�  
               Yita_ref_cell{p,1}=[phi_ref(p,1);Y_ref(p,1)];
    end     
    Yita_ref=cell2mat(Yita_ref_cell);
    error_1=-Yita_ref+PSI*kesi+GAMMA*PHI; %��ƫ��
    UU = ones(Nc,1)*U(1);
    f_cell=cell(1,2);
    f_cell{1,1}=2*error_1'*Q*THETA + 2*UU'*S*MM;
    f_cell{1,2}=0;
    f=cell2mat(f_cell);
    
 %% ����ΪԼ����������
 %������Լ��
    A_t=zeros(Nc,Nc);%��falcone���� P181
    for p=1:1:Nc
        for q=1:1:Nc
            if q<=p 
                A_t(p,q)=1;
            else 
                A_t(p,q)=0;
            end
        end 
    end 
    A_I=kron(A_t,eye(Nu));%������ڿ˻�
    Ut=kron(ones(Nc,1),U(1));
    umin = -0.1745;%ά������Ʊ����ĸ�����ͬ
    umax = 0.1745;
    delta_umin=-0.014*1;
    delta_umax=0.014*1;
    Umin=kron(ones(Nc,1),umin);
    Umax=kron(ones(Nc,1),umax);
    
    
     %��ƫ���������Լ��
    y_cepian_max = [ 0.0489; 0.0489; 0.0489; 0.0489];%%%%%%%%%%%+-2.8��
    y_cepian_min = [ -0.0489; -0.0489; -0.0489; -0.0489];
    Y_cepian_max = kron(ones(Np,1),y_cepian_max);
    Y_cepian_min = kron(ones(Np,1),y_cepian_min);
    
    %����y�ͺ��������ӲԼ��
    ycmax=[0.21;5];
    ycmin=[-0.3;-3];
    Ycmax=kron(ones(Np,1),ycmax);
    Ycmin=kron(ones(Np,1),ycmin);
    
    %��Ͽ�����Լ���������Լ������Լ��
    A_cons_cell={A_I zeros(Nu*Nc,1);-A_I zeros(Nu*Nc,1); THETA_sideslip -1*ones(Ny_cepian*Np,1);-THETA_sideslip ones(Ny_cepian*Np,1);THETA zeros(Ny*Np,1);-THETA zeros(Ny*Np,1)};
    b_cons_cell={Umax-Ut;-Umin+Ut; Y_cepian_max-PSI_sideslip*kesi-GAMMA_sideslip*PHI - PHI_sideslip; -Y_cepian_min+PSI_sideslip*kesi+GAMMA_sideslip*PHI + PHI_sideslip;Ycmax-PSI*kesi-GAMMA*PHI;-Ycmin+PSI*kesi+GAMMA*PHI};
    A_cons=cell2mat(A_cons_cell);%����ⷽ�̣�״̬������ʽԼ���������ת��Ϊ����ֵ��ȡֵ��Χ
    b_cons=cell2mat(b_cons_cell);%����ⷽ�̣�״̬������ʽԼ����ȡֵ
       %״̬��Լ��
    M=0.2*pi/180;   %��Լ�����ֵ �˴�M��������Լ��Լ��ֵ������ֵ̫�����׳���ʧ�ȣ�ʧȥ������Լ�������壬̫С��������޽⣬���ƶ�ʧ��
%       M=0;  %��Լ�����ֵ �˴�M��������Լ��Լ��ֵ������ֵ̫�����׳���ʧ�ȣ�ʧȥ������Լ�������壬̫С��������޽⣬���ƶ�ʧ��
    delta_Umin=kron(ones(Nc,1),delta_umin);
    delta_Umax=kron(ones(Nc,1),delta_umax);
    lb=[delta_Umin;0];%����ⷽ�̣�״̬���½磬��������ʱ���ڿ����������ɳ�����
    ub=[delta_Umax;M];%����ⷽ�̣�״̬���Ͻ磬��������ʱ���ڿ����������ɳ�����
    
    %% ��ʼ������
%        options = optimset('Algorithm','interior-point-convex');
options = optimset('Algorithm','active-set');
       x_start=zeros(Nc+1,1);%����һ����ʼ��
      [X,fval,exitflag]=quadprog(H,f,A_cons,b_cons,[],[],lb,ub,x_start,options);
      fprintf('exitflag=%d\n',exitflag); 
      fprintf('H=%4.2f\n',H(1,1));
      fprintf('f=%4.2f\n',f(1,1));
    %% �������
    u_piao=X(1);%�õ���������
    U(1)=kesi(7,1)+u_piao;%��ǰʱ�̵Ŀ�����Ϊ��һ��ʱ�̿���+��������
   
   %U(2)=Yita_ref(2);%���dphi_ref
    sys(1)= U(1);
    sys(2)=  phi_ref(p,1);
    sys(3)= Y_ref(p,1);
    sys(4)= X_predict(Np,1);
    sys(5)=u_piao;
    toc
% End of mdlOutputs.
  

   
