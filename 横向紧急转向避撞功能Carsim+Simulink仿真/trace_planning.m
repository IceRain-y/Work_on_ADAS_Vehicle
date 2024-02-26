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
U=0;%控制量初始化,这里面加了一个期望轨迹的输出，如果去掉，U为一维的
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
    Nx=6;%状态量的个数
    Nu=1;%控制量的个数
    Ny=2;%输出量的个数
    Np =15;%预测步长
    Nc=1;%控制步长
    Ny_cepian = 4;
    Row=1000;%松弛因子权重
    c = 1;
    fprintf('Update start, t=%6.3f\n',t)
   
    %输入接口转换,x_dot后面加一个非常小的数，是防止出现分母为零的情况
   % y_dot=u(1)/3.6-0.000000001*0.4786;%CarSim输出的是km/h，转换为m/s
    y_dot=u(1)/3.6;
    x_dot=u(2)/3.6+0.0001;%CarSim输出的是km/h，转换为m/s
    phi=u(3)*3.141592654/180; %CarSim输出的为角度，角度转换为弧度
    phi_dot=u(4)*3.141592654/180;
    Y=u(6);%单位为m
    X=u(5);%单位为米
    Y_dot=u(7);
    X_dot=u(8);
%% 车辆参数输入
%syms sf sr;%分别为前后车轮的滑移率,需要提供
    Sf=0.2; Sr=0.2;
%syms lf lr;%前后车轮距离车辆质心的距离，车辆固有参数
    lf=1.232;lr=1.468;
%syms C_cf C_cr C_lf C_lr;%分别为前后车轮的纵横向侧偏刚度，车辆固有参数
    Ccf=66900;Ccr=62700;Clf=66900;Clr=62700;
%syms m g I;%m为车辆质量，g为重力加速度，I为车辆绕Z轴的转动惯量，车辆固有参数
    m=1723;g=9.8;I=4175;
    w=3.75;
    %sf=34;%ssss是保持直线行驶的距离
    %sf=55.8;
    %%
    %重要参数
    sf=33.333;
    %
    d=78.329;
%  load('dlc.mat')
%% 参考轨迹生成
    shape=2.4;%参数名称，用于参考轨迹生成
    dx1=25;dx2=21.95;%没有任何实际意义，只是参数名称
    dy1=4.05;dy2=5.7;%没有任何实际意义，只是参数名称
    Xs1=27.19;Xs2=56.46;%参数名称
    X_predict=zeros(Np,1);%用于保存预测时域内的纵向位置信息，这是计算期望轨迹的基础
    phi_ref=zeros(Np,1);%用于保存预测时域内的期望轨迹
    
    %  以下计算kesi,即状态量与控制量合在一起   
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
    T=0.02;%仿真步长
%     T_all=20;%总的仿真时间，主要功能是防止计算期望轨迹越界
         %权重矩阵设置 
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
    %最基本也最重要的矩阵，是控制器的基础，采用动力学模型，该矩阵与车辆参数密切相关，通过对动力学方程求解雅克比矩阵得到
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
    d_k=zeros(Nx,1);%计算偏差
    state_k1=zeros(Nx,1);%预测的下一时刻状态量，用于计算偏差
    %以下即为根据离散非线性模型预测下一时刻状态量
    %注意，为避免前后轴距的表达式（a,b）与控制器的a,b矩冲突，将前后轴距的表达式改为lf和lr
    state_k1(1,1)=y_dot+T*(-x_dot*phi_dot+2*(Ccf*(delta_f-(y_dot+lf*phi_dot)/x_dot)+Ccr*(lr*phi_dot-y_dot)/x_dot)/m);
    state_k1(2,1)=x_dot+T*(y_dot*phi_dot+2*(Clf*Sf+Clr*Sr+Ccf*delta_f*(delta_f-(y_dot+phi_dot*lf)/x_dot))/m);
    state_k1(3,1)=phi+T*phi_dot;
    state_k1(4,1)=phi_dot+T*((2*lf*Ccf*(delta_f-(y_dot+lf*phi_dot)/x_dot)-2*lr*Ccr*(lr*phi_dot-y_dot)/x_dot)/I);
    state_k1(5,1)=Y+T*(x_dot*sin(phi)+y_dot*cos(phi));
    state_k1(6,1)=X+T*(x_dot*cos(phi)-y_dot*sin(phi));
    d_k=state_k1-a*kesi(1:6,1)-b*kesi(7,1);%根据falcone公式（2.11b）求得d(k,t)
    d_piao_k=zeros(Nx+Nu,1);%d_k的增广形式，参考falcone(B,4c)
    d_piao_k(1:6,1)=d_k;
    d_piao_k(7,1)=0;
    %% 侧偏角e(t)计算 这段主要用来后面的侧偏角软约束 用以实现车辆在极端工况下的安全行驶
    y_dot = state_k1(1,1);  %%  y_dot  x_dot  phi_dot主要用来更新 Alpha_fl r ..等侧偏角
    x_dot = state_k1(2,1);  %%  y_dot  x_dot  phi_dot主要用来更新 Alpha_fl r ..等侧偏角
    phi_dot = state_k1(4,1);%%  y_dot  x_dot  phi_dot主要用来更新 Alpha_fl r ..等侧偏角
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
 -1]; %% 侧偏角输出矩阵 D
C_d =  [ -(1.0*(cos(delta_f)/(cos(delta_f)*(phi_dot - 1.0*x_dot) - 1.0*sin(delta_f)*(1.232*phi_dot + y_dot)) + (sin(delta_f)*(cos(delta_f)*(1.232*phi_dot + y_dot) + sin(delta_f)*(phi_dot - 1.0*x_dot)))/(cos(delta_f)*(phi_dot - 1.0*x_dot) - 1.0*sin(delta_f)*(1.232*phi_dot + y_dot))^2))/((cos(delta_f)*(1.232*phi_dot + y_dot) + sin(delta_f)*(phi_dot - 1.0*x_dot))^2/(cos(delta_f)*(phi_dot - 1.0*x_dot) - 1.0*sin(delta_f)*(1.232*phi_dot + y_dot))^2 + 1.0),             (sin(delta_f)/(cos(delta_f)*(phi_dot - 1.0*x_dot) - 1.0*sin(delta_f)*(1.232*phi_dot + y_dot)) - (1.0*cos(delta_f)*(cos(delta_f)*(1.232*phi_dot + y_dot) + sin(delta_f)*(phi_dot - 1.0*x_dot)))/(cos(delta_f)*(phi_dot - 1.0*x_dot) - 1.0*sin(delta_f)*(1.232*phi_dot + y_dot))^2)/((cos(delta_f)*(1.232*phi_dot + y_dot) + sin(delta_f)*(phi_dot - 1.0*x_dot))^2/(cos(delta_f)*(phi_dot - 1.0*x_dot) - 1.0*sin(delta_f)*(1.232*phi_dot + y_dot))^2 + 1.0), 0,              -(1.0*((1.232*cos(delta_f) + sin(delta_f))/(cos(delta_f)*(phi_dot - 1.0*x_dot) - 1.0*sin(delta_f)*(1.232*phi_dot + y_dot)) - (1.0*(cos(delta_f)*(1.232*phi_dot + y_dot) + sin(delta_f)*(phi_dot - 1.0*x_dot))*(cos(delta_f) - 1.232*sin(delta_f)))/(cos(delta_f)*(phi_dot - 1.0*x_dot) - 1.0*sin(delta_f)*(1.232*phi_dot + y_dot))^2))/((cos(delta_f)*(1.232*phi_dot + y_dot) + sin(delta_f)*(phi_dot - 1.0*x_dot))^2/(cos(delta_f)*(phi_dot - 1.0*x_dot) - 1.0*sin(delta_f)*(1.232*phi_dot + y_dot))^2 + 1.0), 0, 0
        (cos(delta_f)/(sin(delta_f)*(1.232*phi_dot + y_dot) + cos(delta_f)*(phi_dot + x_dot)) - (sin(delta_f)*(1.0*cos(delta_f)*(1.232*phi_dot + y_dot) - sin(delta_f)*(phi_dot + x_dot)))/(sin(delta_f)*(1.232*phi_dot + y_dot) + cos(delta_f)*(phi_dot + x_dot))^2)/((1.0*cos(delta_f)*(1.232*phi_dot + y_dot) - sin(delta_f)*(phi_dot + x_dot))^2/(sin(delta_f)*(1.232*phi_dot + y_dot) + cos(delta_f)*(phi_dot + x_dot))^2 + 1.0),                              -(1.0*(sin(delta_f)/(sin(delta_f)*(1.232*phi_dot + y_dot) + cos(delta_f)*(phi_dot + x_dot)) + (1.0*cos(delta_f)*(1.0*cos(delta_f)*(1.232*phi_dot + y_dot) - sin(delta_f)*(phi_dot + x_dot)))/(sin(delta_f)*(1.232*phi_dot + y_dot) + cos(delta_f)*(phi_dot + x_dot))^2))/((1.0*cos(delta_f)*(1.232*phi_dot + y_dot) - sin(delta_f)*(phi_dot + x_dot))^2/(sin(delta_f)*(1.232*phi_dot + y_dot) + cos(delta_f)*(phi_dot + x_dot))^2 + 1.0), 0,                                             ((1.232*cos(delta_f) - 1.0*sin(delta_f))/(sin(delta_f)*(1.232*phi_dot + y_dot) + cos(delta_f)*(phi_dot + x_dot)) - ((1.0*cos(delta_f)*(1.232*phi_dot + y_dot) - sin(delta_f)*(phi_dot + x_dot))*(cos(delta_f) + 1.232*sin(delta_f)))/(sin(delta_f)*(1.232*phi_dot + y_dot) + cos(delta_f)*(phi_dot + x_dot))^2)/((1.0*cos(delta_f)*(1.232*phi_dot + y_dot) - sin(delta_f)*(phi_dot + x_dot))^2/(sin(delta_f)*(1.232*phi_dot + y_dot) + cos(delta_f)*(phi_dot + x_dot))^2 + 1.0), 0, 0
        -(1.0*(cos(delta_f)/(cos(delta_f)*(phi_dot - 1.0*x_dot) + sin(delta_f)*(1.468*phi_dot - 1.0*y_dot)) + (1.0*sin(delta_f)*(1.0*sin(delta_f)*(phi_dot - 1.0*x_dot) - cos(delta_f)*(1.468*phi_dot - 1.0*y_dot)))/(cos(delta_f)*(phi_dot - 1.0*x_dot) + sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2))/((1.0*sin(delta_f)*(phi_dot - 1.0*x_dot) - cos(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2/(cos(delta_f)*(phi_dot - 1.0*x_dot) + sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2 + 1.0), (sin(delta_f)/(cos(delta_f)*(phi_dot - 1.0*x_dot) + sin(delta_f)*(1.468*phi_dot - 1.0*y_dot)) - (cos(delta_f)*(1.0*sin(delta_f)*(phi_dot - 1.0*x_dot) - cos(delta_f)*(1.468*phi_dot - 1.0*y_dot)))/(cos(delta_f)*(phi_dot - 1.0*x_dot) + sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2)/((1.0*sin(delta_f)*(phi_dot - 1.0*x_dot) - cos(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2/(cos(delta_f)*(phi_dot - 1.0*x_dot) + sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2 + 1.0), 0, ((1.468*cos(delta_f) - 1.0*sin(delta_f))/(cos(delta_f)*(phi_dot - 1.0*x_dot) + sin(delta_f)*(1.468*phi_dot - 1.0*y_dot)) + (1.0*(1.0*sin(delta_f)*(phi_dot - 1.0*x_dot) - cos(delta_f)*(1.468*phi_dot - 1.0*y_dot))*(cos(delta_f) + 1.468*sin(delta_f)))/(cos(delta_f)*(phi_dot - 1.0*x_dot) + sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2)/((1.0*sin(delta_f)*(phi_dot - 1.0*x_dot) - cos(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2/(cos(delta_f)*(phi_dot - 1.0*x_dot) + sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2 + 1.0), 0, 0
        (1.0*(cos(delta_f)/(1.0*cos(delta_f)*(phi_dot + x_dot) - sin(delta_f)*(1.468*phi_dot - 1.0*y_dot)) + (1.0*sin(delta_f)*(sin(delta_f)*(phi_dot + x_dot) + cos(delta_f)*(1.468*phi_dot - 1.0*y_dot)))/(1.0*cos(delta_f)*(phi_dot + x_dot) - sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2))/((sin(delta_f)*(phi_dot + x_dot) + cos(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2/(1.0*cos(delta_f)*(phi_dot + x_dot) - sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2 + 1.0),                -(sin(delta_f)/(1.0*cos(delta_f)*(phi_dot + x_dot) - sin(delta_f)*(1.468*phi_dot - 1.0*y_dot)) - (cos(delta_f)*(sin(delta_f)*(phi_dot + x_dot) + cos(delta_f)*(1.468*phi_dot - 1.0*y_dot)))/(1.0*cos(delta_f)*(phi_dot + x_dot) - sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2)/((sin(delta_f)*(phi_dot + x_dot) + cos(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2/(1.0*cos(delta_f)*(phi_dot + x_dot) - sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2 + 1.0), 0,                        -((1.468*cos(delta_f) + sin(delta_f))/(1.0*cos(delta_f)*(phi_dot + x_dot) - sin(delta_f)*(1.468*phi_dot - 1.0*y_dot)) - ((sin(delta_f)*(phi_dot + x_dot) + cos(delta_f)*(1.468*phi_dot - 1.0*y_dot))*(cos(delta_f) - 1.468*sin(delta_f)))/(1.0*cos(delta_f)*(phi_dot + x_dot) - sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2)/((sin(delta_f)*(phi_dot + x_dot) + cos(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2/(1.0*cos(delta_f)*(phi_dot + x_dot) - sin(delta_f)*(1.468*phi_dot - 1.0*y_dot))^2 + 1.0), 0, 0
     ]; %% 侧偏角输出 = C_d 乘以状态量kesi因子 + D_d *delat_u + e_k
 
    CC_d_cell = cell(1,2);
    CC_d_cell{1, 1} = C_d;
    CC_d_cell{1, 2} = D_d;
    Cc_d = cell2mat(CC_d_cell);
    side_slip_error = state2_k1 - C_d*kesi(1:6,1) - D_d*kesi(7,1);%根据上面更新过来计算新的值 理论上要计算未来每个预测时域内的差值（最好不能强行理解为误差等，此处只是为了变换才这样表示，可看论文基于模型预测控制算法里面的推导公式）
    sideslip_piao_k=zeros(Ny_cepian, 1);%d_k的增广形式，参考falcone(B,4c)
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
        PHI_cell{p,1} = d_piao_k;%理论上来说，这个是要实时更新的，但是为了简便，这里又一次近似
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
    GAMMA=cell2mat(GAMMA_cell);%大写的GAMMA
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

               X_DOT=x_dot*cos(phi)-y_dot*sin(phi);%惯性坐标系下纵向速度,这是大地坐标系下的
               X_predict(Np,1)=X+X_DOT*p*T;%首先计算出未来X的位置，          
               x1=X+X_DOT*p*T-sf;
               x2=x1-d-10;
               Y_ref(p,1)=0.*(x1<0)+(10*w/(d^3)* x1^3-15*w/(d^4)* x1^4+6*w/(d^5)*x1^5).*((0<x1)&&(x1<=d))+3.75.*(x1>d);%期望换道路径
               phi_ref(p,1)=0.*(x1<0)+(atan(30*w/d^3* x1^2-60*w/d^4* x1^3+30*w/d^5*x1^4)).*((0<x1)&&(x1<=d))+0.*(x1>d);%期望横摆角  
               %Y_ref(p,1)=0.*(x1<0)+(10*w/(d^3)* x1^3-15*w/(d^4)* x1^4+6*w/(d^5)*x1^5).*((0<x1)&&(x1<=d))+3.75.*((d<x1)&&(x1<=d+10))+(3.75-((10*w/(d^3)* x2^3-15*w/(d^4)* x2^4+6*w/(d^5)*x2^5))).*((d+10<=x1)&&(x1<=2*d+10))+0.*(x1>=2*d+10);%期望换道路径
%                phi_ref(p,1)=0.*(x1<0)+(atan(30*w/d^3* x1^2-60*w/d^4* x1^3+30*w/d^5*x1^4)).*((0<x1)&&(x1<=d))+0.*(x1>d);%期望横摆角  
               Yita_ref_cell{p,1}=[phi_ref(p,1);Y_ref(p,1)];
    end     
    Yita_ref=cell2mat(Yita_ref_cell);
    error_1=-Yita_ref+PSI*kesi+GAMMA*PHI; %求偏差
    UU = ones(Nc,1)*U(1);
    f_cell=cell(1,2);
    f_cell{1,1}=2*error_1'*Q*THETA + 2*UU'*S*MM;
    f_cell{1,2}=0;
    f=cell2mat(f_cell);
    
 %% 以下为约束生成区域
 %控制量约束
    A_t=zeros(Nc,Nc);%见falcone论文 P181
    for p=1:1:Nc
        for q=1:1:Nc
            if q<=p 
                A_t(p,q)=1;
            else 
                A_t(p,q)=0;
            end
        end 
    end 
    A_I=kron(A_t,eye(Nu));%求克罗内克积
    Ut=kron(ones(Nc,1),U(1));
    umin = -0.1745;%维数与控制变量的个数相同
    umax = 0.1745;
    delta_umin=-0.014*1;
    delta_umax=0.014*1;
    Umin=kron(ones(Nc,1),umin);
    Umax=kron(ones(Nc,1),umax);
    
    
     %侧偏角输出软量约束
    y_cepian_max = [ 0.0489; 0.0489; 0.0489; 0.0489];%%%%%%%%%%%+-2.8°
    y_cepian_min = [ -0.0489; -0.0489; -0.0489; -0.0489];
    Y_cepian_max = kron(ones(Np,1),y_cepian_max);
    Y_cepian_min = kron(ones(Np,1),y_cepian_min);
    
    %横向y和航向输出量硬约束
    ycmax=[0.21;5];
    ycmin=[-0.3;-3];
    Ycmax=kron(ones(Np,1),ycmax);
    Ycmin=kron(ones(Np,1),ycmin);
    
    %结合控制量约束和输出量约束和软约束
    A_cons_cell={A_I zeros(Nu*Nc,1);-A_I zeros(Nu*Nc,1); THETA_sideslip -1*ones(Ny_cepian*Np,1);-THETA_sideslip ones(Ny_cepian*Np,1);THETA zeros(Ny*Np,1);-THETA zeros(Ny*Np,1)};
    b_cons_cell={Umax-Ut;-Umin+Ut; Y_cepian_max-PSI_sideslip*kesi-GAMMA_sideslip*PHI - PHI_sideslip; -Y_cepian_min+PSI_sideslip*kesi+GAMMA_sideslip*PHI + PHI_sideslip;Ycmax-PSI*kesi-GAMMA*PHI;-Ycmin+PSI*kesi+GAMMA*PHI};
    A_cons=cell2mat(A_cons_cell);%（求解方程）状态量不等式约束增益矩阵，转换为绝对值的取值范围
    b_cons=cell2mat(b_cons_cell);%（求解方程）状态量不等式约束的取值
       %状态量约束
    M=0.2*pi/180;   %软约束最大值 此处M是最大的软约束约束值，当此值太大容易车量失稳，失去加入软约束的意义，太小容易造成无解，控制丢失。
%       M=0;  %软约束最大值 此处M是最大的软约束约束值，当此值太大容易车量失稳，失去加入软约束的意义，太小容易造成无解，控制丢失。
    delta_Umin=kron(ones(Nc,1),delta_umin);
    delta_Umax=kron(ones(Nc,1),delta_umax);
    lb=[delta_Umin;0];%（求解方程）状态量下界，包含控制时域内控制增量和松弛因子
    ub=[delta_Umax;M];%（求解方程）状态量上界，包含控制时域内控制增量和松弛因子
    
    %% 开始求解过程
%        options = optimset('Algorithm','interior-point-convex');
options = optimset('Algorithm','active-set');
       x_start=zeros(Nc+1,1);%加入一个起始点
      [X,fval,exitflag]=quadprog(H,f,A_cons,b_cons,[],[],lb,ub,x_start,options);
      fprintf('exitflag=%d\n',exitflag); 
      fprintf('H=%4.2f\n',H(1,1));
      fprintf('f=%4.2f\n',f(1,1));
    %% 计算输出
    u_piao=X(1);%得到控制增量
    U(1)=kesi(7,1)+u_piao;%当前时刻的控制量为上一刻时刻控制+控制增量
   
   %U(2)=Yita_ref(2);%输出dphi_ref
    sys(1)= U(1);
    sys(2)=  phi_ref(p,1);
    sys(3)= Y_ref(p,1);
    sys(4)= X_predict(Np,1);
    sys(5)=u_piao;
    toc
% End of mdlOutputs.
  

   
