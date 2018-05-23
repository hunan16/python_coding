%增强型轨道炮  输入时序 输出电流波形、电枢速度波形  20161219
%发射试验 20161123-1  装弹位置 400mm

%7组  比较计算结果和实验结果

clear
clc
%****************参数******************************
num_Loop=100;
pos_min=0.01;
pos_max=0.4;
m_min=0.008;
m_max=0.015;
for ik=1:num_Loop
N=[9 3 4 3 4 4 3];
Trig=[0 280 400 600 760 960 1160]*10^-6*rand(1);  %确定维数  触发通道数确定 修改后可以扩充
ZDWZ=pos_min+(pos_max-pos_min)*rand(1);         %装弹位置 
m=m_min+(m_max-m_min)*rand(1);                                      %电枢质量 kg
%PFU参数----约束参数-----轨道炮负载参数---电路模型求解参数
E_total=4;                        %总能量 单位MJ
Imax=800*10^3;              %imax  导轨电流最大阈值  单位A
Imin=700*10^3;               %imin   导轨电流最低阈值  单位A

n_total=40 ;                    %PFU总单元数
% U0=round((E_total/n_total*10^3)^0.5*10^3*1000) /1000;                   %PFU电容预充电压 单位V 总能量一定
U0=2.8*10^3;
C=2*10^-3 ;                    %PFU储能电容  单位F
L=13*10^-6 ;                   %PFU脉冲调波电感L  单位 H
R=7*10^-3;                     %PFU电感电阻与引线电阻之和 单位 欧姆

dt=10^-7;                          % 微小时间段步长 0.1us
N_u_step=1;                    %电压接入 分段步数 0
N_capacity=10^6;            %矩阵维数 

Lr_initial=1.5*10^-6;             %导轨初始电感   ！！！为了让程序不出尖刺，缩小微分作用的影响 电感接入 分段步数 
N_Lr0=50;                             %电感接入 分段步数 

%连接发射器汇流排和电源汇流排的粗电缆    1123-2后2根增加至6根
R_cabel=3.8*0.001*1/2;  %单根电缆 除以根数 并联
L_cabel=0.9*10^-6*1/2;   %单根电缆
%发射器参数
Length_rail=2500*10^-3;          %轨道长度，用于计算电阻

Length=Length_rail-ZDWZ;        %加速距离 更改！！
height=0.02;                     %导轨高度9cm
b=0.03;                             %导轨间距9cm

Parm=26/3*10^-3;          %轨道与电枢的接触长度 15g电枢
A_arm=11*28*10^-6;      %电枢的横截面积
K1=0.7;                           %电枢径向力与轴向力之比
uf1=0.9;                           %电枢与轨道 铝和铜间的静摩擦系数
uf2=0.3;                            %电枢与轨道 铝和铜间的动摩擦系数   后面改为 0.02
lc=(11+28)*2*10^-3;
FN0=20*10^3;


Permeability=1.26*10^-6;      %导轨材料磁导率1.26微亨/米  μ
L_gradient=1.1*10^-6;          %导轨电感梯度0.8微亨/米
resistivity10=1.75*10^-8;          %导轨初始等效电阻率   17.5毫欧*m/mm2    ρ10
resistivity20=2.83*10^-8;          %电枢初始等效电阻率 米*欧姆 ρ20
bert=3.6*10^-16;                       %电阻率热效应常数 β   απηρ
kc=0.9;                                      %初始接触电阻常数
Rvc=10^-9;                                %接触速度趋肤效应电阻常数
kVSEC=0.5;                              %接触速度趋肤效应电阻分配比例
% RP=10^-3;                                  %炮口熄弧电阻


af=R/2/L;                   %衰减常数  方程复根 实部
w=(1/L/C-af^2)^0.5;  %方程复根 虚部  振荡频率
tao=L/R;                    %一阶方程时间常数 

%****************参数书写完毕******************************
%  matrix_Isource=zeros(1,10000);  % 导轨电流矩阵 每个阶段都会改变
matrixu_C1=zeros(1,100000);
matrixu_C2=zeros(1,100000);
matrixu_C3=zeros(1,100000);
matrixu_C4=zeros(1,100000);
matrixu_C5=zeros(1,100000);
matrixu_C6=zeros(1,100000);
matrixu_C7=zeros(1,100000);
matrixI_L1 =zeros(1,100000);
matrixI_L2=zeros(1,100000);
matrixI_L3=zeros(1,100000);
matrixI_L4 =zeros(1,100000);
matrixI_L5 =zeros(1,100000);
matrixI_L6 =zeros(1,100000);
matrixI_L7 =zeros(1,100000);
matrixu_load=zeros(1,100000);
matrix_Isource=zeros(1,100000);
matrix_Ia=zeros(1,100000);
matrix_v=zeros(1,100000);
matrix_x=zeros(1,100000);

%初始状态 不改变认为这一组没有模块触发～
flag1=0;
flag2=0;
flag3=0;
flag4=0;
flag5=0;
flag6=0; 
flag7=0;

t=Trig(1);
j=1;
flag1=1;

u_C0=U0;      %电容初始电压 10kV

Lr0=0;    %!!! 逐渐加电感
I_L1=0;          %电感电流初始值为0
u_load=0;        %负载两端电压 
%负载初始状态 
v=0;
x=0;
Isource=0;
Ia=0;
ku1=1;   
while(~(t>Trig(2)))
        matrixu_C1(j)=u_C0;
        matrixI_L1(j) = I_L1;
        matrixu_load(j)=u_load;
        matrix_v(j)=v;
        matrix_x(j)=x;
        matrix_Isource(j)=Isource;
        matrix_Ia(j)=Ia;
        
        if(flag1==1)
        u_C1=u_load+(u_C0-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(u_C0-u_load)-I_L1/C)*exp(-af*dt)*sin(w*dt);
        I_L1=I_L1*exp(-af*dt)*cos(w*dt)+1/w*((u_C0-u_load)/L-af*I_L1)*exp(-af*dt)*sin(w*dt);   
        if(u_C1<1)
            flag1=2;
            u_C1=0;
        end
        else if(flag1==2)
            I_L1=(u_load/R+I_L1)*exp(-dt/tao)-u_load/R;
        if(I_L1<0)
            flag1=3;
            I_L1=0;
        end         
        end
        end
        Isource=N(1)*I_L1;
          
        t=t+dt;            %!!! 考虑分母不为零  v 在时间端点上跳变，电阻公式里t采用“t+dt” 而v采用“t时刻的v”！！！
        
        %负载等效电路模型求解
       resistivity1=resistivity10+bert*Isource/height;
       resistivity2=resistivity20+bert*Ia/height;
       %导轨电阻 外轨和内轨 
       Rra=8*Length_rail/3/height*(Permeability*resistivity1/2/pi/t)^0.5;
       Rrb=8*x/3/height*(Permeability*resistivity1/2/pi/t)^0.5+8*ZDWZ/3/height*(Permeability*resistivity1/2/pi/t)^0.5;
       Rr=Rra+Rrb;
       Rc=kc/height*(Permeability*(resistivity1+resistivity2)/4/pi/t)^0.5;
       RVSEC=Rvc*v^1.5;
       RVSEC1=(1-kVSEC)*RVSEC;
       RVSEC2=kVSEC*RVSEC;
       Ra=b/height*(Permeability*resistivity2/t/pi/2)^0.5;
       uEMF=L_gradient*v*Isource;
        
       if(j>1)
       di=Isource-matrix_Isource(j-1);
             else
        di=(Isource-matrix_Isource(j))*2;    
       end
       
        if(ku1>N_Lr0)
           Lr0=Lr0;
       else
           Lr0=Lr0+Lr_initial/N_Lr0;
        end
       if(ku1>N_Lr0)
           Lr0=Lr0;
       else
           Lr0=Lr0+Lr_initial/N_Lr0;
       end
       Lr=Lr0+L_gradient*x; 
%    u1_load=(Rr+Rc+RVSEC1+(RVSEC2+Ra)*RP/(RVSEC2+Ra+RP)+L_gradient*v)*Isource;%+di/dt*Lr;       
       u1_load=(Rr+Rc+RVSEC1+RVSEC2+Ra+L_gradient*v)*Isource;%+di/dt*Lr;
       u2_load=R_cabel*Isource+L_cabel*di/dt;
       u3_load=di/dt/2*Lr;
       u_load=u1_load+u2_load+u3_load; 
       
       if(u_load>3000)
           u_load=1600;%-u_load;  % 遏制振荡   从它变为负数开始 强制校正为正 
       end
       if(u_load<0)
           u_load=50;%-u_load;  % 遏制振荡   从它变为负数开始 强制校正为正 
       end
%     Ia=RP/(RVSEC2+Ra+RP)*Isource;
       F=L_gradient*Isource^2/2;
       
       FN=FN0+K1*Parm*lc/4/A_arm*F;
       Ff_staic=FN*uf1;
       
       if(v>0)
           Ff=uf2*FN;
       else
           Ff=Ff_staic;
       end
       if(F>Ff)
           a=(F-Ff)/m;
       else
           a=0;
       end
       
       v=v+a*dt;
       x=x+matrix_v(j)*dt+1/2*a*(dt^2+2*(t-dt)*dt);
      
        u_C0=u_C1;
               
        j=j+1;
        ku1=ku1+1;
        end
        
 t=t-dt;
 j=j-1;
 flag2=1;

  u_C0=U0;      %电容初始电压
  
  Lr0=0;    %!!! 逐渐加电感
  I_L2=0;          %电感电流初始值为0
  u_load=matrixu_load(j);        %负载两端电压 
  %负载初始状态 
 v=matrix_v(j);
 x=matrix_x(j);
 Isource=matrix_Isource(j);
 Ia= matrix_Ia(j);
 I_L1=matrixI_L1(j) ; 
 u_C1=matrixu_C1(j);
 ku2=1;
while(~(t>Trig(3)))
        matrixu_C2(j)=u_C0;
        matrixI_L2(j) = I_L2;
        matrixu_load(j)=u_load;
        matrix_v(j)=v;
        matrix_x(j)=x;
        matrix_Isource(j)=Isource;
        matrix_Ia(j)=Ia;
        
        matrixI_L1(j)=I_L1;
        matrixu_C1(j)=u_C1;
        
        if(flag2==1)
             u_C2=u_load+(u_C0-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(u_C0-u_load)-I_L2/C)*exp(-af*dt)*sin(w*dt);
             I_L2=I_L2*exp(-af*dt)*cos(w*dt)+1/w*((u_C0-u_load)/L-af*I_L2)*exp(-af*dt)*sin(w*dt);
        if(u_C2<0)
            flag2=2;
            u_C2=0;
            I_L2=matrixI_L2(j) ;
        end
        else  if(flag2==2)
                   I_L2=(u_load/R+I_L2)*exp(-dt/tao)-u_load/R;
                   if(I_L2<0)
                       flag2=3;
                       I_L2=0;
                   end
            end
        end
     
        
     if(flag1==1)
        m1=u_C1;  %middle variable
        u_C1=u_load+(m1-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m1-u_load)-I_L1/C)*exp(-af*dt)*sin(w*dt);
        I_L1=I_L1*exp(-af*dt)*cos(w*dt)+1/w*((m1-u_load)/L-af*I_L1)*exp(-af*dt)*sin(w*dt);   
        if(u_C1<1)
            flag1=2;
            u_C1=0;
        end
        else if(flag1==2)
            I_L1=(u_load/R+I_L1)*exp(-dt/tao)-u_load/R;
        if(I_L1<0)
            flag1=3;
            I_L1=0;
        end         
     end
     end
        
        Isource=N(2)*I_L2+N(1)*I_L1 ;    %n个模块  第一次计算时的公式！！！ 历史电流？
   
        t=t+dt;            
        
        %负载等效电路模型求解
       resistivity1=resistivity10+bert*Isource/height;
       resistivity2=resistivity20+bert*Ia/height;
      %导轨电阻 外轨和内轨 
       Rra=8*Length_rail/3/height*(Permeability*resistivity1/2/pi/t)^0.5;
       Rrb=8*x/3/height*(Permeability*resistivity1/2/pi/t)^0.5+8*ZDWZ/3/height*(Permeability*resistivity1/2/pi/t)^0.5;
       Rr=Rra+Rrb;
       Rc=kc/height*(Permeability*(resistivity1+resistivity2)/4/pi/t)^0.5;
       RVSEC=Rvc*v^1.5;
       RVSEC1=(1-kVSEC)*RVSEC;
       RVSEC2=kVSEC*RVSEC;
       Ra=b/height*(Permeability*resistivity2/t/pi/2)^0.5;
       uEMF=L_gradient*v*Isource;
       
       if(j>1)
       di=Isource-matrix_Isource(j-1);
             else
        di=(Isource-matrix_Isource(j))*2;    
       end
            
        if(ku2>N_Lr0)
           Lr0=Lr0;
       else
           Lr0=Lr0+Lr_initial/N_Lr0;
        end
       
       Lr=Lr0+L_gradient*x;   
       
       u1_load=(Rr+Rc+RVSEC1+RVSEC2+Ra+L_gradient*v)*Isource;%+di/dt*Lr;
       u2_load=R_cabel*Isource+L_cabel*di/dt;
       u3_load=di/dt/2*Lr;
       u_load=u1_load+u2_load+u3_load; 
       
       if(u_load>4000)
           u_load=2000;%-u_load;  % 遏制振荡   从它变为负数开始 强制校正为正 
       end
       if(u_load<0)
            u_load=50;%-u_load;  % 遏制振荡   从它变为负数开始 强制校正为正 
       end
      
      F=L_gradient*Isource^2/2;
      if(v>0)
           Ff=uf2*FN;
       else
           Ff=Ff_staic;
       end
       if(F>Ff)
           a=(F-Ff)/m;
       else
           a=0;
       end
       
       v=v+a*dt;
       x=x+matrix_v(j)*dt+1/2*a*(dt^2+2*(t-dt)*dt);
       
       if(x>Length)
           break
       end
        u_C0=u_C2;
              
        j=j+1;
        ku2=ku2+1;
end

%第三组放电计算
 t=t-dt;
 j=j-1;
 flag3=1;

  u_C0=U0;      %电容初始电压

  Lr0=0;    %!!! 逐渐加电感
  I_L3=0;          %电感电流初始值为0
  u_load=matrixu_load(j);        %负载两端电压 
  %负载初始状态 
 v=matrix_v(j);
 x=matrix_x(j);
 Isource=matrix_Isource(j);
 Ia= matrix_Ia(j);
 I_L1=matrixI_L1(j) ; 
 u_C1=matrixu_C1(j);
 I_L2=matrixI_L2(j);
u_C2=matrixu_C2(j);
ku3=1;
while(~(t>Trig(4)))
        matrixu_C3(j)=u_C0;
        matrixI_L3(j) = I_L3;
        matrixu_load(j)=u_load;
        matrix_v(j)=v;
        matrix_x(j)=x;
        matrix_Isource(j)=Isource;
        matrix_Ia(j)=Ia;
        
        matrixI_L1(j)=I_L1;
        u_C1=matrixu_C1(j);
        matrixI_L2(j)= I_L2;
        matrixu_C2(j)=u_C2;
        
        if(flag3==1)
             u_C3=u_load+(u_C0-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(u_C0-u_load)-I_L3/C)*exp(-af*dt)*sin(w*dt);
             I_L3=I_L3*exp(-af*dt)*cos(w*dt)+1/w*((u_C0-u_load)/L-af*I_L3)*exp(-af*dt)*sin(w*dt);
        if(u_C3<0)
            flag3=2;
            u_C3=0;
            I_L3=matrixI_L3(j) ;
        end
        else  if(flag3==2)
                   I_L3=(u_load/R+I_L3)*exp(-dt/tao)-u_load/R;
                   if(I_L3<0)
                       flag3=3;
                       I_L3=0;
                   end
            end
        end
        
        if(flag2==1)
            m2=u_C2;
             u_C2=u_load+(m2-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m2-u_load)-I_L2/C)*exp(-af*dt)*sin(w*dt);
             I_L2=I_L2*exp(-af*dt)*cos(w*dt)+1/w*((m2-u_load)/L-af*I_L2)*exp(-af*dt)*sin(w*dt);
        if(u_C2<0)
            flag2=2;
            u_C2=0;
            I_L2=matrixI_L2(j) ;
        end
        else  if(flag2==2)
                   I_L2=(u_load/R+I_L2)*exp(-dt/tao)-u_load/R;
                   if(I_L2<0)
                       flag2=3;
                       I_L2=0;
                   end
            end
        end
     
        
     if(flag1==1)
        m1=u_C1;  %middle variable
        u_C1=u_load+(m1-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m1-u_load)-I_L1/C)*exp(-af*dt)*sin(w*dt);
        I_L1=I_L1*exp(-af*dt)*cos(w*dt)+1/w*((m1-u_load)/L-af*I_L1)*exp(-af*dt)*sin(w*dt);   
        if(u_C1<1)
            flag1=2;
            u_C1=0;
        end
        else if(flag1==2)
            I_L1=(u_load/R+I_L1)*exp(-dt/tao)-u_load/R;
        if(I_L1<0)
            flag1=3;
            I_L1=0;
        end         
     end
     end
        
        Isource=N(2)*I_L2+N(1)*I_L1+N(3)*I_L3;    %n个模块  第一次计算时的公式！！！ 历史电流？
   
        t=t+dt;            
        
        %负载等效电路模型求解
       resistivity1=resistivity10+bert*Isource/height;
       resistivity2=resistivity20+bert*Ia/height;
      %导轨电阻 外轨和内轨 
       Rra=8*Length_rail/3/height*(Permeability*resistivity1/2/pi/t)^0.5;
       Rrb=8*x/3/height*(Permeability*resistivity1/2/pi/t)^0.5+8*ZDWZ/3/height*(Permeability*resistivity1/2/pi/t)^0.5;
       Rr=Rra+Rrb;
       
       Rc=kc/height*(Permeability*(resistivity1+resistivity2)/4/pi/t)^0.5;
       RVSEC=Rvc*v^1.5;
       RVSEC1=(1-kVSEC)*RVSEC;
       RVSEC2=kVSEC*RVSEC;
       Ra=b/height*(Permeability*resistivity2/t/pi/2)^0.5;
       uEMF=L_gradient*v*Isource;
        
       if(j>1)
       di=Isource-matrix_Isource(j-1);
             else
        di=(Isource-matrix_Isource(j))*2;    
       end
       
        if(ku3>N_Lr0)
           Lr0=Lr0;
       else
           Lr0=Lr0+Lr_initial/N_Lr0;
       end
       Lr=Lr0+L_gradient*x;   
%        u_load=(Rr+Rc+RVSEC1+(RVSEC2+Ra)*RP/(RVSEC2+Ra+RP)+L_gradient*v)*Isource;%+di/dt*Lr;       
%        Ia=RP/(RVSEC2+Ra+RP)*Isource;
       u1_load=(2*Rr+Rc+RVSEC1+RVSEC2+Ra+L_gradient*v)*Isource;%+di/dt*Lr;
       u2_load=R_cabel*Isource+L_cabel*di/dt;
       u3_load=di/dt/2*Lr;   %电感有两项  另外一项是在u1-load里 电感梯度*速度*电流
       u_load=u1_load+u2_load+u3_load; 
      
       if(u_load>2000)
           u_load=1000;%-u_load;  % 遏制振荡   从它变为负数开始 强制校正为正 
       end
       if(u_load<0)
           u_load=50;%-u_load;  % 遏制振荡   从它变为负数开始 强制校正为正 
       end    
      
      F=L_gradient*Isource^2/2;
      if(v>0)
           Ff=uf2*FN;
       else
           Ff=Ff_staic;
       end
       if(F>Ff)
           a=(F-Ff)/m;
       else
           a=0;
       end
       
       v=v+a*dt;
       x=x+matrix_v(j)*dt+1/2*a*(dt^2+2*(t-dt)*dt);
       if(x>Length)
           break
       end
        u_C0=u_C3;
              
        j=j+1;
        ku3=ku3+1;
end
%*************
L=13*10^-6 ;                    %PFU脉冲调波电感L  单位 H
R=12*10^-3;                     %PFU电感电阻与引线电阻之和 单位 欧姆
%****************
%第四组放电计算
 t=t-dt;
 j=j-1;
 flag4=1;

  u_C0=U0;      %电容初始电压
  
  Lr0=0;    %!!! 逐渐加电感
  I_L4=0;          %电感电流初始值为0
  u_load=matrixu_load(j);        %负载两端电压 
  %负载初始状态 
 v=matrix_v(j);
 x=matrix_x(j);
 Isource=matrix_Isource(j);
 Ia= matrix_Ia(j);
 I_L1=matrixI_L1(j) ; 
 u_C1=matrixu_C1(j);
 I_L2=matrixI_L2(j);
u_C2=matrixu_C2(j);
 I_L3=matrixI_L3(j);
u_C3=matrixu_C3(j);
ku4=1;
while(~(t>Trig(5)))
        matrixu_C4(j)=u_C0;
        matrixI_L4(j) = I_L4;
        matrixu_load(j)=u_load;
        matrix_v(j)=v;
        matrix_x(j)=x;
        matrix_Isource(j)=Isource;
        matrix_Ia(j)=Ia;
        
        matrixI_L1(j)=I_L1;
        matrixu_C1(j)=u_C1;
        matrixI_L2(j)= I_L2;
        matrixu_C2(j)=u_C2;
        matrixI_L3(j)= I_L3;
        matrixu_C3(j)=u_C3;
        
        
        if(flag4==1)
             u_C4=u_load+(u_C0-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(u_C0-u_load)-I_L4/C)*exp(-af*dt)*sin(w*dt);
             I_L4=I_L4*exp(-af*dt)*cos(w*dt)+1/w*((u_C0-u_load)/L-af*I_L4)*exp(-af*dt)*sin(w*dt);
        if(u_C4<0)
            flag4=2;
            u_C4=0;
            I_L4=matrixI_L4(j) ;
        end
        else if(flag4==2)
                   I_L4=(u_load/R+I_L4)*exp(-dt/tao)-u_load/R;
                   if(I_L4<0)
                       flag4=3;
                       I_L4=0;
                   end
            end
        end
        
           if(flag3==1)
             m3=u_C3;
             u_C3=u_load+(m3-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m3-u_load)-I_L3/C)*exp(-af*dt)*sin(w*dt);
             I_L3=I_L3*exp(-af*dt)*cos(w*dt)+1/w*((m3-u_load)/L-af*I_L3)*exp(-af*dt)*sin(w*dt);
        if(u_C3<0)
            flag3=2;
            u_C3=0;
            I_L3=matrixI_L3(j) ;
        end
       else if(flag3==2)
                   I_L3=(u_load/R+I_L3)*exp(-dt/tao)-u_load/R;
                   if(I_L3<0)
                       flag3=3;
                       I_L3=0;
                   end
           end
           end
        
        if(flag2==1)
            m2=u_C2;
             u_C2=u_load+(m2-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m2-u_load)-I_L2/C)*exp(-af*dt)*sin(w*dt);
             I_L2=I_L2*exp(-af*dt)*cos(w*dt)+1/w*((m2-u_load)/L-af*I_L2)*exp(-af*dt)*sin(w*dt);
        if(u_C2<0)
            flag2=2;
            u_C2=0;
            I_L2=matrixI_L2(j) ;
        end
        else  if(flag2==2)
                   I_L2=(u_load/R+I_L2)*exp(-dt/tao)-u_load/R;
                   if(I_L2<0)
                       flag2=3;
                       I_L2=0;
                   end
            end
        end
     
        
     if(flag1==1)
        m1=u_C1;  %middle variable
        u_C1=u_load+(m1-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m1-u_load)-I_L1/C)*exp(-af*dt)*sin(w*dt);
        I_L1=I_L1*exp(-af*dt)*cos(w*dt)+1/w*((m1-u_load)/L-af*I_L1)*exp(-af*dt)*sin(w*dt);   
        if(u_C1<1)
            flag1=2;
            u_C1=0;
        end
        else if(flag1==2)
            I_L1=(u_load/R+I_L1)*exp(-dt/tao)-u_load/R;
        if(I_L1<0)
            flag1=3;
            I_L1=0;
        end         
     end
     end
        
        Isource=N(2)*I_L2+N(1)*I_L1+N(3)*I_L3+N(4)*I_L4;    %n个模块  第一次计算时的公式！！！ 历史电流？
   
        t=t+dt;            
        
        %负载等效电路模型求解
       resistivity1=resistivity10+bert*Isource/height;
       resistivity2=resistivity20+bert*Ia/height;
      %导轨电阻 外轨和内轨 
       Rra=8*Length_rail/3/height*(Permeability*resistivity1/2/pi/t)^0.5;
       Rrb=8*x/3/height*(Permeability*resistivity1/2/pi/t)^0.5+8*ZDWZ/3/height*(Permeability*resistivity1/2/pi/t)^0.5;
       Rr=Rra+Rrb;
       
       Rc=kc/height*(Permeability*(resistivity1+resistivity2)/4/pi/t)^0.5;
       RVSEC=Rvc*v^1.5;
       RVSEC1=(1-kVSEC)*RVSEC;
       RVSEC2=kVSEC*RVSEC;
       Ra=b/height*(Permeability*resistivity2/t/pi/2)^0.5;
       uEMF=L_gradient*v*Isource;
        
       if(j>1)
       di=Isource-matrix_Isource(j-1);
             else
        di=(Isource-matrix_Isource(j))*2;    
       end
       
        if(ku4>N_Lr0)
           Lr0=Lr0;
       else
           Lr0=Lr0+Lr_initial/N_Lr0;
       end
       Lr=Lr0+L_gradient*x;   
%        u_load=(Rr+Rc+RVSEC1+(RVSEC2+Ra)*RP/(RVSEC2+Ra+RP)+L_gradient*v)*Isource;%+di/dt*Lr;       
%        Ia=RP/(RVSEC2+Ra+RP)*Isource;
       u1_load=(2*Rr+Rc+RVSEC1+RVSEC2+Ra+L_gradient*v)*Isource;%+di/dt*Lr;
       u2_load=R_cabel*Isource+L_cabel*di/dt;
       u3_load=di/dt/2*Lr;   
       u_load=u1_load+u2_load+u3_load;
       
       if(u_load>3000)
           u_load=2000;%-u_load;  % 遏制振荡   从它变为负数开始 强制校正为正 
       end
       if(u_load<0)
           u_load=60;%-u_load;  % 遏制振荡   从它变为负数开始 强制校正为正 
       end
       
      F=L_gradient*Isource^2/2;
      if(v>0)
           Ff=uf2*FN;
       else
           Ff=Ff_staic;
       end
       if(F>Ff)
           a=(F-Ff)/m;
       else
           a=0;
       end
       v=v+a*dt;
       x=x+matrix_v(j)*dt+1/2*a*(dt^2+2*(t-dt)*dt);
       
       if(x>Length)
           break
       end
        u_C0=u_C4;
              
        j=j+1;
        ku4=ku4+1;
end

%修改预紧力 和摩擦系数
FN0=10*10^3;   
uf2=0.02;
%第五组放电计算
 t=t-dt;
 j=j-1;
 flag5=1;

  u_C0=U0;      %电容初始电压
  
  Lr0=0;    %!!! 逐渐加电感
  I_L5=0;          %电感电流初始值为0
  u_load=matrixu_load(j);        %负载两端电压 
  %负载初始状态 
 v=matrix_v(j);
 x=matrix_x(j);
 Isource=matrix_Isource(j);
 Ia= matrix_Ia(j);
 
 I_L1=matrixI_L1(j) ; 
 u_C1=matrixu_C1(j);
 I_L2=matrixI_L2(j);
u_C2=matrixu_C2(j);
 I_L3=matrixI_L3(j);
u_C3=matrixu_C3(j);
I_L4=matrixI_L4(j);
u_C4=matrixu_C4(j);
ku5=1;
while(~(t>Trig(6)))
        matrixu_C5(j)=u_C0;
        matrixI_L5(j) = I_L5;
        matrixu_load(j)=u_load;
        matrix_v(j)=v;
        matrix_x(j)=x;
        matrix_Isource(j)=Isource;
        matrix_Ia(j)=Ia;
        
        matrixI_L1(j)=I_L1;
        matrixu_C1(j)=u_C1;
        matrixI_L2(j)= I_L2;
        matrixu_C2(j)=u_C2;
        matrixI_L3(j)= I_L3;
        matrixu_C3(j)=u_C3;
        matrixI_L3(j)= I_L3;
        matrixu_C3(j)=u_C3;
         matrixI_L4(j)= I_L4;
        matrixu_C4(j)=u_C4;
        
        if(flag5==1)
             u_C5=u_load+(u_C0-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(u_C0-u_load)-I_L5/C)*exp(-af*dt)*sin(w*dt);
             I_L5=I_L5*exp(-af*dt)*cos(w*dt)+1/w*((u_C0-u_load)/L-af*I_L5)*exp(-af*dt)*sin(w*dt);
        if(u_C5<0)
            flag5=2;
            u_C5=0;
            I_L5=matrixI_L5(j) ;
        end
        else if(flag5==2)
                   I_L5=(u_load/R+I_L5)*exp(-dt/tao)-u_load/R;
                   if(I_L5<0)
                       flag5=3;
                       I_L5=0;
                   end
            end
       end
        
        if(flag4==1)
            m4=u_C4;
             u_C4=u_load+(m4-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m4-u_load)-I_L4/C)*exp(-af*dt)*sin(w*dt);
             I_L4=I_L4*exp(-af*dt)*cos(w*dt)+1/w*((m4-u_load)/L-af*I_L4)*exp(-af*dt)*sin(w*dt);
        if(u_C4<0)
            flag4=2;
            u_C4=0;
            I_L4=matrixI_L4(j) ;
        end
        else if(flag4==2)
                   I_L4=(u_load/R+I_L4)*exp(-dt/tao)-u_load/R;
                   if(I_L4<0)
                       flag4=3;
                       I_L4=0;
                   end
            end
        end
        
           if(flag3==1)
             m3=u_C3;
             u_C3=u_load+(m3-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m3-u_load)-I_L3/C)*exp(-af*dt)*sin(w*dt);
             I_L3=I_L3*exp(-af*dt)*cos(w*dt)+1/w*((m3-u_load)/L-af*I_L3)*exp(-af*dt)*sin(w*dt);
        if(u_C3<0)
            flag3=2;
            u_C3=0;
            I_L3=matrixI_L3(j) ;
        end
       else if(flag3==2)
                   I_L3=(u_load/R+I_L3)*exp(-dt/tao)-u_load/R;
                   if(I_L3<0)
                       flag3=3;
                       I_L3=0;
                   end
           end
           end
        
        if(flag2==1)
            m2=u_C2;
             u_C2=u_load+(m2-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m2-u_load)-I_L2/C)*exp(-af*dt)*sin(w*dt);
             I_L2=I_L2*exp(-af*dt)*cos(w*dt)+1/w*((m2-u_load)/L-af*I_L2)*exp(-af*dt)*sin(w*dt);
        if(u_C2<0)
            flag2=2;
            u_C2=0;
            I_L2=matrixI_L2(j) ;
        end
        else  if(flag2==2)
                   I_L2=(u_load/R+I_L2)*exp(-dt/tao)-u_load/R;
                   if(I_L2<0)
                       flag2=3;
                       I_L2=0;
                   end
            end
        end
     
        
     if(flag1==1)
        m1=u_C1;  %middle variable
        u_C1=u_load+(m1-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m1-u_load)-I_L1/C)*exp(-af*dt)*sin(w*dt);
        I_L1=I_L1*exp(-af*dt)*cos(w*dt)+1/w*((m1-u_load)/L-af*I_L1)*exp(-af*dt)*sin(w*dt);   
        if(u_C1<1)
            flag1=2;
            u_C1=0;
        end
        else if(flag1==2)
            I_L1=(u_load/R+I_L1)*exp(-dt/tao)-u_load/R;
        if(I_L1<0)
            flag1=3;
            I_L1=0;
        end         
     end
     end
        
        Isource=N(2)*I_L2+N(1)*I_L1+N(3)*I_L3+N(4)*I_L4+N(5)*I_L5;    %n个模块  第一次计算时的公式！！！ 历史电流？
   
        t=t+dt;            
        
        %负载等效电路模型求解
       resistivity1=resistivity10+bert*Isource/height;
       resistivity2=resistivity20+bert*Ia/height;
      %导轨电阻 外轨和内轨 
       Rra=8*Length_rail/3/height*(Permeability*resistivity1/2/pi/t)^0.5;
       Rrb=8*x/3/height*(Permeability*resistivity1/2/pi/t)^0.5+8*ZDWZ/3/height*(Permeability*resistivity1/2/pi/t)^0.5;
       Rr=Rra+Rrb;
       Rc=kc/height*(Permeability*(resistivity1+resistivity2)/4/pi/t)^0.5;
       RVSEC=Rvc*v^1.5;
       RVSEC1=(1-kVSEC)*RVSEC;
       RVSEC2=kVSEC*RVSEC;
       Ra=b/height*(Permeability*resistivity2/t/pi/2)^0.5;
       uEMF=L_gradient*v*Isource;
        
       if(j>1)
       di=Isource-matrix_Isource(j-1);
             else
        di=(Isource-matrix_Isource(j))*2;    
       end
       
        if(ku5>N_Lr0)
           Lr0=Lr0;
       else
           Lr0=Lr0+Lr_initial/N_Lr0;
       end
       Lr=Lr0+L_gradient*x;   
       u1_load=(Rr+Rc+RVSEC1+RVSEC2+Ra+L_gradient*v)*Isource;%+di/dt*Lr;
       u2_load=R_cabel*Isource+L_cabel*di/dt;
       u3_load=di/dt/2*Lr;   
       u_load=u1_load+u2_load+u3_load; 
       
       if(u_load>2000)
           u_load=1000;%-u_load;  % 遏制振荡   从它变为负数开始 强制校正为正 
       end
       if(u_load<0)
           u_load=60;%-u_load;  % 遏制振荡   从它变为负数开始 强制校正为正 
       end
        
      F=L_gradient*Isource^2/2;
      if(v>0)
           Ff=uf2*FN;
       else
           Ff=Ff_staic;
       end
       if(F>Ff)
           a=(F-Ff)/m;
       else
           a=0;
       end
       v=v+a*dt;
       x=x+matrix_v(j)*dt+1/2*a*(dt^2+2*(t-dt)*dt);
       if(x>Length)
           break
       end
        u_C0=u_C5;
              
        j=j+1;
        ku5=ku5+1;
end
%第六组放电计算
 t=t-dt;
 j=j-1;
 flag6=1;

  u_C0=U0;      %电容初始电压
  
  Lr0=0;    %!!! 逐渐加电感
  I_L6=0;          %电感电流初始值为0
  u_load=matrixu_load(j);        %负载两端电压 
  %负载初始状态 
 v=matrix_v(j);
 x=matrix_x(j);
 Isource=matrix_Isource(j);
 Ia= matrix_Ia(j);
 
 I_L1=matrixI_L1(j) ; 
 u_C1=matrixu_C1(j);
 I_L2=matrixI_L2(j);
u_C2=matrixu_C2(j);
 I_L3=matrixI_L3(j);
u_C3=matrixu_C3(j);
I_L4=matrixI_L4(j);
u_C4=matrixu_C4(j);
I_L5=matrixI_L5(j);
u_C5=matrixu_C5(j);

ku6=1;
while(~(t>Trig(7)))
        matrixu_C6(j)=u_C0;
        matrixI_L6(j) = I_L6;
        matrixu_load(j)=u_load;
        matrix_v(j)=v;
        matrix_x(j)=x;
        matrix_Isource(j)=Isource;
        matrix_Ia(j)=Ia;
        
        matrixI_L1(j)=I_L1;
        matrixu_C1(j)=u_C1;
        matrixI_L2(j)= I_L2;
        matrixu_C2(j)=u_C2;
        matrixI_L3(j)= I_L3;
        matrixu_C3(j)=u_C3;
        matrixI_L3(j)= I_L3;
        matrixu_C3(j)=u_C3;
        matrixI_L4(j)= I_L4;
        matrixu_C4(j)=u_C4;
        matrixI_L5(j)= I_L5;
        matrixu_C5(j)=u_C5;
        
        if(flag6==1)
             u_C6=u_load+(u_C0-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(u_C0-u_load)-I_L6/C)*exp(-af*dt)*sin(w*dt);
             I_L6=I_L6*exp(-af*dt)*cos(w*dt)+1/w*((u_C0-u_load)/L-af*I_L6)*exp(-af*dt)*sin(w*dt);
        if(u_C6<0)
            flag6=2;
            u_C6=0;
            I_L6=matrixI_L6(j) ;
        end
        else if(flag6==2)
                   I_L6=(u_load/R+I_L6)*exp(-dt/tao)-u_load/R;
                   if(I_L6<0)
                       flag6=3;
                       I_L6=0;
                   end
            end
        end
      
        if(flag5==1)
             m5=u_C5;
             u_C5=u_load+(m5-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m5-u_load)-I_L5/C)*exp(-af*dt)*sin(w*dt);
             I_L5=I_L5*exp(-af*dt)*cos(w*dt)+1/w*((m5-u_load)/L-af*I_L5)*exp(-af*dt)*sin(w*dt);
       if(u_C5<0)
            flag5=2;
            u_C5=0;
            I_L5=matrixI_L5(j) ;
       end
       else if(flag5==2)
                   I_L5=(u_load/R+I_L5)*exp(-dt/tao)-u_load/R;
                   if(I_L5<0)
                       flag5=3;
                       I_L5=0;
                   end
            end
       end
       
       if(flag4==1)
            m4=u_C4;
             u_C4=u_load+(m4-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m4-u_load)-I_L4/C)*exp(-af*dt)*sin(w*dt);
             I_L4=I_L4*exp(-af*dt)*cos(w*dt)+1/w*((m4-u_load)/L-af*I_L4)*exp(-af*dt)*sin(w*dt);
       if(u_C4<0)
            flag4=2;
            u_C4=0;
            I_L4=matrixI_L4(j) ;
       end
       else if(flag4==2)
                   I_L4=(u_load/R+I_L4)*exp(-dt/tao)-u_load/R;
                   if(I_L4<0)
                       flag4=3;
                       I_L4=0;
                   end
            end
       end
        
        if(flag3==1)
             m3=u_C3;
             u_C3=u_load+(m3-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m3-u_load)-I_L3/C)*exp(-af*dt)*sin(w*dt);
             I_L3=I_L3*exp(-af*dt)*cos(w*dt)+1/w*((m3-u_load)/L-af*I_L3)*exp(-af*dt)*sin(w*dt);
        if(u_C3<0)
            flag3=2;
            u_C3=0;
            I_L3=matrixI_L3(j) ;
        end
       else if(flag3==2)
                   I_L3=(u_load/R+I_L3)*exp(-dt/tao)-u_load/R;
                   if(I_L3<0)
                       flag3=3;
                       I_L3=0;
                   end
           end
           end
        
        if(flag2==1)
            m2=u_C2;
             u_C2=u_load+(m2-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m2-u_load)-I_L2/C)*exp(-af*dt)*sin(w*dt);
             I_L2=I_L2*exp(-af*dt)*cos(w*dt)+1/w*((m2-u_load)/L-af*I_L2)*exp(-af*dt)*sin(w*dt);
        if(u_C2<0)
            flag2=2;
            u_C2=0;
            I_L2=matrixI_L2(j) ;
        end
        else  if(flag2==2)
                   I_L2=(u_load/R+I_L2)*exp(-dt/tao)-u_load/R;
                   if(I_L2<0)
                       flag2=3;
                       I_L2=0;
                   end
            end
        end
     
        
     if(flag1==1)
        m1=u_C1;  %middle variable
        u_C1=u_load+(m1-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m1-u_load)-I_L1/C)*exp(-af*dt)*sin(w*dt);
        I_L1=I_L1*exp(-af*dt)*cos(w*dt)+1/w*((m1-u_load)/L-af*I_L1)*exp(-af*dt)*sin(w*dt);   
        if(u_C1<1)
            flag1=2;
            u_C1=0;
        end
        else if(flag1==2)
            I_L1=(u_load/R+I_L1)*exp(-dt/tao)-u_load/R;
        if(I_L1<0)
            flag1=3;
            I_L1=0;
        end         
     end
     end
        
        Isource=N(1)*I_L1+N(2)*I_L2+N(3)*I_L3+N(4)*I_L4+N(5)*I_L5+N(6)*I_L6;    
        
        t=t+dt;            
        
        %负载等效电路模型求解
       resistivity1=resistivity10+bert*Isource/height;
       resistivity2=resistivity20+bert*Ia/height;
      %导轨电阻 外轨和内轨 
       Rra=8*Length_rail/3/height*(Permeability*resistivity1/2/pi/t)^0.5;
       Rrb=8*x/3/height*(Permeability*resistivity1/2/pi/t)^0.5+8*ZDWZ/3/height*(Permeability*resistivity1/2/pi/t)^0.5;
       Rr=Rra+Rrb;
       
       Rc=kc/height*(Permeability*(resistivity1+resistivity2)/4/pi/t)^0.5;
       RVSEC=Rvc*v^1.5;
       RVSEC1=(1-kVSEC)*RVSEC;
       RVSEC2=kVSEC*RVSEC;
       Ra=b/height*(Permeability*resistivity2/t/pi/2)^0.5;
       uEMF=L_gradient*v*Isource;
       
       
       if(j>1)
       di=Isource-matrix_Isource(j-1);
             else
        di=(Isource-matrix_Isource(j))*2;    
       end
       
        if(ku6>N_Lr0)
           Lr0=Lr0;
       else
           Lr0=Lr0+Lr_initial/N_Lr0;
       end
       Lr=Lr0+L_gradient*x;   
       u1_load=(Rr+Rc+RVSEC1+RVSEC2+Ra+L_gradient*v)*Isource;%+di/dt*Lr;
       u2_load=R_cabel*Isource+L_cabel*di/dt;
       u3_load=di/dt/2*Lr;   
       u_load=u1_load+u2_load+u3_load; 
       
       if(u_load>1500)
           u_load=800;%-u_load;  % 遏制振荡   从它变为负数开始 强制校正为正 
       end
       if(u_load<0)
           u_load=60;%-u_load;  % 遏制振荡   从它变为负数开始 强制校正为正 
       end
       
       F=L_gradient*Isource^2/2;
      if(v>0)
           Ff=uf2*FN;
       else
           Ff=Ff_staic;
       end
       if(F>Ff)
           a=(F-Ff)/m;
       else
           a=0;
       end
       
       v=v+a*dt;
        x=x+matrix_v(j)*dt+1/2*a*(dt^2+2*(t-dt)*dt);
       if(x>Length)
           break
       end
        u_C0=u_C6;
              
        j=j+1;
        ku6=ku6+1;
end

%最后一组 或最后一个放电阶段 
%第七组放电计算
 t=t-dt;
 j=j-1;
 flag7=1;

  u_C0=U0;      %电容初始电压
  
  Lr0=0;    %!!! 逐渐加电感
  I_L7=0;          %电感电流初始值为0
  u_load=matrixu_load(j);        %负载两端电压 
  %负载初始状态 
 v=matrix_v(j);
 x=matrix_x(j);
 Isource=matrix_Isource(j);
 Ia= matrix_Ia(j);
 
 I_L1=matrixI_L1(j) ; 
 u_C1=matrixu_C1(j);
 I_L2=matrixI_L2(j);
u_C2=matrixu_C2(j);
 I_L3=matrixI_L3(j);
u_C3=matrixu_C3(j);
I_L4=matrixI_L4(j);
u_C4=matrixu_C4(j);
I_L5=matrixI_L5(j);
u_C5=matrixu_C5(j);
I_L6=matrixI_L6(j);
u_C6=matrixu_C6(j);

ku7=1;
while(x<Length)
        matrixu_C7(j)=u_C0;
        matrixI_L7(j) = I_L7;
        matrixu_load(j)=u_load;
        matrix_v(j)=v;
        matrix_x(j)=x;
        matrix_Isource(j)=Isource;
        matrix_Ia(j)=Ia;
        
        matrixI_L1(j)=I_L1;
        matrixu_C1(j)=u_C1;
        matrixI_L2(j)= I_L2;
        matrixu_C2(j)=u_C2;
        matrixI_L3(j)= I_L3;
        matrixu_C3(j)=u_C3;
        matrixI_L3(j)= I_L3;
        matrixu_C3(j)=u_C3;
        matrixI_L4(j)= I_L4;
        matrixu_C4(j)=u_C4;
        matrixI_L5(j)= I_L5;
        matrixu_C5(j)=u_C5;
        matrixI_L6(j)= I_L6;
        matrixu_C6(j)=u_C6;
        
       if(flag7==1)
             u_C7=u_load+(u_C0-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(u_C0-u_load)-I_L7/C)*exp(-af*dt)*sin(w*dt);
             I_L7=I_L7*exp(-af*dt)*cos(w*dt)+1/w*((u_C0-u_load)/L-af*I_L7)*exp(-af*dt)*sin(w*dt);
        if(u_C7<0)
            flag7=2;
            u_C7=0;
            I_L7=matrixI_L7(j) ;
        end
        else if(flag7==2)
                   I_L7=(u_load/R+I_L7)*exp(-dt/tao)-u_load/R;
                   if(I_L7<0)
                       flag7=3;
                       I_L7=0;
                   end
            end
       end
       
        if(flag6==1)
            m6=u_C6;
             u_C6=u_load+(m6-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m6-u_load)-I_L6/C)*exp(-af*dt)*sin(w*dt);
             I_L6=I_L6*exp(-af*dt)*cos(w*dt)+1/w*((m6-u_load)/L-af*I_L6)*exp(-af*dt)*sin(w*dt);
        if(u_C6<0)
            flag6=2;
            u_C6=0;
            I_L6=matrixI_L6(j) ;
        end
        else if(flag6==2)
                   I_L6=(u_load/R+I_L6)*exp(-dt/tao)-u_load/R;
                   if(I_L6<0)
                       flag6=3;
                       I_L6=0;
                   end
            end
        end
       
        if(flag5==1)
            m5=u_C5;
             u_C5=u_load+(m5-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m5-u_load)-I_L5/C)*exp(-af*dt)*sin(w*dt);
             I_L5=I_L5*exp(-af*dt)*cos(w*dt)+1/w*((m5-u_load)/L-af*I_L5)*exp(-af*dt)*sin(w*dt);
        if(u_C5<0)
            flag5=2;
            u_C5=0;
            I_L5=matrixI_L5(j) ;
        end
        else if(flag5==2)
                   I_L5=(u_load/R+I_L5)*exp(-dt/tao)-u_load/R;
                   if(I_L5<0)
                       flag5=3;
                       I_L5=0;
                   end
            end
       end
        
        if(flag4==1)
            m4=u_C4;
             u_C4=u_load+(m4-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m4-u_load)-I_L4/C)*exp(-af*dt)*sin(w*dt);
             I_L4=I_L4*exp(-af*dt)*cos(w*dt)+1/w*((m4-u_load)/L-af*I_L4)*exp(-af*dt)*sin(w*dt);
        if(u_C4<0)
            flag4=2;
            u_C4=0;
            I_L4=matrixI_L4(j) ;
        end
        else if(flag4==2)
                   I_L4=(u_load/R+I_L4)*exp(-dt/tao)-u_load/R;
                   if(I_L4<0)
                       flag4=3;
                       I_L4=0;
                   end
            end
        end
        
           if(flag3==1)
             m3=u_C3;
             u_C3=u_load+(m3-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m3-u_load)-I_L3/C)*exp(-af*dt)*sin(w*dt);
             I_L3=I_L3*exp(-af*dt)*cos(w*dt)+1/w*((m3-u_load)/L-af*I_L3)*exp(-af*dt)*sin(w*dt);
        if(u_C3<0)
            flag3=2;
            u_C3=0;
            I_L3=matrixI_L3(j) ;
        end
       else if(flag3==2)
                   I_L3=(u_load/R+I_L3)*exp(-dt/tao)-u_load/R;
                   if(I_L3<0)
                       flag3=3;
                       I_L3=0;
                   end
           end
           end
        
        if(flag2==1)
            m2=u_C2;
             u_C2=u_load+(m2-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m2-u_load)-I_L2/C)*exp(-af*dt)*sin(w*dt);
             I_L2=I_L2*exp(-af*dt)*cos(w*dt)+1/w*((m2-u_load)/L-af*I_L2)*exp(-af*dt)*sin(w*dt);
        if(u_C2<0)
            flag2=2;
            u_C2=0;
            I_L2=matrixI_L2(j) ;
        end
        else  if(flag2==2)
                   I_L2=(u_load/R+I_L2)*exp(-dt/tao)-u_load/R;
                   if(I_L2<0)
                       flag2=3;
                       I_L2=0;
                   end
            end
        end
     
        
     if(flag1==1)
        m1=u_C1;  %middle variable
        u_C1=u_load+(m1-u_load)*exp(-af*dt)*cos(w*dt)+1/w*(af*(m1-u_load)-I_L1/C)*exp(-af*dt)*sin(w*dt);
        I_L1=I_L1*exp(-af*dt)*cos(w*dt)+1/w*((m1-u_load)/L-af*I_L1)*exp(-af*dt)*sin(w*dt);   
        if(u_C1<1)
            flag1=2;
            u_C1=0;
        end
        else if(flag1==2)
            I_L1=(u_load/R+I_L1)*exp(-dt/tao)-u_load/R;
        if(I_L1<0)
            flag1=3;
            I_L1=0;
        end         
     end
     end
        
        Isource=N(1)*I_L1+N(2)*I_L2+N(3)*I_L3+N(4)*I_L4+N(6)*I_L5+N(6)*I_L6+N(7)*I_L7;    %n个模块  第一次计算时的公式！！！ 历史电流？
   
        t=t+dt;            
        
        %负载等效电路模型求解
       resistivity1=resistivity10+bert*Isource/height;
       resistivity2=resistivity20+bert*Ia/height;
      %导轨电阻 外轨和内轨 
       Rra=8*Length_rail/3/height*(Permeability*resistivity1/2/pi/t)^0.5;
       Rrb=8*x/3/height*(Permeability*resistivity1/2/pi/t)^0.5+8*ZDWZ/3/height*(Permeability*resistivity1/2/pi/t)^0.5;
       Rr=Rra+Rrb;
       
       Rc=kc/height*(Permeability*(resistivity1+resistivity2)/4/pi/t)^0.5;
       RVSEC=Rvc*v^1.5;
       RVSEC1=(1-kVSEC)*RVSEC;
       RVSEC2=kVSEC*RVSEC;
       Ra=b/height*(Permeability*resistivity2/t/pi/2)^0.5;
       uEMF=L_gradient*v*Isource;
        
       if(j>1)
       di=Isource-matrix_Isource(j-1);
             else
       di=(Isource-matrix_Isource(j))*2;    
       end
       
        if(ku7>N_Lr0)
           Lr0=Lr0;
       else
           Lr0=Lr0+Lr_initial/N_Lr0;
       end
       Lr=Lr0+L_gradient*x;   
       u1_load=(Rr+Rc+RVSEC1+RVSEC2+Ra+L_gradient*v)*Isource; %+di/dt*Lr;
       u2_load=R_cabel*Isource+L_cabel*di/dt;
       u3_load=di/dt/2*Lr;   
       u_load=u1_load+u2_load+u3_load; 
            
       if(u_load>1500)
           u_load=1000;%-u_load;  % 遏制振荡   从它变为负数开始 强制校正为正 
       end
       if(u_load<0)
           u_load=60;%-u_load;  % 遏制振荡   从它变为负数开始 强制校正为正 
       end
       
      F=L_gradient*Isource^2/2;
      if(v>0)
           Ff=uf2*FN;
       else
           Ff=Ff_staic;
       end
       if(F>Ff)
           a=(F-Ff)/m;
       else
           a=0;
       end
       
       v=v+a*dt;
       x=x+matrix_v(j)*dt+1/2*a*(dt^2+2*(t-dt)*dt);
    
        u_C0=u_C7;
              
        j=j+1;
        ku7=ku7+1;
end
t=t-dt;
j=j-1;
tim=(0:(j-1))*dt*1000;  %时间向量  单位ms

matrix_Isource=matrix_Isource(1:j)*10^-3;
D_sum_I(ik)=sum(matrix_Isource)*dt*1000;
matrixI_L1=matrixI_L1(1:j)*10^-3;
matrixI_L2=matrixI_L2(1:j)*10^-3;
matrixI_L3=matrixI_L3(1:j)*10^-3;
matrixI_L4=matrixI_L4(1:j)*10^-3;
matrixI_L5=matrixI_L5(1:j)*10^-3;
matrixI_L6=matrixI_L6(1:j)*10^-3;
matrixI_L7=matrixI_L7(1:j)*10^-3;

SumI_L1=matrixI_L1*N(1);
SumI_L2=matrixI_L2*N(2);
SumI_L3=matrixI_L3*N(3);
SumI_L4=matrixI_L4*N(4);
SumI_L5=matrixI_L5*N(5);
SumI_L6=matrixI_L6*N(6);
SumI_L7=matrixI_L7*N(7);

matrixu_C1=matrixu_C1(1:j)*10^-3;
matrixu_C2=matrixu_C2(1:j)*10^-3;
matrixu_C3=matrixu_C3(1:j)*10^-3;
matrixu_C4=matrixu_C4(1:j)*10^-3;
matrixu_C5=matrixu_C5(1:j)*10^-3;
matrixu_C6=matrixu_C6(1:j)*10^-3;
matrixu_C7=matrixu_C7(1:j)*10^-3;
matrixu_load=matrixu_load(1:j);

%画图及数据输出
% x
% t
Loop_v(ik,1)=v; 
Loop_Length(ik,1)=Length;
Loop_m(ik,1)=m;
Loop_I = D_sum_I';



%{
figure(1*ik)
plot(tim,matrixI_L1,tim,matrixI_L2,tim,matrixI_L3,tim,matrixI_L4,tim,matrixI_L5,tim,matrixI_L6,tim,matrixI_L7);
grid on;
xlabel('t/ms');
ylabel('I/kA');
title('1123-1');
legend('matrixI_L1',' matrixI_L2','matrixI_L3','matrixI_L4',' matrixI_L5','matrixI_L6','matrixI_L7');

figure(20*ik);
plot(tim,SumI_L1,'--',tim,SumI_L2,'--',tim,SumI_L3,'--',tim,SumI_L4,'--',tim,SumI_L5,'--',tim,SumI_L6,'--',tim,SumI_L7,'--',tim,matrix_Isource);
grid on;
xlabel('t/ms');
ylabel('I/kA');
legend('SumI_L1',' SumI_L2','SumI_L3','SumI_L4',' SumI_L5','SumI_L6','SumI_L7','matrix_Isource');
%}

end