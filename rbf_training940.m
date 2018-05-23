%% 案例7：RBF网络的回归-非线性函数回归的实现 


%% 清空环境变量
clc
clear all
close all
%% 产生训练样本（训练输入，训练输出）
% ld为样本例数
V_data=xlsread('C:\Users\hn\AnacondaProjects\w-RBFforEML\rbf_training_data940ILmv.xls');
ld=length(V_data(:,1)); 
% 产生2*ld的矩阵 
% x=rand(2,ld); 

% 将x转换到[-1.5 1.5]之间
% x=(x-0.5)*1.5*2; 

% x的第一列为x1，第二列为x2.
% x1=x(1,:);
% x2=x(2,:);
S_num=ld;
train_num = S_num - 19
test_num = 19

x=V_data(:,1:3)';
x1_I=V_data(:,1)';
x2_L=V_data(:,2)';
x3_m=V_data(:,3)';

v=V_data(: ,4)';
%% 数据标准化

data_scal2 = []

for i = 1:4
    mean_data = mean(V_data(:,i), 1)
    std_data = std(V_data(:,i), 0, 1)
    data_scal2(:,i) = ( V_data(:,i) - mean_data ) / std_data
    
end

x_scal = data_scal2(:,1:3)';
v_scal = data_scal2(:,4)';

%% 建立RBF神经网络 
% 采用approximate RBF神经网络。spread为默认值

x_train = x_scal(:,test_num+1:940);
v_train = v_scal(:,test_num+1:940);

net=newrb(x_train ,v_train);

%% 建立测试样本

x_test = x_scal(: , 1:test_num)
y_test = v_scal(: , 1:test_num)

%% 使用建立的RBF网络进行模拟，得出网络输出

y_test_pred = sim(net,x_test);

%% 标准化还原函数

ytest_pred_orig = y_test_pred*std_data+ mean_data
ytest_orig = y_test*std_data+ mean_data

diff = ytest_orig - ytest_pred_orig
    
%% 使用图像，画出3维图

% 真正的函数图像
t1=21;
t2=40;
yanV_data=V_data(t1:t2,1)';
yanM_data=V_data(t1:t2,2)';
yan_data=[yanV_data;yanM_data];
yan_y_data=sim(net,yan_data);
figure(3)
plot(t1-20:t2-20,10^(12)*(yan_y_data(1,:)-V_data(t1:t2,3)'),'r-')
xlabel('样本编号')
ylabel('计算的上升时间与原始上升时间的差')
title('上升时间误差分布')
print -djpeg -r600 'C:\Users\hn\AnacondaProjects\w-RBFforEML\RBF图片\1.jpg';
tEan=mean(abs(10^(12)*(yan_y_data(1,:)-V_data(t1:t2,3)')));
tsd=std2(10^(12)*(yan_y_data(1,:)-V_data(t1:t2,3)'));
figure(4)
plot(t1-20:t2-20,50*(round(yan_y_data(2,:))-V_data(t1:t2,4)'),'r-')
xlabel('样本编号')
ylabel('计算的最大电流与原始最大电流的差')
title('最大电流误差分布')
print -djpeg -r600 'C:\Users\hn\AnacondaProjects\w-RBFforEML\RBF图片\2.jpg';
IEan=mean(abs(50*(round(yan_y_data(2,:))-V_data(t1:t2,4)')));
Isd=std2(50*(round(yan_y_data(2,:))-V_data(t1:t2,4)'));

% 网络得出的函数图像
% v=reshape(ty,row);
% subplot(1,3,2)
% mesh(i,j,v);
% zlim([0,2])
% title('RBF神经网络结果')


% 误差图像
% subplot(1,3,3)
% mesh(x1,x2,F-v);
% zlim([0,60])
% title('误差图像')
% 
% set(gcf,'position',[300 ,250,900,400])

%%
% 
% <html>
% <table align="center" >	<tr>		<td align="center"><font size="2">版权所有：</font><a
% href="http://www.ilovematlab.cn/">Matlab中文论坛</a>&nbsp;&nbsp; <script
% src="http://s3.cnzz.com/stat.php?id=971931&web_id=971931&show=pic" language="JavaScript" ></script>&nbsp;</td>	</tr></table>
% </html>
% 

