%% 案例7：RBF网络的回归-非线性函数回归的实现 


%% 清空环境变量
clc
clear all
close all
%% 产生训练样本（训练输入，训练输出）
% ld为样本例数
V_data=xlsread('C:\Users\hn\AnacondaProjects\w-RBFforEML\rbf_training_data940ILmv.xls');
ld=length(V_data(:,1)); 

S_num=ld;
test_num = 19
train_num = S_num - test_num

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
v_scal = data_scal2(:,4)'


%% 建立RBF神经网络 
% 采用approximate RBF神经网络。spread为默认值

x_train = x_scal(:,test_num+1:300);
v_train = v_scal(:,test_num+1:300);

net=newrb(x_train ,v_train);

%% 建立测试样本

x_test = x_scal(: , 1:test_num);
y_test = v_scal(: , 1:test_num);

%% 使用建立的RBF网络进行模拟，得出网络输出

y_test_pred = sim(net,x_test);

%% 标准化还原函数

ytest_pred_orig = y_test_pred*std_data+ mean_data
ytest_orig = y_test*std_data+ mean_data

diff = ytest_orig - ytest_pred_orig
diff_mean = mean(diff)
diff_std = std(diff)


