%% 案例7：RBF网络的回归-非线性函数回归的实现 


%% 清空环境变量
clc
clear all
close all
%% 产生训练样本（训练输入，训练输出）
% ld为样本例数
V_data=xlsread('C:\Users\hn\AnacondaProjects\python_coding\rbf_training_data940ILmv.xls');
ld=length(V_data(:,1)); 

S_num=ld;
test_num = 20;
train_num = 400 - test_num;
ik_int = 0;


%% 数据标准化

for ik = 0.1:0.1:1

    train_num_ik = train_num*ik;
    data_scal = [];

    for i = 1:4
        mean_data = mean(V_data(1:train_num_ik,i), 1);
        std_data = std(V_data(1:train_num_ik,i), 0, 1);
        data_scal(:,i) = ( V_data(1:train_num_ik,i) - mean_data ) / std_data;
    
    end

x_train = data_scal(test_num+1:train_num_ik , 1:3)';
y_train = data_scal(test_num+1:train_num_ik ,4)';


%% 建立RBF神经网络 
% 采用approximate RBF神经网络。spread为默认值

%x_train = x_scal(:,test_num+1:train_num*ik);
%v_train = v_scal(:,test_num+1:train_num*ik);

net=newrb(x_train ,y_train,0.01,0.93,160,4);

%% 建立测试样本

x_test = data_scal( 1:test_num , 1:3)';
y_test = data_scal( 1:test_num , 4)';

%% 使用建立的RBF网络进行模拟，得出网络输出

y_test_pred = sim(net,x_test);
y_train_pred = sim(net,x_train);

%% 标准化还原函数

ytest_pred_orig = y_test_pred*std_data+ mean_data;
ytest_orig = y_test*std_data+ mean_data;

ytrain_pred_orig = y_train_pred*std_data+ mean_data;
ytrain_orig = y_train*std_data+ mean_data;


ik_int = ik_int+1;

train_num_list(1,ik_int) = train_num_ik;
train_diff = ytrain_orig - ytrain_pred_orig;
train_diffmean(1,ik_int) = mean(train_diff);
train_diffstd(1,ik_int) = std(train_diff);

test_diff = ytest_pred_orig - ytest_orig
test_diff_abs = abs(ytest_orig - ytest_pred_orig);
test_diffmean(1,ik_int) = mean(test_diff_abs);
test_diffstd(1,ik_int) = std(test_diff);

end



















