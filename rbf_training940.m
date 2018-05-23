%% ����7��RBF����Ļع�-�����Ժ����ع��ʵ�� 


%% ��ջ�������
clc
clear all
close all
%% ����ѵ��������ѵ�����룬ѵ�������
% ldΪ��������
V_data=xlsread('C:\Users\hn\AnacondaProjects\w-RBFforEML\rbf_training_data940ILmv.xls');
ld=length(V_data(:,1)); 
% ����2*ld�ľ��� 
% x=rand(2,ld); 

% ��xת����[-1.5 1.5]֮��
% x=(x-0.5)*1.5*2; 

% x�ĵ�һ��Ϊx1���ڶ���Ϊx2.
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
%% ���ݱ�׼��

data_scal2 = []

for i = 1:4
    mean_data = mean(V_data(:,i), 1)
    std_data = std(V_data(:,i), 0, 1)
    data_scal2(:,i) = ( V_data(:,i) - mean_data ) / std_data
    
end

x_scal = data_scal2(:,1:3)';
v_scal = data_scal2(:,4)';

%% ����RBF������ 
% ����approximate RBF�����硣spreadΪĬ��ֵ

x_train = x_scal(:,test_num+1:940);
v_train = v_scal(:,test_num+1:940);

net=newrb(x_train ,v_train);

%% ������������

x_test = x_scal(: , 1:test_num)
y_test = v_scal(: , 1:test_num)

%% ʹ�ý�����RBF�������ģ�⣬�ó��������

y_test_pred = sim(net,x_test);

%% ��׼����ԭ����

ytest_pred_orig = y_test_pred*std_data+ mean_data
ytest_orig = y_test*std_data+ mean_data

diff = ytest_orig - ytest_pred_orig
    
%% ʹ��ͼ�񣬻���3άͼ

% �����ĺ���ͼ��
t1=21;
t2=40;
yanV_data=V_data(t1:t2,1)';
yanM_data=V_data(t1:t2,2)';
yan_data=[yanV_data;yanM_data];
yan_y_data=sim(net,yan_data);
figure(3)
plot(t1-20:t2-20,10^(12)*(yan_y_data(1,:)-V_data(t1:t2,3)'),'r-')
xlabel('�������')
ylabel('���������ʱ����ԭʼ����ʱ��Ĳ�')
title('����ʱ�����ֲ�')
print -djpeg -r600 'C:\Users\hn\AnacondaProjects\w-RBFforEML\RBFͼƬ\1.jpg';
tEan=mean(abs(10^(12)*(yan_y_data(1,:)-V_data(t1:t2,3)')));
tsd=std2(10^(12)*(yan_y_data(1,:)-V_data(t1:t2,3)'));
figure(4)
plot(t1-20:t2-20,50*(round(yan_y_data(2,:))-V_data(t1:t2,4)'),'r-')
xlabel('�������')
ylabel('�������������ԭʼ�������Ĳ�')
title('���������ֲ�')
print -djpeg -r600 'C:\Users\hn\AnacondaProjects\w-RBFforEML\RBFͼƬ\2.jpg';
IEan=mean(abs(50*(round(yan_y_data(2,:))-V_data(t1:t2,4)')));
Isd=std2(50*(round(yan_y_data(2,:))-V_data(t1:t2,4)'));

% ����ó��ĺ���ͼ��
% v=reshape(ty,row);
% subplot(1,3,2)
% mesh(i,j,v);
% zlim([0,2])
% title('RBF��������')


% ���ͼ��
% subplot(1,3,3)
% mesh(x1,x2,F-v);
% zlim([0,60])
% title('���ͼ��')
% 
% set(gcf,'position',[300 ,250,900,400])

%%
% 
% <html>
% <table align="center" >	<tr>		<td align="center"><font size="2">��Ȩ���У�</font><a
% href="http://www.ilovematlab.cn/">Matlab������̳</a>&nbsp;&nbsp; <script
% src="http://s3.cnzz.com/stat.php?id=971931&web_id=971931&show=pic" language="JavaScript" ></script>&nbsp;</td>	</tr></table>
% </html>
% 

