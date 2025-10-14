% % 绘制初始温盐场
clc
clear

%% 读取网格信息
grd_lb=textread('TEST.grd','','headerlines',17); %网格文件
NUM=44561; %网格节点数!!
KB=11; %垂向分层数+1!!
fsm=grd_lb(1:NUM,6); %干湿判断
lon=grd_lb(1:NUM,4); %经度
lat=grd_lb(1:NUM,5); %纬度

%% 读取温盐初始场数据
its=readtable('TEST_its.dat','NumHeaderLines',1);
its=its{2:2:end,:};
its1=its(:,2); %表层
its2=its(:,3); %第2层
its3=its(:,4); %第3层
its4=its(:,5); %第4层
its5=its(:,6); %第5层
its6=its(:,7); %第6层
its7=its(:,8); %第7层
its8=its(:,9); %第8层
its9=its(:,10); %第9层
its10=its(:,11); %底层
its11=its(:,12); %b

%% 剔除干点
lon(fsm<0)=[];
lat(fsm<0)=[];
its1(fsm<0)=[];
its2(fsm<0)=[];
its3(fsm<0)=[];
its4(fsm<0)=[];
its5(fsm<0)=[];
its6(fsm<0)=[];
its7(fsm<0)=[];
its8(fsm<0)=[];
its9(fsm<0)=[];
its10(fsm<0)=[];
its11(fsm<0)=[];

%% 画图
cla
scatter(lon,lat,3,its1,'filled');hold on
% scatter(lon,lat,5,its11,'filled');
h=colorbar;
h.Label.String = 'Sailinity';
axis equal
box on

%设置为0的区域
x=[117,120.1,120.1,117,117];
y=[32.5,32.5,31,31,32.5];
plot(x,y,'r')
legend('','手动设置为0的区域')
fclose('all')
