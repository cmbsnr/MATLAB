%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      网络结构分析                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = net(state, result,pro)
%分析网络结构并绘图
if isnumeric(state) && state == 0
    result = net_init(result,pro);
    return;
elseif isnumeric(state) && state == -1
    result = net_term(result,pro);
    return;
elseif ~isstruct(state)
    disp('Error input state for processing')
    pauss;
end

if ~isfield(result,'const') %不变量
    result.const.TotFunc0 = sum((state.gr == 1 | state.gr == 2),'all');%总官能度，不变
    result.const.gr0 = sum(state.gr > 0, 2);%官能团数
    result.const.MaxDeg0 = max(result.const.gr0,'all');%最大官能团数，不变
    result.const.isNode0 = result.const.gr0 > 2;%筛选多官能团点，不变
    result.const.TotNode0 = sum(result.const.isNode0);%多官能团点总数，不变
end

func = getfunc;
G0 = func.ord2graph(state,pro);%将ord转变为图
[Treeindex,STree] = conncomp(G0); %每个点所属树，每个树的尺寸
Treeindex = Treeindex';
state.G0 = G0;
state.Treeindex = Treeindex;%传递给morph

NTree = length(STree);%总分子数
Nmaxtree = max(STree);%最大分子
Nave_mean = mean(STree); Nave_std = std(STree);%数均分子量及其偏差
Wave_mean = sum(STree.^2)/state.TotN; Wave_std = std(STree.^2)*NTree/state.TotN;%重均分子量及其偏差
PolyDisp = Wave_mean/Nave_mean;%多分散性
ReactFunc = sum((state.gr == 1 | state.gr == 2) & state.ca ~= 0,'all');%反应的官能团数
ReactDeg = ReactFunc/result.const.TotFunc0;%反应程度
isReactNode = sum (state.ca > 0, 2) > 2; %是否为外连大于2的节点
ReactNodeDeg = sum (isReactNode)/result.const.TotNode0; %多官能团点的反应程度（连接度大于2的节点占总多官能团点的比例）
NodeFuncDeg =  mean(sum(state.ca(result.const.isNode0,:) > 0, 2)); % 官能团>2点的平均支化度（平均外连数）

G = ord2simplegraph(state,pro);%将ord转变为简化图,链段变带权重直线
ReactNodeList = findnode(G,string(find(isReactNode)));%图G中支化节点的序号
%BranchRatioMean = Net(1)*2/sum(degree(G,ReactNodeList)); %节点的平均网络率，联通至网络的链占臂的比例，一个net链连两个节点
NetEdgesList = find(G.Edges.Style=='Net');%网络链的序号列表
BranchRatio = zeros(length(ReactNodeList),1);
j = 0;
for i = ReactNodeList'
    j = j + 1;
    e = outedges(G,i);%网络节点的出链序号
    %BranchRatio(j) = sum(G.Edges(e,:).Style == 'Net')/length(e);%每个支化节点的网络率，网络边在所有边中所占的比例，每次判断较慢
    BranchRatio(j) = sum(ismember(e,NetEdgesList))/length(e);%每个支化节点的网络率,网络边在支化点外连所有边中所占的比例
end
BranchRatio_mean = mean(BranchRatio);%网络支化率均值
BranchRatio_std = std(BranchRatio);%网络支化率方差

fxx = @(x) [numel(x),mean(x),std(x)];%统计网络链
[LGroup,LGroupName] = findgroups(G.Edges.Style);
info = splitapply(fxx,G.Edges.Weight,LGroup);

%画简化掉环，重复边，悬挂链和中间节点%
SG = simplify(G,'min');%去掉自环，按权重最小保留重复边
SG = rmnode(SG,find(degree(SG)==1));%去掉悬挂链
SG = chain2edge(SG);%将链变成边
NetStrong = fxx(SG.Edges.Weight);%强网络链

[~,Edgecircle] = cyclebasis(G);%搜索网络环
fx1 = @(x) sum(G.Edges(x,:).Weight);%网络环权重函数
circleweight = cellfun(fx1,Edgecircle);
NetCircle = fxx(circleweight);%网络环的统计


m = max([NetStrong(1);NetCircle(1);info(:,1)]);
if  m > result.maxrow_hist
    result.maxrow_hist = m;
    result.data_hist(m,width(result.data_hist)) = 0; %扩充已有hist数组
end

[i,j]=ismember('Dangle',LGroupName);
if ~i
    Dangle = [0,0,0];
    Dangleweight = zeros(result.maxrow_hist,1);
else
    Dangle = info(j,:);
    Dangleweight = G.Edges.Weight(LGroup==j);
end
[i,j]=ismember('Free',LGroupName);
if ~i
    Free = [0,0,0];
    Freeweight = zeros(result.maxrow_hist,1);
else
    Free = info(j,:);
    Freeweight = G.Edges.Weight(LGroup==j);
end
[i,j]=ismember('Loop',LGroupName);
if ~i
    Loop = [0,0,0];
    Loopweight = zeros(result.maxrow_hist,1);
else
    Loop = info(j,:);
    Loopweight = G.Edges.Weight(LGroup==j);
end
[i,j]=ismember('Net',LGroupName);
if ~i
    Net = [0,0,0];
    Netweight = zeros(result.maxrow_hist,1);
else
    Net = info(j,:);
    Netweight = G.Edges.Weight(LGroup==j);
end

result.data_series(pro.Index,:) = [{NTree,Nmaxtree,Nave_mean,Nave_std,Wave_mean,Wave_std,PolyDisp,...%汇总series数据
    ReactDeg,ReactNodeDeg,NodeFuncDeg,BranchRatio_mean,BranchRatio_std},...
    num2cell([Dangle,Free,Loop,Net,NetStrong,NetCircle])];

fmaxrow = @(x)(vertcat(x,zeros(result.maxrow_hist-height(x),1)));
data_hist = cell2mat(cellfun(fmaxrow,{Dangleweight,Freeweight,Loopweight,Netweight,SG.Edges.Weight,circleweight},'UniformOutput',false));
result.data_hist(:,pro.Index:pro.Cyclenum:6*pro.Cyclenum) = data_hist;%汇总
if pro.Noplot%不画图
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%画图%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure(result.fig(1));%%%%%%%%%%%%%%%%%%%%%%%%%统计图
index = 1:pro.Index;
Cycle = pro.Resultlist(index,1);
a1=nexttile(1);
%反应程度，节点反应程度，节点支化率，节点官能度
yyaxis(a1,'left');
hold on
plot(a1,Cycle,result.data_series.ReactDeg(index),'-o');%ReactDeg
plot(a1,Cycle,result.data_series.ReactNodeDeg(index),'-d');%ReactNodeDeg
errorbar(a1,Cycle,result.data_series.NetBranchRatio_ave(index),result.data_series.NetBranchRatio_err(index),'-s');%BranchRatio_mean, BranchRatio_std
ylabel("Degree")
yyaxis(a1,'right');
plot(a1,Cycle,result.data_series.NodeFunctionDegree(index),'-h');%NodeFuncDeg
ylabel("Functionality")
legend('Degree of reaction','Degree of Node reaction','Degree of Net Branching','Functionality of Nodes');
xlim auto
hold off

a3 = nexttile(3);
%数均分子量，重均分子量，分子数
yyaxis(a3,'left');
hold on
%errorbar(a3,Cycle,result.data_series.Nave(index),result.data_series.Nave_err(index),'-o');%Nave_mean,Nave_std
%errorbar(a3,Cycle,result.data_series.Wave(index),result.data_series.Wave_err(index),'-d');%Wave_mean,Wave_std
plot(a3,Cycle,result.data_series.Nave(index),'-d');
plot(a3,Cycle,result.data_series.Wave(index),'-d');
ylabel("Molecular Weight")
yyaxis(a3,'right');
plot(a3,Cycle,result.data_series.NTree(index),'-s');%NTree
ylabel("Molecular Number")
legend('Number Average of Molecular Weight','Weight Average of Molecular Weight','Number of Trees');
xlim auto
hold off

a5=nexttile(5);
%悬挂链、自由链、环形链、网络链、强网络链、网络环长度
hold on
errorbar(a5,Cycle,result.data_series.Dangle_ave(index),result.data_series.Dangle_err(index),'-o');%Dangle(2),Dangle(3)
errorbar(a5,Cycle,result.data_series.Free_ave(index),result.data_series.Free_err(index),'-d');%Free(2),Free(3)
errorbar(a5,Cycle,result.data_series.Loop_ave(index),result.data_series.Loop_err(index),'-s');%Loop(2),Loop(3)
errorbar(a5,Cycle,result.data_series.Net_ave(index),result.data_series.Net_err(index),'-h');%Net(2),Net(3)
errorbar(a5,Cycle,result.data_series.NetStrong_ave(index),result.data_series.NetStrong_err(index),'-*');%NetStrong(2),NetStrong(3)
errorbar(a5,Cycle,result.data_series.NetCircle_ave(index),result.data_series.NetCircle_err(index),'-+');%NetCircle(2),NetCircle(3)
legend('Dangle','Free','Loop','Net','NetStrong','NetCircle');
ylabel("Length")
xlim auto
hold off

a7 = nexttile(7);
%悬挂链、自由链、环形链、网络链、强网络链、网络环数量
hold on
plot(a7,Cycle,result.data_series.Dangle_num(index),'-o','color','#0072BD');%Dangle(1)
plot(a7,Cycle,result.data_series.Free_num(index),'-d','color','#D95319');%Free(1)
plot(a7,Cycle,result.data_series.Loop_num(index),'-s','color','#EDB120');%Loop(1)
plot(a7,Cycle,result.data_series.Net_num(index),'-h','color','#7E2F8E');%Net(1)
plot(a7,Cycle,result.data_series.NetStrong_num(index),'-*','color','#77AC30');%NetStrong(1)
plot(a7,Cycle,result.data_series.NetCircle_num(index),'-+','color','#4DBEEE');%NetCircle(1)
legend('Dangle','Free','Loop','Net','NetStrong','NetCircle');
xlabel("Cycle");
ylabel("Number")
xlim auto
hold off
linkaxes([a1,a3,a5,a7],'x');

a2 = nexttile(2);
%悬挂链分布
histogram(a2,Dangleweight,'Normalization','probability','Facecolor','#0072BD','FaceAlpha',0.3);
legend('Dangle');
xlim auto
a4 =nexttile(4);
%自由链分布
histogram(a4,Freeweight,'Normalization','probability','Facecolor','#D95319','FaceAlpha',0.3);   
legend('Free');
xlim auto
a6 =nexttile(6);
%环形链分布
histogram(a6,Loopweight,'Normalization','probability','Facecolor','#EDB120','FaceAlpha',0.3);    
legend('Loop');
xlim auto
a8 =nexttile(8);
%网络链分布
hold on
histogram(a8,Netweight,'Normalization','probability','Facecolor','#7E2F8E','FaceAlpha',0.3);
histogram(a8,SG.Edges.Weight,'Normalization','probability','Facecolor','#77AC30','FaceAlpha',0.3);
histogram(a8,circleweight,'Normalization','probability','Facecolor','#4DBEEE','FaceAlpha',0.3);
legend('Net','NetStrong','NetCircle');
xlim auto
hold off
linkaxes([a2,a4,a6,a8],'x');
drawnow limitrate
writeVideo(result.video{1}, getframe(f));
linkaxes([a1,a3,a5,a7],'off');
linkaxes([a2,a4,a6,a8],'off');

f = figure(result.fig(end));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%形貌图
nexttile(1);
cla
hold on
%画基本图，调用morph
result = morph(state,result,pro);
%画网络键
if numel(G.Edges) > 0 %非空图
    Nodecoor = state.tfcoor(double(string(G.Nodes.Name)),:);%节点三维坐标
    LWidths = 5 * G.Edges.Weight/max(G.Edges.Weight);%按距离画粗细
    g = plot(G,'Xdata',Nodecoor(:,1),'Ydata',Nodecoor(:,2),'Zdata',Nodecoor(:,3),'EdgeLabel',{},...G.Edges.Weight,...
        'NodeLabel',{},'LineWidth',LWidths,'LineStyle','none','Marker','none','EdgeAlpha',0.4,'EdgeColor',[0 0.4470 0.7410]);%画蓝色
    GEdges = double(string(G.Edges.EndNodes));%边的首尾
    Edgeindex = all(abs(state.tfcoor(GEdges(:,1),:)-state.tfcoor(GEdges(:,2),:))< pro.Dim/2,2);%不跨越边界的边序号
    if ismember("SimpMol", pro.Plotlist)
    highlight(g,'Edges',Edgeindex,'LineStyle','-');%画不跨界的线
    end
    if ismember("Net", pro.Plotlist)
    highlight(g,'Edges',G.Edges.Style == 'Net' & Edgeindex,'EdgeColor',[0.9 0.3 0.1],'LineStyle','-');%网络键显示红色
    end
    %画网络环键
    if ismember("NetCircle", pro.Plotlist)
    CGEdges = unique(cell2mat(Edgecircle'));%环形边序号
    highlight(g,'Edges',CGEdges(Edgeindex(CGEdges)),'EdgeColor',[0.4940 0.1840 0.5560],'LineWidth',5,'LineStyle','-');%环形边显示紫色 
    end
else
    g = plot(graph([]));
end
%画强网络键
if numel(SG.Edges) > 0 %非空图
    SNodecoor = state.tfcoor(double(string(SG.Nodes.Name)),:);%节点三维坐标
    LWidths = 5 * SG.Edges.Weight/max(SG.Edges.Weight);%按距离画粗细
    sg = plot(SG,'Xdata',SNodecoor(:,1),'Ydata',SNodecoor(:,2),'Zdata',SNodecoor(:,3),'EdgeLabel',{},...
        'NodeLabel',{},'LineWidth',LWidths,'LineStyle','none','EdgeColor',[0.4660 0.6740 0.1880],'EdgeAlpha',0.6,'Marker','none');%画绿色
    if ismember("NetStrong",pro.Plotlist)%画强网络键
    SGEdges = double(string(SG.Edges.EndNodes));%边的首尾
    Edgeindex = all(abs(state.tfcoor(SGEdges(:,1),:)-state.tfcoor(SGEdges(:,2),:))< pro.Dim/2,2);%不跨越边界的边序号
    highlight(sg,'Edges',Edgeindex,'LineStyle','-');%画不跨界的线
    end
else
    sg = plot(graph([]));
end
axis equal
axis tight
view(-40,40);
hold off

a = nexttile(2);
cla
copyobj([result.g0,g,sg],a);%复制到新窗口
axis equal
axis tight
view([0,90])
drawnow limitrate
writeVideo(result.video{end}, getframe(f));%保存形貌图到video
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = net_init(result,pro)
if ~pro.Noplot
result.fig = figure("Name","Network Structure");
result.fig.OuterPosition = [50 50.6000 1400.4000 900.8000];
tiledlayout(result.fig,4,2,'TileSpacing','tight','Padding','compact');
VMorph = VideoWriter([pro.filepath, '\net']);
VMorph.FrameRate = 5;
open(VMorph);
result.video = VMorph;
result = morph(0,result,pro);%加入morph中的fig和video放在最后的图
end

result.stacname = {'NTree','Nmaxtree','Nave','Nave_err','Wave','Wave_err','Polydispersity','ReactDeg'...
    'ReactNodeDeg','NodeFunctionDegree','NetBranchRatio_ave','NetBranchRatio_err',...
    'Dangle_num','Dangle_ave','Dangle_err','Free_num','Free_ave','Free_err',...
    'Loop_num','Loop_ave','Loop_err','Net_num','Net_ave','Net_err',...
    'NetStrong_num','NetStrong_ave','NetStrong_err','NetCircle_num','NetCircle_ave','NetCircle_err'};
result.stacname_hist = ["Dangle","Free","Loop","Net","NetStrong","NetCircle"];
result.data_series = array2table(zeros(pro.Cyclenum,length(result.stacname)),'VariableNames',result.stacname);%series，表格格式
result.maxrow_hist = 100;%初步最大行数
result.data_hist = zeros(result.maxrow_hist,pro.Cyclenum*length(result.stacname_hist));%建立hist数组
end

function result = net_term(result,pro)
%程序终止任务
if ~pro.Noplot
cellfun(@(x)(close(x)),result.video)
end
end

function G = ord2simplegraph(state,pro)
%%%将ord的链替换为边的生成图
cat = state.ca;%临时存储ca
ca0 = sum(cat ~=0 ,2);%反应的官能度
chainmon = find(ca0 == 2);%中间链单元
chainend = find(ca0 ~= 2);%端部链单元

visited = zeros(state.TotN,1);%记录访问过的中间链单元，访问过（1），未访问（0）
Edge = zeros(pro.NBranch * state.TotN, 4);%记录链的两端、链长、链属性
TotL = 0;%总键数

%从中间点向端点或节点方向查找%
for i = 1:length(chainmon)
   n = chainmon(i);%中间点
   if visited(n) == 0
       [s,w1,cat,ca0,visited] = searchend(n,1,cat,ca0,visited);
       if s ~= n % 不是孤立圈
          [t,w2,cat,ca0,visited] = searchend(n,2,cat,ca0,visited);
       else%孤立圈
           t = n; w2 = 0;
       end
       visited(n) = 1;
       TotL = TotL + 1;
       if s == t 
            Edge(TotL,:) = [s, t, w1 + w2, 0];%loop
       elseif ca0(s) == 1 && ca0(t) == 1
            Edge(TotL,:) = [s, t, w1 + w2, 2];% free
       elseif (ca0(s) == 1 && ca0(t) > 2) || (ca0(t) == 1 && ca0(s) > 2)
            Edge(TotL,:) = [s, t, w1 + w2, 1];% dangle
       else
            Edge(TotL,:) = [s, t, w1 + w2, 3];% net
       end
   end
end

%查找端点或节点间的直连%
for i = 1:length(chainend)
    s = chainend(i);%起点
    if any(cat(s,:),2)%存在未清零的连接（不连中间点的连接）
        n = find(cat(s,:)~=0);
        for j = n
            t = cat(s,j);%终点
            TotL = TotL + 1;
            if ca0(s) > 2 && ca0(t) > 2 %两个节点直连
                Edge(TotL,:) = [s,t,1,3]; %short net
            elseif ca0(s) == 1 && ca0(t) == 1%两个端点直连
                Edge(TotL,:) = [s,t,1,2]; %free
            else%节点与端点或端点与端点直连
                Edge(TotL,:) = [s,t,1,1]; %dangle
            end
            cat(s,j) = 0;
            cat(t,cat(t,:) == s) = 0;         
        end
    end
end

%画将链等效为边的网络图，此图为多重图%
%Style = replace(string(Edge(1:TotL,4)),{'0','1','2','3'},{'Loop','Dangle','Free','Net'});%不同数字转为不同类型
Style = replace(string(Edge(1:TotL,4)),["0","1","2","3"],["Loop","Dangle","Free","Net"]);%不同数字转为不同类型
Nodes = unique(Edge(1:TotL,[1,2]),'sorted');%节点序号
EdgeTable = table(string(Edge(1:TotL,[1,2])),Edge(1:TotL,3),Style, ...
    'VariableNames',{'EndNodes' 'Weight' 'Style'});%边的属性，起点终点，长度，类型
NodeTable = table(string(Nodes(:)),'VariableNames',{'Name'});%节点属性，节点字符串
G = graph(EdgeTable,NodeTable);%建立图
end

function [n,w,cat,ca0,visited] = searchend(n,direction,cat,ca0,visited)
%%%沿着中间点向某个方向查找终点（循环点、三度点和端点）
    w = 0; n1 = n;
    while 1
        n0 = n;
        n = cat(n,direction);
        w = w + 1;
        if ca0(n)~=2 || n == n1
            break;
        end
        visited(n) = 1;
    end
    if ca0(n) ~= 2 %节点或端点
        cat(n ,cat(n,:) == n0) = 0;
    else %孤立圈
        visited(n) = 1;
    end
end

function [SG] = chain2edge(SG)
%%%利用邻接矩阵将图中链转化为边，去除图中环、孤立节点，重复边取小值，交联节点不变
Sa = adjacency(SG,'weighted');%利用邻接矩阵去除中间节点
Nodes = SG.Nodes;
Sd = degree(SG);
Sd2 = find(Sd == 2);%中间节点
for i = Sd2'
    [~,n,m] = find(Sa(i,:));
    if length(n) ~= 2 %已经成环
        continue;
    end
    s = n(1); t = n(2); w = sum(m);
    Sa(i, s) = 0; Sa (s, i) = 0;
    Sa(i, t) = 0; Sa (t, i) = 0;    
    if Sa(s, t) ~= 0 %已经存在边
        if Sd(s) > 2 && Sd(t) > 2 %支化点间多重边
           w = min(w, Sa(s, t));%选短边
        else%成环
           w = 0;%删除环
           %w = w + Sa(s,t);%如需保留环则用此三行
           %Sa(s, t) = 0; Sa(t, s) = 0;
           %t = s;
        end
    end
    Sa(s,t) = w; Sa(t, s) = w;
end
SG = graph(Sa);
SG.Nodes.Name = Nodes.Name;
SG = rmnode(SG,find(degree(SG)==0));%删除孤立点
end

