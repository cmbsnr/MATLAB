G = graph([1 2 2 11 10 3 3 5 8 4 4 6],[2 11 10 3 5 5 4 8 9 7 6 7]);
% plot(G)
%%%利用邻接矩阵将图中链转化为边，去除图中环、孤立节点，重复边取小值，交联节点不变
Sa = adjacency(G,'weighted');%利用邻接矩阵去除中间节点
Nodes = G.Nodes;
Sd = degree(G);
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
           %w = 0;%删除环
           w = w + Sa(s,t);%如需保留环则用此三行
           Sa(s, t) = 0; Sa(t, s) = 0;
           t = s;
        end
    end
    Sa(s,t) = w; Sa(t, s) = w;
end
SG = graph(Sa);
% SG.Nodes.Name = Nodes.Name;
SG = rmnode(SG,find(degree(SG)==0));%删除孤立点
h = plot(SG);
weights = SG.Edges.Weight;
labeledge(h,1:numedges(SG),weights)