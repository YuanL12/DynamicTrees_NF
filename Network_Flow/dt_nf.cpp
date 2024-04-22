#include "dt_nf.hpp"

 

std::vector<std::map<int, int>> adj_inv(std::vector<std::map<int, int>> adj){
    int n = adj.size();
    std::vector<std::map<int, int>> inv(n);
    for(int v = 1; v < n; v++){
        for (const auto& [u, c]: adj[v]){
            inv[u][v] = c;
        }
    }
    return inv;
};

int dinicMaxFlow(int s, int t, 
std::vector<std::map<int, int>> adj) {
    int n = adj.size(); //the number of nodes
    GraphManager* g = new  GraphManager(n);
    ST_Tree* tree = new ST_Tree(false, n-1, 0);
    std::vector<std::map<int, int>> inv = adj_inv(adj);
    std::vector<std::map<int, int>> adj_copy = adj_inv(inv);
    int flow = 0;
    while (true){
        int v = tree->root(s);
        if (v==t){
            v = tree->mincost(s);
            int c = tree->cost(v);
            tree->update(s, -c);
            flow += c;
            g->displayCombinedGraph(tree->getAllEdges(), tree->getAllDashEdges(), "flow", 0);
            while(true){
                v = tree->mincost(s);
                if (tree->cost(v) == 0)
                    tree->cut(v);
                else
                    break;
                g->displayCombinedGraph(tree->getAllEdges(), tree->getAllDashEdges(), "flow", 0);
            }
        }
        else if (adj_copy[v].begin() == adj_copy[v].end()){
            if (v==s){
                break;
            }
            else{
                for (const auto& [u, c]: inv[v]){
                    if (v == tree->parent(u)){
                        tree->cut(u);
                        g->displayCombinedGraph(tree->getAllEdges(), tree->getAllDashEdges(), "flow", 0);
                    }
                    adj[u].erase(v);
                }
                inv[v].clear();
            }
        }
        else
        {
            int w;
            int c;
            for (const auto& [a, b] : adj_copy[v])
            {
                if (tree->root(a) != v)
                {
                    w = a;
                    c = b;
                    break;
                }
            }
            adj_copy[v].erase(w);
            tree->link(v, w, c);
            g->displayCombinedGraph(tree->getAllEdges(), tree->getAllDashEdges(), "flow", 0);
        }
    }
    return flow;
}


