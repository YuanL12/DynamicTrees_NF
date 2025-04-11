#pragma once

#include <map>
#include <set>
#include <vector>
#include <iostream>


struct ST_Node // Nodes of underlying tree - can represent vertices or edges of original tree. Their subtree represents a path in original tree. Single Node trees represent trivial/single node paths of original tree.
{
    bool external; // true/false
    int vertex_id; // vertex number
    ST_Node* bparent {nullptr};
    ST_Node *bleft {nullptr}, *bright {nullptr};
    ST_Node *bhead {nullptr}, *btail {nullptr}; // bottom-most and top-most nodes in path 

    // for optimization
    int wt;
    int rank;
    
    bool reversed; //to determine if reverse has occured

    // External node/vertex constructor
    ST_Node(bool ext, int vert) : external{ext}, vertex_id{vert}, bparent{nullptr}, bleft{this}, bright{this}, bhead{this},  btail{this}, wt{1}, rank{0}, reversed{false} {}; 

    // to determine cost of node + min edge
    double netmin = 0; 
    double netcost = 0;

    // Internal node/edge constructor
    ST_Node (bool ext, ST_Node* par, bool rev, double netmi, double netcst) : external{ext}, vertex_id{-1}, bparent{par}, rank{-1}, reversed{rev}, netmin{netmi}, netcost{netcst} {};
};

class ST_Tree
{
    private:
        // Tree Base Attributes
        std::map<int, ST_Node*> vertices; // A map from vertex number to its corresponding node in path tree
        std::map<int, int> dparent; // dashed edges/parents map
        std::map<int, double> dcost; // cost of dashed edge

        // Elementary Path operations
        // Static Operations
        ST_Node* path(int v);  // Return the node representing the path v belongs to
        int head(ST_Node* p);  // Return head/lower-most node of path
        int tail(ST_Node* p);  // Return tail/upper-most node of path

        int before(int v);  // return the vertex right before v in the path, unless it is tail in which case return -1
        int after(int v);  // return the vertex right after v in the path, unless it is head in which case return -1

        double pcost(int v);  // Return the cost of the edge between v, after(v) -> parent in path
        double pmincost(ST_Node* p);  // Return the vertex closest to tail(p) such that (u, after(u)) has minimum cost among edges on p
        void pupdate(ST_Node* p, double x);  // Increase the cost of each in the path by x
        void reverse(ST_Node* p);  // reverse the path

        // Cost Helpers
        double grosscost(ST_Node* v); // compute grosscost of a node (weight of the edge)
        double grossmin(ST_Node* v); // compute grossmin of a path (edge with least weight) represented by a node

        // Dynamic Operations
        // Suboperations
        void construct(ST_Node* v, ST_Node* w, double x);  // Create a new node corresponding to an edge connecting the two paths represented by v and w  
        std::tuple<ST_Node*, ST_Node*, double> destroy (ST_Node* u);  // Destroy the root of the tree, breaking it into two trees/represented paths 
        
        void rotate(ST_Node* v); // rotate a node (representing an edge) upwards
        ST_Node* tilt_left(ST_Node* x); // for node rank management
        ST_Node* tilt_right(ST_Node* x); // for node rank management

        // Operations
        ST_Node* concatenate(ST_Node* p, ST_Node* q, double x); // Connect two paths through an edge of cost x
        std::tuple<ST_Node*, ST_Node*, double, double> split(int v); // Break a path at a node into path before and path after. Selected node becomes a trivial/single node path 

        // Path Partition functions - bold and dashed edges
        ST_Node* splice(ST_Node* p); // Extend current bold path upwards by converting one dashed edge
        ST_Node* expose(int v); // Create bold path from this node to root of tree

        // Helper Functions
        // Functions related to visualisation
        std::vector<ST_Node*> getAllUniquePaths(); // Return all unique paths in the tree
        std::vector<std::vector<int>> getAllGraphs(); // Return all connected components of the tree

        // Miscellaneous operations
        std::pair<int, int> getEdge(ST_Node* eNode); // Get the vertices that this node connects
        bool get_reversal_state(ST_Node* v); // Calculate reversal state of the node
        
        int representation_number; // For numbering visualisation images
        int debug_mode; // Flag to create Visualisations in suboperations
        bool optimized; // Flag for using optimization

    public:
        // Constructors and Destructor
        ST_Tree(bool optim, std::map<int, int>& treePar, int n, int debug);  // Construct based on input tree/forest and number of nodes (named 1 to n)
        ST_Tree(bool optim, int n, int debug);
        
         // Create tree of n unconnected nodes (named 1 to n)
        ~ST_Tree(); // To implement


        // Dynamic Tree Operations
        int parent(int v); // return parent of v (or return null if it has no parent)
        int root(int v); // return root of the tree containing v

        double cost(int v); // Return cost of the edge (v, par(v))
        double mincost(int v); // Return the node with the least cost from v to root(v)
        void update(int v, double x); // Increase weight of all edges on the path from v to root(v) by x

        void link(int v, int w, double x); // Let v be a root of a tree. Connect w to v, effectively joining two trees
        double cut(int v); // Divide the tree into two by breaking at vertex v

        void evert(int v); // Make v the root of the tree

        // Visualisation Helper Functions
        std::vector<std::vector<int>> getAllEdges(); // Return all edges in the tree
        std::vector<std::vector<int>> getAllDashEdges(); // Return all dashed edges in the tree
        // std::vector<std::vector<int>> getAllDashEdges(int n); // Return all dashed edges in the tree
        
        // Internal Graph Visualiser
        void displayInternalGraph(int mode = 0);

        // Getter for debugging
        ST_Node* get_vertex_ptr(int vertex_id); // Return node for a given vertex - handy for debugging
        
        // void bfs(int source); // Perform BFS on the tree
        // void dfs(ST_Node* u, int t, double flow); // Perform DFS on the tree
};