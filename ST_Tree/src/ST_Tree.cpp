#include "ST_Tree.hpp"
#include <algorithm>
#include <fstream>
#include <queue>
#include <cmath>

// Constructors
ST_Tree::ST_Tree(bool optim, std::map<int, int>& treePar, int n, int debug) // Construct based on input tree/forest and number of nodes (named 1 to n)
{
    optimized = optim;
    debug_mode = 0;

    for (int i = 1; i <= n; i++)
    {
        vertices[i] = new ST_Node(true, i);
        dparent[i] = -1;
    }

    for (const auto &[u, v] : treePar)
    {
        link(u, v, 0);
    }

    debug_mode = debug;
    representation_number = 0;
}

ST_Tree::ST_Tree(bool optim, int n, int debug) // Create tree of n unconnected nodes (named 1 to n)
{ 
    optimized = optim;

    for (int i = 1; i <= n; i++)
        vertices[i] = new ST_Node(true, i), dparent[i] = -1;
    
    debug_mode = debug;
    representation_number = 0;
}

// destructor
ST_Tree::~ST_Tree()
{
    // First, collect all internal nodes to delete
    std::vector<ST_Node*> internal_nodes;
    
    // Get all unique paths
    std::vector<ST_Node*> paths = getAllUniquePaths();
    
    // For each path, collect all internal nodes
    for (ST_Node* p : paths)
    {
        std::queue<ST_Node*> bfs;
        bfs.push(p);
        while (!bfs.empty())
        {
            ST_Node* current = bfs.front();
            bfs.pop();
            
            if (!current->external)
            {
                internal_nodes.push_back(current);
                
                // Add children to queue if they're not external
                if (!current->bleft->external)
                    bfs.push(current->bleft);
                if (!current->bright->external)
                    bfs.push(current->bright);
            }
        }
    }
    
    // Delete all internal nodes
    for (ST_Node* node : internal_nodes)
    {
        // Set pointers to nullptr to avoid dangling pointers
        node->bparent = nullptr;
        node->bleft = nullptr;
        node->bright = nullptr;
        node->bhead = nullptr;
        node->btail = nullptr;
        delete node;
    }
    
    // Delete all external nodes (vertices)
    for (auto& pair : vertices)
    {
        // Set pointers to nullptr to avoid dangling pointers
        pair.second->bparent = nullptr;
        pair.second->bleft = nullptr;
        pair.second->bright = nullptr;
        pair.second->bhead = nullptr;
        pair.second->btail = nullptr;
        delete pair.second;
    }
    
    // Clear maps
    vertices.clear();
    dparent.clear();
    dcost.clear();
}

// Elemntary Path operations
// Static Operations
ST_Node* ST_Tree::path(int v) // Return the node representing the path v belongs to
{
    ST_Node* cur = vertices[v];
    while (cur->bparent) 
        cur = cur->bparent;
    return cur;
}

int ST_Tree::head(ST_Node* p) // Return head/lower-most node of path
{
    return (p->reversed) ? p->btail->vertex_id : p->bhead->vertex_id;
}

int ST_Tree::tail(ST_Node* p) // Return tail/upper-most node of path
{
    return (p->reversed) ? p->bhead->vertex_id : p->btail->vertex_id;
}


int ST_Tree::before(int v) { // returns the vertex before v on path(v), if v is the tail return -1
    ST_Node* u = vertices[v];

    // deepest node that is the right child of its parent
    ST_Node* deepest_right = nullptr;
    ST_Node* current = u;
    bool rev = false;
    if (current->bparent)
    {
        rev = get_reversal_state(current);
        while (current->bparent) 
        {
            rev ^= current->reversed;
            if ((current->bparent->bright == current) ^ rev) {  // if i am a right child
                deepest_right = current->bparent;
                break;
            }
            current = current->bparent;
        }
    }

    // rightmost external descendant of left sibling
    if (deepest_right) {
        if (rev)
            return head(deepest_right->bright); 
        else
            return tail(deepest_right->bleft);
    }

    return -1;
}

int ST_Tree::after(int v) { // returns the vertex after v on path(v), if v is the head return -1
    ST_Node* u = vertices[v];

    // deepest node that is the left child of its parent
    ST_Node* deepest_left = nullptr;
    ST_Node* current = u;
    bool rev = false;
    if (current->bparent)
    {
        rev = get_reversal_state(current);
        while (current->bparent) {
            rev ^= current->reversed; // update reversal state while moving up
            if ((current->bparent->bleft == current) ^ rev) {  // if i am a left child
                deepest_left = current->bparent;
                break;
            }
            current = current->bparent;
        }
    }

    // leftmost external descendant of right sibling
    if (deepest_left) {
        if (rev)
            return tail(deepest_left->bleft); 
        else
            return head(deepest_left->bright);
    }

    return -1;
}


double ST_Tree::pcost(int v){ // Return the cost of the edge between v, after(v) -> parent in path
    ST_Node* u = vertices[v];

    // deepest node that is the left child of its parent
    ST_Node* deepest_left = nullptr;
    ST_Node* current = u;
    bool rev = get_reversal_state(current->bparent);
    while (current->bparent) {
        if (!rev){
            if (current->bparent->bleft == current) {  
            deepest_left = current;
            break;
            }
        current = current->bparent;
        rev ^= current->reversed;
        }
        else{
            if (current->bparent->bright == current) {  
            deepest_left = current;
            break;
            }
        current = current->bparent;
        rev ^= current->reversed;
        }

    }

    return grosscost(deepest_left->bparent);
}

double ST_Tree::pmincost(ST_Node* p){ // Return the vertex closest to tail(p) such that (u, after(u)) has minimum cost among edges on p
    ST_Node* u = p;

    bool rev = get_reversal_state(u);
    while (u->netcost !=0 ||
             !((!rev &&
                 (u->bright->external ||
                  u->bright->netmin >0)) ||
                   (rev &&
                    (u->bleft->external ||
                     u->bleft->netmin > 0)))){
        if (!rev){
            if (!u->bright->external && u->bright->netmin==0){
                u = u->bright;
            }
            else if (u->netcost > 0){
                u = u->bleft;
            }
            rev ^= u->reversed;
        }
        else{
            if (!u->bleft->external && u->bleft->netmin==0){
            u = u->bleft;
            }
            else if (u->netcost > 0){
                u = u->bright;
            }
            rev ^= u->reversed;
        }
    }

    return (rev) ? head(u->bright) : tail(u->bleft);
}

void ST_Tree::pupdate(ST_Node* p, double x){ // Increase the cost of each in the path by x
    p->netmin = p->netmin + x;
}

void ST_Tree::reverse(ST_Node* p){ // reverse the path
    p->reversed = !p->reversed;
}

// Cost Helpers
double ST_Tree::grosscost(ST_Node* v){ // compute grosscost of a node (weight of the edge)
    return ((v->netcost) + (grossmin(v)));
}

double ST_Tree::grossmin(ST_Node* v){ // compute grossmin of a path (edge with least weight) represented by a node
  double min_value = v->netmin;

  // traverse upwards till root
  while (v->bparent != nullptr) {
        v = v->bparent;
        min_value += v->netmin;
  }

  return min_value;
}


// Dynamic Operations
// Suboperations
void ST_Tree::construct(ST_Node* v, ST_Node* w, double x) // Create a new node corresponding to an edge connecting the two paths represented by v and w
{
    ST_Node* u = new ST_Node(false, nullptr, false, 0, 0); // internal node, with no parent, not reversed, netmin and netcost 0
    
    // Handle left node
    u->bleft = v;
    v->bparent = u;
    u->bhead = vertices[head(v)];
    
    // Handle right node
    u->bright = w;
    w->bparent = u; 
    u->btail = vertices[tail(w)];

    // set my weight
    u->wt = v->wt + w->wt;

    // update costs
    // update netmin - min is min of myself and both my children
    u->netmin = x;
    if (!u->bleft->external)
        u->netmin = std::min(u->netmin, u->bleft->netmin);
    if (!u->bright->external)
        u->netmin = std::min(u->netmin, u->bright->netmin);

    // set netcost
    u->netcost = x - u->netmin;

    // update children - netmin is difference between original netmin and value assigned to parent
    if (!u->bleft->external)
        u->bleft->netmin -= u->netmin; 
    if (!u->bright->external)
        u->bright->netmin -= u->netmin; 
}

std::tuple<ST_Node*, ST_Node*, double> ST_Tree::destroy (ST_Node* u) // Destroy the root of the tree, breaking it into two trees/represented paths -> return (path1, path2, edge1cost, edge2cost) - edge1, edge2 are the edges broken
{
    std::tuple<ST_Node*, ST_Node*, double> result = {u->bleft, u->bright, 0};
    
    // Update costs of children - reverse of construct operation
    if (!u->bleft->external)
        u->bleft->netmin += u->netmin;
    if (!u->bright->external)
        u->bright->netmin += u->netmin;

    if (u->reversed) // if the node I am destroying is reversed, then reverse both children
    {
        if (!u->bright->external) u->bright->reversed = !u->bright->reversed;
        if (!u->bleft->external) u->bleft->reversed = !u->bleft->reversed;
    }

    // remove references to u and delete u
    u->bleft->bparent = nullptr;
    u->bright->bparent = nullptr;
    delete u;

    return result;
}

void ST_Tree::rotate(ST_Node* v) // rotate a node (representing an edge) upwards
{
    if (v->bparent)
    {
        ST_Node* u = v->bparent;

        if (v->reversed) // handle v's reversal by computing it down 1 stage
        {
            v->reversed = !v->reversed;

            // reverse my pointers
            ST_Node* tmp = v->bhead;
            v->bhead = v->btail;
            v->btail = tmp;
            
            tmp = v->bleft;
            v->bleft = v->bright;
            v->bright = tmp;

            // reverse my non-external children
            if (!v->bleft->external) v->bleft->reversed = !v->bleft->reversed;
            if (!v->bright->external) v->bright->reversed = !v->bright->reversed;
        }

        
        if (u->bright == v) // rotate left
        {
            // give v's left child - reverse already handled so my left child by pointer is my actual left child
            u->bright = v->bleft;
            u->bright->bparent = u;

            // update v's parent
            v->bparent = u->bparent;
            if (v->bparent)
            {
                if (v->bparent->bleft == u)
                    v->bparent->bleft = v;
                else
                    v->bparent->bright = v;
            }

            // make u left child
            v->bleft = u;
            u->bparent = v;

            // update costs
            double orig_v_nmin = v->netmin;
            // update v's right child
            if (!v->bright->external)
                v->bright->netmin += v->netmin; // as v->grossmin reduces by v->netmin
            // update v
            v->netcost += v->netmin; // grossmin reduces by netmin
            v->netmin = u->netmin;  // my grossmin is now u's grossmin and i now have the former parent of u so my netmin is u's netmin

            // update u
            // my netmin is minimum of my netcost and my children's netmin
            u->netmin = u->netcost;
            if (!u->bleft->external)
                u->netmin = std::min(u->netmin, u->bleft->netmin);
            if (!u->bright->external)
                u->netmin = std::min(u->netmin, u->bright->netmin);
            u->netcost -= u->netmin; // netmin represents change in grossmin after rotation

            // update rotated parents children
            if (!u->bleft->external)
                u->bleft->netmin -= u->netmin; // u acts as an intermediate between me and the smallest node
            if (!u->bright->external)
            {
                if (orig_v_nmin) // if my previous parent i.e. v, did not have the min value, so my parent must contain the minimum value - so the difference between me and my parent's min value increases by original netmin of v
                    u->bright->netmin += orig_v_nmin;
                else // if my previous parent had min value then my new parent acts as an intermediate between me and the smallest node
                    u->bright->netmin -= u->netmin;
            }
            
        }
        else // rotate right
        {
            // give v's right child - reverse already handled so my right child by pointer is my actual right child
            u->bleft = v->bright;
            u->bleft->bparent = u;

            // update v's parent
            v->bparent = u->bparent;
            if (v->bparent)
            {
                if (v->bparent->bleft == u)
                    v->bparent->bleft = v;
                else
                    v->bparent->bright = v;
            }
            
            // make u right child
            v->bright = u;
            u->bparent = v;

            // update costs
            double orig_v_nmin = v->netmin;
            // update v's right child
            if (!v->bleft->external) // as v->grossmin reduces by v->netmin
                v->bleft->netmin += v->netmin;
            // update v
            v->netcost += v->netmin; // grossmin reduces by netmin
            v->netmin = u->netmin; // my grossmin is now u's grossmin and i now have the former parent of u so my netmin is u's netmin

            // update u
            // my netmin is minimum of my netcost and my children's netmin
            u->netmin = u->netcost;
            if (!u->bleft->external)
                u->netmin = std::min(u->netmin, u->bleft->netmin);
            if (!u->bright->external)
                u->netmin = std::min(u->netmin, u->bright->netmin);
            u->netcost -= u->netmin; // netmin represents change in grossmin after rotation

            // update rotated parents children
            if (!u->bright->external)
                u->bright->netmin -= u->netmin; // u acts as an intermediate between me and the smallest node
            if (!u->bleft->external)
            {
                if (orig_v_nmin) // if my previous parent i.e. v, did not have the min value, so my parent must contain the minimum value - so the difference between me and my parent's min value increases by original netmin of v
                    u->bleft->netmin += orig_v_nmin;
                else // if my previous parent had min value then my new parent acts as an intermediate between me and the smallest node
                    u->bleft->netmin -= u->netmin;
            }
        }

        // update head and tail
        u->bhead = vertices[head(u->bleft)];
        u->btail = vertices[tail(u->bright)];

        v->bhead = vertices[head(v->bleft)];
        v->btail = vertices[tail(v->bright)];

        bool tmp = u->reversed; // if u was reversed, if i set v to be reversed, this has the same effect as the new subtree of v has the same nodes as the old subtree of u
        u->reversed = v->reversed;
        v->reversed = tmp;

        // update weights
        v->wt = u->wt;
        u->wt = u->bleft->wt + u->bright->wt;
    }
    if (debug_mode == 1)
    {
        displayInternalGraph();
        std::cout << "Rotated " << getEdge(v).first << "-" << getEdge(v).second << " | graph num " << representation_number << std::endl;
    }
}

ST_Node* ST_Tree::tilt_left(ST_Node* x) // for node rank management
{
    if ((x->bleft->rank == x->bright->rank) && (x->bright->rank == x->rank))
    {
        ++x->rank;
    }
    else if (!x->reversed && x->bright->rank == x->rank)
    {
        rotate(x->bright);
        return x->bparent;
    }
    else if (x->reversed && x->bleft->rank == x->rank)
    {
        rotate(x->bleft);
        return x->bparent;
    }
    return x;
}

ST_Node* ST_Tree::tilt_right(ST_Node* x) // for node rank management
{
    if ((x->bleft->rank == x->bright->rank) && (x->bright->rank == x->rank))
    {
        ++x->rank;
    }
    else if (!x->reversed && x->bleft->rank == x->rank)
    {
        rotate(x->bleft);
        return x->bparent;
    }
    else if (x->reversed && x->bright->rank == x->rank)
    {
        rotate(x->bright);
        return x->bparent;
    }
    return x;
}

// Operations
ST_Node* ST_Tree::concatenate(ST_Node* p, ST_Node* q, double x) // Connect two paths through an edge of cost x 
{
    std::vector<int> conc_paths = {getEdge(p).first, getEdge(p).second, getEdge(q).first, getEdge(q).second};

    // unoptimized
    if (!optimized)
    {
        construct(p, q, x);
        p = p->bparent;
    }

    // optimized
    else
    {
        if ((p->rank == q->rank) || (p->rank > q->rank && p->external) || (p->rank < q->rank && q->external)) // Case 1
        {
            construct(p, q, x);
            p->bparent->rank = std::max(p->rank, q->rank) + 1;
            p = p->bparent;
        }
        else if (p->rank > q->rank && !(p->external)) // Case 2
        {
            p = tilt_left(p);
            
            ST_Node* y = (!p->reversed) ? p->bleft : p->bright;
            ST_Node* z = (p->reversed) ? p->bleft : p->bright;
            double c = p->netcost + p->netmin;
            int r = p->rank;

            // delete node and rebuild it later to maintain cost and reversal state
            destroy(p); // delete p

            ST_Node* u = concatenate(z, q, x); // local join (z, q)

            construct(y, u, c); // reconstruct p
            p = y->bparent;
            p->rank = r;
        }
        else if (p->rank < q->rank && !(q->external)) // Case 3
        {
            q = tilt_right(q);
            
            ST_Node* y = (!q->reversed) ? q->bleft : q->bright;
            ST_Node* z = (q->reversed) ? q->bleft : q->bright;
            double c = q->netcost + q->netmin;
            int r = q->rank;

            // delete node and rebuild it later to maintain cost and reversal state
            destroy(q); // delete q

            ST_Node* u = concatenate(p, y, x); // local join (p, y)

            construct(u, z, c); // reconstruct q

            q = z->bparent;
            q->rank = r;
            p = q;
        }
    }

    if (debug_mode == 1 )
    {
        std::cout << "concatenated " << conc_paths[0] << "-" << conc_paths[1] << " with " << conc_paths[2] << "-" << conc_paths[3] << " | graph num : " << representation_number+1 << std::endl;
        displayInternalGraph();
    }
    return p;
}

std::tuple<ST_Node*, ST_Node*, double, double> ST_Tree::split(int v) // Break a path at a node into path before and path after. Selected node becomes a trivial/single node path
{
    ST_Node* vNode = vertices[v];
    ST_Node *p {nullptr}, *q {nullptr}; // left / right or below / above
    double x {-1}, y {-1}; // cost of edge to left / right or below / above

    // unoptimized
    if (!optimized)
    {
        while (vNode->bparent) // while there are edges connected to me in the path
        {
            ST_Node* cur = vNode->bparent; // Find an edge connected to me
            while (cur->bparent) // rotate edge to top to prepare for deletion
                rotate(cur); 

            // store whether left edge or right edge in path
            if ((!cur->reversed && tail(cur->bleft) == v) || (cur->reversed && head(cur->bright) == v)) // if i am in left, then the path being separated is from me to above - after(v) to tail(path)
                q = (cur->reversed) ? cur->bleft : cur->bright, y = cur->netcost + cur->netmin;
            else if ((!cur->reversed && head(cur->bright) == v) || (cur->reversed && tail(cur->bleft) == v)) // if i am in right, then the path being separated is from me to bottom - from head(path) to before(v)
                p = (cur->reversed) ? cur->bright : cur->bleft, x = cur->netcost + cur->netmin;

            destroy(cur); // delete edge
        }
    }

    // optimized
    else
    {
        if (vNode->bparent)
        {
            std::queue<ST_Node*> pList;
            std::queue<double> pCosts;
            std::queue<ST_Node*> qList;
            std::queue<double> qCosts;
            std::vector<ST_Node*> parList;


            parList.push_back(vNode);
            vNode = vNode->bparent; // ignore external node
            double curMin = grossmin(vNode);
            bool rev = get_reversal_state(vNode);
            while (vNode)
            {
                if ((vNode->bright == parList.back()) ^ rev) // if i was a right child
                {
                    if (pList.size())
                        pCosts.push(vNode->netcost + curMin);
                    else
                        x = vNode->netcost + curMin;
                    pList.push((!rev) ? vNode->bleft : vNode->bright);
                }
                else // if i was a left child
                {
                    if (qList.size())
                        qCosts.push(vNode->netcost + curMin);
                    else
                        y = vNode->netcost + curMin;
                    qList.push((!rev) ? vNode->bright : vNode->bleft);

                }

                parList.push_back(vNode);
                curMin -= vNode->netmin;
                rev ^= vNode->reversed;
                vNode = vNode->bparent;
            }

            while (parList.size() > 1) // ignore external node - that's why we stop at 1 
            {
                destroy(parList.back());
                parList.pop_back();
            }

            if (pList.size())
            {
                p = pList.front();
                pList.pop();
                while (pList.size())
                {
                    p = concatenate(pList.front(), p, pCosts.front());
                    pList.pop();
                    pCosts.pop();
                }
            }

            if (qList.size())
            {
                q = qList.front();
                qList.pop();
                while (qList.size())
                {
                    q = concatenate(q, qList.front(), qCosts.front());
                    qList.pop();
                    qCosts.pop();
                }
            }
        }
    }
    
    if (debug_mode == 1)
    {
        displayInternalGraph();
        std::cout << "Split " << v << " completed in " << representation_number << std::endl;
    }
    return {p, q, x, y};
}

// Path Partition functions - bold and dashed edges
ST_Node* ST_Tree::splice(ST_Node* p) // Extend current bold path upwards by converting one dashed edge 
{
    int v = dparent[tail(p)]; // find node above me outside my path
    auto [q, r, x, y] = split(v); // r is path to root
    // note that weight of a path is the same as the size of tail(p) since it is essentially the sum of the weights of all the vertices in the path - this creates a telescoping sum that leaves us with size(tail(p))
    vertices[v]->wt -= p->wt; // subtract size of node below v (tail(p)) 
    vertices[v]->rank = (int) log2(vertices[v]->wt);
    if (q) // q is other downward path
    {
        dparent[tail(q)] = v; // make it dashed
        dcost[tail(q)] = x;
        vertices[v]->wt += q->wt; // add the size of the node previously connected to v
    }

    vertices[v]->rank = (int) log2(vertices[v]->wt); // update rank of v

    dparent[tail(p)] = -1; // dashed edge removed
    p = concatenate(p, path(v), dcost[tail(p)]); // connect me to the node above me
    if (r) // if more nodes on way to root
    {
        if (debug_mode ==1)
        {
            std::cout << "Splice " << getEdge(p).first << "-" << getEdge(p).second << " completed in " << representation_number+1 << std::endl;
        }
        return concatenate(p, r, y); // connect me to the path above v
    }
    else // no path to root
    {
        if (debug_mode ==1)
        {
            std::cout << "Splice " << getEdge(p).first << "-" << getEdge(p).second << " completed in " << representation_number << std::endl;
        }
        return p;
    }
}

ST_Node* ST_Tree::expose(int v) // Create bold path from this node to root of tree 
{
    auto [q, r, x, y] = split(v); // cut me off from anything below
    if (q) // if path afer me
    {
        dparent[tail(q)] = v; // then make me its dashed parent
        dcost[tail(q)] = x;
        vertices[v]->wt += q->wt; // add the size of the node previously connected to v
    }
    vertices[v]->rank = (int) log2(vertices[v]->wt); // update rank of v

    // connect me upwards
    ST_Node* p; 
    if (r)
        p = concatenate(path(v), r, y);
    else
        p = path(v);
        
    // if not connected to root then keep connecting upwards
    while (dparent[tail(p)] != -1)
        p = splice(p);
    
    if (debug_mode == 1)
    {
        std::cout << "exposed " << v << "| graph num : " << representation_number+1 << std::endl;
        displayInternalGraph();
    }
    return p;
}


// Helper Functions
// Functions related to visualisation
std::vector<ST_Node*> ST_Tree::getAllUniquePaths(){
    std::map <ST_Node*, int> paths;
    for ( int i = vertices.size(); i > 0; i--){
        ST_Node* p = path(i);
        if (paths.find(p) == paths.end())
            paths[p] = 1;
    };
    std::vector<ST_Node*> uniquePaths;
    for (auto const& [key, val] : paths)
        uniquePaths.push_back(key);
    return uniquePaths;
}

std::vector<std::vector<int>> ST_Tree::getAllGraphs(){
    std::vector<std::vector<int>> graphs;
    std::vector<ST_Node*> paths = getAllUniquePaths();
    for (ST_Node* p : paths){
        std::vector<int> graph;
        int cur_val = head(p);
        while (cur_val != tail(p))
        {
            graph.push_back(cur_val);
            cur_val = after(cur_val);
        }
        graph.push_back(cur_val);
        std::reverse(graph.begin(), graph.end());
        graphs.push_back(graph);
    }
    return graphs;
}

// Miscellaneous operations
std::pair<int, int> ST_Tree::getEdge(ST_Node* eNode)
{
    return {tail(eNode->bleft), head(eNode->bright)};
}

bool ST_Tree::get_reversal_state(ST_Node* v)
{
    bool rev = v->reversed;
    while (v->bparent)
    {
        v = v->bparent;
        rev = rev ^ v->reversed;
    }
    return rev;
}

// Dynamic Tree Operations
int ST_Tree::parent(int v){ // return parent of v (or return null if it has no parent)
    if (tail(path(v)) == v){ //so no parent or dashed parent
        return dparent[v];
    }

    return after(v);
}

int ST_Tree::root(int v){ // return root of the tree containing v
    return tail(expose(v));
}

double ST_Tree::cost(int v){ // Return cost of the edge (v, par(v))
    if (v==tail(path(v))){
        return dcost[v];
    }
    else{
        return pcost(v);
    }
}


double ST_Tree::mincost(int v){ // Return the node with the least cost from v to root(v)
    return pmincost(expose(v));
}

void ST_Tree::update(int v, double x){ // Increase weight of all edges on the path from v to root(v) by x
    return pupdate(expose(v), x);
}

void ST_Tree::link(int v, int w, double x) // Let v be a root of a tree. Connect w to v, effectively joining two trees  
{
    concatenate(path(v), expose(w), x); // makes path from root of v to w bold
}

double ST_Tree::cut(int v) // Divide the tree into two by breaking at vertex v  
{
    expose(v); // make path from me to root bold and everything below dashed
    auto [p, q, x, y] = split(v);
    dparent[v] = -1; // don't connect me to what I split off from
    return y;
}

void ST_Tree::evert(int v){ // Make v the root of the tree
    reverse(expose(v)); // connect me to root and then reverse the whole connection
    dparent[v] = -1;
}

// Visualisation Helper Functions
std::vector<std::vector<int>> ST_Tree::getAllEdges(){
    std::vector<std::vector<int>> edges;
    std::vector<std::vector<int>> graphs = getAllGraphs();
    // std::cout << "Internal Paths:\n";
    for (std::vector<int> graph : graphs){
        // std::cout << graph[0];
        for (int i = 0; i < (int) graph.size() - 1; i++){
            // std::cout << " " << graph[i+1];
            std::vector<int> edge = {graph[i], graph[i+1], (int) cost(graph[i+1])};
            edges.push_back(edge);
        }
        // std::cout << std::endl;
    }
    return edges;
};

std::vector<std::vector<int>> ST_Tree::getAllDashEdges(){
    std::vector<std::vector<int>> edges;
    for (const auto &[u, v] : dparent){
        if (v!=-1){
            std::vector<int> edge = {v, u, (int) dcost[u]};
            edges.push_back(edge);
        }
    }
    return edges;
};


// Internal Graph Visualiser
void ST_Tree::displayInternalGraph(int mode) {
    std::vector<ST_Node*> path_list = getAllUniquePaths(); 
    std::ofstream dot_file("internal_graph.dot");
    dot_file << "digraph G {\n";
    
    // Add Connections for each path
    for (ST_Node* p : path_list)
    {
        std::queue<ST_Node*> bfs;
        bfs.push(p);
        while (bfs.size())
        {
            ST_Node* cur = bfs.front();
            bfs.pop();
            if (!cur->external)
            {
                std::pair<int, int> edgeName = getEdge(cur);
                dot_file << "\"(" << edgeName.first << "," << edgeName.second <<  ((cur->reversed) ? ")*" : ")") << "\" [shape=box" << ((optimized) ? ",xlabel=" + std::to_string(cur->rank) : "") << "];\n"; // create internal node shape
                // left child connection
                if (cur->bleft->external) // if external
                {
                    dot_file << cur->bleft->vertex_id << " [shape=circle" << ((optimized) ? ",xlabel=" + std::to_string(cur->bleft->rank) : "") << "];\n";
                    dot_file << "\"(" << edgeName.first << "," << edgeName.second <<  ((cur->reversed) ? ")*" : ")") << "\" -> " << cur->bleft->vertex_id << " [label=l];\n";
                }
                else // if internal
                {
                    std::pair<int, int> childEdgeName = getEdge(cur->bleft);
                    dot_file << "\"(" << edgeName.first << "," << edgeName.second <<  ((cur->reversed) ? ")*" : ")") << "\" -> " << "\"(" << childEdgeName.first << "," << childEdgeName.second <<  ((cur->bleft->reversed) ? ")*" : ")") << "\" [label=l];\n";
                    bfs.push(cur->bleft);                    
                }

                // right child connection
                if (cur->bright->external) // if external
                {
                    dot_file << cur->bright->vertex_id << " [shape=circle" << ((optimized) ? ",xlabel=" + std::to_string(cur->bright->rank) : "") << "];\n";
                    dot_file << "\"(" << edgeName.first << "," << edgeName.second <<  ((cur->reversed) ? ")*" : ")") << "\" -> " << cur->bright->vertex_id << " [label=r];\n";
                }
                else // if internal
                {
                    std::pair<int, int> childEdgeName = getEdge(cur->bright);
                    dot_file << "\"(" << edgeName.first << "," << edgeName.second <<  ((cur->reversed) ? ")*" : ")") << "\" -> " << "\"(" << childEdgeName.first << "," << childEdgeName.second <<  ((cur->bright->reversed) ? ")*" : ")") << "\" [label=r];\n";
                    bfs.push(cur->bright);                    
                }
                
                // head and tail
                dot_file << "\"(" << edgeName.first << "," << edgeName.second <<  ((cur->reversed) ? ")*" : ")") << "\" -> " << cur->bhead->vertex_id << " [style=dashed, label = h];\n";
                dot_file << "\"(" << edgeName.first << "," << edgeName.second <<  ((cur->reversed) ? ")*" : ")") << "\" -> " << cur->btail->vertex_id << " [style=dashed, label = t];\n";
            }
            else
                dot_file << cur->vertex_id << " [shape=circle" << ((optimized) ? ",xlabel=" + std::to_string(cur->rank) : "") << "];\n";
        }
    }

    
    dot_file << "}\n";

    // Close the file
    dot_file.close();
    std::string command = "dot -Tpng internal_graph.dot -o internal_graphs/i" + std::to_string(++representation_number) + ".png";
    system(command.c_str());
    if (mode == 2)
    {
        command = "open internal_graphs/i" + std::to_string(representation_number)+ ".png";
        system(command.c_str());
    }
}



// Getter for debugging
ST_Node* ST_Tree::get_vertex_ptr(int vertex_id)
{
    return vertices[vertex_id];
}