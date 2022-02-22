#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <bits/stdc++.h>
#include <iostream>
#include <random>
#include <chrono>
using namespace std;

// Command to Compile the c++ file: c++ -std=gnu++2a -Wall -g -O3 randmst.cc -o randmst

// Acknowledgement: We referenced https://www.geeksforgeeks.org/kruskals-minimum-spanning-tree-using-stl-in-c/
// for advice on the syntax used in representing the disjoint sets data structure.
// We implemented the logic of the kruskal algorithm and the data structures ourselves.
  
// Sets Disjoint data structure
struct SetsDisjoint
{
    int *parent;
    int *rank;
  
    SetsDisjoint(int n)
    {
        parent = new int[n+1];
        rank = new int[n+1];
  
        // Initially, all vertices are in different sets and their ranks are 0.
        for (int i = 0; i <= n; i++)
        {
            parent[i] = i; // Every element is it's own parent
            rank[i] = 0;
        }
    }
  
    // Union by rank heuristic
    void merge(int u, int v)
    {
        u = find(u);
        v = find(v);
  
        // Make the tree with lower rank into a subtree of the tree with higher rank
        if (rank[u] > rank[v])
            parent[v] = u;
        else
            parent[u] = v; // If rank[u] <= rank[v]
  
        if (rank[u] == rank[v])
            ++rank[v];
    }

    // Return the parent of vertex u
    // Implements the Path Compression heuristic
    int find(int u)
    {
        // Make the parent of the nodes in the path from u--> parent[u] point to parent[u]
        if (u != parent[u])
            parent[u] = find(parent[u]);
        return parent[u];
    }
};

struct Graph
{
    int V;
    vector< pair<double, pair<int, int>> > edges;
    
    Graph(int V)
    {
        this->V = V;
    }

    // Add an edge to the graph from node u to node v with weight w
    void addEdge(int u, int v, double w)
    {
        edges.push_back({w, {u, v}});
    }

    // Return the weight of the MST
    double kruskal()
    {
        double mst_weight = 0;
        // Sort edges in increasing order by edge weight
        sort(edges.begin(), edges.end());
        // Initialize disjoint set where each vertex in the graph is its own set
        SetsDisjoint sd(V);
    
        // Iterate through all sorted edges
        vector< pair<double, pair<int, int>> >::iterator it;
        for (it=edges.begin(); it!=edges.end(); it++)
        {
            int u = it->second.first;
            int v = it->second.second;
            int set_u = sd.find(u);
            int set_v = sd.find(v);
            if (set_u != set_v)
            {
                mst_weight += it->first;
                sd.merge(set_u, set_v);
            }
        }
        return mst_weight;
    }
};

Graph getGraph(int numpoints,uint dimension) {
    
    Graph g(numpoints);

    // Reseed the uniform_real_distribution random number generator
    random_device rd; 
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0,1.0);

    // Generate a graph for the 0-dimensional case
    if (dimension==0) {
        double edge_cutoff = 6*pow(numpoints,-0.919);
        for (int i=0; i<numpoints; ++i) {
            for (int j=i+1; j<numpoints; ++j) {
                    double weight = (double) dis(gen);
                    if (weight < edge_cutoff) {
                        g.addEdge(i,j,weight);
                    }   
            }
        }
        return g;
    }
    
    // Generate a 2-D array representing 2,3, or 4 dimensional coordinates for higher-dimensional cases
    double** arr = new double*[numpoints];
    for (int i = 0; i < numpoints; ++i) {
        arr[i] = new double[dimension];
        for (uint j = 0; j < dimension; ++j) {
            arr[i][j] = (double) dis(gen);
        }
    }

    // Add all n choose 2 edges
    double edge_cutoff=10;
    if (dimension==2) {
        edge_cutoff=2.6*pow(numpoints,-.514);
    }
    if (dimension==3) {
        edge_cutoff=2*pow(numpoints,-0.342);
    }
    if (dimension==4) {
        edge_cutoff=1.6*pow(numpoints,-0.243);
    }
    double weight = 0;
    for (int i=0; i<numpoints; ++i) {
        for (int j=i+1; j<numpoints; ++j) {
                weight = 0; // Euclidean distance between arr[i] and arr[j]
                if (dimension==2) {
                    weight = sqrt(pow(arr[i][0]-arr[j][0],2)+pow(arr[i][1]-arr[j][1],2));
                }
                if (dimension==3) {
                    weight = sqrt(pow(arr[i][0]-arr[j][0],2)+pow(arr[i][1]-arr[j][1],2)+pow(arr[i][2]-arr[j][2],2));
                }
                if (dimension==4) {
                    weight = sqrt(pow(arr[i][0]-arr[j][0],2)+pow(arr[i][1]-arr[j][1],2)+pow(arr[i][2]-arr[j][2],2)+pow(arr[i][3]-arr[j][3],2));
                }
                if (weight < edge_cutoff) {
                    g.addEdge(i,j,weight);
                }   
        }
    }
    // Free arr
    for(int i = 0; i < numpoints; ++i) { // Free each sub-array
        delete[] arr[i];   
    }
    delete[] arr; // Free the array of pointers
    return g;
}

double kruskal_trials(int numpoints, uint numtrials, uint dimension) {
    double mst_weight_sum = 0;
    for (uint i = 0; i < numtrials; ++i) {
        Graph g=getGraph(numpoints,dimension);
        
        double mst_weight = g.kruskal();
        mst_weight_sum+=mst_weight;
    }
    cout << mst_weight_sum/numtrials << endl;
    return mst_weight_sum/numtrials;
}

void produce_table_times(uint dimension) {
    ofstream outdata;
    int num_n_values_to_test = 12;
    int n_values[12] = {128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144};

    outdata.open("avg_mst_weight_dim"+to_string(dimension)+".dat", ios_base::app); // opens the file
    if( !outdata ) { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    for (int i = 0; i < num_n_values_to_test; ++i) {
        outdata << n_values[i] << ",";
    }
    outdata << endl;
    for (int i = 0; i < num_n_values_to_test; ++i) {
        outdata << i << "," << kruskal_trials(n_values[i],5,dimension);
    }
    outdata << endl;
    outdata.close();
    return;
    }

// ./kruskal 0 numpoints numtrials dimension
int main(int argc, char **argv) {
    assert(argc==5);
    int numpoints = atoi(argv[2]);
    uint flag=atoi(argv[1]),numtrials=atoi(argv[3]),dimension=atoi(argv[4]);
    (void) numtrials; // avoid uninitialized variable warnings
    if (flag==0) {
        kruskal_trials(numpoints, numtrials, dimension);
    }
    if (flag==1) { // run kruskal for one graph with the given numpoints, dimension
        Graph g=getGraph(numpoints,dimension);
        cout << g.kruskal() << endl;
    }
    if (flag == 2) { // generate all tables of values
        produce_table_times(0);
        produce_table_times(2);
        produce_table_times(3);
        produce_table_times(4);
    }
    if (flag == 5) { // generate just one table of values for the first argument
        produce_table_times(dimension);
    }
    
    return 0;

    
}
