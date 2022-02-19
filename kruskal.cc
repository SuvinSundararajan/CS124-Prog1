#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <bits/stdc++.h>
#include <iostream>
#include <random>
#include <chrono>
using namespace std;

// Compile: c++ -std=gnu++2a -Wall -g -O3 kruskal.cc -o kruskal

// Acknowledgement: Reference https://www.geeksforgeeks.org/kruskals-minimum-spanning-tree-using-stl-in-c/ for how to represent a graph

// Creating shortcut for an integer pair
typedef  pair<int, int> iPair;
  
// Structure to represent a graph
struct Graph
{
    int V, E;
    vector< pair<double, iPair> > edges;
  
    // Constructor
    Graph(int V, int E)
    {
        this->V = V;
        this->E = E;
    }
  
    // Utility function to add an edge
    void addEdge(int u, int v, double w)
    {
        edges.push_back({w, {u, v}});
    }

    double kruskal();
};
  
// To represent Disjoint Sets
struct DisjointSets
{
    int *parent, *rank;
    int n;
  
    DisjointSets(int n)
    {
        this->n = n;
        parent = new int[n+1];
        rank = new int[n+1];
  
        // Initially, all vertices are in
        // different sets and have rank 0.
        for (int i = 0; i <= n; i++)
        {
            rank[i] = 0;
  
            //every element is parent of itself
            parent[i] = i;
        }
    }
  
    // Find the parent of a node 'u'
    // Path Compression
    int find(int u)
    {
        /* Make the parent of the nodes in the path
           from u--> parent[u] point to parent[u] */
        if (u != parent[u])
            parent[u] = find(parent[u]);
        return parent[u];
    }
  
    // Union by rank
    void merge(int x, int y)
    {
        x = find(x), y = find(y);
  
        /* Make tree with smaller height
           a subtree of the other tree  */
        if (rank[x] > rank[y])
            parent[y] = x;
        else // If rank[x] <= rank[y]
            parent[x] = y;
  
        if (rank[x] == rank[y])
            rank[y]++;
    }
};

double Graph::kruskal()
{
    double mst_weight = 0;
    // Sort edges in increasing order by edge weight
    sort(edges.begin(), edges.end());
    // Initialize disjoint set
    DisjointSets ds(V);
  
    // Iterate through all sorted edges
    vector< pair<double, iPair> >::iterator it;
    for (it=edges.begin(); it!=edges.end(); it++)
    {
        int u = it->second.first;
        int v = it->second.second;
  
        int set_u = ds.find(u);
        int set_v = ds.find(v);
        if (set_u != set_v)
        {
            mst_weight += it->first;
            ds.merge(set_u, set_v);
        }
    }
    return mst_weight;
}

Graph getGraph(uint numpoints,uint dimension) {
    
    uint num_edges = numpoints*(numpoints-1)/2;
    Graph g(numpoints, num_edges);

    random_device rd; 
    mt19937 gen(rd()); 
    uniform_real_distribution<> dis(0,1.0);

    double arr[numpoints][dimension];
    if (dimension >= 2) {
        for (uint i = 0; i <numpoints; i++){
            for(uint k = 0; k < dimension; k++){
                arr[i][k] = (double) dis(gen);
            }
        }
    }
    
    // Add all n choose 2 edges
    for (uint i=0; i<numpoints; ++i) {
        for (uint j=i+1; j<numpoints; ++j) {
            if (i!=j) {
                double weight = 0; // Euclidean distance between arr[i] and arr[j]
                if (dimension==0) {
                    weight = (double) dis(gen);
                }
                if (dimension==2) {
                    weight = sqrt(pow(arr[i][0]-arr[j][0],2)+pow(arr[i][1]-arr[j][1],2));
                }
                if (dimension==3) {
                    weight = sqrt(pow(arr[i][0]-arr[j][0],2)+pow(arr[i][1]-arr[j][1],2)+pow(arr[i][2]-arr[j][2],2));
                }
                if (dimension==4) {
                    weight = sqrt(pow(arr[i][0]-arr[j][0],2)+pow(arr[i][1]-arr[j][1],2)+pow(arr[i][2]-arr[j][2],2)+pow(arr[i][3]-arr[j][3],2));
                }
                g.addEdge(i,j,weight);
            }
        }
    }
    return g;
}

double kruskal_trials(uint numpoints, uint numtrials, uint dimension) {
    double mst_weight_sum = 0;
    for (uint i = 0; i < numtrials; ++i) {
        Graph g=getGraph(numpoints,dimension);
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        mst_weight_sum+=g.kruskal();
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << numpoints << " points, dimension " << dimension << ", " << "Time difference: = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    
    }
    cout << "\nMST Weight: " << mst_weight_sum/numtrials << endl;
    return mst_weight_sum/numtrials;
}

void produce_table_times(uint dimension) {
    ofstream outdata; // outdata is like cin
    int num_n_values_to_test = 8;
    int n_values[num_n_values_to_test] = {128, 256, 512, 1024, 2048, 4096, 8192, 16384};
    //int n_values[12] = {128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144};

    outdata.open("runtimes_dim"+to_string(dimension)+".dat"); // opens the file
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

// ./randmst 0 numpoints numtrials dimension
int main(int argc, char **argv) {
    cout << argc << endl;
    assert(argc==5);
    uint flag=atoi(argv[1]),numpoints=atoi(argv[2]),numtrials=atoi(argv[3]),dimension=atoi(argv[4]);

    if (flag==1) {
        Graph g=getGraph(numpoints,dimension);
        mst_weight_sum+=g.kruskal();
    }
    if (flag == 2) {
        produce_table_times(0);
        produce_table_times(2);
        produce_table_times(3);
        produce_table_times(4);
    }
    else {
        produce_table_times(dimension);
    }
    
    return 0;

    
}
