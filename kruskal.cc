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

/**
Note: Time to compute 30,000 points = 80 seconds
Runtime increases by about 4x for each 2x increase in the number of points
I expect ~320 seconds for each trial of 65536, ~1280 seconds for 131072, and ~5120 seconds for 262144
Running 5 trials on 4 different dimensions for 262144 may take approximately 5120/60*5*4/60 = 28 hours
    Add ~7 hours for running on 131072 points, ~2 hours for running on 65536 points, ~1 hours for running on the rest
    Total estimated time to run: ~38 hours
*/

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

    ~Graph() {
        //delete &edges;
    }
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
    double max_edge_weight_in_mst = 0;
  
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
            max_edge_weight_in_mst = std::max(max_edge_weight_in_mst,it->first);
            mst_weight += it->first;
            ds.merge(set_u, set_v);
        }
    }
    cout << "-------------MST's max edge weight: " << max_edge_weight_in_mst << "---------" << endl;
    return mst_weight;
}

Graph getGraphThrowAwaySomeEdges(int numpoints,uint dimension) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    int num_edges = numpoints*(numpoints-1)/2;
    Graph g(numpoints, num_edges);

    random_device rd; 
    mt19937 gen(rd()); 
    uniform_real_distribution<> dis(0,1.0);

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
        // Note: without the edge weight cutoff, kruskal takes ~85 seconds for MST with 32768 points of dimension 3
        // Note: with the edge weight cutoff, kruskal takes ~35 milliseconds for MST with 32768 points of dimension 3
        edge_cutoff=2*pow(numpoints,-0.342);
    }
    if (dimension==4) {
        // Note: without the edge weight cutoff, kruskal takes ~85 seconds for MST with 32768 points of dimension 3
        // Note: with the edge weight cutoff, kruskal takes ~35 milliseconds for MST with 32768 points of dimension 3
        edge_cutoff=1.6*pow(numpoints,-0.243);
    }
    for (int i=0; i<numpoints; ++i) {
        for (int j=i+1; j<numpoints; ++j) {
                double weight = 0; // Euclidean distance between arr[i] and arr[j]
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
    // Free each sub-array
    for(int i = 0; i < numpoints; ++i) {
        delete[] arr[i];   
    }
    // Free the array of pointers
    delete[] arr;

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time to generate graph: = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    return g;
}

Graph getGraph(int numpoints,uint dimension) {
    
    int num_edges = numpoints*(numpoints-1)/2;
    Graph g(numpoints, num_edges);

    random_device rd; 
    mt19937 gen(rd()); 
    uniform_real_distribution<> dis(0,1.0);
    
    // Add all n choose 2 edges
    for (int i=0; i<numpoints; ++i) {
        for (int j=i+1; j<numpoints; ++j) {
                double weight = 0; // Euclidean distance between arr[i] and arr[j]
                if (dimension==0) {
                    weight = (double) dis(gen);
                }
                else {
                    double** arr = new double*[numpoints];
                    for (int i = 0; i < numpoints; ++i) {
                        arr[i] = new double[dimension];
                        for (uint j = 0; j < dimension; ++j) {
                            arr[i][j] = (double) dis(gen);
                        }
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
                    // Free arr
                    // Free each sub-array
                    for(int i = 0; i < numpoints; ++i) {
                        delete[] arr[i];   
                    }
                    // Free the array of pointers
                    delete[] arr;
                }
                
                g.addEdge(i,j,weight);
        }
    }
    return g;
}

double kruskal_trials(int numpoints, uint numtrials, uint dimension) {
    double mst_weight_sum = 0;
    for (uint i = 0; i < numtrials; ++i) {
        // NOTE: throw away edges
        Graph g=getGraphThrowAwaySomeEdges(numpoints,dimension);
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        double mst_weight = g.kruskal();
        mst_weight_sum+=mst_weight;
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        //delete &g;
        std::cout << numpoints << " points, dimension " << dimension << ", weight " << mst_weight << ", Time difference: = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    
    }
    cout << "\nMST Weight: " << mst_weight_sum/numtrials << endl;
    return mst_weight_sum/numtrials;
}

void produce_table_times(uint dimension) {
    ofstream outdata; // outdata is like cin
    int num_n_values_to_test = 12;
    //int n_values[num_n_values_to_test] = {128, 256, 512, 1024, 2048, 4096, 8192, 16384};
    int n_values[12] = {128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144};

    outdata.open("runtimes_dim"+to_string(dimension)+".dat", std::ios_base::app); // opens the file
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

    if (flag==1) { // run kruskal for one graph with the given numpoints, dimension
        Graph g=getGraphThrowAwaySomeEdges(numpoints,dimension);
        //Graph g=getGraph(numpoints,dimension);
        cout << g.kruskal() << endl;
    }
    if (flag == 2) { // generate all tables of values
        produce_table_times(0);
        produce_table_times(2);
        produce_table_times(3);
        produce_table_times(4);
    }
    if (flag == 0) { // generate just one table of values for the first argument
        produce_table_times(dimension);
    }
    if (flag == 4) {
        Graph g_=getGraphThrowAwaySomeEdges(numpoints,dimension);
        cout << "Throw away edges:" << g_.kruskal() << endl;
        Graph g=getGraph(numpoints,dimension);
        cout << "Regular: " << g.kruskal() << endl;
    }
    
    return 0;

    
}
