#include <iostream>
#include <cstdlib> 
#include <queue>
#include <vector>
#include <climits>
#include <ctime>
#include <fstream>
#include <cstring>
#include <chrono>

#include <stdio.h>
#include <omp.h>

#define NIL 0
#define INF INT_MAX

using namespace std;

class BGraph
{
    int m, n;

    std::vector<int> *adj;

public:
	int *pair_u, *pair_v, *dist;
	
	double solution = 0;
	
    BGraph(int m, int n);
	
    void addEdge(int u, int v);

	void addEdge(int u, int v, double val);
	
    bool bfs();

    bool dfs(int u);

    int hopcroftKarpAlgorithm();
	
	//std::vector<double> *values;
};


int BGraph::hopcroftKarpAlgorithm()
{
    pair_u = new int[m + 1];
    pair_v = new int[n + 1];
    dist = new int[m + 1];

    for (int u = 0; u <= m; u++)
        pair_u[u] = NIL;
    for (int v = 0; v <= n; v++)
        pair_v[v] = NIL;

    int result = 0;

    while (bfs())
    {
        for (int u = 1; u <= m; u++)
		{
            if (pair_u[u] == NIL && dfs(u))
			{
				result++;
				//solution += values[u].front();
			}	
			/*
			else
			{
				values[u].pop_front();
			}
			*/	
		}
    }
    return result;
}

bool BGraph::bfs()
{
    std::queue<int> q;

    for (int u = 1; u <= m; u++)
    {
        if (pair_u[u] == NIL)
        {
            dist[u] = 0;
            q.push(u);
        }

        else
            dist[u] = INF;
    }

    dist[NIL] = INF;

	#pragma omp parallel for
	for(int i = q.size(); i>0; i--)
    //while (!q.empty())
    {
        int u = q.front();
        q.pop();

        if (dist[u] < dist[NIL])
        {
            std::vector<int>::iterator it;
            for (it = adj[u].begin(); it != adj[u].end(); ++it)
            {
                int v = *it;

                if (dist[pair_v[v]] == INF)
                {
                    dist[pair_v[v]] = dist[u] + 1;
                    q.push(pair_v[v]);
                }
            }
        }
    }

    return (dist[NIL] != INF);
}

bool BGraph::dfs(int u)
{
    if (u != NIL)
    {
        std::vector<int>::iterator it;
        for (it = adj[u].begin(); it != adj[u].end(); ++it)
        {
            int v = *it;

            if (dist[pair_v[v]] == dist[u] + 1)
            {
                if (dfs(pair_v[v]) == true)
                {
                    pair_v[v] = u;
                    pair_u[u] = v;
                    return true;
                }
            }
        }

        dist[u] = INF;
        return false;
    }
    return true;
}

BGraph::BGraph(int m, int n)
{
    this->m = m;
    this->n = n;
    adj = new std::vector<int>[m + 1];
	//values = new std::vector<double>[m+1];
}

void BGraph::addEdge(int u, int v)
{
    adj[u].push_back(v);
	//values[u].push_back(0);
}

void BGraph::addEdge(int u, int v, double val)
{
    adj[u].push_back(v);
	//values[u].push_back(val);
}


int main(int argc, char* argv[])
{
	//omp_set_num_threads(4);
	
	string filename;
	if (argc > 1)
	{
		filename = argv[1];
	}
	else
	{
	/////////////////////////////////////	
		filename = "guys/3534_cage9.mtx";
	/////////////////////////////////////
	}
	ifstream file(filename);
	
	int n, m, max;
	int u, v;
	double val;

	while(file.peek() == '%') file.ignore(2048, '\n');
	
	file >> n >> m >> max;	
	//printf("%d %d %d\n", n, m, max);

	BGraph g(n, n);
	
	for (int l = 0; l < max; l++)
	{	
		file >> u >> v >> val;
		//printf("%d %d %f\n", u, v, val);
		g.addEdge(u, v, val);
	}

/*
	int val;
	n = 3534;
	BGraph g(n, n);
    FILE* fp;
    fp = fopen("guys/3534.txt", "rt");
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n;j++){
            fscanf(fp, "%d ", &val);
			g.addEdge(i, j, val);
        }
    }
	
	fclose(fp);
*/

	auto start = std::chrono::steady_clock::now();
	
	int res = g.hopcroftKarpAlgorithm();
	std::cout << "\nmax match = " << res;

	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	 
	std::cout << "\nelapsed time: " << elapsed_seconds.count() << "s\n";
	
	//std::cout << "\nsolution: " << g.solution;
	
	//file.close();
		
	//system("pause");
	return 0;
	
}

/*int conslolniymain()
{
    int v, u, e;
	std::cout << "u v num of edges\n";
    std::cin >> v >> u >> e; 
    BGraph g(v, u);
	std::cout << "pairs:\n";
    for (int i = 0; i < e; ++i)
    {
        std::cin >> u >> v;
        g.addEdge(u, v);
    }
  
    int res = g.hopcroftKarpAlgorithm();
    std::cout << "Maximum matching is " << res <<"\n";

    return 0;
}
*/
