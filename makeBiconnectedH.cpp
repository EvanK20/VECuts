#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm>

using namespace std;

vector <pair<int,int>> edgelist;
vector <int>* adj;
int n,m;
pair<int,int> edge;

char is_2v_connected()
{
    int* dfs = new int[n]; 
    int* idfs = new int[n]; 
    int* p = new int[n]; 
    p[0]=-1;
    for(int i=0;i<n;i++){dfs[i]=-1;} //initialize

    int* low = new int[n]; 
    for(int i=0;i<n;i++){low[i]=i;}
 
    int* stack = new int[n];
    vector<int>::iterator* temp_out = new vector<int>::iterator[n];

    int Nr=0; //current DFS label
    dfs[0]=Nr; idfs[Nr++]=0;
    stack[0]=0; temp_out[0]=adj[0].begin();
    int SP=0;

    //perform DFS from 0, and calculate low[v]
    while(SP!=-1)
    {
       int v=stack[SP]; //current vertex
       char down=0;
       vector<int>::iterator i;
       for(i=temp_out[SP]; i<adj[v].end(); i++)
       {
          int u = *i; //the current neighbor of v
          if(dfs[u]==-1)
          {
             dfs[u]=Nr; idfs[Nr++]=u; p[u]=v;
             stack[SP+1]=u; temp_out[SP+1]=adj[u].begin(); temp_out[SP]=i;
             down=1; break;
          }
          if(v==0){continue;}
          if(dfs[u]<dfs[v] && u!=p[v])
          {
             if(dfs[low[v]]>dfs[u]){low[v]=u;}
          }
          else if(dfs[u]>dfs[v] && p[u]==v)
          {
             if(dfs[low[v]]>dfs[low[u]]){low[v]=low[u];}
          }
       }
       if(down){SP++; continue;}
       SP--;
    }
   
    char nc=0; //number of children of the root
    for(int i=1;i<n;i++)
    {
       if(p[i]==0){nc++; if(nc==2){return 0;}}
       else
       {
          if(dfs[low[i]]>=dfs[p[i]]){return 0;}
       }
    }
    return 1;
}

char is_connected()
{
    char* found = new char[n];
    for(int i=0;i<n;i++){found[i]=0;}
    int* stack = new int[n];
    vector<int>::iterator* nextOut = new vector<int>::iterator[n];
    
    found[0]=1; int k=1; stack[0]=0; nextOut[0]=adj[0].begin();
    int SP=0;
    while(SP!=-1)
    {
       int v=stack[SP];
       char down=0;
       vector<int>::iterator i;
       for(i=nextOut[SP]; i<adj[v].end(); i++)
       {
          if(found[*i]){continue;}
          found[*i]=1; k++; 
          stack[SP+1]=*i; nextOut[SP+1]=adj[*i].begin(); nextOut[SP]=i+1;
          down=1; break;
       }
       if(down){SP++;continue;}
       SP--;
    }    

    return k==n;
}

void add_Hcycle()
{
    srand(time(NULL));
    vector<int> vertices;
    for (int i=0; i<n; i++) vertices.push_back(i);
    random_shuffle( vertices.begin(), vertices.end() );
    vertices.push_back(vertices[0]);
    vertices.push_back(-1);

    vector<int>::iterator i;
    for(i=vertices.begin(); *(i+1)!=-1; i++)
    {
       char e_exists=0;
       vector<int>::iterator j;
       for(j=adj[*i].begin(); j<adj[*i].end(); j++)
       {
          if(*j==*(i+1)){e_exists=1; break;}
       } 
       if(e_exists){continue;}
       edge.first=*i; edge.second=*(i+1);
       edgelist.push_back(edge);
       m++;
    }

    random_shuffle( edgelist.begin(), edgelist.end() );
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
       printf("\n usage: %s <input_file> <output_file>\n\n", argv[0]);
       exit(-1);
    }
    
    FILE* fp = fopen(argv[1],"r");
    fscanf(fp,"%d %d",&n,&m);
    adj = new vector<int>[n];
    for(int i=0;i<m;i++)
    {
       fscanf(fp,"%d %d",&edge.first,&edge.second);
       edge.first--; edge.second--;
       edgelist.push_back(edge);
       adj[edge.first].push_back(edge.second);
       adj[edge.second].push_back(edge.first);
    }
    fclose(fp);

    if(!is_connected())
    {
       add_Hcycle();
    }
    else if(!is_2v_connected())
    {
       add_Hcycle();
    }


    fp=fopen(argv[2],"w");
    fprintf(fp,"%d %d\n",n,m);
    for(int i=0;i<m;i++)
    {
       fprintf(fp,"%d %d\n",edgelist[i].first+1,edgelist[i].second+1);
    } 
    fclose(fp);

    return 0;
}
