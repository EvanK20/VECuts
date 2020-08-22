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

int* dfs; int* idfs; int* p; int* low;

char is_2v_connected()
{
    dfs = new int[n]; 
    idfs = new int[n]; 
    p = new int[n]; 
    p[0]=-1;
    for(int i=0;i<n;i++){dfs[i]=-1;} //initialize

    low = new int[n]; 
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

void make_2v_connected()
{
    //construct the dfs tree
    vector<int>* T = new vector<int>[n];
    for(int i=1;i<n;i++)
    {
       T[p[i]].push_back(i);
    }

    //find the articulation points
    char* is_cutvertex = new char[n];
    for(int i=0;i<n;i++){is_cutvertex[i]=0;}
    is_cutvertex[0]=T[0].size()>1;
    for(int i=1;i<n;i++)
    {
       if(p[i]==0){continue;}
       if(dfs[low[i]]>=dfs[p[i]]){is_cutvertex[p[i]]=1;}
    }

    //compute the 2-vertex-connected components in the dfs order
    vector<vector<int>> components;
    int count=0;
   
    int* stack = new int[n]; //stack for vertex
    int* stackC = new int[n]; //stack for component
    vector<int>::iterator* next_child = new vector<int>::iterator[n];

    if(!is_cutvertex[0])
    {
       vector<int> C0;
       C0.push_back(0);
       components.push_back(C0); count=1;
    }
    stack[0]=0; stackC[0]=0; next_child[0]=T[0].begin();
    int SP=0;
    while(SP!=-1)
    {
       int v=stack[SP];
       int k=stackC[SP];
       vector<int>::iterator u = next_child[SP];
       if(u==T[v].end())
       {
          SP--; continue;
       }
       next_child[SP]++; 
       stack[SP+1]=*u; next_child[SP+1]=T[*u].begin();
       if(v!=0 && dfs[low[*u]]>=dfs[v])
       {
          vector<int> newC;
          newC.push_back(v); newC.push_back(*u);
          components.push_back(newC);
          stackC[SP+1]=count++;
          SP++; continue;
       }
       else if(v==0 && is_cutvertex[0])
       {
          vector<int> newC;
          newC.push_back(v); newC.push_back(*u);
          components.push_back(newC);
          stackC[SP+1]=count++;
          SP++; continue;
       }
       components[k].push_back(*u);
       stackC[SP+1]=k;
       SP++;  
    }
    
    //print tree
    //for(int i=0;i<n;i++){printf("%d: ",i+1);for(int j=0;j<T[i].size();j++){printf("%d ",T[i][j]+1);}printf("\n");}

    //printf("cut vertices: ");for(int i=0;i<n;i++){if(is_cutvertex[i]){printf("%d ",i+1);}}printf("\n\n");
  
    //print components
    /*for(int i=0;i<components.size();i++)
    {
       for(int j=0;j<components[i].size();j++)
       {
          printf("%d ",components[i][j]+1);
       }printf("\n");
    }*/

    //find all components that are leaves in the block-graph
    //a component is a leaf if and only if it contains only one articulation point
    vector<int> leaves;
    int L=0;
    for(int i=0;i<count;i++)
    {
       int countAP=0;
       int size=components[i].size();
       for(int j=0;j<size;j++)
       {
          countAP+=is_cutvertex[components[i][j]];
       }
       if(countAP==1){leaves.push_back(i); L++;}
    }
    
    //print leaves
    //for(int i=0;i<L;i++){for(int j=0;j<components[leaves[i]].size();j++){printf("%d ",components[leaves[i]][j]+1);}printf("\n");}

    //connect the leaves in a path structure, by joining vertices that are not articulation points
    for(int i=0;i<L-1;i++)
    {
       int l1=leaves[i];
       int l2=leaves[(i+1)%L];
       int indx1 = rand()%components[l1].size();
       int indx2 = rand()%components[l2].size();
       //none of the end points must be an articulation point
       if(is_cutvertex[components[l1][indx1]]){indx1=(indx1+1)%components[l1].size();}
       if(is_cutvertex[components[l2][indx2]]){indx2=(indx2+1)%components[l2].size();}

       edge.first=components[l1][indx1];
       edge.second=components[l2][indx2];
       edgelist.push_back(edge);
       m++;
    }
    random_shuffle(edgelist.begin(), edgelist.end());
}

void make_connected()
{
    vector <vector<int>> components;
    int k=0; //number of components

    char* found = new char[n]; 
    for(int i=0;i<n;i++){found[i]=0;}
    int* stack = new int[n];
    vector<int>::iterator* next_out = new vector<int>::iterator[n];

    for(int i=0;i<n;i++)
    {
       if(found[i]){continue;}
       vector<int> C;
       C.push_back(i);
       found[i]=1;
       stack[0]=i; next_out[0]=adj[i].begin();
       int SP=0;
       while(SP!=-1)
       {
          int v=stack[SP];
          char down=0;
          vector<int>::iterator j;
          for(j=next_out[SP]; j<adj[v].end(); j++)
          {
             if(found[*j]){continue;}
             C.push_back(*j);
             found[*j]=1;
             stack[SP+1]=*j; next_out[SP+1]=adj[*j].begin(); next_out[SP]=j+1;
             down=1; break;
          }
          if(down){SP++;continue;}
          SP--;
       }
       components.push_back(C);
       k++;
    }

    //connect the components in a path structure
    for(int i=0;i<k-1;i++)
    {
       edge.first=components[i][rand()%components[i].size()];
       edge.second=components[i+1][rand()%components[i+1].size()];
       edgelist.push_back(edge);
       adj[edge.first].push_back(edge.second);
       adj[edge.second].push_back(edge.first);
       m++;       
    }
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
    }
    fclose(fp);

    srand(time(NULL));

   // random_shuffle(edgelist.begin(), edgelist.end());   
 
    for(int i=0;i<m;i++)
    {
       adj[edgelist[i].first].push_back(edgelist[i].second);
       adj[edgelist[i].second].push_back(edgelist[i].first);
    }

    if(!is_connected())
    { 
       make_connected();
    }
    if(!is_2v_connected())
    { 
       make_2v_connected();
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
