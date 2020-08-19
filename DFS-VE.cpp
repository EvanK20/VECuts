#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <chrono>

int n, m; //number of vertices and edges
int *input_graph;
int *l, *low, *high, *M, *Mp, *up, *down; //DFS data
int *invM, *invMfirst, *invMp, *invMpfirst, *invHigh, *invHighfirst; //lists of inverse functions 
int *dfs, *idfs, *p; //DFS numbering, its inverse function, and parent
int r; //the root of DFS tree
int *T, *TfirstChild, *BackEdges; //DFS Tree
int *Gout, *GfirstOut; // adjacency list
int *ufparent, *representative, *ufrank; // union-find structure


//performs a DFS from "root"; it's better for the graph to be connected
//computes: dfs, idfs, p, l, low
//constructs: DFS Tree; the lists of children are sorted in increasing order
void DFS(int root)
{
    r=root;
    T = new int[n]; TfirstChild = new int[n+1];
    for(int i=0;i<=n;i++){TfirstChild[i]=0;}
    BackEdges = new int[2*(m-n+1)];
    int bp=0; //back-edge pointer  

    dfs = new int[n]; idfs = new int[n]; p = new int[n]; 
    p[0]=-1;
    for(int i=0;i<n;i++){dfs[i]=-1;} //initialize

    l = new int[n]; low = new int[n]; 
    for(int i=0;i<n;i++){l[i]=i;low[i]=i;}
 
    int* stack = new int[n];
    int* temp_out = new int[n];

    int Nr=0; //current DFS label
    dfs[r]=Nr; idfs[Nr++]=r;
    stack[0]=r; temp_out[0]=GfirstOut[r];
    int SP=0;

    up = new int[n]; down = new int[n];
    for(int i=0;i<n;i++){up[i]=0; down[i]=0;}

    int* temp_child = new int[n];
 
    //perform DFS
    while(SP!=-1)
    {
       int v=stack[SP]; //current vertex
       char descend=0;
       for(int i=temp_out[SP];i<GfirstOut[v+1];i++)
       {
          int u = Gout[i]; //the current neighbor of v
          if(dfs[u]==-1)
          {
             dfs[u]=Nr; idfs[Nr++]=u; p[u]=v; TfirstChild[v+1]++; temp_child[v]=u;
             stack[SP+1]=u; temp_out[SP+1]=GfirstOut[u]; temp_out[SP]=i;
             descend=1; break;
          }
          if(v==r){continue;}
          if(dfs[u]<dfs[v] && u!=p[v])
          {
             if(dfs[l[v]]>dfs[u]) 
             {
                l[v]=u;
                if(dfs[low[v]]>dfs[u]){low[v]=u;}
             }
             BackEdges[bp++]=v; BackEdges[bp++]=u;
             up[v]++; down[temp_child[u]]++;
          }
          else if(dfs[u]>dfs[v] && p[u]==v)
          {
             if(dfs[low[v]]>dfs[low[u]]){low[v]=low[u];}
          }
       }
       if(descend){SP++; continue;}
       SP--;
    }

    //construct DFS tree
    int* TnextChild = new int [n+1];
    for(int i=1;i<=n;i++)
    {
       TfirstChild[i]+=TfirstChild[i-1]; //the list of children of i starts from TfirstChild[i]
    }
    for(int i=0;i<=n;i++){TnextChild[i]=TfirstChild[i];}
    for(int i=1;i<n;i++)
    {
       int v = idfs[i]; 
       T[TnextChild[p[v]]++]=v;
    }

    //delete[] TnextChild; delete[] stack; delete[] temp_out;
}


//sorts the list of back-edges in increasing order w.r.t. the (dfs of the) lower end of its elements
void sort_BackEdgesLow()
{
    int* indx = new int[n+1];
    for(int i=0;i<=n;i++){indx[i]=0;}

    int N=m-n+1; //number of back-edges
    int* copyBE = new int[2*N];
    for(int i=0;i<2*N;i++){copyBE[i]=BackEdges[i];}

    for(int i=0;i<N;i++)
    {
       indx[dfs[BackEdges[2*i+1]]+1]+=2; //count how many times every value apears (+1, to include the other end-point)   
    }
    for(int i=1;i<=n;i++)
    {
       indx[i]+=indx[i-1];
    }
    for(int i=0;i<N;i++)
    {
       int u=copyBE[2*i]; 
       int v=copyBE[2*i+1];
       BackEdges[indx[dfs[v]]++]=u;
       BackEdges[indx[dfs[v]]++]=v;
    }

    //delete[] indx; delete[] copyBE;
}


// path compression
void ufcompress(int v) {
    int p;
    if ((p = ufparent[v]) != v) {
        ufcompress(p);
        ufparent[v] = ufparent[p];
    }
}

int find(int v) {
    ufcompress(v);
    return ufparent[v];
}

// note: makes w the representative element
void unite(int v, int u, int w) {
    int p = find(v);
    int q = find(u);
    if ( p == q) return;
    if ( ufrank[p] > ufrank[q] )
    {
        // swap p and q
        int t = p;
        p = q;
        q = t;
    }
    else if ( ufrank[p] == ufrank[q] ) ufrank[q]++;

    ufparent[p] = q;
    representative[q] = w;
}

void ufinit(int N) {
    int i;
    for (i=0; i<=N; i++)
    {
        ufparent[i] = representative[i] = i;
        ufrank[i] = 0;
    }
}

//computes high (=highp) and invHigh
void ComputeHigh()
{
    sort_BackEdgesLow(); //sort the back-edges in increasing order w.r.t. their lower end

    //initialize the dsu structure
    ufparent = new int[n];
    representative = new int[n];  
    ufrank = new int[n];
    ufinit(n-1);

    int N=m-n+1;
    high = new int[n]; 
    //process all back-edges in decreasing order
    for (int i = N - 1; i >= 0; i--) 
    {
       int u = BackEdges[2 * i];
       int v = BackEdges[2 * i + 1];
       // (u,v) is a back-edge where v in an ancestor of u
       int x = representative[find(u)];
       while (p[x]!=v) 
       {
          high[x] = v;    
          int next = representative[find(p[x])];
          unite(x, p[x], next);
          x = next;
       }
    }

    //calculate the inverse lists invHigh, and have their elements sorted in increasing order
    invHigh = new int[n];
    invHighfirst = new int[n+1];
    for(int i=0;i<=n;i++){invHighfirst[i]=0;}
    //ignore the root and the child of the root
    for(int i=2;i<n;i++)
    {
       invHighfirst[high[idfs[i]]+1]++;
    } 
    for(int i=1;i<=n;i++){invHighfirst[i]+=invHighfirst[i-1];}
    int* invHighnext = new int[n+1];
    for(int i=0;i<=n;i++){invHighnext[i]=invHighfirst[i];}
    for(int i=2;i<n;i++)
    {
       int v=idfs[i];
       invHigh[invHighnext[high[v]]++]=v;
    }
    //delete[] ufparent; delete[] representative; delete[] ufrank; delete[] invHighnext;
}

//computes M, Mp, invM, and invMp
void ComputeM()
{
    M = new int[n];
    Mp = new int[n];
    for(int i=0;i<n;i++){M[i]=-1;Mp[i]=-1;}
   
    int* L = new int[n];
    int* R = new int[n];
    for(int i=0;i<n;i++){L[i]=-1;}

    //begin computation, in a bottom-up fashion
    //for M, ignore the root; for Mp, ignore the child of the root as well
    for(int i=n-1;i>0;i--)
    {
       int v=idfs[i]; //i = dfs[v]
      
       //fix the pointers L and R
       for(int j=TfirstChild[v];j<TfirstChild[v+1];j++)
       {
          int c=T[j]; //c is a child of v
          if(dfs[low[c]]<i)
          {
             if(L[v]==-1){L[v]=j;}
             R[v]=j;
          }
       }

       //find M[v]
       if(l[v]!=v){M[v]=v;}
       else if(L[v]!=R[v]){M[v]=v;}
       else
       {
          int d=M[T[L[v]]];
          while(1)
          {
             if(dfs[l[d]]<i){M[v]=d;break;}
             while(dfs[low[T[L[d]]]]>=i){L[d]++;}
             while(dfs[low[T[R[d]]]]>=i){R[d]--;}
             if(L[d]!=R[d]){M[v]=d;break;}
             d=M[T[L[d]]];
          }
       }
      
       if(i==1){break;} //v is the child of the root
      
       int dfsP = dfs[p[v]];

       //find Mp[v]
       if(dfs[l[v]]<dfsP){Mp[v]=v;continue;}
       while(dfs[low[T[L[v]]]]>=dfsP){L[v]++;}        
       while(dfs[low[T[R[v]]]]>=dfsP){R[v]--;}
       if(L[v]!=R[v]){Mp[v]=v;continue;}
       else
       {
          int d=Mp[T[L[v]]];
          while(1)
          {
             if(dfs[l[d]]<dfsP){Mp[v]=d;break;}
             while(dfs[low[T[L[d]]]]>=dfsP){L[d]++;}
             while(dfs[low[T[R[d]]]]>=dfsP){R[d]--;}
             if(L[d]!=R[d]){Mp[v]=d;break;}
             d=Mp[T[L[d]]];
          }
       }
    }

    //calculate the inverse lists, and have their elements sorted in decreasing order
    invM = new int[n];
    invMfirst = new int[n+1];
    for(int i=0;i<=n;i++){invMfirst[i]=0;}
    for(int i=1;i<n;i++)
    {
       int v=idfs[i];
       invMfirst[M[v]+1]++;
    }
    for(int i=1;i<=n;i++){invMfirst[i]+=invMfirst[i-1];}
    int* invMnext = new int[n+1];
    for(int i=0;i<=n;i++){invMnext[i]=invMfirst[i];}
    for(int i=n-1;i>0;i--)
    {
       int v=idfs[i];
       invM[invMnext[M[v]]++]=v;
    }   

    invMp = new int[n];
    invMpfirst = new int[n+1];
    for(int i=0;i<=n;i++){invMpfirst[i]=0;}
    for(int i=2;i<n;i++)
    {
       int v=idfs[i];
       invMpfirst[Mp[v]+1]++;
    }
    for(int i=1;i<=n;i++){invMpfirst[i]+=invMpfirst[i-1];}
    for(int i=0;i<=n;i++){invMnext[i]=invMpfirst[i];}
    for(int i=n-1;i>1;i--)
    {
       int v=idfs[i];
       invMp[invMnext[Mp[v]]++]=v;
    } 

    //delete[] L; delete[] R; delete[] invMnext;
}

//sorts every list of children in decreasing order w.r.t. the highp of its elements
void sort_ChildrenHigh()
{
    //the inverse lists invHigh must have been computed, and sorted in increasing order
    
    int* TnextChild = new int[n+1];
    for(int i=0;i<=n;i++){TnextChild[i]=TfirstChild[i];}
    
    for(int i=n-1;i>=0;i--)
    {
       int v=idfs[i];
       for(int j=invHighfirst[v+1]-1;j>=invHighfirst[v];j--)
       {
          int x=invHigh[j];
          T[TnextChild[p[x]]++]=x;
       }  
    }    

/*for(int i=0;i<n;i++)
{
  printf("%d: ",i);
  for(int j=TfirstChild[i];j<TfirstChild[i+1];j++)
   {
     printf("%d(%d) ",T[j],high[T[j]]);
   }printf("\n");
}*/
    //delete[] TnextChild;
}

int* count;
//the main algorithm; computes count[v] in time O(m)
//input graph must be 2-vertex-connected
void compute_count()
{
    count = new int[n];
    for(int i=0;i<n;i++){count[i]=0;}
        
    ///vertex-edge cut-pairs with back-edges
                           
    int* b_edges = new int[n];
    for(int i=0;i<n;i++){b_edges[i]=0;}

    //calculate b_edges[v];
    //ignore the root and the child of the root
    for(int i=n-1;i>1;i--)
    {
       int v=idfs[i];
       b_edges[v]=up[v];
       for(int j=TfirstChild[v];j<TfirstChild[v+1];j++)
       {
          b_edges[v]+=b_edges[T[j]];
       }       
       b_edges[v]-=down[v];

       //check whether there is only one back-edge from T(v) to T(p(v),r]
       if(b_edges[v]==1){count[p[v]]++;}
    }
   

    ///     
    ///vertex-edge cut-pairs with tree-edges        

    ///Get all the data we need              
    ComputeM();
    ComputeHigh();

    //find all vertices u with high(u)=highp(u); they have down(u)=0
    //is_valid[u]=1 if and only if high[u]=highp[u]
    char* is_valid = new char[n];
    for(int i=0;i<n;i++){is_valid[i]=down[i]==0;}

    ///Algorithm: M(u)>v
    for(int x=0;x<n;x++)
    {
       int c = invMpfirst[x];
       int u = invMfirst[x];
       int endMp = invMpfirst[x+1];
       int endM = invMfirst[x+1];  
       while(c!=endMp && u!=endM)
       {
          while(u!=endM && dfs[invM[u]]>=dfs[p[invMp[c]]]){u++;}
          if(u==endM){break;}
          if(dfs[high[invMp[c]]]<dfs[invM[u]])
          {
             int n_edges=0;
             int first = u;
             while(u!=endM && dfs[high[invMp[c]]]<dfs[invM[u]])
             {
                n_edges++;
                u++;
             }
             int last = u-1;
             count[p[invMp[c]]]+=n_edges;
             c++;
             while(c!=endMp && dfs[invM[last]]<dfs[p[invMp[c]]])
             {
                while(dfs[invM[first]]>=dfs[p[invMp[c]]])
                {
                   n_edges--;
                   first++;
                }
                count[p[invMp[c]]]+=n_edges;
                c++;
             }
          }
          else
          {
             c++;
          }
       }
    }


    //Algorithm:high(u)=v
    for(int v=0;v<n;v++)
    {
       int u=invHighfirst[v];
       int c=TfirstChild[v];
       int endH=invHighfirst[v+1];
       int endT=TfirstChild[v+1];
       while(u!=endH)
       {    
          if(!is_valid[invHigh[u]]){u++;continue;} //check if high(u)=highp(u)
          while(!(c==endT-1 || dfs[T[c+1]]>dfs[invHigh[u]]))
          {
             c++;
          }
          if(low[invHigh[u]]==v || dfs[invHigh[u]]<=dfs[Mp[T[c]]])
          {
             count[v]++;
          }
          u++;
       } 
    }

    //Algorithm:high(u)<v
    for(int x=0;x<n;x++)
    {
       int u=invMfirst[x];
       int c=invMpfirst[x];
       int endM=invMfirst[x+1];
       int endMp=invMpfirst[x+1];
       while(u!=endM && c!=endMp)
       {
          while(c!=endMp && dfs[invMp[c]]>=dfs[invM[u]]){c++;}
          if(c==endMp){break;}
          if(is_valid[invM[u]] && dfs[high[invM[u]]]<dfs[p[invMp[c]]])
          {
             int n_edges=0;
             int h=dfs[high[invM[u]]];
             while(c!=endMp && h<dfs[p[invMp[c]]])
             {
                while(u!=endM && dfs[invMp[c]]<dfs[invM[u]])
                {
                   n_edges++;
                   u++;
                }
                count[p[invMp[c]]]+=n_edges;
                c++;
             }
          }
          else
          {
             u++;
          }
       }
    }

    //Algorithm: M(u)=v
    sort_ChildrenHigh(); //sort the lists of children in decreasing order w.r.t. the highp of their elements

    for(int v=0;v<n;v++)
    {
       if(invMfirst[v]==invMfirst[v+1]){continue;}
       int u=invMfirst[v]+1;
       int c=TfirstChild[v];
       int endM=invMfirst[v+1];
       int endT=TfirstChild[v+1];
       int min=v;
       while(u!=endM && c!=endT)
       {
          min=high[T[c]];
          while(u!=endM && dfs[invM[u]]>dfs[min])
          {
             count[v]++;
             u++;
          }
          min=low[T[c]];
          c++;
          while(c!=endT && dfs[high[T[c]]]>=dfs[min])
          {
             if(dfs[low[T[c]]]<dfs[min]) 
             {
                min=low[T[c]]; 
             }
             c++;
          }
          while(u!=endM && dfs[invM[u]]>dfs[min]){u++;}
       }
       while(u!=endM)
       {
          if(dfs[invM[u]]<=dfs[min])
          {
             count[v]++;
          }
          u++;
       }
    }

   // delete[] b_edges; delete[] TlastChild; delete[] is_valid;
}

//create the adjacency list structure from input_graph
void processInput()
{
    Gout = new int[2*m];
    GfirstOut = new int[n+1];
    int* GnextOut = new int[n+1];
    GnextOut[0]=0;
    for(int x=0;x<=n;x++){GfirstOut[x]=0;} //initialize
    
    for(int i=0;i<2*m;i++)
    {
       int x=input_graph[i];
       GfirstOut[x+1]++; //update the degree of x
    } 
    for(int x=1;x<=n;x++)
    {     
       GfirstOut[x] += GfirstOut[x-1]; //the list of neighbours of x starts from GfirstOut[x]
       GnextOut[x] = GfirstOut[x];
    }
    for(int i=0;i<2*m;i+=2)
    {
       int x=input_graph[i];
       int y=input_graph[i+1];
       Gout[GnextOut[x]++]=y;
       Gout[GnextOut[y]++]=x;
    }

    //delete[] GnextOut;
}

void readGraph(char* filename)
{
    FILE* fp = fopen(filename,"r");
    fscanf(fp,"%d %d",&n,&m);
    input_graph = new int[2*m];
    for(int i=0;i<m;i++)
    {
       fscanf(fp,"%d %d",&input_graph[2*i],&input_graph[2*i+1]);
       input_graph[2*i]--; input_graph[2*i+1]--; //vertex numbering starts from 1
    }
    fclose(fp);
}

int main(int argc, char *argv[])
{
    if (argc != 3) 
    {
       printf("\n usage: %s <input file> <output file>\n\n", argv[0]);
       exit(-1);
    }
    char *file = argv[1];

    readGraph(file); //extract input_graph

    processInput(); //create adjacency list

     using namespace std::chrono;
     high_resolution_clock::time_point t1 = high_resolution_clock::now();

    r=0;
    DFS(r); //perform a dfs from r

    compute_count();

     high_resolution_clock::time_point t2 = high_resolution_clock::now();
     duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

   FILE* fp = fopen(argv[2],"w");
   for(int i=0;i<n;i++){fprintf(fp,"%d\n",count[i]);}
   fprintf(fp,"%f\n",time_span.count());
   fclose(fp);

   printf("%f\n",time_span.count());
   return 0;	
}



