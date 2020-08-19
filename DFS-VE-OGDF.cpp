#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <ogdf/basic/Graph.h>

using namespace ogdf;

Graph G;
NodeArray<List<node>> Adj;

int n, m; //number of vertices and edges

//DFS data
NodeArray<node> l;
NodeArray<node> low;
NodeArray<node> high;
NodeArray<node> M;
NodeArray<node> Mp;
NodeArray<int> up;
NodeArray<int> down; 

//lists of inverse functions 
NodeArray<List<node>> invM;
NodeArray<List<node>> invMp;
NodeArray<List<node>> invHigh; 

//DFS numbering, its inverse function, and parent
NodeArray<int> dfs;
node *idfs;
NodeArray<node> p; 

node r; //the root of DFS tree

//DFS Tree
NodeArray<List<node>> T;
node *BackEdges;

// union-find structure
NodeArray<node> ufparent;
NodeArray<node> representative;
NodeArray<int> ufrank; 


//performs a DFS from "root"; it's better for the graph to be connected
//computes: dfs, idfs, p, l, low
//constructs: DFS Tree; the lists of children are sorted in increasing order
void DFS(node root)
{
    r=root;
    T.init(G);
    BackEdges = new node[2*(m-n+1)];
    int bp=0; //back-edge pointer  

    dfs.init(G); 
    idfs = new node[n]; 
    p.init(G); 
    p[r]=nullptr;
    for(node v=G.firstNode(); v!=nullptr; v=v->succ()){dfs[v]=-1;} //initialize

    l.init(G); low.init(G); 
    for(node v=G.firstNode(); v!=nullptr; v=v->succ())
    {
       l[v]=v; low[v]=v;
    }
 
    node* stack = new node[n];
    ListIterator<node>* temp_out = new ListIterator<node>[n];

    int Nr=0; //current DFS label
    dfs[r]=Nr; idfs[Nr++]=r;
    stack[0]=r; temp_out[0]=Adj[r].begin();
    int SP=0;

    up.init(G); down.init(G);
    for(node v=G.firstNode(); v!=nullptr; v=v->succ()){up[v]=0; down[v]=0;}

    NodeArray<node> temp_child;
    temp_child.init(G);

    //perform DFS
    while(SP!=-1)
    {
       node v=stack[SP]; //current vertex
       char descend=0;
       for(ListIterator<node> i=temp_out[SP]; i!=Adj[v].end(); i++)
       {
          node u = *i; //the current neighbor of v
          if(dfs[u]==-1)
          {
             dfs[u]=Nr; idfs[Nr++]=u; p[u]=v; temp_child[v]=u; T[v].pushBack(u);
             stack[SP+1]=u; temp_out[SP+1]=Adj[u].begin(); temp_out[SP]=i;
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

    //delete[] TnextChild; delete[] stack; delete[] temp_out;
}

//sorts the list of back-edges in increasing order w.r.t. the (dfs of the) lower end of its elements
void sort_BackEdgesLow()
{
    int* indx = new int[n+1];
    for(int i=0;i<=n;i++){indx[i]=0;}

    int N=m-n+1; //number of back-edges
    node* copyBE = new node[2*N];
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
       node u=copyBE[2*i]; 
       node v=copyBE[2*i+1];
       BackEdges[indx[dfs[v]]++]=u;
       BackEdges[indx[dfs[v]]++]=v;
    }

    //delete[] indx; delete[] copyBE;
}

// path compression
void ufcompress(node v) {
    node p;
    if ((p = ufparent[v]) != v) {
        ufcompress(p);
        ufparent[v] = ufparent[p];
    }
}

node find(node v) {
    ufcompress(v);
    return ufparent[v];
}

// note: makes w the representative element
void unite(node v, node u, node w) {
    node p = find(v);
    node q = find(u);
    if ( p == q) return;
    if ( ufrank[p] > ufrank[q] )
    {
        // swap p and q
        node t = p;
        p = q;
        q = t;
    }
    else if ( ufrank[p] == ufrank[q] ) ufrank[q]++;

    ufparent[p] = q;
    representative[q] = w;
}

void ufinit() {
    ufparent.init(G); representative.init(G); ufrank.init(G);
    for (node v=G.firstNode(); v!=nullptr; v=v->succ())
    {
        ufparent[v] = representative[v] = v;
        ufrank[v] = 0;
    }
}


//computes high (=highp) and invHigh
void ComputeHigh()
{
    sort_BackEdgesLow(); //sort the back-edges in increasing order w.r.t. their lower end

    //initialize the dsu structure
    ufinit();

    int N=m-n+1;
    high.init(G);
    high[idfs[0]]=nullptr;
    high[idfs[1]]=nullptr;
    //process all back-edges in decreasing order
    for (int i = N - 1; i >= 0; i--) 
    {
       node u = BackEdges[2 * i];
       node v = BackEdges[2 * i + 1];
       // (u,v) is a back-edge where v in an ancestor of u
       node x = representative[find(u)];
       while (p[x]!=v) 
       {
          high[x] = v;    
          node next = representative[find(p[x])];
          unite(x, p[x], next);
          x = next;
       }
    }

    //calculate the inverse lists invHigh, and have their elements sorted in increasing order
    invHigh.init(G);
    for(int i=2;i<n;i++)
    {
       node v=idfs[i];
       invHigh[high[v]].pushBack(v);
    }
    //delete[] ufparent; delete[] representative; delete[] ufrank; delete[] invHighnext;
}


//computes M, Mp, invM, and invMp
void ComputeM()
{
    M.init(G);
    Mp.init(G);
    for(node v=G.firstNode(); v!=nullptr; v=v->succ()){M[v]=nullptr;Mp[v]=nullptr;}
   
    NodeArray<ListIterator<node>> L;
    NodeArray<ListIterator<node>> R;
    L.init(G); R.init(G);
    NodeArray<char> setL;
    setL.init(G);
    for(node v=G.firstNode(); v!=nullptr; v=v->succ()){setL[v]=0;}

    //begin computation, in a bottom-up fashion
    //for M, ignore the root; for Mp, ignore the child of the root as well
    for(int i=n-1;i>0;i--)
    {
       node v=idfs[i]; //i = dfs[v]
      
       //fix the pointers L and R
       for(ListIterator<node> cIndx=T[v].begin(); cIndx!=T[v].end(); cIndx++)
       {
          node c = *cIndx; //c is a child of v
          if(dfs[low[c]]<i)
          {
             if(!setL[v]){L[v]=cIndx; setL[v]=1;}
             R[v]=cIndx;
          }
       }

       //find M[v]
       if(l[v]!=v){M[v]=v;}
       else if(L[v]!=R[v]){M[v]=v;}
       else
       {
          node d=M[*L[v]];
          while(1)
          {
             if(dfs[l[d]]<i){M[v]=d;break;}
             while(dfs[low[*L[d]]]>=i){L[d]++;}
             while(dfs[low[*R[d]]]>=i){R[d]--;}
             if(L[d]!=R[d]){M[v]=d;break;}
             d=M[*L[d]];
          }
       }
      
       if(i==1){break;} //v is the child of the root
      
       int dfsP = dfs[p[v]];

       //find Mp[v]
       if(dfs[l[v]]<dfsP){Mp[v]=v;continue;}
       while(dfs[low[*L[v]]]>=dfsP){L[v]++;}        
       while(dfs[low[*R[v]]]>=dfsP){R[v]--;}
       if(L[v]!=R[v]){Mp[v]=v;continue;}
       else
       {
          node d=Mp[*L[v]];
          while(1)
          {
             if(dfs[l[d]]<dfsP){Mp[v]=d;break;}
             while(dfs[low[*L[d]]]>=dfsP){L[d]++;}
             while(dfs[low[*R[d]]]>=dfsP){R[d]--;}
             if(L[d]!=R[d]){Mp[v]=d;break;}
             d=Mp[*L[d]];
          }
       }
    }

    //calculate the inverse lists, and have their elements sorted in decreasing order
    invM.init(G);
    for(int i=n-1;i>0;i--)
    {
       node v=idfs[i];
       invM[M[v]].pushBack(v);
    }   

    invMp.init(G);
    for(int i=n-1;i>1;i--)
    {
       node v=idfs[i];
       invMp[Mp[v]].pushBack(v);
    } 

    //delete[] L; delete[] R; delete[] invMnext;
}


//sorts every list of children in decreasing order w.r.t. the highp of its elements
void sort_ChildrenHigh()
{
    //the inverse lists invHigh must have been computed, and sorted in increasing order
    
    for(node v=G.firstNode(); v!=nullptr; v=v->succ())
    {
       T[v].clear();
    }

    for(int i=n-1; i>-1; i--)
    {
       node x=idfs[i];
       List<node> temp;
       for(ListIterator<node> c = invHigh[x].begin(); c.valid(); c++)
       {
          temp.pushFront(*c);
       }           
       for(ListIterator<node> c = temp.begin(); c.valid(); c++)
       {
          T[p[*c]].pushBack(*c);
       }      

    }
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
                           
    NodeArray<int> b_edges;
    b_edges.init(G);
    for(node v=G.firstNode(); v!=nullptr; v=v->succ()){b_edges[v]=0;}

    //calculate b_edges[v];
    //ignore the root and the child of the root
    for(int i=n-1;i>1;i--)
    {
       node v=idfs[i];
       b_edges[v]=up[v];
       for(ListIterator<node> j=T[v].begin(); j!=T[v].end(); j++)
       {
          b_edges[v]+=b_edges[*j];
       }       
       b_edges[v]-=down[v];

       //check whether there is only one back-edge from T(v) to T(p(v),r]
       if(b_edges[v]==1){count[p[v]->index()]++;}
    }
   

    ///     
    ///vertex-edge cut-pairs with tree-edges        

    ///Get all the data we need              
    ComputeM();
    ComputeHigh();

    ///Algorithm: M(u)>v
    for(node x=G.firstNode(); x!=nullptr; x=x->succ())
    {
       ListIterator<node> c = invMp[x].begin();
       ListIterator<node> u = invM[x].begin();
       while(c.valid() && u.valid())
       {
          while(u.valid() && dfs[*u]>=dfs[p[*c]]){u++;}
          if(!u.valid()){break;}
          if(dfs[high[*c]]<dfs[*u])
          {
             int n_edges=0;
             ListIterator<node> first = u;
             ListIterator<node> last;
             while(u.valid() && dfs[high[*c]]<dfs[*u])
             {
                n_edges++;
                last=u;
                u++;
             }
             count[p[*c]->index()]+=n_edges;
             c++;
             while(c.valid() && dfs[*last]<dfs[p[*c]])
             {
                while(dfs[*first]>=dfs[p[*c]])
                {
                   n_edges--;
                   first++;
                }
                count[p[*c]->index()]+=n_edges;
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
    for(node v=G.firstNode(); v!=nullptr; v=v->succ())
    {
       ListIterator<node> u = invHigh[v].begin();
       ListIterator<node> c = T[v].begin();
       while(u.valid())
       {    
          if(down[*u]!=0){u++;continue;} //check if high(u)=highp(u)
          while(!(*c==T[v].back() || dfs[*(c.succ())]>dfs[*u]))
          {
             c++;
          }
          if(low[*u]==v || dfs[*u]<=dfs[Mp[*c]] )
          {
             count[v->index()]++;
          }
          u++;
       } 
    }

//for(node v=G.firstNode(); v!=nullptr; v=v->succ()){printf("%d: %d\n",v->index(),count[v->index()]);}

    //Algorithm:high(u)<v
    for(node x=G.firstNode(); x!=nullptr; x=x->succ())
    {
       ListIterator<node> u = invM[x].begin();
       ListIterator<node> c = invMp[x].begin();
       while(u.valid() && c.valid())
       {
          while(c.valid() && dfs[*c]>=dfs[*u]){c++;}
          if(!c.valid()){break;}
          if(down[*u]==0 && dfs[high[*u]]<dfs[p[*c]])
          {
             int n_edges=0;
             int h=dfs[high[*u]];
             while(c.valid() && h<dfs[p[*c]])
             {
                while(u.valid() && dfs[*c]<dfs[*u])
                {
                   n_edges++;
                   u++;
                }
                count[p[*c]->index()]+=n_edges;
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

    for(node v=G.firstNode(); v!=nullptr; v=v->succ())
    {
       if(invM[v].size()==0){continue;}
       ListIterator<node> u = invM[v].begin();
       u++;
       ListIterator<node> c = T[v].begin();
       node min=v;
       while(u.valid() && c.valid())
       {
          min=high[*c];
          while(u.valid() && dfs[*u]>dfs[min])
          {
             count[v->index()]++;
             u++;
          }
          min=low[*c];
          c++;
          while(c.valid() && dfs[high[*c]]>=dfs[min])
          {
             if(dfs[low[*c]]<dfs[min]) 
             {
                min=low[*c]; 
             }
             c++;
          }
          while(u.valid() && dfs[*u]>dfs[min]){u++;}
       }
       while(u.valid())
       {
          if(dfs[*u]<=dfs[min])
          {
             count[v->index()]++;
          }
          u++;
       }
    }
   // delete[] b_edges; delete[] TlastChild; delete[] is_valid;
}

int main(int argc, char *argv[])
{
    if (argc != 3) 
    {
       printf("\n usage: %s <input file> <output file>\n\n", argv[0]);
       exit(-1);
    }
    char *file = argv[1];

   /* read graph from file */
   FILE* fp = fopen(argv[1],"r");
   fscanf(fp,"%d %d",&n,&m);
   node* vertex = new node[n];
   for(int i=0;i<n;i++){vertex[i]=G.newNode();} //create vertices
   for(int i=0;i<m;i++)
   {
      int x,y;
      fscanf(fp,"%d %d",&x,&y);
      G.newEdge(vertex[x-1],vertex[y-1]); //create edges
   }
   fclose(fp);

     using namespace std::chrono;
     high_resolution_clock::time_point t1 = high_resolution_clock::now();

    //construct adjacency list
    Adj.init(G);
    for(edge e = G.firstEdge(); e!=nullptr; e=e->succ())
    {
       node x = e->source();
       node y = e->target();
       Adj[x].pushBack(y);
       Adj[y].pushBack(x);
    }

    r=G.firstNode();
    DFS(r); //perform a dfs from r

    compute_count();

     high_resolution_clock::time_point t2 = high_resolution_clock::now();
     duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

   fp = fopen(argv[2],"w");
   for(int i=0;i<n;i++){fprintf(fp,"%d\n",count[i]);}
   fprintf(fp,"%f\n",time_span.count());
   fclose(fp);

   printf("%f\n",time_span.count());
   return 0;	
}



