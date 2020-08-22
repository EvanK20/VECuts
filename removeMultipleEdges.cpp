#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include <algorithm>

using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
       printf("\n usage: %s <input_file> <output_file>\n\n", argv[0]);
       exit(-1);
    }
    
    int n,m;
    FILE* fp = fopen(argv[1],"r");
    fscanf(fp,"%d %d",&n,&m);
    int* edges = new int[2*m]; 
    for(int i=0;i<m;i++)
    {
       int x,y;
       fscanf(fp,"%d %d",&x,&y);
       if(x>y){int temp=x; x=y; y=temp;}
       edges[2*i]=x; edges[2*i+1]=y; 
    }
    fclose(fp);

    //sort w.r.t. the first end-point
    int* copyE = new int[2*m];
    for(int i=0;i<2*m;i++){copyE[i]=edges[i];}

    int* first = new int[n+2];
    for(int i=0;i<=n+1;i++){first[i]=0;}
    for(int i=0;i<m;i++){first[edges[2*i]+1]+=2;}
    for(int i=1;i<=n+1;i++){first[i]+=first[i-1];}
    int* next_indx = new int[n+2];
    for(int i=0;i<=n+1;i++){next_indx[i]=first[i];}    
    
    for(int i=0;i<m;i++)
    {
       int x=copyE[2*i];
       int y=copyE[2*i+1];
       edges[next_indx[x]++]=x;
       edges[next_indx[x]++]=y; 
    }


    //sort w.r.t. the second end-point
    for(int i=0;i<2*m;i++){copyE[i]=edges[i];}

    first = new int[n+2];
    for(int i=0;i<=n+1;i++){first[i]=0;}
    for(int i=0;i<m;i++){first[edges[2*i+1]+1]+=2;}
    for(int i=1;i<=n+1;i++){first[i]+=first[i-1];}
    next_indx = new int[n+2];
    for(int i=0;i<=n+1;i++){next_indx[i]=first[i];}    
    
    for(int i=0;i<m;i++)
    {
       int x=copyE[2*i];
       int y=copyE[2*i+1];
       edges[next_indx[y]++]=x;
       edges[next_indx[y]++]=y; 
    }
    

    //count the number of single edges
    int mNew=0;
    for(int i=0;i<m;i++)
    {
       mNew++;
       int x=edges[2*i];
       int y=edges[2*i+1];
       int indx=i+1;
       while(indx<m)
       {
          if(edges[2*indx]==x && edges[2*indx+1]==y)
          {
             indx++;
          }
          else
          {
             break;
          }
       }
       i+=indx-i-1;
    }

    //create edge list
    vector <pair<int,int> > edgelist;
    edgelist.reserve(mNew);
    for(int i=0;i<m;i++)
    {
       int x=edges[2*i];
       int y=edges[2*i+1];
       int indx=i+1;
       while(indx<m)
       {
          if(edges[2*indx]==x && edges[2*indx+1]==y)
          {
             indx++;
          }
          else
          {
             break;
          }
       }
       i+=indx-i-1;

       if(rand()%2==0){int temp=x; x=y; y=temp;}
       pair <int,int> e;
       e.first=x; e.second=y;
       edgelist.push_back(e);
    }
   
    srand(time(NULL));
    random_shuffle( edgelist.begin(), edgelist.end() );

    fp=fopen(argv[2],"w");
    fprintf(fp,"%d %d\n",n,mNew);
    vector <pair<int,int> >::iterator k;
    for (k=edgelist.begin(); k<edgelist.end(); k++)
    {
        int x = (*k).first;
        int y = (*k).second;
        fprintf(fp,"%d %d\n", x, y);
    }
    fclose(fp);
    return 0;
}





