/* Create a random graph according to the G_np model */
/* It optionally adds a random Hamiltonian cycle */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <utility>
#include <algorithm>


using namespace std;


int main(int argc, char *argv[]) {
    
    if (argc != 5) {
        printf("\n usage: %s <nodes> <edges> <H> <seed>\n\n", argv[0]);
	exit(-1);
    }

    vector <pair<int,int> > edgelist;
    vector <int> vertices;
    

    int n = atoi(argv[1]);
    int m = atoi(argv[2]);
    int H = atoi(argv[3]);
    int seed = atoi(argv[4]);

    srand(seed);

    vector <int> *adj;
    adj = new vector<int> [n+1];

    pair <int,int> edge;

    int nedges = 0;


    double p = (double) (2*m)/(n*(n-1));

    for(int i=1; i<n; i++)
        for(int j=i+1; j<=n; j++)
            if (rand() < p*RAND_MAX)
            {
                edge.first = i;
                edge.second = j;
                edgelist.push_back(edge);
                adj[i].push_back(j);
                nedges++;   
            }



    if (H)
    {
        for (int i=1; i<=n; i++) vertices.push_back(i);
        random_shuffle( vertices.begin(), vertices.end() );
        vertices.push_back(vertices.front());
        vertices.push_back(0);

        vector<int>::iterator j;

        int x,y;
        for (j=vertices.begin(); *(j+1)!=0; j++)
        {
            char e_exists=0; //does edge (j,j+1) exist?
            vector<int>::iterator i;
            for(i=adj[*j].begin(); i<adj[*j].end(); i++)
            {
               if(*i==*(j+1)){e_exists=1;break;}
            }
            if(e_exists){continue;}
            for(i=adj[*(j+1)].begin(); i<adj[*(j+1)].end(); i++)
            {
               if(*i==*j){e_exists=1;break;}
            }
            if(e_exists){continue;}

            x = *j;
            y = *(j+1);
            edge.first = x;
            edge.second = y;
            edgelist.push_back(edge);
            nedges++;
        }
    }

    printf("nedges = %d, RAND_MAX=%d, p=%g\n", nedges, RAND_MAX, p);

    random_shuffle( edgelist.begin(), edgelist.end() );


    char fname[80];
    strcpy(fname, "randGr.");
    strcat(fname, argv[1]);
    strcat(fname, ".");
    strcat(fname, argv[2]);
    strcat(fname, ".");
    strcat(fname, argv[3]);
    strcat(fname, ".");
    strcat(fname, argv[4]);
    FILE *output = fopen (fname, "w");
    if (!output) {
        fprintf (stderr, "Error creating file \"%s\".\n", fname);
        exit(-1);
    }
    fprintf(output,"%d %d\n", n, nedges);
    
    vector <pair<int,int> >::iterator k;
    for (k=edgelist.begin(); k<edgelist.end(); k++)
    {
        int x = (*k).first;
        int y = (*k).second;
        fprintf(output,"%d %d\n", x, y);
    }

    fclose(output);
    return 0;
}
