//  LSH
//  Created by Bharath Kumar Kande on 11/6/17.
//  Team of 3. Developed together by Bharath Kumar Kande, Manisha Siddartha Nalla and Sunandha Perumalla
//

double* read_file(char *file_path);
void read_binary_file(char *file_path, double *data);
int LSH(int dimensions, int ndata, double *data, int m, double **r, double *b, double w, int clusters, int **H, int *cluster_assign);
int search(int dimensions, int ndata, double *data, int m, double **r, double *b, double w, int clusters, int *cluster_size, int *cluster_start, double *query, int **H);

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define INSTANCES 5000
#define ATTRIBUTES 17

int main(int argc, const char * argv[]) {
    int dimensions, noofpoints, clusters=0, ndata, *cluster_assign, i, j, m, **H;
    double *data, **r, *b, w, sum=0.0;
    
    ndata=INSTANCES;
    dimensions=ATTRIBUTES;
    m=3;
    
    noofpoints=ndata*dimensions;
    
    data=(double *)malloc(noofpoints*sizeof(double));
    cluster_assign=(int *)malloc(ndata*sizeof(int));
    b=(double *)malloc(m*sizeof(double));
    r=(double **)malloc(m*sizeof(double *));
    for(i=0;i<m;i++)
    *(r+i)=(double *)malloc(dimensions*sizeof(double));
    H=(int **)malloc(ndata*sizeof(int *));
    for(i=0;i<ndata;i++)
    *(H+i)=(int *)malloc(m*sizeof(int));
    
    // Importing binary data file
    read_binary_file("/Users/bharathkumarkande/Desktop/bin files/ethylene_CO_bin.bin", data);
    
    // r vectors generation and normalization
    for(i=0;i<m;i++) {
        sum=0.0;
        for(j=0;j<dimensions;j++) {
            r[i][j]=((double)rand()*(101))/(double)RAND_MAX;
            sum+=r[i][j];
        }
        sum=sqrt(sum);
        for(j=0;j<dimensions;j++)
        r[i][j]=r[i][j]/sum;
    }
    
    // b vector(Can be all zeros)
    for(i=0;i<m;i++) {
        b[i]=0.0;
    }
    
    // Initialize w
    w=9000.0;
    
    LSH(dimensions, ndata, data, m, r, b, w, clusters, H, cluster_assign);
    return 0;
}

int LSH(int dimensions, int ndata, double *data, int m, double **r, double *b, double w, int clusters, int **H, int *cluster_assign) {
    int i,j,k,same=0,match=0,left=0,temp,clustersize=0, *cluster_size, *cluster_start;
    double sum=0.0, *query,temp1;
    query=(double *)malloc(dimensions*sizeof(double));
    
    // Find H[ndata][m] and cluster assign
    for(i=0;i<ndata;i++) {
        for(j=0;j<m;j++) {
            sum=0.0;
            for(k=0;k<dimensions;k++)
            sum+=data[(i*dimensions)+k]*r[j][k];
            sum=(sum-b[j])/w;
            sum=floor(sum);
            H[i][j]=(int)sum;
        }
        if(i==0) {
            cluster_assign[i]=0;
            clusters++;
        } else {
            match=0;
            for(j=0;j<i;j++) {
                same=0;
                for(k=0;k<m;k++) {
                    if(H[i][k]==H[j][k]) {
                        same++;
                    } else {
                        break;
                    }
                }
                match++;
                if(same==m) {
                    cluster_assign[i]=cluster_assign[j];
                    break;
                }
            }
            if(match==i) {
                cluster_assign[i]=clusters;
                clusters++;
            }
        }
    }
    // End: Find H[ndata][m] and cluster assign
    
    // Allocate memory for cluster start and cluster size
    cluster_size=(int *)malloc(sizeof(int)*clusters);
    cluster_start=(int *)malloc(sizeof(int)*clusters);
    
    left=0;
    // Sorting Cluster assign and find cluster sizes
    for(i=0;i<clusters;i++) {
        clustersize=0;
        for(j=left;j<ndata;j++) {
            if(cluster_assign[j]==i) {
                clustersize++;
                temp=cluster_assign[j];
                cluster_assign[j]=cluster_assign[left];
                cluster_assign[left]=temp;
                for(k=0;k<dimensions;k++) {
                    temp1=data[(j*dimensions)+k];
                    data[(j*dimensions)+k]=data[(left*dimensions)+k];
                    data[(left*dimensions)+k]=temp1;
                }
                for(k=0;k<m;k++) {
                    temp=H[j][k];
                    H[j][k]=H[left][k];
                    H[left][k]=temp;
                }
                left++;
            }
        }
        cluster_size[i]=clustersize;
    }
    // End: Sorting Cluster assign and find cluster sizes
    
    // Find cluster starts
    cluster_start[0]=0;
    for(i=1;i<clusters;i++)
    cluster_start[i]=cluster_start[i-1]+cluster_size[i-1];
    // End: Find cluster starts
    
    printf("Clusters: %d\n", clusters);
    
    // Search for 10 points for the dataset given
    double *query_total;
    query_total = (double*) malloc((dimensions*10) * sizeof(double));
    query_total = read_file("/Users/bharathkumarkande/Desktop/bin files/ethy_queries.txt");
    for(i=0;i<10;i++){
        for(j=0;j<dimensions;j++)
            query[j] = query_total[(i*dimensions)+j];
        search(dimensions, ndata, data, m, r, b, w, clusters, cluster_size, cluster_start, query, H);
    }
    // End: Search for 10 points for the dataset given

    return 0;
}


int search(int dimensions, int ndata, double *data, int m, double **r, double *b, double w, int clusters, int *cluster_size, int *cluster_start, double *query, int **H) {
    int j,k,match=0,closestcluster=-1, indexofNN=-1, *hofquery;
    double sum=0.0, min=-1.0;
    hofquery=(int *)malloc(sizeof(int)*m);
    
    // Hashing the query point
    for(j=0;j<m;j++) {
        sum=0.0;
        for(k=0;k<dimensions;k++)
        sum+=query[k]*r[j][k];
        sum=(sum-b[j])/w;
        sum=floor(sum);
        hofquery[j]=(int)sum;
    }
    // End: Hashing the query point
    
    
    // Find cluster, hash match
    for(j=0;j<clusters;j++) {
        match=0;
        for(k=0;k<m;k++) {
            if(hofquery[k]==H[cluster_start[j]][k]) {
                match++;
            } else {
                break;
            }
        }
        if(match==m) {
            closestcluster=j;
            break;
        }
    }
    // End: Find cluster, hash match
    
    
    if(closestcluster!=-1) {
        // Exhaustive search in nearest bucket
        min=-1.0;
        for(j=cluster_start[closestcluster];j<(cluster_start[closestcluster]+cluster_size[closestcluster]);j++) {
            sum=0.0;
            for(k=0;k<dimensions;k++) {
                sum+=pow((data[(j*dimensions)+k]-query[k]), 2);
            }
            sum=sqrt(sum);
            if((sum < min) || (min == -1.0)) {
                min=sum;
                indexofNN=j;
            }
        }
        // End: Exhaustive search in nearest bucket
        
        printf("Visits: %d\n", cluster_size[closestcluster]);
        printf("Distance to Nearest Neighbor: %lf\n\n", min);
        
    } else {
        printf("No matching Hash Values, adjust value w.\n");
    }
    
    return 0;
}


void read_binary_file(char *file_path, double *data) {
    FILE *file = fopen(file_path, "rb");
    int i;
    double temp_data = -1.0;
    for(i=0;i<(INSTANCES*ATTRIBUTES);i++){
        fread(&temp_data, sizeof(double), 1, file);
        data[i] = temp_data;
    }
    fclose(file);
}

double* read_file(char* file_path){
    FILE *file = fopen(file_path, "r");
    int i;
    double *data = (double*) malloc(10*ATTRIBUTES * sizeof(double));
    char *split_data = malloc(256 * sizeof(char));
    for(i=0;i<10*ATTRIBUTES;i++){
        fscanf(file, "%s", split_data);
        data[i] = strtod(split_data, NULL);
    }
    fclose(file);
    return data;
}


