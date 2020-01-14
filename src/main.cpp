#include <sys/time.h>
#include "divsufsort.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fm-tree.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unistd.h> /* getopt */
#include <vector>
#include "lib/bitmask.h"
#include "lib/bitmask_bitset.h"
#include "lib/data.h"
#include "lib/rank_select.h"
#include "lib/rb_tree.h"
#include <windows.h>
#include <psapi.h>
#define ARRAY_SIZE 100

using namespace std;

/*
	Creates and returns wavelet tree.
	@Author: Robert Jambrecic
*/
rank_select* create_wavelet_tree(const unsigned char *bwt, int n){

    string content((const char *)bwt);
    function<bitmask *(uint32_t)> c;
    unordered_map<char, uint32_t> counters;
    int size=(n/2)*2+2;
    uint32_t word_size = size;
    c = bitmask_bitset::create;
    bitmask::set_creator(c);
    vector<data*> data_chunks;

    if (content.length()) {
        data_chunks.push_back(data::create(content, word_size, counters));
    }

    return new rb_tree(data_chunks);

}


/*
	Creates and returns binary wavelet tree.
	@Author: Robert Jambrecic
*/
rank_select* create_binary_wavelet_tree(int *B, int n){

    string b_str;
    for(int i=0;i<n;i++){
        b_str += to_string(test_bit(B,i));
    }

    vector<data *> data_chunks_bin;
    function<bitmask *(uint32_t)> c_bin;
    unordered_map<char, uint32_t> counters_bin;
    int size=(n/2)*2+2;
    uint32_t word_size_bin = n;
    c_bin = bitmask_bitset::create;
    bitmask::set_creator(c_bin);

    if (b_str.length()) {
        data_chunks_bin.push_back(data::create(b_str, word_size_bin, counters_bin));
    }

    return new rb_tree(data_chunks_bin);
}


/*
	Writes set of locations R into result file.
	@Author: Andrea Bernat
*/
void make_result_file(set<int> R){

    ofstream outfile ("result.txt");

    outfile << "Number of locations: ";
    outfile << R.size();
    outfile << "\n\nLocations: \n";

    int i = 0;
    for (auto elem : R) {
        i++;
        outfile << elem;

        if(i != R.size()){
            outfile << ", ";
        }
    }
}

/*
	Returns the size of the file.
	@Author: Anel Hadzimuratagic
*/
long size_of_file(FILE *fp, char *fname){

    long n;

    /* Get the file size. */
    if(fseek(fp, 0, SEEK_END) == 0) {
        n = ftell(fp);
        rewind(fp);
        if(n < 0) {
            fprintf(stderr, "%s: Cannot ftell `%s': ", "main", fname);
            perror(NULL);
            exit(EXIT_FAILURE);
        }
    } else {
        fprintf(stderr, "%s: Cannot fseek `%s': ", "main", fname);
        perror(NULL);
        exit(EXIT_FAILURE);
    }

    return n;
}


/*
	Uppercase whole char array P.
	@Author: Andrea Bernat
*/
void upperCase(char *P){

    for (int i=0; i < ARRAY_SIZE; i++){
        P[i] = toupper(P[i]);
    }
}

/*
	Asks the user to provide the pattern file and reads the pattern from the file.
	@Author: Robert Jambrecic
*/
char* make_pattern(){
    FILE *fp;
    char fname[ARRAY_SIZE];
    long n;

    fprintf(stdout, "\n****************************************************************"
                    "\nConsider that if your file content differ from {A,C,G,T}"
                    "\nyou will get wrong results!"
                    "\n****************************************************************\n");
    fprintf(stdout, "\nInput name of pattern file: ");
    scanf("%s", fname);


    if((fp = fopen(fname, "rb")) == NULL) {
        printf(" Cannot open file ");
        perror(NULL);
        exit(EXIT_FAILURE);
    }

    n = size_of_file(fp, fname);

    char *P = (char *)malloc((size_t)n * sizeof(char));

    /* Read n bytes of data. */
    if(fread(P, sizeof(unsigned char), (size_t)n, fp) != (size_t)n) {
        fprintf(stderr, "%s: %s `%s': ",
                "main",
                (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
                fname);
        perror(NULL);
        exit(EXIT_FAILURE);
    }

    fclose(fp);

    return P;
}

/*
	Function for calculation of memory usage of program
	@Author: Anel Hadzimuratagic

*/
void calculate_memory_usage(){
    PROCESS_MEMORY_COUNTERS pmc;
    BOOL result = GetProcessMemoryInfo(GetCurrentProcess(),
                                       &pmc,
                                       sizeof( pmc ));
    if(result){
        printf( "\nThere is %d KB memory in use.\n", pmc.WorkingSetSize/1024 );
    }
}


/*
	Main function for FM tree execution.
	@Author: Anel Hadzimuratagic, Robert Jambrecic, Andrea Bernat
*/
int main() {
    FILE *fp;
    char fname[ARRAY_SIZE];
    long n;
    rank_select *rank;
    rank_select *rank_b;
    int sp=0;
    int ep=0;
    char *P;
    int D;
	
	fprintf(stdout, "\n****************************************************************"
                    "\nConsider that if your file content differ from {A,C,G,T,$}"
                    "\nyou will get wrong results!"
                    "\n****************************************************************\n");   
    fprintf(stdout, "\nInput name of text file: ");
    scanf("%s", fname);

    if((fp = fopen(fname, "rb")) == NULL) {
        printf(" Cannot open file ");
        perror(NULL);
        exit(EXIT_FAILURE);
    }

	fprintf(stdout, "Input sampling distance: ");
	scanf("%d", &D);

	while(D<1 || D>10){
		fprintf(stdout, "\n***Sampling distance must be in range [1,10]***\n");
		fprintf(stdout, "\nInput sampling distance: ");
	    scanf("%d", &D);
	};

    n = size_of_file(fp, fname);

    /* Allocate memory. */
    unsigned char *T = (unsigned char *)malloc((size_t)n * sizeof(unsigned char));
    unsigned char *bwt = (unsigned char *)malloc((size_t)n * sizeof(unsigned char));
    unsigned char *sortedT = (unsigned char *)malloc((size_t)n * sizeof(unsigned char));
    int *SA = (int *)malloc((size_t)n * sizeof(int));
    unsigned int *SSA = (unsigned int *)malloc((size_t)(int(ceil(n/(float)D))) * sizeof(int));
    int *C = (int *)malloc((size_t)5 * sizeof(int));
    int *B = (int *)malloc((size_t)((size_t)int(ceil(n/(float)32)) * sizeof(int)));

    if((T == NULL) ||(bwt == NULL) ||(sortedT == NULL) || (SA == NULL)
        || (SSA == NULL)|| (C == NULL) || (B == NULL)) {
        fprintf(stderr, "%s: Cannot allocate memory.\n", "main");
        exit(EXIT_FAILURE);
    }

    /* Read n bytes of data. */
    if(fread(T, sizeof(unsigned char), (size_t)n, fp) != (size_t)n) {
        fprintf(stderr, "%s: %s `%s': ",
                "main",
                (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
                fname);
        perror(NULL);
        exit(EXIT_FAILURE);
    }

    fclose(fp);

    /* Construct the suffix array. */
    fprintf(stderr, "%s: %ld bytes ... \n", fname, n);
    if(divsufsort(T, SA, (int)n) != 0) {
        fprintf(stderr, "%s: Cannot allocate memory.\n", "main");
        exit(EXIT_FAILURE);
    }

    /* Construct BWT */
    divbwt(T,bwt,NULL,n);

    /* Construct sorted array F */
    getSortedT(T,SA,sortedT,n);

    /* Construct C array */
    getC(sortedT,C,n);

    free(sortedT);

    while(true){
        int mode;
        fprintf(stdout, "\nInput mode:\n1 - make patterns\n2 - build FMIndex \n3 - FMTree search\n\n");
        fprintf(stdout, "Choose a method: ");
        scanf("%d", &mode);

        if(mode == 1){
            P = make_pattern();
            upperCase(P);

        }else if(mode == 2){

            //---------------- create walvet tree for bwt -------------

            rank = create_wavelet_tree(bwt, n);
            printf("\nCount found %d locations.\n",count(C, P, &sp, &ep, n, rank));
            createSSA(SA,bwt,B,SSA,D,n);
            free(SA);

            //---------------- create binary walvet tree for bwt -------------

            rank_b = create_binary_wavelet_tree(B, n);

        }else if(mode == 3) {

			if(rank==NULL || rank_b==NULL){
				fprintf(stdout,"\nPlease build index before search!\n");
				continue;
			}
			
			if(P==NULL){
                fprintf(stdout, "\nPlease make pattern before search!\n");
				continue;
			}
            struct timeval tvalStart, tvalEnd;
            gettimeofday(&tvalStart, NULL);
            std::set<int> R=FM_tree(bwt, C, B, SSA, P, &sp, &ep, n, D, rank, rank_b);
            gettimeofday(&tvalEnd, NULL);

            printf("\n\n*********************Result*********************");
            printf("\nSearching time: %ld microsec.", ((tvalEnd.tv_sec- tvalStart.tv_sec)*1000000L + tvalEnd.tv_usec) - tvalStart.tv_usec);
            printf("\nNumber of matched locations = %ld", R.size());
            printf("\n************************************************\n");

            make_result_file(R);
            break;

        }else{
            fprintf(stdout, "Do not input required options! FMtree will exit ...\n");
            break;
        }

    }

    free(C);
    free(T);
    free(bwt);
    free(B);
    free(SSA);

    calculate_memory_usage();

    return 0;
}
