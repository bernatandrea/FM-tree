#include <sys/time.h>
#include "divsufsort.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fm-tree.h"
#include <math.h>
#include <unordered_map>
#include <unistd.h> /* getopt */
#include <vector>
#include "lib/bitmask.h"
#include "lib/bitmask_bitset.h"
#include "lib/data.h"
#include "lib/rank_select.h"
#include "lib/rb_tree.h"

using namespace std;

rank_select* create_wavelet_tree(const unsigned char *bwt){

    std::string content((const char *)bwt);
    std::function<bitmask *(uint32_t)> c;
    std::unordered_map<char, uint32_t> counters;
    uint32_t word_size = 23508030;
    c = bitmask_bitset::create;
    bitmask::set_creator(c);
    std::vector<data*> data_chunks;

    if (content.length()) {
        data_chunks.push_back(data::create(content, word_size, counters));
    }

    return new rb_tree(data_chunks);

}


rank_select* create_binary_wavelet_tree(int *B, int n){

    string b_str;
    for(int i=0;i<n;i++){
        b_str += to_string(test_bit(B,i));
    }

    std::vector<data *> data_chunks_bin;
    std::function<bitmask *(uint32_t)> c_bin;
    std::unordered_map<char, uint32_t> counters_bin;
    uint32_t word_size_bin = 2358030;
    c_bin = bitmask_bitset::create;
    bitmask::set_creator(c_bin);

    if (b_str.length()) {
        data_chunks_bin.push_back(data::create(b_str, word_size_bin, counters_bin));
    }


    return new rb_tree(data_chunks_bin);
}


int main() {
    FILE *fp;
    const char *fname;
    long n;
    rank_select *t;
    rank_select *tb;
    int sp=0;
    int ep=0;
    char P[100];
    int D;


    /* Open a file for reading. */
    if((fp = fopen(fname = "ecoli.txt", "rb")) == NULL) {
        printf(" Cannot open file ");
        perror(NULL);
        exit(EXIT_FAILURE);
    }

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

    fprintf(stdout, "Input sampling distance: ");
    scanf("%d", &D);

   /* Allocate 5n bytes of memory. */
    unsigned char *T = (unsigned char *)malloc((size_t)n * sizeof(unsigned char));
    unsigned char *bwt = (unsigned char *)malloc((size_t)n * sizeof(unsigned char));
    unsigned char *sortedT = (unsigned char *)malloc((size_t)n * sizeof(unsigned char));
    int *SA = (int *)malloc((size_t)n * sizeof(int));
    int *SSA = (int *)malloc((size_t)(int(ceil(n/(float)D))) * sizeof(int));
    int *C = (int *)malloc((size_t)5 * sizeof(int));
    int *B = (int *)malloc((size_t)n * sizeof(int));

    if((T == NULL) ||(bwt == NULL) ||(sortedT == NULL) || (SA == NULL) || (SSA == NULL)|| (C == NULL)) {
        fprintf(stderr, "%s: Cannot allocate memory.\n", "main");
        exit(EXIT_FAILURE);
    }

    /* Read n bytes of data. */
    if(fread(T, sizeof(unsigned char), (size_t)n, fp) != (size_t)n) {
        fprintf(stderr, "%s: %s `%s': ",
                "main",
                (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
                "test-short.txt");
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
        fprintf(stdout, "\nInput mode:\n1 - create wavelet tree\n2 - FMTree search\n3 - make patterns:");
        scanf("%d", &mode);

        if(mode == 1){
            //---------------- create walvet tree for bwt -------------
            t = create_wavelet_tree(bwt);

            count(bwt,C,P,&sp,&ep,n,t);

	    createSSA(SA,B,SSA,D,n);

    	    free(SA);

            //---------------- create binary walvet tree for bwt -------------
            tb = create_binary_wavelet_tree(B, n);

        }else if(mode == 2){

            struct timeval tvalStart, tvalEnd;
            gettimeofday(&tvalStart, NULL);
            std::set<int> R=FM_tree(bwt,T,C,B,SSA,P,&sp,&ep,n,t,tb);
            gettimeofday(&tvalEnd, NULL);

            printf("\n\n\n*********************Result*********************");
            printf("\nSearching time: %ld microsec.", ((tvalEnd.tv_sec- tvalStart.tv_sec)*1000000L + tvalEnd.tv_usec) - tvalStart.tv_usec);
            printf("\nNumber of matched locations = %ld", R.size());
            printf("\n************************************************\n");
            break;

        }else if(mode == 3) {
            fprintf(stdout, "Please input the pattern for search:\n");
            scanf("%s", P);
        }else{
            fprintf(stdout, "Do not input required options! FMtree will exit ...\n");
            break;
        }

    }

    /*struct rusage *usage;
    int who = RUSAGE_SELF;
    
    getrusage(who, usage);
    printf("\n %ld", usage->ru_maxrss);*/
    

    free(C);
    free(T);
    free(bwt);
    free(B);
    free(SSA);
    

    return 0;
}
