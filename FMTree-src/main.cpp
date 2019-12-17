#include <iostream>
#include <sys/time.h>
#include "divsufsort.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "divsufsort.h"
#include "fm-tree.h"
#include <queue>
#include <math.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unistd.h> /* getopt */
#include <vector>
#include <chrono>
#include "lib/bitmask.h"
#include "lib/bitmask_bitset.h"
#include "lib/bitmask_vector.h"
#include "lib/data.h"
#include "lib/rank_select.h"
#include "lib/lookup_list.h"
#include "lib/rb_tree.h"

using namespace std;


int main() {
    FILE *fp;
    const char *fname;
    long n;
    int D=4;
    clock_t start, finish;
    clock_t start_rank, finish_rank;

    char P[]="ACTGGAGAAAGCAGCCCGGA";
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
    fprintf(stderr, "%s: %ld bytes ... ", fname, n);
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



    //---------------- create walven tree for bwt -------------

    start=clock();

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

    rank_select *t;
    t = new rb_tree(data_chunks);

    //-----------------------------------------
    /* Testing count method */
    int sp=0;
    int ep=0;
    int rez=count(bwt,C,P,&sp,&ep,n,t);
    printf("\nCount = %d",rez);

    createSSA(SA,B,SSA,D,n);
    free(SA);

    // ------------- create walven tree for binary ------
    string b_str;
    for(int i=0;i<n;i++){
        b_str+=to_string(test_bit(B,i));
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

    rank_select *tb;
    tb = new rb_tree(data_chunks_bin);

    finish=clock();
    fprintf(stderr, "\nStvaranje walveta: %.4f milisec\n", ((double)(finish - start) / (double)CLOCKS_PER_SEC)*1000);

    //----------------------------------------------------

    printf("\nFM tree\n");

    struct timeval tvalStart, tvalEnd;
    gettimeofday(&tvalStart, NULL);
    //start_rank = clock();
    std::set<int> R=FM_tree(bwt,T,C,B,SSA,P,&sp,&ep,n,t,tb);
    gettimeofday(&tvalEnd, NULL);
    //finish_rank = clock();

    //fprintf(stderr, "\nFM tree: %.4f milisec\n", ((double)(finish_rank - start_rank) / (double)CLOCKS_PER_SEC)*1000);

    printf("\\nFM tree: %ld", ((tvalEnd.tv_sec- tvalStart.tv_sec)*1000000L + tvalEnd.tv_usec) - tvalStart.tv_usec);
    printf("\nFINISH\n");

    // print locations
    /*for (auto it=R.begin(); it != R.end(); ++it)
        cout << ' ' << *it;*/

    free(C);
    free(T);
    free(bwt);
    free(B);
    free(SSA);

    return 0;
}