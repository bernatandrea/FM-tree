#ifndef FMTREE_FM_TREE_H
#define FMTREE_FM_TREE_H

#include <set>
#include "lib/bitmask.h"
#include "lib/bitmask_bitset.h"
#include "lib/bitmask_vector.h"
#include "lib/data.h"
#include "lib/rank_select.h"
#include "lib/lookup_list.h"
#include "lib/rb_tree.h"

class Node
{
private:
    int sp;
    int ep;
    int layer;

public:

    Node(int sp, int ep, int layer) { // Constructor with parameters
        this->sp = sp;
        this->ep = ep;
        this->layer = layer;
    }

    int getSp() { return sp; }
    int getEp() { return ep; }
    int getLayer()  { return layer; }

    void setSp(int sp) { this->sp=sp; }
    void setEp(int ep) { this->ep=ep; }
    void setLayer(int layer)  { this->layer=layer; }
};

int test_bit(int* B, int k);
void getSortedT(const unsigned char *T,int *SA, unsigned char *sortedT, int n);
void getC(unsigned char *sortedT,int *C, int n);
int occ(const unsigned char *T, char c, int n);
int rankOp(char s,const unsigned char *bwt,int index);
int LFoperation(int *C,char s,const unsigned char *bwt,int index);
int count(int *C,char* P, int *sp, int *ep ,int n, rank_select *rank);
void createSSA(int *SA, unsigned char *bwt, int *B,unsigned int *SSA,int D, int n);
std::set<int> early_leaf_node(int *C, char P[], unsigned int *SSA, int *B, int n, rank_select* rank, rank_select *rank_b);
std::set<int> FM_tree(const unsigned char *bwt,int *C,int *B,unsigned int* SSA, char P[], int *sp, int *ep ,int n,int D, rank_select *rank,rank_select *rank_b);
std::set<int> locate(const unsigned char *bwt,int *C,int *B,int *SSA, char P[],int *sp, int *ep,rank_select *rank, rank_select *rank_b);

#endif //FMTREE_FM_TREE_H