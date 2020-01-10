#include <string.h>
#include "fm-tree.h"
#include <queue>
#include <bitset>
#include <set>
#include "lib/rank_select.h"

using namespace std;

void set_bit(int A[], int k){
    A[k/32] |= 1 << (k%32);
}

void clear_bit(int A[],int k){
    A[k/32] &= ~(1 << (k%32));
    //intValue &= ~(1 << bitPosition);
}

int test_bit(int A[],int k){
    return ( (A[k/32] & (1 << (k%32))) !=0 );
}

void getSortedT(const unsigned char *T, int *SA,unsigned char *sortedT, int n){

    for(int i=0;i<n;i++){
        sortedT[i]=T[SA[i]];
    }
}

/**
 returning C
**/
void getC(unsigned char *sortedT,int *C,int n){
    //aphabet acgt + $ = 5
    int occ_end=1; //$
    int occ_a=0; //a
    int occ_c=0;  //c
    int occ_g=0; //g
    int occ_t=0;  //t

    for(int i=1;i<n;i++){
        if(sortedT[i]=='A'){
            occ_a++;
        }else if(sortedT[i]=='C'){
            occ_c++;
        }else if(sortedT[i]=='G'){
            occ_g++;
        }else if(sortedT[i]=='T'){
            if(occ_t!=0){
                break;
            }
            occ_t++;
        }
    }

    C[0]=0; //$
    C[1]=occ_a!=0?occ_end:0;
    C[2]=occ_c!=0?(occ_end+occ_a):0;
    C[3]=occ_g!=0?(occ_end+occ_a+occ_c):0;
    C[4]=occ_t!=0?(occ_end+occ_a+occ_c+occ_g):0;

}

/**
Returning occ
**/
int occ(const unsigned char *T, char c, int n){
    int count=0;

    for(int i=0;i<=n;i++){
        if(tolower((char)T[i])==c){
            count+=1;
        }
    }
    return count;
}

/**
Returning rank
**/
int rankOp(char s,const unsigned char *bwt,int index){
    return occ(bwt,s,index-1);
}

int acgtToInt(char c){

    if(c=='$'){
        return 0;
    }else if(c=='A'){
        return 1;
    }else if(c=='C'){
        return 2;
    }else if(c=='G'){
        return 3;
    }else if(c=='T'){
        return 4;
    }
    return -1;
}

int LFoperation(int *C,char s,const unsigned char *bwt,int index){
    return C[acgtToInt(s)]+rankOp(s,bwt,index);
}

int occBin(int* B,int index){
    int occ=0;
    for(int i=0;i<=index;i++){
        if(test_bit(B,i)==1){
            occ++;
        }
    }
    return occ;
}


int rankBin(int* B, int index){
    return occBin(B,index-1);
}

void createSSA(int *SA, int *B,int *SSA,int D, int n){
    int j=0;
    for(int i=0; i<n;i++){
        if(SA[i]%D==0){
            set_bit(B,i);
            SSA[j++]=SA[i];
        } /*else{
            B[i]=0;
        }*/
    }
}

std::set<int> early_leaf_node(const unsigned char *T, const unsigned char *bwt,int *C, char P[], int *SSA, int *B, int n, rank_select* t,rank_select *tb){
    int sp = 0;
    int ep = 0;
    int i = strlen(P)-1;
    char P_1[i-1];
    std::set<int> R;
    string s = P;
    strcpy(P_1, s.substr(1,i).c_str());

    count(bwt, C, P_1, &sp, &ep, n, t);

    for(int j = sp; j<=ep; j++){
        if(test_bit(B,i)==1){
            int index = tb->rank('1', i-1);
            if(P[0] == T[SSA[index]-1]){
                R.insert(SSA[index]-1);
            }
        }
    }

    return R;
}


int count(const unsigned char *bwt,int *C,char P[], int *sp, int *ep ,int n, rank_select *t){
    int i=strlen(P)-1;
    int s=P[i--];
    int si=acgtToInt(s);

    //ako nema pojave slova izlazi
    if(C[si]==0 && si!=0){
        return -1;
    }
    *sp=C[si];

    //ako sljedece abecedno slovo ima 0 pojava trazi prvo slovo koje se ponavlja bar 1x
    while(C[si+1]==0){
        si++;
    }

    //ako ne postoji sljedece znaci da je s zadnje abecedno slovo i ide do kraja
    // ili??
    if((si+1)>=5){
        *ep=n-1;
    }else{
        *ep=C[si+1]-1;
    }

    while(i>=0 && *sp<*ep ){
        s=P[i--];
        *sp=C[acgtToInt(s)]+ t->rank(toupper(s), *sp-1);
        *ep=C[acgtToInt(s)]+t->rank(toupper(s), *ep)-1;
    }

    if(*sp>*ep){
        return 0;
    }else{
        return *ep-*sp+1;
    }
}

/**
 * Returning num of located objects
 * **/
std::set<int> locate(const unsigned char *bwt,int *C,int *B,int *SSA, int sp, int ep,set<int> R, rank_select *ran,rank_select *tb) {
    int i=0;
    int j=0;
    int m=0;
    for(i=sp;i<=ep;i++){
        j=i;m=0;
        while(B[j]!=1){
            //LF operation
            j=C[acgtToInt(bwt[j])]+ran->rank(bwt[j],j-1);
            m+=1;
        }
        R.insert(SSA[tb->rank('1',j-1)]+m);
    }
    return R;
}


#define D 4
#define treshold 0
#define s_a 1
#define s_c 2
#define s_g 3
#define s_t 4
std::set<int> FM_tree(const unsigned char *bwt,const unsigned char *T,int *C,int *B,int* SSA, char P[], int *sp, int *ep ,int n, rank_select *ran,rank_select *tb) {
    int total_num = *ep - *sp + 1;
    int num = 0;
    int layer=0;
    int R_size=0;
    int sp_child[4];
    int ep_child[4];
    int t=0;
    int ssp=0;
    int sep=0;

    if(total_num == 0){
        return std::set<int>();
    }
    std::set<int> R = early_leaf_node(T, bwt, C, P, SSA, B, n, ran, tb);
    R_size=R.size();
    num+=R_size;

    int tree_height = D - 1;
    Node node(*sp, *ep, layer);
    queue <Node> nodeQueue;
    nodeQueue.push(node);

    while(!nodeQueue.empty() && num<total_num){
        node=nodeQueue.front();
        *sp=node.getSp();
        *ep=node.getEp();
        layer=node.getLayer();

        ssp=rankBin(B,*sp);
        sep=rankBin(B,*ep);
        num=num+sep-ssp+1;
        if(ep-sp+1<treshold){
            /**
             * klasicni locate u R upisuje lokacije i vraca broj upisanih
             * **/
            int old_size=R.size();
            R = locate(bwt,C,B,SSA,*sp,*ep,R,ran,tb);
            int new_size=R.size();
            num+=(new_size-old_size);
        }else{
        //
        ssp=tb->rank('1',*sp-1);
        sep=tb->rank('1',*ep-1);
        num=num+sep-ssp+1;

        for(int k=ssp;k<=sep;k++){
            R.insert(SSA[k]+layer);
            R_size++;
            if(R_size>total_num){
                return R;
            }
        }

        if((layer+1)<tree_height){
            if(C[s_a]!=0){
                sp_child[0]=C[s_a]+ran->rank('A',*sp-1);
                ep_child[0]=C[s_a]+ran->rank('A',*ep)-1;
            }
            if(C[s_c]!=0){
                sp_child[1]=C[s_c]+ran->rank('C',*sp-1);
                ep_child[1]=C[s_c]+ran->rank('C',*ep)-1;
            }
            if(C[s_g]!=0){
                sp_child[2]=C[s_g]+ran->rank('G',*sp-1);
                ep_child[2]=C[s_g]+ran->rank('G',*ep)-1;
            }
            if(C[s_t]!=0){
                sp_child[3]=C[s_t]+ran->rank('T',*sp-1);
                ep_child[3]=C[s_t]+ran->rank('T',*ep)-1;
            }

            for(t=0;t<=3;t++){
                if(sp_child[t]!=-1 && ep_child[t]!=-1 && sp_child[t]<=ep_child[t]) {
                    Node nodeNew(sp_child[t], ep_child[t], layer + 1);
                    nodeQueue.push(nodeNew);
                }
            }
        }
    }
}
return R;
}
