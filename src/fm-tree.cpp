#include <string.h>
#include "fm-tree.h"
#include <queue>
#include <bitset>
#include <set>
#include "lib/rank_select.h"
using namespace std;

unsigned int modeBWT=3758096384; // bit mask 111000..00(28 zeros)
unsigned int modeSA=536870911;  // bit mask 000111..11(28 ones)

/*
	Setting bit in array A on position k
	@Author:Andrea Bernat
*/
void set_bit(int A[], int k){
    A[k/32] |= 1 << (k%32);
}

/*
	Remove/clear bit in array A on position k
	@Author:Andrea Bernat
*/
void clear_bit(int A[],int k){
    A[k/32] &= ~(1 << (k%32));
    //intValue &= ~(1 << bitPosition);
}

/*
	Check bit in array A on position k,
	if bit is set returns 1
	else returns 0
	@Author:Andrea Bernat
*/
int test_bit(int A[],int k){
    return ( (A[k/32] & (1 << (k%32))) !=0 );
}

/*
	Create sorted char array sortedT by sorting
	char array T with help of SA array
	@Author:Anel Hadzimuratagic
*/
void getSortedT(const unsigned char *T, int *SA,unsigned char *sortedT, int n){
    for(int i=0;i<n;i++){
        sortedT[i]=T[SA[i]];
    }
}

/*
	Create array C with 5 elements:
		- occurrences of A,C,G,T and $ elements
	@Author: Robert Jambrecic
*/
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

/*
	Calculates occurrences of char c in char array T
	from start of the array to index n.
	Returns number of occurrences.
	@Author: Robert Jambrecic
*/
int occ(const unsigned char *T, char c, int n){
    int count=0;

    for(int i=0;i<=n;i++){
        if(tolower((char)T[i])==c){
            count+=1;
        }
    }
    return count;
}

/*
	Calculates linear rank operation of char c
	in BWT char array on position index.
	Returns rank.
	@Author: Robert Jambrecic
*/
int rankOp(char s,const unsigned char *bwt,int index){
    return occ(bwt,s,index-1);
}

/*
	Converting char to integer.
	Returns integer presented by char c.
	@Author: Andrea Bernat
*/
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

/*
	Calculates LF operation.
	Returns result of operation.
	@Author: Anel Hadzimuratagic
*/
int LFoperation(int *C,char s,const unsigned char *bwt,int index){
    return C[acgtToInt(s)]+rankOp(s,bwt,index);
}

/*
	Counting number of set bits in B array
	from start to index.
	Returns number of occurrences.
	@Author: Andrea Bernat
*/
int occBin(int* B,int index){
    int occ=0;
    for(int i=0;i<=index;i++){
        if(test_bit(B,i)==1){
            occ++;
        }
    }
    return occ;
}

/*
	Linear binary rank operation.
	Returns rank of set 1 in array B to position index.
	@Author: Andrea Bernat
*/
int rankBin(int* B, int index){
    return occBin(B,index-1);
}

/*
	Creates SSA and bit array B with sampled elements from SA array.
	Where D is sampling distance.
	@Author: Anel Hadzimuratagic
*/

void createSSA(int *SA, unsigned char *bwt, int *B,unsigned int *SSA,int D, int n){
    int j=0;
    for(int i=0; i<n;i++){
        if(SA[i]%D==0){
            set_bit(B,i);
            SSA[j]=SA[i];

            // setting on first 3 bits bwt[i] so that can be used in early leaf calculation
            if(bwt[i]=='$'){  // 100 2^31 - 4
                SSA[j]=SSA[j] | 2147483648;
            }else if(bwt[i]=='C'){ // 001 2^29 - 1
                SSA[j]=SSA[j] | 536870912;
            }else if(bwt[i]=='G'){ // 010 2^30 - 2
                SSA[j]=SSA[j] | 1073741824;
            } else if(bwt[i]=='T'){ // 011 (2^29 + 2^30) - 3
                SSA[j]=SSA[j] | 1610612736;
            }
            j++;
        }else{
            clear_bit(B,i);
        }
    }
}

/*
	Implementation of early leaf node calculation to avoid
	expensive D-1 step in FM-tree.
	Returns found set of locations of pattern P in target array T
	in D-1 step of FM tree.
	@Author: Anel Hadzimuratagic
*/
set<int> early_leaf_node(int *C, char P[], unsigned int *SSA, int *B, int n, rank_select* rank, rank_select *rank_b){
    int s;
    int sp = 0;
    int ep = 0;
    int i = strlen(P);
    char P_1[i-1];
    set<int> R;

    if(i==1){
        return R;
    }

    string str = P;
    strcpy(P_1, str.substr(1,i).c_str());

    count(C, P_1, &sp, &ep, n, rank);

    // preparation to compare with first 3 bits in SSA
    if(P[0]=='A'){
        s=0; // 000
    } else if(P[0]=='C'){
        s=1; // 001
    } else if(P[0]=='G'){
        s=2; // 010
    } else if(P[0]=='T'){
        s=3; // 011
    } else if(P[0]=='$'){
        s=4; // 100
    }

    int ssp=rank_b->rank('1', sp - 1);
    int sep=rank_b->rank('1', ep);

    unsigned int SA_element;
    unsigned int BWT_element;

    if(ssp == sep && test_bit(B, sp) == 1){
        sep+=1;
    }
    for(int i=ssp; i < sep; i++){
        SA_element=(SSA[i]&modeSA);
        BWT_element=(SSA[i]&modeBWT)>>29;

        if(SA_element>=0 && BWT_element==s){
            R.insert(SA_element-1);
        }
    }
    return R;
}

/*
	Counts number of pattern P occurrences
	in target array T, from bwt array.
	Returns occurrence number and range [sp,ep] such that
	sortedT[sp,ep] includes all sortedT rows prefixed by P.
	@Author: Robert Jambrecic
*/
int count(int *C,char* P, int *sp, int *ep ,int n, rank_select *rank){
    int len= strlen(P) - 1;
    char s=P[len--];
    int si=acgtToInt(s);

    //if there is no occurrence or bad char return 0
    if(C[si]==0){
		
		//set sp and ep that result *ep-*sp+1=0
		*sp=1;
		*ep=0;
        return 0;
    }
    *sp=C[si];

    //look for next letter that appears at least once.
    while(C[si+1]==0){
        si++;
    }

    if((si+1)>=5){
        *ep=n-1;
    }else{
        *ep=C[si+1]-1;
    }

    while(len >= 0 && *sp <= *ep ){
        s=P[len--];
        *sp= C[acgtToInt(s)] + rank->rank(toupper(s), *sp - 1);
        *ep= C[acgtToInt(s)] + rank->rank(toupper(s), *ep) - 1;
    }

    if(*sp>*ep){
        return 0;
    }else{
        return *ep-*sp+1;
    }
}

/*
	FM index locate function.
	Calculates locations of pattern P occurrences.
	Returns set of found locations.
	@Author: Andrea Bernat
*/
set<int> locate(const unsigned char *bwt, int *C, int *B, unsigned int *SSA, int sp, int ep, set<int> R, rank_select *rank, rank_select *rank_b) {
    int i=0;
    int j=0;
    int m=0;
    int SA_element;

    for(i=sp;i<=ep;i++){
        m=0;
        j=i;
        while(test_bit(B,j)!=1){
            //LF operation
            j= C[acgtToInt(bwt[j])] + rank->rank(bwt[j], j - 1);
            m+=1;
        }
        // obtain SSA element
        SA_element= SSA[rank_b->rank('1', j - 1)] & modeSA;
        R.insert(SA_element+m);
    }
    return R;
}

// trashold can be changed
#define treshold 10
#define s_a 1
#define s_c 2
#define s_g 3
#define s_t 4


/*
	Implementation of FM-tree function calculation.
	Calculates locations of pattern P occurrences.
	Returns set of found locations.
	Optimized FM index locate function.
	@Author: Anel Hadzimuratagic, Andrea Bernat
*/
set<int> FM_tree(const unsigned char *bwt,int *C,int *B,unsigned int* SSA, char P[], int *sp, int *ep ,int n,int D, rank_select *rank,rank_select *rank_b) {
    int sp_child[4];
    int ep_child[4];
    int nextNodeSuccess[4];
    int total_num = *ep - *sp + 1;
    int tree_height = D - 1;
    int layer=0;
    int ssp=0;
    int sep=0;
    unsigned int SA_element;
    int i=0;

    if(total_num == 0){
        return set<int>();
    }

    set<int> R = early_leaf_node(C, P, SSA, B, n, rank, rank_b);

    Node node(*sp, *ep, layer);
    queue <Node> nodeQueue;
    nodeQueue.push(node);

    while(!nodeQueue.empty() && R.size()<total_num){

        node=nodeQueue.front();
        nodeQueue.pop();

        *sp=node.getSp();
        *ep=node.getEp();
        layer=node.getLayer();

        if(*ep-*sp+1<treshold){
            R = locate(bwt,C,B,SSA,*sp,*ep,R,rank,rank_b);
        }else{
            ssp=rank_b->rank('1',*sp-1);
            sep=rank_b->rank('1',*ep);

            if((ssp == sep) && test_bit(B, *sp)==1){
                sep++;
            }

            for(int k=ssp;k<sep;k++){
                SA_element=(SSA[k]&modeSA);
                R.insert(SA_element+layer);
            }

            if((layer+1)<=tree_height){
                if(C[s_a]!=0){
                    nextNodeSuccess[0]=1;
                    sp_child[0]=C[s_a]+rank->rank('A',*sp-1);
                    ep_child[0]=C[s_a]+rank->rank('A',*ep)-1;
                }
                if(C[s_c]!=0){
                    nextNodeSuccess[1]=1;
                    sp_child[1]=C[s_c]+rank->rank('C',*sp-1);
                    ep_child[1]=C[s_c]+rank->rank('C',*ep)-1;
                }
                if(C[s_g]!=0){
                    nextNodeSuccess[2]=1;
                    sp_child[2]=C[s_g]+rank->rank('G',*sp-1);
                    ep_child[2]=C[s_g]+rank->rank('G',*ep)-1;
                }
                if(C[s_t]!=0){
                    nextNodeSuccess[3]=1;
                    sp_child[3]=C[s_t]+rank->rank('T',*sp-1);
                    ep_child[3]=C[s_t]+rank->rank('T',*ep)-1;
                }

                for(i=0; i <= 3; i++){
                    if(nextNodeSuccess[i] == 1 && sp_child[i] != -1 && ep_child[i] != -1 && sp_child[i] <= ep_child[i]) {
                        Node nodeNew(sp_child[i], ep_child[i], layer + 1);
                        nodeQueue.push(nodeNew);
                    }
                    nextNodeSuccess[i]=0;
                }
            }
        }
    }
    return R;
}
