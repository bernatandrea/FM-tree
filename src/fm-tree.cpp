#include <string.h>
#include "fm-tree.h"
#include <queue>
#include <bitset>
#include <set>
#include "lib/rank_select.h"

using namespace std;

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
void createSSA(int *SA, int *B,int *SSA,int D, int n){
    int j=0;
    for(int i=0; i<n;i++){
        if(SA[i]%D==0){
            set_bit(B,i);
            SSA[j++]=SA[i];
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
std::set<int> early_leaf_node(const unsigned char *T, const unsigned char *bwt,int *C, char P[], int *SSA, int *B, int n, rank_select* t,rank_select *tb){
    int sp = 0;
    int ep = 0;
    int i = strlen(P);
    char P_1[i-1];
    std::set<int> R;

    if(i==1){
        return R;
    }

    string s = P;
    strcpy(P_1, s.substr(1,i).c_str());

    count(bwt, C, P_1, &sp, &ep, n, t);

    for(int j = sp; j<=ep; j++){
        if(test_bit(B,j)==1){
            int index = tb->rank('1', j-1);
            if(P[0] == T[SSA[index]-1]){
                R.insert(SSA[index]-1);
            }
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
int count(const unsigned char *bwt,int *C,char P[], int *sp, int *ep ,int n, rank_select *t){
    int i=strlen(P)-1;
    int s=P[i--];
    int si=acgtToInt(s);

    //if there is no occurrence return -1
    if(C[si]==0 && si!=0){
        return -1;
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

/*
	FM index locate function.
	Calculates locations of pattern P occurrences.
	Returns set of found locations.
	@Author: Andrea Bernat
*/
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

/*
	Implementation of FM-tree function calculation.
	Calculates locations of pattern P occurrences.
	Returns set of found locations.
	Optimized FM index locate function.
	@Author: Anel Hadzimuratagic, Andrea Bernat
*/
std::set<int> FM_tree(const unsigned char *bwt,const unsigned char *T,int *C,int *B,int* SSA, char P[], int *sp, int *ep ,int n, rank_select *ran,rank_select *tb) {
    int total_num = *ep - *sp + 1;
    int layer=0;
    int sp_child[4];
    int ep_child[4];
    int t=0;
    int ssp=0;
    int sep=0;

    int boolChar[4];

    int v_sp=*sp;
    int v_ep=*ep;

    if(total_num == 0){
        return std::set<int>();
    }
    std::set<int> R = early_leaf_node(T, bwt, C, P, SSA, B, n, ran, tb);

    int tree_height = D - 1;
    Node node(v_sp, v_ep, layer);
    queue <Node> nodeQueue;
    nodeQueue.push(node);

    while(!nodeQueue.empty() && R.size()<total_num){

        boolChar[0]=0;
        boolChar[1]=0;
        boolChar[2]=0;
        boolChar[3]=0;

        node=nodeQueue.front();
        nodeQueue.pop();

        v_sp=node.getSp();
        v_ep=node.getEp();
        layer=node.getLayer();

        //ssp=rankBin(B,*sp);
        //sep=rankBin(B,*ep);
        //num=num+sep-ssp+1;
        if(v_ep-v_sp+1<treshold){
            R = locate(bwt,C,B,SSA,v_sp,v_ep,R,ran,tb);

        }else{
            ssp=tb->rank('1',v_sp-1);
            sep=tb->rank('1',v_ep);

            if((ssp == sep) && test_bit(B, v_sp)==1){
                sep++;
            }

            for(int k=ssp;k<sep;k++){
                R.insert(SSA[k]+layer);
            }

            if((layer+1)<=tree_height){
                if(C[s_a]!=0){
                    boolChar[0]=1;
                    sp_child[0]=C[s_a]+ran->rank('A',v_sp-1);
                    ep_child[0]=C[s_a]+ran->rank('A',v_ep)-1;
                }
                if(C[s_c]!=0){
                    boolChar[1]=1;
                    sp_child[1]=C[s_c]+ran->rank('C',v_sp-1);
                    ep_child[1]=C[s_c]+ran->rank('C',v_ep)-1;
                }
                if(C[s_g]!=0){
                    boolChar[2]=1;
                    sp_child[2]=C[s_g]+ran->rank('G',v_sp-1);
                    ep_child[2]=C[s_g]+ran->rank('G',v_ep)-1;
                }
                if(C[s_t]!=0){
                    boolChar[3]=1;
                    sp_child[3]=C[s_t]+ran->rank('T',v_sp-1);
                    ep_child[3]=C[s_t]+ran->rank('T',v_ep)-1;
                }

                for(t=0;t<=3;t++){
                    if(boolChar[t]==1 && sp_child[t]!=-1 && ep_child[t]!=-1 && sp_child[t]<=ep_child[t]) {
                        Node nodeNew(sp_child[t], ep_child[t], layer + 1);
                        nodeQueue.push(nodeNew);
                    }
                }
            }
        }
    }
    return R;
}
