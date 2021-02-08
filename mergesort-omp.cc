/**
 *  \file mergesort.cc
 *
 *  \brief Implement your mergesort in this file.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <algorithm>
#include <iostream>
#include "sort.hh"
void merge(keytype* X, int n, keytype* tmp) {
   int i = 0;
   int j = n/2;
   int ti = 0;

   while (i<n/2 && j<n) {
      if (X[i] < X[j]) {
         tmp[ti] = X[i];
         ti++; i++;
      } else {
         tmp[ti] = X[j];
         ti++; j++;
      }
   }
   while (i<n/2) { 
      tmp[ti] = X[i];
      ti++; i++;
   }
      while (j<n) { 
         tmp[ti] = X[j];
         ti++; j++;
   }
   memcpy(X, tmp, n*sizeof(keytype));

} // end of merge()

void mergesort(keytype* X, int n, keytype* tmp)
{
   if (n < 2) return;

   mergesort(X, n/2, tmp);

   mergesort(X+(n/2), n-(n/2), tmp);

   
   merge(X, n, tmp);
}


int BSearch(int x, keytype* T, int p, int r)
{
        int low = p;
        int high = std::max(p,r+1);
        while (low < high) {
                int mid = (low + high) / 2;
                if (x <= T[mid])
                        high = mid;
                else low = mid + 1;
        }
        return high;
}



void PMerge(keytype* T,int l_bound1,int r_bound1,int l_bound2,int r_bound2, keytype* A,int left_bound3, int depth)
{
	int n1 = r_bound1 - l_bound1 +1;
	int n2 = r_bound2 - l_bound2 +1;
	
	if (n1<n2)
	{
		int temp = 0;
		temp = l_bound1;
		l_bound1 = l_bound2;
		l_bound2 = temp;
		
		temp = r_bound1;
		r_bound1 = r_bound2;
		r_bound2 = temp;
		
		temp = n1;
		n1 = n2;
		n2 = temp;
	}
	
	if (n1 == 0)
		return;
		
	else
	{	
			int q1 = (l_bound1 + r_bound1)/2;
			int q2 = BSearch(T[q1],T,l_bound2,r_bound2);
			int q3 = left_bound3 + (q1 - l_bound1) + (q2 - l_bound2);
			A[q3] = T[q1];
		
		if (depth>1){
#pragma omp task
					PMerge(T,l_bound1,q1 - 1,l_bound2,q2 - 1, A, left_bound3,depth/2);
					PMerge(T, q1 + 1, r_bound1,q2,r_bound2,A,q3+1,depth/2);
#pragma omp taskwait
		}

		else
		{
			
			PMerge(T,l_bound1,q1 - 1,l_bound2,q2 - 1, A, left_bound3,0);
			PMerge(T, q1 + 1, r_bound1,q2,r_bound2,A,q3+1,0);
		}
	}
}

void parallel_mergesort(keytype* A, int left_bound , int right_bound, keytype* B, int left_bound2, int depth)
{
	int size = right_bound - left_bound + 1;
	if (size < 1)
		return;
	if (size == 1) {
		B[left_bound2] = A[left_bound];
	}
	
	else 
	{ 
		int new_right = (right_bound + left_bound)/2; //upper bound for first call lower bound for right call
		int  new_left = new_right - left_bound + 1;//lower bound for temp array
		if(depth>1) {
			keytype* T = newKeys(size);
#pragma omp task
					parallel_mergesort(A,left_bound, new_right, T, 0, depth/2); 
					parallel_mergesort(A, new_right + 1,right_bound,T,new_left,depth/2);
#pragma omp taskwait
		PMerge(T, 0, new_left -1, new_left, size-1,B,left_bound2,depth);
		free(T);
				}
		else
		{
			
			keytype* temp = newKeys(size);
			mergesort(A+left_bound,size/2,temp);
			mergesort(A+new_right+1, size-size/2,temp+(size - size/2));
			PMerge(temp,0,new_left-1,new_left,size-1,B,left_bound2,0);			
		}
	}
	

}
void
mySort (int N, keytype* A)
{

	keytype* keys = newKeys(N);
	//std::cout<<A[0]<<"\n";
	//A = A+1;
	//std::cout<< A[0]<<"\n";
#pragma omp parallel
        {
		int threads = omp_get_num_threads();
#pragma omp single
			parallel_mergesort(A, 0, N - 1, A, 0, threads);
        }
	//mergesort(A,N,N);

	//mergesort(A,N,keys);
//	for(int i=0;i<N;i++){
//		std::cout<<A[i]<<"\n";
//				}	
}

/* eof */
