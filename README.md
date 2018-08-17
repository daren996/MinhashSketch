# MinhashSketch

It is a project when I was internship at the University of Washington in St. Louis under the guidance of [Prof. Buhler](https://www.cse.wustl.edu/~jbuhler/).

## What is MinhashSketch

It provides a complete program which implements a parallel algorithm to get  minhash sketch with DNA sequence as input under [CUDA](http://supercomputingblog.com/cuda-tutorials/) and [CUB](https://nvlabs.github.io/cub/).  

[Minhash](https://en.wikipedia.org/wiki/MinHash) is a technique for quickly estimating how similar two sets are. For More information of similarity and minhash sketch, you can see [Class Note 24](https://classes.engineering.wustl.edu/cse584a/notes/class24.pdf) and [Class Note 25](https://classes.engineering.wustl.edu/cse584a/notes/class25.pdf)

### Steps

In general, the process divides the long DNA into sevelral CHUNKs and cope with them sequentially. For each CHUNK, it will get hash velues, radix sort, mark, scan and merge in GPU. 

<p>
<img src="./git_picture/Steps1.png" width="400" align=center />

### Input

1. **DNA sequence**. (e.g. we used [E.coli K12 genome](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP010438.1?report=fasta) as experimental data)
2. **k**: The length of subsequences (or [k-mers](https://en.wikipedia.org/wiki/K-mer)). 
3. **m**: The number of sketches stored. 
4. **t**: The number of hash functions. 

### Output

The output should be the minhash sketch of input DNA sequences and the similarity of them. 

<pr>
<img src="./git_picture/running_example1.png" width="1000" align=center />

## How do MinhashSketch work

### Subsequence

We use a pointer to access sequence, which need noly to record the last base read, and directly assign four bases(A, C, T, G) to four kinds of 2-bit values in binary. 
If the length of the subsequence exceeds 32 (i.e. one uint64 is not enough to indicate the length of subsequences), we need to use more than one ([k / 32] + 1) uint64 to store subsequences.

For example, k is 100 here:

<img src="./git_picture/Subsequence1.png" width="1000" align=center />


## Implemented on the GPU

### Process Description

#### Divide DNA sequence into several chunks

We pre-defined parameters BLOCKS\_NUM, BLOCK\_THREADS and ITEMS\_PER\_THREAD. Then we can process BLOCKS\_NUM * BLOCK\_THREADS * ITEMS\_PER\_THREAD values at a time, which is also the length of a chunk.

If the length of a chunk is (BLOCKS\_NUM * BLOCK\_THREADS * ITEMS\_PER\_THREAD), the length of a sebsequence of DNA is (BLOCKS\_NUM * BLOCK\_THREADS * ITEMS\_PER\_THREAD + k - 1). And the number of chunks is 

	if (length % (BLOCKS_NUM * BLOCK_THREADS * ITEMS_PER_THREAD) == 0)
        CHUNKS_NUM = (length - k + 1) / (BLOCKS_NUM * BLOCK_THREADS * ITEMS_PER_THREAD);
    else
        CHUNKS_NUM = (length - k + 1) / (BLOCKS_NUM * BLOCK_THREADS * ITEMS_PER_THREAD) + 1;



### Merge Between Blocks

In general, merge two ordered lists into one need O(m + n) time complexity. But we could use a efficient algorithm to merge two ordered lists with O(log(n)). 

	Given two ordered Lists A and B
	rank(x | A) : The number of elements in A that are not greater than x.
	rank(B | A) : An array (r1,r2,..., rn), where ri = rank(B[i]:A).
	rank(x | AUB) : The number of elements in the AUB that are not greater than x, so rank(x | AUB) = rank(X | A) + rank(x | B).
	Further : rank(A | AUB)= rank(A | A) + rank(A | B).

rank(A | AUB) is offset vector of elements in A.
We could get rank(A | AUB) and rank(B | AUB) with time complexity O(log(n)), where n is number of elements in A (or B). 
Then write back m smallest values according to offset.

And I considered the duplicating situation, that is, there are values in the two lists are duplicates. In this case, we can give it a mark, do not write back and then reduce the offset of values which bigger than that the mark value by one.  



-------------------
It is a project when I was internship at the University of Washington in St. Louis under the guidance of [Prof. Buhler](https://www.cse.wustl.edu/~jbuhler/).

 