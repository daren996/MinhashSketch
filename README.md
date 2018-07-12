# MinhashSketch

## Some improvements when implementing the algorithm
When meeting with Prof. Buhler, I got several tips about the implemention.

### Simplify the process of taking k-mers subsequence

I used hash template struct of STL to transfer sequence of DNA to integer before. 
This will cost O[k*(l-k+1)] complexity to get each base and O(l-k+1) to caculate hash values totally which is terrible inefficient. 

We can use a pointer to access this sequence, and directly assign four bases(A, C, T, G) to four kinds of 2-bit values in binary.

Such as:
<img src="./git_picture/process_of_taking_k-mers1.png" width="350" align=center /> 
k is 5 here.

You will get the first k-mer subsequence, ATTGC, then transfer to 0011111001 in binary. Then move pointer to next pot. Then you can get the second k-mer, TTGCC. What you need to do now is left shift 2 bits of last k-mer in binary and add the 2-bit values in binary of new base.

<img src="./git_picture/process_of_taking_k-mers2.png" width="343" align=center />

As a result, you only need O(l-k+1) take-operations and no hash operations to transfer sequence to integer.

-------------------
This is a project when I was internship at the University of Washington in St. Louis under the guidance of [Prof. Buhler](https://www.cse.wustl.edu/~jbuhler/).

 