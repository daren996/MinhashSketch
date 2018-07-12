# MinhashSketch

## Some improvements when implementing the algorithm
When meeting with Prof. Buhler, I got several tips about the code.

### Simplify the process of taking k-mers

I used hash template struct of STL to transfer sequence of DNA to integer before. 
This will cost O(k*(l-k+1)) complexity totally which is terrible inefficient. 

Directly assign four bases to 0-3 four digits, and use the first two pointers to access the sequence.


-------------------
This is a project when I was internship at the University of Washington in St. Louis under the guidance of Prof. [Buhler](https://www.cse.wustl.edu/~jbuhler/).

 