# count

This example demonstrates how to use `Counter` classes to calculate the number of the following combinatorial objects.

* Binary matrices with prescribed row and column sums
* Binary matrices whose row and column sums lie in prescribed intervals
* Perfect and near-perfect matchings in bipartite graphs

## Usage

    count TYPE INSTANCE
    
The `INSTANCE` format depends on the associated problem `TYPE`, which must be substituted by one of the following strings.

* `matching`: Count the number of perfect and near perfect matchings in the bipartite graph `G=(U,V,E)` that is specified by the `INSTANCE` encoding. 

  An `INSTANCE` has the form `([0|1]^n)^n`, where n is a positive integer. Such a 0-1 string is interpreted as the n times n bi-adjacency matrix `M = (m_ij)` of a bipartite graph `G=(U,V,E)` with `|U|=|V|`, such that `m_ij = 1`, if and only if `(i,j)` is in `E`.

  For example, the input string  '110101011' corresponds to

                                 u1  u2  u3
                 1 1 0           |\ / \ /|
            M =  1 0 1      G =  | X   X |
                 0 1 1           |/ \ / \|
                                 v1  v2  v3

* `fixed`: Count the number of binary matrices whose row and column sums match
   prescribed integers. 
   
   An `INSTANCE` has the form `r(,r)*;c(,c)*`, 
   where `r` and `c` are positive integers. While the i-th occurrence of `r`
   defines the sum of row i, the j-th occurrence of `c` is the sum of 
   column j.  
  
   For example, the instance `2,2,2;1,2,1,2` corresponds to the vector pair

            row sums:    (2,2,2)
            column sums: (1,2,1,2)


* `interval`:
  Count the number of binary matrices whose row and column sums lie in
  the intervals prescribed by the encoded `INSTANCE`. 
  
  The parameter
  `INSTANCE` is a string of the form `l-u(,l-u)*;l-u(,l-u)*`. 
  The semicolon
  separates the row sums from the column sums. The i-th pair on the left
  side of the semicolon prescribes a lower and upper bound on the sum of
  row i. In contrast, the j-th pair on the right side of the semicolon
  defines the sum of column j. For convenience, a pair `l-u` can be
  replaced by an integer `l` if `l=u`.
  
  For example, the string `1-2,2,2-3;0-2,0-1,1-1,1-3` corresponds to the

            lower bounds on row sums:    (1,2,2)
            upper bounds on row sums:    (2,2,3)
            lower bounds on column sums: (0,0,1,1)
            upper bounds on column sums: (2,1,1,3)
