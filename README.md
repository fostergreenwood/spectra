# Spectra of Cayley Graphs of Complex Reflection Groups
This repository contains code used to compute eigenvalues of the adjacency, distance, and codimension matrices associated with the Cayley graph of a complex reflection group.  The code accompanies the following paper:

Foster-Greenwood, B. and Kriloff, C.: Spectra of Cayley graphs of complex reflection groups. Journal of Algebraic Combinatorics **44**(1), 33--57 (2016) [journal](http://link.springer.com/article/10.1007/s10801-015-0652-8) | [arXiv](http://arxiv.org/abs/1502.07392)


## Background
Let G be a complex reflection group, and let T be the set of all reflections in G.  The left Cayley graph of G with respect to T has a vertex for each element of the group and an edge joining h to g if gh<sup>-1</sup> is a reflection.  Various matrices associated to the group and Cayley graph arise from a general construction:  Given a class function f on G, define a matrix with rows and columns indexed by the elements of G and (g,h)-entry the value f(gh<sup>-1</sup>).  

Specifically, the *adjacency matrix* of the Cayley graph arises by using the reflection indicator function 
```
g -> 1 if g is a reflection and g -> 0 otherwise,
``` 
the *distance matrix* of the Cayley graph arises by using the absolute reflection length function 
```
g -> min number of factors to express g as a product of reflections,
``` 
and the *codimension matrix* (with respect to the reflection representation) arises from the codimension function 
```
g -> codimension of the fixed point space of g.
``` 

The eigenvalues of the adjacency, distance, and codimension matrices can be calculated via a character theoretic formula, and if f is integer-valued and constant on rational conjugacy classes, then the eigenvalues are guaranteed to be integers.  (See the [article](http://arxiv.org/abs/1502.07392) for details and references.)

## Running the code

The code in this repository is written for use in GAP 3 with the CHEVIE package.  The CHEVIE package contains functionality for working with complex reflection groups and has not yet been ported to GAP 4.  For more information, see [GAP](http://www.gap-system.org) and the [GAP3 distribution](http://webusers.imj-prg.fr/~jean.michel/gap3/) prepared by Jean Michel.

The purpose of the code in this repository is to calculate eigenvalues of the adjacency, distance, and codimension matrices via the character formula and to determine if absolute reflection length is constant on rational conjugacy classes.  The files 
```
spectraprompt.gap
constantprompt.gap
```
prompt the user with how to call the relevant functions from the file
```
spectra.gap
```
for the specific reflection groups of interest.
