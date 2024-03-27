# dna-melting-defect

This repository hosts the Matlab codes used in the paper
[Equilibrium melting probabilities of a DNA molecule with a defect: An exact solution of the Poland–Scheraga model](https://pubs.aip.org/aip/jcp/article/159/14/145102/2916016/Equilibrium-melting-probabilities-of-a-DNA)
on DNA melting with a defect basepair. 

fig2.m produces the melting map, i.e. the probability that a basepair is open as a function of the position of the basepair along the DNA strand. The theoretical predicition is compared to the numerical solution of the Poland-Scheraga algorithm. 

fig3.m produces the melting curves, i.e. the probability that a basepair is open as a function of temperature, for three different basepairs: the defect, a neighboring basepair, and a background basepair; and for three different values of the loop exponent: c<1, 1<c<2 and c>2. 
