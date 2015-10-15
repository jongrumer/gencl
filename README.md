
# G E N C L

A code for GEN-eration of C-onfiguration state function L-ists

* Modern fortran version of gencsl(f), based on latest ATSP2K version [3].
* Development version -- use with caution!

This program computes all possible couplings for each member of
the reference set, generates all unique, possible configurations
and their couplings  from the active set, and for each replacement,
generates configurations and their couplings, then outputs the
configurations and couplings with given final term(s).


### Version History (pre-git)

Authors                        | Year | Version description
------------------------------ | ---- | ----------------------------------------
C. Froese Fischer and Bin Liu  | 1983 | Original Version, published 1991 [1,2]
G. Gaigalas                    | 1997 | Added support for f-shells
J. Grumer                      | 2015 | Cleaned and re-implemented in modern Fortran

### Description (from [1])

#### Nature of problem:
In configuration interaction calculations as well as variational MCHF 
calculations, the wave function is expanded in terms of configuration state 
functions. In complex atoms, the configuration state lists are best generated 
systematically according to rules. The present program embodies several rules 
which have been found useful in dealing with the problem of correlation.

#### Solution method:
The notion of a complex [4] and single and double replacements are two 
important concepts that have evolved from the study of correlation in 
many-electron systems. By systematically generating all couplings of a set of 
electrons or by considering all possible replacements, this program generates 
configuration state lists in the "clist" format required by the MCHF atomic 
structure package [2].

#### Restrictions:
A maximum of 5 (five) subshells (in addition to the common closed shells) is 
allowed in any given configuration state and no more than 15 electrons in each 
of the different sets that are created. These restrictions may be removed by 
changing dimension statements and some format statements and, in the case of 
the active set, the code will have to be extended.

### I/O 

#### Input: (Interactive)
  1. Header
  2. List of closed shells
  3. The "reference" set
  4. The "active" set
  5. Replacements from the reference set
    * Virtual set if replacement is 's' or 'd' or 'sd'
  6. Final Term

#### Output:
  1.   Header
  2.   List of closed shells
  3.   Configurations      : FORMAT(8(1X,A3,'(',I2,')'))
  4.   and their couplings : FORMAT(15(1X,A3)) 

### References
* [1] C. F. Fischer and Bin Liu, Computer Physics Communications, 64,  406 (1991)
* [2] C. F. Fischer et al, Computer Physics Communications, 64,  369 (1991)
* [3] C. F. Fischer et al, Computer Physics Communications, 176, 559 (2007)
* [4] D. Layzer, Ann. Phys. 8, 271 (1959)

