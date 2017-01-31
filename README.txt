Original code by professor Andrew Scott, 2014.
Modified by Michael Hoang, 2016.

FindSIC is a number crunching algorithm based on SIC-POVMs, vectors of complex numbers that each describe mutually perpendicular lines - coordinate systems for large dimensions.

It was accepted that at least one SIC existed for every dimension greater than 1, but this was mere conjecture. FindSIC sets out to prove that, for each dimension, there is indeed at least one SIC.

It does so by selecting a dimension (referred to as 'd'), forming a vector of d complex numbers, then optimizing the numbers in accordance to the SIC function. After a few attempts, it will evaluate if the current set of complex numbers passes the SIC function within a threshhold (1e-6).

Should the current set fail, the results are replaced with a new random set of complex numbers, which are then optimized.

Should the current set pass, the results are saved as a strong argument that a SIC exists in dimension d. Then the algorithm restarts, working on dimension d+1.