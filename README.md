# README

## Info

To compute (an approximation of) the top eignenvalues and associated eigenvectors of a large sparse matrix using the pwer iteration method: https://en.wikipedia.org/wiki/Power_iteration


## To compile

gcc topeigen.c -o topeigen -O9 -lm

## To execute

./topeigen net.txt k res.txt
- k = number of top eignenvectors/eigenvalues to compute
- net.txt should contain on each line: "i j value". The matrix has to be symetric (if "i j" is here then "j i" is not here)
- will print the result in res.txt: on column k the k highest eigenvalue (in absolute value) followed by the entries of its eigenvector.


## Initial contributors

Maximilien Danisch  
September 2017  
http://bit.ly/danisch  
maximilien.danisch@gmail.com
