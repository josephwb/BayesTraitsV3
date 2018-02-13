BayesTraits version 3.x
---------------
Forked from [Mark Pagel's Software page](http://www.evolution.rdg.ac.uk/BayesTraitsV3.0.1/BayesTraitsV3.0.1.html).

This is simply a convenient way for me to play with the code. If you want the canonical version please visit the website above.

*NOTE*: this version is not guaranteed to change in step with the original source, or to even execute correctly. Beware!

If you are brave enough to use this version, compilation is straightforward (on linux, anyway). 

### Dependencies
BayesTraits requires both the GNU Scientific Library (GSL) and the Linear Algebra Package (LAPACK). Install them with:

```
sudo apt-get install libgsl-dev gsl-bin 
sudo apt-get install libblas-dev liblapack-dev
```

### Compile
I made a super simple Makefile. To compile the serial version, just do:

```
make
```

To compile the threaded version, do:

```
make threaded
```

