# asxp

### What asxp does

It runs under Linux (Ubuntu 16.04 LTS on my machine) and generates perspective views of algebraic surfaces via raycasting:

You define a certain polynomial f(x, y, z, a, b) in the variables x, y, z and parameters a and b in the configuration file of *asxp*. 


Then you can go into a visualizing window for the surface, where the position of the mouse determines the parameters a and b as a0 and b0 and the surface plotted is the (two dimensional projection of the) zero-set
```
f(x, y, z, a0, b0) = 0
```
in 3d-space with coordinates x, y, z.


### Look at the video

The current version runs with CUDA support and uses the parallel processing power of the graphic card:

[![asxp demo video](https://img.youtube.com/vi/hFiTgNpNDK8/0.jpg)](https://www.youtube.com/watch?v=hFiTgNpNDK8)

### Building asxp

#### Introduction

At the moment I give only a sketchy explanation for experts who know how to extract the necessary dependency informations from the makefiles and can change the makefiles accordingly. If you need more detailed information how to build asxp please write me an e-mail.

#### Steps

Suppose you want to build asxp in your directory ```software``` anywhere in your filesystem.

```
cd software
git clone https://www.github.com/juergenboehm/asxp_npr2.git
cd asxp_npr2
```

Modify paths and variable values in ```make.pro``` and in ```cuda/Makefile``` if necessary.

```
./qmakeit
cd cuda
make
cd ..
make
```
At the moment, the cuda library ```cuda/lib/libasxp.a``` must always be built separately.


### Author

Jürgen Böhm - see [www.aviduratas.de](http://www.aviduratas.de)





