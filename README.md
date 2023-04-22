## `SN-VM`: Sleptsov Net VM

# Sleptsov net Virtual Machine(SN-VM) 

Sleptsov net Virtual Machine(SN-VM) - usual multicore software for running low-level Sleptsov nets (LSN).

Sleptsov net virtual machine is implemented as a software interpreter of SN behavior, which can load the transitions, arcs and tokens of LSN. This processor implements the concurrent firing of transitions in multiple instances in a single step.

How to use `SN-VM-GPU` as a part of experimental `SNC IDE&VM`:
--------------------------------------------------------------

We list references to components in "Compatibility" section.

1) Use `Tina` `nd` as graphical editor and its labels with special syntax (section "Transition substitution label") to specify transition substitution of `HSN`.

2) Use `NDRtoSN` to convert `NDR` file of `Tina` into `HSN` or `LSN`. 

3) Use `HSNtoLSN` to compile and link HSN file and mentioned in it `LSN` files into a single `LSN` file.

4) Run `LSN` file on `SN-VM` or `SN-VM-GPU`.


Compatibility: 
-------------- 

`Tina`, `nd`, and `NDR` file format according to https://projects.laas.fr/tina/index.php

`NDRtoSN` and Transition substitution labels according to https://github.com/dazeorgacm/NDRtoSN

`SN-VM` and `LSN` file format according to https://github.com/zhangq9919/Sleptsov-net-processor

`HSNtoLSN` and `HSN` file format according to https://github.com/HfZhao1998/Compiler-and-Linker-of-Sleptsov-net-Program

`SN-VM-GPU` and `MSN` file format according to https://github.com/tishtri/SN-VM-GPU

Command line formats:
------------
SNVM lsn_file.lsn output_file.txt 


Input (file) format:
------------
.txt / .lsn

m n k l NST
arcs
initial marking: p mu 

Please, consult LSN_format.txt files for detail

Input files can be obtained from *NDRtoSN.c

*NDRtoSN.c converts .ndr file into .lsn file and stores tables of names for places and transitions.

Output (file) format:
------------
.txt

Sleptsov net behavior

Program options:
-----------
-d debug_level: Level of output information     

(a): 'debug_level = 0'-----only final marking

(b): 'debug_level = 1'-----matrices of incoming and outgoing arcs and final marking

(c): 'debug_level = 2'-----all details

-nth nth: number of thread

-pm pm: choose format for printing marking (0 - usual or 1 - sparse vector)

-smax step: specified steps

-wm: print raw matices and vectors

Parameters:
-----------

m: number of places

n: number of transitions

k: number of arcs

l: number of nonzero markings

NST: number of substitution transition

Generator
------------
gen_pol: generator of .hsn for computing polynomials  n=k

>gen_pol 2 pol2.hsn

k=2: a2x^2+a1x+a0

*HSN can be converted to LSN using *Compiler and linker of Sleptsov net program.

Examples
------------
>NDRtoSN sn_add.ndr sn_add.lsn

-sn_add.ndr is the model for calculating z=x+y; converts sn_add.ndr file into an_add.lsn file and stores tables of names for places and transitions.

>SNVM sn_add.lsn add_result.txt

-Run sleptsov net for computing z=x+y

-Watch Sleptsov net behavior in add_result.txt

Tools to generate models:
------------
Tina Toolbox for analysis of Petri nets and Time Petri nets http://www.laas.fr/tina

References:
------------
1. Zaitsev DA (2016) Sleptsov Nets Run Fast. IEEE Trans Syst Man Cybern Syst 46:682–693. https://doi.org/10.1109/TSMC.2015.2444414

2. Zaitsev D.A., Jürjens J. Programming in the Sleptsov net language for systems control, Advances in Mechanical Engineering, 2016, Vol. 8(4), 1-11. https://doi.org/10.1177%2F1687814016640159

3. Zaitsev D.A. Universal Sleptsov Net, International Journal of Computer Mathematics, 94(12) 2017, 2396-2408. http://dx.doi.org/10.1080/00207160.2017.1283410

4. Qing Zhang, Ding Liu, Yifan Hou, Sleptsov Net Processor, International Conference ”Problems of Infocommunications. Science and Technology” (PICST2022), 10-12 October, 2022, Kyiv, Ukraine.

5. Hongfei Zhao, Ding Liu, Yifan Hou, Compiler and Linker of Sleptsov Net Program,International Conference ”Problems of Infocommunications. Science and Technology” (PICST2022), 10-12 October, 2022, Kyiv, Ukraine.
