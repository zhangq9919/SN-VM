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


Input file (.lsn) format:
------------
m n k l NST
arcs
initial marking: p mu 

Please, consult LSN_format.txt files for detail

Input files can be obtained from *NDRtoSN.c

*NDRtoSN.c converts .ndr file into .lsn file and stores tables of names for places and transitions.

Output file (.txt) format:
------------
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

-rm: print raw matices and vectors

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

`HSN` file can be converted to `LSN` file using *Compiler and linker of Sleptsov net program.

Examples
------------
>NDRtoSN sn_add.ndr sn_add.lsn

-sn_add.ndr is the model for calculating z=x+y; converts sn_add.ndr file into sn_add.lsn and stores tables of names for places and transitions.

>SNVM sn_add.lsn add_result.txt

-Run sleptsov net for computing z=x+y

-Watch Sleptsov net behavior in add_result.txt

The more examples of `LSN` files and `HSN` files:

add2.hsn & add2.lsn: Add two numbers

matrix2.hsn & matrix2.lsn: Two-dimensional matrix multiplication

mul3.hsn & mul3.lsn: Multiply three numbers

pol2.hsn & pol2.lsn: Computational polynomial n=2

fdiv.lsn: Division

depnz3.lsn: Exact double exponent counters 2^2^k k=3

References:
------------
1. Zaitsev DA (2016) Sleptsov Nets Run Fast. IEEE Trans Syst Man Cybern Syst 46:682–693. https://doi.org/10.1109/TSMC.2015.2444414

2. Zaitsev D.A., Jürjens J. Programming in the Sleptsov net language for systems control, Advances in Mechanical Engineering, 2016, Vol. 8(4), 1-11. https://doi.org/10.1177%2F1687814016640159

3. Zaitsev D.A. Universal Sleptsov Net, International Journal of Computer Mathematics, 94(12) 2017, 2396-2408. http://dx.doi.org/10.1080/00207160.2017.1283410

4. Tatiana R. Shmeleva, Jan W. Owsiński, Abdulmalik Ahmad Lawan (2021) Deep learning on Sleptsov nets, International Journal of Parallel, Emergent and Distributed Systems, 36:6, 535-548, https://doi.org/10.1080/17445760.2021.1945055

5. Qing Zhang, Ding Liu, Yifan Hou, Sleptsov Net Processor, International Conference ”Problems of Infocommunications. Science and Technology” (PICST2022), 10-12 October, 2022, Kyiv, Ukraine.

6. Hongfei Zhao, Ding Liu, Yifan Hou, Compiler and Linker of Sleptsov Net Program,International Conference ”Problems of Infocommunications. Science and Technology” (PICST2022), 10-12 October, 2022, Kyiv, Ukraine.
