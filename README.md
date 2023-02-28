# Description
Sleptsov net processor(SNP) - usual multicore software for running low-level Sleptsov nets.

Sleptsov net processor is implemented as a software interpreter of SN behavior, which can load the transitions, arcs and tokens of SNs. This processor implements the concurrent firing of transitions in multiple instances in a single step.

Command line formats:
------------
SNP -d debug_level -nth nth  <input_file > output_file


Input (file) format:
------------
.txt / .lsn

m n k      
initial marking  
p t w

Please, consult SN_format.txt files for detail

Input files can be obtained from NDRtoLSN.c

NDRtoLSN.c converts .ndr file into .lsn file and stores tables of names for places and transitions.

Output (file) format:
------------
.txt

Sleptsov net behavior

Parameters:
-----------
-d debug_level: Level of output information     

(a): 'debug_level = 0'-----only final marking

(b): 'debug_level = 1'-----matrices of incoming and outgoing arcs and final marking

(c): 'debug_level = 2'-----all details

-nth nth: number of thread

m: number of places

n: number of transitions

k: number of arcs

An example
------------
>NDRtoLSN sn_add.ndr sn_add.txt

-sn_add.ndr is the model for calculating z=x+y; converts sn_add.ndr file into an_add.lsn file; sn_add.txt.nmp and sn_add.txt.nmt stores tables of names for places and transitions.

>SNP -d 0 -nth 4  <sn_add.lsn/sn_add.txt > add_result.txt

-Run sleptsov net for computing z=x+y; debug_level is 0, only the output of the final marking; nth is 4, 4 threads.

-Wathch Sleptsov net behavior in add_result.txt

Tools to generate models:
------------
Tina Toolbox for analysis of Petri nets and Time Petri nets http://www.laas.fr/tina

References:
------------
Zaitsev DA (2016) Sleptsov Nets Run Fast. IEEE Trans Syst Man Cybern Syst 46:682–693. https://doi.org/10.1109/TSMC.2015.2444414
