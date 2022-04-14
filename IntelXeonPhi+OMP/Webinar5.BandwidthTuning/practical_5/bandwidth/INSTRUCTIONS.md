# Developing Understanding of High Bandwidth Memory (MCDRAM)


MCDRAM is One of the biggest conceptual differences which challenge
existing HPC programmers on the new Knights Landing architecture.

What makes MCDRAM unique compared to memory for normal CPUs is that
it resides on the Knights Landing chip itself separate from 
any existing RAM that may also be on the system. Another 
novelty of MCDRAM is the exceptionally high bandwidth
it offers applications, provided it is used correctly
(this is why the term "High Bandwidth Memory" is often
used).

Yet another difficulty in understanding MCDRAM comes
from exactly how programmers can use it. Since MCDRAM
is separate from RAM, the programmer must generally
take extra steps to use it effectively. What steps
the programmer must take depends on what "memory mode"
the Knights Landing chip is in. Knights landing 
can operate in one of three modes:

* Cache
* Flat
* Hybrid


In Cache mode the MCDRAM is used as a last-level cache and so
its use is completely transparent to the programmer (but
special programming strategies may be employed to improve
its performance).

In Flat mode the MCDRAM is enabled as a separate NUMA node
which the programmer must *explicitly* indicate to use in code.

Hybrid mode takes some fixed percentage of MCDRAM to act as 
cache and leaves the remaining to be explicitly programmable.


Generally speaking if an application can run in Flat mode it will
achieve the best bandwidth, but this requires extra programming effort.
On the other hand cache mode will still achieve respectable bandwidth with minimal 
programming effort and so it is often recommended that applications start in
cache mode and move into flat mode as an additional tuning step.


Here we will investigate the impact of MCDRAM both in cache and in flat mode
on a memory bound application. First to get a baseline we will run the 
STREAM benchmark authored by John D. McCalphin, Ph.D and is hosted at
the website: https://www.cs.virginia.edu/stream/. For simplicity I have
already copied the relevant files and made appropriate changes to the makefile.


## MCDRAM Mode: Cache

The TACC KNL machines are by default given in Cache mode, but it is
possible to request flat mode machines as well. We will begin
with the default.

### Running the Stream Benchmark

Go into the `stream` directory and type `make`. Next run
`sh runstream.sh` This will give you some numbers which we will
use as a baseline going forward. Here is an example output
from my own run:


```
-------------------------------------------------------------
Function    Best Rate MB/s  Avg time     Min time     Max time
Copy:          259320.3     0.004972     0.004936     0.005014
Scale:         261174.8     0.004958     0.004901     0.005013
Add:           286443.2     0.007568     0.006703     0.008066
Triad:         317312.1     0.006104     0.006051     0.006166
-------------------------------------------------------------
Solution Validates: avg error less than 1.000000e-13 on all three arrays
-------------------------------------------------------------
```

We will be especially interested in the vector add case, as we have 
used that example a lot throughout these exercises and lectures.



### Preparing Our Own Stream Benchmark

Here we will compute our own vector add and see what steps it
takes to achieve a decent fraction of peak bandwidth from MCDRAM.
As we are currently in Cache mode this will not require very much work yet,
but it will be a good point of comparison for when we move into flat mode.



#### Step 1 Parallelizing + Vectorizing the vector add example

At this point you should be fairly comfortable with OpenMP.
Find the vector add code in either main.cpp or main.f90
and add the appropriate OpenMP directives to parallelize
and vectorize this code. Don't forget to add
`num_threads(68)`.


This code should already contain timing logic and logic to compute
the effective bandwidth of this code, assuming that the cost
of arithmetic can be ignored. Run this code. It may not match
exactly with the STREAM benchmark result for a variety of reasons,
but we should be able to get at least 200 GB/s.


## MCDRAM mode: Flat
Now we will investigate flat mode as a way to improve performance of
bandwidth-limited code. Note that in order to request a flat-mode
KNL on TACC you must issue the following modified `idev` command:

```shell
idev -p Flat-Quadrant
```

### Running STREAM
In this mode we will actually need to run stream twice to get two bandwidth values,
one running purely with RAM and the other running purely with MCDRAM. This can be achieved
with `numactl` as follows:


Pure RAM run:
```
sh run_stream.sh
```
The peak bandwidth from this on flat mode should be around 90 or 100 GB/s.

Pure MCDRAM run:
```
numactl -m1 sh run_stream.sh
```

The peak bandwidth from this on flat mode should be around 400 GB/s.

### Writing our own STREAM
The files for this are located under the C or Fortran directories.

#### Step 1: Parallel+Vectorizing

With this step we can use the code already written for Cache mode and run it twice. Again as in the STREAM case we use numactl to 
indicate which memory we want.

Pure RAM:
```
./main
```

Pure MCDRAM:
```
numactl -m1 ./main
```


#### Step 2: Using Libmemkind

Using numactl -m1 to bind all memory allocations to MCDRAM will work fine if the required memory does not exceed its capacity (16GB), but
this is not a realistic assumption for many applications. The way to get the benefit of flat mode while mixing both MCDRAM and RAM allocations
is to use libmemkind. The way libmemkind is used differs on C++ and Fortran, so I will start by explaining how it is done in C++ and then
follow with Fortran.

For C++ the libmemkind library gives a special allocator to indicate which memory a C++ container should use. A short example is given below:

```c++
#include <vector>
#include <hbw_allocator.h>

std::vector<float,hbw::allocator<float> > x(500); // Allocate 500 floats on MCDRAM
std::vector<float                       > y(500); // Allocate 500 floats on normal RAM
```

It will also be necessary to link libmemkind. An example (minimal) build command is given below

```
icpc --std=c++11 main.cpp -o main -lmemkind
```


Since Fortran does not have the same notion of allocators as C++, memory must be targeted differently.


```
real, allocatable :: x (:)
real, allocatable :: y (:)
!dir$ attributes fastmem :: x

allocate (x(500)) !allocate 500 reals in MCDRAM
allocate (y(500)) !allocate 500 reals in normal RAM
```

again it is necessary to link libmemkind here, which I show by example:


```shell
ifort main.f90 -o main -lmemkind
```



In our vector add example, simply replace all major allocations (the three arrays) 
with allocations that specifically target MCDRAM with libmemkind. 
Run and record the result.


# Troubleshooting

You can tell whch memory mode you are in by running the following command:

```shell
numactl -H
```
Note that we have not discussed *cluster* modes yet, but that will also impact the output of numactl -H. For now we are assuming that 
all the KNL machines are in *quadrant* mode.

## Cache+Quadrant Mode Output of Numactl -H

```
available: 1 nodes (0)
node 0 cpus: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271
node 0 size: 98207 MB
node 0 free: 89768 MB
node distances:
node   0
  0:  10
```

## Flat+Quadrant Mode Output of Numactl -H
```
available: 2 nodes (0-1)
node 0 cpus: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271
node 0 size: 98207 MB
node 0 free: 90456 MB
node 1 cpus:
node 1 size: 16384 MB
node 1 free: 15758 MB
node distances:
node   0   1
  0:  10  31
  1:  31  10
```
Note the main change here is that now there are *two* NUMA memory nodes instead of just one.


