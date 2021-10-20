# findSimilarity

Finding a similar signal in a longer signal is a common problem in signal processing. [Waveform Similarity Overlap Add (WSOLA)](https://github.com/flatmax/gtkiostream/blob/master/applications/WSOLA.C#L82) is an example of an algorithm which uses signal similarity to speed up and slow down audio without changing its pitch.

The brute force method to find a similar signal between a reference (s) and a longer waveform (y) is based on RMS value of the distance :
```matlab
Ns=size(s,1);
res=buffer(y,Ns,Ns-1,'nodelay');
res=s-res;
errM=rms(res);
[err,nI]=min(errM);
```
This brute force method is slow because it uses a huge amount of memory in buffering the signal with one sample overlap and the computational complexity of the RMS function.

Methods based on the DFT [for example Mueen's implementation](https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html) break the nonlinear distance measure into a linear component which is convolved using the DFT :
```
(s-y)^2 = s^2 - 2s*y + y^2
```
The s*y convolution is performed in the DFT domain which is a speed up for large vectors. However for smaller s and y vectors the computation is still not efficient enough.

The solution proposed here reduces the order of the computation by breaking it into two steps. The first step computes the distance measure on a smaller signal space looking for the possible location of the global minimum. The second iteration does a complete search of the global minimum's region looking for the exact result.
