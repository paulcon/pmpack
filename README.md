Parameterized Matrix Pack (PMPACK) Version 1.0
==============================================

In short, a MATLAB(R) package for working with parameterized marices.

Parameterized matrices or parametric matrices are those that depend analytically on one or more parameters.  PMPack is a suite of Matlab functions for applying polynomial spectral methods  (Galerkin and pseudospectral) to linear equations with parameterized matrrices.

Synposis
--------

### Simple parameterized matrix

	epsilon=0.5;
	A=@(s) [1+epsilon s; s 1];
	Ax=@(s,x) A(s)*x;
	b=@(s) [2; 1];
	iAb=@(s) A(s)\b(s);
	N=2;
	s=parameter();

	%% Standard pseudospectral
	p_order=15;
	X=pseudospectral(iAb,s,p_order);
	
### More examples to come	

History
-------

Current version: Pre 1.0

### Pre 1.0
Paul and I are working on this code to release with a SISC publication.

License
-------

Copyright 2009-2010, Paul G. Constantine and David F. Gleich. 
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY Paul G. Constantine and David F. Gleich ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Paul G. Constantine and David F. Gleich OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the authors and should not be interpreted as representing official policies, either expressed or implied, of Paul G. Constantine and David F. Gleich.

  
