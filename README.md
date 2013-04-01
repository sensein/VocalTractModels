The contents of this project include the Vocal Tract model simulator and
area function synthesizer developed by Dr. Shinji Maeda and refined for
use in Matlab by Satrajit Ghosh.

# VTCalcs and VTsynth for Matlab

VTCalcs is a vocal tract model written by Dr. Shinji Maeda. It was
originally available for DOS only. This Matlab version can be compiled
and run on any platform supported by Mathworks. In addition, a newer
version of the synthesizer is available (VTsynth). This synthesizer is
used in VTCalcs, but has to be compiled separately from the mex files in
VTcalcs.

The projects contain both Windows and Linux mex files (although a little
old at this point). The mex file from VTsynth has already been put at
the appropriate location inside the VTCalcs distribution. Also, the
files are a mixture of C and C++ code and the compiler has to be changed
to g++ on Linux to compile correctly. 

Please contact Satrajit Ghosh for any clarifications or questions.

# LICENSE

Licensed under the Apache License, Version 2.0 (the "License");

    http://www.apache.org/licenses/LICENSE-2.0

See LICENSE file for details

# Publications

  @incollection{
  year={1990},
  isbn={978-94-010-7414-8},
  booktitle={Speech Production and Speech Modelling},
  volume={55},
  series={NATO ASI Series},
  editor={Hardcastle, WilliamJ. and Marchal, Alain},
  doi={10.1007/978-94-009-2037-8_6},
  title={Compensatory Articulation During Speech: Evidence from the Analysis and Synthesis of Vocal-Tract Shapes Using an Articulatory Model},
  url={http://dx.doi.org/10.1007/978-94-009-2037-8_6},
  publisher={Springer Netherlands},
  author={Maeda, Shinji},
  pages={131-149},
  language={English}
  }


