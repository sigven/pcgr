# Third-party licenses

PCGR is distributed under the MIT License (see `LICENSE.md`). Parts of PCGR are
derived from third-party software, which is distributed under the licenses
reproduced below.

## scarHRD

The genomic instability ("HRD") scores in PCGR - HRD-LOH, large-scale state
transitions (LST) and telomeric allelic imbalance (TAI) - are computed by
`pcgr/hrd.py`, a Python port of the algorithms in the
[scarHRD](https://github.com/sztup/scarHRD) R package.

- Source: https://github.com/sztup/scarHRD
- Reference: Sztupinszki Z, Diossy M, Krzystanek M, et al. Migrating the SNP
  array-based homologous recombination deficiency measures to next generation
  sequencing data of breast cancer. npj Breast Cancer 2018;4:16.
  https://doi.org/10.1038/s41523-018-0066-6
- License: MIT

```
MIT License

Copyright (c) 2020 Zsofia Sztupinszki

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
