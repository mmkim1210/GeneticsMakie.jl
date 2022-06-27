# GeneticsMakie

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mmkim1210.github.io/GeneticsMakie.jl/dev)
[![Build Status](https://github.com/mmkim1210/GeneticsMakie.jl/workflows/CI/badge.svg)](https://github.com/mmkim1210/GeneticsMakie.jl/actions)
[![Coverage](https://codecov.io/gh/mmkim1210/GeneticsMakie.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mmkim1210/GeneticsMakie.jl)

## Installation
GeneticsMakie.jl supports Julia v1.6 or later.
```julia
julia> ]
pkg> add GeneticsMakie
```

## Examples
Visualize 1 Mb window around KMT2E:
<p align="center"><img width="70%" style="border-radius: 5px;" src="figs/KMT2E-gene.png"></p>

Highlight KMT2E and several nearby genes:
<p align="center"><img width="70%" style="border-radius: 5px;" src="figs/KMT2E-gene-highlight.png"></p>

Visualize KMT2E isoforms:
<p align="center"><img width="70%" style="border-radius: 5px;" src="figs/KMT2E-isoform.png"></p>

Visualize KMT2E and nearby isoforms simultaneously:
<p align="center"><img width="90%" style="border-radius: 5px;" src="figs/KMT2E-isoform-others.png"></p>

Visualize KMT2E isoforms w/ randomly generated expression
<p align="center"><img width="85%" style="border-radius: 5px;" src="figs/KMT2E-expression.png"></p>

Visualize LD for KMT2E locus:
<p align="center"><img width="80%" style="border-radius: 5px;" src="figs/KMT2E-LD.png"></p>

Visualize LD as a square for KMT2E locus:
<p align="center"><img width="60%" style="border-radius: 5px;" src="figs/KMT2E-LD-square.png"></p>

Visualize GWAS results for KMT2E locus:
<p align="center"><img width="70%" style="border-radius: 5px;" src="figs/KMT2E-locuszoom.png"></p>

Visualize GWAS results for KMT2E locus with the same reference SNPs for LD:
<p align="center"><img width="100%" style="border-radius: 5px;" src="figs/KMT2E-locuszoom-index.png"></p>

Visualize Manhattan plot:
<p align="center"><img width="75%" style="border-radius: 5px;" src="figs/manhattan.png"></p>

Visualize Miami/Hudson plot:
<p align="center"><img width="75%" style="border-radius: 5px;" src="figs/miami.png"></p>

Visualize QQ plot of P values:
<p align="center"><img width="80%" style="border-radius: 5px;" src="figs/QQ.png"></p>

## Further examples
Patterns of association and LD:
<p align="center"><img width="80%" style="border-radius: 5px;" src="figs/effect-LD.png"></p>
Colocalization of GWAS signals:
<p align="center"><img width="90%" style="border-radius: 5px;" src="figs/coloc.png"></p>
Association results across multiple phenotypes:
<p align="center"><img width="100%" style="border-radius: 5px;" src="figs/KMT2E-phewas2.png"></p>
MHC association for schizophrenia with increasing sample size:
<p align="center"><img width="80%" style="border-radius: 5px;" src="figs/C4A-locuszoom1.png"></p>
<p align="center"><img width="80%" style="border-radius: 5px;" src="figs/C4A-locuszoom2.png"></p>
LD structure for ~66,000 SNPs in MHC region:
<p align="center"><img width="60%" style="border-radius: 5px;" src="figs/MHC-LD-square.png"></p>