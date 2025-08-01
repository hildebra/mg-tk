# MG-TK

## Introduction 

MG-TK is a pipeline developed to 
- Assemble metagenomes, profile miTags, profile functions, profile taxonomy, build MAGs (metagenomic asssembled genomes) using a variety of approaches (MG-TK.pl)
- Build a gene catalog based on these assemblies and predicted genes, build abundance matrices from these and annotate the genes functionally (secScripts/geneCat.pl)
- Merge MAGs into MGS (metagenomic species), build inter- and intra-speciies phylogenies (secScripts/MGS.pl)
- Additional functionalities are available, such as building phylogenetic tree automatically for many genomes and calcualting population genetic statistics on these (secScripts/phylo/buildTree5.pl)

MG-TK is implemented in Perl, C++ and uses R and python scripts. 

Author: Falk Hildebrand <Falk.Hildebrand@gmail.com> 


## Detailed information

- [Installing MG-TK](helpers/documentation/install.md)
	- [I/O considerations](helpers/documentation/install.md#I/O-considerations)
	- [Known Issues](helpers/documentation/install.md#Known-issues)
- [MG-TK manual](helpers/documentation/manual.md)
	- [MG-TK.pl flags](helpers/documentation/manual.md#MG-TKpl-flags)
	- [Genecat.pl flags](helpers/documentation/manual.md#genecatpl-flags)
	- [MGS.pl flags](helpers/documentation/manual.md#mgspl-flags)
	- [buildTree5.pl flags](helpers/documentation/manual.md#buildtree5pl-flags)
- [Outputs](helpers/documentation/manual.md#outputs)
	- [Abundance matrices](helpers/documentation/manual.md#abundance-matrices)
	- [Gene function & MAG/MGS gene content](helpers/documentation/manual.md#gene-function--magmgs-gene-content)
- [Running MG-TK](helpers/documentation/usage.md)
	- [Temporary and output files](helpers/documentation/usage.md#temporary-and-output-files)
	- [Mapping file](helpers/documentation/usage.md#mapping-file-for-MG-TK)
	- [Additional usage scenarios](helpers/documentation/usage.md#additional-usage-scenarios)
		- [map2tar mode](helpers/documentation/usage.md#map2tar-mode)
		- [Building phylogenetic trees with MG-TK](helpers/documentation/usage.md#building-phylogenetic-trees-with-MG-TK)
- [Examples](helpers/documentation/examples.md)
	- [MG-TK metagenomic assembly and gene catalog](helpers/documentation/examples.md#MG-TK-metagenomic-assembly-and-gene-catalog)
	- [Assembly-independent MG-TK mode](helpers/documentation/examples.md#assembly-independent-MG-TK-mode)
	- [Hybrid assemblies](helpers/documentation/examples.md#hybrid-assemblies)
- [FAQ](helpers/documentation/FAQ.md)




## License, citations etc

<details>
  <summary> Expand section </summary>

### Used software

plenty.. please refer to helpers/install/\*.yml for software that is available on Conda.

Other software used that was adapted and/or developed specifically for MG-TK (all implemented in C++):
- [sdm](https://github.com/hildebra/sdm), [LCA](https://github.com/hildebra/LCA): for 1) read qual filtering 2) least common ancestor calculations in tax assignments. Both are part of our amplicons sequencing pipeline [LotuS2](https://lotus2.earlham.ac.uk/)
- [clusterMAGs](https://github.com/hildebra/clusterMAGs): cluster MAGs into MGS (metagenomic species) based on conserved marker genes
- [camopy2](https://github.com/hildebra/canopy2): canopy clustering of gene catalogues, a much more efficient implementation of the [original algorithm](http://www.nature.com/articles/nbt.2939)
- [MSAfix](https://github.com/hildebra/MSAfix): fixes frameshifts in MSA that sometimes occur, important for high senstive intraspecific phylogenies
- [rdCover](https://github.com/hildebra/rdCover): calculates read coverage of genes, genomes etc based on mappings.

### Citing MG-TK

**Please cite MG-TK with:**
- Assembly mode: Hildebrand, F. et al. Antibiotics-induced monodominance of a novel gut bacterial order. Gut 68, 1781–1790 (2019). 
- Strain mode: Hildebrand, F. et al. Dispersal strategies shape persistence and evolution of human gut bacteria. Cell Host & Microbe 29, 1167-1176.e9 (2021). 
- Assembly-independent mode: Bahram, M. et al. Metagenomic assessment of the global diversity and distribution of bacteria and fungi. Environmental Microbiology 23, 316–326 (2021).
- sdm, LCA: Özkurt, E. et al. Microbiome (2022).

Falk Hildebrand <Falk.Hildebrand@gmail.com>

### License

 Copyright (c) 2017-2024 Falk Hildebrand

 MG-TK is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 2 of the License, or
 (at your option) any later version.

 MG-TK is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 See the file LICENSE for more details.

 You should have received a copy of the GNU General Public License
 along with the source code.  If not, see <http://www.gnu.org/licenses/>.

</details>
