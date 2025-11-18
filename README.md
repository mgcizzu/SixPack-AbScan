# SixPack-AbScan
A tool to predict antibody cross-reactivity on non-target species

Large collections of novel antibodies are currently being generated as part of various initiatives including, but not limited to, the Human Protein Atlas. These extensive collections facilitate proteome-wide analysis of cell diversity and are normally targeted towards human tissues, although the antibodies are sometimes tested and validated in other biomedically-relevant mammal species (Mouse, Rat), or a bit less commonly against standard research model animals (Zebrafish, Drosophila).

Researchers working on less common model or non-model species who wish to analyse via immunostaining the distribution of a given protein (or a cell type marked by a given protein) need typically to invest considerable time and resources in the re-screening of commercial antibodies in their target species. This repo aims to make their life a bit easier.

This is a very simple workflow that allows, given a genome or a transcriptome, to pre-screen in silico the cross reactivity of existing epitope-mapped antibodies. 
The workflow starts by generating a 6-frame translation of the genome/transcriptome, producing by brute-force all the possible sets of protein fragments potentially encoded by the input sequence(s). 
The user also provides a list with one or more epitopes as aminoacid stretches, and the code will search for the epitope stretches in the predicted proteome from the previous step.
The tool finally reports the potential hits, also suggesting potential multiple reactivity (ie. similar proteins encoded by paralogue genes).

I don't think the idea behind this tool is particularly novel, I heard about this strategy a decade ago from my Post-doc advisor Michalis Averof (www.averof-lab.org), but I think it's a very old concept.

Now that spatial proteomics via multiplexed-IF is getting traction, I often get the question of how realistic it is to implement in non-model species. 
That's a hard question, but I think a good starting point is considering whether the available molecular tools will work at all.

The tool is intrinsically limited in a few ways:
It’s only as good as the coverage of the genome/transcriptome. 
It only checks for perfect matches (full cross-reactivity) and not for partial ones (to be extended in future releases).
It only checks for linear epitope matches (doesn’t take into account the protein conformation).



