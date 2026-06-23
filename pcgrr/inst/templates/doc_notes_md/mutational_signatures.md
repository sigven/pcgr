The set of somatic mutations observed in a tumor reflects the varied mutational 
processes that have been active during its life history, providing insights into 
the routes taken to carcinogenesis. Exogenous mutagens, such as tobacco smoke 
and ultraviolet light, and endogenous processes, such as APOBEC enzymatic 
family functional activity or DNA mismatch repair deficiency, result in 
characteristic patterns of mutation. There is growing evidence that mutational 
signatures can explain therapeutic response [@Brady2022-hr; @Levatic2022-hk].

In PCGR, we apply the **MutationalPatterns** package [@Blokzijl2018-nc] to 
identify the contribution of _known mutational signatures_ in a single tumor 
sample. Specifically, we apply [MutationalPatterns](https://github.com/UMCUGenetics/MutationalPatterns) 
to optimally reconstruct the observed spectrum of mutations through a 
[reference collection of known mutational signatures]({pubmed_url}32025018). 
By default, we restrict the signatures in the reference collection to those 
already observed in the tumor type in question (i.e. from large-scale <i>de novo</i> 
signature extraction on ICGC tumor samples).
