#Plot.flowchart.At_Reseq.R

library(DiagrammeR)
grViz("digraph flowchart {
      # node definitions with substituted label text
      node [fontname = Helvetica]        
      tab1 [label = '@@1', shape = box, fontsize = 28]
      tab2 [label = '@@2', shape = box, fontsize = 28]
      tab3 [label = '@@3', shape = box, fontsize = 28]
      tab4 [label = '@@4', shape = box, fontsize = 28]
      tab5 [label = '@@5', shape = box, fontsize = 28]
      tab6 [label = '@@6', shape = box, fontsize = 28]
      tab7 [label = '@@7', shape = box, fontsize = 28]
      tab8 [label = '@@8', shape = box, fontsize = 28]
      tab9 [label = '@@9', shape = box, fontsize = 28]
      tab10 [label = '@@10', shape = box, fontsize = 28]
      tab11 [label = '@@11', shape = box, fontsize = 28]
      tab12 [label = '@@12', shape = box, fontsize = 28]
      tab13 [label = '@@13', shape = box, fontsize = 28]
      tab14 [label = '@@14', shape = diamond, fontsize = 20]
      tab15 [label = '@@15', shape = diamond, fontsize = 20]
      tab16 [label = '@@16', shape = diamond, fontsize = 20]
      tab17 [label = '@@17', shape = diamond, fontsize = 20]
      tab18 [label = '@@18', shape = diamond, fontsize = 20]
      tab19 [label = '@@19', shape = diamond, fontsize = 20]
      tab20 [label = '@@20', shape = diamond, fontsize = 20]
      tab21 [label = '@@21', shape = diamond, fontsize = 20]
      tab22 [label = '@@22', shape = diamond, fontsize = 20]
      tab23 [label = '@@23', shape = diamond, fontsize = 20]
      tab24 [label = '@@24', shape = diamond, fontsize = 20]
      tab25 [label = '@@25', shape = diamond, fontsize = 20]
      tab26 [label = '@@26', shape = diamond, fontsize = 20]

      # edge definitions with the node IDs
      tab1 -> tab2 -> tab3 -> tab4 -> tab5 -> tab6 -> tab7 -> tab8 -> tab9 -> tab10 -> tab12 -> tab13;
      tab4 -> tab11 -> tab12;
      tab14 -> tab15 -> tab16 -> tab17 -> tab18 -> tab19 -> tab20 -> tab21  -> tab22 -> tab23 -> tab25 -> tab26;
      tab17 -> tab24 -> tab25;      
      }

      [1]: 'Quality-trimmed Reads'
      [2]: 'Map Reads to Reference'
      [3]: 'Mark Duplicates'
      [4]: 'Base Recalibration'
      [5]: 'Call Variants'
      [6]: 'Consolidate GVCFs'
      [7]: 'Joint-Call Cohort'
      [8]: 'Splitting SNPs + Indels'
      [9]: 'Filter Variants'
      [10]: 'Mutation Identification'
      [11]: 'Call Indels using Pindel'
      [12]: 'Unify variants'
      [13]: 'Annotation'
      [14]: 'fastp'
      [15]: 'Pipe 1'
      [16]: 'Pipe 2'
      [17]: 'Pipe 3,4'      
      [18]: 'Pipe 5'
      [19]: 'Pipe 6'
      [20]: 'Pipe 7'
      [21]: 'Pipe 8'
      [22]: 'Pipe 9'
      [23]: 'Pipe 10'
      [24]: 'Pipe 11'
      [25]: 'Pipe 12'
      [26]: 'Pipe 13'
")
