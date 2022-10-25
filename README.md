# LowpassingHM
Improve your single cell gene expression heatmap with one function

## Example
```
LowPassingHM(seurat_object,beta = 10,markers = top10$gene) # the inputted seurat object need build snn graph at first
```

### Raw heatmap
<img src="https://github.com/WWXkenmo/LowpassingHM/blob/main/hm_raw.png" alt="raw" width="600" />

### New heatmap
<img src="https://github.com/WWXkenmo/LowpassingHM/blob/main/hm_new.png" alt="new" width="600" />
