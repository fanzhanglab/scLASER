
<!-- rnb-text-begin -->

---
title: "R Notebook"
output: html_notebook
---



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-output-begin eyJkYXRhIjoiXG48IS0tIHJuYi1zb3VyY2UtYmVnaW4gZXlKa1lYUmhJam9pWUdCZ2NseHVYRzVuWlhSM1pDZ3BYRzVjYm5ObGRIZGtLRndpUkc5M2JteHZZV1J6WENJcFhHNWdZR0FpZlE9PSAtLT5cblxuYGBgclxuXG5nZXR3ZCgpXG5cbnNldHdkKFxcRG93bmxvYWRzXFwpXG5gYGBcblxuPCEtLSBybmItc291cmNlLWVuZCAtLT5cbiJ9 -->

````

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuXG5nZXR3ZCgpXG5cbnNldHdkKFwiRG93bmxvYWRzXCIpXG5gYGAifQ== -->

```r

getwd()

setwd(\Downloads\)
```

<!-- rnb-source-end -->
````



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYGBgclxuXG5nZXR3ZCgpXG5cbnNldHdkKFxcRG93bmxvYWRzXFwpXG5gYGBcbmBgYCJ9 -->

```r
```r

getwd()

setwd(\Downloads\)

<!-- rnb-source-end -->


<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->






<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-output-begin eyJkYXRhIjoiXG48IS0tIHJuYi1zb3VyY2UtYmVnaW4gZXlKa1lYUmhJam9pWUdCZ2NseHViR2xpY21GeWVTaHpZMHhCVTBWU0tWeHVZR0JnSW4wPSAtLT5cblxuYGBgclxubGlicmFyeShzY0xBU0VSKVxuYGBgXG5cbjwhLS0gcm5iLXNvdXJjZS1lbmQgLS0+XG4ifQ== -->

````

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubGlicmFyeShzY0xBU0VSKVxuYGBgIn0= -->

```r
library(scLASER)
```

<!-- rnb-source-end -->
````



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYGBgclxubGlicmFyeShzY0xBU0VSKVxuYGBgXG5gYGAifQ== -->

```r
```r
library(scLASER)

<!-- rnb-source-end -->


<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-output-begin eyJkYXRhIjoiXG48IS0tIHJuYi1zb3VyY2UtYmVnaW4gZXlKa1lYUmhJam9pWUdCZ2NseHViMkpxSUQwZ2NtVmhaRkpFVXloY0luTmpURUZUUlZKZmIySnFMVEl1Y21SelhDSXBYRzVnWUdBaWZRPT0gLS0+XG5cbmBgYHJcbm9iaiA9IHJlYWRSRFMoXFxzY0xBU0VSX29iai0yLnJkc1xcKVxuYGBgXG5cbjwhLS0gcm5iLXNvdXJjZS1lbmQgLS0+XG4ifQ== -->

````

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxub2JqID0gcmVhZFJEUyhcInNjTEFTRVJfb2JqLTIucmRzXCIpXG5gYGAifQ== -->

```r
obj = readRDS(\scLASER_obj-2.rds\)
```

<!-- rnb-source-end -->
````



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYGBgclxub2JqID0gcmVhZFJEUyhcXHNjTEFTRVJfb2JqLTIucmRzXFwpXG5gYGBcbmBgYCJ9 -->

```r
```r
obj = readRDS(\scLASER_obj-2.rds\)

<!-- rnb-source-end -->


<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-output-begin eyJkYXRhIjoiXG48IS0tIHJuYi1zb3VyY2UtYmVnaW4gZXlKa1lYUmhJam9pWUdCZ2NseHViMkpxTVNBOUlITmpURUZUUlZJb2IySnFRRzFsZEdGa1lYUmhMRzlpYWtCd1kzTXBYRzVnWUdBaWZRPT0gLS0+XG5cbmBgYHJcbm9iajEgPSBzY0xBU0VSKG9iakBtZXRhZGF0YSxvYmpAcGNzKVxuYGBgXG5cbjwhLS0gcm5iLXNvdXJjZS1lbmQgLS0+XG4ifQ== -->

````

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxub2JqMSA9IHNjTEFTRVIob2JqQG1ldGFkYXRhLG9iakBwY3MpXG5gYGAifQ== -->

```r
obj1 = scLASER(obj@metadata,obj@pcs)
```

<!-- rnb-source-end -->
````



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYGBgclxub2JqMSA9IHNjTEFTRVIob2JqQG1ldGFkYXRhLG9iakBwY3MpXG5gYGBcbmBgYCJ9 -->

```r
```r
obj1 = scLASER(obj@metadata,obj@pcs)

<!-- rnb-source-end -->


<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->





<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-output-begin eyJkYXRhIjoiXG48IS0tIHJuYi1zb3VyY2UtYmVnaW4gZXlKa1lYUmhJam9pWUdCZ2NseHViMkpxTVNBOExTQnlkVzVmYUdGeWJXOXVlU2hjYmlBZ2IySnFNU3hjYmlBZ1ltRjBZMmhmWTI5c0lEMGdYQ0p6WVcxd2JHVmZhV1JjSWl4Y2JpQWdkR2hsZEdFZ1BTQXlMRnh1SUNCbGNITnBiRzl1TG1Oc2RYTjBaWElnUFNBdFNXNW1MRnh1SUNCbGNITnBiRzl1TG1oaGNtMXZibmtnUFNBdFNXNW1MRnh1SUNCdFlYZ3VhWFJsY2k1amJIVnpkR1Z5SUQwZ016QXNYRzRnSUcxaGVDNXBkR1Z5TG1oaGNtMXZibmtnUFNBeE1DeGNiaUFnY0d4dmRGOWpiMjUyWlhKblpXNWpaU0E5SUVaQlRGTkZYRzRwWEc1Z1lHQWlmUT09IC0tPlxuXG5gYGByXG5vYmoxIDwtIHJ1bl9oYXJtb255KFxuICBvYmoxLFxuICBiYXRjaF9jb2wgPSBcXHNhbXBsZV9pZFxcLFxuICB0aGV0YSA9IDIsXG4gIGVwc2lsb24uY2x1c3RlciA9IC1JbmYsXG4gIGVwc2lsb24uaGFybW9ueSA9IC1JbmYsXG4gIG1heC5pdGVyLmNsdXN0ZXIgPSAzMCxcbiAgbWF4Lml0ZXIuaGFybW9ueSA9IDEwLFxuICBwbG90X2NvbnZlcmdlbmNlID0gRkFMU0VcbilcbmBgYFxuXG48IS0tIHJuYi1zb3VyY2UtZW5kIC0tPlxuIn0= -->

````

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxub2JqMSA8LSBydW5faGFybW9ueShcbiAgb2JqMSxcbiAgYmF0Y2hfY29sID0gXCJzYW1wbGVfaWRcIixcbiAgdGhldGEgPSAyLFxuICBlcHNpbG9uLmNsdXN0ZXIgPSAtSW5mLFxuICBlcHNpbG9uLmhhcm1vbnkgPSAtSW5mLFxuICBtYXguaXRlci5jbHVzdGVyID0gMzAsXG4gIG1heC5pdGVyLmhhcm1vbnkgPSAxMCxcbiAgcGxvdF9jb252ZXJnZW5jZSA9IEZBTFNFXG4pXG5gYGAifQ== -->

```r
obj1 <- run_harmony(
  obj1,
  batch_col = \sample_id\,
  theta = 2,
  epsilon.cluster = -Inf,
  epsilon.harmony = -Inf,
  max.iter.cluster = 30,
  max.iter.harmony = 10,
  plot_convergence = FALSE
)
```

<!-- rnb-source-end -->
````



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYGBgclxub2JqMSA8LSBydW5faGFybW9ueShcbiAgb2JqMSxcbiAgYmF0Y2hfY29sID0gXFxzYW1wbGVfaWRcXCxcbiAgdGhldGEgPSAyLFxuICBlcHNpbG9uLmNsdXN0ZXIgPSAtSW5mLFxuICBlcHNpbG9uLmhhcm1vbnkgPSAtSW5mLFxuICBtYXguaXRlci5jbHVzdGVyID0gMzAsXG4gIG1heC5pdGVyLmhhcm1vbnkgPSAxMCxcbiAgcGxvdF9jb252ZXJnZW5jZSA9IEZBTFNFXG4pXG5gYGBcbmBgYCJ9 -->

```r
```r
obj1 <- run_harmony(
  obj1,
  batch_col = \sample_id\,
  theta = 2,
  epsilon.cluster = -Inf,
  epsilon.harmony = -Inf,
  max.iter.cluster = 30,
  max.iter.harmony = 10,
  plot_convergence = FALSE
)

<!-- rnb-source-end -->


<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->





<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuXG5cblxuXG5vYmoxQG1ldGFkYXRhJGRpc2Vhc2UgPSBvYmoxQG1ldGFkYXRhJGRpc2Vhc2UgJT4lIGFzLm51bWVyaWMoKVxub2JqMSA8LSBhc3NvY2lhdGlvbl9uYW1fTFYoXG4gIG9iaiA9IG9iajEsXG4gIHRlc3RfdmFyID0gXCJkaXNlYXNlXCIsXG4gIHNhbXBsZW1fa2V5ID0gXCJzYW1wbGVfaWRcIixcbiAgZ3JhcGhfdXNlID0gXCJSTkFfbm5cIixcbiAgYmF0Y2hlcyA9IFwiYmF0Y2hcIixcbiAgY292cyA9IGMoXCJhZ2VcIiwgXCJzZXhcIiksXG4gIE5udWxsID0gMTAwMCxcbiAgbG9jYWxfdGVzdCA9IFRSVUUsbl9wY3MgPSAyMCxcbiAgc2VlZCA9IDEyMzQsIHJldHVybl9uYW0gPSBULCBrZXkgPSBcIk5BTV9cIlxuKVxuXG5cbmBgYCJ9 -->

```r




obj1@metadata$disease = obj1@metadata$disease %>% as.numeric()
obj1 <- association_nam_LV(
  obj = obj1,
  test_var = "disease",
  samplem_key = "sample_id",
  graph_use = "RNA_nn",
  batches = "batch",
  covs = c("age", "sex"),
  Nnull = 1000,
  local_test = TRUE,n_pcs = 20,
  seed = 1234, return_nam = T, key = "NAM_"
)

```

<!-- rnb-source-end -->


<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->





assosiation function long analysis




<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuXG5cbmZvckFuYWx5c2lzIDwtIGRhdGEuZnJhbWUobWV0YSwgbmFtX3BjcylcbmZvckFuYWx5c2lzJHZpc2l0X2ZhY3RvciA8LSBwYXN0ZTAoXCJWXCIsIGZvckFuYWx5c2lzJHZpc2l0KVxub3V0Y29tZXMgPC0gY29sbmFtZXMobmFtX3BjcylcblxuIyAyIHRpbWUgcG9pbnQsIG5vdCBzZWxlY3RpbmdcbiMgZnVuY3Rpb24gZm9yIGFzc29jaWF0aW9uXG5nZXRMTU0gPC0gZnVuY3Rpb24ocGMpe1xuICB0bXAgPC0gZm9yQW5hbHlzaXNcbiAgdG1wJG5hbV9QQyA8LSB0bXBbLHdoaWNoKGNvbG5hbWVzKHRtcCkgPT0gcGMpXVxuXG4gIG1vZCA8LSBsbWUoZml4ZWQgPSBuYW1fUEMgfiBkaXNlYXNlKnZpc2l0ICsgYWdlICsgc2V4LCBcbiAgICAgICAgICAgICByYW5kb20gPSB+IDEgfCBzdWJqZWN0X2lkL3NhbXBsZV9pZCwgXG4gICAgICAgICAgICAgZGF0YSA9IHRtcCwgbWV0aG9kPVwiTUxcIilcbiAgbnVsbCA8LSBsbWUoZml4ZWQgPSBuYW1fUEMgfiBkaXNlYXNlICsgdmlzaXQgKyBhZ2UgKyBzZXgsIFxuICAgICAgICAgICAgIHJhbmRvbSA9IH4gMSB8IHN1YmplY3RfaWQvc2FtcGxlX2lkLCBcbiAgICAgICAgICAgICBkYXRhID0gdG1wLCBtZXRob2Q9XCJNTFwiKVxuICBcbiAgZnVsbF90VGFiIDwtIHN1bW1hcnkobW9kKSR0VGFibGVcbiAgbHJ0IDwtIGFub3ZhKG51bGwsIG1vZClcbiAgbHJ0X3B2YWwgPC0gbHJ0JGBwLXZhbHVlYFsyXVxuICB3YW50ID0gZGF0YS5mcmFtZShOQU1zY29yZSA9IHBjLCBjb2VmZmljaWVudCA9IHJvd25hbWVzKGZ1bGxfdFRhYiksIGZ1bGxfdFRhYiwgbHJ0X3B2YWwgPSBscnRfcHZhbClcbn1cblxuXG4jc2VsZWN0aW5nXG5cblxuXG5cblxuXG5gYGAifQ== -->

```r


forAnalysis <- data.frame(meta, nam_pcs)
forAnalysis$visit_factor <- paste0("V", forAnalysis$visit)
outcomes <- colnames(nam_pcs)

# 2 time point, not selecting
# function for association
getLMM <- function(pc){
  tmp <- forAnalysis
  tmp$nam_PC <- tmp[,which(colnames(tmp) == pc)]

  mod <- lme(fixed = nam_PC ~ disease*visit + age + sex, 
             random = ~ 1 | subject_id/sample_id, 
             data = tmp, method="ML")
  null <- lme(fixed = nam_PC ~ disease + visit + age + sex, 
             random = ~ 1 | subject_id/sample_id, 
             data = tmp, method="ML")
  
  full_tTab <- summary(mod)$tTable
  lrt <- anova(null, mod)
  lrt_pval <- lrt$`p-value`[2]
  want = data.frame(NAMscore = pc, coefficient = rownames(full_tTab), full_tTab, lrt_pval = lrt_pval)
}


#selecting

```

<!-- rnb-source-end -->


<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuXG5gYGAifQ== -->

```r

```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiV2FybmluZyBtZXNzYWdlOlxuUiBncmFwaGljcyBlbmdpbmUgdmVyc2lvbiAxNCBpcyBub3Qgc3VwcG9ydGVkIGJ5IHRoaXMgdmVyc2lvbiBvZiBSU3R1ZGlvLiBUaGUgUGxvdHMgdGFiIHdpbGwgYmUgZGlzYWJsZWQgdW50aWwgYSBuZXdlciB2ZXJzaW9uIG9mIFJTdHVkaW8gaXMgaW5zdGFsbGVkLiBcbiJ9 -->

```
Warning message:
R graphics engine version 14 is not supported by this version of RStudio. The Plots tab will be disabled until a newer version of RStudio is installed. 
```



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZGV2dG9vbHM6Omluc3RhbGxfZ2l0aHViKFwiZmFuemhhbmdsYWIvc2NMQVNFUlwiLGF1dGhfdG9rZW4gID0gXCJnaHBfdngxenFiQ21rMkphUmhnRGlhdnp5SnhHVFBVbnZsMmtGamVjXCIsZm9yY2UgPSBUKVxuYGBgIn0= -->

```r
devtools::install_github("fanzhanglab/scLASER",auth_token  = "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",force = T)
```

<!-- rnb-source-end -->

<!-- rnb-output-begin eyJkYXRhIjoiUmVnaXN0ZXJlZCBTMyBtZXRob2RzIG92ZXJ3cml0dGVuIGJ5ICdodG1sdG9vbHMnOlxuICBtZXRob2QgICAgICAgICAgICAgICBmcm9tICAgICAgICAgXG4gIHByaW50Lmh0bWwgICAgICAgICAgIHRvb2xzOnJzdHVkaW9cbiAgcHJpbnQuc2hpbnkudGFnICAgICAgdG9vbHM6cnN0dWRpb1xuICBwcmludC5zaGlueS50YWcubGlzdCB0b29sczpyc3R1ZGlvXG5SZWdpc3RlcmVkIFMzIG1ldGhvZCBvdmVyd3JpdHRlbiBieSAnaHRtbHdpZGdldHMnOlxuICBtZXRob2QgICAgICAgICAgIGZyb20gICAgICAgICBcbiAgcHJpbnQuaHRtbHdpZGdldCB0b29sczpyc3R1ZGlvXG5Eb3dubG9hZGluZyBHaXRIdWIgcmVwbyBmYW56aGFuZ2xhYi9zY0xBU0VSQEhFQURcbiJ9 -->

```
Registered S3 methods overwritten by 'htmltools':
  method               from         
  print.html           tools:rstudio
  print.shiny.tag      tools:rstudio
  print.shiny.tag.list tools:rstudio
Registered S3 method overwritten by 'htmlwidgets':
  method           from         
  print.htmlwidget tools:rstudio
Downloading GitHub repo fanzhanglab/scLASER@HEAD
```



<!-- rnb-output-end -->

<!-- rnb-output-begin eyJkYXRhIjoiVGhlc2UgcGFja2FnZXMgaGF2ZSBtb3JlIHJlY2VudCB2ZXJzaW9ucyBhdmFpbGFibGUuXG5JdCBpcyByZWNvbW1lbmRlZCB0byB1cGRhdGUgYWxsIG9mIHRoZW0uXG5XaGljaCB3b3VsZCB5b3UgbGlrZSB0byB1cGRhdGU/XG5cbiAgMTogQWxsICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFxuICAyOiBDUkFOIHBhY2thZ2VzIG9ubHkgICAgICAgICAgICAgICAgICAgICAgICAgICAgXG4gIDM6IE5vbmUgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBcbiAgNDogQkggICAgICAgICAgICgxLjgxLjAtMSAgIC0+IDEuOTAuMC0xICApIFtDUkFOXVxuICA1OiBSY3BwICAgICAgICAgKDEuMC4xMCAgICAgLT4gMS4xLjAgICAgICkgW0NSQU5dXG4gIDY6IHZjdHJzICAgICAgICAoMC42LjMgICAgICAtPiAwLjYuNSAgICAgKSBbQ1JBTl1cbiAgNzogc3RyaW5naSAgICAgICgxLjcuMTIgICAgIC0+IDEuOC43ICAgICApIFtDUkFOXVxuICA4OiBybGFuZyAgICAgICAgKDEuMS4xICAgICAgLT4gMS4xLjYgICAgICkgW0NSQU5dXG4gIDk6IG1hZ3JpdHRyICAgICAoMi4wLjMgICAgICAtPiAyLjAuNCAgICAgKSBbQ1JBTl1cbiAxMDogbGlmZWN5Y2xlICAgICgxLjAuMyAgICAgIC0+IDEuMC40ICAgICApIFtDUkFOXVxuIDExOiBnbHVlICAgICAgICAgKDEuNi4yICAgICAgLT4gMS44LjAgICAgICkgW0NSQU5dXG4gMTI6IGNsaSAgICAgICAgICAoMy42LjEgICAgICAtPiAzLjYuNSAgICAgKSBbQ1JBTl1cbiAxMzogc3RyaW5nciAgICAgICgxLjUuMCAgICAgIC0+IDEuNi4wICAgICApIFtDUkFOXVxuIDE0OiBwbHlyICAgICAgICAgKDEuOC44ICAgICAgLT4gMS44LjkgICAgICkgW0NSQU5dXG4gMTU6IHBhcmFsbGVsbHkgICAoMS4zNi4wICAgICAtPiAxLjQ2LjAgICAgKSBbQ1JBTl1cbiAxNjogbGlzdGVudiAgICAgICgwLjkuMCAgICAgIC0+IDAuMTAuMCAgICApIFtDUkFOXVxuIDE3OiBmYXJ2ZXIgICAgICAgKDIuMS4xICAgICAgLT4gMi4xLjIgICAgICkgW0NSQU5dXG4gMTg6IGRpZ2VzdCAgICAgICAoMC42LjMxICAgICAtPiAwLjYuMzkgICAgKSBbQ1JBTl1cbiAxOTogZ2xvYmFscyAgICAgICgwLjE2LjIgICAgIC0+IDAuMTguMCAgICApIFtDUkFOXVxuIDIwOiBmdXR1cmUgICAgICAgKDEuMzMuMCAgICAgLT4gMS42OC4wICAgICkgW0NSQU5dXG4gMjE6IHByb2dyZXNzciAgICAoMC4xNC4wICAgICAtPiAwLjE4LjAgICAgKSBbQ1JBTl1cbiAyMjogZnV0dXJlLmFwcGx5ICgxLjExLjAgICAgIC0+IDEuMjAuMSAgICApIFtDUkFOXVxuIDIzOiBTNyAgICAgICAgICAgKDAuMi4wICAgICAgLT4gMC4yLjEgICAgICkgW0NSQU5dXG4gMjQ6IGlzb2JhbmQgICAgICAoMC4yLjcgICAgICAtPiAwLjMuMCAgICAgKSBbQ1JBTl1cbiAyNTogc2hhcGUgICAgICAgICgxLjQuNiAgICAgIC0+IDEuNC42LjEgICApIFtDUkFOXVxuIDI2OiBsYXZhICAgICAgICAgKDEuNy4yLjEgICAgLT4gMS44LjIgICAgICkgW0NSQU5dXG4gMjc6IGdncGxvdDIgICAgICAoNC4wLjAgICAgICAtPiA0LjAuMSAgICAgKSBbQ1JBTl1cbiAyODogZGF0YS50YWJsZSAgICgxLjE0LjggICAgIC0+IDEuMTcuOCAgICApIFtDUkFOXVxuIDI5OiB1dGY4ICAgICAgICAgKDEuMi4zICAgICAgLT4gMS4yLjYgICAgICkgW0NSQU5dXG4gMzA6IHRpbWVjaGFuZ2UgICAoMC4yLjAgICAgICAtPiAwLjMuMCAgICAgKSBbQ1JBTl1cbiAzMTogcHJvZGxpbSAgICAgICgyMDIzLjAzLjMxIC0+IDIwMjUuMDQuMjgpIFtDUkFOXVxuIDMyOiBSNiAgICAgICAgICAgKDIuNS4xICAgICAgLT4gMi42LjEgICAgICkgW0NSQU5dXG4gMzM6IHBpbGxhciAgICAgICAoMS45LjAgICAgICAtPiAxLjExLjEgICAgKSBbQ1JBTl1cbiAzNDogY3BwMTEgICAgICAgICgwLjQuNiAgICAgIC0+IDAuNS4yICAgICApIFtDUkFOXVxuIDM1OiB0emRiICAgICAgICAgKDAuMy4wICAgICAgLT4gMC41LjAgICAgICkgW0NSQU5dXG4gMzY6IHdpdGhyICAgICAgICAoMi41LjAgICAgICAtPiAzLjAuMiAgICAgKSBbQ1JBTl1cbiAzNzogdGltZURhdGUgICAgICg0MDIyLjEwOCAgIC0+IDQwNTEuMTExICApIFtDUkFOXVxuIDM4OiB0aWR5c2VsZWN0ICAgKDEuMi4wICAgICAgLT4gMS4yLjEgICAgICkgW0NSQU5dXG4gMzk6IHRpZHlyICAgICAgICAoMS4zLjAgICAgICAtPiAxLjMuMSAgICAgKSBbQ1JBTl1cbiA0MDogdGliYmxlICAgICAgICgzLjIuMSAgICAgIC0+IDMuMy4wICAgICApIFtDUkFOXVxuIDQxOiBzcGFyc2V2Y3RycyAgKDAuMy40ICAgICAgLT4gMC4zLjUgICAgICkgW0NSQU5dXG4gNDI6IHB1cnJyICAgICAgICAoMS4wLjIgICAgICAtPiAxLjIuMCAgICAgKSBbQ1JBTl1cbiA0MzogbHVicmlkYXRlICAgICgxLjkuMiAgICAgIC0+IDEuOS40ICAgICApIFtDUkFOXVxuIDQ0OiBpcHJlZCAgICAgICAgKDAuOS0xNCAgICAgLT4gMC45LTE1ICAgICkgW0NSQU5dXG4gNDU6IGhhcmRoYXQgICAgICAoMS4zLjAgICAgICAtPiAxLjQuMiAgICAgKSBbQ1JBTl1cbiA0NjogZ293ZXIgICAgICAgICgxLjAuMSAgICAgIC0+IDEuMC4yICAgICApIFtDUkFOXVxuIDQ3OiBnZW5lcmljcyAgICAgKDAuMS4zICAgICAgLT4gMC4xLjQgICAgICkgW0NSQU5dXG4gNDg6IGNsb2NrICAgICAgICAoMC42LjEgICAgICAtPiAwLjcuMyAgICAgKSBbQ1JBTl1cbiA0OTogZHBseXIgICAgICAgICgxLjEuMiAgICAgIC0+IDEuMS40ICAgICApIFtDUkFOXVxuIDUwOiBwcm94eSAgICAgICAgKDAuNC0yNyAgICAgLT4gMC40LTI4ICAgICkgW0NSQU5dXG4gNTE6IGJhY2twb3J0cyAgICAoMS40LjEgICAgICAtPiAxLjUuMCAgICAgKSBbQ1JBTl1cbiA1MjogUmNwcEVpZ2VuICAgICgwLjMuMy45LjMgIC0+IDAuMy40LjAuMiApIFtDUkFOXVxuIDUzOiBkcXJuZyAgICAgICAgKDAuMy4wICAgICAgLT4gMC40LjEgICAgICkgW0NSQU5dXG4gNTQ6IFJTcGVjdHJhICAgICAoMC4xNi0xICAgICAtPiAwLjE2LTIgICAgKSBbQ1JBTl1cbiA1NTogUmNwcEFubm95ICAgICgwLjAuMjAgICAgIC0+IDAuMC4yMiAgICApIFtDUkFOXVxuIDU2OiBGTk4gICAgICAgICAgKDEuMS4zLjIgICAgLT4gMS4xLjQuMSAgICkgW0NSQU5dXG4gNTc6IHBvbHljbGlwICAgICAoMS4xMC00ICAgICAtPiAxLjEwLTcgICAgKSBbQ1JBTl1cbiA1ODogZGVsZGlyICAgICAgICgxLjAtNiAgICAgIC0+IDIuMC00ICAgICApIFtDUkFOXVxuIDU5OiBzcGF0c3RhdC4uLi4gKDMuMC0yICAgICAgLT4gMy4yLTAgICAgICkgW0NSQU5dXG4gNjA6IHNwYXRzdGF0Li4uLiAoMy4wLTEgICAgICAtPiAzLjEtOSAgICAgKSBbQ1JBTl1cbiA2MTogdGVuc29yICAgICAgICgxLjUgICAgICAgIC0+IDEuNS4xICAgICApIFtDUkFOXVxuIDYyOiBhYmluZCAgICAgICAgKDEuNC01ICAgICAgLT4gMS40LTggICAgICkgW0NSQU5dXG4gNjM6IHNwYXRzdGF0Li4uLiAoMy4wLTEgICAgICAtPiAzLjEtMCAgICAgKSBbQ1JBTl1cbiA2NDogc3BhdHN0YXQuLi4uICgzLjEtNCAgICAgIC0+IDMuNC0zICAgICApIFtDUkFOXVxuIDY1OiBzcGF0c3RhdC4uLi4gKDMuMS0wICAgICAgLT4gMy42LTEgICAgICkgW0NSQU5dXG4gNjY6IGZzICAgICAgICAgICAoMS41LjIgICAgICAtPiAxLjYuNiAgICAgKSBbQ1JBTl1cbiA2Nzogc2FzcyAgICAgICAgICgwLjQuNSAgICAgIC0+IDAuNC4xMCAgICApIFtDUkFOXVxuIDY4OiBwcm9taXNlcyAgICAgKDEuMi4wLjEgICAgLT4gMS41LjAgICAgICkgW0NSQU5dXG4gNjk6IG1pbWUgICAgICAgICAoMC4xMiAgICAgICAtPiAwLjEzICAgICAgKSBbQ1JBTl1cbiA3MDogbGF0ZXIgICAgICAgICgxLjMuMCAgICAgIC0+IDEuNC40ICAgICApIFtDUkFOXVxuIDcxOiBqc29ubGl0ZSAgICAgKDEuOC40ICAgICAgLT4gMi4wLjAgICAgICkgW0NSQU5dXG4gNzI6IGh0dHB1diAgICAgICAoMS42LjkgICAgICAtPiAxLjYuMTYgICAgKSBbQ1JBTl1cbiA3MzogaHRtbHRvb2xzICAgICgwLjUuNSAgICAgIC0+IDAuNS45ICAgICApIFtDUkFOXVxuIDc0OiBmb250YXdlc29tZSAgKDAuNS4yICAgICAgLT4gMC41LjMgICAgICkgW0NSQU5dXG4gNzU6IGZhc3RtYXAgICAgICAoMS4xLjEgICAgICAtPiAxLjIuMCAgICAgKSBbQ1JBTl1cbiA3NjogY29tbW9ubWFyayAgICgxLjkuMCAgICAgIC0+IDIuMC4wICAgICApIFtDUkFOXVxuIDc3OiBjYWNoZW0gICAgICAgKDEuMC43ICAgICAgLT4gMS4xLjAgICAgICkgW0NSQU5dXG4gNzg6IGJzbGliICAgICAgICAoMC41LjEgICAgICAtPiAwLjkuMCAgICAgKSBbQ1JBTl1cbiA3OTogUmNwcEFybWFkLi4uICgwLjEyLjIuMC4wIC0+IDE1LjIuMi0xICApIFtDUkFOXVxuIDgwOiBtYXRyaXhTdGF0cyAgKDAuNjMuMCAgICAgLT4gMS41LjAgICAgICkgW0NSQU5dXG4gODE6IHJlc2hhcGUyICAgICAoMS40LjQgICAgICAtPiAxLjQuNSAgICAgKSBbQ1JBTl1cbiA4MjogYml0b3BzICAgICAgICgxLjAtNyAgICAgIC0+IDEuMC05ICAgICApIFtDUkFOXVxuIDgzOiBjYVRvb2xzICAgICAgKDEuMTguMiAgICAgLT4gMS4xOC4zICAgICkgW0NSQU5dXG4gODQ6IGd0b29scyAgICAgICAoMy45LjQgICAgICAtPiAzLjkuNSAgICAgKSBbQ1JBTl1cbiA4NTogZ3Bsb3RzICAgICAgICgzLjEuMyAgICAgIC0+IDMuMy4wICAgICApIFtDUkFOXVxuIDg2OiBycHJvanJvb3QgICAgKDIuMC4zICAgICAgLT4gMi4xLjEgICAgICkgW0NSQU5dXG4gODc6IGhlcmUgICAgICAgICAoMS4wLjEgICAgICAtPiAxLjAuMiAgICAgKSBbQ1JBTl1cbiA4ODogUmNwcFRPTUwgICAgICgwLjIuMiAgICAgIC0+IDAuMi4zICAgICApIFtDUkFOXVxuIDg5OiBzeXMgICAgICAgICAgKDMuNC4xICAgICAgLT4gMy40LjMgICAgICkgW0NSQU5dXG4gOTA6IHRpbnl0ZXggICAgICAoMC40NiAgICAgICAtPiAwLjU4ICAgICAgKSBbQ1JBTl1cbiA5MTogYXNrcGFzcyAgICAgICgxLjEgICAgICAgIC0+IDEuMi4xICAgICApIFtDUkFOXVxuIDkyOiB4ZnVuICAgICAgICAgKDAuMzkgICAgICAgLT4gMC41NCAgICAgICkgW0NSQU5dXG4gOTM6IGhpZ2hyICAgICAgICAoMC4xMCAgICAgICAtPiAwLjExICAgICAgKSBbQ1JBTl1cbiA5NDogZXZhbHVhdGUgICAgICgwLjIxICAgICAgIC0+IDEuMC41ICAgICApIFtDUkFOXVxuIDk1OiBvcGVuc3NsICAgICAgKDIuMC42ICAgICAgLT4gMi4zLjQgICAgICkgW0NSQU5dXG4gOTY6IGN1cmwgICAgICAgICAoNS4wLjAgICAgICAtPiA3LjAuMCAgICAgKSBbQ1JBTl1cbiA5NzogeWFtbCAgICAgICAgICgyLjMuNyAgICAgIC0+IDIuMy4xMiAgICApIFtDUkFOXVxuIDk4OiBybWFya2Rvd24gICAgKDIuMjUgICAgICAgLT4gMi4zMCAgICAgICkgW0NSQU5dXG4gOTk6IGtuaXRyICAgICAgICAoMS40NCAgICAgICAtPiAxLjUwICAgICAgKSBbQ1JBTl1cbjEwMDogY3Jvc3N0YWxrICAgICgxLjIuMCAgICAgIC0+IDEuMi4yICAgICApIFtDUkFOXVxuMTAxOiBodG1sd2lkZ2V0cyAgKDEuNi4yICAgICAgLT4gMS42LjQgICAgICkgW0NSQU5dXG4xMDI6IHNoaW55ICAgICAgICAoMS43LjUgICAgICAtPiAxLjEyLjEgICAgKSBbQ1JBTl1cbjEwMzogem9vICAgICAgICAgICgxLjgtMTIgICAgIC0+IDEuOC0xNCAgICApIFtDUkFOXVxuMTA0OiBkb3RDYWxsNjQgICAgKDEuMC0yICAgICAgLT4gMS4yICAgICAgICkgW0NSQU5dXG4xMDU6IHNwYW0gICAgICAgICAoMi45LTEgICAgICAtPiAyLjExLTEgICAgKSBbQ1JBTl1cbjEwNjogc3AgICAgICAgICAgICgxLjYtMCAgICAgIC0+IDIuMi0wICAgICApIFtDUkFOXVxuMTA3OiByYmlidXRpbHMgICAgKDIuMi4xMyAgICAgLT4gMi40ICAgICAgICkgW0NSQU5dXG4xMDg6IHJlY2lwZXMgICAgICAoMS4wLjggICAgICAtPiAxLjMuMSAgICAgKSBbQ1JBTl1cbjEwOTogcFJPQyAgICAgICAgICgxLjE4LjAgICAgIC0+IDEuMTkuMC4xICApIFtDUkFOXVxuMTEwOiBlMTA3MSAgICAgICAgKDEuNy0xMyAgICAgLT4gMS43LTE2ICAgICkgW0NSQU5dXG4xMTE6IGZvcmNhdHMgICAgICAoMS4wLjAgICAgICAtPiAxLjAuMSAgICAgKSBbQ1JBTl1cbjExMjogY29kYSAgICAgICAgICgwLjE5LTQgICAgIC0+IDAuMTktNC4xICApIFtDUkFOXVxuMTEzOiBicm9vbSAgICAgICAgKDEuMC41ICAgICAgLT4gMS4wLjExICAgICkgW0NSQU5dXG4xMTQ6IHV3b3QgICAgICAgICAoMC4xLjE0ICAgICAtPiAwLjIuNCAgICAgKSBbQ1JBTl1cbjExNTogc3BhdHN0YXQuLi4uICgzLjEtMCAgICAgIC0+IDMuNi0wICAgICApIFtDUkFOXVxuMTE2OiBzY3RyYW5zZm9ybSAgKDAuNC4xICAgICAgLT4gMC40LjIgICAgICkgW0NSQU5dXG4xMTc6IFJ0c25lICAgICAgICAoMC4xNiAgICAgICAtPiAwLjE3ICAgICAgKSBbQ1JBTl1cbjExODogcmV0aWN1bGF0ZSAgICgxLjI4ICAgICAgIC0+IDEuNDQuMSAgICApIFtDUkFOXVxuMTE5OiBSY3BwSE5TVyAgICAgKDAuNC4xICAgICAgLT4gMC42LjAgICAgICkgW0NSQU5dXG4xMjA6IFJBTk4gICAgICAgICAoMi42LjEgICAgICAtPiAyLjYuMiAgICAgKSBbQ1JBTl1cbjEyMTogcGxvdGx5ICAgICAgICg0LjEwLjIgICAgIC0+IDQuMTEuMCAgICApIFtDUkFOXVxuMTIyOiBwYmFwcGx5ICAgICAgKDEuNy0yICAgICAgLT4gMS43LTQgICAgICkgW0NSQU5dXG4xMjM6IHBhdGNod29yayAgICAoMS4xLjMgICAgICAtPiAxLjMuMiAgICAgKSBbQ1JBTl1cbjEyNDogbWluaVVJICAgICAgICgwLjEuMS4xICAgIC0+IDAuMS4yICAgICApIFtDUkFOXVxuMTI1OiBpZ3JhcGggICAgICAgKDEuNC4yICAgICAgLT4gMi4yLjEgICAgICkgW0NSQU5dXG4xMjY6IGdncmlkZ2VzICAgICAoMC41LjQgICAgICAtPiAwLjUuNyAgICAgKSBbQ1JBTl1cbjEyNzogZ2dyZXBlbCAgICAgICgwLjkuMyAgICAgIC0+IDAuOS42ICAgICApIFtDUkFOXVxuMTI4OiBmaXRkaXN0cnBsdXMgKDEuMS0xMSAgICAgLT4gMS4yLTQgICAgICkgW0NSQU5dXG4xMjk6IGZhc3REdW1taWVzICAoMS43LjMgICAgICAtPiAxLjcuNSAgICAgKSBbQ1JBTl1cbjEzMDogU2V1cmF0T2JqZWN0ICg0LjEuMyAgICAgIC0+IDUuMy4wICAgICApIFtDUkFOXVxuMTMxOiBubG9wdHIgICAgICAgKDIuMC4zICAgICAgLT4gMi4yLjEgICAgICkgW0NSQU5dXG4xMzI6IG1pbnFhICAgICAgICAoMS4yLjUgICAgICAtPiAxLjIuOCAgICAgKSBbQ1JBTl1cbjEzMzogaGFybW9ueSAgICAgICgwLjEuMSAgICAgIC0+IDEuMi40ICAgICApIFtDUkFOXVxuMTM0OiBjYXJldCAgICAgICAgKDYuMC05NCAgICAgLT4gNy4wLTEgICAgICkgW0NSQU5dXG4xMzU6IFNldXJhdCAgICAgICAoNC4zLjAgICAgICAtPiA1LjQuMCAgICAgKSBbQ1JBTl1cbjEzNjogbG1lNCAgICAgICAgICgxLjEtMzIgICAgIC0+IDEuMS0zOCAgICApIFtDUkFOXVxuIn0= -->

```
These packages have more recent versions available.
It is recommended to update all of them.
Which would you like to update?

  1: All                                           
  2: CRAN packages only                            
  3: None                                          
  4: BH           (1.81.0-1   -> 1.90.0-1  ) [CRAN]
  5: Rcpp         (1.0.10     -> 1.1.0     ) [CRAN]
  6: vctrs        (0.6.3      -> 0.6.5     ) [CRAN]
  7: stringi      (1.7.12     -> 1.8.7     ) [CRAN]
  8: rlang        (1.1.1      -> 1.1.6     ) [CRAN]
  9: magrittr     (2.0.3      -> 2.0.4     ) [CRAN]
 10: lifecycle    (1.0.3      -> 1.0.4     ) [CRAN]
 11: glue         (1.6.2      -> 1.8.0     ) [CRAN]
 12: cli          (3.6.1      -> 3.6.5     ) [CRAN]
 13: stringr      (1.5.0      -> 1.6.0     ) [CRAN]
 14: plyr         (1.8.8      -> 1.8.9     ) [CRAN]
 15: parallelly   (1.36.0     -> 1.46.0    ) [CRAN]
 16: listenv      (0.9.0      -> 0.10.0    ) [CRAN]
 17: farver       (2.1.1      -> 2.1.2     ) [CRAN]
 18: digest       (0.6.31     -> 0.6.39    ) [CRAN]
 19: globals      (0.16.2     -> 0.18.0    ) [CRAN]
 20: future       (1.33.0     -> 1.68.0    ) [CRAN]
 21: progressr    (0.14.0     -> 0.18.0    ) [CRAN]
 22: future.apply (1.11.0     -> 1.20.1    ) [CRAN]
 23: S7           (0.2.0      -> 0.2.1     ) [CRAN]
 24: isoband      (0.2.7      -> 0.3.0     ) [CRAN]
 25: shape        (1.4.6      -> 1.4.6.1   ) [CRAN]
 26: lava         (1.7.2.1    -> 1.8.2     ) [CRAN]
 27: ggplot2      (4.0.0      -> 4.0.1     ) [CRAN]
 28: data.table   (1.14.8     -> 1.17.8    ) [CRAN]
 29: utf8         (1.2.3      -> 1.2.6     ) [CRAN]
 30: timechange   (0.2.0      -> 0.3.0     ) [CRAN]
 31: prodlim      (2023.03.31 -> 2025.04.28) [CRAN]
 32: R6           (2.5.1      -> 2.6.1     ) [CRAN]
 33: pillar       (1.9.0      -> 1.11.1    ) [CRAN]
 34: cpp11        (0.4.6      -> 0.5.2     ) [CRAN]
 35: tzdb         (0.3.0      -> 0.5.0     ) [CRAN]
 36: withr        (2.5.0      -> 3.0.2     ) [CRAN]
 37: timeDate     (4022.108   -> 4051.111  ) [CRAN]
 38: tidyselect   (1.2.0      -> 1.2.1     ) [CRAN]
 39: tidyr        (1.3.0      -> 1.3.1     ) [CRAN]
 40: tibble       (3.2.1      -> 3.3.0     ) [CRAN]
 41: sparsevctrs  (0.3.4      -> 0.3.5     ) [CRAN]
 42: purrr        (1.0.2      -> 1.2.0     ) [CRAN]
 43: lubridate    (1.9.2      -> 1.9.4     ) [CRAN]
 44: ipred        (0.9-14     -> 0.9-15    ) [CRAN]
 45: hardhat      (1.3.0      -> 1.4.2     ) [CRAN]
 46: gower        (1.0.1      -> 1.0.2     ) [CRAN]
 47: generics     (0.1.3      -> 0.1.4     ) [CRAN]
 48: clock        (0.6.1      -> 0.7.3     ) [CRAN]
 49: dplyr        (1.1.2      -> 1.1.4     ) [CRAN]
 50: proxy        (0.4-27     -> 0.4-28    ) [CRAN]
 51: backports    (1.4.1      -> 1.5.0     ) [CRAN]
 52: RcppEigen    (0.3.3.9.3  -> 0.3.4.0.2 ) [CRAN]
 53: dqrng        (0.3.0      -> 0.4.1     ) [CRAN]
 54: RSpectra     (0.16-1     -> 0.16-2    ) [CRAN]
 55: RcppAnnoy    (0.0.20     -> 0.0.22    ) [CRAN]
 56: FNN          (1.1.3.2    -> 1.1.4.1   ) [CRAN]
 57: polyclip     (1.10-4     -> 1.10-7    ) [CRAN]
 58: deldir       (1.0-6      -> 2.0-4     ) [CRAN]
 59: spatstat.... (3.0-2      -> 3.2-0     ) [CRAN]
 60: spatstat.... (3.0-1      -> 3.1-9     ) [CRAN]
 61: tensor       (1.5        -> 1.5.1     ) [CRAN]
 62: abind        (1.4-5      -> 1.4-8     ) [CRAN]
 63: spatstat.... (3.0-1      -> 3.1-0     ) [CRAN]
 64: spatstat.... (3.1-4      -> 3.4-3     ) [CRAN]
 65: spatstat.... (3.1-0      -> 3.6-1     ) [CRAN]
 66: fs           (1.5.2      -> 1.6.6     ) [CRAN]
 67: sass         (0.4.5      -> 0.4.10    ) [CRAN]
 68: promises     (1.2.0.1    -> 1.5.0     ) [CRAN]
 69: mime         (0.12       -> 0.13      ) [CRAN]
 70: later        (1.3.0      -> 1.4.4     ) [CRAN]
 71: jsonlite     (1.8.4      -> 2.0.0     ) [CRAN]
 72: httpuv       (1.6.9      -> 1.6.16    ) [CRAN]
 73: htmltools    (0.5.5      -> 0.5.9     ) [CRAN]
 74: fontawesome  (0.5.2      -> 0.5.3     ) [CRAN]
 75: fastmap      (1.1.1      -> 1.2.0     ) [CRAN]
 76: commonmark   (1.9.0      -> 2.0.0     ) [CRAN]
 77: cachem       (1.0.7      -> 1.1.0     ) [CRAN]
 78: bslib        (0.5.1      -> 0.9.0     ) [CRAN]
 79: RcppArmad... (0.12.2.0.0 -> 15.2.2-1  ) [CRAN]
 80: matrixStats  (0.63.0     -> 1.5.0     ) [CRAN]
 81: reshape2     (1.4.4      -> 1.4.5     ) [CRAN]
 82: bitops       (1.0-7      -> 1.0-9     ) [CRAN]
 83: caTools      (1.18.2     -> 1.18.3    ) [CRAN]
 84: gtools       (3.9.4      -> 3.9.5     ) [CRAN]
 85: gplots       (3.1.3      -> 3.3.0     ) [CRAN]
 86: rprojroot    (2.0.3      -> 2.1.1     ) [CRAN]
 87: here         (1.0.1      -> 1.0.2     ) [CRAN]
 88: RcppTOML     (0.2.2      -> 0.2.3     ) [CRAN]
 89: sys          (3.4.1      -> 3.4.3     ) [CRAN]
 90: tinytex      (0.46       -> 0.58      ) [CRAN]
 91: askpass      (1.1        -> 1.2.1     ) [CRAN]
 92: xfun         (0.39       -> 0.54      ) [CRAN]
 93: highr        (0.10       -> 0.11      ) [CRAN]
 94: evaluate     (0.21       -> 1.0.5     ) [CRAN]
 95: openssl      (2.0.6      -> 2.3.4     ) [CRAN]
 96: curl         (5.0.0      -> 7.0.0     ) [CRAN]
 97: yaml         (2.3.7      -> 2.3.12    ) [CRAN]
 98: rmarkdown    (2.25       -> 2.30      ) [CRAN]
 99: knitr        (1.44       -> 1.50      ) [CRAN]
100: crosstalk    (1.2.0      -> 1.2.2     ) [CRAN]
101: htmlwidgets  (1.6.2      -> 1.6.4     ) [CRAN]
102: shiny        (1.7.5      -> 1.12.1    ) [CRAN]
103: zoo          (1.8-12     -> 1.8-14    ) [CRAN]
104: dotCall64    (1.0-2      -> 1.2       ) [CRAN]
105: spam         (2.9-1      -> 2.11-1    ) [CRAN]
106: sp           (1.6-0      -> 2.2-0     ) [CRAN]
107: rbibutils    (2.2.13     -> 2.4       ) [CRAN]
108: recipes      (1.0.8      -> 1.3.1     ) [CRAN]
109: pROC         (1.18.0     -> 1.19.0.1  ) [CRAN]
110: e1071        (1.7-13     -> 1.7-16    ) [CRAN]
111: forcats      (1.0.0      -> 1.0.1     ) [CRAN]
112: coda         (0.19-4     -> 0.19-4.1  ) [CRAN]
113: broom        (1.0.5      -> 1.0.11    ) [CRAN]
114: uwot         (0.1.14     -> 0.2.4     ) [CRAN]
115: spatstat.... (3.1-0      -> 3.6-0     ) [CRAN]
116: sctransform  (0.4.1      -> 0.4.2     ) [CRAN]
117: Rtsne        (0.16       -> 0.17      ) [CRAN]
118: reticulate   (1.28       -> 1.44.1    ) [CRAN]
119: RcppHNSW     (0.4.1      -> 0.6.0     ) [CRAN]
120: RANN         (2.6.1      -> 2.6.2     ) [CRAN]
121: plotly       (4.10.2     -> 4.11.0    ) [CRAN]
122: pbapply      (1.7-2      -> 1.7-4     ) [CRAN]
123: patchwork    (1.1.3      -> 1.3.2     ) [CRAN]
124: miniUI       (0.1.1.1    -> 0.1.2     ) [CRAN]
125: igraph       (1.4.2      -> 2.2.1     ) [CRAN]
126: ggridges     (0.5.4      -> 0.5.7     ) [CRAN]
127: ggrepel      (0.9.3      -> 0.9.6     ) [CRAN]
128: fitdistrplus (1.1-11     -> 1.2-4     ) [CRAN]
129: fastDummies  (1.7.3      -> 1.7.5     ) [CRAN]
130: SeuratObject (4.1.3      -> 5.3.0     ) [CRAN]
131: nloptr       (2.0.3      -> 2.2.1     ) [CRAN]
132: minqa        (1.2.5      -> 1.2.8     ) [CRAN]
133: harmony      (0.1.1      -> 1.2.4     ) [CRAN]
134: caret        (6.0-94     -> 7.0-1     ) [CRAN]
135: Seurat       (4.3.0      -> 5.4.0     ) [CRAN]
136: lme4         (1.1-32     -> 1.1-38    ) [CRAN]
```



<!-- rnb-output-end -->

<!-- rnb-output-begin eyJkYXRhIjoiSW5zdGFsbGluZyAxIHBhY2thZ2VzOiBzcGF0c3RhdC51bml2YXJcbkluc3RhbGxpbmcgcGFja2FnZSBpbnRvIJFDOi9Vc2Vycy9KdWFuL0RvY3VtZW50cy9SL3dpbi1saWJyYXJ5LzQuMZJcbihhcyCRbGlikiBpcyB1bnNwZWNpZmllZClcbmFsc28gaW5zdGFsbGluZyB0aGUgZGVwZW5kZW5jeSCRc3BhdHN0YXQudXRpbHOSXG5cblJlZ2lzdGVyZWQgUzMgbWV0aG9kIG92ZXJ3cml0dGVuIGJ5ICdkYXRhLnRhYmxlJzpcbiAgbWV0aG9kICAgICAgICAgICBmcm9tXG4gIHByaW50LmRhdGEudGFibGUgICAgIFxuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCDigJhkYXRhLnRhYmxlOjpsYXN04oCZIGJ5IOKAmGRwbHlyOjpsYXN04oCZIHdoZW4gbG9hZGluZyDigJhzY0xBU0VS4oCZXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0IOKAmGRhdGEudGFibGU6OmZpcnN04oCZIGJ5IOKAmGRwbHlyOjpmaXJzdOKAmSB3aGVuIGxvYWRpbmcg4oCYc2NMQVNFUuKAmVxuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCDigJhNQVNTOjpzZWxlY3TigJkgYnkg4oCYZHBseXI6OnNlbGVjdOKAmSB3aGVuIGxvYWRpbmcg4oCYc2NMQVNFUuKAmVxuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCDigJhkYXRhLnRhYmxlOjpiZXR3ZWVu4oCZIGJ5IOKAmGRwbHlyOjpiZXR3ZWVu4oCZIHdoZW4gbG9hZGluZyDigJhzY0xBU0VS4oCZXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0IOKAmGNhcmV0OjpwbHNkYeKAmSBieSDigJhtaXhPbWljczo6cGxzZGHigJkgd2hlbiBsb2FkaW5nIOKAmHNjTEFTRVLigJlcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQg4oCYY2FyZXQ6OnNwbHNkYeKAmSBieSDigJhtaXhPbWljczo6c3Bsc2Rh4oCZIHdoZW4gbG9hZGluZyDigJhzY0xBU0VS4oCZXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0IOKAmGNhcmV0OjpuZWFyWmVyb1ZhcuKAmSBieSDigJhtaXhPbWljczo6bmVhclplcm9WYXLigJkgd2hlbiBsb2FkaW5nIOKAmHNjTEFTRVLigJlcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQg4oCYZHBseXI6OmNvbGxhcHNl4oCZIGJ5IOKAmG5sbWU6OmNvbGxhcHNl4oCZIHdoZW4gbG9hZGluZyDigJhzY0xBU0VS4oCZXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0IOKAmGxtZTQ6OmxtTGlzdOKAmSBieSDigJhubG1lOjpsbUxpc3TigJkgd2hlbiBsb2FkaW5nIOKAmHNjTEFTRVLigJlcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQg4oCYZm9yZWFjaDo6d2hlbuKAmSBieSDigJhwdXJycjo6d2hlbuKAmSB3aGVuIGxvYWRpbmcg4oCYc2NMQVNFUuKAmVxuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCDigJhtaXhPbWljczo6bWFw4oCZIGJ5IOKAmHB1cnJyOjptYXDigJkgd2hlbiBsb2FkaW5nIOKAmHNjTEFTRVLigJlcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQg4oCYY2FyZXQ6OmxpZnTigJkgYnkg4oCYcHVycnI6OmxpZnTigJkgd2hlbiBsb2FkaW5nIOKAmHNjTEFTRVLigJlcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQg4oCYZGF0YS50YWJsZTo6dHJhbnNwb3Nl4oCZIGJ5IOKAmHB1cnJyOjp0cmFuc3Bvc2XigJkgd2hlbiBsb2FkaW5nIOKAmHNjTEFTRVLigJlcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQg4oCYZm9yZWFjaDo6YWNjdW11bGF0ZeKAmSBieSDigJhwdXJycjo6YWNjdW11bGF0ZeKAmSB3aGVuIGxvYWRpbmcg4oCYc2NMQVNFUuKAmVxuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCDigJhNYXRyaXg6OmNvdjJjb3LigJkgYnkg4oCYc3RhdHM6OmNvdjJjb3LigJkgd2hlbiBsb2FkaW5nIOKAmHNjTEFTRVLigJlcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQg4oCYZHBseXI6OmZpbHRlcuKAmSBieSDigJhzdGF0czo6ZmlsdGVy4oCZIHdoZW4gbG9hZGluZyDigJhzY0xBU0VS4oCZXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0IOKAmGRwbHlyOjpsYWfigJkgYnkg4oCYc3RhdHM6OmxhZ+KAmSB3aGVuIGxvYWRpbmcg4oCYc2NMQVNFUuKAmVxuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCDigJhNYXRyaXg6OnRvZXBsaXR64oCZIGJ5IOKAmHN0YXRzOjp0b2VwbGl0euKAmSB3aGVuIGxvYWRpbmcg4oCYc2NMQVNFUuKAmVxuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCDigJhNYXRyaXg6OnVwZGF0ZeKAmSBieSDigJhzdGF0czo6dXBkYXRl4oCZIHdoZW4gbG9hZGluZyDigJhzY0xBU0VS4oCZXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0IOKAmE1hdHJpeDo6ZXhwYW5k4oCZIGJ5IOKAmHRpZHlyOjpleHBhbmTigJkgd2hlbiBsb2FkaW5nIOKAmHNjTEFTRVLigJlcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQg4oCYTWF0cml4OjpwYWNr4oCZIGJ5IOKAmHRpZHlyOjpwYWNr4oCZIHdoZW4gbG9hZGluZyDigJhzY0xBU0VS4oCZXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0IOKAmE1hdHJpeDo6dW5wYWNr4oCZIGJ5IOKAmHRpZHlyOjp1bnBhY2vigJkgd2hlbiBsb2FkaW5nIOKAmHNjTEFTRVLigJlcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQg4oCYTWF0cml4Ojp0YWls4oCZIGJ5IOKAmHV0aWxzOjp0YWls4oCZIHdoZW4gbG9hZGluZyDigJhzY0xBU0VS4oCZXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0IOKAmE1hdHJpeDo6aGVhZOKAmSBieSDigJh1dGlsczo6aGVhZOKAmSB3aGVuIGxvYWRpbmcg4oCYc2NMQVNFUuKAmVxuIn0= -->

```
Installing 1 packages: spatstat.univar
Installing package into 㤼㸱C:/Users/Juan/Documents/R/win-library/4.1㤼㸲
(as 㤼㸱lib㤼㸲 is unspecified)
also installing the dependency 㤼㸱spatstat.utils㤼㸲

Registered S3 method overwritten by 'data.table':
  method           from
  print.data.table     
Warning: replacing previous import ‘data.table::last’ by ‘dplyr::last’ when loading ‘scLASER’
Warning: replacing previous import ‘data.table::first’ by ‘dplyr::first’ when loading ‘scLASER’
Warning: replacing previous import ‘MASS::select’ by ‘dplyr::select’ when loading ‘scLASER’
Warning: replacing previous import ‘data.table::between’ by ‘dplyr::between’ when loading ‘scLASER’
Warning: replacing previous import ‘caret::plsda’ by ‘mixOmics::plsda’ when loading ‘scLASER’
Warning: replacing previous import ‘caret::splsda’ by ‘mixOmics::splsda’ when loading ‘scLASER’
Warning: replacing previous import ‘caret::nearZeroVar’ by ‘mixOmics::nearZeroVar’ when loading ‘scLASER’
Warning: replacing previous import ‘dplyr::collapse’ by ‘nlme::collapse’ when loading ‘scLASER’
Warning: replacing previous import ‘lme4::lmList’ by ‘nlme::lmList’ when loading ‘scLASER’
Warning: replacing previous import ‘foreach::when’ by ‘purrr::when’ when loading ‘scLASER’
Warning: replacing previous import ‘mixOmics::map’ by ‘purrr::map’ when loading ‘scLASER’
Warning: replacing previous import ‘caret::lift’ by ‘purrr::lift’ when loading ‘scLASER’
Warning: replacing previous import ‘data.table::transpose’ by ‘purrr::transpose’ when loading ‘scLASER’
Warning: replacing previous import ‘foreach::accumulate’ by ‘purrr::accumulate’ when loading ‘scLASER’
Warning: replacing previous import ‘Matrix::cov2cor’ by ‘stats::cov2cor’ when loading ‘scLASER’
Warning: replacing previous import ‘dplyr::filter’ by ‘stats::filter’ when loading ‘scLASER’
Warning: replacing previous import ‘dplyr::lag’ by ‘stats::lag’ when loading ‘scLASER’
Warning: replacing previous import ‘Matrix::toeplitz’ by ‘stats::toeplitz’ when loading ‘scLASER’
Warning: replacing previous import ‘Matrix::update’ by ‘stats::update’ when loading ‘scLASER’
Warning: replacing previous import ‘Matrix::expand’ by ‘tidyr::expand’ when loading ‘scLASER’
Warning: replacing previous import ‘Matrix::pack’ by ‘tidyr::pack’ when loading ‘scLASER’
Warning: replacing previous import ‘Matrix::unpack’ by ‘tidyr::unpack’ when loading ‘scLASER’
Warning: replacing previous import ‘Matrix::tail’ by ‘utils::tail’ when loading ‘scLASER’
Warning: replacing previous import ‘Matrix::head’ by ‘utils::head’ when loading ‘scLASER’
```



<!-- rnb-output-end -->

<!-- rnb-output-begin eyJkYXRhIjoiXG4gIFRoZXJlIGlzIGEgYmluYXJ5IHZlcnNpb24gYXZhaWxhYmxlIGJ1dCB0aGUgc291cmNlIHZlcnNpb25cbiAgaXMgbGF0ZXI6XG4ifQ== -->

```

  There is a binary version available but the source version
  is later:
```



<!-- rnb-output-end -->

<!-- rnb-frame-begin eyJtZXRhZGF0YSI6eyJjbGFzc2VzIjpbImRhdGEuZnJhbWUiXSwibnJvdyI6MSwibmNvbCI6M30sInJkZiI6Ikg0c0lBQUFBQUFBQUJtMVEyMDRESVJDZHZWQlRFcXZHNzFqU2J1TVh0QjlnZk9xYm1iS29SQW9iWUZQOStTcUxTOUppSHdnelorYWNtVE12MjkyYTdpZ0FWRkFYSlZRa2hERGJQSy9hcHhhZ0xrTldRQTN6OGY4S1hZK3hGZUQrckVEV2JObTAvOEcyV1lhQVJqQytDem1pOFNEY3hLb21jTGFYR3UxM3lwd1pMQmRUOXFDRjZOd3JONGRlS3ZUUzZFeHhiczJSbmF1bXdzTDE2SjFIendZdmxjc1g0UXBkVHFFZGVtUnZOcWlGN0pSUmJrdy96ZytrY2p3SXljaUZ6WUM3UVk5N2RRMy9HUFJuczdvOEM5eGVpYXUva2VYUEpFWFNUWVIrbHpyZGhDamNDNVZzQnYvUlB1dXQxRDQ1Q2FoajNuaE1mWlFibFpEb0RVNi9hUlAyc1FnQ0FBQT0ifQ== -->

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["binary"],"name":[1],"type":["chr"],"align":["left"]},{"label":["source"],"name":[2],"type":["chr"],"align":["left"]},{"label":["needs_compilation"],"name":[3],"type":["lgl"],"align":["right"]}],"data":[{"1":"3.0-2","2":"3.2-0","3":"TRUE","_rn_":"spatstat.utils"}],"options":{"columns":{"min":{},"max":[10],"total":[3]},"rows":{"min":[10],"max":[10],"total":[1]},"pages":{}}}
  </script>
</div>

<!-- rnb-frame-end -->


<!-- rnb-output-begin eyJkYXRhIjoiUGFja2FnZSB3aGljaCBpcyBvbmx5IGF2YWlsYWJsZSBpbiBzb3VyY2UgZm9ybSwgYW5kIG1heSBuZWVkXG4gIGNvbXBpbGF0aW9uIG9mIEMvQysrL0ZvcnRyYW46IOOkvOO4sXNwYXRzdGF0LnVuaXZhcuOkvOO4slxudHJ5aW5nIFVSTCAnaHR0cHM6Ly9jcmFuLnJzdHVkaW8uY29tL2Jpbi93aW5kb3dzL2NvbnRyaWIvNC4xL3NwYXRzdGF0LnV0aWxzXzMuMC0yLnppcCdcbkNvbnRlbnQgdHlwZSAnYXBwbGljYXRpb24vemlwJyBsZW5ndGggMzkwMzA3IGJ5dGVzICgzODEgS0IpXG5kb3dubG9hZGVkIDM4MSBLQlxuIn0= -->

```
Package which is only available in source form, and may need
  compilation of C/C++/Fortran: 㤼㸱spatstat.univar㤼㸲
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.1/spatstat.utils_3.0-2.zip'
Content type 'application/zip' length 390307 bytes (381 KB)
downloaded 381 KB
```



<!-- rnb-output-end -->

<!-- rnb-output-begin eyJkYXRhIjoicGFja2FnZSDigJhzcGF0c3RhdC51dGlsc+KAmSBzdWNjZXNzZnVsbHkgdW5wYWNrZWQgYW5kIE1ENSBzdW1zIGNoZWNrZWRcbiJ9 -->

```
package ‘spatstat.utils’ successfully unpacked and MD5 sums checked
```



<!-- rnb-output-end -->

<!-- rnb-output-begin eyJkYXRhIjoiY2Fubm90IHJlbW92ZSBwcmlvciBpbnN0YWxsYXRpb24gb2YgcGFja2FnZSDjpLzjuLFzcGF0c3RhdC51dGlsc+OkvOO4snByb2JsZW0gY29weWluZyBDOlxcVXNlcnNcXEp1YW5cXERvY3VtZW50c1xcUlxcd2luLWxpYnJhcnlcXDQuMVxcMDBMT0NLXFxzcGF0c3RhdC51dGlsc1xcbGlic1xceDY0XFxzcGF0c3RhdC51dGlscy5kbGwgdG8gQzpcXFVzZXJzXFxKdWFuXFxEb2N1bWVudHNcXFJcXHdpbi1saWJyYXJ5XFw0LjFcXHNwYXRzdGF0LnV0aWxzXFxsaWJzXFx4NjRcXHNwYXRzdGF0LnV0aWxzLmRsbDogUGVybWlzc2lvbiBkZW5pZWRyZXN0b3JlZCDjpLzjuLFzcGF0c3RhdC51dGlsc+OkvOO4slxuIn0= -->

```
cannot remove prior installation of package 㤼㸱spatstat.utils㤼㸲problem copying C:\Users\Juan\Documents\R\win-library\4.1\00LOCK\spatstat.utils\libs\x64\spatstat.utils.dll to C:\Users\Juan\Documents\R\win-library\4.1\spatstat.utils\libs\x64\spatstat.utils.dll: Permission deniedrestored 㤼㸱spatstat.utils㤼㸲
```



<!-- rnb-output-end -->

<!-- rnb-output-begin eyJkYXRhIjoiXG5UaGUgZG93bmxvYWRlZCBiaW5hcnkgcGFja2FnZXMgYXJlIGluXG5cdEM6XFxVc2Vyc1xcSnVhblxcQXBwRGF0YVxcTG9jYWxcXFRlbXBcXFJ0bXBDVUxBNmJcXGRvd25sb2FkZWRfcGFja2FnZXNcbiJ9 -->

```

The downloaded binary packages are in
	C:\Users\Juan\AppData\Local\Temp\RtmpCULA6b\downloaded_packages
```



<!-- rnb-output-end -->

<!-- rnb-output-begin eyJkYXRhIjoiaW5zdGFsbGluZyB0aGUgc291cmNlIHBhY2thZ2Ug46S847ixc3BhdHN0YXQudW5pdmFy46S847iyXG5cbnRyeWluZyBVUkwgJ2h0dHBzOi8vY3Jhbi5yc3R1ZGlvLmNvbS9zcmMvY29udHJpYi9zcGF0c3RhdC51bml2YXJfMy4xLTUudGFyLmd6J1xuQ29udGVudCB0eXBlICdhcHBsaWNhdGlvbi94LWd6aXAnIGxlbmd0aCA3MjI5OSBieXRlcyAoNzAgS0IpXG5kb3dubG9hZGVkIDcwIEtCXG4ifQ== -->

```
installing the source package 㤼㸱spatstat.univar㤼㸲

trying URL 'https://cran.rstudio.com/src/contrib/spatstat.univar_3.1-5.tar.gz'
Content type 'application/x-gzip' length 72299 bytes (70 KB)
downloaded 70 KB
```



<!-- rnb-output-end -->

<!-- rnb-output-begin eyJkYXRhIjoiKiBpbnN0YWxsaW5nICpzb3VyY2UqIHBhY2thZ2UgJ3NwYXRzdGF0LnVuaXZhcicgLi4uXG4qKiBwYWNrYWdlICdzcGF0c3RhdC51bml2YXInIHN1Y2Nlc3NmdWxseSB1bnBhY2tlZCBhbmQgTUQ1IHN1bXMgY2hlY2tlZFxuKiogdXNpbmcgc3RhZ2VkIGluc3RhbGxhdGlvblxuKiogbGlic1xuXG4qKiogYXJjaCAtIGkzODZcblwiQzovUkJ1aWxkVG9vbHMvNC4wL21pbmd3MzIvYmluL1wiZ2NjICAtSVwiQzovUFJPR1JBfjEvUi9SLTQxfjEuMi9pbmNsdWRlXCIgLUROREVCVUcgICAgICAgICAgLU8yIC1XYWxsICAtc3RkPWdudTk5IC1tZnBtYXRoPXNzZSAtbXNzZTIgLW1zdGFja3JlYWxpZ24gIC1jIGFjY2Vzcy5jIC1vIGFjY2Vzcy5vXG5cIkM6L1JCdWlsZFRvb2xzLzQuMC9taW5ndzMyL2Jpbi9cImdjYyAgLUlcIkM6L1BST0dSQX4xL1IvUi00MX4xLjIvaW5jbHVkZVwiIC1ETkRFQlVHICAgICAgICAgIC1PMiAtV2FsbCAgLXN0ZD1nbnU5OSAtbWZwbWF0aD1zc2UgLW1zc2UyIC1tc3RhY2tyZWFsaWduICAtYyBhZGFwdGl2ZS5jIC1vIGFkYXB0aXZlLm9cblwiQzovUkJ1aWxkVG9vbHMvNC4wL21pbmd3MzIvYmluL1wiZ2NjICAtSVwiQzovUFJPR1JBfjEvUi9SLTQxfjEuMi9pbmNsdWRlXCIgLUROREVCVUcgICAgICAgICAgLU8yIC1XYWxsICAtc3RkPWdudTk5IC1tZnBtYXRoPXNzZSAtbXNzZTIgLW1zdGFja3JlYWxpZ24gIC1jIGNvbG9uZWwuYyAtbyBjb2xvbmVsLm9cbmNvbG9uZWwuYzogSW4gZnVuY3Rpb24gJ2NvbG9uZWwnOlxuY29sb25lbC5jOjU5OjM3OiB3YXJuaW5nOiB1bnVzZWQgdmFyaWFibGUgJ3Jvb3QycGknIFstV3VudXNlZC12YXJpYWJsZV1cbiAgIGRvdWJsZSB4aSwgd2ksIHVpaiwgdGVtcCwga3ZhbHVlLCByb290MnBpO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF5+fn5+fn5cbmNvbG9uZWwuYzo1ODoyMTogd2FybmluZzogdW51c2VkIHZhcmlhYmxlICdiYicgWy1XdW51c2VkLXZhcmlhYmxlXVxuICAgaW50IGksIGosIE54LCBOciwgYmI7XG4gICAgICAgICAgICAgICAgICAgICBeflxuY29sb25lbC5jOiBJbiBmdW5jdGlvbiAnZmNvbG9uZWwnOlxuY29sb25lbC5jOjIxNzo0MTogd2FybmluZzogdW51c2VkIHZhcmlhYmxlICdyb290MnBpJyBbLVd1bnVzZWQtdmFyaWFibGVdXG4gICBkb3VibGUgZHIsIHhpLCB3aSwgdmlqLCB0ZW1wLCBrdmFsdWUsIHJvb3QycGk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF5+fn5+fn5cbmNvbG9uZWwuYzoyMTY6MjQ6IHdhcm5pbmc6IHVudXNlZCB2YXJpYWJsZSAnYmInIFstV3VudXNlZC12YXJpYWJsZV1cbiAgIGludCBpLCBqLCBrLCBOeCwgTnIsIGJiO1xuICAgICAgICAgICAgICAgICAgICAgICAgXn5cbmNvbG9uZWwuYzogSW4gZnVuY3Rpb24gJ2Jjb2xvbmVsJzpcbmNvbG9uZWwuYzo0MTg6NDA6IHdhcm5pbmc6IHVudXNlZCB2YXJpYWJsZSAncm9vdDJwaScgWy1XdW51c2VkLXZhcmlhYmxlXVxuICAgZG91YmxlIGt2YWx1ZSwgZGVub21qLCB0ZW1wLCB0aHJlc2gsIHJvb3QycGk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXn5+fn5+flxuY29sb25lbC5jOjQxNjoyMTogd2FybmluZzogdW51c2VkIHZhcmlhYmxlICdiYicgWy1XdW51c2VkLXZhcmlhYmxlXVxuICAgaW50IGksIGosIE54LCBOciwgYmIsIGpiZHJ5O1xuICAgICAgICAgICAgICAgICAgICAgXn5cbmNvbG9uZWwuYzogSW4gZnVuY3Rpb24gJ2ZiY29sb25lbCc6XG5jb2xvbmVsLmM6NjUzOjQwOiB3YXJuaW5nOiB1bnVzZWQgdmFyaWFibGUgJ3Jvb3QycGknIFstV3VudXNlZC12YXJpYWJsZV1cbiAgIGRvdWJsZSBrdmFsdWUsIGRlbm9taiwgdGVtcCwgdGhyZXNoLCByb290MnBpO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF5+fn5+fn5cbmNvbG9uZWwuYzo2NTE6MzI6IHdhcm5pbmc6IHVudXNlZCB2YXJpYWJsZSAnanVwcGVyaScgWy1XdW51c2VkLXZhcmlhYmxlXVxuICAgaW50IGksIGosIE54LCBOciwgYmIsIGpiZHJ5LCBqdXBwZXJpO1xuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBefn5+fn5+XG5jb2xvbmVsLmM6NjUxOjIxOiB3YXJuaW5nOiB1bnVzZWQgdmFyaWFibGUgJ2JiJyBbLVd1bnVzZWQtdmFyaWFibGVdXG4gICBpbnQgaSwgaiwgTngsIE5yLCBiYiwgamJkcnksIGp1cHBlcmk7XG4gICAgICAgICAgICAgICAgICAgICBeflxuXCJDOi9SQnVpbGRUb29scy80LjAvbWluZ3czMi9iaW4vXCJnY2MgIC1JXCJDOi9QUk9HUkF+MS9SL1ItNDF+MS4yL2luY2x1ZGVcIiAtRE5ERUJVRyAgICAgICAgICAtTzIgLVdhbGwgIC1zdGQ9Z251OTkgLW1mcG1hdGg9c3NlIC1tc3NlMiAtbXN0YWNrcmVhbGlnbiAgLWMgaG90cm9kLmMgLW8gaG90cm9kLm9cblwiQzovUkJ1aWxkVG9vbHMvNC4wL21pbmd3MzIvYmluL1wiZ2NjICAtSVwiQzovUFJPR1JBfjEvUi9SLTQxfjEuMi9pbmNsdWRlXCIgLUROREVCVUcgICAgICAgICAgLU8yIC1XYWxsICAtc3RkPWdudTk5IC1tZnBtYXRoPXNzZSAtbXNzZTIgLW1zdGFja3JlYWxpZ24gIC1jIGluaXQuYyAtbyBpbml0Lm9cblwiQzovUkJ1aWxkVG9vbHMvNC4wL21pbmd3MzIvYmluL1wiZ2NjICAtSVwiQzovUFJPR1JBfjEvUi9SLTQxfjEuMi9pbmNsdWRlXCIgLUROREVCVUcgICAgICAgICAgLU8yIC1XYWxsICAtc3RkPWdudTk5IC1tZnBtYXRoPXNzZSAtbXNzZTIgLW1zdGFja3JlYWxpZ24gIC1jIGtlcm5lbHMuYyAtbyBrZXJuZWxzLm9cblwiQzovUkJ1aWxkVG9vbHMvNC4wL21pbmd3MzIvYmluL1wiZ2NjICAtSVwiQzovUFJPR1JBfjEvUi9SLTQxfjEuMi9pbmNsdWRlXCIgLUROREVCVUcgICAgICAgICAgLU8yIC1XYWxsICAtc3RkPWdudTk5IC1tZnBtYXRoPXNzZSAtbXNzZTIgLW1zdGFja3JlYWxpZ24gIC1jIHRhYm51bS5jIC1vIHRhYm51bS5vXG5cIkM6L1JCdWlsZFRvb2xzLzQuMC9taW5ndzMyL2Jpbi9cImdjYyAgLUlcIkM6L1BST0dSQX4xL1IvUi00MX4xLjIvaW5jbHVkZVwiIC1ETkRFQlVHICAgICAgICAgIC1PMiAtV2FsbCAgLXN0ZD1nbnU5OSAtbWZwbWF0aD1zc2UgLW1zc2UyIC1tc3RhY2tyZWFsaWduICAtYyB0YXlsb3Jib290LmMgLW8gdGF5bG9yYm9vdC5vXG5cIkM6L1JCdWlsZFRvb2xzLzQuMC9taW5ndzMyL2Jpbi9cImdjYyAgLUlcIkM6L1BST0dSQX4xL1IvUi00MX4xLjIvaW5jbHVkZVwiIC1ETkRFQlVHICAgICAgICAgIC1PMiAtV2FsbCAgLXN0ZD1nbnU5OSAtbWZwbWF0aD1zc2UgLW1zc2UyIC1tc3RhY2tyZWFsaWduICAtYyB3aGlzdC5jIC1vIHdoaXN0Lm9cbkM6L1JCdWlsZFRvb2xzLzQuMC9taW5ndzMyL2Jpbi9nY2MgLXNoYXJlZCAtcyAtc3RhdGljLWxpYmdjYyAtbyBzcGF0c3RhdC51bml2YXIuZGxsIHRtcC5kZWYgYWNjZXNzLm8gYWRhcHRpdmUubyBjb2xvbmVsLm8gaG90cm9kLm8gaW5pdC5vIGtlcm5lbHMubyB0YWJudW0ubyB0YXlsb3Jib290Lm8gd2hpc3QubyAtTEM6L1BST0dSQX4xL1IvUi00MX4xLjIvYmluL2kzODYgLWxSXG5pbnN0YWxsaW5nIHRvIEM6L1VzZXJzL0p1YW4vRG9jdW1lbnRzL1Ivd2luLWxpYnJhcnkvNC4xLzAwTE9DSy1zcGF0c3RhdC51bml2YXIvMDBuZXcvc3BhdHN0YXQudW5pdmFyL2xpYnMvaTM4NlxuXG4qKiogYXJjaCAtIHg2NFxuXCJDOi9SQnVpbGRUb29scy80LjAvbWluZ3c2NC9iaW4vXCJnY2MgIC1JXCJDOi9QUk9HUkF+MS9SL1ItNDF+MS4yL2luY2x1ZGVcIiAtRE5ERUJVRyAgICAgICAgICAtTzIgLVdhbGwgIC1zdGQ9Z251OTkgLW1mcG1hdGg9c3NlIC1tc3NlMiAtbXN0YWNrcmVhbGlnbiAgLWMgYWNjZXNzLmMgLW8gYWNjZXNzLm9cblwiQzovUkJ1aWxkVG9vbHMvNC4wL21pbmd3NjQvYmluL1wiZ2NjICAtSVwiQzovUFJPR1JBfjEvUi9SLTQxfjEuMi9pbmNsdWRlXCIgLUROREVCVUcgICAgICAgICAgLU8yIC1XYWxsICAtc3RkPWdudTk5IC1tZnBtYXRoPXNzZSAtbXNzZTIgLW1zdGFja3JlYWxpZ24gIC1jIGFkYXB0aXZlLmMgLW8gYWRhcHRpdmUub1xuXCJDOi9SQnVpbGRUb29scy80LjAvbWluZ3c2NC9iaW4vXCJnY2MgIC1JXCJDOi9QUk9HUkF+MS9SL1ItNDF+MS4yL2luY2x1ZGVcIiAtRE5ERUJVRyAgICAgICAgICAtTzIgLVdhbGwgIC1zdGQ9Z251OTkgLW1mcG1hdGg9c3NlIC1tc3NlMiAtbXN0YWNrcmVhbGlnbiAgLWMgY29sb25lbC5jIC1vIGNvbG9uZWwub1xuY29sb25lbC5jOiBJbiBmdW5jdGlvbiAnY29sb25lbCc6XG5jb2xvbmVsLmM6NTk6Mzc6IHdhcm5pbmc6IHVudXNlZCB2YXJpYWJsZSAncm9vdDJwaScgWy1XdW51c2VkLXZhcmlhYmxlXVxuICAgZG91YmxlIHhpLCB3aSwgdWlqLCB0ZW1wLCBrdmFsdWUsIHJvb3QycGk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXn5+fn5+flxuY29sb25lbC5jOjU4OjIxOiB3YXJuaW5nOiB1bnVzZWQgdmFyaWFibGUgJ2JiJyBbLVd1bnVzZWQtdmFyaWFibGVdXG4gICBpbnQgaSwgaiwgTngsIE5yLCBiYjtcbiAgICAgICAgICAgICAgICAgICAgIF5+XG5jb2xvbmVsLmM6IEluIGZ1bmN0aW9uICdmY29sb25lbCc6XG5jb2xvbmVsLmM6MjE3OjQxOiB3YXJuaW5nOiB1bnVzZWQgdmFyaWFibGUgJ3Jvb3QycGknIFstV3VudXNlZC12YXJpYWJsZV1cbiAgIGRvdWJsZSBkciwgeGksIHdpLCB2aWosIHRlbXAsIGt2YWx1ZSwgcm9vdDJwaTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXn5+fn5+flxuY29sb25lbC5jOjIxNjoyNDogd2FybmluZzogdW51c2VkIHZhcmlhYmxlICdiYicgWy1XdW51c2VkLXZhcmlhYmxlXVxuICAgaW50IGksIGosIGssIE54LCBOciwgYmI7XG4gICAgICAgICAgICAgICAgICAgICAgICBeflxuY29sb25lbC5jOiBJbiBmdW5jdGlvbiAnYmNvbG9uZWwnOlxuY29sb25lbC5jOjQxODo0MDogd2FybmluZzogdW51c2VkIHZhcmlhYmxlICdyb290MnBpJyBbLVd1bnVzZWQtdmFyaWFibGVdXG4gICBkb3VibGUga3ZhbHVlLCBkZW5vbWosIHRlbXAsIHRocmVzaCwgcm9vdDJwaTtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBefn5+fn5+XG5jb2xvbmVsLmM6NDE2OjIxOiB3YXJuaW5nOiB1bnVzZWQgdmFyaWFibGUgJ2JiJyBbLVd1bnVzZWQtdmFyaWFibGVdXG4gICBpbnQgaSwgaiwgTngsIE5yLCBiYiwgamJkcnk7XG4gICAgICAgICAgICAgICAgICAgICBeflxuY29sb25lbC5jOiBJbiBmdW5jdGlvbiAnZmJjb2xvbmVsJzpcbmNvbG9uZWwuYzo2NTM6NDA6IHdhcm5pbmc6IHVudXNlZCB2YXJpYWJsZSAncm9vdDJwaScgWy1XdW51c2VkLXZhcmlhYmxlXVxuICAgZG91YmxlIGt2YWx1ZSwgZGVub21qLCB0ZW1wLCB0aHJlc2gsIHJvb3QycGk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXn5+fn5+flxuY29sb25lbC5jOjY1MTozMjogd2FybmluZzogdW51c2VkIHZhcmlhYmxlICdqdXBwZXJpJyBbLVd1bnVzZWQtdmFyaWFibGVdXG4gICBpbnQgaSwgaiwgTngsIE5yLCBiYiwgamJkcnksIGp1cHBlcmk7XG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF5+fn5+fn5cbmNvbG9uZWwuYzo2NTE6MjE6IHdhcm5pbmc6IHVudXNlZCB2YXJpYWJsZSAnYmInIFstV3VudXNlZC12YXJpYWJsZV1cbiAgIGludCBpLCBqLCBOeCwgTnIsIGJiLCBqYmRyeSwganVwcGVyaTtcbiAgICAgICAgICAgICAgICAgICAgIF5+XG5cIkM6L1JCdWlsZFRvb2xzLzQuMC9taW5ndzY0L2Jpbi9cImdjYyAgLUlcIkM6L1BST0dSQX4xL1IvUi00MX4xLjIvaW5jbHVkZVwiIC1ETkRFQlVHICAgICAgICAgIC1PMiAtV2FsbCAgLXN0ZD1nbnU5OSAtbWZwbWF0aD1zc2UgLW1zc2UyIC1tc3RhY2tyZWFsaWduICAtYyBob3Ryb2QuYyAtbyBob3Ryb2Qub1xuXCJDOi9SQnVpbGRUb29scy80LjAvbWluZ3c2NC9iaW4vXCJnY2MgIC1JXCJDOi9QUk9HUkF+MS9SL1ItNDF+MS4yL2luY2x1ZGVcIiAtRE5ERUJVRyAgICAgICAgICAtTzIgLVdhbGwgIC1zdGQ9Z251OTkgLW1mcG1hdGg9c3NlIC1tc3NlMiAtbXN0YWNrcmVhbGlnbiAgLWMgaW5pdC5jIC1vIGluaXQub1xuXCJDOi9SQnVpbGRUb29scy80LjAvbWluZ3c2NC9iaW4vXCJnY2MgIC1JXCJDOi9QUk9HUkF+MS9SL1ItNDF+MS4yL2luY2x1ZGVcIiAtRE5ERUJVRyAgICAgICAgICAtTzIgLVdhbGwgIC1zdGQ9Z251OTkgLW1mcG1hdGg9c3NlIC1tc3NlMiAtbXN0YWNrcmVhbGlnbiAgLWMga2VybmVscy5jIC1vIGtlcm5lbHMub1xuXCJDOi9SQnVpbGRUb29scy80LjAvbWluZ3c2NC9iaW4vXCJnY2MgIC1JXCJDOi9QUk9HUkF+MS9SL1ItNDF+MS4yL2luY2x1ZGVcIiAtRE5ERUJVRyAgICAgICAgICAtTzIgLVdhbGwgIC1zdGQ9Z251OTkgLW1mcG1hdGg9c3NlIC1tc3NlMiAtbXN0YWNrcmVhbGlnbiAgLWMgdGFibnVtLmMgLW8gdGFibnVtLm9cblwiQzovUkJ1aWxkVG9vbHMvNC4wL21pbmd3NjQvYmluL1wiZ2NjICAtSVwiQzovUFJPR1JBfjEvUi9SLTQxfjEuMi9pbmNsdWRlXCIgLUROREVCVUcgICAgICAgICAgLU8yIC1XYWxsICAtc3RkPWdudTk5IC1tZnBtYXRoPXNzZSAtbXNzZTIgLW1zdGFja3JlYWxpZ24gIC1jIHRheWxvcmJvb3QuYyAtbyB0YXlsb3Jib290Lm9cblwiQzovUkJ1aWxkVG9vbHMvNC4wL21pbmd3NjQvYmluL1wiZ2NjICAtSVwiQzovUFJPR1JBfjEvUi9SLTQxfjEuMi9pbmNsdWRlXCIgLUROREVCVUcgICAgICAgICAgLU8yIC1XYWxsICAtc3RkPWdudTk5IC1tZnBtYXRoPXNzZSAtbXNzZTIgLW1zdGFja3JlYWxpZ24gIC1jIHdoaXN0LmMgLW8gd2hpc3Qub1xuQzovUkJ1aWxkVG9vbHMvNC4wL21pbmd3NjQvYmluL2djYyAtc2hhcmVkIC1zIC1zdGF0aWMtbGliZ2NjIC1vIHNwYXRzdGF0LnVuaXZhci5kbGwgdG1wLmRlZiBhY2Nlc3MubyBhZGFwdGl2ZS5vIGNvbG9uZWwubyBob3Ryb2QubyBpbml0Lm8ga2VybmVscy5vIHRhYm51bS5vIHRheWxvcmJvb3QubyB3aGlzdC5vIC1MQzovUFJPR1JBfjEvUi9SLTQxfjEuMi9iaW4veDY0IC1sUlxuaW5zdGFsbGluZyB0byBDOi9Vc2Vycy9KdWFuL0RvY3VtZW50cy9SL3dpbi1saWJyYXJ5LzQuMS8wMExPQ0stc3BhdHN0YXQudW5pdmFyLzAwbmV3L3NwYXRzdGF0LnVuaXZhci9saWJzL3g2NFxuKiogUlxuKiogaW5zdFxuKiogYnl0ZS1jb21waWxlIGFuZCBwcmVwYXJlIHBhY2thZ2UgZm9yIGxhenkgbG9hZGluZ1xuRXJyb3IgaW4gbG9hZE5hbWVzcGFjZShpLCBjKGxpYi5sb2MsIC5saWJQYXRocygpKSwgdmVyc2lvbkNoZWNrID0gdklbW2ldXSkgOiBcbiAgbmFtZXNwYWNlICdzcGF0c3RhdC51dGlscycgMy4wLTIgaXMgYmVpbmcgbG9hZGVkLCBidXQgPj0gMy4xLjIgaXMgcmVxdWlyZWRcbkNhbGxzOiA8QW5vbnltb3VzPiAuLi4gd2l0aENhbGxpbmdIYW5kbGVycyAtPiBsb2FkTmFtZXNwYWNlIC0+IG5hbWVzcGFjZUltcG9ydCAtPiBsb2FkTmFtZXNwYWNlXG5FeGVjdXRpb24gaGFsdGVkXG5FUlJPUjogbGF6eSBsb2FkaW5nIGZhaWxlZCBmb3IgcGFja2FnZSAnc3BhdHN0YXQudW5pdmFyJ1xuKiByZW1vdmluZyAnQzovVXNlcnMvSnVhbi9Eb2N1bWVudHMvUi93aW4tbGlicmFyeS80LjEvc3BhdHN0YXQudW5pdmFyJ1xuIn0= -->

```
* installing *source* package 'spatstat.univar' ...
** package 'spatstat.univar' successfully unpacked and MD5 sums checked
** using staged installation
** libs

*** arch - i386
"C:/RBuildTools/4.0/mingw32/bin/"gcc  -I"C:/PROGRA~1/R/R-41~1.2/include" -DNDEBUG          -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -c access.c -o access.o
"C:/RBuildTools/4.0/mingw32/bin/"gcc  -I"C:/PROGRA~1/R/R-41~1.2/include" -DNDEBUG          -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -c adaptive.c -o adaptive.o
"C:/RBuildTools/4.0/mingw32/bin/"gcc  -I"C:/PROGRA~1/R/R-41~1.2/include" -DNDEBUG          -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -c colonel.c -o colonel.o
colonel.c: In function 'colonel':
colonel.c:59:37: warning: unused variable 'root2pi' [-Wunused-variable]
   double xi, wi, uij, temp, kvalue, root2pi;
                                     ^~~~~~~
colonel.c:58:21: warning: unused variable 'bb' [-Wunused-variable]
   int i, j, Nx, Nr, bb;
                     ^~
colonel.c: In function 'fcolonel':
colonel.c:217:41: warning: unused variable 'root2pi' [-Wunused-variable]
   double dr, xi, wi, vij, temp, kvalue, root2pi;
                                         ^~~~~~~
colonel.c:216:24: warning: unused variable 'bb' [-Wunused-variable]
   int i, j, k, Nx, Nr, bb;
                        ^~
colonel.c: In function 'bcolonel':
colonel.c:418:40: warning: unused variable 'root2pi' [-Wunused-variable]
   double kvalue, denomj, temp, thresh, root2pi;
                                        ^~~~~~~
colonel.c:416:21: warning: unused variable 'bb' [-Wunused-variable]
   int i, j, Nx, Nr, bb, jbdry;
                     ^~
colonel.c: In function 'fbcolonel':
colonel.c:653:40: warning: unused variable 'root2pi' [-Wunused-variable]
   double kvalue, denomj, temp, thresh, root2pi;
                                        ^~~~~~~
colonel.c:651:32: warning: unused variable 'jupperi' [-Wunused-variable]
   int i, j, Nx, Nr, bb, jbdry, jupperi;
                                ^~~~~~~
colonel.c:651:21: warning: unused variable 'bb' [-Wunused-variable]
   int i, j, Nx, Nr, bb, jbdry, jupperi;
                     ^~
"C:/RBuildTools/4.0/mingw32/bin/"gcc  -I"C:/PROGRA~1/R/R-41~1.2/include" -DNDEBUG          -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -c hotrod.c -o hotrod.o
"C:/RBuildTools/4.0/mingw32/bin/"gcc  -I"C:/PROGRA~1/R/R-41~1.2/include" -DNDEBUG          -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -c init.c -o init.o
"C:/RBuildTools/4.0/mingw32/bin/"gcc  -I"C:/PROGRA~1/R/R-41~1.2/include" -DNDEBUG          -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -c kernels.c -o kernels.o
"C:/RBuildTools/4.0/mingw32/bin/"gcc  -I"C:/PROGRA~1/R/R-41~1.2/include" -DNDEBUG          -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -c tabnum.c -o tabnum.o
"C:/RBuildTools/4.0/mingw32/bin/"gcc  -I"C:/PROGRA~1/R/R-41~1.2/include" -DNDEBUG          -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -c taylorboot.c -o taylorboot.o
"C:/RBuildTools/4.0/mingw32/bin/"gcc  -I"C:/PROGRA~1/R/R-41~1.2/include" -DNDEBUG          -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -c whist.c -o whist.o
C:/RBuildTools/4.0/mingw32/bin/gcc -shared -s -static-libgcc -o spatstat.univar.dll tmp.def access.o adaptive.o colonel.o hotrod.o init.o kernels.o tabnum.o taylorboot.o whist.o -LC:/PROGRA~1/R/R-41~1.2/bin/i386 -lR
installing to C:/Users/Juan/Documents/R/win-library/4.1/00LOCK-spatstat.univar/00new/spatstat.univar/libs/i386

*** arch - x64
"C:/RBuildTools/4.0/mingw64/bin/"gcc  -I"C:/PROGRA~1/R/R-41~1.2/include" -DNDEBUG          -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -c access.c -o access.o
"C:/RBuildTools/4.0/mingw64/bin/"gcc  -I"C:/PROGRA~1/R/R-41~1.2/include" -DNDEBUG          -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -c adaptive.c -o adaptive.o
"C:/RBuildTools/4.0/mingw64/bin/"gcc  -I"C:/PROGRA~1/R/R-41~1.2/include" -DNDEBUG          -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -c colonel.c -o colonel.o
colonel.c: In function 'colonel':
colonel.c:59:37: warning: unused variable 'root2pi' [-Wunused-variable]
   double xi, wi, uij, temp, kvalue, root2pi;
                                     ^~~~~~~
colonel.c:58:21: warning: unused variable 'bb' [-Wunused-variable]
   int i, j, Nx, Nr, bb;
                     ^~
colonel.c: In function 'fcolonel':
colonel.c:217:41: warning: unused variable 'root2pi' [-Wunused-variable]
   double dr, xi, wi, vij, temp, kvalue, root2pi;
                                         ^~~~~~~
colonel.c:216:24: warning: unused variable 'bb' [-Wunused-variable]
   int i, j, k, Nx, Nr, bb;
                        ^~
colonel.c: In function 'bcolonel':
colonel.c:418:40: warning: unused variable 'root2pi' [-Wunused-variable]
   double kvalue, denomj, temp, thresh, root2pi;
                                        ^~~~~~~
colonel.c:416:21: warning: unused variable 'bb' [-Wunused-variable]
   int i, j, Nx, Nr, bb, jbdry;
                     ^~
colonel.c: In function 'fbcolonel':
colonel.c:653:40: warning: unused variable 'root2pi' [-Wunused-variable]
   double kvalue, denomj, temp, thresh, root2pi;
                                        ^~~~~~~
colonel.c:651:32: warning: unused variable 'jupperi' [-Wunused-variable]
   int i, j, Nx, Nr, bb, jbdry, jupperi;
                                ^~~~~~~
colonel.c:651:21: warning: unused variable 'bb' [-Wunused-variable]
   int i, j, Nx, Nr, bb, jbdry, jupperi;
                     ^~
"C:/RBuildTools/4.0/mingw64/bin/"gcc  -I"C:/PROGRA~1/R/R-41~1.2/include" -DNDEBUG          -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -c hotrod.c -o hotrod.o
"C:/RBuildTools/4.0/mingw64/bin/"gcc  -I"C:/PROGRA~1/R/R-41~1.2/include" -DNDEBUG          -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -c init.c -o init.o
"C:/RBuildTools/4.0/mingw64/bin/"gcc  -I"C:/PROGRA~1/R/R-41~1.2/include" -DNDEBUG          -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -c kernels.c -o kernels.o
"C:/RBuildTools/4.0/mingw64/bin/"gcc  -I"C:/PROGRA~1/R/R-41~1.2/include" -DNDEBUG          -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -c tabnum.c -o tabnum.o
"C:/RBuildTools/4.0/mingw64/bin/"gcc  -I"C:/PROGRA~1/R/R-41~1.2/include" -DNDEBUG          -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -c taylorboot.c -o taylorboot.o
"C:/RBuildTools/4.0/mingw64/bin/"gcc  -I"C:/PROGRA~1/R/R-41~1.2/include" -DNDEBUG          -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -c whist.c -o whist.o
C:/RBuildTools/4.0/mingw64/bin/gcc -shared -s -static-libgcc -o spatstat.univar.dll tmp.def access.o adaptive.o colonel.o hotrod.o init.o kernels.o tabnum.o taylorboot.o whist.o -LC:/PROGRA~1/R/R-41~1.2/bin/x64 -lR
installing to C:/Users/Juan/Documents/R/win-library/4.1/00LOCK-spatstat.univar/00new/spatstat.univar/libs/x64
** R
** inst
** byte-compile and prepare package for lazy loading
Error in loadNamespace(i, c(lib.loc, .libPaths()), versionCheck = vI[[i]]) : 
  namespace 'spatstat.utils' 3.0-2 is being loaded, but >= 3.1.2 is required
Calls: <Anonymous> ... withCallingHandlers -> loadNamespace -> namespaceImport -> loadNamespace
Execution halted
ERROR: lazy loading failed for package 'spatstat.univar'
* removing 'C:/Users/Juan/Documents/R/win-library/4.1/spatstat.univar'
```



<!-- rnb-output-end -->

<!-- rnb-output-begin eyJkYXRhIjoiaW5zdGFsbGF0aW9uIG9mIHBhY2thZ2UgkXNwYXRzdGF0LnVuaXZhcpIgaGFkIG5vbi16ZXJvIGV4aXQgc3RhdHVzXG5UaGUgZG93bmxvYWRlZCBzb3VyY2UgcGFja2FnZXMgYXJlIGluXG5cdOKAmEM6XFxVc2Vyc1xcSnVhblxcQXBwRGF0YVxcTG9jYWxcXFRlbXBcXFJ0bXBDVUxBNmJcXGRvd25sb2FkZWRfcGFja2FnZXPigJlcbiJ9 -->

```
installation of package 㤼㸱spatstat.univar㤼㸲 had non-zero exit status
The downloaded source packages are in
	‘C:\Users\Juan\AppData\Local\Temp\RtmpCULA6b\downloaded_packages’
```



<!-- rnb-output-end -->

<!-- rnb-output-begin eyJkYXRhIjoiLS0gUiBDTUQgYnVpbGQgLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLS0tLVxuICBcbiAgXG4gIFxudiAgY2hlY2tpbmcgZm9yIGZpbGUgJ0M6XFxVc2Vyc1xcSnVhblxcQXBwRGF0YVxcTG9jYWxcXFRlbXBcXFJ0bXBDVUxBNmJcXHJlbW90ZXM0MGM0YzQxMzE1ZFxcZmFuemhhbmdsYWItc2NMQVNFUi05ZWQwOGQ0ZjZkNTQ0MWM0NmVlYjhkMjVkMWYzMzJkYzMwNDAzYzI2L0RFU0NSSVBUSU9OJ1xuXG4gIFxuICBcbiAgXG4tICBwcmVwYXJpbmcgJ3NjTEFTRVInOlxuICAgY2hlY2tpbmcgREVTQ1JJUFRJT04gbWV0YS1pbmZvcm1hdGlvbiAuLi5cbiAgXG4gICBjaGVja2luZyBERVNDUklQVElPTiBtZXRhLWluZm9ybWF0aW9uIC4uLiBcbiAgXG52ICBjaGVja2luZyBERVNDUklQVElPTiBtZXRhLWluZm9ybWF0aW9uXG5cbiAgXG4gIFxuICBcbi0gIGNoZWNraW5nIGZvciBMRiBsaW5lLWVuZGluZ3MgaW4gc291cmNlIGFuZCBtYWtlIGZpbGVzIGFuZCBzaGVsbCBzY3JpcHRzXG5cbiAgXG4tICBjaGVja2luZyBmb3IgZW1wdHkgb3IgdW5uZWVkZWQgZGlyZWN0b3JpZXNcblxuICBcbi0gIGJ1aWxkaW5nICdzY0xBU0VSXzAuMC4wLjkwMDAudGFyLmd6J1xuXG4gIFxuICAgXG4ifQ== -->

```
-- R CMD build -------------------------------------------------------
  
  
  
v  checking for file 'C:\Users\Juan\AppData\Local\Temp\RtmpCULA6b\remotes40c4c41315d\fanzhanglab-scLASER-9ed08d4f6d5441c46eeb8d25d1f332dc30403c26/DESCRIPTION'

  
  
  
-  preparing 'scLASER':
   checking DESCRIPTION meta-information ...
  
   checking DESCRIPTION meta-information ... 
  
v  checking DESCRIPTION meta-information

  
  
  
-  checking for LF line-endings in source and make files and shell scripts

  
-  checking for empty or unneeded directories

  
-  building 'scLASER_0.0.0.9000.tar.gz'

  
   
```



<!-- rnb-output-end -->

<!-- rnb-output-begin eyJkYXRhIjoiSW5zdGFsbGluZyBwYWNrYWdlIGludG8gkUM6L1VzZXJzL0p1YW4vRG9jdW1lbnRzL1Ivd2luLWxpYnJhcnkvNC4xklxuKGFzIJFsaWKSIGlzIHVuc3BlY2lmaWVkKVxuIn0= -->

```
Installing package into 㤼㸱C:/Users/Juan/Documents/R/win-library/4.1㤼㸲
(as 㤼㸱lib㤼㸲 is unspecified)
```



<!-- rnb-output-end -->

<!-- rnb-output-begin eyJkYXRhIjoiKiBpbnN0YWxsaW5nICpzb3VyY2UqIHBhY2thZ2UgJ3NjTEFTRVInIC4uLlxuKiogdXNpbmcgc3RhZ2VkIGluc3RhbGxhdGlvblxuKiogUlxuKiogYnl0ZS1jb21waWxlIGFuZCBwcmVwYXJlIHBhY2thZ2UgZm9yIGxhenkgbG9hZGluZ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnZGF0YS50YWJsZTo6bGFzdCcgYnkgJ2RwbHlyOjpsYXN0JyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdkYXRhLnRhYmxlOjpmaXJzdCcgYnkgJ2RwbHlyOjpmaXJzdCcgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnTUFTUzo6c2VsZWN0JyBieSAnZHBseXI6OnNlbGVjdCcgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnZGF0YS50YWJsZTo6YmV0d2VlbicgYnkgJ2RwbHlyOjpiZXR3ZWVuJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdjYXJldDo6cGxzZGEnIGJ5ICdtaXhPbWljczo6cGxzZGEnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2NhcmV0OjpzcGxzZGEnIGJ5ICdtaXhPbWljczo6c3Bsc2RhJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdjYXJldDo6bmVhclplcm9WYXInIGJ5ICdtaXhPbWljczo6bmVhclplcm9WYXInIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2RwbHlyOjpjb2xsYXBzZScgYnkgJ25sbWU6OmNvbGxhcHNlJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdsbWU0OjpsbUxpc3QnIGJ5ICdubG1lOjpsbUxpc3QnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2ZvcmVhY2g6OndoZW4nIGJ5ICdwdXJycjo6d2hlbicgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnbWl4T21pY3M6Om1hcCcgYnkgJ3B1cnJyOjptYXAnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2NhcmV0OjpsaWZ0JyBieSAncHVycnI6OmxpZnQnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2RhdGEudGFibGU6OnRyYW5zcG9zZScgYnkgJ3B1cnJyOjp0cmFuc3Bvc2UnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2ZvcmVhY2g6OmFjY3VtdWxhdGUnIGJ5ICdwdXJycjo6YWNjdW11bGF0ZScgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnTWF0cml4Ojpjb3YyY29yJyBieSAnc3RhdHM6OmNvdjJjb3InIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2RwbHlyOjpmaWx0ZXInIGJ5ICdzdGF0czo6ZmlsdGVyJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdkcGx5cjo6bGFnJyBieSAnc3RhdHM6OmxhZycgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnTWF0cml4Ojp0b2VwbGl0eicgYnkgJ3N0YXRzOjp0b2VwbGl0eicgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnTWF0cml4Ojp1cGRhdGUnIGJ5ICdzdGF0czo6dXBkYXRlJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdNYXRyaXg6OmV4cGFuZCcgYnkgJ3RpZHlyOjpleHBhbmQnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ01hdHJpeDo6cGFjaycgYnkgJ3RpZHlyOjpwYWNrJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdNYXRyaXg6OnVucGFjaycgYnkgJ3RpZHlyOjp1bnBhY2snIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ01hdHJpeDo6dGFpbCcgYnkgJ3V0aWxzOjp0YWlsJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdNYXRyaXg6OmhlYWQnIGJ5ICd1dGlsczo6aGVhZCcgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuaW4gbWV0aG9kIGZvciAnbG9uZ0FuYWx5c2VzX3BsdXNNb2RlbFNlbGVjdGlvbicgd2l0aCBzaWduYXR1cmUgJ29iamVjdD1cInNjTEFTRVJcIic6IG5vIGRlZmluaXRpb24gZm9yIGNsYXNzIFwic2NMQVNFUlwiXG5pbiBtZXRob2QgZm9yICdzaG93JyB3aXRoIHNpZ25hdHVyZSAnXCJzY0xBU0VSXCInOiBubyBkZWZpbml0aW9uIGZvciBjbGFzcyBcInNjTEFTRVJcIlxuKiogaGVscFxuKioqIGluc3RhbGxpbmcgaGVscCBpbmRpY2VzXG4gIGNvbnZlcnRpbmcgaGVscCBmb3IgcGFja2FnZSAnc2NMQVNFUidcbiAgICBmaW5kaW5nIEhUTUwgbGlua3MgLi4uIGRvbmVcbiAgICBhc3NvY2lhdGlvbl9uYW1fc2NMQVNFUiAgICAgICAgICAgICAgICAgaHRtbCAgXG4gICAgY29tcHV0ZV91bWFwICAgICAgICAgICAgICAgICAgICAgICAgICAgIGh0bWwgIFxuICAgIGdlbmVyYXRlX2R1bW15X2RhdGEgICAgICAgICAgICAgICAgICAgICBodG1sICBcbiAgICBnZW5lcmF0ZV9kdW1teV9kYXRhXzJ0aW1lICAgICAgICAgICAgICAgaHRtbCAgXG4gICAgZ2VuZXJhdGVfcHNldWRvX3Bjc190aW1lICAgICAgICAgICAgICAgIGh0bWwgIFxuICAgIGxvbmdBbmFseXNlc19wbHVzTW9kZWxTZWxlY3Rpb24gICAgICAgICBodG1sICBcbiAgICBtaXhPbWljc19wbHNkYSAgICAgICAgICAgICAgICAgICAgICAgICAgaHRtbCAgXG4gICAgcGxvdF9jZWxsdHlwZV9wcm9wb3J0aW9ucyAgICAgICAgICAgICAgIGh0bWwgIFxuICAgIHJ1bl9oYXJtb255ICAgICAgICAgICAgICAgICAgICAgICAgICAgICBodG1sICBcbiAgICBzY0xBU0VSLWNsYXNzICAgICAgICAgICAgICAgICAgICAgICAgICAgaHRtbCAgXG4gICAgc2NMQVNFUiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGh0bWwgIFxuICAgIHNjTEFTRVJfbW9kZWxpbmcgICAgICAgICAgICAgICAgICAgICAgICBodG1sICBcbiAgICBzaG93LXNjTEFTRVItbWV0aG9kICAgICAgICAgICAgICAgICAgICAgaHRtbCAgXG4qKiBidWlsZGluZyBwYWNrYWdlIGluZGljZXNcbioqIGluc3RhbGxpbmcgdmlnbmV0dGVzXG4qKiB0ZXN0aW5nIGlmIGluc3RhbGxlZCBwYWNrYWdlIGNhbiBiZSBsb2FkZWQgZnJvbSB0ZW1wb3JhcnkgbG9jYXRpb25cbioqKiBhcmNoIC0gaTM4NlxuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnZGF0YS50YWJsZTo6bGFzdCcgYnkgJ2RwbHlyOjpsYXN0JyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdkYXRhLnRhYmxlOjpmaXJzdCcgYnkgJ2RwbHlyOjpmaXJzdCcgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnTUFTUzo6c2VsZWN0JyBieSAnZHBseXI6OnNlbGVjdCcgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnZGF0YS50YWJsZTo6YmV0d2VlbicgYnkgJ2RwbHlyOjpiZXR3ZWVuJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdjYXJldDo6cGxzZGEnIGJ5ICdtaXhPbWljczo6cGxzZGEnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2NhcmV0OjpzcGxzZGEnIGJ5ICdtaXhPbWljczo6c3Bsc2RhJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdjYXJldDo6bmVhclplcm9WYXInIGJ5ICdtaXhPbWljczo6bmVhclplcm9WYXInIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2RwbHlyOjpjb2xsYXBzZScgYnkgJ25sbWU6OmNvbGxhcHNlJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdsbWU0OjpsbUxpc3QnIGJ5ICdubG1lOjpsbUxpc3QnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2ZvcmVhY2g6OndoZW4nIGJ5ICdwdXJycjo6d2hlbicgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnbWl4T21pY3M6Om1hcCcgYnkgJ3B1cnJyOjptYXAnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2NhcmV0OjpsaWZ0JyBieSAncHVycnI6OmxpZnQnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2RhdGEudGFibGU6OnRyYW5zcG9zZScgYnkgJ3B1cnJyOjp0cmFuc3Bvc2UnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2ZvcmVhY2g6OmFjY3VtdWxhdGUnIGJ5ICdwdXJycjo6YWNjdW11bGF0ZScgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnTWF0cml4Ojpjb3YyY29yJyBieSAnc3RhdHM6OmNvdjJjb3InIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2RwbHlyOjpmaWx0ZXInIGJ5ICdzdGF0czo6ZmlsdGVyJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdkcGx5cjo6bGFnJyBieSAnc3RhdHM6OmxhZycgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnTWF0cml4Ojp0b2VwbGl0eicgYnkgJ3N0YXRzOjp0b2VwbGl0eicgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnTWF0cml4Ojp1cGRhdGUnIGJ5ICdzdGF0czo6dXBkYXRlJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdNYXRyaXg6OmV4cGFuZCcgYnkgJ3RpZHlyOjpleHBhbmQnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ01hdHJpeDo6cGFjaycgYnkgJ3RpZHlyOjpwYWNrJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdNYXRyaXg6OnVucGFjaycgYnkgJ3RpZHlyOjp1bnBhY2snIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ01hdHJpeDo6dGFpbCcgYnkgJ3V0aWxzOjp0YWlsJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdNYXRyaXg6OmhlYWQnIGJ5ICd1dGlsczo6aGVhZCcgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuKioqIGFyY2ggLSB4NjRcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2RhdGEudGFibGU6Omxhc3QnIGJ5ICdkcGx5cjo6bGFzdCcgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnZGF0YS50YWJsZTo6Zmlyc3QnIGJ5ICdkcGx5cjo6Zmlyc3QnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ01BU1M6OnNlbGVjdCcgYnkgJ2RwbHlyOjpzZWxlY3QnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2RhdGEudGFibGU6OmJldHdlZW4nIGJ5ICdkcGx5cjo6YmV0d2Vlbicgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnY2FyZXQ6OnBsc2RhJyBieSAnbWl4T21pY3M6OnBsc2RhJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdjYXJldDo6c3Bsc2RhJyBieSAnbWl4T21pY3M6OnNwbHNkYScgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnY2FyZXQ6Om5lYXJaZXJvVmFyJyBieSAnbWl4T21pY3M6Om5lYXJaZXJvVmFyJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdkcGx5cjo6Y29sbGFwc2UnIGJ5ICdubG1lOjpjb2xsYXBzZScgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnbG1lNDo6bG1MaXN0JyBieSAnbmxtZTo6bG1MaXN0JyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdmb3JlYWNoOjp3aGVuJyBieSAncHVycnI6OndoZW4nIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ21peE9taWNzOjptYXAnIGJ5ICdwdXJycjo6bWFwJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdjYXJldDo6bGlmdCcgYnkgJ3B1cnJyOjpsaWZ0JyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdkYXRhLnRhYmxlOjp0cmFuc3Bvc2UnIGJ5ICdwdXJycjo6dHJhbnNwb3NlJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdmb3JlYWNoOjphY2N1bXVsYXRlJyBieSAncHVycnI6OmFjY3VtdWxhdGUnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ01hdHJpeDo6Y292MmNvcicgYnkgJ3N0YXRzOjpjb3YyY29yJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdkcGx5cjo6ZmlsdGVyJyBieSAnc3RhdHM6OmZpbHRlcicgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnZHBseXI6OmxhZycgYnkgJ3N0YXRzOjpsYWcnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ01hdHJpeDo6dG9lcGxpdHonIGJ5ICdzdGF0czo6dG9lcGxpdHonIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ01hdHJpeDo6dXBkYXRlJyBieSAnc3RhdHM6OnVwZGF0ZScgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnTWF0cml4OjpleHBhbmQnIGJ5ICd0aWR5cjo6ZXhwYW5kJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdNYXRyaXg6OnBhY2snIGJ5ICd0aWR5cjo6cGFjaycgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnTWF0cml4Ojp1bnBhY2snIGJ5ICd0aWR5cjo6dW5wYWNrJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdNYXRyaXg6OnRhaWwnIGJ5ICd1dGlsczo6dGFpbCcgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnTWF0cml4OjpoZWFkJyBieSAndXRpbHM6OmhlYWQnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbioqIHRlc3RpbmcgaWYgaW5zdGFsbGVkIHBhY2thZ2UgY2FuIGJlIGxvYWRlZCBmcm9tIGZpbmFsIGxvY2F0aW9uXG4qKiogYXJjaCAtIGkzODZcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2RhdGEudGFibGU6Omxhc3QnIGJ5ICdkcGx5cjo6bGFzdCcgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnZGF0YS50YWJsZTo6Zmlyc3QnIGJ5ICdkcGx5cjo6Zmlyc3QnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ01BU1M6OnNlbGVjdCcgYnkgJ2RwbHlyOjpzZWxlY3QnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2RhdGEudGFibGU6OmJldHdlZW4nIGJ5ICdkcGx5cjo6YmV0d2Vlbicgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnY2FyZXQ6OnBsc2RhJyBieSAnbWl4T21pY3M6OnBsc2RhJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdjYXJldDo6c3Bsc2RhJyBieSAnbWl4T21pY3M6OnNwbHNkYScgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnY2FyZXQ6Om5lYXJaZXJvVmFyJyBieSAnbWl4T21pY3M6Om5lYXJaZXJvVmFyJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdkcGx5cjo6Y29sbGFwc2UnIGJ5ICdubG1lOjpjb2xsYXBzZScgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnbG1lNDo6bG1MaXN0JyBieSAnbmxtZTo6bG1MaXN0JyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdmb3JlYWNoOjp3aGVuJyBieSAncHVycnI6OndoZW4nIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ21peE9taWNzOjptYXAnIGJ5ICdwdXJycjo6bWFwJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdjYXJldDo6bGlmdCcgYnkgJ3B1cnJyOjpsaWZ0JyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdkYXRhLnRhYmxlOjp0cmFuc3Bvc2UnIGJ5ICdwdXJycjo6dHJhbnNwb3NlJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdmb3JlYWNoOjphY2N1bXVsYXRlJyBieSAncHVycnI6OmFjY3VtdWxhdGUnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ01hdHJpeDo6Y292MmNvcicgYnkgJ3N0YXRzOjpjb3YyY29yJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdkcGx5cjo6ZmlsdGVyJyBieSAnc3RhdHM6OmZpbHRlcicgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnZHBseXI6OmxhZycgYnkgJ3N0YXRzOjpsYWcnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ01hdHJpeDo6dG9lcGxpdHonIGJ5ICdzdGF0czo6dG9lcGxpdHonIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ01hdHJpeDo6dXBkYXRlJyBieSAnc3RhdHM6OnVwZGF0ZScgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnTWF0cml4OjpleHBhbmQnIGJ5ICd0aWR5cjo6ZXhwYW5kJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdNYXRyaXg6OnBhY2snIGJ5ICd0aWR5cjo6cGFjaycgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnTWF0cml4Ojp1bnBhY2snIGJ5ICd0aWR5cjo6dW5wYWNrJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdNYXRyaXg6OnRhaWwnIGJ5ICd1dGlsczo6dGFpbCcgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnTWF0cml4OjpoZWFkJyBieSAndXRpbHM6OmhlYWQnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbioqKiBhcmNoIC0geDY0XG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdkYXRhLnRhYmxlOjpsYXN0JyBieSAnZHBseXI6Omxhc3QnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2RhdGEudGFibGU6OmZpcnN0JyBieSAnZHBseXI6OmZpcnN0JyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdNQVNTOjpzZWxlY3QnIGJ5ICdkcGx5cjo6c2VsZWN0JyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdkYXRhLnRhYmxlOjpiZXR3ZWVuJyBieSAnZHBseXI6OmJldHdlZW4nIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2NhcmV0OjpwbHNkYScgYnkgJ21peE9taWNzOjpwbHNkYScgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnY2FyZXQ6OnNwbHNkYScgYnkgJ21peE9taWNzOjpzcGxzZGEnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2NhcmV0OjpuZWFyWmVyb1ZhcicgYnkgJ21peE9taWNzOjpuZWFyWmVyb1Zhcicgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnZHBseXI6OmNvbGxhcHNlJyBieSAnbmxtZTo6Y29sbGFwc2UnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2xtZTQ6OmxtTGlzdCcgYnkgJ25sbWU6OmxtTGlzdCcgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnZm9yZWFjaDo6d2hlbicgYnkgJ3B1cnJyOjp3aGVuJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdtaXhPbWljczo6bWFwJyBieSAncHVycnI6Om1hcCcgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnY2FyZXQ6OmxpZnQnIGJ5ICdwdXJycjo6bGlmdCcgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnZGF0YS50YWJsZTo6dHJhbnNwb3NlJyBieSAncHVycnI6OnRyYW5zcG9zZScgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnZm9yZWFjaDo6YWNjdW11bGF0ZScgYnkgJ3B1cnJyOjphY2N1bXVsYXRlJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdNYXRyaXg6OmNvdjJjb3InIGJ5ICdzdGF0czo6Y292MmNvcicgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnZHBseXI6OmZpbHRlcicgYnkgJ3N0YXRzOjpmaWx0ZXInIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ2RwbHlyOjpsYWcnIGJ5ICdzdGF0czo6bGFnJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdNYXRyaXg6OnRvZXBsaXR6JyBieSAnc3RhdHM6OnRvZXBsaXR6JyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG5XYXJuaW5nOiByZXBsYWNpbmcgcHJldmlvdXMgaW1wb3J0ICdNYXRyaXg6OnVwZGF0ZScgYnkgJ3N0YXRzOjp1cGRhdGUnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ01hdHJpeDo6ZXhwYW5kJyBieSAndGlkeXI6OmV4cGFuZCcgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnTWF0cml4OjpwYWNrJyBieSAndGlkeXI6OnBhY2snIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ01hdHJpeDo6dW5wYWNrJyBieSAndGlkeXI6OnVucGFjaycgd2hlbiBsb2FkaW5nICdzY0xBU0VSJ1xuV2FybmluZzogcmVwbGFjaW5nIHByZXZpb3VzIGltcG9ydCAnTWF0cml4Ojp0YWlsJyBieSAndXRpbHM6OnRhaWwnIHdoZW4gbG9hZGluZyAnc2NMQVNFUidcbldhcm5pbmc6IHJlcGxhY2luZyBwcmV2aW91cyBpbXBvcnQgJ01hdHJpeDo6aGVhZCcgYnkgJ3V0aWxzOjpoZWFkJyB3aGVuIGxvYWRpbmcgJ3NjTEFTRVInXG4qKiB0ZXN0aW5nIGlmIGluc3RhbGxlZCBwYWNrYWdlIGtlZXBzIGEgcmVjb3JkIG9mIHRlbXBvcmFyeSBpbnN0YWxsYXRpb24gcGF0aFxuKiBET05FIChzY0xBU0VSKVxuIn0= -->

```
* installing *source* package 'scLASER' ...
** using staged installation
** R
** byte-compile and prepare package for lazy loading
Warning: replacing previous import 'data.table::last' by 'dplyr::last' when loading 'scLASER'
Warning: replacing previous import 'data.table::first' by 'dplyr::first' when loading 'scLASER'
Warning: replacing previous import 'MASS::select' by 'dplyr::select' when loading 'scLASER'
Warning: replacing previous import 'data.table::between' by 'dplyr::between' when loading 'scLASER'
Warning: replacing previous import 'caret::plsda' by 'mixOmics::plsda' when loading 'scLASER'
Warning: replacing previous import 'caret::splsda' by 'mixOmics::splsda' when loading 'scLASER'
Warning: replacing previous import 'caret::nearZeroVar' by 'mixOmics::nearZeroVar' when loading 'scLASER'
Warning: replacing previous import 'dplyr::collapse' by 'nlme::collapse' when loading 'scLASER'
Warning: replacing previous import 'lme4::lmList' by 'nlme::lmList' when loading 'scLASER'
Warning: replacing previous import 'foreach::when' by 'purrr::when' when loading 'scLASER'
Warning: replacing previous import 'mixOmics::map' by 'purrr::map' when loading 'scLASER'
Warning: replacing previous import 'caret::lift' by 'purrr::lift' when loading 'scLASER'
Warning: replacing previous import 'data.table::transpose' by 'purrr::transpose' when loading 'scLASER'
Warning: replacing previous import 'foreach::accumulate' by 'purrr::accumulate' when loading 'scLASER'
Warning: replacing previous import 'Matrix::cov2cor' by 'stats::cov2cor' when loading 'scLASER'
Warning: replacing previous import 'dplyr::filter' by 'stats::filter' when loading 'scLASER'
Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading 'scLASER'
Warning: replacing previous import 'Matrix::toeplitz' by 'stats::toeplitz' when loading 'scLASER'
Warning: replacing previous import 'Matrix::update' by 'stats::update' when loading 'scLASER'
Warning: replacing previous import 'Matrix::expand' by 'tidyr::expand' when loading 'scLASER'
Warning: replacing previous import 'Matrix::pack' by 'tidyr::pack' when loading 'scLASER'
Warning: replacing previous import 'Matrix::unpack' by 'tidyr::unpack' when loading 'scLASER'
Warning: replacing previous import 'Matrix::tail' by 'utils::tail' when loading 'scLASER'
Warning: replacing previous import 'Matrix::head' by 'utils::head' when loading 'scLASER'
in method for 'longAnalyses_plusModelSelection' with signature 'object="scLASER"': no definition for class "scLASER"
in method for 'show' with signature '"scLASER"': no definition for class "scLASER"
** help
*** installing help indices
  converting help for package 'scLASER'
    finding HTML links ... done
    association_nam_scLASER                 html  
    compute_umap                            html  
    generate_dummy_data                     html  
    generate_dummy_data_2time               html  
    generate_pseudo_pcs_time                html  
    longAnalyses_plusModelSelection         html  
    mixOmics_plsda                          html  
    plot_celltype_proportions               html  
    run_harmony                             html  
    scLASER-class                           html  
    scLASER                                 html  
    scLASER_modeling                        html  
    show-scLASER-method                     html  
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
*** arch - i386
Warning: replacing previous import 'data.table::last' by 'dplyr::last' when loading 'scLASER'
Warning: replacing previous import 'data.table::first' by 'dplyr::first' when loading 'scLASER'
Warning: replacing previous import 'MASS::select' by 'dplyr::select' when loading 'scLASER'
Warning: replacing previous import 'data.table::between' by 'dplyr::between' when loading 'scLASER'
Warning: replacing previous import 'caret::plsda' by 'mixOmics::plsda' when loading 'scLASER'
Warning: replacing previous import 'caret::splsda' by 'mixOmics::splsda' when loading 'scLASER'
Warning: replacing previous import 'caret::nearZeroVar' by 'mixOmics::nearZeroVar' when loading 'scLASER'
Warning: replacing previous import 'dplyr::collapse' by 'nlme::collapse' when loading 'scLASER'
Warning: replacing previous import 'lme4::lmList' by 'nlme::lmList' when loading 'scLASER'
Warning: replacing previous import 'foreach::when' by 'purrr::when' when loading 'scLASER'
Warning: replacing previous import 'mixOmics::map' by 'purrr::map' when loading 'scLASER'
Warning: replacing previous import 'caret::lift' by 'purrr::lift' when loading 'scLASER'
Warning: replacing previous import 'data.table::transpose' by 'purrr::transpose' when loading 'scLASER'
Warning: replacing previous import 'foreach::accumulate' by 'purrr::accumulate' when loading 'scLASER'
Warning: replacing previous import 'Matrix::cov2cor' by 'stats::cov2cor' when loading 'scLASER'
Warning: replacing previous import 'dplyr::filter' by 'stats::filter' when loading 'scLASER'
Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading 'scLASER'
Warning: replacing previous import 'Matrix::toeplitz' by 'stats::toeplitz' when loading 'scLASER'
Warning: replacing previous import 'Matrix::update' by 'stats::update' when loading 'scLASER'
Warning: replacing previous import 'Matrix::expand' by 'tidyr::expand' when loading 'scLASER'
Warning: replacing previous import 'Matrix::pack' by 'tidyr::pack' when loading 'scLASER'
Warning: replacing previous import 'Matrix::unpack' by 'tidyr::unpack' when loading 'scLASER'
Warning: replacing previous import 'Matrix::tail' by 'utils::tail' when loading 'scLASER'
Warning: replacing previous import 'Matrix::head' by 'utils::head' when loading 'scLASER'
*** arch - x64
Warning: replacing previous import 'data.table::last' by 'dplyr::last' when loading 'scLASER'
Warning: replacing previous import 'data.table::first' by 'dplyr::first' when loading 'scLASER'
Warning: replacing previous import 'MASS::select' by 'dplyr::select' when loading 'scLASER'
Warning: replacing previous import 'data.table::between' by 'dplyr::between' when loading 'scLASER'
Warning: replacing previous import 'caret::plsda' by 'mixOmics::plsda' when loading 'scLASER'
Warning: replacing previous import 'caret::splsda' by 'mixOmics::splsda' when loading 'scLASER'
Warning: replacing previous import 'caret::nearZeroVar' by 'mixOmics::nearZeroVar' when loading 'scLASER'
Warning: replacing previous import 'dplyr::collapse' by 'nlme::collapse' when loading 'scLASER'
Warning: replacing previous import 'lme4::lmList' by 'nlme::lmList' when loading 'scLASER'
Warning: replacing previous import 'foreach::when' by 'purrr::when' when loading 'scLASER'
Warning: replacing previous import 'mixOmics::map' by 'purrr::map' when loading 'scLASER'
Warning: replacing previous import 'caret::lift' by 'purrr::lift' when loading 'scLASER'
Warning: replacing previous import 'data.table::transpose' by 'purrr::transpose' when loading 'scLASER'
Warning: replacing previous import 'foreach::accumulate' by 'purrr::accumulate' when loading 'scLASER'
Warning: replacing previous import 'Matrix::cov2cor' by 'stats::cov2cor' when loading 'scLASER'
Warning: replacing previous import 'dplyr::filter' by 'stats::filter' when loading 'scLASER'
Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading 'scLASER'
Warning: replacing previous import 'Matrix::toeplitz' by 'stats::toeplitz' when loading 'scLASER'
Warning: replacing previous import 'Matrix::update' by 'stats::update' when loading 'scLASER'
Warning: replacing previous import 'Matrix::expand' by 'tidyr::expand' when loading 'scLASER'
Warning: replacing previous import 'Matrix::pack' by 'tidyr::pack' when loading 'scLASER'
Warning: replacing previous import 'Matrix::unpack' by 'tidyr::unpack' when loading 'scLASER'
Warning: replacing previous import 'Matrix::tail' by 'utils::tail' when loading 'scLASER'
Warning: replacing previous import 'Matrix::head' by 'utils::head' when loading 'scLASER'
** testing if installed package can be loaded from final location
*** arch - i386
Warning: replacing previous import 'data.table::last' by 'dplyr::last' when loading 'scLASER'
Warning: replacing previous import 'data.table::first' by 'dplyr::first' when loading 'scLASER'
Warning: replacing previous import 'MASS::select' by 'dplyr::select' when loading 'scLASER'
Warning: replacing previous import 'data.table::between' by 'dplyr::between' when loading 'scLASER'
Warning: replacing previous import 'caret::plsda' by 'mixOmics::plsda' when loading 'scLASER'
Warning: replacing previous import 'caret::splsda' by 'mixOmics::splsda' when loading 'scLASER'
Warning: replacing previous import 'caret::nearZeroVar' by 'mixOmics::nearZeroVar' when loading 'scLASER'
Warning: replacing previous import 'dplyr::collapse' by 'nlme::collapse' when loading 'scLASER'
Warning: replacing previous import 'lme4::lmList' by 'nlme::lmList' when loading 'scLASER'
Warning: replacing previous import 'foreach::when' by 'purrr::when' when loading 'scLASER'
Warning: replacing previous import 'mixOmics::map' by 'purrr::map' when loading 'scLASER'
Warning: replacing previous import 'caret::lift' by 'purrr::lift' when loading 'scLASER'
Warning: replacing previous import 'data.table::transpose' by 'purrr::transpose' when loading 'scLASER'
Warning: replacing previous import 'foreach::accumulate' by 'purrr::accumulate' when loading 'scLASER'
Warning: replacing previous import 'Matrix::cov2cor' by 'stats::cov2cor' when loading 'scLASER'
Warning: replacing previous import 'dplyr::filter' by 'stats::filter' when loading 'scLASER'
Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading 'scLASER'
Warning: replacing previous import 'Matrix::toeplitz' by 'stats::toeplitz' when loading 'scLASER'
Warning: replacing previous import 'Matrix::update' by 'stats::update' when loading 'scLASER'
Warning: replacing previous import 'Matrix::expand' by 'tidyr::expand' when loading 'scLASER'
Warning: replacing previous import 'Matrix::pack' by 'tidyr::pack' when loading 'scLASER'
Warning: replacing previous import 'Matrix::unpack' by 'tidyr::unpack' when loading 'scLASER'
Warning: replacing previous import 'Matrix::tail' by 'utils::tail' when loading 'scLASER'
Warning: replacing previous import 'Matrix::head' by 'utils::head' when loading 'scLASER'
*** arch - x64
Warning: replacing previous import 'data.table::last' by 'dplyr::last' when loading 'scLASER'
Warning: replacing previous import 'data.table::first' by 'dplyr::first' when loading 'scLASER'
Warning: replacing previous import 'MASS::select' by 'dplyr::select' when loading 'scLASER'
Warning: replacing previous import 'data.table::between' by 'dplyr::between' when loading 'scLASER'
Warning: replacing previous import 'caret::plsda' by 'mixOmics::plsda' when loading 'scLASER'
Warning: replacing previous import 'caret::splsda' by 'mixOmics::splsda' when loading 'scLASER'
Warning: replacing previous import 'caret::nearZeroVar' by 'mixOmics::nearZeroVar' when loading 'scLASER'
Warning: replacing previous import 'dplyr::collapse' by 'nlme::collapse' when loading 'scLASER'
Warning: replacing previous import 'lme4::lmList' by 'nlme::lmList' when loading 'scLASER'
Warning: replacing previous import 'foreach::when' by 'purrr::when' when loading 'scLASER'
Warning: replacing previous import 'mixOmics::map' by 'purrr::map' when loading 'scLASER'
Warning: replacing previous import 'caret::lift' by 'purrr::lift' when loading 'scLASER'
Warning: replacing previous import 'data.table::transpose' by 'purrr::transpose' when loading 'scLASER'
Warning: replacing previous import 'foreach::accumulate' by 'purrr::accumulate' when loading 'scLASER'
Warning: replacing previous import 'Matrix::cov2cor' by 'stats::cov2cor' when loading 'scLASER'
Warning: replacing previous import 'dplyr::filter' by 'stats::filter' when loading 'scLASER'
Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading 'scLASER'
Warning: replacing previous import 'Matrix::toeplitz' by 'stats::toeplitz' when loading 'scLASER'
Warning: replacing previous import 'Matrix::update' by 'stats::update' when loading 'scLASER'
Warning: replacing previous import 'Matrix::expand' by 'tidyr::expand' when loading 'scLASER'
Warning: replacing previous import 'Matrix::pack' by 'tidyr::pack' when loading 'scLASER'
Warning: replacing previous import 'Matrix::unpack' by 'tidyr::unpack' when loading 'scLASER'
Warning: replacing previous import 'Matrix::tail' by 'utils::tail' when loading 'scLASER'
Warning: replacing previous import 'Matrix::head' by 'utils::head' when loading 'scLASER'
** testing if installed package keeps a record of temporary installation path
* DONE (scLASER)
```



<!-- rnb-output-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




## fixing the data simulation 


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-output-begin eyJkYXRhIjoiXG48IS0tIHJuYi1zb3VyY2UtYmVnaW4gZXlKa1lYUmhJam9pWUdCZ2NseHVYRzRqSUc5c1pDQm1kVzVqZEdsdmJseHVjMlYwTG5ObFpXUW9NakF5TUNsY2JuUmxjM1JmYjJ4a0lEd3RJSE5qVEVGVFJWSTZPbWRsYm1WeVlYUmxYMlIxYlcxNVgyUmhkR0VvWEc0Z0lHNWZZMlZzYkhNZ1BTQXhOVEFzWEc0Z0lITmtYMk5sYkd4MGVYQmxjeUE5SURBdU1TeGNiaUFnYmw5dFlXcHZjbDlqWld4c1gzUjVjR1Z6SUQwZ055eGNiaUFnYmw5dGFXNXZjbDlqWld4c1gzUjVjR1Z6SUQwZ015eGNiaUFnY21Wc1lYUnBkbVZmWVdKMWJtUmhibU5sSUQwZ01DNDBMRnh1SUNCdVgyMWhhbTl5WDJsdWRHVnlZV04wWDJObGJHeDBlWEJsY3lBOUlESXNYRzRnSUc1ZmJXbHViM0pmYVc1MFpYSmhZM1JmWTJWc2JIUjVjR1Z6SUQwZ01peGNiaUFnYmw5cGJtUnBkbWxrZFdGc2N5QTlJRE13TEZ4dUlDQnVYMkpoZEdOb2N5QTlJRFFzWEc0Z0lHbHVkR1Z5WVdOMGFXOXVYMlpsWVhSMWNtVWdQU0JjSW5acGMybDBYQ0lzWEc0Z0lIUnBiV1ZmY0c5cGJuUnpJRDBnTXl4Y2JpQWdkR1Z6ZEY5MllYSWdQU0JjSW1ScGMyVmhjMlZjSWl4Y2JpQWdjSEp2Y0Y5a2FYTmxZWE5sSUQwZ01DNDFMRnh1SUNCbVkxOXBiblJsY21GamRDQTlJREF1Tml4Y2JpQWdhVzUwWlhKaFkzUnBiMjVmZEhsd1pTQTlJRndpYzNCbFkybG1hV05jSWl4Y2JpQWdjMlZsWkNBOUlESXdNakFzWEc0Z0lIWnBjMmwwWDJWbVptVmpkSE5mY0hKdlozSmxjM052Y2lBOUlHTW9NQzR3TUN3Z01DNDJNQ3dnTUM0ek1Da3NYRzRnSUhacGMybDBYMlZtWm1WamRITmZZMjl1ZEhKdmJDQWdJQ0E5SUdNb01Dd2dNQ3dnTUNrc1hHNGdJR1JwY21WamRHbHZibDlpZVY5amJIVnpkR1Z5SUNBZ0lDQTlJR01vS3pFc0lDMHhMQ0FyTVN3Z0xURXBYRzRwWEc1Z1lHQWlmUT09IC0tPlxuXG5gYGByXG5cbiMgb2xkIGZ1bmN0aW9uXG5zZXQuc2VlZCgyMDIwKVxudGVzdF9vbGQgPC0gc2NMQVNFUjo6Z2VuZXJhdGVfZHVtbXlfZGF0YShcbiAgbl9jZWxscyA9IDE1MCxcbiAgc2RfY2VsbHR5cGVzID0gMC4xLFxuICBuX21ham9yX2NlbGxfdHlwZXMgPSA3LFxuICBuX21pbm9yX2NlbGxfdHlwZXMgPSAzLFxuICByZWxhdGl2ZV9hYnVuZGFuY2UgPSAwLjQsXG4gIG5fbWFqb3JfaW50ZXJhY3RfY2VsbHR5cGVzID0gMixcbiAgbl9taW5vcl9pbnRlcmFjdF9jZWxsdHlwZXMgPSAyLFxuICBuX2luZGl2aWR1YWxzID0gMzAsXG4gIG5fYmF0Y2hzID0gNCxcbiAgaW50ZXJhY3Rpb25fZmVhdHVyZSA9IFxcdmlzaXRcXCxcbiAgdGltZV9wb2ludHMgPSAzLFxuICB0ZXN0X3ZhciA9IFxcZGlzZWFzZVxcLFxuICBwcm9wX2Rpc2Vhc2UgPSAwLjUsXG4gIGZjX2ludGVyYWN0ID0gMC42LFxuICBpbnRlcmFjdGlvbl90eXBlID0gXFxzcGVjaWZpY1xcLFxuICBzZWVkID0gMjAyMCxcbiAgdmlzaXRfZWZmZWN0c19wcm9ncmVzc29yID0gYygwLjAwLCAwLjYwLCAwLjMwKSxcbiAgdmlzaXRfZWZmZWN0c19jb250cm9sICAgID0gYygwLCAwLCAwKSxcbiAgZGlyZWN0aW9uX2J5X2NsdXN0ZXIgICAgID0gYygrMSwgLTEsICsxLCAtMSlcbilcbmBgYFxuXG48IS0tIHJuYi1zb3VyY2UtZW5kIC0tPlxuIn0= -->

````

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuXG4jIG9sZCBmdW5jdGlvblxuc2V0LnNlZWQoMjAyMClcbnRlc3Rfb2xkIDwtIHNjTEFTRVI6OmdlbmVyYXRlX2R1bW15X2RhdGEoXG4gIG5fY2VsbHMgPSAxNTAsXG4gIHNkX2NlbGx0eXBlcyA9IDAuMSxcbiAgbl9tYWpvcl9jZWxsX3R5cGVzID0gNyxcbiAgbl9taW5vcl9jZWxsX3R5cGVzID0gMyxcbiAgcmVsYXRpdmVfYWJ1bmRhbmNlID0gMC40LFxuICBuX21ham9yX2ludGVyYWN0X2NlbGx0eXBlcyA9IDIsXG4gIG5fbWlub3JfaW50ZXJhY3RfY2VsbHR5cGVzID0gMixcbiAgbl9pbmRpdmlkdWFscyA9IDMwLFxuICBuX2JhdGNocyA9IDQsXG4gIGludGVyYWN0aW9uX2ZlYXR1cmUgPSBcInZpc2l0XCIsXG4gIHRpbWVfcG9pbnRzID0gMyxcbiAgdGVzdF92YXIgPSBcImRpc2Vhc2VcIixcbiAgcHJvcF9kaXNlYXNlID0gMC41LFxuICBmY19pbnRlcmFjdCA9IDAuNixcbiAgaW50ZXJhY3Rpb25fdHlwZSA9IFwic3BlY2lmaWNcIixcbiAgc2VlZCA9IDIwMjAsXG4gIHZpc2l0X2VmZmVjdHNfcHJvZ3Jlc3NvciA9IGMoMC4wMCwgMC42MCwgMC4zMCksXG4gIHZpc2l0X2VmZmVjdHNfY29udHJvbCAgICA9IGMoMCwgMCwgMCksXG4gIGRpcmVjdGlvbl9ieV9jbHVzdGVyICAgICA9IGMoKzEsIC0xLCArMSwgLTEpXG4pXG5gYGAifQ== -->

```r

# old function
set.seed(2020)
test_old <- scLASER::generate_dummy_data(
  n_cells = 150,
  sd_celltypes = 0.1,
  n_major_cell_types = 7,
  n_minor_cell_types = 3,
  relative_abundance = 0.4,
  n_major_interact_celltypes = 2,
  n_minor_interact_celltypes = 2,
  n_individuals = 30,
  n_batchs = 4,
  interaction_feature = \visit\,
  time_points = 3,
  test_var = \disease\,
  prop_disease = 0.5,
  fc_interact = 0.6,
  interaction_type = \specific\,
  seed = 2020,
  visit_effects_progressor = c(0.00, 0.60, 0.30),
  visit_effects_control    = c(0, 0, 0),
  direction_by_cluster     = c(+1, -1, +1, -1)
)
```

<!-- rnb-source-end -->
````



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYGBgclxuXG4jIG9sZCBmdW5jdGlvblxuc2V0LnNlZWQoMjAyMClcbnRlc3Rfb2xkIDwtIHNjTEFTRVI6OmdlbmVyYXRlX2R1bW15X2RhdGEoXG4gIG5fY2VsbHMgPSAxNTAsXG4gIHNkX2NlbGx0eXBlcyA9IDAuMSxcbiAgbl9tYWpvcl9jZWxsX3R5cGVzID0gNyxcbiAgbl9taW5vcl9jZWxsX3R5cGVzID0gMyxcbiAgcmVsYXRpdmVfYWJ1bmRhbmNlID0gMC40LFxuICBuX21ham9yX2ludGVyYWN0X2NlbGx0eXBlcyA9IDIsXG4gIG5fbWlub3JfaW50ZXJhY3RfY2VsbHR5cGVzID0gMixcbiAgbl9pbmRpdmlkdWFscyA9IDMwLFxuICBuX2JhdGNocyA9IDQsXG4gIGludGVyYWN0aW9uX2ZlYXR1cmUgPSBcXHZpc2l0XFwsXG4gIHRpbWVfcG9pbnRzID0gMyxcbiAgdGVzdF92YXIgPSBcXGRpc2Vhc2VcXCxcbiAgcHJvcF9kaXNlYXNlID0gMC41LFxuICBmY19pbnRlcmFjdCA9IDAuNixcbiAgaW50ZXJhY3Rpb25fdHlwZSA9IFxcc3BlY2lmaWNcXCxcbiAgc2VlZCA9IDIwMjAsXG4gIHZpc2l0X2VmZmVjdHNfcHJvZ3Jlc3NvciA9IGMoMC4wMCwgMC42MCwgMC4zMCksXG4gIHZpc2l0X2VmZmVjdHNfY29udHJvbCAgICA9IGMoMCwgMCwgMCksXG4gIGRpcmVjdGlvbl9ieV9jbHVzdGVyICAgICA9IGMoKzEsIC0xLCArMSwgLTEpXG4pXG5gYGBcbmBgYCJ9 -->

```r
```r

# old function
set.seed(2020)
test_old <- scLASER::generate_dummy_data(
  n_cells = 150,
  sd_celltypes = 0.1,
  n_major_cell_types = 7,
  n_minor_cell_types = 3,
  relative_abundance = 0.4,
  n_major_interact_celltypes = 2,
  n_minor_interact_celltypes = 2,
  n_individuals = 30,
  n_batchs = 4,
  interaction_feature = \visit\,
  time_points = 3,
  test_var = \disease\,
  prop_disease = 0.5,
  fc_interact = 0.6,
  interaction_type = \specific\,
  seed = 2020,
  visit_effects_progressor = c(0.00, 0.60, 0.30),
  visit_effects_control    = c(0, 0, 0),
  direction_by_cluster     = c(+1, -1, +1, -1)
)

<!-- rnb-source-end -->


<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-output-begin eyJkYXRhIjoiXG48IS0tIHJuYi1zb3VyY2UtYmVnaW4gZXlKa1lYUmhJam9pWUdCZ2NseHVJMjVsZHlCbWRXNWpkR2x2Ymx4dVhHNGpKeUJVYVhSc1pWeHVJeWRjYmlNbklFQndZWEpoYlNCdVgyTmxiR3h6WEc0akp5QkFjR0Z5WVcwZ2MyUmZZMlZzYkhSNWNHVnpYRzRqSnlCQWNHRnlZVzBnYmw5dFlXcHZjbDlqWld4c1gzUjVjR1Z6WEc0akp5QkFjR0Z5WVcwZ2JsOXRhVzV2Y2w5alpXeHNYM1I1Y0dWelhHNGpKeUJBY0dGeVlXMGdjbVZzWVhScGRtVmZZV0oxYm1SaGJtTmxYRzRqSnlCQWNHRnlZVzBnYmw5dFlXcHZjbDlwYm5SbGNtRmpkRjlqWld4c2RIbHdaWE5jYmlNbklFQndZWEpoYlNCdVgyMXBibTl5WDJsdWRHVnlZV04wWDJObGJHeDBlWEJsYzF4dUl5Y2dRSEJoY21GdElHNWZhVzVrYVhacFpIVmhiSE5jYmlNbklFQndZWEpoYlNCdVgySmhkR05vYzF4dUl5Y2dRSEJoY21GdElHbHVkR1Z5WVdOMGFXOXVYMlpsWVhSMWNtVmNiaU1uSUVCd1lYSmhiU0IwYVcxbFgzQnZhVzUwYzF4dUl5Y2dRSEJoY21GdElIUmxjM1JmZG1GeVhHNGpKeUJBY0dGeVlXMGdjSEp2Y0Y5a2FYTmxZWE5sWEc0akp5QkFjR0Z5WVcwZ1ptTmZhVzUwWlhKaFkzUmNiaU1uSUVCd1lYSmhiU0JwYm5SbGNtRmpkR2x2Ymw5MGVYQmxYRzRqSnlCQWNHRnlZVzBnYzJWbFpGeHVJeWNnUUhCaGNtRnRJSFpwYzJsMFgyVm1abVZqZEhOZmNISnZaM0psYzNOdmNseHVJeWNnUUhCaGNtRnRJSFpwYzJsMFgyVm1abVZqZEhOZlkyOXVkSEp2YkZ4dUl5Y2dRSEJoY21GdElHUnBjbVZqZEdsdmJsOWllVjlqYkhWemRHVnlYRzRqSjF4dUl5Y2dRSEpsZEhWeWJseHVJeWNnUUdWNGNHOXlkRnh1SXlkY2JpTW5JRUJsZUdGdGNHeGxjMXh1WjJWdVpYSmhkR1ZmWkhWdGJYbGZaR0YwWVY5dVpYY2dQQzBnWm5WdVkzUnBiMjRvWEc0Z0lHNWZZMlZzYkhNZ1BTQXpNREF3TENBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNNZ1ltRnpaV3hwYm1VZ1kyVnNiSE1nY0dWeUlHMWhhbTl5SUdObGJHd2dkSGx3WlNCd1pYSWdjMkZ0Y0d4bFhHNGdJSE5rWDJObGJHeDBlWEJsY3lBOUlEQXVNVEFzSUNBZ0lDQWdJQ0FnSUNBZ0lDTWdjbVZzWVhScGRtVWdjMlFnWm05eUlHTnZkVzUwYzF4dUlDQnVYMjFoYW05eVgyTmxiR3hmZEhsd1pYTWdQU0EzTEZ4dUlDQnVYMjFwYm05eVgyTmxiR3hmZEhsd1pYTWdQU0F6TEZ4dUlDQnlaV3hoZEdsMlpWOWhZblZ1WkdGdVkyVWdQU0F3TGpFd0xDQWdJQ0FnSUNBaklHMXBibTl5SUhaeklHMWhhbTl5SUdKaGMyVnNhVzVsSUhKaGRHbHZYRzRnSUc1ZmJXRnFiM0pmYVc1MFpYSmhZM1JmWTJWc2JIUjVjR1Z6SUQwZ01Td2dJQ01nYUc5M0lHMWhibmtnYldGcWIzSnpJR0Z5WlNCY0ltbHVkR1Z5WVdOMGFXNW5YQ0pjYmlBZ2JsOXRhVzV2Y2w5cGJuUmxjbUZqZEY5alpXeHNkSGx3WlhNZ1BTQXhMQ0FnSXlCb2IzY2diV0Z1ZVNCdGFXNXZjbk1nWVhKbElGd2lhVzUwWlhKaFkzUnBibWRjSWx4dUlDQnVYMmx1WkdsMmFXUjFZV3h6SUQwZ016QXNYRzRnSUc1ZlltRjBZMmh6SUQwZ05DeGNibHh1SUNCcGJuUmxjbUZqZEdsdmJsOW1aV0YwZFhKbElEMGdYQ0oyYVhOcGRGd2lMQ0FnSUNNZ2EyVndkQ0JtYjNJZ2JHRmlaV3hwYm1kY2JpQWdkR2x0WlY5d2IybHVkSE1nUFNBMExDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0l5QStQU0F5T3lCM2IzSnJjeUJtYjNJZ015dGNiaUFnZEdWemRGOTJZWElnUFNCY0ltUnBjMlZoYzJWY0lpeGNiaUFnY0hKdmNGOWthWE5sWVhObElEMGdNQzQxTUN4Y2JseHVJQ0JtWTE5cGJuUmxjbUZqZENBOUlEQXVNVEFzSUNBZ0lDQWdJQ0FnSUNBZ0lDQWpJR1ZtWm1WamRDQnRZV2R1YVhSMVpHVWdkWE5sWkNCaWVTQmtaV1poZFd4MGN5QmlaV3h2ZDF4dUlDQnBiblJsY21GamRHbHZibDkwZVhCbElEMGdZeWhjSW5Od1pXTnBabWxqWENJc1hDSmthV1ptWlhKbGJuUnBZV3hjSWl4Y0ltOXdjRzl6YVhSbFhDSXBMRnh1SUNCelpXVmtJRDBnTVRJek5DeGNibHh1SUNCMmFYTnBkRjlsWm1abFkzUnpYM0J5YjJkeVpYTnpiM0lnUFNCT1ZVeE1MQ0FqSUcxMWJIUnBjR3hwWTJGMGFYWmxJR05vWVc1blpUb2dLekF1TVNCdFpXRnVjeUFyTVRBbElIWnpJR0poYzJWc2FXNWxYRzRnSUhacGMybDBYMlZtWm1WamRITmZZMjl1ZEhKdmJDQWdJQ0E5SUU1VlRFd3NJQ01nWkdWbVlYVnNkQ0F3YzF4dUlDQmthWEpsWTNScGIyNWZZbmxmWTJ4MWMzUmxjaUFnSUNBZ1BTQk9WVXhNSUNBaklDc3hMeTB4SUdadmNpQnBiblJsY21GamRHbHVaeUJqYkhWemRHVnljeXdnY21WamVXTnNaV1FnWVhNZ2JtVmxaR1ZrWEc0cElIdGNiaUFnYzJWMExuTmxaV1FvYzJWbFpDbGNiaUFnYVc1MFpYSmhZM1JwYjI1ZmRIbHdaU0E4TFNCdFlYUmphQzVoY21jb2FXNTBaWEpoWTNScGIyNWZkSGx3WlNsY2JseHVJQ0FqSUMwdExTQkNZWE5wWXlCelpYUjFjQ0F0TFMwdFhHNGdJRzVmWTJWc2JGOTBlWEJsY3lBOExTQnVYMjFoYW05eVgyTmxiR3hmZEhsd1pYTWdLeUJ1WDIxcGJtOXlYMk5sYkd4ZmRIbHdaWE5jYmlBZ2MzUnZjR2xtYm05MEtIUnBiV1ZmY0c5cGJuUnpJRDQ5SURJc0lHNWZZMlZzYkY5MGVYQmxjeUErUFNBeEtWeHVYRzRnSUdObGJHeGZkSGx3WlhNZ1BDMGdURVZVVkVWU1UxdHpaWEZmYkdWdUtHNWZZMlZzYkY5MGVYQmxjeWxkWEc0Z0lHMWhhbTl5WDJsa2VDQWdQQzBnYzJWeFgyeGxiaWh1WDIxaGFtOXlYMk5sYkd4ZmRIbHdaWE1wWEc0Z0lHMXBibTl5WDJsa2VDQWdQQzBnYVdZZ0tHNWZiV2x1YjNKZlkyVnNiRjkwZVhCbGN5QStJREFwSUNodVgyMWhhbTl5WDJObGJHeGZkSGx3WlhNZ0t5QnpaWEZmYkdWdUtHNWZiV2x1YjNKZlkyVnNiRjkwZVhCbGN5a3BJR1ZzYzJVZ2FXNTBaV2RsY2lnd0tWeHVYRzRnSUNNZ2QyaHBZMmdnWTJWc2JDQjBlWEJsY3lCaGNtVWdYQ0pwYm5SbGNtRmpkR2x1WjF3aUlDaG1hWEp6ZENCemIyMWxJRzFoYW05eWN5d2diR0Z6ZENCemIyMWxJRzFwYm05eWN5bGNiaUFnYVc1MFpYSmhZM1JmYVdSNElEd3RJR01vWEc0Z0lDQWdhR1ZoWkNodFlXcHZjbDlwWkhnc0lHNWZiV0ZxYjNKZmFXNTBaWEpoWTNSZlkyVnNiSFI1Y0dWektTeGNiaUFnSUNCMFlXbHNLSE5sY1Y5c1pXNG9ibDlqWld4c1gzUjVjR1Z6S1N3Z2JsOXRhVzV2Y2w5cGJuUmxjbUZqZEY5alpXeHNkSGx3WlhNcFhHNGdJQ2xjYmlBZ2FXNTBaWEpoWTNSZmFXUjRJRHd0SUdsdWRHVnljMlZqZENocGJuUmxjbUZqZEY5cFpIZ3NJSE5sY1Y5c1pXNG9ibDlqWld4c1gzUjVjR1Z6S1NrZ0lDTWdaM1ZoY21RZ2NtRnBiSE5jYmlBZ2FXNTBaWEpoWTNSZlkyVnNiRjkwZVhCbGN5QThMU0JqWld4c1gzUjVjR1Z6VzJsdWRHVnlZV04wWDJsa2VGMWNibHh1SUNBaklDMHRMU0JUZFdKcVpXTjBjeUJoYm1RZ2RtbHphWFJ6SUMwdExTMWNiaUFnYzNWaWFtVmpkRjlwWkNBOExTQndZWE4wWlRBb1hDSlRWVUpmWENJc0lITmxjVjlzWlc0b2JsOXBibVJwZG1sa2RXRnNjeWtwWEc0Z0lHUnBjMlZoYzJWZmRtVmpJRHd0SUdNb2NtVndLREZNTENCeWIzVnVaQ2h1WDJsdVpHbDJhV1IxWVd4eklDb2djSEp2Y0Y5a2FYTmxZWE5sS1Nrc1hHNGdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdjbVZ3S0RCTUxDQnVYMmx1WkdsMmFXUjFZV3h6SUMwZ2NtOTFibVFvYmw5cGJtUnBkbWxrZFdGc2N5QXFJSEJ5YjNCZlpHbHpaV0Z6WlNrcEtWeHVJQ0JrYVhObFlYTmxYM1psWXlBOExTQnpZVzF3YkdVb1pHbHpaV0Z6WlY5MlpXTXNJRzVmYVc1a2FYWnBaSFZoYkhNcFhHNWNiaUFnYzJWNFgzWmxZeUE4TFNCellXMXdiR1VvWXlnd1RDd2dNVXdwTENCdVgybHVaR2wyYVdSMVlXeHpMQ0J5WlhCc1lXTmxJRDBnVkZKVlJTa2dJQ0FnSUNBZ0lDQWdJeUF3THpGY2JpQWdZV2RsWDNabFl5QThMU0J6WVcxd2JHVW9NVGc2TmpBc0lHNWZhVzVrYVhacFpIVmhiSE1zSUhKbGNHeGhZMlVnUFNCVVVsVkZLVnh1SUNCaWJXbGZkbVZqSUR3dElITmhiWEJzWlNneE5Ub3pOU3dnYmw5cGJtUnBkbWxrZFdGc2N5d2djbVZ3YkdGalpTQTlJRlJTVlVVcFhHNGdJR0poZEdOb1gzWmxZeUE4TFNCeVpYQW9jMlZ4WDJ4bGJpaHVYMkpoZEdOb2N5a3NJR3hsYm1kMGFDNXZkWFFnUFNCdVgybHVaR2wyYVdSMVlXeHpLVnh1WEc0Z0lITjFZbXBsWTNSeklEd3RJR1JoZEdFdVpuSmhiV1VvWEc0Z0lDQWdjM1ZpYW1WamRGOXBaQ0E5SUhOMVltcGxZM1JmYVdRc1hHNGdJQ0FnYzJWNElDQWdQU0J6WlhoZmRtVmpMRnh1SUNBZ0lHUnBjMlZoYzJVZ1BTQmthWE5sWVhObFgzWmxZeXdnSUNBZ0lDQWdJQ0FnSUNNZ1kyRnViMjVwWTJGc0lHVjRjRzl6ZFhKbElITjBiM0poWjJWY2JpQWdJQ0JoWjJVZ0lDQTlJR0ZuWlY5MlpXTXNYRzRnSUNBZ1ltMXBJQ0FnUFNCaWJXbGZkbVZqTEZ4dUlDQWdJR0poZEdOb0lEMGdabUZqZEc5eUtHSmhkR05vWDNabFl5a3NYRzRnSUNBZ2MzUnlhVzVuYzBGelJtRmpkRzl5Y3lBOUlFWkJURk5GWEc0Z0lDbGNibHh1SUNCMmFYTnBkSE1nUEMwZ1pHRjBZUzVtY21GdFpTaGNiaUFnSUNCemRXSnFaV04wWDJsa0lEMGdjbVZ3S0hOMVltcGxZM1JmYVdRc0lHVmhZMmdnUFNCMGFXMWxYM0J2YVc1MGN5a3NYRzRnSUNBZ2RtbHphWFFnSUNBZ0lDQTlJSEpsY0Nnd09paDBhVzFsWDNCdmFXNTBjeUF0SURFcExDQjBhVzFsY3lBOUlHNWZhVzVrYVhacFpIVmhiSE1wTEZ4dUlDQWdJSE4wY21sdVozTkJjMFpoWTNSdmNuTWdQU0JHUVV4VFJWeHVJQ0FwWEc0Z0lIWnBjMmwwY3lSellXMXdiR1ZmYVdRZ1BDMGdjR0Z6ZEdVd0tIWnBjMmwwY3lSemRXSnFaV04wWDJsa0xDQmNJbDlXWENJc0lIWnBjMmwwY3lSMmFYTnBkQ2xjYmx4dUlDQnRaWFJoSUR3dElHMWxjbWRsS0hacGMybDBjeXdnYzNWaWFtVmpkSE1zSUdKNUlEMGdYQ0p6ZFdKcVpXTjBYMmxrWENJc0lITnZjblFnUFNCR1FVeFRSU2xjYmx4dUlDQWpJRTFoYTJVZ2MzVnlaU0JoSUdOdmJIVnRiaUJ1WVcxbFpDQmdkR1Z6ZEY5MllYSmdJR1Y0YVhOMGN5QW9aWFpsYmlCcFppQjBaWE4wWDNaaGNpQWhQU0JjSW1ScGMyVmhjMlZjSWlsY2JpQWdhV1lnS0NGcFpHVnVkR2xqWVd3b2RHVnpkRjkyWVhJc0lGd2laR2x6WldGelpWd2lLU2tnZTF4dUlDQWdJRzFsZEdGYlczUmxjM1JmZG1GeVhWMGdQQzBnYldWMFlWdGJYQ0prYVhObFlYTmxYQ0pkWFZ4dUlDQjlYRzVjYmlBZ0l5Qk1ZV0psYkNCaGJtUWdTVTVVUlZKQlExUkpUMDRnVkVWU1RTQW9jR1Z5YzJsemRDQjBieUJ2ZFhSd2RYUXBYRzRnSUcxbGRHRWthVzUwWlhKaFkzUnBiMjRnSUNBZ1BDMGdjR0Z6ZEdVd0tHbHVkR1Z5WVdOMGFXOXVYMlpsWVhSMWNtVXNJRndpT2x3aUxDQjBaWE4wWDNaaGNpbGNiaUFnYldWMFlTUnBiblJsY21GamRGOTBaWEp0SUNBOExTQmhjeTVwYm5SbFoyVnlLRzFsZEdGYlcybHVkR1Z5WVdOMGFXOXVYMlpsWVhSMWNtVmRYU2tnS2lCaGN5NXBiblJsWjJWeUtHMWxkR0ZiVzNSbGMzUmZkbUZ5WFYwcFhHNWNiaUFnSXlBdExTMGdSR1ZtWVhWc2RDQjJhWE5wZENCbFptWmxZM1J6SUNoc1pXNW5kR2dnUFNCMGFXMWxYM0J2YVc1MGN5a2dMUzB0TFZ4dUlDQnBaaUFvYVhNdWJuVnNiQ2gyYVhOcGRGOWxabVpsWTNSelgzQnliMmR5WlhOemIzSXBLU0I3WEc0Z0lDQWdhV1lnS0dsdWRHVnlZV04wYVc5dVgzUjVjR1VnUFQwZ1hDSnpjR1ZqYVdacFkxd2lLU0I3WEc0Z0lDQWdJQ0IyWlNBOExTQnlaWEFvTUN3Z2RHbHRaVjl3YjJsdWRITXBYRzRnSUNBZ0lDQnBaaUFvZEdsdFpWOXdiMmx1ZEhNZ1BqMGdNaWtnZG1WYk1sMGdQQzBnWm1OZmFXNTBaWEpoWTNRZ0lDTWdZblZ0Y0NCV01TQnZibXg1WEc0Z0lDQWdJQ0IyYVhOcGRGOWxabVpsWTNSelgzQnliMmR5WlhOemIzSWdQQzBnZG1WY2JpQWdJQ0I5SUdWc2MyVWdhV1lnS0dsdWRHVnlZV04wYVc5dVgzUjVjR1VnUFQwZ1hDSmthV1ptWlhKbGJuUnBZV3hjSWlrZ2UxeHVJQ0FnSUNBZ2RtbHphWFJmWldabVpXTjBjMTl3Y205bmNtVnpjMjl5SUR3dElISmxjQ2htWTE5cGJuUmxjbUZqZEN3Z2RHbHRaVjl3YjJsdWRITXBYRzRnSUNBZ2ZTQmxiSE5sSUhzZ0l5QmNJbTl3Y0c5emFYUmxYQ0k2SUdGc2RHVnlibUYwWlNBckx5MGdjM1JoY25ScGJtY2dZWFFnVmpCY2JpQWdJQ0FnSUhabElEd3RJSEpsY0Nnd0xDQjBhVzFsWDNCdmFXNTBjeWxjYmlBZ0lDQWdJSFpsVzNObGNTZ3hMQ0IwYVcxbFgzQnZhVzUwY3l3Z1lua2dQU0F5S1YwZ1BDMGdLMlpqWDJsdWRHVnlZV04wSUNBaklGWXdMQ0JXTWl3Z0xpNHVYRzRnSUNBZ0lDQnBaaUFvZEdsdFpWOXdiMmx1ZEhNZ1BqMGdNaWtnZG1WYmMyVnhLRElzSUhScGJXVmZjRzlwYm5SekxDQmllU0E5SURJcFhTQThMU0F0Wm1OZmFXNTBaWEpoWTNRZ0l5QldNU3dnVmpNc0lDNHVMbHh1SUNBZ0lDQWdkbWx6YVhSZlpXWm1aV04wYzE5d2NtOW5jbVZ6YzI5eUlEd3RJSFpsWEc0Z0lDQWdmVnh1SUNCOVhHNGdJR2xtSUNocGN5NXVkV3hzS0hacGMybDBYMlZtWm1WamRITmZZMjl1ZEhKdmJDa3BJSHRjYmlBZ0lDQjJhWE5wZEY5bFptWmxZM1J6WDJOdmJuUnliMndnUEMwZ2NtVndLREFzSUhScGJXVmZjRzlwYm5SektWeHVJQ0I5WEc0Z0lITjBiM0JwWm01dmRDaHNaVzVuZEdnb2RtbHphWFJmWldabVpXTjBjMTl3Y205bmNtVnpjMjl5S1NBOVBTQjBhVzFsWDNCdmFXNTBjeXhjYmlBZ0lDQWdJQ0FnSUNBZ0lHeGxibWQwYUNoMmFYTnBkRjlsWm1abFkzUnpYMk52Ym5SeWIyd3BJQ0FnSUQwOUlIUnBiV1ZmY0c5cGJuUnpLVnh1WEc0Z0lDTWdSR2x5WldOMGFXOXVJSEJsY2lCcGJuUmxjbUZqZEdsdVp5QmpiSFZ6ZEdWeVhHNGdJR2xtSUNocGN5NXVkV3hzS0dScGNtVmpkR2x2Ymw5aWVWOWpiSFZ6ZEdWeUtTa2dlMXh1SUNBZ0lDTWdiMnhrSUdSbFptRjFiSFE2SUdGc2JDQnBiblJsY21GamRHbHVaeUJqWld4c0lIUjVjR1Z6SUdkdklHbHVJSFJvWlNCellXMWxJR1JwY21WamRHbHZiaUFvS3pFcFhHNGdJQ0FnWkdseVpXTjBhVzl1WDJKNVgyTnNkWE4wWlhJZ1BDMGdjbVZ3S0RGTUxDQnRZWGdvTVN3Z2JHVnVaM1JvS0dsdWRHVnlZV04wWDJObGJHeGZkSGx3WlhNcEtTbGNiaUFnSUNCdVlXMWxjeWhrYVhKbFkzUnBiMjVmWW5sZlkyeDFjM1JsY2lrZ1BDMGdhVzUwWlhKaFkzUmZZMlZzYkY5MGVYQmxjMXh1SUNCOUlHVnNjMlVnZTF4dUlDQWdJR2xtSUNnaGFYTXViblZzYkNodVlXMWxjeWhrYVhKbFkzUnBiMjVmWW5sZlkyeDFjM1JsY2lrcEtTQjdYRzRnSUNBZ0lDQWpJRTVGVnpvZ2RYTmxjaUJ3WVhOelpXUWdZU0FxYm1GdFpXUXFJSFpsWTNSdmNpd2daUzVuTGlCaktFRWdQU0F4TENCRElEMGdMVEVzSUVvZ1BTQXhLVnh1SUNBZ0lDQWdJeUJYWlNCaGJHbG5iaUIwYnlCcGJuUmxjbUZqZEY5alpXeHNYM1I1Y0dWeklHRnVaQ0JrWldaaGRXeDBJRzFwYzNOcGJtY2diMjVsY3lCMGJ5QXJNUzVjYmlBZ0lDQWdJSFJ0Y0NBOExTQnlaWEFvTVV3c0lHMWhlQ2d4TENCc1pXNW5kR2dvYVc1MFpYSmhZM1JmWTJWc2JGOTBlWEJsY3lrcEtWeHVJQ0FnSUNBZ2JtRnRaWE1vZEcxd0tTQThMU0JwYm5SbGNtRmpkRjlqWld4c1gzUjVjR1Z6WEc1Y2JpQWdJQ0FnSUcxaGRHTm9aV1FnUEMwZ2FXNTBaWEp6WldOMEtHNWhiV1Z6S0dScGNtVmpkR2x2Ymw5aWVWOWpiSFZ6ZEdWeUtTd2dhVzUwWlhKaFkzUmZZMlZzYkY5MGVYQmxjeWxjYmlBZ0lDQWdJSFJ0Y0Z0dFlYUmphR1ZrWFNBOExTQmthWEpsWTNScGIyNWZZbmxmWTJ4MWMzUmxjbHR0WVhSamFHVmtYVnh1WEc0Z0lDQWdJQ0JrYVhKbFkzUnBiMjVmWW5sZlkyeDFjM1JsY2lBOExTQjBiWEJjYmlBZ0lDQjlJR1ZzYzJVZ2UxeHVJQ0FnSUNBZ0l5QnZiR1FnWW1Wb1lYWnBiM0lnWm05eUlIVnVibUZ0WldRZ2RtVmpkRzl5T2lCeVpXTjVZMnhsSUc5MlpYSWdhVzUwWlhKaFkzUnBibWNnWTJWc2JDQjBlWEJsYzF4dUlDQWdJQ0FnWkdseVpXTjBhVzl1WDJKNVgyTnNkWE4wWlhJZ1BDMGdjbVZ3S0dScGNtVmpkR2x2Ymw5aWVWOWpiSFZ6ZEdWeUxGeHVJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lHeGxibWQwYUM1dmRYUWdQU0JzWlc1bmRHZ29hVzUwWlhKaFkzUmZZMlZzYkY5MGVYQmxjeWtwWEc0Z0lDQWdJQ0J1WVcxbGN5aGthWEpsWTNScGIyNWZZbmxmWTJ4MWMzUmxjaWtnUEMwZ2FXNTBaWEpoWTNSZlkyVnNiRjkwZVhCbGMxeHVJQ0FnSUgxY2JpQWdmVnh1WEc0Z0lDTWdMUzB0SUVKaGMyVnNhVzVsSUdOdmRXNTBjeUJ3WlhJZ0tITmhiWEJzWlN3Z1kyVnNiRjkwZVhCbEtTQXRMUzB0WEc0Z0lHOXVaVjl6WVcxd2JHVmZZMjkxYm5SeklEd3RJR1oxYm1OMGFXOXVLQ2tnZTF4dUlDQWdJRzFoYW05eVgyTnZkVzUwY3lBOExTQnliM1Z1WkNoeWRXNXBaaWh1WDIxaGFtOXlYMk5sYkd4ZmRIbHdaWE1zWEc0Z0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJRzFwYmlBOUlHNWZZMlZzYkhNZ0tpQW9NU0F0SUhOa1gyTmxiR3gwZVhCbGN5a3NYRzRnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lHMWhlQ0E5SUc1ZlkyVnNiSE1nS2lBb01TQXJJSE5rWDJObGJHeDBlWEJsY3lrcEtWeHVJQ0FnSUcxcGJtOXlYMk52ZFc1MGN5QThMU0JwWmlBb2JsOXRhVzV2Y2w5alpXeHNYM1I1Y0dWeklENGdNQ2tnZTF4dUlDQWdJQ0FnY205MWJtUW9jblZ1YVdZb2JsOXRhVzV2Y2w5alpXeHNYM1I1Y0dWekxGeHVJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdiV2x1SUQwZ2JsOWpaV3hzY3lBcUlISmxiR0YwYVhabFgyRmlkVzVrWVc1alpTQXFJQ2d4SUMwZ2MyUmZZMlZzYkhSNWNHVnpLU3hjYmlBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUcxaGVDQTlJRzVmWTJWc2JITWdLaUJ5Wld4aGRHbDJaVjloWW5WdVpHRnVZMlVnS2lBb01TQXJJSE5rWDJObGJHeDBlWEJsY3lrcEtWeHVJQ0FnSUgwZ1pXeHpaU0JwYm5SbFoyVnlLREFwWEc0Z0lDQWdZeWh0WVdwdmNsOWpiM1Z1ZEhNc0lHMXBibTl5WDJOdmRXNTBjeWxjYmlBZ2ZWeHVYRzRnSUdOdmRXNTBjMTlzYVhOMElEd3RJSEpsY0d4cFkyRjBaU2h1Y205M0tHMWxkR0VwTENCdmJtVmZjMkZ0Y0d4bFgyTnZkVzUwY3lncExDQnphVzF3YkdsbWVTQTlJRVpCVEZORktWeHVJQ0JqYjNWdWRITmZaR1lnUEMwZ1pHOHVZMkZzYkNoeVltbHVaQ3dnWTI5MWJuUnpYMnhwYzNRcFhHNGdJR052Ykc1aGJXVnpLR052ZFc1MGMxOWtaaWtnUEMwZ1kyVnNiRjkwZVhCbGMxeHVYRzRnSUdOdmRXNTBjMTlzYjI1bklEd3RJSEpsYzJoaGNHVW9YRzRnSUNBZ1pHRjBZUzVtY21GdFpTaHpZVzF3YkdWZmFXUWdQU0J0WlhSaEpITmhiWEJzWlY5cFpDd2dZMjkxYm5SelgyUm1MQ0JqYUdWamF5NXVZVzFsY3lBOUlFWkJURk5GS1N4Y2JpQWdJQ0IyWVhKNWFXNW5JRDBnWTJWc2JGOTBlWEJsY3l3Z2RpNXVZVzFsY3lBOUlGd2lZMjkxYm5SY0lpd2dkR2x0WlhaaGNpQTlJRndpWTJWc2JGOTBlWEJsWENJc1hHNGdJQ0FnZEdsdFpYTWdQU0JqWld4c1gzUjVjR1Z6TENCa2FYSmxZM1JwYjI0Z1BTQmNJbXh2Ym1kY0lseHVJQ0FwWEc0Z0lISnZkMjVoYldWektHTnZkVzUwYzE5c2IyNW5LU0E4TFNCT1ZVeE1YRzVjYmlBZ0l5Qk5aWEpuWlNCa2FYTmxZWE5sTDNacGMybDBJSE52SUhkbElHTmhiaUJoY0hCc2VTQmxabVpsWTNSelhHNGdJR052ZFc1MGMxOXNiMjVuSUR3dElHMWxjbWRsS0dOdmRXNTBjMTlzYjI1bkxGeHVJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNCdFpYUmhXeXdnWXloY0luTmhiWEJzWlY5cFpGd2lMQ0JjSW5acGMybDBYQ0lzSUZ3aVpHbHpaV0Z6WlZ3aUtWMHNYRzRnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lHSjVJRDBnWENKellXMXdiR1ZmYVdSY0lpd2djMjl5ZENBOUlFWkJURk5GS1Z4dVhHNGdJQ01nTFMwdElFRndjR3g1SUdWbVptVmpkSE1nYjI1c2VTQjBieUJwYm5SbGNtRmpkR2x1WnlCalpXeHNJSFI1Y0dWeklDMHRMUzFjYmlBZ2FYTmZhVzUwWlhKaFkzUnBibWNnUEMwZ1kyOTFiblJ6WDJ4dmJtY2tZMlZzYkY5MGVYQmxJQ1ZwYmlVZ2FXNTBaWEpoWTNSZlkyVnNiRjkwZVhCbGMxeHVJQ0JsWm1aZmRtVmpJRHd0SUc1MWJXVnlhV01vYm5KdmR5aGpiM1Z1ZEhOZmJHOXVaeWtwWEc0Z0lHbG1JQ2hoYm5rb2FYTmZhVzUwWlhKaFkzUnBibWNwS1NCN1hHNGdJQ0FnSXlCdFlYQWdaR2x5WldOMGFXOXVJSEJsY2lCalpXeHNJSFI1Y0dWY2JpQWdJQ0JrYVhKZmJXRndJRHd0SUhObGRFNWhiV1Z6S0dScGNtVmpkR2x2Ymw5aWVWOWpiSFZ6ZEdWeUxDQnViU0E5SUc1aGJXVnpLR1JwY21WamRHbHZibDlpZVY5amJIVnpkR1Z5S1NsY2JpQWdJQ0JrYVhKZlkzUWdJRHd0SUhWdWJtRnRaU2hrYVhKZmJXRndXMk52ZFc1MGMxOXNiMjVuSkdObGJHeGZkSGx3WlZ0cGMxOXBiblJsY21GamRHbHVaMTFkS1Z4dUlDQWdJQ01nZG1semFYUWdaV1ptWldOMGN5QmllU0JuY205MWNGeHVJQ0FnSUhaZmFXUjRJQ0FnUEMwZ1kyOTFiblJ6WDJ4dmJtY2tkbWx6YVhSYmFYTmZhVzUwWlhKaFkzUnBibWRkSUNzZ01VeGNiaUFnSUNCcGMxOXdjbTluSUR3dElHTnZkVzUwYzE5c2IyNW5KR1JwYzJWaGMyVmJhWE5mYVc1MFpYSmhZM1JwYm1kZElEMDlJREZNWEc0Z0lDQWdkbVVnSUNBZ0lDQThMU0JwWm1Wc2MyVW9hWE5mY0hKdlp5d2dkbWx6YVhSZlpXWm1aV04wYzE5d2NtOW5jbVZ6YzI5eVczWmZhV1I0WFN3Z2RtbHphWFJmWldabVpXTjBjMTlqYjI1MGNtOXNXM1pmYVdSNFhTbGNiaUFnSUNCbFptWmZkbVZqVzJselgybHVkR1Z5WVdOMGFXNW5YU0E4TFNCa2FYSmZZM1FnS2lCMlpWeHVJQ0I5WEc1Y2JpQWdZMjkxYm5SelgyeHZibWNrWVdScVgyTnZkVzUwSUR3dElIQnRZWGdvTUV3c0lISnZkVzVrS0dOdmRXNTBjMTlzYjI1bkpHTnZkVzUwSUNvZ0tERWdLeUJsWm1aZmRtVmpLU2twWEc1Y2JpQWdJeUF0TFMwZ1JYaHdZVzVrSUhSdklIQmxjaTFqWld4c0lISnZkM01nWVc1a0lHRjBkR0ZqYUNCdFpYUmhaR0YwWVNBb1NVNURURlZFU1U1SElHbHVkR1Z5WVdOMFgzUmxjbTBwSUMwdExTMWNiaUFnY21Wd1gyVmhZMmdnUEMwZ1puVnVZM1JwYjI0b2VDd2dkR2x0WlhNcElHbG1JQ2hzWlc1bmRHZ29lQ2tnUFQwZ01Da2dlQ0JsYkhObElISmxjQ2g0TENCMGFXMWxjeUE5SUhScGJXVnpLVnh1SUNCbGVIQmhibVJsWkNBOExTQmtZWFJoTG1aeVlXMWxLRnh1SUNBZ0lITmhiWEJzWlY5cFpDQTlJSEpsY0Y5bFlXTm9LR052ZFc1MGMxOXNiMjVuSkhOaGJYQnNaVjlwWkN3Z1kyOTFiblJ6WDJ4dmJtY2tZV1JxWDJOdmRXNTBLU3hjYmlBZ0lDQmpaV3hzWDNSNWNHVWdQU0J5WlhCZlpXRmphQ2hqYjNWdWRITmZiRzl1WnlSalpXeHNYM1I1Y0dVc0lHTnZkVzUwYzE5c2IyNW5KR0ZrYWw5amIzVnVkQ2tzWEc0Z0lDQWdjM1J5YVc1bmMwRnpSbUZqZEc5eWN5QTlJRVpCVEZORlhHNGdJQ2xjYmx4dUlDQWpJRTFsY21kbElHMWxkR0ZrWVhSaE95QnJaV1Z3SUdsdWRHVnlZV04wYVc5dUlDc2dhVzUwWlhKaFkzUmZkR1Z5YlZ4dUlDQnJaV1Z3WDJOdmJITWdQQzBnWXloY0luTmhiWEJzWlY5cFpGd2lMRndpYzNWaWFtVmpkRjlwWkZ3aUxGd2lkbWx6YVhSY0lpeGNJbk5sZUZ3aUxGd2laR2x6WldGelpWd2lMRndpWVdkbFhDSXNYQ0ppYldsY0lpeGNJbUpoZEdOb1hDSXNYRzRnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJRndpYVc1MFpYSmhZM1JwYjI1Y0lpeGNJbWx1ZEdWeVlXTjBYM1JsY20xY0lpbGNiaUFnWkhWdGJYbGZaR0YwWVNBOExTQnRaWEpuWlNobGVIQmhibVJsWkN3Z2JXVjBZVnNzSUd0bFpYQmZZMjlzYzEwc0lHSjVJRDBnWENKellXMXdiR1ZmYVdSY0lpd2djMjl5ZENBOUlFWkJURk5GS1Z4dVhHNGdJQ01nVTJoMVptWnNaU0J5YjNkeklHWnZjaUJ5WldGc2FYTnRYRzRnSUdsbUlDaHVjbTkzS0dSMWJXMTVYMlJoZEdFcElENGdNU2tnZTF4dUlDQWdJR1IxYlcxNVgyUmhkR0VnUEMwZ1pIVnRiWGxmWkdGMFlWdHpZVzF3YkdVdWFXNTBLRzV5YjNjb1pIVnRiWGxmWkdGMFlTa3BMQ0FzSUdSeWIzQWdQU0JHUVV4VFJWMWNiaUFnSUNCeWIzZHVZVzFsY3loa2RXMXRlVjlrWVhSaEtTQThMU0JPVlV4TVhHNGdJSDFjYmx4dUlDQmtkVzF0ZVY5a1lYUmhYRzU5WEc1Y2JseHVYRzVjYmx4dVlHQmdJbjA9IC0tPlxuXG5gYGByXG4jbmV3IGZ1bmN0aW9uXG5cbiMnIFRpdGxlXG4jJ1xuIycgQHBhcmFtIG5fY2VsbHNcbiMnIEBwYXJhbSBzZF9jZWxsdHlwZXNcbiMnIEBwYXJhbSBuX21ham9yX2NlbGxfdHlwZXNcbiMnIEBwYXJhbSBuX21pbm9yX2NlbGxfdHlwZXNcbiMnIEBwYXJhbSByZWxhdGl2ZV9hYnVuZGFuY2VcbiMnIEBwYXJhbSBuX21ham9yX2ludGVyYWN0X2NlbGx0eXBlc1xuIycgQHBhcmFtIG5fbWlub3JfaW50ZXJhY3RfY2VsbHR5cGVzXG4jJyBAcGFyYW0gbl9pbmRpdmlkdWFsc1xuIycgQHBhcmFtIG5fYmF0Y2hzXG4jJyBAcGFyYW0gaW50ZXJhY3Rpb25fZmVhdHVyZVxuIycgQHBhcmFtIHRpbWVfcG9pbnRzXG4jJyBAcGFyYW0gdGVzdF92YXJcbiMnIEBwYXJhbSBwcm9wX2Rpc2Vhc2VcbiMnIEBwYXJhbSBmY19pbnRlcmFjdFxuIycgQHBhcmFtIGludGVyYWN0aW9uX3R5cGVcbiMnIEBwYXJhbSBzZWVkXG4jJyBAcGFyYW0gdmlzaXRfZWZmZWN0c19wcm9ncmVzc29yXG4jJyBAcGFyYW0gdmlzaXRfZWZmZWN0c19jb250cm9sXG4jJyBAcGFyYW0gZGlyZWN0aW9uX2J5X2NsdXN0ZXJcbiMnXG4jJyBAcmV0dXJuXG4jJyBAZXhwb3J0XG4jJ1xuIycgQGV4YW1wbGVzXG5nZW5lcmF0ZV9kdW1teV9kYXRhX25ldyA8LSBmdW5jdGlvbihcbiAgbl9jZWxscyA9IDMwMDAsICAgICAgICAgICAgICAgICAgIyBiYXNlbGluZSBjZWxscyBwZXIgbWFqb3IgY2VsbCB0eXBlIHBlciBzYW1wbGVcbiAgc2RfY2VsbHR5cGVzID0gMC4xMCwgICAgICAgICAgICAgIyByZWxhdGl2ZSBzZCBmb3IgY291bnRzXG4gIG5fbWFqb3JfY2VsbF90eXBlcyA9IDcsXG4gIG5fbWlub3JfY2VsbF90eXBlcyA9IDMsXG4gIHJlbGF0aXZlX2FidW5kYW5jZSA9IDAuMTAsICAgICAgICMgbWlub3IgdnMgbWFqb3IgYmFzZWxpbmUgcmF0aW9cbiAgbl9tYWpvcl9pbnRlcmFjdF9jZWxsdHlwZXMgPSAxLCAgIyBob3cgbWFueSBtYWpvcnMgYXJlIFxcaW50ZXJhY3RpbmdcXFxuICBuX21pbm9yX2ludGVyYWN0X2NlbGx0eXBlcyA9IDEsICAjIGhvdyBtYW55IG1pbm9ycyBhcmUgXFxpbnRlcmFjdGluZ1xcXG4gIG5faW5kaXZpZHVhbHMgPSAzMCxcbiAgbl9iYXRjaHMgPSA0LFxuXG4gIGludGVyYWN0aW9uX2ZlYXR1cmUgPSBcXHZpc2l0XFwsICAgIyBrZXB0IGZvciBsYWJlbGluZ1xuICB0aW1lX3BvaW50cyA9IDQsICAgICAgICAgICAgICAgICAjID49IDI7IHdvcmtzIGZvciAzK1xuICB0ZXN0X3ZhciA9IFxcZGlzZWFzZVxcLFxuICBwcm9wX2Rpc2Vhc2UgPSAwLjUwLFxuXG4gIGZjX2ludGVyYWN0ID0gMC4xMCwgICAgICAgICAgICAgICMgZWZmZWN0IG1hZ25pdHVkZSB1c2VkIGJ5IGRlZmF1bHRzIGJlbG93XG4gIGludGVyYWN0aW9uX3R5cGUgPSBjKFxcc3BlY2lmaWNcXCxcXGRpZmZlcmVudGlhbFxcLFxcb3Bwb3NpdGVcXCksXG4gIHNlZWQgPSAxMjM0LFxuXG4gIHZpc2l0X2VmZmVjdHNfcHJvZ3Jlc3NvciA9IE5VTEwsICMgbXVsdGlwbGljYXRpdmUgY2hhbmdlOiArMC4xIG1lYW5zICsxMCUgdnMgYmFzZWxpbmVcbiAgdmlzaXRfZWZmZWN0c19jb250cm9sICAgID0gTlVMTCwgIyBkZWZhdWx0IDBzXG4gIGRpcmVjdGlvbl9ieV9jbHVzdGVyICAgICA9IE5VTEwgICMgKzEvLTEgZm9yIGludGVyYWN0aW5nIGNsdXN0ZXJzLCByZWN5Y2xlZCBhcyBuZWVkZWRcbikge1xuICBzZXQuc2VlZChzZWVkKVxuICBpbnRlcmFjdGlvbl90eXBlIDwtIG1hdGNoLmFyZyhpbnRlcmFjdGlvbl90eXBlKVxuXG4gICMgLS0tIEJhc2ljIHNldHVwIC0tLS1cbiAgbl9jZWxsX3R5cGVzIDwtIG5fbWFqb3JfY2VsbF90eXBlcyArIG5fbWlub3JfY2VsbF90eXBlc1xuICBzdG9waWZub3QodGltZV9wb2ludHMgPj0gMiwgbl9jZWxsX3R5cGVzID49IDEpXG5cbiAgY2VsbF90eXBlcyA8LSBMRVRURVJTW3NlcV9sZW4obl9jZWxsX3R5cGVzKV1cbiAgbWFqb3JfaWR4ICA8LSBzZXFfbGVuKG5fbWFqb3JfY2VsbF90eXBlcylcbiAgbWlub3JfaWR4ICA8LSBpZiAobl9taW5vcl9jZWxsX3R5cGVzID4gMCkgKG5fbWFqb3JfY2VsbF90eXBlcyArIHNlcV9sZW4obl9taW5vcl9jZWxsX3R5cGVzKSkgZWxzZSBpbnRlZ2VyKDApXG5cbiAgIyB3aGljaCBjZWxsIHR5cGVzIGFyZSBcXGludGVyYWN0aW5nXFwgKGZpcnN0IHNvbWUgbWFqb3JzLCBsYXN0IHNvbWUgbWlub3JzKVxuICBpbnRlcmFjdF9pZHggPC0gYyhcbiAgICBoZWFkKG1ham9yX2lkeCwgbl9tYWpvcl9pbnRlcmFjdF9jZWxsdHlwZXMpLFxuICAgIHRhaWwoc2VxX2xlbihuX2NlbGxfdHlwZXMpLCBuX21pbm9yX2ludGVyYWN0X2NlbGx0eXBlcylcbiAgKVxuICBpbnRlcmFjdF9pZHggPC0gaW50ZXJzZWN0KGludGVyYWN0X2lkeCwgc2VxX2xlbihuX2NlbGxfdHlwZXMpKSAgIyBndWFyZCByYWlsc1xuICBpbnRlcmFjdF9jZWxsX3R5cGVzIDwtIGNlbGxfdHlwZXNbaW50ZXJhY3RfaWR4XVxuXG4gICMgLS0tIFN1YmplY3RzIGFuZCB2aXNpdHMgLS0tLVxuICBzdWJqZWN0X2lkIDwtIHBhc3RlMChcXFNVQl9cXCwgc2VxX2xlbihuX2luZGl2aWR1YWxzKSlcbiAgZGlzZWFzZV92ZWMgPC0gYyhyZXAoMUwsIHJvdW5kKG5faW5kaXZpZHVhbHMgKiBwcm9wX2Rpc2Vhc2UpKSxcbiAgICAgICAgICAgICAgICAgICByZXAoMEwsIG5faW5kaXZpZHVhbHMgLSByb3VuZChuX2luZGl2aWR1YWxzICogcHJvcF9kaXNlYXNlKSkpXG4gIGRpc2Vhc2VfdmVjIDwtIHNhbXBsZShkaXNlYXNlX3ZlYywgbl9pbmRpdmlkdWFscylcblxuICBzZXhfdmVjIDwtIHNhbXBsZShjKDBMLCAxTCksIG5faW5kaXZpZHVhbHMsIHJlcGxhY2UgPSBUUlVFKSAgICAgICAgICAjIDAvMVxuICBhZ2VfdmVjIDwtIHNhbXBsZSgxODo2MCwgbl9pbmRpdmlkdWFscywgcmVwbGFjZSA9IFRSVUUpXG4gIGJtaV92ZWMgPC0gc2FtcGxlKDE1OjM1LCBuX2luZGl2aWR1YWxzLCByZXBsYWNlID0gVFJVRSlcbiAgYmF0Y2hfdmVjIDwtIHJlcChzZXFfbGVuKG5fYmF0Y2hzKSwgbGVuZ3RoLm91dCA9IG5faW5kaXZpZHVhbHMpXG5cbiAgc3ViamVjdHMgPC0gZGF0YS5mcmFtZShcbiAgICBzdWJqZWN0X2lkID0gc3ViamVjdF9pZCxcbiAgICBzZXggICA9IHNleF92ZWMsXG4gICAgZGlzZWFzZSA9IGRpc2Vhc2VfdmVjLCAgICAgICAgICAgIyBjYW5vbmljYWwgZXhwb3N1cmUgc3RvcmFnZVxuICAgIGFnZSAgID0gYWdlX3ZlYyxcbiAgICBibWkgICA9IGJtaV92ZWMsXG4gICAgYmF0Y2ggPSBmYWN0b3IoYmF0Y2hfdmVjKSxcbiAgICBzdHJpbmdzQXNGYWN0b3JzID0gRkFMU0VcbiAgKVxuXG4gIHZpc2l0cyA8LSBkYXRhLmZyYW1lKFxuICAgIHN1YmplY3RfaWQgPSByZXAoc3ViamVjdF9pZCwgZWFjaCA9IHRpbWVfcG9pbnRzKSxcbiAgICB2aXNpdCAgICAgID0gcmVwKDA6KHRpbWVfcG9pbnRzIC0gMSksIHRpbWVzID0gbl9pbmRpdmlkdWFscyksXG4gICAgc3RyaW5nc0FzRmFjdG9ycyA9IEZBTFNFXG4gIClcbiAgdmlzaXRzJHNhbXBsZV9pZCA8LSBwYXN0ZTAodmlzaXRzJHN1YmplY3RfaWQsIFxcX1ZcXCwgdmlzaXRzJHZpc2l0KVxuXG4gIG1ldGEgPC0gbWVyZ2UodmlzaXRzLCBzdWJqZWN0cywgYnkgPSBcXHN1YmplY3RfaWRcXCwgc29ydCA9IEZBTFNFKVxuXG4gICMgTWFrZSBzdXJlIGEgY29sdW1uIG5hbWVkIGB0ZXN0X3ZhcmAgZXhpc3RzIChldmVuIGlmIHRlc3RfdmFyICE9IFxcZGlzZWFzZVxcKVxuICBpZiAoIWlkZW50aWNhbCh0ZXN0X3ZhciwgXFxkaXNlYXNlXFwpKSB7XG4gICAgbWV0YVtbdGVzdF92YXJdXSA8LSBtZXRhW1tcXGRpc2Vhc2VcXF1dXG4gIH1cblxuICAjIExhYmVsIGFuZCBJTlRFUkFDVElPTiBURVJNIChwZXJzaXN0IHRvIG91dHB1dClcbiAgbWV0YSRpbnRlcmFjdGlvbiAgICA8LSBwYXN0ZTAoaW50ZXJhY3Rpb25fZmVhdHVyZSwgXFw6XFwsIHRlc3RfdmFyKVxuICBtZXRhJGludGVyYWN0X3Rlcm0gIDwtIGFzLmludGVnZXIobWV0YVtbaW50ZXJhY3Rpb25fZmVhdHVyZV1dKSAqIGFzLmludGVnZXIobWV0YVtbdGVzdF92YXJdXSlcblxuICAjIC0tLSBEZWZhdWx0IHZpc2l0IGVmZmVjdHMgKGxlbmd0aCA9IHRpbWVfcG9pbnRzKSAtLS0tXG4gIGlmIChpcy5udWxsKHZpc2l0X2VmZmVjdHNfcHJvZ3Jlc3NvcikpIHtcbiAgICBpZiAoaW50ZXJhY3Rpb25fdHlwZSA9PSBcXHNwZWNpZmljXFwpIHtcbiAgICAgIHZlIDwtIHJlcCgwLCB0aW1lX3BvaW50cylcbiAgICAgIGlmICh0aW1lX3BvaW50cyA+PSAyKSB2ZVsyXSA8LSBmY19pbnRlcmFjdCAgIyBidW1wIFYxIG9ubHlcbiAgICAgIHZpc2l0X2VmZmVjdHNfcHJvZ3Jlc3NvciA8LSB2ZVxuICAgIH0gZWxzZSBpZiAoaW50ZXJhY3Rpb25fdHlwZSA9PSBcXGRpZmZlcmVudGlhbFxcKSB7XG4gICAgICB2aXNpdF9lZmZlY3RzX3Byb2dyZXNzb3IgPC0gcmVwKGZjX2ludGVyYWN0LCB0aW1lX3BvaW50cylcbiAgICB9IGVsc2UgeyAjIFxcb3Bwb3NpdGVcXDogYWx0ZXJuYXRlICsvLSBzdGFydGluZyBhdCBWMFxuICAgICAgdmUgPC0gcmVwKDAsIHRpbWVfcG9pbnRzKVxuICAgICAgdmVbc2VxKDEsIHRpbWVfcG9pbnRzLCBieSA9IDIpXSA8LSArZmNfaW50ZXJhY3QgICMgVjAsIFYyLCAuLi5cbiAgICAgIGlmICh0aW1lX3BvaW50cyA+PSAyKSB2ZVtzZXEoMiwgdGltZV9wb2ludHMsIGJ5ID0gMildIDwtIC1mY19pbnRlcmFjdCAjIFYxLCBWMywgLi4uXG4gICAgICB2aXNpdF9lZmZlY3RzX3Byb2dyZXNzb3IgPC0gdmVcbiAgICB9XG4gIH1cbiAgaWYgKGlzLm51bGwodmlzaXRfZWZmZWN0c19jb250cm9sKSkge1xuICAgIHZpc2l0X2VmZmVjdHNfY29udHJvbCA8LSByZXAoMCwgdGltZV9wb2ludHMpXG4gIH1cbiAgc3RvcGlmbm90KGxlbmd0aCh2aXNpdF9lZmZlY3RzX3Byb2dyZXNzb3IpID09IHRpbWVfcG9pbnRzLFxuICAgICAgICAgICAgbGVuZ3RoKHZpc2l0X2VmZmVjdHNfY29udHJvbCkgICAgPT0gdGltZV9wb2ludHMpXG5cbiAgIyBEaXJlY3Rpb24gcGVyIGludGVyYWN0aW5nIGNsdXN0ZXJcbiAgaWYgKGlzLm51bGwoZGlyZWN0aW9uX2J5X2NsdXN0ZXIpKSB7XG4gICAgIyBvbGQgZGVmYXVsdDogYWxsIGludGVyYWN0aW5nIGNlbGwgdHlwZXMgZ28gaW4gdGhlIHNhbWUgZGlyZWN0aW9uICgrMSlcbiAgICBkaXJlY3Rpb25fYnlfY2x1c3RlciA8LSByZXAoMUwsIG1heCgxLCBsZW5ndGgoaW50ZXJhY3RfY2VsbF90eXBlcykpKVxuICAgIG5hbWVzKGRpcmVjdGlvbl9ieV9jbHVzdGVyKSA8LSBpbnRlcmFjdF9jZWxsX3R5cGVzXG4gIH0gZWxzZSB7XG4gICAgaWYgKCFpcy5udWxsKG5hbWVzKGRpcmVjdGlvbl9ieV9jbHVzdGVyKSkpIHtcbiAgICAgICMgTkVXOiB1c2VyIHBhc3NlZCBhICpuYW1lZCogdmVjdG9yLCBlLmcuIGMoQSA9IDEsIEMgPSAtMSwgSiA9IDEpXG4gICAgICAjIFdlIGFsaWduIHRvIGludGVyYWN0X2NlbGxfdHlwZXMgYW5kIGRlZmF1bHQgbWlzc2luZyBvbmVzIHRvICsxLlxuICAgICAgdG1wIDwtIHJlcCgxTCwgbWF4KDEsIGxlbmd0aChpbnRlcmFjdF9jZWxsX3R5cGVzKSkpXG4gICAgICBuYW1lcyh0bXApIDwtIGludGVyYWN0X2NlbGxfdHlwZXNcblxuICAgICAgbWF0Y2hlZCA8LSBpbnRlcnNlY3QobmFtZXMoZGlyZWN0aW9uX2J5X2NsdXN0ZXIpLCBpbnRlcmFjdF9jZWxsX3R5cGVzKVxuICAgICAgdG1wW21hdGNoZWRdIDwtIGRpcmVjdGlvbl9ieV9jbHVzdGVyW21hdGNoZWRdXG5cbiAgICAgIGRpcmVjdGlvbl9ieV9jbHVzdGVyIDwtIHRtcFxuICAgIH0gZWxzZSB7XG4gICAgICAjIG9sZCBiZWhhdmlvciBmb3IgdW5uYW1lZCB2ZWN0b3I6IHJlY3ljbGUgb3ZlciBpbnRlcmFjdGluZyBjZWxsIHR5cGVzXG4gICAgICBkaXJlY3Rpb25fYnlfY2x1c3RlciA8LSByZXAoZGlyZWN0aW9uX2J5X2NsdXN0ZXIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgbGVuZ3RoLm91dCA9IGxlbmd0aChpbnRlcmFjdF9jZWxsX3R5cGVzKSlcbiAgICAgIG5hbWVzKGRpcmVjdGlvbl9ieV9jbHVzdGVyKSA8LSBpbnRlcmFjdF9jZWxsX3R5cGVzXG4gICAgfVxuICB9XG5cbiAgIyAtLS0gQmFzZWxpbmUgY291bnRzIHBlciAoc2FtcGxlLCBjZWxsX3R5cGUpIC0tLS1cbiAgb25lX3NhbXBsZV9jb3VudHMgPC0gZnVuY3Rpb24oKSB7XG4gICAgbWFqb3JfY291bnRzIDwtIHJvdW5kKHJ1bmlmKG5fbWFqb3JfY2VsbF90eXBlcyxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgbWluID0gbl9jZWxscyAqICgxIC0gc2RfY2VsbHR5cGVzKSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgbWF4ID0gbl9jZWxscyAqICgxICsgc2RfY2VsbHR5cGVzKSkpXG4gICAgbWlub3JfY291bnRzIDwtIGlmIChuX21pbm9yX2NlbGxfdHlwZXMgPiAwKSB7XG4gICAgICByb3VuZChydW5pZihuX21pbm9yX2NlbGxfdHlwZXMsXG4gICAgICAgICAgICAgICAgICBtaW4gPSBuX2NlbGxzICogcmVsYXRpdmVfYWJ1bmRhbmNlICogKDEgLSBzZF9jZWxsdHlwZXMpLFxuICAgICAgICAgICAgICAgICAgbWF4ID0gbl9jZWxscyAqIHJlbGF0aXZlX2FidW5kYW5jZSAqICgxICsgc2RfY2VsbHR5cGVzKSkpXG4gICAgfSBlbHNlIGludGVnZXIoMClcbiAgICBjKG1ham9yX2NvdW50cywgbWlub3JfY291bnRzKVxuICB9XG5cbiAgY291bnRzX2xpc3QgPC0gcmVwbGljYXRlKG5yb3cobWV0YSksIG9uZV9zYW1wbGVfY291bnRzKCksIHNpbXBsaWZ5ID0gRkFMU0UpXG4gIGNvdW50c19kZiA8LSBkby5jYWxsKHJiaW5kLCBjb3VudHNfbGlzdClcbiAgY29sbmFtZXMoY291bnRzX2RmKSA8LSBjZWxsX3R5cGVzXG5cbiAgY291bnRzX2xvbmcgPC0gcmVzaGFwZShcbiAgICBkYXRhLmZyYW1lKHNhbXBsZV9pZCA9IG1ldGEkc2FtcGxlX2lkLCBjb3VudHNfZGYsIGNoZWNrLm5hbWVzID0gRkFMU0UpLFxuICAgIHZhcnlpbmcgPSBjZWxsX3R5cGVzLCB2Lm5hbWVzID0gXFxjb3VudFxcLCB0aW1ldmFyID0gXFxjZWxsX3R5cGVcXCxcbiAgICB0aW1lcyA9IGNlbGxfdHlwZXMsIGRpcmVjdGlvbiA9IFxcbG9uZ1xcXG4gIClcbiAgcm93bmFtZXMoY291bnRzX2xvbmcpIDwtIE5VTExcblxuICAjIE1lcmdlIGRpc2Vhc2UvdmlzaXQgc28gd2UgY2FuIGFwcGx5IGVmZmVjdHNcbiAgY291bnRzX2xvbmcgPC0gbWVyZ2UoY291bnRzX2xvbmcsXG4gICAgICAgICAgICAgICAgICAgICAgIG1ldGFbLCBjKFxcc2FtcGxlX2lkXFwsIFxcdmlzaXRcXCwgXFxkaXNlYXNlXFwpXSxcbiAgICAgICAgICAgICAgICAgICAgICAgYnkgPSBcXHNhbXBsZV9pZFxcLCBzb3J0ID0gRkFMU0UpXG5cbiAgIyAtLS0gQXBwbHkgZWZmZWN0cyBvbmx5IHRvIGludGVyYWN0aW5nIGNlbGwgdHlwZXMgLS0tLVxuICBpc19pbnRlcmFjdGluZyA8LSBjb3VudHNfbG9uZyRjZWxsX3R5cGUgJWluJSBpbnRlcmFjdF9jZWxsX3R5cGVzXG4gIGVmZl92ZWMgPC0gbnVtZXJpYyhucm93KGNvdW50c19sb25nKSlcbiAgaWYgKGFueShpc19pbnRlcmFjdGluZykpIHtcbiAgICAjIG1hcCBkaXJlY3Rpb24gcGVyIGNlbGwgdHlwZVxuICAgIGRpcl9tYXAgPC0gc2V0TmFtZXMoZGlyZWN0aW9uX2J5X2NsdXN0ZXIsIG5tID0gbmFtZXMoZGlyZWN0aW9uX2J5X2NsdXN0ZXIpKVxuICAgIGRpcl9jdCAgPC0gdW5uYW1lKGRpcl9tYXBbY291bnRzX2xvbmckY2VsbF90eXBlW2lzX2ludGVyYWN0aW5nXV0pXG4gICAgIyB2aXNpdCBlZmZlY3RzIGJ5IGdyb3VwXG4gICAgdl9pZHggICA8LSBjb3VudHNfbG9uZyR2aXNpdFtpc19pbnRlcmFjdGluZ10gKyAxTFxuICAgIGlzX3Byb2cgPC0gY291bnRzX2xvbmckZGlzZWFzZVtpc19pbnRlcmFjdGluZ10gPT0gMUxcbiAgICB2ZSAgICAgIDwtIGlmZWxzZShpc19wcm9nLCB2aXNpdF9lZmZlY3RzX3Byb2dyZXNzb3Jbdl9pZHhdLCB2aXNpdF9lZmZlY3RzX2NvbnRyb2xbdl9pZHhdKVxuICAgIGVmZl92ZWNbaXNfaW50ZXJhY3RpbmddIDwtIGRpcl9jdCAqIHZlXG4gIH1cblxuICBjb3VudHNfbG9uZyRhZGpfY291bnQgPC0gcG1heCgwTCwgcm91bmQoY291bnRzX2xvbmckY291bnQgKiAoMSArIGVmZl92ZWMpKSlcblxuICAjIC0tLSBFeHBhbmQgdG8gcGVyLWNlbGwgcm93cyBhbmQgYXR0YWNoIG1ldGFkYXRhIChJTkNMVURJTkcgaW50ZXJhY3RfdGVybSkgLS0tLVxuICByZXBfZWFjaCA8LSBmdW5jdGlvbih4LCB0aW1lcykgaWYgKGxlbmd0aCh4KSA9PSAwKSB4IGVsc2UgcmVwKHgsIHRpbWVzID0gdGltZXMpXG4gIGV4cGFuZGVkIDwtIGRhdGEuZnJhbWUoXG4gICAgc2FtcGxlX2lkID0gcmVwX2VhY2goY291bnRzX2xvbmckc2FtcGxlX2lkLCBjb3VudHNfbG9uZyRhZGpfY291bnQpLFxuICAgIGNlbGxfdHlwZSA9IHJlcF9lYWNoKGNvdW50c19sb25nJGNlbGxfdHlwZSwgY291bnRzX2xvbmckYWRqX2NvdW50KSxcbiAgICBzdHJpbmdzQXNGYWN0b3JzID0gRkFMU0VcbiAgKVxuXG4gICMgTWVyZ2UgbWV0YWRhdGE7IGtlZXAgaW50ZXJhY3Rpb24gKyBpbnRlcmFjdF90ZXJtXG4gIGtlZXBfY29scyA8LSBjKFxcc2FtcGxlX2lkXFwsXFxzdWJqZWN0X2lkXFwsXFx2aXNpdFxcLFxcc2V4XFwsXFxkaXNlYXNlXFwsXFxhZ2VcXCxcXGJtaVxcLFxcYmF0Y2hcXCxcbiAgICAgICAgICAgICAgICAgXFxpbnRlcmFjdGlvblxcLFxcaW50ZXJhY3RfdGVybVxcKVxuICBkdW1teV9kYXRhIDwtIG1lcmdlKGV4cGFuZGVkLCBtZXRhWywga2VlcF9jb2xzXSwgYnkgPSBcXHNhbXBsZV9pZFxcLCBzb3J0ID0gRkFMU0UpXG5cbiAgIyBTaHVmZmxlIHJvd3MgZm9yIHJlYWxpc21cbiAgaWYgKG5yb3coZHVtbXlfZGF0YSkgPiAxKSB7XG4gICAgZHVtbXlfZGF0YSA8LSBkdW1teV9kYXRhW3NhbXBsZS5pbnQobnJvdyhkdW1teV9kYXRhKSksICwgZHJvcCA9IEZBTFNFXVxuICAgIHJvd25hbWVzKGR1bW15X2RhdGEpIDwtIE5VTExcbiAgfVxuXG4gIGR1bW15X2RhdGFcbn1cblxuXG5cblxuXG5gYGBcblxuPCEtLSBybmItc291cmNlLWVuZCAtLT5cbiJ9 -->

````

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuI25ldyBmdW5jdGlvblxuXG4jJyBUaXRsZVxuIydcbiMnIEBwYXJhbSBuX2NlbGxzXG4jJyBAcGFyYW0gc2RfY2VsbHR5cGVzXG4jJyBAcGFyYW0gbl9tYWpvcl9jZWxsX3R5cGVzXG4jJyBAcGFyYW0gbl9taW5vcl9jZWxsX3R5cGVzXG4jJyBAcGFyYW0gcmVsYXRpdmVfYWJ1bmRhbmNlXG4jJyBAcGFyYW0gbl9tYWpvcl9pbnRlcmFjdF9jZWxsdHlwZXNcbiMnIEBwYXJhbSBuX21pbm9yX2ludGVyYWN0X2NlbGx0eXBlc1xuIycgQHBhcmFtIG5faW5kaXZpZHVhbHNcbiMnIEBwYXJhbSBuX2JhdGNoc1xuIycgQHBhcmFtIGludGVyYWN0aW9uX2ZlYXR1cmVcbiMnIEBwYXJhbSB0aW1lX3BvaW50c1xuIycgQHBhcmFtIHRlc3RfdmFyXG4jJyBAcGFyYW0gcHJvcF9kaXNlYXNlXG4jJyBAcGFyYW0gZmNfaW50ZXJhY3RcbiMnIEBwYXJhbSBpbnRlcmFjdGlvbl90eXBlXG4jJyBAcGFyYW0gc2VlZFxuIycgQHBhcmFtIHZpc2l0X2VmZmVjdHNfcHJvZ3Jlc3NvclxuIycgQHBhcmFtIHZpc2l0X2VmZmVjdHNfY29udHJvbFxuIycgQHBhcmFtIGRpcmVjdGlvbl9ieV9jbHVzdGVyXG4jJ1xuIycgQHJldHVyblxuIycgQGV4cG9ydFxuIydcbiMnIEBleGFtcGxlc1xuZ2VuZXJhdGVfZHVtbXlfZGF0YV9uZXcgPC0gZnVuY3Rpb24oXG4gIG5fY2VsbHMgPSAzMDAwLCAgICAgICAgICAgICAgICAgICMgYmFzZWxpbmUgY2VsbHMgcGVyIG1ham9yIGNlbGwgdHlwZSBwZXIgc2FtcGxlXG4gIHNkX2NlbGx0eXBlcyA9IDAuMTAsICAgICAgICAgICAgICMgcmVsYXRpdmUgc2QgZm9yIGNvdW50c1xuICBuX21ham9yX2NlbGxfdHlwZXMgPSA3LFxuICBuX21pbm9yX2NlbGxfdHlwZXMgPSAzLFxuICByZWxhdGl2ZV9hYnVuZGFuY2UgPSAwLjEwLCAgICAgICAjIG1pbm9yIHZzIG1ham9yIGJhc2VsaW5lIHJhdGlvXG4gIG5fbWFqb3JfaW50ZXJhY3RfY2VsbHR5cGVzID0gMSwgICMgaG93IG1hbnkgbWFqb3JzIGFyZSBcImludGVyYWN0aW5nXCJcbiAgbl9taW5vcl9pbnRlcmFjdF9jZWxsdHlwZXMgPSAxLCAgIyBob3cgbWFueSBtaW5vcnMgYXJlIFwiaW50ZXJhY3RpbmdcIlxuICBuX2luZGl2aWR1YWxzID0gMzAsXG4gIG5fYmF0Y2hzID0gNCxcblxuICBpbnRlcmFjdGlvbl9mZWF0dXJlID0gXCJ2aXNpdFwiLCAgICMga2VwdCBmb3IgbGFiZWxpbmdcbiAgdGltZV9wb2ludHMgPSA0LCAgICAgICAgICAgICAgICAgIyA+PSAyOyB3b3JrcyBmb3IgMytcbiAgdGVzdF92YXIgPSBcImRpc2Vhc2VcIixcbiAgcHJvcF9kaXNlYXNlID0gMC41MCxcblxuICBmY19pbnRlcmFjdCA9IDAuMTAsICAgICAgICAgICAgICAjIGVmZmVjdCBtYWduaXR1ZGUgdXNlZCBieSBkZWZhdWx0cyBiZWxvd1xuICBpbnRlcmFjdGlvbl90eXBlID0gYyhcInNwZWNpZmljXCIsXCJkaWZmZXJlbnRpYWxcIixcIm9wcG9zaXRlXCIpLFxuICBzZWVkID0gMTIzNCxcblxuICB2aXNpdF9lZmZlY3RzX3Byb2dyZXNzb3IgPSBOVUxMLCAjIG11bHRpcGxpY2F0aXZlIGNoYW5nZTogKzAuMSBtZWFucyArMTAlIHZzIGJhc2VsaW5lXG4gIHZpc2l0X2VmZmVjdHNfY29udHJvbCAgICA9IE5VTEwsICMgZGVmYXVsdCAwc1xuICBkaXJlY3Rpb25fYnlfY2x1c3RlciAgICAgPSBOVUxMICAjICsxLy0xIGZvciBpbnRlcmFjdGluZyBjbHVzdGVycywgcmVjeWNsZWQgYXMgbmVlZGVkXG4pIHtcbiAgc2V0LnNlZWQoc2VlZClcbiAgaW50ZXJhY3Rpb25fdHlwZSA8LSBtYXRjaC5hcmcoaW50ZXJhY3Rpb25fdHlwZSlcblxuICAjIC0tLSBCYXNpYyBzZXR1cCAtLS0tXG4gIG5fY2VsbF90eXBlcyA8LSBuX21ham9yX2NlbGxfdHlwZXMgKyBuX21pbm9yX2NlbGxfdHlwZXNcbiAgc3RvcGlmbm90KHRpbWVfcG9pbnRzID49IDIsIG5fY2VsbF90eXBlcyA+PSAxKVxuXG4gIGNlbGxfdHlwZXMgPC0gTEVUVEVSU1tzZXFfbGVuKG5fY2VsbF90eXBlcyldXG4gIG1ham9yX2lkeCAgPC0gc2VxX2xlbihuX21ham9yX2NlbGxfdHlwZXMpXG4gIG1pbm9yX2lkeCAgPC0gaWYgKG5fbWlub3JfY2VsbF90eXBlcyA+IDApIChuX21ham9yX2NlbGxfdHlwZXMgKyBzZXFfbGVuKG5fbWlub3JfY2VsbF90eXBlcykpIGVsc2UgaW50ZWdlcigwKVxuXG4gICMgd2hpY2ggY2VsbCB0eXBlcyBhcmUgXCJpbnRlcmFjdGluZ1wiIChmaXJzdCBzb21lIG1ham9ycywgbGFzdCBzb21lIG1pbm9ycylcbiAgaW50ZXJhY3RfaWR4IDwtIGMoXG4gICAgaGVhZChtYWpvcl9pZHgsIG5fbWFqb3JfaW50ZXJhY3RfY2VsbHR5cGVzKSxcbiAgICB0YWlsKHNlcV9sZW4obl9jZWxsX3R5cGVzKSwgbl9taW5vcl9pbnRlcmFjdF9jZWxsdHlwZXMpXG4gIClcbiAgaW50ZXJhY3RfaWR4IDwtIGludGVyc2VjdChpbnRlcmFjdF9pZHgsIHNlcV9sZW4obl9jZWxsX3R5cGVzKSkgICMgZ3VhcmQgcmFpbHNcbiAgaW50ZXJhY3RfY2VsbF90eXBlcyA8LSBjZWxsX3R5cGVzW2ludGVyYWN0X2lkeF1cblxuICAjIC0tLSBTdWJqZWN0cyBhbmQgdmlzaXRzIC0tLS1cbiAgc3ViamVjdF9pZCA8LSBwYXN0ZTAoXCJTVUJfXCIsIHNlcV9sZW4obl9pbmRpdmlkdWFscykpXG4gIGRpc2Vhc2VfdmVjIDwtIGMocmVwKDFMLCByb3VuZChuX2luZGl2aWR1YWxzICogcHJvcF9kaXNlYXNlKSksXG4gICAgICAgICAgICAgICAgICAgcmVwKDBMLCBuX2luZGl2aWR1YWxzIC0gcm91bmQobl9pbmRpdmlkdWFscyAqIHByb3BfZGlzZWFzZSkpKVxuICBkaXNlYXNlX3ZlYyA8LSBzYW1wbGUoZGlzZWFzZV92ZWMsIG5faW5kaXZpZHVhbHMpXG5cbiAgc2V4X3ZlYyA8LSBzYW1wbGUoYygwTCwgMUwpLCBuX2luZGl2aWR1YWxzLCByZXBsYWNlID0gVFJVRSkgICAgICAgICAgIyAwLzFcbiAgYWdlX3ZlYyA8LSBzYW1wbGUoMTg6NjAsIG5faW5kaXZpZHVhbHMsIHJlcGxhY2UgPSBUUlVFKVxuICBibWlfdmVjIDwtIHNhbXBsZSgxNTozNSwgbl9pbmRpdmlkdWFscywgcmVwbGFjZSA9IFRSVUUpXG4gIGJhdGNoX3ZlYyA8LSByZXAoc2VxX2xlbihuX2JhdGNocyksIGxlbmd0aC5vdXQgPSBuX2luZGl2aWR1YWxzKVxuXG4gIHN1YmplY3RzIDwtIGRhdGEuZnJhbWUoXG4gICAgc3ViamVjdF9pZCA9IHN1YmplY3RfaWQsXG4gICAgc2V4ICAgPSBzZXhfdmVjLFxuICAgIGRpc2Vhc2UgPSBkaXNlYXNlX3ZlYywgICAgICAgICAgICMgY2Fub25pY2FsIGV4cG9zdXJlIHN0b3JhZ2VcbiAgICBhZ2UgICA9IGFnZV92ZWMsXG4gICAgYm1pICAgPSBibWlfdmVjLFxuICAgIGJhdGNoID0gZmFjdG9yKGJhdGNoX3ZlYyksXG4gICAgc3RyaW5nc0FzRmFjdG9ycyA9IEZBTFNFXG4gIClcblxuICB2aXNpdHMgPC0gZGF0YS5mcmFtZShcbiAgICBzdWJqZWN0X2lkID0gcmVwKHN1YmplY3RfaWQsIGVhY2ggPSB0aW1lX3BvaW50cyksXG4gICAgdmlzaXQgICAgICA9IHJlcCgwOih0aW1lX3BvaW50cyAtIDEpLCB0aW1lcyA9IG5faW5kaXZpZHVhbHMpLFxuICAgIHN0cmluZ3NBc0ZhY3RvcnMgPSBGQUxTRVxuICApXG4gIHZpc2l0cyRzYW1wbGVfaWQgPC0gcGFzdGUwKHZpc2l0cyRzdWJqZWN0X2lkLCBcIl9WXCIsIHZpc2l0cyR2aXNpdClcblxuICBtZXRhIDwtIG1lcmdlKHZpc2l0cywgc3ViamVjdHMsIGJ5ID0gXCJzdWJqZWN0X2lkXCIsIHNvcnQgPSBGQUxTRSlcblxuICAjIE1ha2Ugc3VyZSBhIGNvbHVtbiBuYW1lZCBgdGVzdF92YXJgIGV4aXN0cyAoZXZlbiBpZiB0ZXN0X3ZhciAhPSBcImRpc2Vhc2VcIilcbiAgaWYgKCFpZGVudGljYWwodGVzdF92YXIsIFwiZGlzZWFzZVwiKSkge1xuICAgIG1ldGFbW3Rlc3RfdmFyXV0gPC0gbWV0YVtbXCJkaXNlYXNlXCJdXVxuICB9XG5cbiAgIyBMYWJlbCBhbmQgSU5URVJBQ1RJT04gVEVSTSAocGVyc2lzdCB0byBvdXRwdXQpXG4gIG1ldGEkaW50ZXJhY3Rpb24gICAgPC0gcGFzdGUwKGludGVyYWN0aW9uX2ZlYXR1cmUsIFwiOlwiLCB0ZXN0X3ZhcilcbiAgbWV0YSRpbnRlcmFjdF90ZXJtICA8LSBhcy5pbnRlZ2VyKG1ldGFbW2ludGVyYWN0aW9uX2ZlYXR1cmVdXSkgKiBhcy5pbnRlZ2VyKG1ldGFbW3Rlc3RfdmFyXV0pXG5cbiAgIyAtLS0gRGVmYXVsdCB2aXNpdCBlZmZlY3RzIChsZW5ndGggPSB0aW1lX3BvaW50cykgLS0tLVxuICBpZiAoaXMubnVsbCh2aXNpdF9lZmZlY3RzX3Byb2dyZXNzb3IpKSB7XG4gICAgaWYgKGludGVyYWN0aW9uX3R5cGUgPT0gXCJzcGVjaWZpY1wiKSB7XG4gICAgICB2ZSA8LSByZXAoMCwgdGltZV9wb2ludHMpXG4gICAgICBpZiAodGltZV9wb2ludHMgPj0gMikgdmVbMl0gPC0gZmNfaW50ZXJhY3QgICMgYnVtcCBWMSBvbmx5XG4gICAgICB2aXNpdF9lZmZlY3RzX3Byb2dyZXNzb3IgPC0gdmVcbiAgICB9IGVsc2UgaWYgKGludGVyYWN0aW9uX3R5cGUgPT0gXCJkaWZmZXJlbnRpYWxcIikge1xuICAgICAgdmlzaXRfZWZmZWN0c19wcm9ncmVzc29yIDwtIHJlcChmY19pbnRlcmFjdCwgdGltZV9wb2ludHMpXG4gICAgfSBlbHNlIHsgIyBcIm9wcG9zaXRlXCI6IGFsdGVybmF0ZSArLy0gc3RhcnRpbmcgYXQgVjBcbiAgICAgIHZlIDwtIHJlcCgwLCB0aW1lX3BvaW50cylcbiAgICAgIHZlW3NlcSgxLCB0aW1lX3BvaW50cywgYnkgPSAyKV0gPC0gK2ZjX2ludGVyYWN0ICAjIFYwLCBWMiwgLi4uXG4gICAgICBpZiAodGltZV9wb2ludHMgPj0gMikgdmVbc2VxKDIsIHRpbWVfcG9pbnRzLCBieSA9IDIpXSA8LSAtZmNfaW50ZXJhY3QgIyBWMSwgVjMsIC4uLlxuICAgICAgdmlzaXRfZWZmZWN0c19wcm9ncmVzc29yIDwtIHZlXG4gICAgfVxuICB9XG4gIGlmIChpcy5udWxsKHZpc2l0X2VmZmVjdHNfY29udHJvbCkpIHtcbiAgICB2aXNpdF9lZmZlY3RzX2NvbnRyb2wgPC0gcmVwKDAsIHRpbWVfcG9pbnRzKVxuICB9XG4gIHN0b3BpZm5vdChsZW5ndGgodmlzaXRfZWZmZWN0c19wcm9ncmVzc29yKSA9PSB0aW1lX3BvaW50cyxcbiAgICAgICAgICAgIGxlbmd0aCh2aXNpdF9lZmZlY3RzX2NvbnRyb2wpICAgID09IHRpbWVfcG9pbnRzKVxuXG4gICMgRGlyZWN0aW9uIHBlciBpbnRlcmFjdGluZyBjbHVzdGVyXG4gIGlmIChpcy5udWxsKGRpcmVjdGlvbl9ieV9jbHVzdGVyKSkge1xuICAgICMgb2xkIGRlZmF1bHQ6IGFsbCBpbnRlcmFjdGluZyBjZWxsIHR5cGVzIGdvIGluIHRoZSBzYW1lIGRpcmVjdGlvbiAoKzEpXG4gICAgZGlyZWN0aW9uX2J5X2NsdXN0ZXIgPC0gcmVwKDFMLCBtYXgoMSwgbGVuZ3RoKGludGVyYWN0X2NlbGxfdHlwZXMpKSlcbiAgICBuYW1lcyhkaXJlY3Rpb25fYnlfY2x1c3RlcikgPC0gaW50ZXJhY3RfY2VsbF90eXBlc1xuICB9IGVsc2Uge1xuICAgIGlmICghaXMubnVsbChuYW1lcyhkaXJlY3Rpb25fYnlfY2x1c3RlcikpKSB7XG4gICAgICAjIE5FVzogdXNlciBwYXNzZWQgYSAqbmFtZWQqIHZlY3RvciwgZS5nLiBjKEEgPSAxLCBDID0gLTEsIEogPSAxKVxuICAgICAgIyBXZSBhbGlnbiB0byBpbnRlcmFjdF9jZWxsX3R5cGVzIGFuZCBkZWZhdWx0IG1pc3Npbmcgb25lcyB0byArMS5cbiAgICAgIHRtcCA8LSByZXAoMUwsIG1heCgxLCBsZW5ndGgoaW50ZXJhY3RfY2VsbF90eXBlcykpKVxuICAgICAgbmFtZXModG1wKSA8LSBpbnRlcmFjdF9jZWxsX3R5cGVzXG5cbiAgICAgIG1hdGNoZWQgPC0gaW50ZXJzZWN0KG5hbWVzKGRpcmVjdGlvbl9ieV9jbHVzdGVyKSwgaW50ZXJhY3RfY2VsbF90eXBlcylcbiAgICAgIHRtcFttYXRjaGVkXSA8LSBkaXJlY3Rpb25fYnlfY2x1c3RlclttYXRjaGVkXVxuXG4gICAgICBkaXJlY3Rpb25fYnlfY2x1c3RlciA8LSB0bXBcbiAgICB9IGVsc2Uge1xuICAgICAgIyBvbGQgYmVoYXZpb3IgZm9yIHVubmFtZWQgdmVjdG9yOiByZWN5Y2xlIG92ZXIgaW50ZXJhY3RpbmcgY2VsbCB0eXBlc1xuICAgICAgZGlyZWN0aW9uX2J5X2NsdXN0ZXIgPC0gcmVwKGRpcmVjdGlvbl9ieV9jbHVzdGVyLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGxlbmd0aC5vdXQgPSBsZW5ndGgoaW50ZXJhY3RfY2VsbF90eXBlcykpXG4gICAgICBuYW1lcyhkaXJlY3Rpb25fYnlfY2x1c3RlcikgPC0gaW50ZXJhY3RfY2VsbF90eXBlc1xuICAgIH1cbiAgfVxuXG4gICMgLS0tIEJhc2VsaW5lIGNvdW50cyBwZXIgKHNhbXBsZSwgY2VsbF90eXBlKSAtLS0tXG4gIG9uZV9zYW1wbGVfY291bnRzIDwtIGZ1bmN0aW9uKCkge1xuICAgIG1ham9yX2NvdW50cyA8LSByb3VuZChydW5pZihuX21ham9yX2NlbGxfdHlwZXMsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIG1pbiA9IG5fY2VsbHMgKiAoMSAtIHNkX2NlbGx0eXBlcyksXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIG1heCA9IG5fY2VsbHMgKiAoMSArIHNkX2NlbGx0eXBlcykpKVxuICAgIG1pbm9yX2NvdW50cyA8LSBpZiAobl9taW5vcl9jZWxsX3R5cGVzID4gMCkge1xuICAgICAgcm91bmQocnVuaWYobl9taW5vcl9jZWxsX3R5cGVzLFxuICAgICAgICAgICAgICAgICAgbWluID0gbl9jZWxscyAqIHJlbGF0aXZlX2FidW5kYW5jZSAqICgxIC0gc2RfY2VsbHR5cGVzKSxcbiAgICAgICAgICAgICAgICAgIG1heCA9IG5fY2VsbHMgKiByZWxhdGl2ZV9hYnVuZGFuY2UgKiAoMSArIHNkX2NlbGx0eXBlcykpKVxuICAgIH0gZWxzZSBpbnRlZ2VyKDApXG4gICAgYyhtYWpvcl9jb3VudHMsIG1pbm9yX2NvdW50cylcbiAgfVxuXG4gIGNvdW50c19saXN0IDwtIHJlcGxpY2F0ZShucm93KG1ldGEpLCBvbmVfc2FtcGxlX2NvdW50cygpLCBzaW1wbGlmeSA9IEZBTFNFKVxuICBjb3VudHNfZGYgPC0gZG8uY2FsbChyYmluZCwgY291bnRzX2xpc3QpXG4gIGNvbG5hbWVzKGNvdW50c19kZikgPC0gY2VsbF90eXBlc1xuXG4gIGNvdW50c19sb25nIDwtIHJlc2hhcGUoXG4gICAgZGF0YS5mcmFtZShzYW1wbGVfaWQgPSBtZXRhJHNhbXBsZV9pZCwgY291bnRzX2RmLCBjaGVjay5uYW1lcyA9IEZBTFNFKSxcbiAgICB2YXJ5aW5nID0gY2VsbF90eXBlcywgdi5uYW1lcyA9IFwiY291bnRcIiwgdGltZXZhciA9IFwiY2VsbF90eXBlXCIsXG4gICAgdGltZXMgPSBjZWxsX3R5cGVzLCBkaXJlY3Rpb24gPSBcImxvbmdcIlxuICApXG4gIHJvd25hbWVzKGNvdW50c19sb25nKSA8LSBOVUxMXG5cbiAgIyBNZXJnZSBkaXNlYXNlL3Zpc2l0IHNvIHdlIGNhbiBhcHBseSBlZmZlY3RzXG4gIGNvdW50c19sb25nIDwtIG1lcmdlKGNvdW50c19sb25nLFxuICAgICAgICAgICAgICAgICAgICAgICBtZXRhWywgYyhcInNhbXBsZV9pZFwiLCBcInZpc2l0XCIsIFwiZGlzZWFzZVwiKV0sXG4gICAgICAgICAgICAgICAgICAgICAgIGJ5ID0gXCJzYW1wbGVfaWRcIiwgc29ydCA9IEZBTFNFKVxuXG4gICMgLS0tIEFwcGx5IGVmZmVjdHMgb25seSB0byBpbnRlcmFjdGluZyBjZWxsIHR5cGVzIC0tLS1cbiAgaXNfaW50ZXJhY3RpbmcgPC0gY291bnRzX2xvbmckY2VsbF90eXBlICVpbiUgaW50ZXJhY3RfY2VsbF90eXBlc1xuICBlZmZfdmVjIDwtIG51bWVyaWMobnJvdyhjb3VudHNfbG9uZykpXG4gIGlmIChhbnkoaXNfaW50ZXJhY3RpbmcpKSB7XG4gICAgIyBtYXAgZGlyZWN0aW9uIHBlciBjZWxsIHR5cGVcbiAgICBkaXJfbWFwIDwtIHNldE5hbWVzKGRpcmVjdGlvbl9ieV9jbHVzdGVyLCBubSA9IG5hbWVzKGRpcmVjdGlvbl9ieV9jbHVzdGVyKSlcbiAgICBkaXJfY3QgIDwtIHVubmFtZShkaXJfbWFwW2NvdW50c19sb25nJGNlbGxfdHlwZVtpc19pbnRlcmFjdGluZ11dKVxuICAgICMgdmlzaXQgZWZmZWN0cyBieSBncm91cFxuICAgIHZfaWR4ICAgPC0gY291bnRzX2xvbmckdmlzaXRbaXNfaW50ZXJhY3RpbmddICsgMUxcbiAgICBpc19wcm9nIDwtIGNvdW50c19sb25nJGRpc2Vhc2VbaXNfaW50ZXJhY3RpbmddID09IDFMXG4gICAgdmUgICAgICA8LSBpZmVsc2UoaXNfcHJvZywgdmlzaXRfZWZmZWN0c19wcm9ncmVzc29yW3ZfaWR4XSwgdmlzaXRfZWZmZWN0c19jb250cm9sW3ZfaWR4XSlcbiAgICBlZmZfdmVjW2lzX2ludGVyYWN0aW5nXSA8LSBkaXJfY3QgKiB2ZVxuICB9XG5cbiAgY291bnRzX2xvbmckYWRqX2NvdW50IDwtIHBtYXgoMEwsIHJvdW5kKGNvdW50c19sb25nJGNvdW50ICogKDEgKyBlZmZfdmVjKSkpXG5cbiAgIyAtLS0gRXhwYW5kIHRvIHBlci1jZWxsIHJvd3MgYW5kIGF0dGFjaCBtZXRhZGF0YSAoSU5DTFVESU5HIGludGVyYWN0X3Rlcm0pIC0tLS1cbiAgcmVwX2VhY2ggPC0gZnVuY3Rpb24oeCwgdGltZXMpIGlmIChsZW5ndGgoeCkgPT0gMCkgeCBlbHNlIHJlcCh4LCB0aW1lcyA9IHRpbWVzKVxuICBleHBhbmRlZCA8LSBkYXRhLmZyYW1lKFxuICAgIHNhbXBsZV9pZCA9IHJlcF9lYWNoKGNvdW50c19sb25nJHNhbXBsZV9pZCwgY291bnRzX2xvbmckYWRqX2NvdW50KSxcbiAgICBjZWxsX3R5cGUgPSByZXBfZWFjaChjb3VudHNfbG9uZyRjZWxsX3R5cGUsIGNvdW50c19sb25nJGFkal9jb3VudCksXG4gICAgc3RyaW5nc0FzRmFjdG9ycyA9IEZBTFNFXG4gIClcblxuICAjIE1lcmdlIG1ldGFkYXRhOyBrZWVwIGludGVyYWN0aW9uICsgaW50ZXJhY3RfdGVybVxuICBrZWVwX2NvbHMgPC0gYyhcInNhbXBsZV9pZFwiLFwic3ViamVjdF9pZFwiLFwidmlzaXRcIixcInNleFwiLFwiZGlzZWFzZVwiLFwiYWdlXCIsXCJibWlcIixcImJhdGNoXCIsXG4gICAgICAgICAgICAgICAgIFwiaW50ZXJhY3Rpb25cIixcImludGVyYWN0X3Rlcm1cIilcbiAgZHVtbXlfZGF0YSA8LSBtZXJnZShleHBhbmRlZCwgbWV0YVssIGtlZXBfY29sc10sIGJ5ID0gXCJzYW1wbGVfaWRcIiwgc29ydCA9IEZBTFNFKVxuXG4gICMgU2h1ZmZsZSByb3dzIGZvciByZWFsaXNtXG4gIGlmIChucm93KGR1bW15X2RhdGEpID4gMSkge1xuICAgIGR1bW15X2RhdGEgPC0gZHVtbXlfZGF0YVtzYW1wbGUuaW50KG5yb3coZHVtbXlfZGF0YSkpLCAsIGRyb3AgPSBGQUxTRV1cbiAgICByb3duYW1lcyhkdW1teV9kYXRhKSA8LSBOVUxMXG4gIH1cblxuICBkdW1teV9kYXRhXG59XG5cblxuXG5cblxuYGBgIn0= -->

```r
#new function

#' Title
#'
#' @param n_cells
#' @param sd_celltypes
#' @param n_major_cell_types
#' @param n_minor_cell_types
#' @param relative_abundance
#' @param n_major_interact_celltypes
#' @param n_minor_interact_celltypes
#' @param n_individuals
#' @param n_batchs
#' @param interaction_feature
#' @param time_points
#' @param test_var
#' @param prop_disease
#' @param fc_interact
#' @param interaction_type
#' @param seed
#' @param visit_effects_progressor
#' @param visit_effects_control
#' @param direction_by_cluster
#'
#' @return
#' @export
#'
#' @examples
generate_dummy_data_new <- function(
  n_cells = 3000,                  # baseline cells per major cell type per sample
  sd_celltypes = 0.10,             # relative sd for counts
  n_major_cell_types = 7,
  n_minor_cell_types = 3,
  relative_abundance = 0.10,       # minor vs major baseline ratio
  n_major_interact_celltypes = 1,  # how many majors are \interacting\
  n_minor_interact_celltypes = 1,  # how many minors are \interacting\
  n_individuals = 30,
  n_batchs = 4,

  interaction_feature = \visit\,   # kept for labeling
  time_points = 4,                 # >= 2; works for 3+
  test_var = \disease\,
  prop_disease = 0.50,

  fc_interact = 0.10,              # effect magnitude used by defaults below
  interaction_type = c(\specific\,\differential\,\opposite\),
  seed = 1234,

  visit_effects_progressor = NULL, # multiplicative change: +0.1 means +10% vs baseline
  visit_effects_control    = NULL, # default 0s
  direction_by_cluster     = NULL  # +1/-1 for interacting clusters, recycled as needed
) {
  set.seed(seed)
  interaction_type <- match.arg(interaction_type)

  # --- Basic setup ----
  n_cell_types <- n_major_cell_types + n_minor_cell_types
  stopifnot(time_points >= 2, n_cell_types >= 1)

  cell_types <- LETTERS[seq_len(n_cell_types)]
  major_idx  <- seq_len(n_major_cell_types)
  minor_idx  <- if (n_minor_cell_types > 0) (n_major_cell_types + seq_len(n_minor_cell_types)) else integer(0)

  # which cell types are \interacting\ (first some majors, last some minors)
  interact_idx <- c(
    head(major_idx, n_major_interact_celltypes),
    tail(seq_len(n_cell_types), n_minor_interact_celltypes)
  )
  interact_idx <- intersect(interact_idx, seq_len(n_cell_types))  # guard rails
  interact_cell_types <- cell_types[interact_idx]

  # --- Subjects and visits ----
  subject_id <- paste0(\SUB_\, seq_len(n_individuals))
  disease_vec <- c(rep(1L, round(n_individuals * prop_disease)),
                   rep(0L, n_individuals - round(n_individuals * prop_disease)))
  disease_vec <- sample(disease_vec, n_individuals)

  sex_vec <- sample(c(0L, 1L), n_individuals, replace = TRUE)          # 0/1
  age_vec <- sample(18:60, n_individuals, replace = TRUE)
  bmi_vec <- sample(15:35, n_individuals, replace = TRUE)
  batch_vec <- rep(seq_len(n_batchs), length.out = n_individuals)

  subjects <- data.frame(
    subject_id = subject_id,
    sex   = sex_vec,
    disease = disease_vec,           # canonical exposure storage
    age   = age_vec,
    bmi   = bmi_vec,
    batch = factor(batch_vec),
    stringsAsFactors = FALSE
  )

  visits <- data.frame(
    subject_id = rep(subject_id, each = time_points),
    visit      = rep(0:(time_points - 1), times = n_individuals),
    stringsAsFactors = FALSE
  )
  visits$sample_id <- paste0(visits$subject_id, \_V\, visits$visit)

  meta <- merge(visits, subjects, by = \subject_id\, sort = FALSE)

  # Make sure a column named `test_var` exists (even if test_var != \disease\)
  if (!identical(test_var, \disease\)) {
    meta[[test_var]] <- meta[[\disease\]]
  }

  # Label and INTERACTION TERM (persist to output)
  meta$interaction    <- paste0(interaction_feature, \:\, test_var)
  meta$interact_term  <- as.integer(meta[[interaction_feature]]) * as.integer(meta[[test_var]])

  # --- Default visit effects (length = time_points) ----
  if (is.null(visit_effects_progressor)) {
    if (interaction_type == \specific\) {
      ve <- rep(0, time_points)
      if (time_points >= 2) ve[2] <- fc_interact  # bump V1 only
      visit_effects_progressor <- ve
    } else if (interaction_type == \differential\) {
      visit_effects_progressor <- rep(fc_interact, time_points)
    } else { # \opposite\: alternate +/- starting at V0
      ve <- rep(0, time_points)
      ve[seq(1, time_points, by = 2)] <- +fc_interact  # V0, V2, ...
      if (time_points >= 2) ve[seq(2, time_points, by = 2)] <- -fc_interact # V1, V3, ...
      visit_effects_progressor <- ve
    }
  }
  if (is.null(visit_effects_control)) {
    visit_effects_control <- rep(0, time_points)
  }
  stopifnot(length(visit_effects_progressor) == time_points,
            length(visit_effects_control)    == time_points)

  # Direction per interacting cluster
  if (is.null(direction_by_cluster)) {
    # old default: all interacting cell types go in the same direction (+1)
    direction_by_cluster <- rep(1L, max(1, length(interact_cell_types)))
    names(direction_by_cluster) <- interact_cell_types
  } else {
    if (!is.null(names(direction_by_cluster))) {
      # NEW: user passed a *named* vector, e.g. c(A = 1, C = -1, J = 1)
      # We align to interact_cell_types and default missing ones to +1.
      tmp <- rep(1L, max(1, length(interact_cell_types)))
      names(tmp) <- interact_cell_types

      matched <- intersect(names(direction_by_cluster), interact_cell_types)
      tmp[matched] <- direction_by_cluster[matched]

      direction_by_cluster <- tmp
    } else {
      # old behavior for unnamed vector: recycle over interacting cell types
      direction_by_cluster <- rep(direction_by_cluster,
                                  length.out = length(interact_cell_types))
      names(direction_by_cluster) <- interact_cell_types
    }
  }

  # --- Baseline counts per (sample, cell_type) ----
  one_sample_counts <- function() {
    major_counts <- round(runif(n_major_cell_types,
                                min = n_cells * (1 - sd_celltypes),
                                max = n_cells * (1 + sd_celltypes)))
    minor_counts <- if (n_minor_cell_types > 0) {
      round(runif(n_minor_cell_types,
                  min = n_cells * relative_abundance * (1 - sd_celltypes),
                  max = n_cells * relative_abundance * (1 + sd_celltypes)))
    } else integer(0)
    c(major_counts, minor_counts)
  }

  counts_list <- replicate(nrow(meta), one_sample_counts(), simplify = FALSE)
  counts_df <- do.call(rbind, counts_list)
  colnames(counts_df) <- cell_types

  counts_long <- reshape(
    data.frame(sample_id = meta$sample_id, counts_df, check.names = FALSE),
    varying = cell_types, v.names = \count\, timevar = \cell_type\,
    times = cell_types, direction = \long\
  )
  rownames(counts_long) <- NULL

  # Merge disease/visit so we can apply effects
  counts_long <- merge(counts_long,
                       meta[, c(\sample_id\, \visit\, \disease\)],
                       by = \sample_id\, sort = FALSE)

  # --- Apply effects only to interacting cell types ----
  is_interacting <- counts_long$cell_type %in% interact_cell_types
  eff_vec <- numeric(nrow(counts_long))
  if (any(is_interacting)) {
    # map direction per cell type
    dir_map <- setNames(direction_by_cluster, nm = names(direction_by_cluster))
    dir_ct  <- unname(dir_map[counts_long$cell_type[is_interacting]])
    # visit effects by group
    v_idx   <- counts_long$visit[is_interacting] + 1L
    is_prog <- counts_long$disease[is_interacting] == 1L
    ve      <- ifelse(is_prog, visit_effects_progressor[v_idx], visit_effects_control[v_idx])
    eff_vec[is_interacting] <- dir_ct * ve
  }

  counts_long$adj_count <- pmax(0L, round(counts_long$count * (1 + eff_vec)))

  # --- Expand to per-cell rows and attach metadata (INCLUDING interact_term) ----
  rep_each <- function(x, times) if (length(x) == 0) x else rep(x, times = times)
  expanded <- data.frame(
    sample_id = rep_each(counts_long$sample_id, counts_long$adj_count),
    cell_type = rep_each(counts_long$cell_type, counts_long$adj_count),
    stringsAsFactors = FALSE
  )

  # Merge metadata; keep interaction + interact_term
  keep_cols <- c(\sample_id\,\subject_id\,\visit\,\sex\,\disease\,\age\,\bmi\,\batch\,
                 \interaction\,\interact_term\)
  dummy_data <- merge(expanded, meta[, keep_cols], by = \sample_id\, sort = FALSE)

  # Shuffle rows for realism
  if (nrow(dummy_data) > 1) {
    dummy_data <- dummy_data[sample.int(nrow(dummy_data)), , drop = FALSE]
    rownames(dummy_data) <- NULL
  }

  dummy_data
}

```

<!-- rnb-source-end -->
````



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYGBgclxuI25ldyBmdW5jdGlvblxuXG4jJyBUaXRsZVxuIydcbiMnIEBwYXJhbSBuX2NlbGxzXG4jJyBAcGFyYW0gc2RfY2VsbHR5cGVzXG4jJyBAcGFyYW0gbl9tYWpvcl9jZWxsX3R5cGVzXG4jJyBAcGFyYW0gbl9taW5vcl9jZWxsX3R5cGVzXG4jJyBAcGFyYW0gcmVsYXRpdmVfYWJ1bmRhbmNlXG4jJyBAcGFyYW0gbl9tYWpvcl9pbnRlcmFjdF9jZWxsdHlwZXNcbiMnIEBwYXJhbSBuX21pbm9yX2ludGVyYWN0X2NlbGx0eXBlc1xuIycgQHBhcmFtIG5faW5kaXZpZHVhbHNcbiMnIEBwYXJhbSBuX2JhdGNoc1xuIycgQHBhcmFtIGludGVyYWN0aW9uX2ZlYXR1cmVcbiMnIEBwYXJhbSB0aW1lX3BvaW50c1xuIycgQHBhcmFtIHRlc3RfdmFyXG4jJyBAcGFyYW0gcHJvcF9kaXNlYXNlXG4jJyBAcGFyYW0gZmNfaW50ZXJhY3RcbiMnIEBwYXJhbSBpbnRlcmFjdGlvbl90eXBlXG4jJyBAcGFyYW0gc2VlZFxuIycgQHBhcmFtIHZpc2l0X2VmZmVjdHNfcHJvZ3Jlc3NvclxuIycgQHBhcmFtIHZpc2l0X2VmZmVjdHNfY29udHJvbFxuIycgQHBhcmFtIGRpcmVjdGlvbl9ieV9jbHVzdGVyXG4jJ1xuIycgQHJldHVyblxuIycgQGV4cG9ydFxuIydcbiMnIEBleGFtcGxlc1xuZ2VuZXJhdGVfZHVtbXlfZGF0YV9uZXcgPC0gZnVuY3Rpb24oXG4gIG5fY2VsbHMgPSAzMDAwLCAgICAgICAgICAgICAgICAgICMgYmFzZWxpbmUgY2VsbHMgcGVyIG1ham9yIGNlbGwgdHlwZSBwZXIgc2FtcGxlXG4gIHNkX2NlbGx0eXBlcyA9IDAuMTAsICAgICAgICAgICAgICMgcmVsYXRpdmUgc2QgZm9yIGNvdW50c1xuICBuX21ham9yX2NlbGxfdHlwZXMgPSA3LFxuICBuX21pbm9yX2NlbGxfdHlwZXMgPSAzLFxuICByZWxhdGl2ZV9hYnVuZGFuY2UgPSAwLjEwLCAgICAgICAjIG1pbm9yIHZzIG1ham9yIGJhc2VsaW5lIHJhdGlvXG4gIG5fbWFqb3JfaW50ZXJhY3RfY2VsbHR5cGVzID0gMSwgICMgaG93IG1hbnkgbWFqb3JzIGFyZSBcXGludGVyYWN0aW5nXFxcbiAgbl9taW5vcl9pbnRlcmFjdF9jZWxsdHlwZXMgPSAxLCAgIyBob3cgbWFueSBtaW5vcnMgYXJlIFxcaW50ZXJhY3RpbmdcXFxuICBuX2luZGl2aWR1YWxzID0gMzAsXG4gIG5fYmF0Y2hzID0gNCxcblxuICBpbnRlcmFjdGlvbl9mZWF0dXJlID0gXFx2aXNpdFxcLCAgICMga2VwdCBmb3IgbGFiZWxpbmdcbiAgdGltZV9wb2ludHMgPSA0LCAgICAgICAgICAgICAgICAgIyA+PSAyOyB3b3JrcyBmb3IgMytcbiAgdGVzdF92YXIgPSBcXGRpc2Vhc2VcXCxcbiAgcHJvcF9kaXNlYXNlID0gMC41MCxcblxuICBmY19pbnRlcmFjdCA9IDAuMTAsICAgICAgICAgICAgICAjIGVmZmVjdCBtYWduaXR1ZGUgdXNlZCBieSBkZWZhdWx0cyBiZWxvd1xuICBpbnRlcmFjdGlvbl90eXBlID0gYyhcXHNwZWNpZmljXFwsXFxkaWZmZXJlbnRpYWxcXCxcXG9wcG9zaXRlXFwpLFxuICBzZWVkID0gMTIzNCxcblxuICB2aXNpdF9lZmZlY3RzX3Byb2dyZXNzb3IgPSBOVUxMLCAjIG11bHRpcGxpY2F0aXZlIGNoYW5nZTogKzAuMSBtZWFucyArMTAlIHZzIGJhc2VsaW5lXG4gIHZpc2l0X2VmZmVjdHNfY29udHJvbCAgICA9IE5VTEwsICMgZGVmYXVsdCAwc1xuICBkaXJlY3Rpb25fYnlfY2x1c3RlciAgICAgPSBOVUxMICAjICsxLy0xIGZvciBpbnRlcmFjdGluZyBjbHVzdGVycywgcmVjeWNsZWQgYXMgbmVlZGVkXG4pIHtcbiAgc2V0LnNlZWQoc2VlZClcbiAgaW50ZXJhY3Rpb25fdHlwZSA8LSBtYXRjaC5hcmcoaW50ZXJhY3Rpb25fdHlwZSlcblxuICAjIC0tLSBCYXNpYyBzZXR1cCAtLS0tXG4gIG5fY2VsbF90eXBlcyA8LSBuX21ham9yX2NlbGxfdHlwZXMgKyBuX21pbm9yX2NlbGxfdHlwZXNcbiAgc3RvcGlmbm90KHRpbWVfcG9pbnRzID49IDIsIG5fY2VsbF90eXBlcyA+PSAxKVxuXG4gIGNlbGxfdHlwZXMgPC0gTEVUVEVSU1tzZXFfbGVuKG5fY2VsbF90eXBlcyldXG4gIG1ham9yX2lkeCAgPC0gc2VxX2xlbihuX21ham9yX2NlbGxfdHlwZXMpXG4gIG1pbm9yX2lkeCAgPC0gaWYgKG5fbWlub3JfY2VsbF90eXBlcyA+IDApIChuX21ham9yX2NlbGxfdHlwZXMgKyBzZXFfbGVuKG5fbWlub3JfY2VsbF90eXBlcykpIGVsc2UgaW50ZWdlcigwKVxuXG4gICMgd2hpY2ggY2VsbCB0eXBlcyBhcmUgXFxpbnRlcmFjdGluZ1xcIChmaXJzdCBzb21lIG1ham9ycywgbGFzdCBzb21lIG1pbm9ycylcbiAgaW50ZXJhY3RfaWR4IDwtIGMoXG4gICAgaGVhZChtYWpvcl9pZHgsIG5fbWFqb3JfaW50ZXJhY3RfY2VsbHR5cGVzKSxcbiAgICB0YWlsKHNlcV9sZW4obl9jZWxsX3R5cGVzKSwgbl9taW5vcl9pbnRlcmFjdF9jZWxsdHlwZXMpXG4gIClcbiAgaW50ZXJhY3RfaWR4IDwtIGludGVyc2VjdChpbnRlcmFjdF9pZHgsIHNlcV9sZW4obl9jZWxsX3R5cGVzKSkgICMgZ3VhcmQgcmFpbHNcbiAgaW50ZXJhY3RfY2VsbF90eXBlcyA8LSBjZWxsX3R5cGVzW2ludGVyYWN0X2lkeF1cblxuICAjIC0tLSBTdWJqZWN0cyBhbmQgdmlzaXRzIC0tLS1cbiAgc3ViamVjdF9pZCA8LSBwYXN0ZTAoXFxTVUJfXFwsIHNlcV9sZW4obl9pbmRpdmlkdWFscykpXG4gIGRpc2Vhc2VfdmVjIDwtIGMocmVwKDFMLCByb3VuZChuX2luZGl2aWR1YWxzICogcHJvcF9kaXNlYXNlKSksXG4gICAgICAgICAgICAgICAgICAgcmVwKDBMLCBuX2luZGl2aWR1YWxzIC0gcm91bmQobl9pbmRpdmlkdWFscyAqIHByb3BfZGlzZWFzZSkpKVxuICBkaXNlYXNlX3ZlYyA8LSBzYW1wbGUoZGlzZWFzZV92ZWMsIG5faW5kaXZpZHVhbHMpXG5cbiAgc2V4X3ZlYyA8LSBzYW1wbGUoYygwTCwgMUwpLCBuX2luZGl2aWR1YWxzLCByZXBsYWNlID0gVFJVRSkgICAgICAgICAgIyAwLzFcbiAgYWdlX3ZlYyA8LSBzYW1wbGUoMTg6NjAsIG5faW5kaXZpZHVhbHMsIHJlcGxhY2UgPSBUUlVFKVxuICBibWlfdmVjIDwtIHNhbXBsZSgxNTozNSwgbl9pbmRpdmlkdWFscywgcmVwbGFjZSA9IFRSVUUpXG4gIGJhdGNoX3ZlYyA8LSByZXAoc2VxX2xlbihuX2JhdGNocyksIGxlbmd0aC5vdXQgPSBuX2luZGl2aWR1YWxzKVxuXG4gIHN1YmplY3RzIDwtIGRhdGEuZnJhbWUoXG4gICAgc3ViamVjdF9pZCA9IHN1YmplY3RfaWQsXG4gICAgc2V4ICAgPSBzZXhfdmVjLFxuICAgIGRpc2Vhc2UgPSBkaXNlYXNlX3ZlYywgICAgICAgICAgICMgY2Fub25pY2FsIGV4cG9zdXJlIHN0b3JhZ2VcbiAgICBhZ2UgICA9IGFnZV92ZWMsXG4gICAgYm1pICAgPSBibWlfdmVjLFxuICAgIGJhdGNoID0gZmFjdG9yKGJhdGNoX3ZlYyksXG4gICAgc3RyaW5nc0FzRmFjdG9ycyA9IEZBTFNFXG4gIClcblxuICB2aXNpdHMgPC0gZGF0YS5mcmFtZShcbiAgICBzdWJqZWN0X2lkID0gcmVwKHN1YmplY3RfaWQsIGVhY2ggPSB0aW1lX3BvaW50cyksXG4gICAgdmlzaXQgICAgICA9IHJlcCgwOih0aW1lX3BvaW50cyAtIDEpLCB0aW1lcyA9IG5faW5kaXZpZHVhbHMpLFxuICAgIHN0cmluZ3NBc0ZhY3RvcnMgPSBGQUxTRVxuICApXG4gIHZpc2l0cyRzYW1wbGVfaWQgPC0gcGFzdGUwKHZpc2l0cyRzdWJqZWN0X2lkLCBcXF9WXFwsIHZpc2l0cyR2aXNpdClcblxuICBtZXRhIDwtIG1lcmdlKHZpc2l0cywgc3ViamVjdHMsIGJ5ID0gXFxzdWJqZWN0X2lkXFwsIHNvcnQgPSBGQUxTRSlcblxuICAjIE1ha2Ugc3VyZSBhIGNvbHVtbiBuYW1lZCBgdGVzdF92YXJgIGV4aXN0cyAoZXZlbiBpZiB0ZXN0X3ZhciAhPSBcXGRpc2Vhc2VcXClcbiAgaWYgKCFpZGVudGljYWwodGVzdF92YXIsIFxcZGlzZWFzZVxcKSkge1xuICAgIG1ldGFbW3Rlc3RfdmFyXV0gPC0gbWV0YVtbXFxkaXNlYXNlXFxdXVxuICB9XG5cbiAgIyBMYWJlbCBhbmQgSU5URVJBQ1RJT04gVEVSTSAocGVyc2lzdCB0byBvdXRwdXQpXG4gIG1ldGEkaW50ZXJhY3Rpb24gICAgPC0gcGFzdGUwKGludGVyYWN0aW9uX2ZlYXR1cmUsIFxcOlxcLCB0ZXN0X3ZhcilcbiAgbWV0YSRpbnRlcmFjdF90ZXJtICA8LSBhcy5pbnRlZ2VyKG1ldGFbW2ludGVyYWN0aW9uX2ZlYXR1cmVdXSkgKiBhcy5pbnRlZ2VyKG1ldGFbW3Rlc3RfdmFyXV0pXG5cbiAgIyAtLS0gRGVmYXVsdCB2aXNpdCBlZmZlY3RzIChsZW5ndGggPSB0aW1lX3BvaW50cykgLS0tLVxuICBpZiAoaXMubnVsbCh2aXNpdF9lZmZlY3RzX3Byb2dyZXNzb3IpKSB7XG4gICAgaWYgKGludGVyYWN0aW9uX3R5cGUgPT0gXFxzcGVjaWZpY1xcKSB7XG4gICAgICB2ZSA8LSByZXAoMCwgdGltZV9wb2ludHMpXG4gICAgICBpZiAodGltZV9wb2ludHMgPj0gMikgdmVbMl0gPC0gZmNfaW50ZXJhY3QgICMgYnVtcCBWMSBvbmx5XG4gICAgICB2aXNpdF9lZmZlY3RzX3Byb2dyZXNzb3IgPC0gdmVcbiAgICB9IGVsc2UgaWYgKGludGVyYWN0aW9uX3R5cGUgPT0gXFxkaWZmZXJlbnRpYWxcXCkge1xuICAgICAgdmlzaXRfZWZmZWN0c19wcm9ncmVzc29yIDwtIHJlcChmY19pbnRlcmFjdCwgdGltZV9wb2ludHMpXG4gICAgfSBlbHNlIHsgIyBcXG9wcG9zaXRlXFw6IGFsdGVybmF0ZSArLy0gc3RhcnRpbmcgYXQgVjBcbiAgICAgIHZlIDwtIHJlcCgwLCB0aW1lX3BvaW50cylcbiAgICAgIHZlW3NlcSgxLCB0aW1lX3BvaW50cywgYnkgPSAyKV0gPC0gK2ZjX2ludGVyYWN0ICAjIFYwLCBWMiwgLi4uXG4gICAgICBpZiAodGltZV9wb2ludHMgPj0gMikgdmVbc2VxKDIsIHRpbWVfcG9pbnRzLCBieSA9IDIpXSA8LSAtZmNfaW50ZXJhY3QgIyBWMSwgVjMsIC4uLlxuICAgICAgdmlzaXRfZWZmZWN0c19wcm9ncmVzc29yIDwtIHZlXG4gICAgfVxuICB9XG4gIGlmIChpcy5udWxsKHZpc2l0X2VmZmVjdHNfY29udHJvbCkpIHtcbiAgICB2aXNpdF9lZmZlY3RzX2NvbnRyb2wgPC0gcmVwKDAsIHRpbWVfcG9pbnRzKVxuICB9XG4gIHN0b3BpZm5vdChsZW5ndGgodmlzaXRfZWZmZWN0c19wcm9ncmVzc29yKSA9PSB0aW1lX3BvaW50cyxcbiAgICAgICAgICAgIGxlbmd0aCh2aXNpdF9lZmZlY3RzX2NvbnRyb2wpICAgID09IHRpbWVfcG9pbnRzKVxuXG4gICMgRGlyZWN0aW9uIHBlciBpbnRlcmFjdGluZyBjbHVzdGVyXG4gIGlmIChpcy5udWxsKGRpcmVjdGlvbl9ieV9jbHVzdGVyKSkge1xuICAgICMgb2xkIGRlZmF1bHQ6IGFsbCBpbnRlcmFjdGluZyBjZWxsIHR5cGVzIGdvIGluIHRoZSBzYW1lIGRpcmVjdGlvbiAoKzEpXG4gICAgZGlyZWN0aW9uX2J5X2NsdXN0ZXIgPC0gcmVwKDFMLCBtYXgoMSwgbGVuZ3RoKGludGVyYWN0X2NlbGxfdHlwZXMpKSlcbiAgICBuYW1lcyhkaXJlY3Rpb25fYnlfY2x1c3RlcikgPC0gaW50ZXJhY3RfY2VsbF90eXBlc1xuICB9IGVsc2Uge1xuICAgIGlmICghaXMubnVsbChuYW1lcyhkaXJlY3Rpb25fYnlfY2x1c3RlcikpKSB7XG4gICAgICAjIE5FVzogdXNlciBwYXNzZWQgYSAqbmFtZWQqIHZlY3RvciwgZS5nLiBjKEEgPSAxLCBDID0gLTEsIEogPSAxKVxuICAgICAgIyBXZSBhbGlnbiB0byBpbnRlcmFjdF9jZWxsX3R5cGVzIGFuZCBkZWZhdWx0IG1pc3Npbmcgb25lcyB0byArMS5cbiAgICAgIHRtcCA8LSByZXAoMUwsIG1heCgxLCBsZW5ndGgoaW50ZXJhY3RfY2VsbF90eXBlcykpKVxuICAgICAgbmFtZXModG1wKSA8LSBpbnRlcmFjdF9jZWxsX3R5cGVzXG5cbiAgICAgIG1hdGNoZWQgPC0gaW50ZXJzZWN0KG5hbWVzKGRpcmVjdGlvbl9ieV9jbHVzdGVyKSwgaW50ZXJhY3RfY2VsbF90eXBlcylcbiAgICAgIHRtcFttYXRjaGVkXSA8LSBkaXJlY3Rpb25fYnlfY2x1c3RlclttYXRjaGVkXVxuXG4gICAgICBkaXJlY3Rpb25fYnlfY2x1c3RlciA8LSB0bXBcbiAgICB9IGVsc2Uge1xuICAgICAgIyBvbGQgYmVoYXZpb3IgZm9yIHVubmFtZWQgdmVjdG9yOiByZWN5Y2xlIG92ZXIgaW50ZXJhY3RpbmcgY2VsbCB0eXBlc1xuICAgICAgZGlyZWN0aW9uX2J5X2NsdXN0ZXIgPC0gcmVwKGRpcmVjdGlvbl9ieV9jbHVzdGVyLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIGxlbmd0aC5vdXQgPSBsZW5ndGgoaW50ZXJhY3RfY2VsbF90eXBlcykpXG4gICAgICBuYW1lcyhkaXJlY3Rpb25fYnlfY2x1c3RlcikgPC0gaW50ZXJhY3RfY2VsbF90eXBlc1xuICAgIH1cbiAgfVxuXG4gICMgLS0tIEJhc2VsaW5lIGNvdW50cyBwZXIgKHNhbXBsZSwgY2VsbF90eXBlKSAtLS0tXG4gIG9uZV9zYW1wbGVfY291bnRzIDwtIGZ1bmN0aW9uKCkge1xuICAgIG1ham9yX2NvdW50cyA8LSByb3VuZChydW5pZihuX21ham9yX2NlbGxfdHlwZXMsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIG1pbiA9IG5fY2VsbHMgKiAoMSAtIHNkX2NlbGx0eXBlcyksXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIG1heCA9IG5fY2VsbHMgKiAoMSArIHNkX2NlbGx0eXBlcykpKVxuICAgIG1pbm9yX2NvdW50cyA8LSBpZiAobl9taW5vcl9jZWxsX3R5cGVzID4gMCkge1xuICAgICAgcm91bmQocnVuaWYobl9taW5vcl9jZWxsX3R5cGVzLFxuICAgICAgICAgICAgICAgICAgbWluID0gbl9jZWxscyAqIHJlbGF0aXZlX2FidW5kYW5jZSAqICgxIC0gc2RfY2VsbHR5cGVzKSxcbiAgICAgICAgICAgICAgICAgIG1heCA9IG5fY2VsbHMgKiByZWxhdGl2ZV9hYnVuZGFuY2UgKiAoMSArIHNkX2NlbGx0eXBlcykpKVxuICAgIH0gZWxzZSBpbnRlZ2VyKDApXG4gICAgYyhtYWpvcl9jb3VudHMsIG1pbm9yX2NvdW50cylcbiAgfVxuXG4gIGNvdW50c19saXN0IDwtIHJlcGxpY2F0ZShucm93KG1ldGEpLCBvbmVfc2FtcGxlX2NvdW50cygpLCBzaW1wbGlmeSA9IEZBTFNFKVxuICBjb3VudHNfZGYgPC0gZG8uY2FsbChyYmluZCwgY291bnRzX2xpc3QpXG4gIGNvbG5hbWVzKGNvdW50c19kZikgPC0gY2VsbF90eXBlc1xuXG4gIGNvdW50c19sb25nIDwtIHJlc2hhcGUoXG4gICAgZGF0YS5mcmFtZShzYW1wbGVfaWQgPSBtZXRhJHNhbXBsZV9pZCwgY291bnRzX2RmLCBjaGVjay5uYW1lcyA9IEZBTFNFKSxcbiAgICB2YXJ5aW5nID0gY2VsbF90eXBlcywgdi5uYW1lcyA9IFxcY291bnRcXCwgdGltZXZhciA9IFxcY2VsbF90eXBlXFwsXG4gICAgdGltZXMgPSBjZWxsX3R5cGVzLCBkaXJlY3Rpb24gPSBcXGxvbmdcXFxuICApXG4gIHJvd25hbWVzKGNvdW50c19sb25nKSA8LSBOVUxMXG5cbiAgIyBNZXJnZSBkaXNlYXNlL3Zpc2l0IHNvIHdlIGNhbiBhcHBseSBlZmZlY3RzXG4gIGNvdW50c19sb25nIDwtIG1lcmdlKGNvdW50c19sb25nLFxuICAgICAgICAgICAgICAgICAgICAgICBtZXRhWywgYyhcXHNhbXBsZV9pZFxcLCBcXHZpc2l0XFwsIFxcZGlzZWFzZVxcKV0sXG4gICAgICAgICAgICAgICAgICAgICAgIGJ5ID0gXFxzYW1wbGVfaWRcXCwgc29ydCA9IEZBTFNFKVxuXG4gICMgLS0tIEFwcGx5IGVmZmVjdHMgb25seSB0byBpbnRlcmFjdGluZyBjZWxsIHR5cGVzIC0tLS1cbiAgaXNfaW50ZXJhY3RpbmcgPC0gY291bnRzX2xvbmckY2VsbF90eXBlICVpbiUgaW50ZXJhY3RfY2VsbF90eXBlc1xuICBlZmZfdmVjIDwtIG51bWVyaWMobnJvdyhjb3VudHNfbG9uZykpXG4gIGlmIChhbnkoaXNfaW50ZXJhY3RpbmcpKSB7XG4gICAgIyBtYXAgZGlyZWN0aW9uIHBlciBjZWxsIHR5cGVcbiAgICBkaXJfbWFwIDwtIHNldE5hbWVzKGRpcmVjdGlvbl9ieV9jbHVzdGVyLCBubSA9IG5hbWVzKGRpcmVjdGlvbl9ieV9jbHVzdGVyKSlcbiAgICBkaXJfY3QgIDwtIHVubmFtZShkaXJfbWFwW2NvdW50c19sb25nJGNlbGxfdHlwZVtpc19pbnRlcmFjdGluZ11dKVxuICAgICMgdmlzaXQgZWZmZWN0cyBieSBncm91cFxuICAgIHZfaWR4ICAgPC0gY291bnRzX2xvbmckdmlzaXRbaXNfaW50ZXJhY3RpbmddICsgMUxcbiAgICBpc19wcm9nIDwtIGNvdW50c19sb25nJGRpc2Vhc2VbaXNfaW50ZXJhY3RpbmddID09IDFMXG4gICAgdmUgICAgICA8LSBpZmVsc2UoaXNfcHJvZywgdmlzaXRfZWZmZWN0c19wcm9ncmVzc29yW3ZfaWR4XSwgdmlzaXRfZWZmZWN0c19jb250cm9sW3ZfaWR4XSlcbiAgICBlZmZfdmVjW2lzX2ludGVyYWN0aW5nXSA8LSBkaXJfY3QgKiB2ZVxuICB9XG5cbiAgY291bnRzX2xvbmckYWRqX2NvdW50IDwtIHBtYXgoMEwsIHJvdW5kKGNvdW50c19sb25nJGNvdW50ICogKDEgKyBlZmZfdmVjKSkpXG5cbiAgIyAtLS0gRXhwYW5kIHRvIHBlci1jZWxsIHJvd3MgYW5kIGF0dGFjaCBtZXRhZGF0YSAoSU5DTFVESU5HIGludGVyYWN0X3Rlcm0pIC0tLS1cbiAgcmVwX2VhY2ggPC0gZnVuY3Rpb24oeCwgdGltZXMpIGlmIChsZW5ndGgoeCkgPT0gMCkgeCBlbHNlIHJlcCh4LCB0aW1lcyA9IHRpbWVzKVxuICBleHBhbmRlZCA8LSBkYXRhLmZyYW1lKFxuICAgIHNhbXBsZV9pZCA9IHJlcF9lYWNoKGNvdW50c19sb25nJHNhbXBsZV9pZCwgY291bnRzX2xvbmckYWRqX2NvdW50KSxcbiAgICBjZWxsX3R5cGUgPSByZXBfZWFjaChjb3VudHNfbG9uZyRjZWxsX3R5cGUsIGNvdW50c19sb25nJGFkal9jb3VudCksXG4gICAgc3RyaW5nc0FzRmFjdG9ycyA9IEZBTFNFXG4gIClcblxuICAjIE1lcmdlIG1ldGFkYXRhOyBrZWVwIGludGVyYWN0aW9uICsgaW50ZXJhY3RfdGVybVxuICBrZWVwX2NvbHMgPC0gYyhcXHNhbXBsZV9pZFxcLFxcc3ViamVjdF9pZFxcLFxcdmlzaXRcXCxcXHNleFxcLFxcZGlzZWFzZVxcLFxcYWdlXFwsXFxibWlcXCxcXGJhdGNoXFwsXG4gICAgICAgICAgICAgICAgIFxcaW50ZXJhY3Rpb25cXCxcXGludGVyYWN0X3Rlcm1cXClcbiAgZHVtbXlfZGF0YSA8LSBtZXJnZShleHBhbmRlZCwgbWV0YVssIGtlZXBfY29sc10sIGJ5ID0gXFxzYW1wbGVfaWRcXCwgc29ydCA9IEZBTFNFKVxuXG4gICMgU2h1ZmZsZSByb3dzIGZvciByZWFsaXNtXG4gIGlmIChucm93KGR1bW15X2RhdGEpID4gMSkge1xuICAgIGR1bW15X2RhdGEgPC0gZHVtbXlfZGF0YVtzYW1wbGUuaW50KG5yb3coZHVtbXlfZGF0YSkpLCAsIGRyb3AgPSBGQUxTRV1cbiAgICByb3duYW1lcyhkdW1teV9kYXRhKSA8LSBOVUxMXG4gIH1cblxuICBkdW1teV9kYXRhXG59XG5cblxuXG5cblxuYGBgXG5gYGAifQ== -->

```r
```r
#new function

#' Title
#'
#' @param n_cells
#' @param sd_celltypes
#' @param n_major_cell_types
#' @param n_minor_cell_types
#' @param relative_abundance
#' @param n_major_interact_celltypes
#' @param n_minor_interact_celltypes
#' @param n_individuals
#' @param n_batchs
#' @param interaction_feature
#' @param time_points
#' @param test_var
#' @param prop_disease
#' @param fc_interact
#' @param interaction_type
#' @param seed
#' @param visit_effects_progressor
#' @param visit_effects_control
#' @param direction_by_cluster
#'
#' @return
#' @export
#'
#' @examples
generate_dummy_data_new <- function(
  n_cells = 3000,                  # baseline cells per major cell type per sample
  sd_celltypes = 0.10,             # relative sd for counts
  n_major_cell_types = 7,
  n_minor_cell_types = 3,
  relative_abundance = 0.10,       # minor vs major baseline ratio
  n_major_interact_celltypes = 1,  # how many majors are \interacting\
  n_minor_interact_celltypes = 1,  # how many minors are \interacting\
  n_individuals = 30,
  n_batchs = 4,

  interaction_feature = \visit\,   # kept for labeling
  time_points = 4,                 # >= 2; works for 3+
  test_var = \disease\,
  prop_disease = 0.50,

  fc_interact = 0.10,              # effect magnitude used by defaults below
  interaction_type = c(\specific\,\differential\,\opposite\),
  seed = 1234,

  visit_effects_progressor = NULL, # multiplicative change: +0.1 means +10% vs baseline
  visit_effects_control    = NULL, # default 0s
  direction_by_cluster     = NULL  # +1/-1 for interacting clusters, recycled as needed
) {
  set.seed(seed)
  interaction_type <- match.arg(interaction_type)

  # --- Basic setup ----
  n_cell_types <- n_major_cell_types + n_minor_cell_types
  stopifnot(time_points >= 2, n_cell_types >= 1)

  cell_types <- LETTERS[seq_len(n_cell_types)]
  major_idx  <- seq_len(n_major_cell_types)
  minor_idx  <- if (n_minor_cell_types > 0) (n_major_cell_types + seq_len(n_minor_cell_types)) else integer(0)

  # which cell types are \interacting\ (first some majors, last some minors)
  interact_idx <- c(
    head(major_idx, n_major_interact_celltypes),
    tail(seq_len(n_cell_types), n_minor_interact_celltypes)
  )
  interact_idx <- intersect(interact_idx, seq_len(n_cell_types))  # guard rails
  interact_cell_types <- cell_types[interact_idx]

  # --- Subjects and visits ----
  subject_id <- paste0(\SUB_\, seq_len(n_individuals))
  disease_vec <- c(rep(1L, round(n_individuals * prop_disease)),
                   rep(0L, n_individuals - round(n_individuals * prop_disease)))
  disease_vec <- sample(disease_vec, n_individuals)

  sex_vec <- sample(c(0L, 1L), n_individuals, replace = TRUE)          # 0/1
  age_vec <- sample(18:60, n_individuals, replace = TRUE)
  bmi_vec <- sample(15:35, n_individuals, replace = TRUE)
  batch_vec <- rep(seq_len(n_batchs), length.out = n_individuals)

  subjects <- data.frame(
    subject_id = subject_id,
    sex   = sex_vec,
    disease = disease_vec,           # canonical exposure storage
    age   = age_vec,
    bmi   = bmi_vec,
    batch = factor(batch_vec),
    stringsAsFactors = FALSE
  )

  visits <- data.frame(
    subject_id = rep(subject_id, each = time_points),
    visit      = rep(0:(time_points - 1), times = n_individuals),
    stringsAsFactors = FALSE
  )
  visits$sample_id <- paste0(visits$subject_id, \_V\, visits$visit)

  meta <- merge(visits, subjects, by = \subject_id\, sort = FALSE)

  # Make sure a column named `test_var` exists (even if test_var != \disease\)
  if (!identical(test_var, \disease\)) {
    meta[[test_var]] <- meta[[\disease\]]
  }

  # Label and INTERACTION TERM (persist to output)
  meta$interaction    <- paste0(interaction_feature, \:\, test_var)
  meta$interact_term  <- as.integer(meta[[interaction_feature]]) * as.integer(meta[[test_var]])

  # --- Default visit effects (length = time_points) ----
  if (is.null(visit_effects_progressor)) {
    if (interaction_type == \specific\) {
      ve <- rep(0, time_points)
      if (time_points >= 2) ve[2] <- fc_interact  # bump V1 only
      visit_effects_progressor <- ve
    } else if (interaction_type == \differential\) {
      visit_effects_progressor <- rep(fc_interact, time_points)
    } else { # \opposite\: alternate +/- starting at V0
      ve <- rep(0, time_points)
      ve[seq(1, time_points, by = 2)] <- +fc_interact  # V0, V2, ...
      if (time_points >= 2) ve[seq(2, time_points, by = 2)] <- -fc_interact # V1, V3, ...
      visit_effects_progressor <- ve
    }
  }
  if (is.null(visit_effects_control)) {
    visit_effects_control <- rep(0, time_points)
  }
  stopifnot(length(visit_effects_progressor) == time_points,
            length(visit_effects_control)    == time_points)

  # Direction per interacting cluster
  if (is.null(direction_by_cluster)) {
    # old default: all interacting cell types go in the same direction (+1)
    direction_by_cluster <- rep(1L, max(1, length(interact_cell_types)))
    names(direction_by_cluster) <- interact_cell_types
  } else {
    if (!is.null(names(direction_by_cluster))) {
      # NEW: user passed a *named* vector, e.g. c(A = 1, C = -1, J = 1)
      # We align to interact_cell_types and default missing ones to +1.
      tmp <- rep(1L, max(1, length(interact_cell_types)))
      names(tmp) <- interact_cell_types

      matched <- intersect(names(direction_by_cluster), interact_cell_types)
      tmp[matched] <- direction_by_cluster[matched]

      direction_by_cluster <- tmp
    } else {
      # old behavior for unnamed vector: recycle over interacting cell types
      direction_by_cluster <- rep(direction_by_cluster,
                                  length.out = length(interact_cell_types))
      names(direction_by_cluster) <- interact_cell_types
    }
  }

  # --- Baseline counts per (sample, cell_type) ----
  one_sample_counts <- function() {
    major_counts <- round(runif(n_major_cell_types,
                                min = n_cells * (1 - sd_celltypes),
                                max = n_cells * (1 + sd_celltypes)))
    minor_counts <- if (n_minor_cell_types > 0) {
      round(runif(n_minor_cell_types,
                  min = n_cells * relative_abundance * (1 - sd_celltypes),
                  max = n_cells * relative_abundance * (1 + sd_celltypes)))
    } else integer(0)
    c(major_counts, minor_counts)
  }

  counts_list <- replicate(nrow(meta), one_sample_counts(), simplify = FALSE)
  counts_df <- do.call(rbind, counts_list)
  colnames(counts_df) <- cell_types

  counts_long <- reshape(
    data.frame(sample_id = meta$sample_id, counts_df, check.names = FALSE),
    varying = cell_types, v.names = \count\, timevar = \cell_type\,
    times = cell_types, direction = \long\
  )
  rownames(counts_long) <- NULL

  # Merge disease/visit so we can apply effects
  counts_long <- merge(counts_long,
                       meta[, c(\sample_id\, \visit\, \disease\)],
                       by = \sample_id\, sort = FALSE)

  # --- Apply effects only to interacting cell types ----
  is_interacting <- counts_long$cell_type %in% interact_cell_types
  eff_vec <- numeric(nrow(counts_long))
  if (any(is_interacting)) {
    # map direction per cell type
    dir_map <- setNames(direction_by_cluster, nm = names(direction_by_cluster))
    dir_ct  <- unname(dir_map[counts_long$cell_type[is_interacting]])
    # visit effects by group
    v_idx   <- counts_long$visit[is_interacting] + 1L
    is_prog <- counts_long$disease[is_interacting] == 1L
    ve      <- ifelse(is_prog, visit_effects_progressor[v_idx], visit_effects_control[v_idx])
    eff_vec[is_interacting] <- dir_ct * ve
  }

  counts_long$adj_count <- pmax(0L, round(counts_long$count * (1 + eff_vec)))

  # --- Expand to per-cell rows and attach metadata (INCLUDING interact_term) ----
  rep_each <- function(x, times) if (length(x) == 0) x else rep(x, times = times)
  expanded <- data.frame(
    sample_id = rep_each(counts_long$sample_id, counts_long$adj_count),
    cell_type = rep_each(counts_long$cell_type, counts_long$adj_count),
    stringsAsFactors = FALSE
  )

  # Merge metadata; keep interaction + interact_term
  keep_cols <- c(\sample_id\,\subject_id\,\visit\,\sex\,\disease\,\age\,\bmi\,\batch\,
                 \interaction\,\interact_term\)
  dummy_data <- merge(expanded, meta[, keep_cols], by = \sample_id\, sort = FALSE)

  # Shuffle rows for realism
  if (nrow(dummy_data) > 1) {
    dummy_data <- dummy_data[sample.int(nrow(dummy_data)), , drop = FALSE]
    rownames(dummy_data) <- NULL
  }

  dummy_data
}


<!-- rnb-source-end -->


<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->





<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-output-begin eyJkYXRhIjoiXG48IS0tIHJuYi1zb3VyY2UtYmVnaW4gZXlKa1lYUmhJam9pWUdCZ2NseHVjMlYwTG5ObFpXUW9NakF5TUNsY2JuUmxjM1JmYm1WM0lEd3RJR2RsYm1WeVlYUmxYMlIxYlcxNVgyUmhkR0ZmYm1WM0tGeHVJQ0J1WDJObGJHeHpJRDBnTVRVd0xGeHVJQ0J6WkY5alpXeHNkSGx3WlhNZ1BTQXdMakVzWEc0Z0lHNWZiV0ZxYjNKZlkyVnNiRjkwZVhCbGN5QTlJRGNzWEc0Z0lHNWZiV2x1YjNKZlkyVnNiRjkwZVhCbGN5QTlJRE1zWEc0Z0lISmxiR0YwYVhabFgyRmlkVzVrWVc1alpTQTlJREF1TkN4Y2JpQWdibDl0WVdwdmNsOXBiblJsY21GamRGOWpaV3hzZEhsd1pYTWdQU0F5TEZ4dUlDQnVYMjFwYm05eVgybHVkR1Z5WVdOMFgyTmxiR3gwZVhCbGN5QTlJRElzWEc0Z0lHNWZhVzVrYVhacFpIVmhiSE1nUFNBek1DeGNiaUFnYmw5aVlYUmphSE1nUFNBMExGeHVJQ0JwYm5SbGNtRmpkR2x2Ymw5bVpXRjBkWEpsSUQwZ1hDSjJhWE5wZEZ3aUxGeHVJQ0IwYVcxbFgzQnZhVzUwY3lBOUlETXNYRzRnSUhSbGMzUmZkbUZ5SUQwZ1hDSmthWE5sWVhObFhDSXNYRzRnSUhCeWIzQmZaR2x6WldGelpTQTlJREF1TlN4Y2JpQWdhVzUwWlhKaFkzUnBiMjVmZEhsd1pTQTlJRndpYzNCbFkybG1hV05jSWl4Y2JpQWdjMlZsWkNBOUlESXdNakFzWEc0Z0lIWnBjMmwwWDJWbVptVmpkSE5mY0hKdlozSmxjM052Y2lBOUlHTW9NQzR3TUN3Z01DNDJNQ3dnTUM0ek1Da3NYRzRnSUhacGMybDBYMlZtWm1WamRITmZZMjl1ZEhKdmJDQWdJQ0E5SUdNb01Dd2dNQ3dnTUNrc1hHNWNiaUFnSXlBZ1pHbHlaV04wYVc5dWN5Qm1iM0lnZEdobElHbHVkR1Z5WVdOMGFXNW5JR05sYkd3Z2RIbHdaWE1nUVN3Z1Fpd2dTU3dnU2x4dUlDQmthWEpsWTNScGIyNWZZbmxmWTJ4MWMzUmxjaUFnSUNBZ1BTQmpLRUVnUFNBck1Td2dJQ01nWlc1eWFXTm9JRUZjYmlBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQkNJRDBnTFRFc0lDQWpJR1JsY0d4bGRHVWdRbHh1SUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lFa2dQU0FyTVN3Z0lDTWdaVzV5YVdOb0lFbGNiaUFnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNBZ0lDQWdJQ0FnSUNCS0lEMGdMVEVwSUNBaklHUmxjR3hsZEdVZ1NseHVLVnh1WUdCZ0luMD0gLS0+XG5cbmBgYHJcbnNldC5zZWVkKDIwMjApXG50ZXN0X25ldyA8LSBnZW5lcmF0ZV9kdW1teV9kYXRhX25ldyhcbiAgbl9jZWxscyA9IDE1MCxcbiAgc2RfY2VsbHR5cGVzID0gMC4xLFxuICBuX21ham9yX2NlbGxfdHlwZXMgPSA3LFxuICBuX21pbm9yX2NlbGxfdHlwZXMgPSAzLFxuICByZWxhdGl2ZV9hYnVuZGFuY2UgPSAwLjQsXG4gIG5fbWFqb3JfaW50ZXJhY3RfY2VsbHR5cGVzID0gMixcbiAgbl9taW5vcl9pbnRlcmFjdF9jZWxsdHlwZXMgPSAyLFxuICBuX2luZGl2aWR1YWxzID0gMzAsXG4gIG5fYmF0Y2hzID0gNCxcbiAgaW50ZXJhY3Rpb25fZmVhdHVyZSA9IFxcdmlzaXRcXCxcbiAgdGltZV9wb2ludHMgPSAzLFxuICB0ZXN0X3ZhciA9IFxcZGlzZWFzZVxcLFxuICBwcm9wX2Rpc2Vhc2UgPSAwLjUsXG4gIGludGVyYWN0aW9uX3R5cGUgPSBcXHNwZWNpZmljXFwsXG4gIHNlZWQgPSAyMDIwLFxuICB2aXNpdF9lZmZlY3RzX3Byb2dyZXNzb3IgPSBjKDAuMDAsIDAuNjAsIDAuMzApLFxuICB2aXNpdF9lZmZlY3RzX2NvbnRyb2wgICAgPSBjKDAsIDAsIDApLFxuXG4gICMgIGRpcmVjdGlvbnMgZm9yIHRoZSBpbnRlcmFjdGluZyBjZWxsIHR5cGVzIEEsIEIsIEksIEpcbiAgZGlyZWN0aW9uX2J5X2NsdXN0ZXIgICAgID0gYyhBID0gKzEsICAjIGVucmljaCBBXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgQiA9IC0xLCAgIyBkZXBsZXRlIEJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBJID0gKzEsICAjIGVucmljaCBJXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgSiA9IC0xKSAgIyBkZXBsZXRlIEpcbilcbmBgYFxuXG48IS0tIHJuYi1zb3VyY2UtZW5kIC0tPlxuIn0= -->

````

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuc2V0LnNlZWQoMjAyMClcbnRlc3RfbmV3IDwtIGdlbmVyYXRlX2R1bW15X2RhdGFfbmV3KFxuICBuX2NlbGxzID0gMTUwLFxuICBzZF9jZWxsdHlwZXMgPSAwLjEsXG4gIG5fbWFqb3JfY2VsbF90eXBlcyA9IDcsXG4gIG5fbWlub3JfY2VsbF90eXBlcyA9IDMsXG4gIHJlbGF0aXZlX2FidW5kYW5jZSA9IDAuNCxcbiAgbl9tYWpvcl9pbnRlcmFjdF9jZWxsdHlwZXMgPSAyLFxuICBuX21pbm9yX2ludGVyYWN0X2NlbGx0eXBlcyA9IDIsXG4gIG5faW5kaXZpZHVhbHMgPSAzMCxcbiAgbl9iYXRjaHMgPSA0LFxuICBpbnRlcmFjdGlvbl9mZWF0dXJlID0gXCJ2aXNpdFwiLFxuICB0aW1lX3BvaW50cyA9IDMsXG4gIHRlc3RfdmFyID0gXCJkaXNlYXNlXCIsXG4gIHByb3BfZGlzZWFzZSA9IDAuNSxcbiAgaW50ZXJhY3Rpb25fdHlwZSA9IFwic3BlY2lmaWNcIixcbiAgc2VlZCA9IDIwMjAsXG4gIHZpc2l0X2VmZmVjdHNfcHJvZ3Jlc3NvciA9IGMoMC4wMCwgMC42MCwgMC4zMCksXG4gIHZpc2l0X2VmZmVjdHNfY29udHJvbCAgICA9IGMoMCwgMCwgMCksXG5cbiAgIyAgZGlyZWN0aW9ucyBmb3IgdGhlIGludGVyYWN0aW5nIGNlbGwgdHlwZXMgQSwgQiwgSSwgSlxuICBkaXJlY3Rpb25fYnlfY2x1c3RlciAgICAgPSBjKEEgPSArMSwgICMgZW5yaWNoIEFcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBCID0gLTEsICAjIGRlcGxldGUgQlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIEkgPSArMSwgICMgZW5yaWNoIElcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBKID0gLTEpICAjIGRlcGxldGUgSlxuKVxuYGBgIn0= -->

```r
set.seed(2020)
test_new <- generate_dummy_data_new(
  n_cells = 150,
  sd_celltypes = 0.1,
  n_major_cell_types = 7,
  n_minor_cell_types = 3,
  relative_abundance = 0.4,
  n_major_interact_celltypes = 2,
  n_minor_interact_celltypes = 2,
  n_individuals = 30,
  n_batchs = 4,
  interaction_feature = \visit\,
  time_points = 3,
  test_var = \disease\,
  prop_disease = 0.5,
  interaction_type = \specific\,
  seed = 2020,
  visit_effects_progressor = c(0.00, 0.60, 0.30),
  visit_effects_control    = c(0, 0, 0),

  #  directions for the interacting cell types A, B, I, J
  direction_by_cluster     = c(A = +1,  # enrich A
                               B = -1,  # deplete B
                               I = +1,  # enrich I
                               J = -1)  # deplete J
)
```

<!-- rnb-source-end -->
````



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYGBgclxuc2V0LnNlZWQoMjAyMClcbnRlc3RfbmV3IDwtIGdlbmVyYXRlX2R1bW15X2RhdGFfbmV3KFxuICBuX2NlbGxzID0gMTUwLFxuICBzZF9jZWxsdHlwZXMgPSAwLjEsXG4gIG5fbWFqb3JfY2VsbF90eXBlcyA9IDcsXG4gIG5fbWlub3JfY2VsbF90eXBlcyA9IDMsXG4gIHJlbGF0aXZlX2FidW5kYW5jZSA9IDAuNCxcbiAgbl9tYWpvcl9pbnRlcmFjdF9jZWxsdHlwZXMgPSAyLFxuICBuX21pbm9yX2ludGVyYWN0X2NlbGx0eXBlcyA9IDIsXG4gIG5faW5kaXZpZHVhbHMgPSAzMCxcbiAgbl9iYXRjaHMgPSA0LFxuICBpbnRlcmFjdGlvbl9mZWF0dXJlID0gXFx2aXNpdFxcLFxuICB0aW1lX3BvaW50cyA9IDMsXG4gIHRlc3RfdmFyID0gXFxkaXNlYXNlXFwsXG4gIHByb3BfZGlzZWFzZSA9IDAuNSxcbiAgaW50ZXJhY3Rpb25fdHlwZSA9IFxcc3BlY2lmaWNcXCxcbiAgc2VlZCA9IDIwMjAsXG4gIHZpc2l0X2VmZmVjdHNfcHJvZ3Jlc3NvciA9IGMoMC4wMCwgMC42MCwgMC4zMCksXG4gIHZpc2l0X2VmZmVjdHNfY29udHJvbCAgICA9IGMoMCwgMCwgMCksXG5cbiAgIyAgZGlyZWN0aW9ucyBmb3IgdGhlIGludGVyYWN0aW5nIGNlbGwgdHlwZXMgQSwgQiwgSSwgSlxuICBkaXJlY3Rpb25fYnlfY2x1c3RlciAgICAgPSBjKEEgPSArMSwgICMgZW5yaWNoIEFcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBCID0gLTEsICAjIGRlcGxldGUgQlxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIEkgPSArMSwgICMgZW5yaWNoIElcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBKID0gLTEpICAjIGRlcGxldGUgSlxuKVxuYGBgXG5gYGAifQ== -->

```r
```r
set.seed(2020)
test_new <- generate_dummy_data_new(
  n_cells = 150,
  sd_celltypes = 0.1,
  n_major_cell_types = 7,
  n_minor_cell_types = 3,
  relative_abundance = 0.4,
  n_major_interact_celltypes = 2,
  n_minor_interact_celltypes = 2,
  n_individuals = 30,
  n_batchs = 4,
  interaction_feature = \visit\,
  time_points = 3,
  test_var = \disease\,
  prop_disease = 0.5,
  interaction_type = \specific\,
  seed = 2020,
  visit_effects_progressor = c(0.00, 0.60, 0.30),
  visit_effects_control    = c(0, 0, 0),

  #  directions for the interacting cell types A, B, I, J
  direction_by_cluster     = c(A = +1,  # enrich A
                               B = -1,  # deplete B
                               I = +1,  # enrich I
                               J = -1)  # deplete J
)

<!-- rnb-source-end -->


<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->





<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-output-begin eyJkYXRhIjoiXG48IS0tIHJuYi1zb3VyY2UtYmVnaW4gZXlKa1lYUmhJam9pWUdCZ2NseHVkR1Z6ZEY5dVpYZGNiblJsYzNSZmIyeGtYRzVnWUdBaWZRPT0gLS0+XG5cbmBgYHJcbnRlc3RfbmV3XG50ZXN0X29sZFxuYGBgXG5cbjwhLS0gcm5iLXNvdXJjZS1lbmQgLS0+XG4ifQ== -->

````

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudGVzdF9uZXdcbnRlc3Rfb2xkXG5gYGAifQ== -->

```r
test_new
test_old
```

<!-- rnb-source-end -->
````



<!-- rnb-output-end -->

<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYGBgclxudGVzdF9uZXdcbnRlc3Rfb2xkXG5gYGBcbmBgYCJ9 -->

```r
```r
test_new
test_old

<!-- rnb-source-end -->


<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



The function below is the data production but this one will need a design matrix


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIycgU2ltdWxhdGUgbG9uZ2l0dWRpbmFsIHNpbmdsZS1jZWxsIGRhdGEgd2l0aCBpbnRlcmFjdGluZyBjZWxsIHR5cGVzXG4jJ1xuIycgQHBhcmFtIG5fY2VsbHMgQmFzZWxpbmUgY2VsbHMgcGVyIG1ham9yIGNlbGwgdHlwZSBwZXIgc2FtcGxlLlxuIycgQHBhcmFtIHNkX2NlbGx0eXBlcyBSZWxhdGl2ZSBTRCBmb3IgY291bnRzLlxuIycgQHBhcmFtIG5fbWFqb3JfY2VsbF90eXBlcyBOdW1iZXIgb2YgbWFqb3IgY2VsbCB0eXBlcy5cbiMnIEBwYXJhbSBuX21pbm9yX2NlbGxfdHlwZXMgTnVtYmVyIG9mIG1pbm9yIGNlbGwgdHlwZXMuXG4jJyBAcGFyYW0gcmVsYXRpdmVfYWJ1bmRhbmNlIEJhc2VsaW5lIG1pbm9yOm1ham9yIHJhdGlvLlxuIycgQHBhcmFtIG5fbWFqb3JfaW50ZXJhY3RfY2VsbHR5cGVzIEhvdyBtYW55IG1ham9yIHR5cGVzIGFyZSBcImludGVyYWN0aW5nXCIuXG4jJyBAcGFyYW0gbl9taW5vcl9pbnRlcmFjdF9jZWxsdHlwZXMgSG93IG1hbnkgbWlub3IgdHlwZXMgYXJlIFwiaW50ZXJhY3RpbmdcIi5cbiMnIEBwYXJhbSBuX2luZGl2aWR1YWxzIE51bWJlciBvZiBzdWJqZWN0cy5cbiMnIEBwYXJhbSBuX2JhdGNocyBOdW1iZXIgb2YgYmF0Y2hlcy5cbiMnIEBwYXJhbSBpbnRlcmFjdGlvbl9mZWF0dXJlIE5hbWUgb2YgdGhlIHRpbWUvdmlzaXQgdmFyaWFibGUgKGUuZy4gXCJ2aXNpdFwiKS5cbiMnIEBwYXJhbSB0aW1lX3BvaW50cyBOdW1iZXIgb2YgbG9uZ2l0dWRpbmFsIHRpbWUgcG9pbnRzICg+PSAyKS5cbiMnIEBwYXJhbSB0ZXN0X3ZhciBOYW1lIG9mIHRoZSBkaXNlYXNlL2V4cG9zdXJlIHZhcmlhYmxlIChlLmcuIFwiZGlzZWFzZVwiKS5cbiMnIEBwYXJhbSBwcm9wX2Rpc2Vhc2UgUHJvcG9ydGlvbiBvZiBkaXNlYXNlZCBzdWJqZWN0cy5cbiMnIEBwYXJhbSBmY19pbnRlcmFjdCBFZmZlY3QgbWFnbml0dWRlIHVzZWQgYnkgdGhlIGRlZmF1bHQgZWZmZWN0IHBhdHRlcm5zLlxuIycgQHBhcmFtIGludGVyYWN0aW9uX3R5cGUgT25lIG9mIFwic3BlY2lmaWNcIiwgXCJkaWZmZXJlbnRpYWxcIiwgb3IgXCJvcHBvc2l0ZVwiXG4jJyAgIGZvciB0aGUgZGVmYXVsdCBsb25naXR1ZGluYWwgcGF0dGVybiAod2hlbiBkZXNpZ24gbWF0cmljZXMgYXJlIG5vdCBnaXZlbikuXG4jJyBAcGFyYW0gc2VlZCBSYW5kb20gc2VlZC5cbiMnIEBwYXJhbSB2aXNpdF9lZmZlY3RzX3Byb2dyZXNzb3IgT3B0aW9uYWwgdmVjdG9yIG9mIGxlbmd0aCBgdGltZV9wb2ludHNgXG4jJyAgIGZvciBkZWZhdWx0IHByb2dyZXNzb3IgdmlzaXQgZWZmZWN0cyAoaWdub3JlZCBpZiBkZXNpZ24gbWF0cmljZXMgZ2l2ZW4pLlxuIycgQHBhcmFtIHZpc2l0X2VmZmVjdHNfY29udHJvbCBPcHRpb25hbCB2ZWN0b3Igb2YgbGVuZ3RoIGB0aW1lX3BvaW50c2BcbiMnICAgZm9yIGRlZmF1bHQgY29udHJvbCB2aXNpdCBlZmZlY3RzIChpZ25vcmVkIGlmIGRlc2lnbiBtYXRyaWNlcyBnaXZlbikuXG4jJyBAcGFyYW0gZGlyZWN0aW9uX2J5X2NsdXN0ZXIgT3B0aW9uYWwgdmVjdG9yIG9mICsxIC8gLTEgcGVyIGludGVyYWN0aW5nXG4jJyAgIGNlbGwgdHlwZSAoaWdub3JlZCBpZiBkZXNpZ24gbWF0cmljZXMgZ2l2ZW4pLlxuIycgQHBhcmFtIGVmZmVjdF9tYXRfcHJvZ3Jlc3NvciBPcHRpb25hbCBudW1lcmljIG1hdHJpeCBvZiBkaW1lbnNpb25cbiMnICAgYG5fY2VsbF90eXBlcyB4IHRpbWVfcG9pbnRzYCBnaXZpbmcgZWZmZWN0IHNpemVzIGZvciBwcm9ncmVzc29ycy5cbiMnICAgUm93cyBjb3JyZXNwb25kIHRvIGNlbGwgdHlwZXMgKGJ5IG5hbWUgb3Igb3JkZXIpLCBjb2x1bW5zIHRvIHZpc2l0c1xuIycgICAoMC4udGltZV9wb2ludHMtMSkuIElmIHByb3ZpZGVkIChhbG9uZyB3aXRoIGBlZmZlY3RfbWF0X2NvbnRyb2xgKSxcbiMnICAgb3ZlcnJpZGVzIGB2aXNpdF9lZmZlY3RzXypgIGFuZCBgZGlyZWN0aW9uX2J5X2NsdXN0ZXJgLlxuIycgQHBhcmFtIGVmZmVjdF9tYXRfY29udHJvbCBPcHRpb25hbCBudW1lcmljIG1hdHJpeCBvZiBzYW1lIGRpbWVuc2lvbiBhc1xuIycgICBgZWZmZWN0X21hdF9wcm9ncmVzc29yYCBnaXZpbmcgZWZmZWN0cyBmb3IgY29udHJvbHMuXG4jJ1xuIycgQHJldHVybiBBIGRhdGEuZnJhbWUgb2Ygc2ltdWxhdGVkIHNpbmdsZS1jZWxsIG1ldGFkYXRhLlxuIycgQGV4cG9ydFxuZ2VuZXJhdGVfZHVtbXlfZGF0YSA8LSBmdW5jdGlvbihcbiAgbl9jZWxscyA9IDMwMDAsICAgICAgICAgICAgICAgICAgIyBiYXNlbGluZSBjZWxscyBwZXIgbWFqb3IgY2VsbCB0eXBlIHBlciBzYW1wbGVcbiAgc2RfY2VsbHR5cGVzID0gMC4xMCwgICAgICAgICAgICAgIyByZWxhdGl2ZSBzZCBmb3IgY291bnRzXG4gIG5fbWFqb3JfY2VsbF90eXBlcyA9IDcsXG4gIG5fbWlub3JfY2VsbF90eXBlcyA9IDMsXG4gIHJlbGF0aXZlX2FidW5kYW5jZSA9IDAuMTAsICAgICAgICMgbWlub3IgdnMgbWFqb3IgYmFzZWxpbmUgcmF0aW9cbiAgbl9tYWpvcl9pbnRlcmFjdF9jZWxsdHlwZXMgPSAxLCAgIyBob3cgbWFueSBtYWpvcnMgYXJlIFwiaW50ZXJhY3RpbmdcIlxuICBuX21pbm9yX2ludGVyYWN0X2NlbGx0eXBlcyA9IDEsICAjIGhvdyBtYW55IG1pbm9ycyBhcmUgXCJpbnRlcmFjdGluZ1wiXG4gIG5faW5kaXZpZHVhbHMgPSAzMCxcbiAgbl9iYXRjaHMgPSA0LFxuXG4gIGludGVyYWN0aW9uX2ZlYXR1cmUgPSBcInZpc2l0XCIsICAgIyBrZXB0IGZvciBsYWJlbGluZ1xuICB0aW1lX3BvaW50cyA9IDQsICAgICAgICAgICAgICAgICAjID49IDI7IHdvcmtzIGZvciAzK1xuICB0ZXN0X3ZhciA9IFwiZGlzZWFzZVwiLFxuICBwcm9wX2Rpc2Vhc2UgPSAwLjUwLFxuXG4gIGZjX2ludGVyYWN0ID0gMC4xMCwgICAgICAgICAgICAgICMgZWZmZWN0IG1hZ25pdHVkZSB1c2VkIGJ5IGRlZmF1bHRzIGJlbG93XG4gIGludGVyYWN0aW9uX3R5cGUgPSBjKFwic3BlY2lmaWNcIixcImRpZmZlcmVudGlhbFwiLFwib3Bwb3NpdGVcIiksXG4gIHNlZWQgPSAxMjM0LFxuXG4gIHZpc2l0X2VmZmVjdHNfcHJvZ3Jlc3NvciA9IE5VTEwsICMgbXVsdGlwbGljYXRpdmUgY2hhbmdlOiArMC4xIG1lYW5zICsxMCUgdnMgYmFzZWxpbmVcbiAgdmlzaXRfZWZmZWN0c19jb250cm9sICAgID0gTlVMTCwgIyBkZWZhdWx0IDBzXG4gIGRpcmVjdGlvbl9ieV9jbHVzdGVyICAgICA9IE5VTEwsICMgKzEvLTEgZm9yIGludGVyYWN0aW5nIGNsdXN0ZXJzLCByZWN5Y2xlZCBhcyBuZWVkZWRcblxuICBlZmZlY3RfbWF0X3Byb2dyZXNzb3IgICAgPSBOVUxMLCAjIE5FVzogY2VsbF90eXBlIHggdGltZV9wb2ludHMgbWF0cml4IGZvciBwcm9ncmVzc29yc1xuICBlZmZlY3RfbWF0X2NvbnRyb2wgICAgICAgPSBOVUxMICAjIE5FVzogc2FtZSBmb3IgY29udHJvbHNcbikge1xuICBzZXQuc2VlZChzZWVkKVxuICBpbnRlcmFjdGlvbl90eXBlIDwtIG1hdGNoLmFyZyhpbnRlcmFjdGlvbl90eXBlKVxuXG4gICMgLS0tIEJhc2ljIHNldHVwIC0tLS1cbiAgbl9jZWxsX3R5cGVzIDwtIG5fbWFqb3JfY2VsbF90eXBlcyArIG5fbWlub3JfY2VsbF90eXBlc1xuICBzdG9waWZub3QodGltZV9wb2ludHMgPj0gMiwgbl9jZWxsX3R5cGVzID49IDEpXG5cbiAgY2VsbF90eXBlcyA8LSBMRVRURVJTW3NlcV9sZW4obl9jZWxsX3R5cGVzKV1cbiAgbWFqb3JfaWR4ICA8LSBzZXFfbGVuKG5fbWFqb3JfY2VsbF90eXBlcylcbiAgbWlub3JfaWR4ICA8LSBpZiAobl9taW5vcl9jZWxsX3R5cGVzID4gMCkgKG5fbWFqb3JfY2VsbF90eXBlcyArIHNlcV9sZW4obl9taW5vcl9jZWxsX3R5cGVzKSkgZWxzZSBpbnRlZ2VyKDApXG5cbiAgIyB3aGljaCBjZWxsIHR5cGVzIGFyZSBcImludGVyYWN0aW5nXCIgKGZpcnN0IHNvbWUgbWFqb3JzLCBsYXN0IHNvbWUgbWlub3JzKVxuICBpbnRlcmFjdF9pZHggPC0gYyhcbiAgICBoZWFkKG1ham9yX2lkeCwgbl9tYWpvcl9pbnRlcmFjdF9jZWxsdHlwZXMpLFxuICAgIHRhaWwoc2VxX2xlbihuX2NlbGxfdHlwZXMpLCBuX21pbm9yX2ludGVyYWN0X2NlbGx0eXBlcylcbiAgKVxuICBpbnRlcmFjdF9pZHggPC0gaW50ZXJzZWN0KGludGVyYWN0X2lkeCwgc2VxX2xlbihuX2NlbGxfdHlwZXMpKSAgIyBndWFyZCByYWlsc1xuICBpbnRlcmFjdF9jZWxsX3R5cGVzIDwtIGNlbGxfdHlwZXNbaW50ZXJhY3RfaWR4XVxuXG4gICMgLS0tIFN1YmplY3RzIGFuZCB2aXNpdHMgLS0tLVxuICBzdWJqZWN0X2lkIDwtIHBhc3RlMChcIlNVQl9cIiwgc2VxX2xlbihuX2luZGl2aWR1YWxzKSlcbiAgZGlzZWFzZV92ZWMgPC0gYyhyZXAoMUwsIHJvdW5kKG5faW5kaXZpZHVhbHMgKiBwcm9wX2Rpc2Vhc2UpKSxcbiAgICAgICAgICAgICAgICAgICByZXAoMEwsIG5faW5kaXZpZHVhbHMgLSByb3VuZChuX2luZGl2aWR1YWxzICogcHJvcF9kaXNlYXNlKSkpXG4gIGRpc2Vhc2VfdmVjIDwtIHNhbXBsZShkaXNlYXNlX3ZlYywgbl9pbmRpdmlkdWFscylcblxuICBzZXhfdmVjIDwtIHNhbXBsZShjKDBMLCAxTCksIG5faW5kaXZpZHVhbHMsIHJlcGxhY2UgPSBUUlVFKSAgICAgICAgICAjIDAvMVxuICBhZ2VfdmVjIDwtIHNhbXBsZSgxODo2MCwgbl9pbmRpdmlkdWFscywgcmVwbGFjZSA9IFRSVUUpXG4gIGJtaV92ZWMgPC0gc2FtcGxlKDE1OjM1LCBuX2luZGl2aWR1YWxzLCByZXBsYWNlID0gVFJVRSlcbiAgYmF0Y2hfdmVjIDwtIHJlcChzZXFfbGVuKG5fYmF0Y2hzKSwgbGVuZ3RoLm91dCA9IG5faW5kaXZpZHVhbHMpXG5cbiAgc3ViamVjdHMgPC0gZGF0YS5mcmFtZShcbiAgICBzdWJqZWN0X2lkID0gc3ViamVjdF9pZCxcbiAgICBzZXggICA9IHNleF92ZWMsXG4gICAgZGlzZWFzZSA9IGRpc2Vhc2VfdmVjLCAgICAgICAgICAgIyBjYW5vbmljYWwgZXhwb3N1cmUgc3RvcmFnZVxuICAgIGFnZSAgID0gYWdlX3ZlYyxcbiAgICBibWkgICA9IGJtaV92ZWMsXG4gICAgYmF0Y2ggPSBmYWN0b3IoYmF0Y2hfdmVjKSxcbiAgICBzdHJpbmdzQXNGYWN0b3JzID0gRkFMU0VcbiAgKVxuXG4gIHZpc2l0cyA8LSBkYXRhLmZyYW1lKFxuICAgIHN1YmplY3RfaWQgPSByZXAoc3ViamVjdF9pZCwgZWFjaCA9IHRpbWVfcG9pbnRzKSxcbiAgICB2aXNpdCAgICAgID0gcmVwKDA6KHRpbWVfcG9pbnRzIC0gMSksIHRpbWVzID0gbl9pbmRpdmlkdWFscyksXG4gICAgc3RyaW5nc0FzRmFjdG9ycyA9IEZBTFNFXG4gIClcbiAgdmlzaXRzJHNhbXBsZV9pZCA8LSBwYXN0ZTAodmlzaXRzJHN1YmplY3RfaWQsIFwiX1ZcIiwgdmlzaXRzJHZpc2l0KVxuXG4gIG1ldGEgPC0gbWVyZ2UodmlzaXRzLCBzdWJqZWN0cywgYnkgPSBcInN1YmplY3RfaWRcIiwgc29ydCA9IEZBTFNFKVxuXG4gICMgTWFrZSBzdXJlIGEgY29sdW1uIG5hbWVkIGB0ZXN0X3ZhcmAgZXhpc3RzIChldmVuIGlmIHRlc3RfdmFyICE9IFwiZGlzZWFzZVwiKVxuICBpZiAoIWlkZW50aWNhbCh0ZXN0X3ZhciwgXCJkaXNlYXNlXCIpKSB7XG4gICAgbWV0YVtbdGVzdF92YXJdXSA8LSBtZXRhW1tcImRpc2Vhc2VcIl1dXG4gIH1cblxuICAjIExhYmVsIGFuZCBJTlRFUkFDVElPTiBURVJNIChwZXJzaXN0IHRvIG91dHB1dClcbiAgbWV0YSRpbnRlcmFjdGlvbiAgICA8LSBwYXN0ZTAoaW50ZXJhY3Rpb25fZmVhdHVyZSwgXCI6XCIsIHRlc3RfdmFyKVxuICBtZXRhJGludGVyYWN0X3Rlcm0gIDwtIGFzLmludGVnZXIobWV0YVtbaW50ZXJhY3Rpb25fZmVhdHVyZV1dKSAqIGFzLmludGVnZXIobWV0YVtbdGVzdF92YXJdXSlcblxuICAjIC0tLSBCYXNlbGluZSBjb3VudHMgcGVyIChzYW1wbGUsIGNlbGxfdHlwZSkgLS0tLVxuICBvbmVfc2FtcGxlX2NvdW50cyA8LSBmdW5jdGlvbigpIHtcbiAgICBtYWpvcl9jb3VudHMgPC0gcm91bmQocnVuaWYobl9tYWpvcl9jZWxsX3R5cGVzLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBtaW4gPSBuX2NlbGxzICogKDEgLSBzZF9jZWxsdHlwZXMpLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBtYXggPSBuX2NlbGxzICogKDEgKyBzZF9jZWxsdHlwZXMpKSlcbiAgICBtaW5vcl9jb3VudHMgPC0gaWYgKG5fbWlub3JfY2VsbF90eXBlcyA+IDApIHtcbiAgICAgIHJvdW5kKHJ1bmlmKG5fbWlub3JfY2VsbF90eXBlcyxcbiAgICAgICAgICAgICAgICAgIG1pbiA9IG5fY2VsbHMgKiByZWxhdGl2ZV9hYnVuZGFuY2UgKiAoMSAtIHNkX2NlbGx0eXBlcyksXG4gICAgICAgICAgICAgICAgICBtYXggPSBuX2NlbGxzICogcmVsYXRpdmVfYWJ1bmRhbmNlICogKDEgKyBzZF9jZWxsdHlwZXMpKSlcbiAgICB9IGVsc2UgaW50ZWdlcigwKVxuICAgIGMobWFqb3JfY291bnRzLCBtaW5vcl9jb3VudHMpXG4gIH1cblxuICBjb3VudHNfbGlzdCA8LSByZXBsaWNhdGUobnJvdyhtZXRhKSwgb25lX3NhbXBsZV9jb3VudHMoKSwgc2ltcGxpZnkgPSBGQUxTRSlcbiAgY291bnRzX2RmIDwtIGRvLmNhbGwocmJpbmQsIGNvdW50c19saXN0KVxuICBjb2xuYW1lcyhjb3VudHNfZGYpIDwtIGNlbGxfdHlwZXNcblxuICBjb3VudHNfbG9uZyA8LSByZXNoYXBlKFxuICAgIGRhdGEuZnJhbWUoc2FtcGxlX2lkID0gbWV0YSRzYW1wbGVfaWQsIGNvdW50c19kZiwgY2hlY2submFtZXMgPSBGQUxTRSksXG4gICAgdmFyeWluZyA9IGNlbGxfdHlwZXMsIHYubmFtZXMgPSBcImNvdW50XCIsIHRpbWV2YXIgPSBcImNlbGxfdHlwZVwiLFxuICAgIHRpbWVzID0gY2VsbF90eXBlcywgZGlyZWN0aW9uID0gXCJsb25nXCJcbiAgKVxuICByb3duYW1lcyhjb3VudHNfbG9uZykgPC0gTlVMTFxuXG4gICMgTWVyZ2UgZGlzZWFzZS92aXNpdCBzbyB3ZSBjYW4gYXBwbHkgZWZmZWN0c1xuICBjb3VudHNfbG9uZyA8LSBtZXJnZShjb3VudHNfbG9uZyxcbiAgICAgICAgICAgICAgICAgICAgICAgbWV0YVssIGMoXCJzYW1wbGVfaWRcIiwgXCJ2aXNpdFwiLCBcImRpc2Vhc2VcIildLFxuICAgICAgICAgICAgICAgICAgICAgICBieSA9IFwic2FtcGxlX2lkXCIsIHNvcnQgPSBGQUxTRSlcblxuICAjIC0tLSBCdWlsZCBlZmZlY3Qgc3RydWN0dXJlIC0tLS1cbiAgZWZmX3ZlYyA8LSBudW1lcmljKG5yb3coY291bnRzX2xvbmcpKVxuXG4gIGlmICghaXMubnVsbChlZmZlY3RfbWF0X3Byb2dyZXNzb3IpIHx8ICFpcy5udWxsKGVmZmVjdF9tYXRfY29udHJvbCkpIHtcbiAgICAjIyAtLS0tIE5FVzogdXNlIGRlc2lnbiBtYXRyaWNlcyBpZiBwcm92aWRlZCAtLS0tXG4gICAgaWYgKGlzLm51bGwoZWZmZWN0X21hdF9wcm9ncmVzc29yKSB8fCBpcy5udWxsKGVmZmVjdF9tYXRfY29udHJvbCkpIHtcbiAgICAgIHN0b3AoXCJCb3RoIGVmZmVjdF9tYXRfcHJvZ3Jlc3NvciBhbmQgZWZmZWN0X21hdF9jb250cm9sIG11c3QgYmUgcHJvdmlkZWQgaWYgdXNpbmcgZGVzaWduIG1hdHJpY2VzLlwiKVxuICAgIH1cblxuICAgIGVmZmVjdF9tYXRfcHJvZ3Jlc3NvciA8LSBhcy5tYXRyaXgoZWZmZWN0X21hdF9wcm9ncmVzc29yKVxuICAgIGVmZmVjdF9tYXRfY29udHJvbCAgICA8LSBhcy5tYXRyaXgoZWZmZWN0X21hdF9jb250cm9sKVxuXG4gICAgaWYgKG5jb2woZWZmZWN0X21hdF9wcm9ncmVzc29yKSAhPSB0aW1lX3BvaW50cyB8fFxuICAgICAgICBuY29sKGVmZmVjdF9tYXRfY29udHJvbCkgICAgIT0gdGltZV9wb2ludHMpIHtcbiAgICAgIHN0b3AoXCJlZmZlY3RfbWF0XyogbXVzdCBoYXZlIG5jb2wgPSB0aW1lX3BvaW50cy5cIilcbiAgICB9XG4gICAgaWYgKG5yb3coZWZmZWN0X21hdF9wcm9ncmVzc29yKSAhPSBuX2NlbGxfdHlwZXMgfHxcbiAgICAgICAgbnJvdyhlZmZlY3RfbWF0X2NvbnRyb2wpICAgICE9IG5fY2VsbF90eXBlcykge1xuICAgICAgc3RvcChcImVmZmVjdF9tYXRfKiBtdXN0IGhhdmUgbnJvdyA9IG5fbWFqb3JfY2VsbF90eXBlcyArIG5fbWlub3JfY2VsbF90eXBlcy5cIilcbiAgICB9XG5cbiAgICAjIEFsaWduIHJvd3MgdG8gY2VsbF90eXBlcyB2aWEgcm93bmFtZXMgaWYgcHJlc2VudDsgb3RoZXJ3aXNlIGFzc3VtZSBpbiBvcmRlclxuICAgIGlmICghaXMubnVsbChyb3duYW1lcyhlZmZlY3RfbWF0X3Byb2dyZXNzb3IpKSkge1xuICAgICAgaWYgKCFhbGwoc29ydChyb3duYW1lcyhlZmZlY3RfbWF0X3Byb2dyZXNzb3IpKSA9PSBzb3J0KGNlbGxfdHlwZXMpKSkge1xuICAgICAgICBzdG9wKFwiUm93IG5hbWVzIG9mIGVmZmVjdF9tYXRfcHJvZ3Jlc3NvciBtdXN0IG1hdGNoIGNlbGwgdHlwZXM6IFwiLFxuICAgICAgICAgICAgIHBhc3RlKGNlbGxfdHlwZXMsIGNvbGxhcHNlID0gXCIsIFwiKSlcbiAgICAgIH1cbiAgICAgIGVmZmVjdF9tYXRfcHJvZ3Jlc3NvciA8LSBlZmZlY3RfbWF0X3Byb2dyZXNzb3JbY2VsbF90eXBlcywgLCBkcm9wID0gRkFMU0VdXG4gICAgfSBlbHNlIHtcbiAgICAgIHJvd25hbWVzKGVmZmVjdF9tYXRfcHJvZ3Jlc3NvcikgPC0gY2VsbF90eXBlc1xuICAgIH1cblxuICAgIGlmICghaXMubnVsbChyb3duYW1lcyhlZmZlY3RfbWF0X2NvbnRyb2wpKSkge1xuICAgICAgaWYgKCFhbGwoc29ydChyb3duYW1lcyhlZmZlY3RfbWF0X2NvbnRyb2wpKSA9PSBzb3J0KGNlbGxfdHlwZXMpKSkge1xuICAgICAgICBzdG9wKFwiUm93IG5hbWVzIG9mIGVmZmVjdF9tYXRfY29udHJvbCBtdXN0IG1hdGNoIGNlbGwgdHlwZXM6IFwiLFxuICAgICAgICAgICAgIHBhc3RlKGNlbGxfdHlwZXMsIGNvbGxhcHNlID0gXCIsIFwiKSlcbiAgICAgIH1cbiAgICAgIGVmZmVjdF9tYXRfY29udHJvbCA8LSBlZmZlY3RfbWF0X2NvbnRyb2xbY2VsbF90eXBlcywgLCBkcm9wID0gRkFMU0VdXG4gICAgfSBlbHNlIHtcbiAgICAgIHJvd25hbWVzKGVmZmVjdF9tYXRfY29udHJvbCkgPC0gY2VsbF90eXBlc1xuICAgIH1cblxuICAgICMgTWFwIGVhY2ggcm93IGluIGNvdW50c19sb25nIHRvIChjZWxsX3R5cGUsIHZpc2l0KSAtPiBlZmZlY3RcbiAgICByb3dfaWR4IDwtIG1hdGNoKGNvdW50c19sb25nJGNlbGxfdHlwZSwgY2VsbF90eXBlcylcbiAgICBjb2xfaWR4IDwtIGNvdW50c19sb25nJHZpc2l0ICsgMUwgICMgdmlzaXRzIGFyZSAwLi50aW1lX3BvaW50cy0xXG5cbiAgICBpc19wcm9nIDwtIGNvdW50c19sb25nJGRpc2Vhc2UgPT0gMUxcblxuICAgIGVmZl92ZWNbaXNfcHJvZ10gIDwtIGVmZmVjdF9tYXRfcHJvZ3Jlc3NvcltjYmluZChyb3dfaWR4W2lzX3Byb2ddLCAgY29sX2lkeFtpc19wcm9nXSldXG4gICAgZWZmX3ZlY1shaXNfcHJvZ10gPC0gZWZmZWN0X21hdF9jb250cm9sWyAgIGNiaW5kKHJvd19pZHhbIWlzX3Byb2ddLCBjb2xfaWR4WyFpc19wcm9nXSldXG5cbiAgfSBlbHNlIHtcbiAgICAjIyAtLS0tIE9MRCBCRUhBVklPUjogdmVjdG9yIHZpc2l0IGVmZmVjdHMgKyBkaXJlY3Rpb25fYnlfY2x1c3RlciAtLS0tXG5cbiAgICAjIC0tLSBEZWZhdWx0IHZpc2l0IGVmZmVjdHMgKGxlbmd0aCA9IHRpbWVfcG9pbnRzKSAtLS0tXG4gICAgaWYgKGlzLm51bGwodmlzaXRfZWZmZWN0c19wcm9ncmVzc29yKSkge1xuICAgICAgaWYgKGludGVyYWN0aW9uX3R5cGUgPT0gXCJzcGVjaWZpY1wiKSB7XG4gICAgICAgIHZlIDwtIHJlcCgwLCB0aW1lX3BvaW50cylcbiAgICAgICAgaWYgKHRpbWVfcG9pbnRzID49IDIpIHZlWzJdIDwtIGZjX2ludGVyYWN0ICAjIGJ1bXAgVjEgb25seVxuICAgICAgICB2aXNpdF9lZmZlY3RzX3Byb2dyZXNzb3IgPC0gdmVcbiAgICAgIH0gZWxzZSBpZiAoaW50ZXJhY3Rpb25fdHlwZSA9PSBcImRpZmZlcmVudGlhbFwiKSB7XG4gICAgICAgIHZpc2l0X2VmZmVjdHNfcHJvZ3Jlc3NvciA8LSByZXAoZmNfaW50ZXJhY3QsIHRpbWVfcG9pbnRzKVxuICAgICAgfSBlbHNlIHsgIyBcIm9wcG9zaXRlXCI6IGFsdGVybmF0ZSArLy0gc3RhcnRpbmcgYXQgVjBcbiAgICAgICAgdmUgPC0gcmVwKDAsIHRpbWVfcG9pbnRzKVxuICAgICAgICB2ZVtzZXEoMSwgdGltZV9wb2ludHMsIGJ5ID0gMildIDwtICtmY19pbnRlcmFjdCAgIyBWMCwgVjIsIC4uLlxuICAgICAgICBpZiAodGltZV9wb2ludHMgPj0gMikgdmVbc2VxKDIsIHRpbWVfcG9pbnRzLCBieSA9IDIpXSA8LSAtZmNfaW50ZXJhY3QgIyBWMSwgVjMsIC4uLlxuICAgICAgICB2aXNpdF9lZmZlY3RzX3Byb2dyZXNzb3IgPC0gdmVcbiAgICAgIH1cbiAgICB9XG4gICAgaWYgKGlzLm51bGwodmlzaXRfZWZmZWN0c19jb250cm9sKSkge1xuICAgICAgdmlzaXRfZWZmZWN0c19jb250cm9sIDwtIHJlcCgwLCB0aW1lX3BvaW50cylcbiAgICB9XG4gICAgc3RvcGlmbm90KGxlbmd0aCh2aXNpdF9lZmZlY3RzX3Byb2dyZXNzb3IpID09IHRpbWVfcG9pbnRzLFxuICAgICAgICAgICAgICBsZW5ndGgodmlzaXRfZWZmZWN0c19jb250cm9sKSAgICA9PSB0aW1lX3BvaW50cylcblxuICAgICMgRGlyZWN0aW9uIHBlciBpbnRlcmFjdGluZyBjbHVzdGVyXG4gICAgaWYgKGlzLm51bGwoZGlyZWN0aW9uX2J5X2NsdXN0ZXIpKSB7XG4gICAgICBkaXJlY3Rpb25fYnlfY2x1c3RlciA8LSByZXAoMUwsIG1heCgxLCBsZW5ndGgoaW50ZXJhY3RfY2VsbF90eXBlcykpKVxuICAgIH1cbiAgICAjIHJlY3ljbGUgdG8gaW50ZXJhY3Rpbmcgc2V0IGxlbmd0aFxuICAgIGRpcmVjdGlvbl9ieV9jbHVzdGVyIDwtIHJlcChkaXJlY3Rpb25fYnlfY2x1c3RlciwgbGVuZ3RoLm91dCA9IGxlbmd0aChpbnRlcmFjdF9jZWxsX3R5cGVzKSlcbiAgICBuYW1lcyhkaXJlY3Rpb25fYnlfY2x1c3RlcikgPC0gaW50ZXJhY3RfY2VsbF90eXBlc1xuXG4gICAgIyBBcHBseSBlZmZlY3RzIG9ubHkgdG8gaW50ZXJhY3RpbmcgY2VsbCB0eXBlc1xuICAgIGlzX2ludGVyYWN0aW5nIDwtIGNvdW50c19sb25nJGNlbGxfdHlwZSAlaW4lIGludGVyYWN0X2NlbGxfdHlwZXNcbiAgICBpZiAoYW55KGlzX2ludGVyYWN0aW5nKSkge1xuICAgICAgIyBtYXAgZGlyZWN0aW9uIHBlciBjZWxsIHR5cGVcbiAgICAgIGRpcl9tYXAgPC0gc2V0TmFtZXMoZGlyZWN0aW9uX2J5X2NsdXN0ZXIsIG5tID0gbmFtZXMoZGlyZWN0aW9uX2J5X2NsdXN0ZXIpKVxuICAgICAgZGlyX2N0ICA8LSB1bm5hbWUoZGlyX21hcFtjb3VudHNfbG9uZyRjZWxsX3R5cGVbaXNfaW50ZXJhY3RpbmddXSlcbiAgICAgICMgdmlzaXQgZWZmZWN0cyBieSBncm91cFxuICAgICAgdl9pZHggICA8LSBjb3VudHNfbG9uZyR2aXNpdFtpc19pbnRlcmFjdGluZ10gKyAxTFxuICAgICAgaXNfcHJvZyA8LSBjb3VudHNfbG9uZyRkaXNlYXNlW2lzX2ludGVyYWN0aW5nXSA9PSAxTFxuICAgICAgdmUgICAgICA8LSBpZmVsc2UoaXNfcHJvZywgdmlzaXRfZWZmZWN0c19wcm9ncmVzc29yW3ZfaWR4XSwgdmlzaXRfZWZmZWN0c19jb250cm9sW3ZfaWR4XSlcbiAgICAgIGVmZl92ZWNbaXNfaW50ZXJhY3RpbmddIDwtIGRpcl9jdCAqIHZlXG4gICAgfVxuICB9XG5cbiAgIyAtLS0gQXBwbHkgZWZmZWN0cyB0byBjb3VudHMgLS0tLVxuICBjb3VudHNfbG9uZyRhZGpfY291bnQgPC0gcG1heCgwTCwgcm91bmQoY291bnRzX2xvbmckY291bnQgKiAoMSArIGVmZl92ZWMpKSlcblxuICAjIC0tLSBFeHBhbmQgdG8gcGVyLWNlbGwgcm93cyBhbmQgYXR0YWNoIG1ldGFkYXRhIChJTkNMVURJTkcgaW50ZXJhY3RfdGVybSkgLS0tLVxuICByZXBfZWFjaCA8LSBmdW5jdGlvbih4LCB0aW1lcykgaWYgKGxlbmd0aCh4KSA9PSAwKSB4IGVsc2UgcmVwKHgsIHRpbWVzID0gdGltZXMpXG4gIGV4cGFuZGVkIDwtIGRhdGEuZnJhbWUoXG4gICAgc2FtcGxlX2lkID0gcmVwX2VhY2goY291bnRzX2xvbmckc2FtcGxlX2lkLCBjb3VudHNfbG9uZyRhZGpfY291bnQpLFxuICAgIGNlbGxfdHlwZSA9IHJlcF9lYWNoKGNvdW50c19sb25nJGNlbGxfdHlwZSwgY291bnRzX2xvbmckYWRqX2NvdW50KSxcbiAgICBzdHJpbmdzQXNGYWN0b3JzID0gRkFMU0VcbiAgKVxuXG4gICMgTWVyZ2UgbWV0YWRhdGE7IGtlZXAgaW50ZXJhY3Rpb24gKyBpbnRlcmFjdF90ZXJtXG4gIGtlZXBfY29scyA8LSBjKFwic2FtcGxlX2lkXCIsXCJzdWJqZWN0X2lkXCIsXCJ2aXNpdFwiLFwic2V4XCIsXCJkaXNlYXNlXCIsXCJhZ2VcIixcImJtaVwiLFwiYmF0Y2hcIixcbiAgICAgICAgICAgICAgICAgXCJpbnRlcmFjdGlvblwiLFwiaW50ZXJhY3RfdGVybVwiKVxuICBkdW1teV9kYXRhIDwtIG1lcmdlKGV4cGFuZGVkLCBtZXRhWywga2VlcF9jb2xzXSwgYnkgPSBcInNhbXBsZV9pZFwiLCBzb3J0ID0gRkFMU0UpXG5cbiAgIyBTaHVmZmxlIHJvd3MgZm9yIHJlYWxpc21cbiAgaWYgKG5yb3coZHVtbXlfZGF0YSkgPiAxKSB7XG4gICAgZHVtbXlfZGF0YSA8LSBkdW1teV9kYXRhW3NhbXBsZS5pbnQobnJvdyhkdW1teV9kYXRhKSksICwgZHJvcCA9IEZBTFNFXVxuICAgIHJvd25hbWVzKGR1bW15X2RhdGEpIDwtIE5VTExcbiAgfVxuXG4gIGR1bW15X2RhdGFcbn1cblxuYGBgIn0= -->

```r
#' Simulate longitudinal single-cell data with interacting cell types
#'
#' @param n_cells Baseline cells per major cell type per sample.
#' @param sd_celltypes Relative SD for counts.
#' @param n_major_cell_types Number of major cell types.
#' @param n_minor_cell_types Number of minor cell types.
#' @param relative_abundance Baseline minor:major ratio.
#' @param n_major_interact_celltypes How many major types are "interacting".
#' @param n_minor_interact_celltypes How many minor types are "interacting".
#' @param n_individuals Number of subjects.
#' @param n_batchs Number of batches.
#' @param interaction_feature Name of the time/visit variable (e.g. "visit").
#' @param time_points Number of longitudinal time points (>= 2).
#' @param test_var Name of the disease/exposure variable (e.g. "disease").
#' @param prop_disease Proportion of diseased subjects.
#' @param fc_interact Effect magnitude used by the default effect patterns.
#' @param interaction_type One of "specific", "differential", or "opposite"
#'   for the default longitudinal pattern (when design matrices are not given).
#' @param seed Random seed.
#' @param visit_effects_progressor Optional vector of length `time_points`
#'   for default progressor visit effects (ignored if design matrices given).
#' @param visit_effects_control Optional vector of length `time_points`
#'   for default control visit effects (ignored if design matrices given).
#' @param direction_by_cluster Optional vector of +1 / -1 per interacting
#'   cell type (ignored if design matrices given).
#' @param effect_mat_progressor Optional numeric matrix of dimension
#'   `n_cell_types x time_points` giving effect sizes for progressors.
#'   Rows correspond to cell types (by name or order), columns to visits
#'   (0..time_points-1). If provided (along with `effect_mat_control`),
#'   overrides `visit_effects_*` and `direction_by_cluster`.
#' @param effect_mat_control Optional numeric matrix of same dimension as
#'   `effect_mat_progressor` giving effects for controls.
#'
#' @return A data.frame of simulated single-cell metadata.
#' @export
generate_dummy_data <- function(
  n_cells = 3000,                  # baseline cells per major cell type per sample
  sd_celltypes = 0.10,             # relative sd for counts
  n_major_cell_types = 7,
  n_minor_cell_types = 3,
  relative_abundance = 0.10,       # minor vs major baseline ratio
  n_major_interact_celltypes = 1,  # how many majors are "interacting"
  n_minor_interact_celltypes = 1,  # how many minors are "interacting"
  n_individuals = 30,
  n_batchs = 4,

  interaction_feature = "visit",   # kept for labeling
  time_points = 4,                 # >= 2; works for 3+
  test_var = "disease",
  prop_disease = 0.50,

  fc_interact = 0.10,              # effect magnitude used by defaults below
  interaction_type = c("specific","differential","opposite"),
  seed = 1234,

  visit_effects_progressor = NULL, # multiplicative change: +0.1 means +10% vs baseline
  visit_effects_control    = NULL, # default 0s
  direction_by_cluster     = NULL, # +1/-1 for interacting clusters, recycled as needed

  effect_mat_progressor    = NULL, # NEW: cell_type x time_points matrix for progressors
  effect_mat_control       = NULL  # NEW: same for controls
) {
  set.seed(seed)
  interaction_type <- match.arg(interaction_type)

  # --- Basic setup ----
  n_cell_types <- n_major_cell_types + n_minor_cell_types
  stopifnot(time_points >= 2, n_cell_types >= 1)

  cell_types <- LETTERS[seq_len(n_cell_types)]
  major_idx  <- seq_len(n_major_cell_types)
  minor_idx  <- if (n_minor_cell_types > 0) (n_major_cell_types + seq_len(n_minor_cell_types)) else integer(0)

  # which cell types are "interacting" (first some majors, last some minors)
  interact_idx <- c(
    head(major_idx, n_major_interact_celltypes),
    tail(seq_len(n_cell_types), n_minor_interact_celltypes)
  )
  interact_idx <- intersect(interact_idx, seq_len(n_cell_types))  # guard rails
  interact_cell_types <- cell_types[interact_idx]

  # --- Subjects and visits ----
  subject_id <- paste0("SUB_", seq_len(n_individuals))
  disease_vec <- c(rep(1L, round(n_individuals * prop_disease)),
                   rep(0L, n_individuals - round(n_individuals * prop_disease)))
  disease_vec <- sample(disease_vec, n_individuals)

  sex_vec <- sample(c(0L, 1L), n_individuals, replace = TRUE)          # 0/1
  age_vec <- sample(18:60, n_individuals, replace = TRUE)
  bmi_vec <- sample(15:35, n_individuals, replace = TRUE)
  batch_vec <- rep(seq_len(n_batchs), length.out = n_individuals)

  subjects <- data.frame(
    subject_id = subject_id,
    sex   = sex_vec,
    disease = disease_vec,           # canonical exposure storage
    age   = age_vec,
    bmi   = bmi_vec,
    batch = factor(batch_vec),
    stringsAsFactors = FALSE
  )

  visits <- data.frame(
    subject_id = rep(subject_id, each = time_points),
    visit      = rep(0:(time_points - 1), times = n_individuals),
    stringsAsFactors = FALSE
  )
  visits$sample_id <- paste0(visits$subject_id, "_V", visits$visit)

  meta <- merge(visits, subjects, by = "subject_id", sort = FALSE)

  # Make sure a column named `test_var` exists (even if test_var != "disease")
  if (!identical(test_var, "disease")) {
    meta[[test_var]] <- meta[["disease"]]
  }

  # Label and INTERACTION TERM (persist to output)
  meta$interaction    <- paste0(interaction_feature, ":", test_var)
  meta$interact_term  <- as.integer(meta[[interaction_feature]]) * as.integer(meta[[test_var]])

  # --- Baseline counts per (sample, cell_type) ----
  one_sample_counts <- function() {
    major_counts <- round(runif(n_major_cell_types,
                                min = n_cells * (1 - sd_celltypes),
                                max = n_cells * (1 + sd_celltypes)))
    minor_counts <- if (n_minor_cell_types > 0) {
      round(runif(n_minor_cell_types,
                  min = n_cells * relative_abundance * (1 - sd_celltypes),
                  max = n_cells * relative_abundance * (1 + sd_celltypes)))
    } else integer(0)
    c(major_counts, minor_counts)
  }

  counts_list <- replicate(nrow(meta), one_sample_counts(), simplify = FALSE)
  counts_df <- do.call(rbind, counts_list)
  colnames(counts_df) <- cell_types

  counts_long <- reshape(
    data.frame(sample_id = meta$sample_id, counts_df, check.names = FALSE),
    varying = cell_types, v.names = "count", timevar = "cell_type",
    times = cell_types, direction = "long"
  )
  rownames(counts_long) <- NULL

  # Merge disease/visit so we can apply effects
  counts_long <- merge(counts_long,
                       meta[, c("sample_id", "visit", "disease")],
                       by = "sample_id", sort = FALSE)

  # --- Build effect structure ----
  eff_vec <- numeric(nrow(counts_long))

  if (!is.null(effect_mat_progressor) || !is.null(effect_mat_control)) {
    ## ---- NEW: use design matrices if provided ----
    if (is.null(effect_mat_progressor) || is.null(effect_mat_control)) {
      stop("Both effect_mat_progressor and effect_mat_control must be provided if using design matrices.")
    }

    effect_mat_progressor <- as.matrix(effect_mat_progressor)
    effect_mat_control    <- as.matrix(effect_mat_control)

    if (ncol(effect_mat_progressor) != time_points ||
        ncol(effect_mat_control)    != time_points) {
      stop("effect_mat_* must have ncol = time_points.")
    }
    if (nrow(effect_mat_progressor) != n_cell_types ||
        nrow(effect_mat_control)    != n_cell_types) {
      stop("effect_mat_* must have nrow = n_major_cell_types + n_minor_cell_types.")
    }

    # Align rows to cell_types via rownames if present; otherwise assume in order
    if (!is.null(rownames(effect_mat_progressor))) {
      if (!all(sort(rownames(effect_mat_progressor)) == sort(cell_types))) {
        stop("Row names of effect_mat_progressor must match cell types: ",
             paste(cell_types, collapse = ", "))
      }
      effect_mat_progressor <- effect_mat_progressor[cell_types, , drop = FALSE]
    } else {
      rownames(effect_mat_progressor) <- cell_types
    }

    if (!is.null(rownames(effect_mat_control))) {
      if (!all(sort(rownames(effect_mat_control)) == sort(cell_types))) {
        stop("Row names of effect_mat_control must match cell types: ",
             paste(cell_types, collapse = ", "))
      }
      effect_mat_control <- effect_mat_control[cell_types, , drop = FALSE]
    } else {
      rownames(effect_mat_control) <- cell_types
    }

    # Map each row in counts_long to (cell_type, visit) -> effect
    row_idx <- match(counts_long$cell_type, cell_types)
    col_idx <- counts_long$visit + 1L  # visits are 0..time_points-1

    is_prog <- counts_long$disease == 1L

    eff_vec[is_prog]  <- effect_mat_progressor[cbind(row_idx[is_prog],  col_idx[is_prog])]
    eff_vec[!is_prog] <- effect_mat_control[   cbind(row_idx[!is_prog], col_idx[!is_prog])]

  } else {
    ## ---- OLD BEHAVIOR: vector visit effects + direction_by_cluster ----

    # --- Default visit effects (length = time_points) ----
    if (is.null(visit_effects_progressor)) {
      if (interaction_type == "specific") {
        ve <- rep(0, time_points)
        if (time_points >= 2) ve[2] <- fc_interact  # bump V1 only
        visit_effects_progressor <- ve
      } else if (interaction_type == "differential") {
        visit_effects_progressor <- rep(fc_interact, time_points)
      } else { # "opposite": alternate +/- starting at V0
        ve <- rep(0, time_points)
        ve[seq(1, time_points, by = 2)] <- +fc_interact  # V0, V2, ...
        if (time_points >= 2) ve[seq(2, time_points, by = 2)] <- -fc_interact # V1, V3, ...
        visit_effects_progressor <- ve
      }
    }
    if (is.null(visit_effects_control)) {
      visit_effects_control <- rep(0, time_points)
    }
    stopifnot(length(visit_effects_progressor) == time_points,
              length(visit_effects_control)    == time_points)

    # Direction per interacting cluster
    if (is.null(direction_by_cluster)) {
      direction_by_cluster <- rep(1L, max(1, length(interact_cell_types)))
    }
    # recycle to interacting set length
    direction_by_cluster <- rep(direction_by_cluster, length.out = length(interact_cell_types))
    names(direction_by_cluster) <- interact_cell_types

    # Apply effects only to interacting cell types
    is_interacting <- counts_long$cell_type %in% interact_cell_types
    if (any(is_interacting)) {
      # map direction per cell type
      dir_map <- setNames(direction_by_cluster, nm = names(direction_by_cluster))
      dir_ct  <- unname(dir_map[counts_long$cell_type[is_interacting]])
      # visit effects by group
      v_idx   <- counts_long$visit[is_interacting] + 1L
      is_prog <- counts_long$disease[is_interacting] == 1L
      ve      <- ifelse(is_prog, visit_effects_progressor[v_idx], visit_effects_control[v_idx])
      eff_vec[is_interacting] <- dir_ct * ve
    }
  }

  # --- Apply effects to counts ----
  counts_long$adj_count <- pmax(0L, round(counts_long$count * (1 + eff_vec)))

  # --- Expand to per-cell rows and attach metadata (INCLUDING interact_term) ----
  rep_each <- function(x, times) if (length(x) == 0) x else rep(x, times = times)
  expanded <- data.frame(
    sample_id = rep_each(counts_long$sample_id, counts_long$adj_count),
    cell_type = rep_each(counts_long$cell_type, counts_long$adj_count),
    stringsAsFactors = FALSE
  )

  # Merge metadata; keep interaction + interact_term
  keep_cols <- c("sample_id","subject_id","visit","sex","disease","age","bmi","batch",
                 "interaction","interact_term")
  dummy_data <- merge(expanded, meta[, keep_cols], by = "sample_id", sort = FALSE)

  # Shuffle rows for realism
  if (nrow(dummy_data) > 1) {
    dummy_data <- dummy_data[sample.int(nrow(dummy_data)), , drop = FALSE]
    rownames(dummy_data) <- NULL
  }

  dummy_data
}

```

<!-- rnb-source-end -->


<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->






<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuY2VsbF90eXBlcyA8LSBMRVRURVJTWzE6MTBdXG50aW1lX3BvaW50cyA8LSAzXG5cbiMgcHJvZ3Jlc3NvciBlZmZlY3RzOiByb3dzID0gY2VsbCB0eXBlcywgY29scyA9IHZpc2l0c1xuRV9wcm9nIDwtIG1hdHJpeCgwLCBucm93ID0gbGVuZ3RoKGNlbGxfdHlwZXMpLCBuY29sID0gdGltZV9wb2ludHMsXG4gICAgICAgICAgICAgICAgIGRpbW5hbWVzID0gbGlzdChjZWxsX3R5cGVzLCBwYXN0ZTAoXCJWXCIsIDA6KHRpbWVfcG9pbnRzIC0gMSkpKSlcblxuIyBleGFtcGxlOiBlbnJpY2ggQSAmIEkgYXQgVjEgYW5kIFYyLCBkZXBsZXRlIEIgJiBKIGF0IFYxXG5FX3Byb2dbXCJBXCIsIF0gPC0gYygwLjAsIDAuNiwgMC4zKVxuRV9wcm9nW1wiSVwiLCBdIDwtIGMoMC4wLCAwLjYsIDAuMylcbkVfcHJvZ1tcIkJcIiwgXSA8LSBjKDAuMCwgLTAuNiwgMC4wKVxuRV9wcm9nW1wiSlwiLCBdIDwtIGMoMC4wLCAtMC42LCAwLjApXG5cbiMgY29udHJvbCBlZmZlY3RzIChhbGwgemVybyBpbiB0aGlzIGV4YW1wbGUpXG5FX2N0cmwgPC0gbWF0cml4KDAsIG5yb3cgPSBsZW5ndGgoY2VsbF90eXBlcyksIG5jb2wgPSB0aW1lX3BvaW50cyxcbiAgICAgICAgICAgICAgICAgZGltbmFtZXMgPSBsaXN0KGNlbGxfdHlwZXMsIHBhc3RlMChcIlZcIiwgMDoodGltZV9wb2ludHMgLSAxKSkpKVxuXG5kMyA8LSBnZW5lcmF0ZV9kdW1teV9kYXRhKFxuICBuX2NlbGxzID0gMTUwLFxuICBzZF9jZWxsdHlwZXMgPSAwLjEsXG4gIG5fbWFqb3JfY2VsbF90eXBlcyA9IDcsXG4gbl9taW5vcl9jZWxsX3R5cGVzID0gMyxcbiAgcmVsYXRpdmVfYWJ1bmRhbmNlID0gMC40LFxuICBuX21ham9yX2ludGVyYWN0X2NlbGx0eXBlcyA9IDIsXG4gIG5fbWlub3JfaW50ZXJhY3RfY2VsbHR5cGVzID0gMixcbiAgbl9pbmRpdmlkdWFscyA9IDMwLFxuICBuX2JhdGNocyA9IDQsXG4gIGludGVyYWN0aW9uX2ZlYXR1cmUgPSBcInZpc2l0XCIsXG4gIHRpbWVfcG9pbnRzID0gMyxcbiAgdGVzdF92YXIgPSBcImRpc2Vhc2VcIixcbiAgcHJvcF9kaXNlYXNlID0gMC41LFxuICBmY19pbnRlcmFjdCA9IDAuNixcbiAgaW50ZXJhY3Rpb25fdHlwZSA9IFwic3BlY2lmaWNcIixcbiAgc2VlZCA9IDIwMjAsXG4gIGVmZmVjdF9tYXRfcHJvZ3Jlc3NvciA9IEVfcHJvZyxcbiAgZWZmZWN0X21hdF9jb250cm9sICAgID0gRV9jdHJsXG4pXG5cbmBgYCJ9 -->

```r
cell_types <- LETTERS[1:10]
time_points <- 3

# progressor effects: rows = cell types, cols = visits
E_prog <- matrix(0, nrow = length(cell_types), ncol = time_points,
                 dimnames = list(cell_types, paste0("V", 0:(time_points - 1))))

# example: enrich A & I at V1 and V2, deplete B & J at V1
E_prog["A", ] <- c(0.0, 0.6, 0.3)
E_prog["I", ] <- c(0.0, 0.6, 0.3)
E_prog["B", ] <- c(0.0, -0.6, 0.0)
E_prog["J", ] <- c(0.0, -0.6, 0.0)

# control effects (all zero in this example)
E_ctrl <- matrix(0, nrow = length(cell_types), ncol = time_points,
                 dimnames = list(cell_types, paste0("V", 0:(time_points - 1))))

d3 <- generate_dummy_data(
  n_cells = 150,
  sd_celltypes = 0.1,
  n_major_cell_types = 7,
 n_minor_cell_types = 3,
  relative_abundance = 0.4,
  n_major_interact_celltypes = 2,
  n_minor_interact_celltypes = 2,
  n_individuals = 30,
  n_batchs = 4,
  interaction_feature = "visit",
  time_points = 3,
  test_var = "disease",
  prop_disease = 0.5,
  fc_interact = 0.6,
  interaction_type = "specific",
  seed = 2020,
  effect_mat_progressor = E_prog,
  effect_mat_control    = E_ctrl
)

```

<!-- rnb-source-end -->


<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYSA9IHNjTEFTRVI6OnNjTEFTRVIoZDMpXG5cbnNjTEFTRVI6OnBsb3RfY2VsbHR5cGVfcHJvcG9ydGlvbnMoYSxoaWdobGlnaHRfY2VsbF90eXBlcyA9IGMoXCJBXCIsXCJCXCIpKVxuYGBgIn0= -->

```r
a = scLASER::scLASER(d3)

scLASER::plot_celltype_proportions(a,highlight_cell_types = c("A","B"))
```

<!-- rnb-source-end -->


<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->




### PLSDA function modification station


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIycgUExTLURBIHN0ZXAgZm9yIHNjTEFTRVIgb2JqZWN0cyB1c2luZyBtaXhPbWljc1xuIydcbiMnIFdyYXBwZXIgYXJvdW5kIFttaXhPbWljczo6cGxzZGEoKV0gdGhhdCBwdWxscyB0aGUgZmVhdHVyZSBtYXRyaXggZnJvbSBhXG4jJyBgc2NMQVNFUmAgb2JqZWN0LCBjb25zdHJ1Y3RzIGNsYXNzIGxhYmVscyBmcm9tIG1ldGFkYXRhLCBvcHRpb25hbGx5IHBlcmZvcm1zXG4jJyBtdWx0aWxldmVsIGRlY29tcG9zaXRpb24sIGZpdHMgUExTLURBLCBhbmQgc2F2ZXMgdGhlIGxhdGVudCB2YXJpYWJsZXNcbiMnIChjb21wb25lbnQgc2NvcmVzKSBpbnRvIGBvYmplY3RAcGxzZGFfTFZgLlxuIydcbiMnIE9wdGlvbmFsbHksIGEgc3ViamVjdC13aXNlIGxvbmdpdHVkaW5hbCBjcm9zcy12YWxpZGF0aW9uIChDVikgaXMgcGVyZm9ybWVkLFxuIycgd2hlcmUgZm9sZHMgYXJlIGNyZWF0ZWQgYXQgdGhlIHN1YmplY3QgbGV2ZWwgKGFsbCB2aXNpdHMgZm9yIGEgc3ViamVjdCBhcmVcbiMnIGtlcHQgaW4gdGhlIHNhbWUgZm9sZCkgYW5kIHRoZSBCYWxhbmNlZCBFcnJvciBSYXRlIChCRVIpIGlzIGNvbXB1dGVkIGZvclxuIycgMS4uYGN2X21heF9uY29tcGAgY29tcG9uZW50cy4gQ1YgcmVzdWx0cyBhcmUgYXR0YWNoZWQgYXNcbiMnIGBhdHRyKG9iamVjdCwgXCJwbHNkYV9jdlwiKWAuXG4jJ1xuIycgQHBhcmFtIG9iamVjdCAgICAgICAgIEEgYHNjTEFTRVJgIG9iamVjdC5cbiMnIEBwYXJhbSByZXNwb25zZV92YXIgICBDaGFyYWN0ZXIuIENvbHVtbiBpbiBgb2JqZWN0QG1ldGFkYXRhYCB1c2VkIGFzIGNsYXNzXG4jJyAgIGxhYmVsIChkZWZhdWx0IGBcImNlbGxfdHlwZVwiYCkuXG4jJyBAcGFyYW0gbXVsdGlsZXZlbF92YXIgT3B0aW9uYWwgY2hhcmFjdGVyLiBDb2x1bW4gaW4gYG9iamVjdEBtZXRhZGF0YWBcbiMnICAgaWRlbnRpZnlpbmcgcmVwZWF0ZWQtbWVhc3VyZXMgdW5pdHMgKGUuZy4gYFwic2FtcGxlX2lkXCJgKS4gSWYgbm90IGBOVUxMYCxcbiMnICAgbXVsdGlsZXZlbCBQTFMtREEgaXMgcGVyZm9ybWVkIHZpYSBbbWl4T21pY3M6OndpdGhpblZhcmlhdGlvbigpXSBwcmlvclxuIycgICB0byBmaXR0aW5nLlxuIycgQHBhcmFtIG5jb21wICAgICAgICAgIEludGVnZXIuIE51bWJlciBvZiBsYXRlbnQgY29tcG9uZW50cyB0byBjb21wdXRlIGZvclxuIycgICB0aGUgbWFpbiBQTFMtREEgZml0IGFuZCBmb3IgQ1YgKHVwcGVyIGJvdW5kKS5cbiMnIEBwYXJhbSBjdl9sb25naXR1ZGluYWwgTG9naWNhbC4gSWYgYFRSVUVgLCBwZXJmb3JtIHN1YmplY3Qtd2lzZSBsb25naXR1ZGluYWxcbiMnICAgQ1YgYXMgZGVzY3JpYmVkIGFib3ZlLiBEZWZhdWx0IGlzIGBGQUxTRWAuXG4jJyBAcGFyYW0gY3Zfc3ViamVjdF92YXIgQ2hhcmFjdGVyLiBDb2x1bW4gaW4gYG9iamVjdEBtZXRhZGF0YWAgdGhhdCBkZWZpbmVzXG4jJyAgIHRoZSBzdWJqZWN0IElEIGZvciBDViBmb2xkcyAoZGVmYXVsdCBgXCJzdWJqZWN0X2lkXCJgKS5cbiMnIEBwYXJhbSBjdl9LICAgICAgICAgICBJbnRlZ2VyLiBOdW1iZXIgb2YgZm9sZHMgZm9yIGxvbmdpdHVkaW5hbCBDVlxuIycgICAoZGVmYXVsdCBgNWApLlxuIycgQHBhcmFtIGN2X21heF9uY29tcCAgIEludGVnZXIuIE1heGltdW0gbnVtYmVyIG9mIGNvbXBvbmVudHMgdG8gZXZhbHVhdGUgaW5cbiMnICAgbG9uZ2l0dWRpbmFsIENWLiBEZWZhdWx0IGlzIGBuY29tcGAuIFZhbHVlcyBncmVhdGVyIHRoYW4gYG5jb21wYCBhcmVcbiMnICAgdHJ1bmNhdGVkIHRvIGBuY29tcGAuXG4jJyBAcGFyYW0gY3Zfc2VlZCAgICAgICAgT3B0aW9uYWwgaW50ZWdlciBzZWVkIGZvciByZXByb2R1Y2libGUgZm9sZCBjcmVhdGlvbi5cbiMnICAgSWYgYE5VTExgLCB0aGUgY3VycmVudCBSTkcgc3RhdGUgaXMgdXNlZC5cbiMnIEBwYXJhbSAuLi4gICAgICAgICAgICBBZGRpdGlvbmFsIGFyZ3VtZW50cyBwYXNzZWQgdG8gW21peE9taWNzOjpwbHNkYSgpXS5cbiMnXG4jJyBAcmV0dXJuIEEgYHNjTEFTRVJgIG9iamVjdCB3aXRoOlxuIycgICBcXGl0ZW1pemV7XG4jJyAgICAgXFxpdGVtIGBAcGxzZGFfTFZgIHBvcHVsYXRlZCB3aXRoIHRoZSBQTFMtREEgbGF0ZW50IHZhcmlhYmxlc1xuIycgICAgICAgKHNjb3JlczsgYG5fY2VsbHMgeCBuY29tcGApLlxuIycgICAgIFxcaXRlbSBJZiBgY3ZfbG9uZ2l0dWRpbmFsID0gVFJVRWAsIGFuIGF0dHJpYnV0ZSBgcGxzZGFfY3ZgIGF0dGFjaGVkOlxuIycgICAgICAgYGF0dHIob2JqZWN0LCBcInBsc2RhX2N2XCIpYCwgYSBsaXN0IHdpdGggZWxlbWVudHNcbiMnICAgICAgIGBiZXJgIChudW1lcmljIEJFUiBwZXIgY29tcG9uZW50KSxcbiMnICAgICAgIGBvcHRfbmNvbXBgIChvcHRpbWFsIG51bWJlciBvZiBjb21wb25lbnRzKSxcbiMnICAgICAgIGFuZCBgbmNvbXBfc2VxYCAoc2VxdWVuY2Ugb2YgY29tcG9uZW50cyBldmFsdWF0ZWQpLlxuIycgICB9XG4jJ1xuIycgQGV4cG9ydFxuIycgQGltcG9ydEZyb20gbWl4T21pY3MgcGxzZGEgd2l0aGluVmFyaWF0aW9uXG5taXhPbWljc19wbHNkYSA8LSBmdW5jdGlvbihvYmplY3QsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICByZXNwb25zZV92YXIgICAgPSBcImNlbGxfdHlwZVwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgbXVsdGlsZXZlbF92YXIgID0gTlVMTCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgIG5jb21wICAgICAgICAgICA9IDIsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICBjdl9sb25naXR1ZGluYWwgPSBGQUxTRSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgIGN2X3N1YmplY3RfdmFyICA9IFwic3ViamVjdF9pZFwiLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgY3ZfSyAgICAgICAgICAgID0gNSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgIGN2X21heF9uY29tcCAgICA9IG5jb21wLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgY3Zfc2VlZCAgICAgICAgID0gTlVMTCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgIC4uLikge1xuXG4gIHN0b3BpZm5vdChpbmhlcml0cyhvYmplY3QsIFwic2NMQVNFUlwiKSlcblxuICAjIyAxLiBGZWF0dXJlIG1hdHJpeCBmcm9tIE5BTV9tYXRyaXggLS0tLVxuICBYIDwtIG9iamVjdEBOQU1fbWF0cml4XG4gIGlmIChpcy5udWxsKFgpKSB7XG4gICAgc3RvcChcbiAgICAgIFwib2JqZWN0QE5BTV9tYXRyaXggaXMgTlVMTC4gXCIsXG4gICAgICBcIlBsZWFzZSBjb21wdXRlIE5BTV9tYXRyaXggYmVmb3JlIGNhbGxpbmcgbWl4T21pY3NfcGxzZGEoKS5cIlxuICAgIClcbiAgfVxuICBYIDwtIGFzLm1hdHJpeChYKVxuXG4gICMjIDIuIFJlc3BvbnNlIGZyb20gbWV0YWRhdGEgLS0tLVxuICBtZXRhIDwtIG9iamVjdEBtZXRhZGF0YVxuICBpZiAoIXJlc3BvbnNlX3ZhciAlaW4lIGNvbG5hbWVzKG1ldGEpKSB7XG4gICAgc3RvcChcInJlc3BvbnNlX3ZhciAnXCIsIHJlc3BvbnNlX3ZhciwgXCInIG5vdCBmb3VuZCBpbiBvYmplY3RAbWV0YWRhdGEuXCIpXG4gIH1cbiAgWSA8LSBmYWN0b3IobWV0YVtbcmVzcG9uc2VfdmFyXV0pXG5cbiAgaWYgKG5yb3coWCkgIT0gbnJvdyhtZXRhKSkge1xuICAgIHN0b3AoXG4gICAgICBcIk51bWJlciBvZiByb3dzIGluIE5BTV9tYXRyaXggKFwiLCBucm93KFgpLFxuICAgICAgXCIpIGRvZXMgbm90IG1hdGNoIG51bWJlciBvZiByb3dzIGluIG1ldGFkYXRhIChcIiwgbnJvdyhtZXRhKSwgXCIpLlwiXG4gICAgKVxuICB9XG5cbiAgIyMgMy4gT3B0aW9uYWwgbXVsdGlsZXZlbCBkZXNpZ24gKHdpdGhpblZhcmlhdGlvbikgLS0tLVxuICBpZiAoIWlzLm51bGwobXVsdGlsZXZlbF92YXIpKSB7XG4gICAgaWYgKCFtdWx0aWxldmVsX3ZhciAlaW4lIGNvbG5hbWVzKG1ldGEpKSB7XG4gICAgICBzdG9wKFwibXVsdGlsZXZlbF92YXIgJ1wiLCBtdWx0aWxldmVsX3ZhcixcbiAgICAgICAgICAgXCInIG5vdCBmb3VuZCBpbiBvYmplY3RAbWV0YWRhdGEuXCIpXG4gICAgfVxuICAgIGRlc2lnbiA8LSBkYXRhLmZyYW1lKHNhbXBsZSA9IG1ldGFbW211bHRpbGV2ZWxfdmFyXV0pXG4gICAgWCA8LSB3aXRoaW5WYXJpYXRpb24oWCwgZGVzaWduID0gZGVzaWduKVxuICB9XG5cbiAgIyMgNC4gRml0IG1haW4gUExTLURBIG1vZGVsIC0tLS1cbiAgIyBVc2UgbWl4T21pY3M6OnBsc2RhOyBuY29tcCBpcyB0aGUgbWF4IG51bWJlciBvZiBjb21wb25lbnRzIHdlIGtlZXBcbiAjIGZpdCA8LSBwbHNkYShYLCBZLCBuY29tcCA9IG5jb21wLCAuLi4pXG5cbiAgIyBMYXRlbnQgdmFyaWFibGVzIChzY29yZXMpOiBuX2NlbGxzIHggbmNvbXBcbiAgI0xWIDwtIGZpdCR2YXJpYXRlcyRYXG4gICNvYmplY3RAcGxzZGFfTFYgPC0gYXMubWF0cml4KExWKVxuXG4gICMjIDUuIE9wdGlvbmFsIHN1YmplY3Qtd2lzZSBsb25naXR1ZGluYWwgQ1YgLS0tLVxuICBpZiAoaXNUUlVFKGN2X2xvbmdpdHVkaW5hbCkpIHtcbiAgICBpZiAoIWN2X3N1YmplY3RfdmFyICVpbiUgY29sbmFtZXMobWV0YSkpIHtcbiAgICAgIHN0b3AoXCJjdl9zdWJqZWN0X3ZhciAnXCIsIGN2X3N1YmplY3RfdmFyLFxuICAgICAgICAgICBcIicgbm90IGZvdW5kIGluIG9iamVjdEBtZXRhZGF0YS5cIilcbiAgICB9XG5cbiAgICAjIFByZXBhcmUgc3ViamVjdC1sZXZlbCB0YWJsZSAodW5pcXVlIChzdWJqZWN0LCBjbGFzcykgY29tYmluYXRpb25zKVxuICAgIHN1YmpfZGYgPC0gbWV0YVssIGMoY3Zfc3ViamVjdF92YXIsIHJlc3BvbnNlX3ZhcildXG4gICAgc3Vial9kZiA8LSBzdWJqX2RmWyFkdXBsaWNhdGVkKHN1YmpfZGYpLCAsIGRyb3AgPSBGQUxTRV1cblxuICAgICMgUmFuZG9taXplIG9yZGVyIG9mIHN1YmplY3RzIGZvciBmb2xkIGFzc2lnbm1lbnRcbiAgICBpZiAoIWlzLm51bGwoY3Zfc2VlZCkpIHtcbiAgICAgIHNldC5zZWVkKGN2X3NlZWQpXG4gICAgfVxuICAgIHN1YmpfZGYgPC0gc3Vial9kZltzYW1wbGUobnJvdyhzdWJqX2RmKSksICwgZHJvcCA9IEZBTFNFXVxuXG4gICAgIyBCdWlsZCBmb2xkIElEcywgc3RyYXRpZmllZCBieSByZXNwb25zZSBjbGFzc1xuICAgIEsgPC0gY3ZfS1xuICAgIGZvbGRfaWQgPC0gaW50ZWdlcihucm93KHN1YmpfZGYpKVxuICAgIGZvciAoY2xzIGluIHVuaXF1ZShzdWJqX2RmW1tyZXNwb25zZV92YXJdXSkpIHtcbiAgICAgIGlkeF9jbHMgPC0gd2hpY2goc3Vial9kZltbcmVzcG9uc2VfdmFyXV0gPT0gY2xzKVxuICAgICAgZm9sZF9pZFtpZHhfY2xzXSA8LSByZXAoc2VxX2xlbihLKSwgbGVuZ3RoLm91dCA9IGxlbmd0aChpZHhfY2xzKSlcbiAgICB9XG5cbiAgICAjIExpc3Qgb2Ygc3ViamVjdCBJRHMgcGVyIGZvbGRcbiAgICBzdWJqZWN0X2ZvbGRzIDwtIHNwbGl0KHN1YmpfZGZbW2N2X3N1YmplY3RfdmFyXV0sIGZvbGRfaWQpXG5cbiAgICAjIENvbnZlcnQgc3ViamVjdCBmb2xkcyB0byByb3cgaW5kZXggZm9sZHMgKGxvbmdpdHVkaW5hbCwgYWxsIHZpc2l0cylcbiAgICBpbmRleF9mb2xkcyA8LSBsYXBwbHkoc3ViamVjdF9mb2xkcywgZnVuY3Rpb24oc3Vial9ncnApIHtcbiAgICAgIHdoaWNoKG1ldGFbW2N2X3N1YmplY3RfdmFyXV0gJWluJSBzdWJqX2dycClcbiAgICB9KVxuXG4gICAgIyBCYWxhbmNlZCBFcnJvciBSYXRlIENWLCBzaW1pbGFyIHRvIHlvdXIgbWFudWFsIGNvZGVcbiAgICBtYXhfbmNvbXAgPC0gbWluKGN2X21heF9uY29tcCwgbmNvbXApXG4gICAgbiAgICAgICAgIDwtIG5yb3coWClcbiAgICBjbGFzc2VzICAgPC0gbGV2ZWxzKFkpXG4gICAgYmVyX2N2ICAgIDwtIG51bWVyaWMobWF4X25jb21wKVxuXG4gICAgZm9yIChuYyBpbiBzZXFfbGVuKG1heF9uY29tcCkpIHtcbiAgICAgIGNsYXNzX2Vycl9zdW0gPC0gc2V0TmFtZXMobnVtZXJpYyhsZW5ndGgoY2xhc3NlcykpLCBjbGFzc2VzKVxuICAgICAgY2xhc3Nfbl9zdW0gICA8LSBzZXROYW1lcyhudW1lcmljKGxlbmd0aChjbGFzc2VzKSksIGNsYXNzZXMpXG5cbiAgICAgIGZvciAoZm9sZCBpbiBzZXFfYWxvbmcoaW5kZXhfZm9sZHMpKSB7XG4gICAgICAgIHRlc3RfaWR4ICA8LSBpbmRleF9mb2xkc1tbZm9sZF1dXG4gICAgICAgIHRyYWluX2lkeCA8LSBzZXRkaWZmKHNlcV9sZW4obiksIHRlc3RfaWR4KVxuXG4gICAgICAgIGZpdF9mb2xkIDwtIHBsc2RhKFxuICAgICAgICAgIFhbdHJhaW5faWR4LCAsIGRyb3AgPSBGQUxTRV0sXG4gICAgICAgICAgWVt0cmFpbl9pZHhdLFxuICAgICAgICAgIG5jb21wID0gbmNcbiAgICAgICAgKVxuXG4gICAgICAgIHByZWQgPC0gbWl4T21pY3M6OnByZWRpY3QoXG4gICAgICAgICAgZml0X2ZvbGQsXG4gICAgICAgICAgWFt0ZXN0X2lkeCwgLCBkcm9wID0gRkFMU0VdLFxuICAgICAgICAgIGRpc3QgPSBcIm1heC5kaXN0XCJcbiAgICAgICAgKVxuXG4gICAgICAgICMgU2FmZWx5IHNlbGVjdCBhdmFpbGFibGUgY29tcG9uZW50XG4gICAgICAgIG5jb21wX2ZpdCA8LSBsZW5ndGgocHJlZCRjbGFzcylcbiAgICAgICAgbmNvbXBfdXNlIDwtIG1pbihuYywgbmNvbXBfZml0KVxuXG4gICAgICAgIHlfcHJlZCA8LSBwcmVkJGNsYXNzW1tuY29tcF91c2VdXVxuICAgICAgICB5X3RydWUgPC0gWVt0ZXN0X2lkeF1cblxuICAgICAgICBmb3IgKGNscyBpbiBjbGFzc2VzKSB7XG4gICAgICAgICAgaWR4X2Nsc19mb2xkIDwtIHdoaWNoKHlfdHJ1ZSA9PSBjbHMpXG4gICAgICAgICAgaWYgKGxlbmd0aChpZHhfY2xzX2ZvbGQpID09IDApIG5leHRcblxuICAgICAgICAgIGVyciA8LSBtZWFuKHlfcHJlZFtpZHhfY2xzX2ZvbGRdICE9IHlfdHJ1ZVtpZHhfY2xzX2ZvbGRdKVxuICAgICAgICAgIGNsYXNzX2Vycl9zdW1bY2xzXSA8LSBjbGFzc19lcnJfc3VtW2Nsc10gKyBlcnIgKiBsZW5ndGgoaWR4X2Nsc19mb2xkKVxuICAgICAgICAgIGNsYXNzX25fc3VtW2Nsc10gICA8LSBjbGFzc19uX3N1bVtjbHNdICAgKyBsZW5ndGgoaWR4X2Nsc19mb2xkKVxuICAgICAgICB9XG4gICAgICB9XG5cbiAgICAgIGNsYXNzX2VyciAgIDwtIGNsYXNzX2Vycl9zdW0gLyBjbGFzc19uX3N1bVxuICAgICAgYmVyX2N2W25jXSAgPC0gbWVhbihjbGFzc19lcnIsIG5hLnJtID0gVFJVRSlcbiAgICB9XG5cbiAgICBvcHRfbmNvbXAgPC0gd2hpY2gubWluKGJlcl9jdilcblxuIGZpdCA8LSBwbHNkYShYLCBZLCBuY29tcCA9IG9wdF9uY29tcCwgLi4uKVxuXG4gICMgTGF0ZW50IHZhcmlhYmxlcyAoc2NvcmVzKTogbl9jZWxscyB4IG5jb21wXG4gIExWIDwtIGZpdCR2YXJpYXRlcyRYXG4gIG9iamVjdEBwbHNkYV9MViA8LSBhcy5tYXRyaXgoTFYpXG4gIH1cblxuICBvYmplY3Rcbn1cblxuYGBgIn0= -->

```r
#' PLS-DA step for scLASER objects using mixOmics
#'
#' Wrapper around [mixOmics::plsda()] that pulls the feature matrix from a
#' `scLASER` object, constructs class labels from metadata, optionally performs
#' multilevel decomposition, fits PLS-DA, and saves the latent variables
#' (component scores) into `object@plsda_LV`.
#'
#' Optionally, a subject-wise longitudinal cross-validation (CV) is performed,
#' where folds are created at the subject level (all visits for a subject are
#' kept in the same fold) and the Balanced Error Rate (BER) is computed for
#' 1..`cv_max_ncomp` components. CV results are attached as
#' `attr(object, "plsda_cv")`.
#'
#' @param object         A `scLASER` object.
#' @param response_var   Character. Column in `object@metadata` used as class
#'   label (default `"cell_type"`).
#' @param multilevel_var Optional character. Column in `object@metadata`
#'   identifying repeated-measures units (e.g. `"sample_id"`). If not `NULL`,
#'   multilevel PLS-DA is performed via [mixOmics::withinVariation()] prior
#'   to fitting.
#' @param ncomp          Integer. Number of latent components to compute for
#'   the main PLS-DA fit and for CV (upper bound).
#' @param cv_longitudinal Logical. If `TRUE`, perform subject-wise longitudinal
#'   CV as described above. Default is `FALSE`.
#' @param cv_subject_var Character. Column in `object@metadata` that defines
#'   the subject ID for CV folds (default `"subject_id"`).
#' @param cv_K           Integer. Number of folds for longitudinal CV
#'   (default `5`).
#' @param cv_max_ncomp   Integer. Maximum number of components to evaluate in
#'   longitudinal CV. Default is `ncomp`. Values greater than `ncomp` are
#'   truncated to `ncomp`.
#' @param cv_seed        Optional integer seed for reproducible fold creation.
#'   If `NULL`, the current RNG state is used.
#' @param ...            Additional arguments passed to [mixOmics::plsda()].
#'
#' @return A `scLASER` object with:
#'   \itemize{
#'     \item `@plsda_LV` populated with the PLS-DA latent variables
#'       (scores; `n_cells x ncomp`).
#'     \item If `cv_longitudinal = TRUE`, an attribute `plsda_cv` attached:
#'       `attr(object, "plsda_cv")`, a list with elements
#'       `ber` (numeric BER per component),
#'       `opt_ncomp` (optimal number of components),
#'       and `ncomp_seq` (sequence of components evaluated).
#'   }
#'
#' @export
#' @importFrom mixOmics plsda withinVariation
mixOmics_plsda <- function(object,
                           response_var    = "cell_type",
                           multilevel_var  = NULL,
                           ncomp           = 2,
                           cv_longitudinal = FALSE,
                           cv_subject_var  = "subject_id",
                           cv_K            = 5,
                           cv_max_ncomp    = ncomp,
                           cv_seed         = NULL,
                           ...) {

  stopifnot(inherits(object, "scLASER"))

  ## 1. Feature matrix from NAM_matrix ----
  X <- object@NAM_matrix
  if (is.null(X)) {
    stop(
      "object@NAM_matrix is NULL. ",
      "Please compute NAM_matrix before calling mixOmics_plsda()."
    )
  }
  X <- as.matrix(X)

  ## 2. Response from metadata ----
  meta <- object@metadata
  if (!response_var %in% colnames(meta)) {
    stop("response_var '", response_var, "' not found in object@metadata.")
  }
  Y <- factor(meta[[response_var]])

  if (nrow(X) != nrow(meta)) {
    stop(
      "Number of rows in NAM_matrix (", nrow(X),
      ") does not match number of rows in metadata (", nrow(meta), ")."
    )
  }

  ## 3. Optional multilevel design (withinVariation) ----
  if (!is.null(multilevel_var)) {
    if (!multilevel_var %in% colnames(meta)) {
      stop("multilevel_var '", multilevel_var,
           "' not found in object@metadata.")
    }
    design <- data.frame(sample = meta[[multilevel_var]])
    X <- withinVariation(X, design = design)
  }

  ## 4. Fit main PLS-DA model ----
  # Use mixOmics::plsda; ncomp is the max number of components we keep
 # fit <- plsda(X, Y, ncomp = ncomp, ...)

  # Latent variables (scores): n_cells x ncomp
  #LV <- fit$variates$X
  #object@plsda_LV <- as.matrix(LV)

  ## 5. Optional subject-wise longitudinal CV ----
  if (isTRUE(cv_longitudinal)) {
    if (!cv_subject_var %in% colnames(meta)) {
      stop("cv_subject_var '", cv_subject_var,
           "' not found in object@metadata.")
    }

    # Prepare subject-level table (unique (subject, class) combinations)
    subj_df <- meta[, c(cv_subject_var, response_var)]
    subj_df <- subj_df[!duplicated(subj_df), , drop = FALSE]

    # Randomize order of subjects for fold assignment
    if (!is.null(cv_seed)) {
      set.seed(cv_seed)
    }
    subj_df <- subj_df[sample(nrow(subj_df)), , drop = FALSE]

    # Build fold IDs, stratified by response class
    K <- cv_K
    fold_id <- integer(nrow(subj_df))
    for (cls in unique(subj_df[[response_var]])) {
      idx_cls <- which(subj_df[[response_var]] == cls)
      fold_id[idx_cls] <- rep(seq_len(K), length.out = length(idx_cls))
    }

    # List of subject IDs per fold
    subject_folds <- split(subj_df[[cv_subject_var]], fold_id)

    # Convert subject folds to row index folds (longitudinal, all visits)
    index_folds <- lapply(subject_folds, function(subj_grp) {
      which(meta[[cv_subject_var]] %in% subj_grp)
    })

    # Balanced Error Rate CV, similar to your manual code
    max_ncomp <- min(cv_max_ncomp, ncomp)
    n         <- nrow(X)
    classes   <- levels(Y)
    ber_cv    <- numeric(max_ncomp)

    for (nc in seq_len(max_ncomp)) {
      class_err_sum <- setNames(numeric(length(classes)), classes)
      class_n_sum   <- setNames(numeric(length(classes)), classes)

      for (fold in seq_along(index_folds)) {
        test_idx  <- index_folds[[fold]]
        train_idx <- setdiff(seq_len(n), test_idx)

        fit_fold <- plsda(
          X[train_idx, , drop = FALSE],
          Y[train_idx],
          ncomp = nc
        )

        pred <- mixOmics::predict(
          fit_fold,
          X[test_idx, , drop = FALSE],
          dist = "max.dist"
        )

        # Safely select available component
        ncomp_fit <- length(pred$class)
        ncomp_use <- min(nc, ncomp_fit)

        y_pred <- pred$class[[ncomp_use]]
        y_true <- Y[test_idx]

        for (cls in classes) {
          idx_cls_fold <- which(y_true == cls)
          if (length(idx_cls_fold) == 0) next

          err <- mean(y_pred[idx_cls_fold] != y_true[idx_cls_fold])
          class_err_sum[cls] <- class_err_sum[cls] + err * length(idx_cls_fold)
          class_n_sum[cls]   <- class_n_sum[cls]   + length(idx_cls_fold)
        }
      }

      class_err   <- class_err_sum / class_n_sum
      ber_cv[nc]  <- mean(class_err, na.rm = TRUE)
    }

    opt_ncomp <- which.min(ber_cv)

 fit <- plsda(X, Y, ncomp = opt_ncomp, ...)

  # Latent variables (scores): n_cells x ncomp
  LV <- fit$variates$X
  object@plsda_LV <- as.matrix(LV)
  }

  object
}

```

<!-- rnb-source-end -->


<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



#testing show method



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYSA9IHNjTEFTRVIoKVxuYGBgIn0= -->

```r
a = scLASER()
```

<!-- rnb-source-end -->


<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->





## 

<!-- rnb-text-end -->

