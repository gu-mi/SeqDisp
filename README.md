SeqDisp
=======

### The Level of Residual Dispersion Variation and the Power of Differential Expression Tests for RNA-Seq Data

`SeqDisp` is an R package for evaluating the negative binomial dispersions methods commonly used in RNA-Seq data analysis. The methodologies are discussed in the manuscript **The Level of Residual Dispersion Variation and the Power of Differential Expression Tests for RNA-Seq Data** (by Gu Mi and Yanming Di, PLOS ONE, under review). Functions in `SeqDisp` have been used to generate all figures and tables displayed in the manuscript, plus some additional analysis tools for further investigations.

******

We provide (Dropbox) links below to download R source codes and related supporting files for preparing the datasets and reproducing figures/tables in the manuscript. Some intermediate key results are also provided when necessary.

### Data Preparations

* R codes for preparing all datasets: [Download R file](https://www.dropbox.com/s/ypdpedvi1i1unbt/Data_Preparations.R?dl=0)
* Human data supporting infomation folder: [Download files](https://www.dropbox.com/sh/7xycce8lg2fq81p/AADz8N6p3Bev8s3sIMl-aXMYa?dl=0)
* Mouse data supporting infomation folder: [Download files](https://www.dropbox.com/sh/1jsgwebcz1jqq52/AACb_ULh7HXEfQL8OlyOUUdma?dl=0)
* Zebrafish data supporting infomation folder: [Download files](https://www.dropbox.com/sh/xtq6rmxmjpqop87/AADAO32t_uauN96mdH83ZkkIa?dl=0)
* Arabidopsis data supporting infomation folder: [Download files](https://www.dropbox.com/sh/33lakumo2u1094m/AABRa8dbGRyCBe0nZ2hcFlg1a?dl=0)
* Fruit fly data supporting infomation folder: [Download files](https://www.dropbox.com/sh/geu5s45wxwrg03q/AACyBcEoIDTlmNPMThKcUIIDa?dl=0)

### Figures

* Figure 1 -- Mean-dispersion plots for the human RNA-Seq dataset: [Download R file and results](https://www.dropbox.com/sh/p3q9wtfxa1oz44q/AACTqnNZPSYsvQrwvjr3vF6La?dl=0)
* Figure 2 -- True Positive Rate (TPR) vs. False Discovery Rate (FDR) plots for the six DE test methods performed on RNA-Seq datasets simulated to mimic real datasets (fold change estimated from real data): [Download R files and results](https://www.dropbox.com/sh/id5kcn43w2nyiul/AAD2H0dDbTaukg3GG1yWq9Ara?dl=0)
* Figure 3 -- True Positive Rate (TPR) vs. False Discovery Rate (FDR) plots for the six DE test methods performed on RNA-Seq datasets simulated to mimic real datasets (fold change = 1.2): [Download R files and results](https://www.dropbox.com/sh/2qh0xd71nejsuvy/AABVSQFnLbZfLesnunUkr544a?dl=0)
* Figure 4 -- True Positive Rate (TPR) vs. False Discovery Rate (FDR) plots for the six DE test methods performed on RNA-Seq datasets simulated to mimic real datasets (fold change = 1.5): [Download R files and results](https://www.dropbox.com/sh/q57kpf6o97d8w39/AAAYd_vMe3aQgAxTgPicCF1-a?dl=0)
* Figure 5 -- True Positive Rate (TPR) vs. False Discovery Rate (FDR) plots for the six DE test methods performed on RNA-Seq dataset simulated to mimic the human dataset: [Download R files and results](https://www.dropbox.com/sh/8ryud7hdo51gx4f/AADYePWYDiF4kQs6qvS4u1oEa?dl=0)
* Figure 6 -- Histograms of p-values for the non-DE genes from the six DE test methods: [Download R file and results](https://www.dropbox.com/sh/4llg6wwuzap8yc4/AAA4xl5va4PQ-QwuZK_L7YGga?dl=0)
* Figure 7 -- Histograms of p-values for the non-DE genes from the six DE test methods: [Download R file and results](https://www.dropbox.com/sh/1k7c3xmv9l40fcv/AAAlIlMqWrYpTdTJ3PEdirxla?dl=0)
* Figure 8 -- MA plots for the trended, genewise, tagwise-trend and QLSpline methods performed on the mouse dataset: [Download R files and results](https://www.dropbox.com/sh/yu8re67rf0xg3zs/AAAW14SuhAqbHEciHlwf7NuEa?dl=0)
* Figure 9 -- Estimation accuracy of sigma.hat: [Download R file and results](https://www.dropbox.com/sh/qwefw68cgoyxkto/AABrNZaxl4yuumfq_0NByMyqa?dl=0)
* Figure 10 -- The calibration plot for estimating residual dispersion variation sigma for the mouse dataset: [Download R file and results](https://www.dropbox.com/sh/8w32v7yr7llr8bw/AACFQ9HRsO7wYIJlcGX3Ey6La?dl=0)

### Tables
* Table 1 -- The proportion of variation in log(phi_MOM) explained by the fitted gamma log-linear, quadratic and cubic regression models: [Download R file](https://www.dropbox.com/s/xpln6q881lf35z2/Table1.R?dl=0)
* Table 2 -- Estimated level of residual dispersion variation in five real RNA-Seq datasets: [Download R file and results](https://www.dropbox.com/sh/iab5id4xabvanlf/AABMi4f88fRUby0pfc01CqxGa?dl=0)
* Table 3 -- (descriptive information; no R programming)
* Table 4 -- Actual FDR for a nominal FDR of 0.1 (fold change estimated from data): [Download R file and results](https://www.dropbox.com/sh/miju0uwmh8sdmop/AAAMxlM8_pO1PgLxwYBW9ONKa?dl=0)
* Table 4 -- Actual FDR for a nominal FDR of 0.1 (fold change = 1.2): [Download R file and results](https://www.dropbox.com/sh/e1z8li11p04xdeh/AACEMbthByz6no63AHPbfJgKa?dl=0)
* Table 4 -- Actual FDR for a nominal FDR of 0.1 (fold change = 1.5): [Download R file and results](https://www.dropbox.com/sh/mjfzfw14ly0qtlw/AADweZyvZ7ILTGPAq7y9FOjua?dl=0)
* Table 5 -- (descriptive information; no R programming)
* Table 6 -- Calibrated sigma.tilde values for the five real datasets: [Download R file and results](https://www.dropbox.com/sh/jvj4m34od6fxqpd/AAByRUqkFsPqiSCztlGBtbLFa?dl=0)

******

If you have any questions, please do not hesitate to email the repository maintainer (Gu Mi) at neo.migu@gmail.com. Thank you for your interests in our research work.