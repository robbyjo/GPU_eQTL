# GPU eQTL

## What it is
To compute trillions of association typically required for quantitative trait locus (QTL) analysis. The statistical model is linear model. This program uses graphical processing unit (GPU) to achieve dramatic speedups. Created in 2009.

## Requirements
* Java version 7 or above
* NVIDIA or AMD graphic card, with the corresponding drivers installed

## How to run through command line
Invoke the following:

    java -jar gpgpu.jar:qgenerics.jar:javacl-1.0.0-RC1-shaded.jar:jdistlib-0.3.5-bin.jar:javacsv.jar:commons-compress-1.3.jar:jgrapht.jar gov.nih.eqtl.QeQTLAnalysis your_file.ini

The structure of `your_file.ini` is the following:

    genotype_file = your_genotype_file.csv
    expression_file = your_expression_file.csv
    output_file = your_output_file.csv
    covariate_file = your_covariate_file.csv
    covariate_fixed = Sex Age Batch PC1 PC2 PC3 PC4 PC5 PCP1 PCP2 PCP3 PCP4 PCP5
    covariate_factor = Sex Batch
    genotype_format = csv
    expression_format = csv
    df_offset = 0
    block_size = 10240
    num_threads = 8
    lambda = 0.75
    library_path = /usr/lib:/usr/lib64
    threshold = pval 1e-4

### Explanation:

`your_genotype_file.csv` is a CSV file with the number of genotypes (row) x number of samples (column). The extra first column must be an identifier of each row. Cannot have any missing values. Currently no other format is accepted.

`your_expression_file.csv` is a CSV file with the number of expression (row) x number of samples (column). The extra first column must be an identifier of each row. Cannot have any missing values. Currently no other format is accepted.

`your_covariate_file.csv` is a CSV file with the number of samples (row) x number of covariates (column). Cannot have any missing values on covariates that are used in the model. Currently no other format is accepted.

**NOTE:** Please make sure that the ordering of the samples match 100% between the genotype, expression, and covariate files because this program did NOT perform additional check other than the number of samples.

In reality, you can put in any Omics file modalities in either `genotype_file` or `expression_file` to match your research objective. For example, if you want to perform methylation QTL analysis, then you put in methylation file as the `expression_file`. If you want to perform splicing-protein analysis, then you put in splicing file as the `genotype_file` and proteomics data as the `expression_file`, and so on.

`covariate_fixed` is the covariate names you want to include in the model. Cannot have any missing values. The names must be an exact match in the covariate file.

`covariate_factor` is part of `covariate_fixed` that this program needs to treat as factors, not as numeric values. Covariate with more than two factors will be encoded using [One-hot encoding](https://en.wikipedia.org/wiki/One-hot).

`df_offset` is useful if you want to further reduce the degrees of freedom (perhaps because you already regress out several problematic covariates) and you want to account those factors into the DF computation. If you did not do any of that, leave it to zero.

`block_size` is the matrix size you want the GPU to handle at a time. Depends on the size of the GPU. 10240 is a good size for old GPUs.

`num_threads` is how many CPU thread you would like to use to perform the computation. Some of the computation (e.g., P-value) is done in CPU. So, values at least two would be good.

`lambda` leave it at 0.75. The intent was to use that much proportion of GPU memory, but evidently it could be unpredictable.

`library_path` is to specify where libOpenCL.so may reside. Linux only. Remove this line for Windows. This is to account for the various ways Linux distributions store the appropriate library.

`threshold` is to record anything below or above a certain threshold. The `pval 1e-4` is to record all the results with p<1e-4. You can have `rsq 0.1`, which means record all results with partial R^2 of at least 0.1.

`your_output_file.csv` is the CSV file of the output format generated. The column description is as follows:
1. The first column will contain the ID of the genotype file.
2. ID of the expression file
3. Partial R^2 of the association
4. The regression (beta) coefficient
5. T statistics
6. log10 of the P-value

## Running in Windows
This program should automatically detect the appropriate drivers.

## Running in Linux (Ubuntu)
You need to install the appropriate drivers. For NVIDIA cards, you need to install `nvidia-driver-<version>`, `nvidia-dkms-<version>`, `nvidia-cuda-toolkit`, and `clinfo`. Check your Linux manual for the appropriate driver version. The package names may vary for other Linux distributions.

For AMD driver, please consult your system administrator.

**IMPORTANT**: You **MUST** reboot your Linux computer upon installing these drivers.



