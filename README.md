# Genetic-priority-score

<h2>Summary</h2>
<p>We created an <i>in-silico</i> genetic priority score that can inform drug target prioritization and validation. This score is constructed as a weighted sum of the effects of eleven genetic features (described below) on drug indication using Firth logistic regression with five-fold cross-validation and was applied to 19,365 protein-coding genes and 348 phenotypes. In our study, we discovered that drugs with a high priority score are much more likely to progress through clinical trials, and we established score thresholds that corresponded to an odds ratio (OR) >=4, OR >=5, and OR >=6 of being a therapeutic target when compared to drugs with lower ranked priority scores.</p>

<p>This web application can be used to search for putative therapeutic targets with a high genetic priority score corresponding to an OR >=4. We included drug indication data from the Open Targets Platform<sup>1</sup> and the Side effect Resource (SIDER) 4.1<sup>2</sup> and mapped all drug indications and genetic phenotypes to phecode integer terms. We have provided evidence of each genetic association, and users can search this website by a gene target or phenotype of interest. We note that these results should be a starting point for users to conduct follow-up analyses, and if a gene and phenotype are not present in this table, this should not be interpreted as a drug target with genetic support against it.</p>

<p>For further details about the methods and analysis behind this study, please see our recently published paper: insert ref. </p>
   

 <h2>Instructions</h2>
 
<p>Users can use this website to explore associations by gene or by phenotype. All associations have a genetic priority score that corresponds to an OR >= 4. This OR cutoff can be refined to OR >= 5 or an OR >= 6 using the OR cutoff slider.</p>
<p>When the Gene tab is selected, users can search the priority scores for their gene of choice in the <b>‘Search Gene Target‘</b> tab. The user can choose to show all phenotype results or select a phenotype of interest in the <b>‘Choose Phenotypes:’</b> dropdown bar. When the phenotype tab is selected, users can select their phenotype of interest from the dropdown <b>‘Select phenotype:’</b> dropdown bar and further restrict by the gene target in the the <b>‘Choose Gene:’</b> dropdown bar.

<p>Once a gene or phenotype is searched on the left side panel, this will generate a table of results in the <b>Table</b> tab. Users can click on the <b>Evidence</b> tab to generate a heatmap with the genetic predictors on the x-axis and the associated gene/phenotypes on the y-axis. Genetic evidence for a predictor is colored in purple, and the user can hover over the cell to get the phenotype description from the phenotype data source. For genes/phenotypes with more than 100 rows of evidence, only the top 100 high-scored genes/phenotypes are shown.</p>

<h2>Genetic Evidence</h2>
    
<p>Using publicly available data sources, we collected eleven genetic features from six types of genetic evidence (clinical variants, coding variants, genome-wide association studies (GWAS) phenotypes, tissue phenotypes, tissue specificity, and constraint). We mapped these features to phecodes and restricted these to phecode integer terms. Phecodes that mapped to phecode categories ‘neoplasms’, ‘infectious diseases’, ‘pregnancy complications’, ‘injuries and poisonings’, and ‘null’ were excluded from the analysis. Each of these genetic features is described below.</p>

<h3>Clinical variants </h3>
              
              
<p><b>OMIM:</b> Mendelian genes from the Online Mendelian Inheritance in Man (OMIM)<sup>3</sup> database.</p>
    
<p><b>HGMD:</b> Disease-causing and likely disease-causing genes for human inherited diseases from the human mutation database (HGMD)<sup>4</sup>.</p>

<p><b>EVA-ClinVar:</b> Clinically relevant genetic variants and diseases from ClinVar<sup>5</sup>. This evidence is obtained from the Open Target platform and was curated by European Variation Archive (EVA)<sup>6</sup>. ClinVar evidence was filtered on clinical significance: likely pathogenic', 'association', 'confers sensitivity', 'drug response', 'protective', and 'pathogenic' and the confidence of the submission was filtered on: 'criteria provided, single submitter', 'criteria provided, conflicting interpretations', 'criteria provided, multiple submitters, no conflicts', 'reviewed by expert panel', and 'practice guideline'.</p>

<h3>Coding variants</h3>
              
<p><b>Genebass:</b> Exome-based association statistics from the UK Biobank hosted by Genebass<sup>7</sup>. We restricted to traits labeled ‘ICD first occurrence’ and extracted genome-wide significant gene associations annotated as either predicted loss of function (pLOF) or missense.</p>
 <h3>Genome-wide association study phenotypes</h3>

<p><b>Coloc phenotype:</b> Genes with a GWA phenotype driven by gene expression regulation through colocalization, using a posterior probability (PPH4)>0.8<sup>8,9</sup>. These phenotypes comprised 1,071 case-control phenotypes from phecodes<sup>10</sup>, and 779 phenotypes from the UK Biobank Neale dataset (Round 1 - http://www.nealelab.is/uk-biobank).</p>

<p><b>eQTL phenotype:</b> Genes with a GWA phenotype driven by gene expression regulation through shared variants. These phenotypes were obtained from the Pan-UK Biobank<sup>11</sup> and consisted of 484 ICD-10 traits, 729 phecodes, and 39 continuous traits that we had manually mapped to phecodes from the UK Biobank Neale.</p>

<p><b>Locus2gene:</b> Genes with a GWA phenotype identified by the Locus2Gene machine learning model from the Open Targets Genetics Portal which prioritizes likely causal genes at each GWAS locus by integrating fine-mapped associations with functional genomic features<sup>12</sup>.</p>

<h3>Tissue phenotypes </h3>
<p><b>Coloc2 and eQTL tissue:</b> Genes with gene expression regulation in a given tissue that underlies a significant genome-wide association through a shared GWA and eQTL signal. We used eQTL summary statistics from GTEX v7 for coloc2 tissue and eQTL summary statistics from GTEX v8 for eQTL tissue. To link the coloc2 and eQTL tissues to phecodes, we used LD score regression applied to specifically expressed genes (LDSC-SEG)<sup>13</sup>, to identify phenotypic relevant tissues. This was applied for 54 tissues from GTEX V8 using GWAS summary statistics for 776 traits from the Neale GWAS resource and 1,112 traits from the Saige GWAS resource.</p>

<h3>Tissue specificity </h3>
              
<p><b>Tissue specificity of gene expression:</b> Tissue specific genes were defined using two sources: TPM, a measure of gene expression per tissue<sup>14</sup>, and Tau, a gene tissue specificity metric<sup>15</sup>. We dichotomized levels of TPM and used a threshold > 0.5 to indicate the presence of gene expression and used a threshold of > 0.8 to dichotomize Tau and indicate tissue specificity<sup>15</sup>. Median TPM values for each gene across all tissue samples and gene counts were obtained from GTEX V8.</p>

<h3>Constraint</h3>
              
<p><b>OE:</b> Gene constraint was measured by the OE (observed/expected) score using the suggested cutoff of < 0.35 of the upper bound of the OE confidence interval<sup>16</sup>.</p>
 

<h2>Citation</h2>

<p></p>


<h2>References</h2>

<p>1.	Ochoa, D. et al. Open Targets Platform: supporting systematic drug-target identification and prioritisation. Nucleic Acids Res 49, D1302-D1310 (2021).</p>
<p>2.	Kuhn, M., Letunic, I., Jensen, L.J. & Bork, P. The SIDER database of drugs and side effects. Nucleic acids research 44, D1075-D1079 (2016).</p>
<p>3.	Hamosh, A. et al. Online Mendelian Inheritance in Man (OMIM), a knowledgebase of human genes and genetic disorders. Nucleic Acids Res 30, 52-5 (2002).</p>
<p>4.	Stenson, P.D. et al. The Human Gene Mutation Database (HGMD(®)): optimizing its use in a clinical diagnostic or research setting. Hum Genet 139, 1197-1207 (2020).</p>
<p>5.	Landrum, M.J. et al. ClinVar: improving access to variant interpretations and supporting evidence. Nucleic Acids Research 46, D1062-D1067 (2017).</p>
<p>6.	Cook, C.E. et al. The European Bioinformatics Institute in 2016: Data growth and integration. Nucleic Acids Research 44, D20-D26 (2015).</p>
<p>7.	Karczewski, K.J. et al. Systematic single-variant and gene-based association testing of thousands of phenotypes in 394,841 UK Biobank exomes. Cell Genomics 2, 100168 (2022).</p>
<p>8.	Dobbyn, A. et al. Landscape of Conditional eQTL in Dorsolateral Prefrontal Cortex and Co-localization with Schizophrenia GWAS. Am J Hum Genet 102, 1169-1184 (2018).</p>
<p>9.	Giambartolomei, C. et al. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. PLoS Genet 10, e1004383 (2014).</p>
<p>10.	Zhou, W. et al. Efficiently controlling for case-control imbalance and sample relatedness in large-scale genetic association studies. Nat Genet 50, 1335-1341 (2018).</p>
<p>11.	Pan-UKB team https://pan.ukbb.broadinstitute.org. (2020).</p>
<p>12.	Mountjoy, E. et al. An open approach to systematically prioritize causal variants and genes at all published human GWAS trait-associated loci. Nat Genet 53, 1527-1533 (2021).</p>
<p>13.	Finucane, H.K. et al. Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types. Nature Genetics 50, 621-629 (2018).</p>
<p>14.	GTEx Portal https://www.gtexportal.org/home/datasets.</p>
<p>15.	Kryuchkova-Mostacci, N. & Robinson-Rechavi, M. A benchmark of gene expression tissue-specificity metrics. Brief Bioinform 18, 205-214 (2017).</p>
<p>16.	Karczewski, K.J. et al. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature 581, 434-443 (2020).</p>

 
