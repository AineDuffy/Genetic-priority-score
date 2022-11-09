# Genetic-priority-score

<h2>Summary</h2>
<p>We created an <i>in-silico</i> genetic priority score that can inform drug target prioritization and validation. This score is constructed as a weighted sum of the effects of eleven genetic features (described below) on drug indication using Firth logistic regression with five-fold cross-validation and was applied to 19,365 genes and 348 phenotypes. In our study, we discovered that drugs with a high priority score are much more likely to progress through clinical trials, and we established score thresholds that corresponded to an odds ratio (OR) >=4, OR >=5, and OR >=6 of being a therapeutic target when compared to drugs with lower ranked priority scores.</p>

<p>This web application can be used to search for putative therapeutic targets with a high genetic priority score corresponding to an OR >=4. We included drug indication data from the Open Targets Platform<sup>1</sup> and the Side effect Resource (SIDER) 4.1<sup>2</sup> and mapped all drug indications and genetic phenotypes to phecode integer terms. We have provided evidence of each genetic association, and users can search this website by a gene target or phenotype of interest. We note that these results should be a starting point for users to conduct follow-up analyses, and if a gene and phenotype are not present in this table, this should not be interpreted as a drug target with genetic support against it.</p>

<p>For further details about the methods and analysis behind this study, please see our recently published paper: insert ref. </p>
   

 <h2>Instructions</h2>
 
<p>Users can use this website to explore associations by gene or by phenotype. All associations have a genetic priority score that corresponds to an OR >= 4. This OR cutoff can be refined to OR >= 5 or an OR >= 6 using the OR cutoff slider.</p>
<p>When the Gene tab is selected, users can search the priority scores for their gene of choice in the <b>‘Search Gene Target‘</b> tab. The user can choose to show all phenotype results or select a phenotype of interest in the <b>‘Choose Phenotypes:’</b> dropdown bar. When the phenotype tab is selected, users can select their phenotype of interest from the dropdown <b>‘Select phenotype:’</b> dropdown bar and further restrict by the gene target in the the <b>‘Choose Gene:’</b> dropdown bar.

<p>Once a gene or phenotype is searched on the left side panel, this will generate a table of results in the", <b>Table</b> tab. Users can click on the <b>Evidence</b> tab to generate a heatmap with the genetic predictors on the x-axis and the associated gene/phenotypes on the y-axis. Genetic evidence for a predictor is colored in purple, and the user can hover over the cell to get the phenotype description from the phenotype data source. For genes/phenotypes with more than 100 rows of evidence, only the top 100 high-scored genes/phenotypes are shown.</p>

<h2>Genetic Evidence</h2>
    
<p>We collected genetic features from six types of genetic evidence (clinical variants, coding variants, genome-wide association studies (GWAS) phenotypes, tissue phenotypes, tissue specificity, and constraint) to investigate the association between human genetic variation in drug target genes and drug indications. We collected this data from several publicly available sources. To allow comparison between the drug and genetic datasets, we mapped all genetic phenotypes to phecodes and restricted these to phecode integer terms. All phecodes mapped to the phecode categories ‘neoplasms’, ‘infectious diseases’, ‘pregnancy complications’, ‘injuries and poisonings’, and ‘null’ were excluded from the analysis. We restricted each gene set to protein-coding genes, for which we obtained a list of 19,365 protein-coding genes from Ensembl<sup>3</sup>. Each of these genetic features is described below.</p>

<h3>Clinical variants </h3>
              
              
<p><b>OMIM:</b>The Online Mendelian Inheritance in Man (OMIM)<sup>4</sup> database contains genes involved in rare mendelian diseases. We extracted 15,846 genes associated with 7,053 Mendelian traits, mapped to 7,216 Human phenotype ontology (HPO) terms from the OMIM database. We used a curated mapping file from the phenome-wide association studies (PheWAS) catalog<sup>5</sup>to map the HPO terms to phecodes. In total, 5,792 phenotypes were mapped to 289 phecode integers associated with 4,147 genes.</p>
    
<p><b>HGMD:</b>The Human Gene Mutation Database (HGMD)<sup>6</sup> contains gene mutations underlying human inherited diseases, curated from the scientific literature. We extracted 12,795 disease-causing, and likely disease-causing genes associated with 16,810 phenotypes mapped to 1,722 Unified Medical Language System (UMLS) codes. We mapped the UMLS codes to phecodes via International Classification of Diseases (ICD) codes using the available UMLS Metathesaurus file and ICD10CM and ICD10 files from the PheWAS resource<sup>5</sup>. There were 271 phecode integers associated with 11,642 genes.</p>

<p><b>EVA-ClinVar:</b>The Open Target Platform summarizes genetic evidence from multiple datasets containing associations between targets and disease, described by Experimental Factor Ontology (EFO) terms.<sup>1</sup> One source is ClinVar<sup>7</sup>, a curated database of associations between clinically relevant genetic variants, diseases, and the variant’s clinical significance. This data is processed and curated by European Variation Archive (EVA)<sup>8</sup> before being added to the Open Target platform. We extracted ClinVar evidence based on two steps. First, evidence was filtered on clinical significance terms: likely pathogenic', 'association', 'confers sensitivity', 'drug response', 'protective', and 'pathogenic'. Secondly, we filtered based on the confidence of the submission assigned as: 'criteria provided, single submitter', 'criteria provided, conflicting interpretations', 'criteria provided, multiple submitters, no conflicts', 'reviewed by expert panel', and 'practice guideline'. We mapped the EFO terms to phecodes using the previously mentioned disease mapping files from the Open Targets Platform, EBISPOT OLS, UMLS, and PheWAS resource. In total, there were 2,772 genes associated with 3,737 phenotypes that mapped to 222 phecode integers.</p>

<h3>Coding variants</h3>
              
<p><b>Genebass:</b>Gene-based Association Summary Statistics (Genebass)<sup>9</sup> is a resource of exome-based association statistics from the UK Biobank, encompassing 4,529 phenotypes across 426,370 individuals with gene-based and single-variant testing. We restricted to traits labeled ‘ICD first occurrence’ and extracted genome-wide significant gene associations annotated as either predicted loss of function (pLOF) or missense. We extracted 18,448 genes associated with 413 ICD codes that mapped to 234 integer phecodes.</p>
 <h3>Genome-wide association study phenotypes</h3>



<p><b>Coloc phenotype:</b>We defined coloc phenotype as genes with a GWA phenotype driven by gene expression regulation through colocalization, using a posterior probability (PPH4) > 0.8<sup>10,11</sup>. We used the results from our coloc2 analysis from our previous study<sup>12</sup>, which consists of putative causal genes and their corresponding phenotype and tissue from co-localization analysis, performed across 1,850 phenotypes from the UK Biobank<sup>13</sup> and expression quantitative trait loci (eQTL) data from 48 tissues from the Genotype-Tissue Expression (GTEx) v7 project<sup>14,15</sup> using the coloc2 method. The description of this approach has been previously described<sup>12</sup>. These phenotypes consisted of 1,071 case-control phenotypes from phecodes<sup>16</sup>, and 779 phenotypes from the UK Biobank Neale dataset (Round 1 - http://www.nealelab.is/uk-biobank). Of these Neale phenotypes, we mapped 112 ICD codes to phecodes and manually mapped 387 continuous traits to phecodes. There was a total of 15,990 coloc2 genes and 313 integer phecodes.</p>

<p><b>eQTL phenotype:</b>We further generated an eQTL phenotype feature, using the updated Pan-UK Biobank<sup>17</sup>, as genes with a GWA phenotype driven by gene expression regulation through shared variants. The Pan-UK Biobank is a multi-ancestry genetic analysis of the UK Biobank across six continental ancestry groups and 7,228 phenotypes. All variants were obtained on hg19 coordinates, and we used LiftOver<sup>18</sup> to convert these to hg38 for compatibility and extracted all genome-wide significant single nucleotide variants ( (span(em("P"))),span(" < 5 x 10-8). We restricted the trait types and obtained summary statistics for 1,252 phenotypes: 484 ICD-10 traits, 729 phecodes, and 39 continuous traits that we had manually mapped to phecodes from the UK Biobank Neale. To identify the putative causal genes to these traits, we used eQTLs from GTEX v8<sup>14</sup> across 49 tissues and extracted significant variants for each tissue with a nominal P less than the gene level threshold. We intersected these variants with the Pan-UK Biobank significant variants to obtain the linked gene for each trait. We mapped all ICD traits to phecodes, and across the 49 tissues, there were 4,446 genes associated with 184 integer phecodes.</p>

<p><b>Locus2gene:</b>The Locus2Gene model is a machine learning method from the Open Targets Genetics Portal which prioritizes likely causal genes at each GWAS locus by integrating fine-mapped associations with functional genomic features<sup>19</sup>. A set of gold standard positive (GSP) genes are used in this method where the causal gene has been confidently assigned. Many of these GSP genes are selected based on known drug–target disease pairs, so we removed these from the final output. To define the Locus2Gene predictor, we restricted to genes prioritized by the model with a score greater than 0.5. We mapped the disease traits to phecodes using the disease mapping files from the Open Targets Platform, EBISPOT OLS, UMLS, and PheWAS resources. We extracted 3,709 genes classified as casual for 1,730 traits mapped to 270 phecode integers.</p>

<h3>Tissue phenotypes </h3>
<p><b>Coloc2 and eQTL tissue:</b>The coloc2 and eQTL tissue features define genes with gene expression regulation in a given tissue that underlies a significant genome-wide association through a shared GWA and eQTL signal. To link the coloc2 and eQTL tissues to phecodes, we used LD score regression applied to specifically expressed genes (LDSC-SEG)<sup>20</sup> across 776 traits from the Neale GWAS resource and 1,112 traits from the Saige GWAS resource, using GWAS summary statistics that were mapped to phecodes. We used LDSC-SEG to estimate the heritability enrichment of specifically expressed genes across 54 GTEX V8 tissues to identify phenotypic relevant tissues. In brief, this method first classifies the top 10% of ‘specifically expressed genes’ in each tissue using a t-statistic test and extends each transcribed region by a 100kb window to create tissue-specific annotations. Stratified LD score regression is then applied to the GWAS summary statistics to test the contribution of the tissue-specific annotations to trait heritability, conditional on the baseline model, which includes 52 annotations including genic, enhancer, and conserved regions, and all genes. We used a P threshold of 0.01 to identify phenotype-relevant tissues, for which there were 421 phecodes across 54 tissues (mean =15.31 phecodes per tissue, sd=6.40). These phecodes were used to represent the phenotype-relevant tissues from our coloc2 and eQTL results.</p>

<h3>Tissue specificity </h3>
              
<p><b>Tissue specificity of gene expression:</b>We constructed a tissue specificity score of gene expression for each gene-phenotype. We used transcripts per million (TPM), a measure of gene expression per tissue, and Tau, a gene tissue specificity metric. We obtained median TPM values for each gene across all tissue samples from GTEX V8. We used LDSC-SEG to map the TPM values per tissue to phecodes and took the median TPM when multiple tissues were mapped. We dichotomized levels of median TPM and used a TPM threshold > 0.5 to indicate the presence of gene expression. We calculated Tau for each gene using gene counts from GTEX V8, as described by Kryuchkova Mostacci and Robinson Rechavi<sup>21</sup>. We used a threshold of > 0.8 to dichotomize Tau and indicate tissue specificity. Using these dichotomized values, we created a genetic feature for tissue specificity of gene expression for each gene-phecode with a TPM > 0.5 and a Tau > 0.8. We identified 6,344 genes as tissue-specific across 207 phecode integers (mean=86.42 phecode integers per gene, sd=74.21).</p>

<h3>Constraint</h3>
              
<p><b>OE:</b>We used the observed/expected (OE) score to measure the constraint per gene and its intolerance to loss of function mutations. We used the suggested cutoff of < 0.35 of the upper bound of the OE confidence interval to identify genes under high constraint<sup>22</sup> and identified 2,854 genes.</p>
 

<h2>Citation</h2>

<p></p>


<h2>References</h2>

<p>1.	Ochoa, D. et al. Open Targets Platform: supporting systematic drug-target identification and prioritisation. Nucleic Acids Res 49, D1302-D1310 (2021).</p>
<p>2.	Kuhn, M., Letunic, I., Jensen, L.J. & Bork, P. The SIDER database of drugs and side effects. Nucleic acids research 44, D1075-D1079 (2016).</p>
<p>3.	Cunningham, F. et al. Ensembl 2022. Nucleic Acids Research 50, D988-D995 (2021).</p>
<p>4.	Hamosh, A. et al. Online Mendelian Inheritance in Man (OMIM), a knowledgebase of human genes and genetic disorders. Nucleic Acids Res 30, 52-5 (2002).</p>
<p>5.	Denny, J.C. et al. Systematic comparison of phenome-wide association study of electronic medical record data and genome-wide association study data. Nature Biotechnology 31, 1102-1111 (2013).</p>
<p>6.	Stenson, P.D. et al. The Human Gene Mutation Database (HGMD(®)): optimizing its use in a clinical diagnostic or research setting. Hum Genet 139, 1197-1207 (2020).</p>
<p>7.	Landrum, M.J. et al. ClinVar: improving access to variant interpretations and supporting evidence. Nucleic Acids Research 46, D1062-D1067 (2017).</p>
<p>8.	Cook, C.E. et al. The European Bioinformatics Institute in 2016: Data growth and integration. Nucleic Acids Research 44, D20-D26 (2015).</p>
<p>9.	Karczewski, K.J. et al. Systematic single-variant and gene-based association testing of thousands of phenotypes in 394,841 UK Biobank exomes. Cell Genomics 2, 100168 (2022).</p>
<p>10.	Dobbyn, A. et al. Landscape of Conditional eQTL in Dorsolateral Prefrontal Cortex and Co-localization with Schizophrenia GWAS. Am J Hum Genet 102, 1169-1184 (2018).</p>
<p>11.	Giambartolomei, C. et al. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. PLoS Genet 10, e1004383 (2014).</p>
<p>12.	Duffy, A. et al. Tissue-specific genetic features inform prediction of drug side effects in clinical trials. Sci Adv 6(2020).</p>
<p>13.	Bycroft, C. et al. The UK Biobank resource with deep phenotyping and genomic data. Nature 562, 203-209 (2018).</p>
<p>14.	Aguet, F. et al. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 369, 1318-1330 (2020).</p>
<p>15.	GTEx Portal https://www.gtexportal.org/home/datasets.</p>
<p>16.	Zhou, W. et al. Efficiently controlling for case-control imbalance and sample relatedness in large-scale genetic association studies. Nat Genet 50, 1335-1341 (2018).</p>
<p>17.	Pan-UKB team https://pan.ukbb.broadinstitute.org. (2020).</p>
<p>18.	Kuhn, R.M., Haussler, D. & Kent, W.J. The UCSC genome browser and associated tools. Briefings in Bioinformatics 14, 144-161 (2012).</p>
<p>19.	Mountjoy, E. et al. An open approach to systematically prioritize causal variants and genes at all published human GWAS trait-associated loci. Nat Genet 53, 1527-1533 (2021).</p>
<p>20.	Finucane, H.K. et al. Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types. Nature Genetics 50, 621-629 (2018).</p>
<p>21.	Kryuchkova-Mostacci, N. & Robinson-Rechavi, M. A benchmark of gene expression tissue-specificity metrics. Brief Bioinform 18, 205-214 (2017).</p>
<p>22.	Karczewski, K.J. et al. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature 581, 434-443 (2020).</p>
 
