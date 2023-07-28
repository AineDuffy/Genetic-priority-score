# Genetic-priority-score

<h2>Summary</h2>
<p>We created an <i>in-silico</i> genetic priority score that can inform drug target prioritization and validation. This score is constructed as a weighted sum of the effects of eight phenotypic specific features (described below) on drug indication using Firth logistic regression with five-fold cross-validation and was applied to 19,365 protein-coding genes and 399 phenotypes. We further incorporated the direction of genetic effect for each predictor using predictions of loss-of-function (LOF) and gain-of-function (GOF) from LoGoFunc<sup>1</sup> for clinical variants and quantitative trait loci estimates for genome-wide association studies (GWAS) phenotypes and created a complementary genetic priority score with the direction of effect (GPS-D). In our study, we discovered that drugs with a high genetic priority score are much more likely to progress through clinical trials and we recommend using a score cut-off > 1.5 that corresponds to an odds ratio (OR) >=9 of being a therapeutic target when compared to drugs with no genetic evidence.</p>

<p>This web application can be used to search for putative therapeutic targets with genetic priority scores and evidence for direction of genetic effect. We included drug indication data from the Open Targets Platform<sup>2</sup> and the Side effect Resource (SIDER) 4.1<sup>3</sup> and mapped all drug indications and genetic phenotypes to phecode integer terms. We have provided evidence of each genetic association, and users can search this website by a gene target or phenotype of interest. We note that these results should be a starting point for users to conduct follow-up analyses, and if a gene and phenotype are not present in this table, this should not be interpreted as a drug target with genetic support against it.</p>

<p>For further details about the methods and analysis behind this study, please see our paper: Duffy, A et al. Development of a human genetics-guided priority score for 19,365 genes and 399 drug indications. Submitted. </p>
   
<p>The Genetic Priority Score R shiny app is accessible at: https://rstudio-connect.hpc.mssm.edu/geneticpriorityscore/ </p>


 <h2>Instructions</h2>
 
<p>Users can use this website to explore associations by gene or by phenotype. All associations have evidence from at least one genetic association and this evidence can be refined using the score cutoff slider at 0.3 increment thresholds of the genetic priority score.</p>
<p>When the Gene tab is selected, users can search the priority scores for their gene of choice in the <b>‘Search Gene Target‘</b> tab. When the phenotype tab is selected, users can select their phenotype of interest from the dropdown <b>‘Select phenotype:’</b> dropdown bar. 

<p>Once a gene or phenotype is searched on the left side panel, this will generate a table of results in the <b>Table</b> tab. Users can click on the <b>Evidence</b>  and the <b>DOE Evidence</b> tab to generate heatmaps of the associations with the genetic predictors on the X-axis and the associated gene/phenotypes on the Y-axis. The <b>Evidence</b> tab generates a heatmap with the genetic evidence for a predictor colored in purple and the user can hover over the cell to get the phenotype description from the phenotype data source. The <b>DOE Evidence</b> tab generates a heatmap with genetic evidence predicted as loss-of-function (LOF) colored in pink and predictions of gain-of-function (GOF) colored in blue. For genes/phenotypes with more than 100 rows of evidence, only the top 100 high-scored genes/phenotypes are shown.</p>

<h2>Genetic Evidence</h2>
    
<p>Using publicly available data sources, we collected eight genetic features from three types of genetic evidence (clinical variants, coding variants and GWAS phenotypes). We mapped these features to phecodes and restricted these to phecode integer terms. Phecodes that mapped to phecode categories: ‘infectious diseases’, ‘pregnancy complications’, ‘injuries and poisonings’ and ‘null’ were excluded from the analysis. Each of these genetic features is described below.</p>

<h3>Clinical variants </h3>
              
              
<p><b>OMIM:</b> Mendelian genes from the Online Mendelian Inheritance in Man (OMIM)<sup>4</sup> database.</p>
    
<p><b>HGMD:</b> Disease-causing and likely disease-causing genes for human inherited diseases from the human mutation database (HGMD)<sup>5</sup>.</p>

<p><b>EVA-ClinVar:</b> Clinically relevant genetic variants and diseases from ClinVar<sup>6</sup>. This evidence is obtained from the Open Target platform and was curated by European Variation Archive (EVA)<sup>7</sup>. ClinVar evidence was filtered on clinical significance: 'likely pathogenic', 'association', 'confers sensitivity', 'drug response', 'protective', and 'pathogenic' and the confidence of the submission was filtered on: 'criteria provided, single submitter', 'criteria provided, conflicting interpretations', 'criteria provided, multiple submitters, no conflicts', 'reviewed by expert panel', and 'practice guideline'.</p>

<h3>Coding variants</h3>
              
<p><b>Gene burden:</b> Exome-based gene burden association statistics from the UK Biobank hosted by Genebass<sup>8</sup>. We extracted traits labeled ‘ICD first occurrence’ and restricted to missense and pLOF burden set ( <i>P</i> < 4.3 x 10<sup>-7</sup> ). </p>

<p><b>Single-Variant:</b> Exome-based single variant association statistics from the UK Biobank hosted by Genebass<sup>8</sup>. We extracted traits labeled ‘ICD first occurrence’ and restricted to missense and pLOF variants ( <i>P</i> < 4.3 x 10<sup>-7</sup> ). </p>
 

<h3>Genome-wide association study phenotypes</h3>

<p><b>eQTL phenotype:</b> Genes with a GWA phenotype driven by gene expression regulation through shared variants. These phenotypes were obtained from the Pan-UK Biobank<sup>9</sup> and consisted of 484 ICD-10 traits, 729 phecodes, and 39 continuous traits that we had manually mapped to phecodes from the UK Biobank Neale.</p>

<p><b>Locus2gene:</b> Genes with a GWA phenotype identified by the Locus2Gene machine learning model from the Open Targets Genetics Portal which prioritizes likely causal genes at each GWAS locus by integrating fine-mapped associations with functional genomic features<sup>10</sup>.</p>

<p><b>pQTL phenotype:</b> Genes with a GWA phenotype that are in high linkage disequilibrium (r<sup>2</sup>>0.80) with at least one cis or trans pQTL from supplementary table 12 from the large-scale plasma proteome study conducted by Ferkingstad et al<sup>11</sup>.</p>
 

<h2>Citation</h2>

<p></p>


<h2>References</h2>

<p>1. David, S. et al. Genome-wide prediction of pathogenic gain- and loss-of-function variants from ensemble learning of a diverse feature set. bioRxiv 2022.06.08.495288 (2022).</p>
<p>2.	Ochoa, D. et al. Open Targets Platform: supporting systematic drug-target identification and prioritisation. Nucleic Acids Res 49, D1302-D1310 (2021).</p>
<p>3.	Kuhn, M., Letunic, I., Jensen, L.J. & Bork, P. The SIDER database of drugs and side effects. Nucleic acids research 44, D1075-D1079 (2016).</p>
<p>4.	Hamosh, A. et al. Online Mendelian Inheritance in Man (OMIM), a knowledgebase of human genes and genetic disorders. Nucleic Acids Res 30, 52-5 (2002).</p>
<p>5.	Stenson, P.D. et al. The Human Gene Mutation Database (HGMD(®)): optimizing its use in a clinical diagnostic or research setting. Hum Genet 139, 1197-1207 (2020).</p>
<p>6.	Landrum, M.J. et al. ClinVar: improving access to variant interpretations and supporting evidence. Nucleic Acids Research 46, D1062-D1067 (2017).</p>
<p>7.	Cook, C.E. et al. The European Bioinformatics Institute in 2016: Data growth and integration. Nucleic Acids Research 44, D20-D26 (2015).</p>
<p>8.	Karczewski, K.J. et al. Systematic single-variant and gene-based association testing of thousands of phenotypes in 394,841 UK Biobank exomes. Cell Genomics 2, 100168 (2022).</p>
<p>8. Dobbyn, A. et al. Landscape of Conditional eQTL in Dorsolateral Prefrontal Cortex and Co-localization with Schizophrenia GWAS. Am J Hum Genet 102, 1169-1184 (2018).</p>
<p>9. Pan-UKB team https://pan.ukbb.broadinstitute.org. (2020).</p>
<p>10. Mountjoy, E. et al. An open approach to systematically prioritize causal variants and genes at all published human GWAS trait-associated loci. Nat Genet 53, 1527-1533 (2021).</p>
<p>11. Ferkingstad, E. et al. Large-scale integration of the plasma proteome with genetics and disease. Nat Genet 53 1712-1721 (2021).</p>	
 
