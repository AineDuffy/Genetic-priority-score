abouttitle='About'
Aboutinfo_page <- tabPanel(
    theme=shinythemes::shinytheme('flatly'),
    title='About',
    h3("Summary"),
    wellPanel(
        p(paste("We created an"),HTML(paste0(tags$em("in-silico "))),("genetic priority score that can inform drug target prioritization and validation for"),paste(comma(length(unique(Drugfile$Gene)))), paste("genes and"), paste(length(unique(Drugfile$`Phecode Integer`))), paste("phenotypes. This score is constructed as a weighted sum of the effects of eight phenotypic specific features (described below) on drug indication using Firth logistic regression with five-fold cross-validation and was applied to"), comma(length(unique(allgenes$gene))),("protein-coding genes and"), length(unique(allgenes$parentterm)), ("phenotypes. We further incorporated the direction of genetic effect for each predictor using predictions of loss-of-function (LOF) and gain-of-function (GOF) from LoGoFunc"),HTML(paste0(tags$sup("1"))),span("for clinical variants and quantitative trait loci estimates for genome-wide association studies (GWAS) phenotypes and created a complementary genetic priority score with the direction of effect (GPS-D). In our study, we discovered that drugs with a high priority score are much more likely to progress through clinical trials and we recommend using a score cut-off > 1.5 that corresponds to an odds ratio (OR) >=9 of being a therapeutic target when compared to drugs with no genetic evidence.")),
        p(paste0("This web application can be used to search for putative therapeutic targets with genetic priority scores and evidence for direction of genetic effect. We included drug indication data from the Open Targets Platform"),HTML(paste0(tags$sup("2"))), ( "and the Side effect Resource (SIDER) 4.1"),HTML(paste0(tags$sup("3"))) ,(" and mapped all drug indications and genetic phenotypes to phecode integer terms. We have provided evidence of each genetic association and predicted direction of genetic effect and users can search this website by a gene target or phenotype of interest. We note that these results should be a starting point for users to conduct follow-up analyses and if a gene and phenotype are not present in this table, this should not be interpreted as a drug target with genetic support against it.")),
        p("For further details about the methods and analysis behind this study, please see our paper: Duffy, A et al. Development of a human genetics-guided priority score for 19,365 genes and 347 drug indications. Submitted.")
    ),

    h3("Instructions"),
    wellPanel(
        p("Users can use this website to explore associations by gene or by phenotype. All associations have evidence from at least one genetic association and this evidence can be refined using the score cutoff slider at 0.3 increment thresholds of the genetic priority score"),
        p(span("When the Gene tab is selected, users can search the priority scores for their gene of choice in the"),
          strong("‘Search Gene Target’"),
          span("tab. When the phenotype tab is selected, users can select their phenotype of interest from the dropdown"), strong("‘Select phenotype:’"),
          span("dropdown bar."),
          br()),
        p(span("Once a gene or phenotype is searched on the left side panel, this will generate a table of results in the", strong("Table"),span("tab. Users can click on the"), strong("Evidence"),span("and the"), strong("DOE Evidence"), span("tab to generate heatmaps of the associations with the genetic predictors on the X-axis and the associated gene/phenotypes on the Y-axis. The"), strong("Evidence"), span("tab generates a heatmap with the genetic evidence for a predictor colored in purple and the user can hover over the cell to get the phenotype description from the phenotype data source. The"),
         strong("DOE Evidence"), span("tab generates a heatmap with genetic evidence predicted as loss-of-function (LOF) colored in pink and predictions of gain-of-function (GOF) colored in blue."),
        span("For genes/phenotypes with more than 100 rows of evidence, only the top 100 high-scored genes/phenotypes are shown.")))
    ),
    h3("Genetic Evidence"),
    wellPanel(
        p(paste0("Using publicly available data sources, we collected eight genetic features from three types of genetic evidence (clinical variants, coding variants and GWAS phenotypes). We mapped these features to phecodes and restricted these to phecode integer terms. Phecodes that mapped to phecode categories ‘neoplasms’, ‘infectious diseases’, ‘pregnancy complications’, ‘injuries and poisonings’ and ‘null’ were excluded from the analysis. Each of these genetic features is described below.")),
        h4("Clinical variants"),
        p(span(strong("OMIM:"), span("Mendelian genes from the Online Mendelian Inheritance in Man (OMIM) database"),HTML(paste0(tags$sup("4"))),span('.'),
        )),
        p(span(strong("HGMD:"), span("Disease-causing and likely disease-causing genes for human inherited diseases from the human mutation database (HGMD) professional"),HTML(paste0(tags$sup("5"))),span('.'),
        )),
        p(span(strong("EVA-ClinVar:"), span("Clinically relevant genetic variants and diseases from ClinVar"),HTML(paste0(tags$sup("6"))), span(". This evidence is obtained from the Open Target platform and was curated by European Variation Archive (EVA)"),HTML(paste0(tags$sup("7"))), span("ClinVar evidence was filtered on clinical significance: likely pathogenic', 'association', 'confers sensitivity', 'drug response', 'protective' and 'pathogenic' and the confidence of the submission was filtered on: 'criteria provided, single submitter'; 'criteria provided, conflicting interpretations'; 'criteria provided, multiple submitters, no conflicts'; 'reviewed by expert panel' and 'practice guideline'."),
        )),
        h4("Coding variants"),

        p(span(strong("Gene burden:"), span("Exome-based gene burden association statistics from the UK Biobank hosted by Genebass"),HTML(paste0(tags$sup("8"))), span(". We extracted traits labeled ‘ICD first occurrence’ and restricted to missense and pLOF burden sets ("),HTML(paste0(tags$em("P"))),span(" < 4.3 x 10"),HTML(paste0(tags$sup("-7"))),span(")."),
        )),
             p(span(strong("Single-Variant:"), span("Exome-based single variant association statistics from the UK Biobank hosted by Genebass"),HTML(paste0(tags$sup("8"))), span(". We extracted traits labeled ‘ICD first occurrence’ and restricted to missense and pLOF variants ("),HTML(paste0(tags$em("P"))),span(" < 4.3 x 10"),HTML(paste0(tags$sup("-7"))),span(")."),
        )),
               h4("Genome-wide association study phenotypes"),
        
        p(span(strong("eQTL phenotype:"), span("Genes with a GWA phenotype driven by gene expression regulation through shared variants. These phenotypes were obtained from the Pan-UK Biobank"),HTML(paste0(tags$sup("9"))), span("and consisted of 484 ICD-10 traits, 729 phecodes and 39 continuous traits that we had manually mapped to phecodes from the UK Biobank Neale."),
        )),
        p(span(strong("Locus2gene:"), span(" Genes with a GWA phenotype identified by the Locus2Gene machine learning model from the Open Targets Genetics Portal which prioritizes likely causal genes at each GWAS locus by integrating fine-mapped associations with functional genomic features"),HTML(paste0(tags$sup("10"))),span('.'),
        )),
          p(span(strong("pQTL phenotype:"), span(" Genes with a GWA phenotype that are in high linkage disequilibrium (r"),HTML(paste0(tags$sup("2"))),span(">0.80) with at least one "),HTML(paste0(tags$em("cis "))),span("or"), HTML(paste0(tags$em("trans "))) , span("pQTL from supplementary table 12 from the large-scale plasma proteome study conducted by Ferkingstad et al"),HTML(paste0(tags$sup("11"))),span('.'),
        ))),
    h3("Citation"),
    wellPanel(
        p('Duffy, A et al. Development of a human genetics-guided priority score for 19,365 genes and 347 drug indications. Submitted.'),
    ),
    h3("Contact"),
    wellPanel(
        p('For any queries please contact:'),
        p(mailtoR(email = "Ron.do@mssm.edu",
                  text = "Ron.do@mssm.edu"),

          use_mailtoR()),
        tags$a(href='https://labs.icahn.mssm.edu/dolab/people/',
               'The Do Lab'),

    ),
    br(),
    h3("References"),
    wellPanel(
        p(span('1.  David, S. et al. Genome-wide prediction of pathogenic gain- and loss-of-function variants from ensemble learning of a diverse feature set. '), em('bioRxiv'),span('2022.06.08.495288 (2022).')),
        p(span('2.	Ochoa, D. et al. Open Targets Platform: supporting systematic drug-target identification and prioritisation. '), em('Nucleic Acids Res '), strong(' 49'), span(', D1302-D1310 (2021).')),
        p(span('3.	Kuhn, M., Letunic, I., Jensen, L.J. & Bork, P. The SIDER database of drugs and side effects. '), em('Nucleic Acids Res '), strong('44'), span(', D1075-D1079 (2016).')),
        p(span('4.	Hamosh, A. et al. Online Mendelian Inheritance in Man (OMIM), a knowledgebase of human genes and genetic disorders. '), em(' Nucleic Acids Res '), strong('30'), span(', 52-5 (2002).')),
        p(span('5.	Stenson, P.D. et al. The Human Gene Mutation Database (HGMD(®)): optimizing its use in a clinical diagnostic or research setting. '), em('Hum Genet '), strong('139'), span(', 1197-1207 (2020).')),
        p(span('6.	Landrum, M.J. et al. ClinVar: improving access to variant interpretations and supporting evidence. '), em('Nucleic Acids Res '), strong('46'), span(', D1062-D1067 (2017).')),
        p(span('7.	Cook, C.E. et al. The European Bioinformatics Institute in 2016: Data growth and integration. '), em('Nucleic Acids Res '), strong('44'), span(', D20-D26 (2015).')),
        p(span('8.	Karczewski, K.J. et al. Systematic single-variant and gene-based association testing of thousands of phenotypes in 394,841 UK Biobank exomes.'), em(' Cell Genomics '), strong('2'), span(', 100168 (2022).')),
         p('9.	Pan-UKB team https://pan.ukbb.broadinstitute.org. (2020).'),
        p(span('10.	Mountjoy, E. et al. An open approach to systematically prioritize causal variants and genes at all published human GWAS trait-associated loci.'), em(' Nat Genet '), strong('53'), span(', 1527-1533 (2021).')),
          p(span('11. Ferkingstad, E. et al. Large-scale integration of the plasma proteome with genetics and disease.'),em(' Nat Genet '),strong('53'), span( '1712-1721 (2021).'))

    ))
