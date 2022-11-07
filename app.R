
library(R.utils)
library(shiny)
library(data.table)
library(DT)
library(dplyr)
library(bslib)
library(packcircles)
library(plotly)
library(reshape2)
library(shinydashboard)
library(shinyWidgets)
library(mailtoR)
library(stringr)
library(scales)


####### load drug - genetic data #######
# file 1 - all phenotypes and genes
allgenes=fread('data/Allgenes_andparentterms.txt', data.table=F)
# file 2 - all phenotypes and genes with OR > 4 and drug, genetic data
Drugfile<-fread('data/Allgenes_or4_prioritytable_miphenotypes_withdescription_updated.txt.gz', data.table =F, quote="")

  ## QC drug data 
names(Drugfile)=c('Gene','Phecode_Integer','SIDER_Indication','SIDER_drug','Open_Targets_Indication','Open_Targets_Clinical_Phase',
                  'Open_Targets_drug','OMIM','Coloc_phenotype','Eqtl_phenotype','Coloc_tissue','HGMD','Genebass','Clinvar',
                  'Open_Targets_Genetics_Portal','Tissue_specific','Eqtl_tissue','OE','Genescore_sum','Number_of_predictors','Order','Percent',
                  'Category','Phecode_description','Clinvar_description','Coloc_phenotype_description','Coloc_tissue_description',
                  'Eqtl_phenotype_description','Eqtl_tissue_description','Genebass_description','HGMD_description','OMIM_description','Open_Targets_Genetics_Portal_description','OR_cutoff')

Drugfile[c(3:5,7,23:33)]=apply(Drugfile[c(3:5,7,23:33)],2, function(X) str_to_sentence(X))
Drugfile$Genescore_sum=round(Drugfile$Genescore_sum,3)
Drugfile$Percent=round(Drugfile$Percent,3)

druggene=Drugfile %>% group_by(Gene) %>% mutate( max=max(OR_cutoff)) %>% distinct(Gene,max) 
druggene1<-druggene %>% arrange((Gene))

genes_phenotype=as.data.frame(rbind(cbind(filteron=unique(Drugfile$Gene), from='Gene'), cbind(filteron=unique(Drugfile$Phecode_Integer),from='Phenotype')))
druggene_pheno<-Drugfile %>% group_by(Gene, Phecode_description) %>% mutate( max=max(OR_cutoff)) %>% distinct(Gene, Phecode_description, max) 


druggene_pheno3 <-druggene_pheno %>% mutate(Phecode_description='All') %>%  distinct()
druggene_pheno<-rbind(druggene_pheno3,druggene_pheno) %>%arrange(Gene)



# Define UI ----

about_page <- tabPanel(
    title = "Home",
    setBackgroundColor(
        color = c("#F5FFFA", "#E6E6FA"),
        gradient = "linear",#E6E6FA"
        direction = "bottom"
    ),

    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    fluidRow(
        column(width=3),
        column(width=6,
               
               div(style = "text-align:center;,
      padding-left: 4px; padding-right: 4px; padding-top: 4px;
      padding-bottom: 4px;padding:4px;
      dborder-color:black;background-color: white;",
                   h1("Genetic Priority Score for Therapeutic Targets", align='center'),
                   br(),
                   div(
                       style = "text-align:center;",
                       actionButton("Search_by_gene", "Search by Gene",class="btn btn-info"),
                       actionButton("Search_by_phenotype", "Search by Indication",class="btn btn-info"),
                       br(),
                       br()),
                   br(), width=8,align="center"),
        )),
    br(),
    br(),
    br(),
    h4(paste('A human genetics-guided priority tool for informing therapeutic target selection for',comma(length(unique((Drugfile$Gene)))), 'genes and', length(unique((Drugfile$Phecode_Integer))), 'phenotypes.'), align='center'),
    h4(paste('Last update:',Sys.Date()), align='center'),
    br(),
    br(),
    br(),
    br(),
    fluidRow(
        column(2),
        column(2),
        column(1),
        column(1,
               div(  align="center",  
                     tags$a(href='https://labs.icahn.mssm.edu/dolab/people/',
                            tags$img(src="Mountsinailogo.jpeg", height='70', width='70'),
                            ''))),
        column(1,
               tags$a(href='https://github.com/rondolab/',
                      tags$img(src="github.png", height='70', width='70')
               ))
    ),
    br(),
    br(),
    br())


## About page 
abouttitle='About'
Aboutinfo_page <- tabPanel(
    theme=shinythemes::shinytheme('flatly'),
    title='About',

    h3("Summary"),
    wellPanel(   
      
      p(paste("We created an in-silico genetic priority score that can inform drug target prioritization and validation. This score is constructed as a weighted sum of the effects of eleven genetic features (described below) on drug indication using Firth logistic regression with five-fold cross-validation and was applied to"), comma(length(unique(allgenes$gene))),("genes and"), length(unique(allgenes$parentterm)), ("phenotypes. In our study, we discovered that drugs with a high priority score are much more likely to progress through clinical trials, and we established score thresholds that corresponded to an odds ratio (OR) >=4, OR >=5, and OR >=6 of being a therapeutic target when compared to drugs with lower ranked priority scores.")),
      p(paste0("This web application can be used to search for putative therapeutic targets with a high genetic priority score corresponding to an OR >=4. We included drug indication data from the Open Targets Platform"),HTML(paste0(tags$sup("1"))), ( "and the Side effect Resource (SIDER) 4.1"),HTML(paste0(tags$sup("2"))) ,(" and mapped all drug indications and genetic phenotypes to phecode integer terms. We have provided evidence of each genetic association, and users can search this website by a gene target or phenotype of interest. We note that these results should be a starting point for users to conduct follow-up analyses, and if a gene and phenotype are not present in this table, this should not be interpreted as a drug target with genetic support against it.")),
      p("For further details about the methods and analysis behind this study, please see our recently published paper: insert ref. ")
    ),
    
    h3("Instructions"),
    wellPanel(   
      p("Users can use this website to explore associations by gene or by phenotype. All associations have a genetic priority score that corresponds to an OR >= 4. This OR cutoff can be refined to OR >= 5 or an OR >= 6 using the OR cutoff slider."),
      p(span("When the Gene tab is selected, users can search the priority scores for their gene of choice in the"),
        strong("‘Search Gene Target‘"),
        span("tab. The user can choose to show all phenotype results or select a phenotype of interest in the"), 
        strong("‘Choose Phenotypes:’"),
        span("dropdown bar. When the phenotype tab is selected, users can select their phenotype of interest from the dropdown"), strong("‘Select phenotype:’"), 
        span("dropdown bar and further restrict by the gene target in the the"), strong("‘Choose Gene:’"), span("dropdown bar."),
        br()),
      p(span("Once a gene or phenotype is searched on the left side panel, this will generate a table of results in the", strong("Table"),span("tab. Users can click on the"), strong("Evidence"), span("tab to generate a heatmap with the genetic predictors on the x-axis and the associated gene/phenotypes on the y-axis. Genetic evidence for a predictor is colored in purple, and the user can hover over the cell to get the phenotype description from the phenotype data source. For genes/phenotypes with more than 100 rows of evidence, only the top 100 high-scored genes/phenotypes are shown.")))
    ),
    h3("Genetic Evidence"),
    
    wellPanel(
      p(paste0("We collected genetic features from six types of genetic evidence (clinical variants, coding variants, genome-wide association studies (GWAS) phenotypes, tissue phenotypes, tissue specificity, and constraint) to investigate the association between human genetic variation in drug target genes and drug indications. We collected this data from several publicly available sources. To allow comparison between the drug and genetic datasets, we mapped all genetic phenotypes to phecodes and restricted these to phecode integer terms. All phecodes mapped to the phecode categories ‘neoplasms’, ‘infectious diseases’, ‘pregnancy complications’, ‘injuries and poisonings’, and ‘null’ were excluded from the analysis. We restricted each gene set to protein-coding genes, for which we obtained a list of "), comma(length(unique(allgenes$gene))), (" protein-coding genes from Ensembl"),HTML(paste0(tags$sup("3"))),(". Each of these genetic features is described below.")),	 
    h4("Clinical variants"),
              
              p(span(strong("OMIM:"), span("The Online Mendelian Inheritance in Man (OMIM)"),HTML(paste0(tags$sup("4"))), span(" database contains genes involved in rare mendelian diseases. We extracted 15,846 genes associated with 7,053 Mendelian traits, mapped to 7,216 Human phenotype ontology (HPO) terms from the OMIM database. We used a curated mapping file from the phenome-wide association studies (PheWAS) catalog"),HTML(paste0(tags$sup("5"))), span("to map the HPO terms to phecodes. In total, 5,792 phenotypes were mapped to 289 phecode integers associated with 4,147 genes."),
              )),
              p(span(strong("HGMD:"), span("The Human Gene Mutation Database (HGMD)"),HTML(paste0(tags$sup("6"))), span(" contains gene mutations underlying human inherited diseases, curated from the scientific literature. We extracted 12,795 disease-causing, and likely disease-causing genes associated with 16,810 phenotypes mapped to 1,722 Unified Medical Language System (UMLS) codes. We mapped the UMLS codes to phecodes via International Classification of Diseases (ICD) codes using the available UMLS Metathesaurus file and ICD10CM and ICD10 files from the PheWAS resource"),HTML(paste0(tags$sup("5"))), span(". There were 271 phecode integers associated with 11,642 genes."),
              )),
              p(span(strong("EVA-ClinVar:"), span("The Open Target Platform summarizes genetic evidence from multiple datasets containing associations between targets and disease, described by Experimental Factor Ontology (EFO) terms."),HTML(paste0(tags$sup("1"))), span(" One source is ClinVar"),HTML(paste0(tags$sup("7"))), span(", a curated database of associations between clinically relevant genetic variants, diseases, and the variant’s clinical significance. This data is processed and curated by European Variation Archive (EVA)"),HTML(paste0(tags$sup("8"))), span(" before being added to the Open Target platform. We extracted ClinVar evidence based on two steps. First, evidence was filtered on clinical significance terms: likely pathogenic', 'association', 'confers sensitivity', 'drug response', 'protective', and 'pathogenic'. Secondly, we filtered based on the confidence of the submission assigned as: 'criteria provided, single submitter', 'criteria provided, conflicting interpretations', 'criteria provided, multiple submitters, no conflicts', 'reviewed by expert panel', and 'practice guideline'. We mapped the EFO terms to phecodes using the previously mentioned disease mapping files from the Open Targets Platform, EBISPOT OLS, UMLS, and PheWAS resource. In total, there were 2,772 genes associated with 3,737 phenotypes that mapped to 222 phecode integers."),
              )),
              h4("Coding variants"),
              
              p(span(strong("Genebass:"), span("Gene-based Association Summary Statistics (Genebass)"),HTML(paste0(tags$sup("9"))), span(" is a resource of exome-based association statistics from the UK Biobank, encompassing 4,529 phenotypes across 426,370 individuals with gene-based and single-variant testing. We restricted to traits labeled ‘ICD first occurrence’ and extracted genome-wide significant gene associations annotated as either predicted loss of function (pLOF) or missense. We extracted 18,448 genes associated with 413 ICD codes that mapped to 234 integer phecodes."),
                     h4("Genome-wide association study phenotypes"),
              )),
              p(span(strong("Coloc phenotype:"), span("We defined coloc phenotype as genes with a GWA phenotype driven by gene expression regulation through colocalization, using a posterior probability (PPH4) > 0.8"),HTML(paste0(tags$sup("10,11"))), span(". We used the results from our coloc2 analysis from our previous study"),HTML(paste0(tags$sup("12"))), span(", which consists of putative causal genes and their corresponding phenotype and tissue from co-localization analysis, performed across 1,850 phenotypes from the UK Biobank"),HTML(paste0(tags$sup("13"))), span(" and expression quantitative trait loci (eQTL) data from 48 tissues from the Genotype-Tissue Expression (GTEx) v7 project"),HTML(paste0(tags$sup("14,15"))), span(" using the coloc2 method. The description of this approach has been previously described"),HTML(paste0(tags$sup("12"))), span(". These phenotypes consisted of 1,071 case-control phenotypes from phecodes"),HTML(paste0(tags$sup("16"))), span(", and 779 phenotypes from the UK Biobank Neale dataset (Round 1 - http://www.nealelab.is/uk-biobank). Of these Neale phenotypes, we mapped 112 ICD codes to phecodes and manually mapped 387 continuous traits to phecodes. There was a total of 15,990 coloc2 genes and 313 integer phecodes."),
              )),
              p(span(strong("eQTL phenotype:"), span("We further generated an eQTL phenotype feature, using the updated Pan-UK Biobank"),HTML(paste0(tags$sup("17"))), span(", as genes with a GWA phenotype driven by gene expression regulation through shared variants. The Pan-UK Biobank is a multi-ancestry genetic analysis of the UK Biobank across six continental ancestry groups and 7,228 phenotypes. All variants were obtained on hg19 coordinates, and we used LiftOver"),HTML(paste0(tags$sup("18"))), span(" to convert these to hg38 for compatibility and extracted all genome-wide significant single nucleotide variants ("), (span(em("P"))),span(" < 5 x 10-8). We restricted the trait types and obtained summary statistics for 1,252 phenotypes: 484 ICD-10 traits, 729 phecodes, and 39 continuous traits that we had manually mapped to phecodes from the UK Biobank Neale. To identify the putative causal genes to these traits, we used eQTLs from GTEX v8"),HTML(paste0(tags$sup("14"))), span(" across 49 tissues and extracted significant variants for each tissue with a nominal P less than the gene level threshold. We intersected these variants with the Pan-UK Biobank significant variants to obtain the linked gene for each trait. We mapped all ICD traits to phecodes, and across the 49 tissues, there were 4,446 genes associated with 184 integer phecodes."),
              )),
              p(span(strong("Locus2gene:"), span("The Locus2Gene model is a machine learning method from the Open Targets Genetics Portal which prioritizes likely causal genes at each GWAS locus by integrating fine-mapped associations with functional genomic features"),HTML(paste0(tags$sup("19"))), span(". A set of gold standard positive (GSP) genes are used in this method where the causal gene has been confidently assigned. Many of these GSP genes are selected based on known drug–target disease pairs, so we removed these from the final output. To define the Locus2Gene predictor, we restricted to genes prioritized by the model with a score greater than 0.5. We mapped the disease traits to phecodes using the disease mapping files from the Open Targets Platform, EBISPOT OLS, UMLS, and PheWAS resources. We extracted 3,709 genes classified as casual for 1,730 traits mapped to 270 phecode integers."),
              )),
              h4("Tissue phenotypes"),
              p(span(strong("Coloc2 and eQTL tissue:"), span("The coloc2 and eQTL tissue features define genes with gene expression regulation in a given tissue that underlies a significant genome-wide association through a shared GWA and eQTL signal. To link the coloc2 and eQTL tissues to phecodes, we used LD score regression applied to specifically expressed genes (LDSC-SEG)"),HTML(paste0(tags$sup("20"))), span(" across 776 traits from the Neale GWAS resource and 1,112 traits from the Saige GWAS resource, using GWAS summary statistics that were mapped to phecodes. We used LDSC-SEG to estimate the heritability enrichment of specifically expressed genes across 54 GTEX V8 tissues to identify phenotypic relevant tissues. In brief, this method first classifies the top 10% of ‘specifically expressed genes’ in each tissue using a t-statistic test and extends each transcribed region by a 100kb window to create tissue-specific annotations. Stratified LD score regression is then applied to the GWAS summary statistics to test the contribution of the tissue-specific annotations to trait heritability, conditional on the baseline model, which includes 52 annotations including genic, enhancer, and conserved regions, and all genes. We used a P threshold of 0.01 to identify phenotype-relevant tissues, for which there were 421 phecodes across 54 tissues (mean =15.31 phecodes per tissue, sd=6.40). These phecodes were used to represent the phenotype-relevant tissues from our coloc2 and eQTL results."),
              )),
              h4("Tissue specificity"),
              
              p(span(strong("Tissue specificity of gene expression:"), span("We constructed a tissue specificity score of gene expression for each gene-phenotype. We used transcripts per million (TPM), a measure of gene expression per tissue, and Tau, a gene tissue specificity metric. We obtained median TPM values for each gene across all tissue samples from GTEX V8. We used LDSC-SEG to map the TPM values per tissue to phecodes and took the median TPM when multiple tissues were mapped. We dichotomized levels of median TPM and used a TPM threshold > 0.5 to indicate the presence of gene expression. We calculated Tau for each gene using gene counts from GTEX V8, as described by Kryuchkova Mostacci and Robinson Rechavi"),HTML(paste0(tags$sup("21"))), span(". We used a threshold of > 0.8 to dichotomize Tau and indicate tissue specificity. Using these dichotomized values, we created a genetic feature for tissue specificity of gene expression for each gene-phecode with a TPM > 0.5 and a Tau > 0.8. We identified 6,344 genes as tissue-specific across 207 phecode integers (mean=86.42 phecode integers per gene, sd=74.21)."), 
              )),
              h4("Constraint"),
              
              p(span(strong("OE:"), span("We used the observed/expected (OE) score to measure the constraint per gene and its intolerance to loss of function mutations. We used the suggested cutoff of < 0.35 of the upper bound of the OE confidence interval to identify genes under high constraint"),HTML(paste0(tags$sup("22"))), span(" and identified 2,854 genes.")
              ))), 

    h3("Citation"),
    wellPanel( 
        br(),
        br()),
    
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
      p('1.	Ochoa, D. et al. Open Targets Platform: supporting systematic drug-target identification and prioritisation. Nucleic Acids Res 49, D1302-D1310 (2021).'),
      p('2.	Kuhn, M., Letunic, I., Jensen, L.J. & Bork, P. The SIDER database of drugs and side effects. Nucleic acids research 44, D1075-D1079 (2016).'),
      p('3.	Cunningham, F. et al. Ensembl 2022. Nucleic Acids Research 50, D988-D995 (2021).'),
      p('4.	Hamosh, A. et al. Online Mendelian Inheritance in Man (OMIM), a knowledgebase of human genes and genetic disorders. Nucleic Acids Res 30, 52-5 (2002).'),
      p('5.	Denny, J.C. et al. Systematic comparison of phenome-wide association study of electronic medical record data and genome-wide association study data. Nature Biotechnology 31, 1102-1111 (2013).'),
      p('6.	Stenson, P.D. et al. The Human Gene Mutation Database (HGMD(®)): optimizing its use in a clinical diagnostic or research setting. Hum Genet 139, 1197-1207 (2020).'),
      p('7.	Landrum, M.J. et al. ClinVar: improving access to variant interpretations and supporting evidence. Nucleic Acids Research 46, D1062-D1067 (2017).'),
      p('8.	Cook, C.E. et al. The European Bioinformatics Institute in 2016: Data growth and integration. Nucleic Acids Research 44, D20-D26 (2015).'),
      p('9.	Karczewski, K.J. et al. Systematic single-variant and gene-based association testing of thousands of phenotypes in 394,841 UK Biobank exomes. Cell Genomics 2, 100168 (2022).'),
      p('10.	Dobbyn, A. et al. Landscape of Conditional eQTL in Dorsolateral Prefrontal Cortex and Co-localization with Schizophrenia GWAS. Am J Hum Genet 102, 1169-1184 (2018).'),
      p('11.	Giambartolomei, C. et al. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. PLoS Genet 10, e1004383 (2014).'),
      p('12.	Duffy, A. et al. Tissue-specific genetic features inform prediction of drug side effects in clinical trials. Sci Adv 6(2020).'),
      p('13.	Bycroft, C. et al. The UK Biobank resource with deep phenotyping and genomic data. Nature 562, 203-209 (2018).'),
      p('14.	Aguet, F. et al. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 369, 1318-1330 (2020).'),
      p('15.	GTEx Portal https://www.gtexportal.org/home/datasets.'),
      p('16.	Zhou, W. et al. Efficiently controlling for case-control imbalance and sample relatedness in large-scale genetic association studies. Nat Genet 50, 1335-1341 (2018).'),
      p('17.	Pan-UKB team https://pan.ukbb.broadinstitute.org. (2020).'),
      p('18.	Kuhn, R.M., Haussler, D. & Kent, W.J. The UCSC genome browser and associated tools. Briefings in Bioinformatics 14, 144-161 (2012).'),
      p('19.	Mountjoy, E. et al. An open approach to systematically prioritize causal variants and genes at all published human GWAS trait-associated loci. Nat Genet 53, 1527-1533 (2021).'),
   p('20.	Finucane, H.K. et al. Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types. Nature Genetics 50, 621-629 (2018).'),
   p('21.	Kryuchkova-Mostacci, N. & Robinson-Rechavi, M. A benchmark of gene expression tissue-specificity metrics. Brief Bioinform 18, 205-214 (2017).'),
   p('22.	Karczewski, K.J. et al. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature 581, 434-443 (2020).')
      )
)



Gene_page <- tabPanel(

    title='Gene',value="Gene",
    titlePanel("Genetic Priority Score"),
    theme=shinythemes::shinytheme('flatly'),
    
    sidebarLayout(
        sidebarPanel(

            # Add OR cutoff slider             
            sliderInput('ORselection', 'Use OR cutoff', min=4, max=6, value=5, step=1),
             # Select Gene and then phenotypes          
            searchInput("selectgene","Search Gene Target:",placeholder = 'Search Gene',  btnSearch = shiny::icon("search"),btnReset = shiny::icon("remove"),width = "100%"),
            helpText("Gene examples: ABCC8, JAK2, PCSK9"),
            selectizeInput("selectphenotype","Choose Phenotypes:", choices=c(unique(druggene_pheno$Phecode_description)),selected=NULL, options = list(create = TRUE, placeholder = 'Phenotypes'))
        ),
        mainPanel(
            tabsetPanel(
                tabPanel('Table',
                         
                         DT::dataTableOutput("table")),
                tabPanel('Evidence', 
                         plotly::plotlyOutput('plot_evidence'))
                
            )                                    
            
        ))
)

Phenotype_page <- tabPanel(

    title='Phenotype', value="Phenotype",
    titlePanel("Genetic Priority Score"),

    sidebarLayout(
        sidebarPanel(
            # Add OR cutoff slider             
            sliderInput('ORselection_phenotype', 'Use OR cutoff', min=4, max=6, value=5, step=1),
            selectizeInput("choosephenotype","Select Phenotype:", choices=NULL,  multiple = FALSE,  selected=character(0),options = list(create = TRUE, placeholder = 'Select Phenotype')),
            selectizeInput("choosegene","Choose Genes:", choices=NULL,selected=character(0), options = list(create = TRUE, placeholder = 'Genes'))
        ),
        
        mainPanel(
            tabsetPanel(
                tabPanel('Table',
                         DT::dataTableOutput("table_pheno")),
                tabPanel('Evidence', 
                         plotly::plotlyOutput('plot_evidence_phenotype'))
                
            ))            
        
    ))


ui <- fluidPage(
    navbarPage(
        id='tabset',
        title = 
            div(
                div(
                    img(src="druggene.png",
                        height = "100px",width = "100px",style = "position: relative; margin:-35px -35px; display:right-align;")),''),
        position = 'fixed-top',
        windowTitle="Genetic Priority Score",
        header=tags$style(type="text/css", "body {padding-top: 70px;},
  
  "),
        
        theme=shinythemes::shinytheme('cosmo'),
        about_page,
        Gene_page, 
        Phenotype_page,
        Aboutinfo_page
    ))

server <- function(input, output,session) {
    

    observeEvent(input$Search_by_gene, {
        updateTabsetPanel(session, "tabset" , selected="Gene")
    })
    
    observeEvent(input$Search_by_phenotype, {
        updateTabsetPanel(session, "tabset" , selected="Phenotype")
    })
    
    observe({
        listchoic_pheno=unique(druggene_pheno$Phecode_description[druggene_pheno$Gene==input$selectgene & druggene_pheno$max>=input$ORselection ])
        updateSelectizeInput(session, "selectphenotype", choices=listchoic_pheno,selected = NULL) 
        
    })
    
    observe({
        or_gene_cutoff_pheno=druggene_pheno %>% filter(max>=input$ORselection_phenotype & Phecode_description!='All') %>% ungroup() %>% distinct(Phecode_description)%>% arrange(Phecode_description)
        or_gene_cutoff_pheno1=or_gene_cutoff_pheno$Phecode_description
        updateSelectizeInput(session, "choosephenotype", choices=or_gene_cutoff_pheno1, selected=character(0)) 
    })
    
    
    observe({
        listchoic_pheno=druggene_pheno$Gene[druggene_pheno$Phecode_description==input$choosephenotype & druggene_pheno$max>=input$ORselection_phenotype]
        listchoic_pheno1=unique(c('All',listchoic_pheno))
        updateSelectizeInput(session, "choosegene", choices=listchoic_pheno1,selected=character(0)) 
        
    })
    
    
    datafile_react<- reactive({
        if ('All' %in% input$selectphenotype) {
            Drugfile1=Drugfile %>% filter(OR_cutoff >=input$ORselection)   %>% 
                filter(Gene %in% input$selectgene) %>% select(Gene, Phecode_Integer,Phecode_description,Open_Targets_Indication, Open_Targets_Clinical_Phase,SIDER_Indication,Category, Number_of_predictors,Genescore_sum,OR_cutoff,Percent) %>% distinct() %>% arrange(desc(Genescore_sum))
            colnames(Drugfile1)=gsub('_',' ', colnames(Drugfile1))
            return(Drugfile1)
            
        } else {
            
            Drugfile1=Drugfile %>%filter(OR_cutoff >=input$ORselection) %>%
                filter(Gene %in% input$selectgene)  %>%  filter(Phecode_description %in% input$selectphenotype)  %>% select(Gene, Phecode_Integer,Phecode_description,Open_Targets_Indication, Open_Targets_Clinical_Phase,SIDER_Indication,Category, Number_of_predictors,Genescore_sum,OR_cutoff, Percent)  %>% distinct()%>% arrange(desc(Genescore_sum))
            colnames(Drugfile1)=gsub('_',' ', colnames(Drugfile1))
            return(Drugfile1)
            
        }
        
    })
    
    output$table <- DT::renderDataTable({ 
        datafile_react()
    })
    
    
    datafile_react_pheno<- reactive({
        if ('All' %in% input$choosegene) {
            Drugfile1=Drugfile %>% filter(OR_cutoff >=input$ORselection_phenotype) %>%
                filter(Phecode_description %in% input$choosephenotype)%>% select( Phecode_Integer,Phecode_description,Gene,Open_Targets_Indication, Open_Targets_Clinical_Phase,SIDER_Indication,Category, Number_of_predictors,Genescore_sum,OR_cutoff,Percent) %>% distinct()%>% arrange(desc(Genescore_sum))
            colnames(Drugfile1)=gsub('_',' ', colnames(Drugfile1))
            return(Drugfile1)
        } else {
            
            Drugfile1=Drugfile %>%filter(OR_cutoff >=input$ORselection_phenotype) %>%
                filter(Gene %in% input$choosegene)  %>%  filter(Phecode_description %in% input$choosephenotype) %>% select( Phecode_Integer,Phecode_description,Gene,Open_Targets_Indication, Open_Targets_Clinical_Phase,SIDER_Indication,Category, Number_of_predictors,Genescore_sum,OR_cutoff,Percent) %>% distinct()%>% arrange(desc(Genescore_sum))
            colnames(Drugfile1)=gsub('_',' ', colnames(Drugfile1))
            return(Drugfile1)
        }
        
    })
    
    output$table_pheno <- DT::renderDataTable({ 
        datafile_react_pheno()
    })
    
    ##Evidence tables
    
    datafile_react_evidence<- reactive({
        
        if ('All' %in% input$selectphenotype) {
            
            Drugfile1=Drugfile %>%filter(OR_cutoff >=input$ORselection) %>%
                filter(Gene %in% input$selectgene)%>% arrange(desc(Genescore_sum)) %>% 
                select(Phecode_description,OMIM,Coloc_phenotype,Eqtl_phenotype,Coloc_tissue,
                       HGMD,Genebass,
                       Clinvar,Open_Targets_Genetics_Portal,Tissue_specific,Eqtl_tissue, OE)
            
            Drugfile2<-reshape2::melt(Drugfile1, id.vars=c('Phecode_description'))
            
            labels=Drugfile %>%filter(OR_cutoff >=input$ORselection) %>%
                filter(Gene %in% input$selectgene) %>% arrange(desc(Genescore_sum))%>% 
                select(Phecode_description,OMIM_description,Coloc_phenotype_description,Eqtl_phenotype_description,Coloc_tissue_description,
                       HGMD_description,Genebass_description,Clinvar_description,
                       Open_Targets_Genetics_Portal_description,Tissue_specific, Eqtl_tissue_description,OE)
            
            colnames(labels)=colnames(Drugfile1)
            labels2<-reshape2::melt(labels, id.vars=c('Phecode_description'),value.name = 'label')
            Drugfile2_labels=inner_join(Drugfile2,labels2, by = c("Phecode_description", "variable") )
            Drugfile2_labels$value<-as.character(Drugfile2_labels$value)
            Drugfile2_labels=Drugfile2_labels%>% distinct() 
            Drugfile2_labels$variable=gsub('_',' ', Drugfile2_labels$variable)
            Drugfile2_labels$label=gsub(";", "; ",     Drugfile2_labels$label)
            
            Drugfile2_labels$label=gsub("([a-z0-9]* [a-z0-9]* [a-z0-9]*) ", "\\1\n",     Drugfile2_labels$label)
            
            return(Drugfile2_labels)
        } else { 
            
            Drugfile1=Drugfile %>%filter(OR_cutoff >=input$ORselection) %>%
                filter(Gene %in% input$selectgene)  %>%  filter(Phecode_description %in% input$selectphenotype) %>% arrange(desc(Genescore_sum))%>% 
                select(Phecode_description,OMIM,Coloc_phenotype,Eqtl_phenotype,Coloc_tissue,
                       HGMD,Genebass,
                       Clinvar,Open_Targets_Genetics_Portal,Tissue_specific,Eqtl_tissue,OE)
            
            Drugfile2<-reshape2::melt(Drugfile1, id.vars=c('Phecode_description'))
            
            labels=Drugfile %>%filter(OR_cutoff >=input$ORselection) %>%
                filter(Gene %in% input$selectgene) %>% arrange(desc(Genescore_sum)) %>% 
                select(Phecode_description,OMIM_description,Coloc_phenotype_description,Eqtl_phenotype_description,Coloc_tissue_description,
                       HGMD_description,Genebass_description,Clinvar_description,
                       Open_Targets_Genetics_Portal_description,Tissue_specific, Eqtl_tissue_description,OE)
            
            colnames(labels)=colnames(Drugfile1)
            labels2<-reshape2::melt(labels, id.vars=c('Phecode_description'),value.name = 'label')
            Drugfile2_labels=inner_join(Drugfile2,labels2, by = c("Phecode_description", "variable") )
            #      Drugfile2_labels$Phecode_description<- factor(Drugfile2_labels$Phecode_description, levels=unique(Drugfile2_labels$Phecode_description))
            Drugfile2_labels$value<-as.character(Drugfile2_labels$value)
            Drugfile2_labels=Drugfile2_labels%>% distinct() 
            Drugfile2_labels$variable=gsub('_',' ', Drugfile2_labels$variable)
            Drugfile2_labels$label=gsub(";", "; ",     Drugfile2_labels$label)
            
            Drugfile2_labels$label=gsub("([a-z0-9]* [a-z0-9]* [a-z0-9]*) ", "\\1\n",     Drugfile2_labels$label)
            
            return(Drugfile2_labels)
            
        }
        
    })  
    
    
    output$plot_evidence <- plotly::renderPlotly({
        evid=datafile_react_evidence() %>%   
            ggplot(aes(x=variable, y=`Phecode_description`, fill= value, text=label)) + 
            geom_tile(colour = "white", lwd = 1.5,
                      linetype = 1) + theme(axis.text.x = element_text(angle = 45), axis.title.x=element_blank(),
                                            axis.title.y=element_blank()) +
            #scale_y_discrete(limits = rev(levels(Drugfile2_labels$Phecode_description))) +
            scale_fill_manual(values=c("azure3","blueviolet")) +labs(fill='Evidence')
        
        
        heightplot=length(unique(evid$data$Phecode_Integer))
        
        if(heightplot<20)  {
            ggplotly(evid, tooltip="text", height=800) 
        }
        else{
            
            ggplotly(evid, tooltip="text", height=2000) }
    })
    
    
    # Evidence phenotype
    
    ##Evidence tables
    
    datafile_react_evidence_phenotype<- reactive({
        
        
        if ('All' %in% input$choosegene) {
            
            Drugfile1=Drugfile %>%filter(OR_cutoff >=input$ORselection_phenotype) %>%
                filter(Phecode_description %in%  input$choosephenotype)%>% 
                select(Gene,Genescore_sum,OMIM,Coloc_phenotype,Eqtl_phenotype,Coloc_tissue,
                       HGMD,Genebass,
                       Clinvar,Open_Targets_Genetics_Portal,Tissue_specific,Eqtl_tissue,OE)
            
            Drugfile2<-reshape2::melt(Drugfile1, id.vars=c('Gene','Genescore_sum'))
            
            labels=Drugfile %>%filter(OR_cutoff >=input$ORselection_phenotype) %>%
                filter(Phecode_description %in%  input$choosephenotype) %>% 
                select(Gene,Genescore_sum,OMIM_description,Coloc_phenotype_description,Eqtl_phenotype_description,Coloc_tissue_description,
                       HGMD_description,Genebass_description,Clinvar_description,
                       Open_Targets_Genetics_Portal_description,Tissue_specific, Eqtl_tissue_description,OE)
            
            
            colnames(labels)=colnames(Drugfile1)
            labels2<-reshape2::melt(labels, id.vars=c('Gene','Genescore_sum'),value.name = 'label')
            Drugfile2_labels=inner_join(Drugfile2,labels2, by = c("Gene", "Genescore_sum", "variable") )

            Drugfile2_labels=Drugfile2_labels %>% distinct()
            Drugfile2_labels$variable=gsub('_',' ', Drugfile2_labels$variable)
            Drugfile2_labels$label=gsub(";", "; ",     Drugfile2_labels$label)
            
            Drugfile2_labels$label=gsub("([a-z0-9]* [a-z0-9]* [a-z0-9]*) ", "\\1\n",     Drugfile2_labels$label)
            if( length(unique(Drugfile2_labels$Gene)) >100) {
                
                topgenepheno=head(Drugfile2_labels  %>% distinct(Gene,Genescore_sum),100)
                
                Drugfile2_labels_topgenes=inner_join(Drugfile2_labels,topgenepheno)
                Drugfile2_labels_topgenes=Drugfile2_labels_topgenes %>% select(-Genescore_sum)
                
                Drugfile2_labels_topgenes$Gene<- factor(Drugfile2_labels_topgenes$Gene, levels=unique(Drugfile2_labels_topgenes$Gene))
                Drugfile2_labels_topgenes$value<-as.character(Drugfile2_labels_topgenes$value)
                
            }else {
                
                Drugfile2_labels_topgenes=Drugfile2_labels %>% select(-Genescore_sum)
                
                Drugfile2_labels_topgenes$Gene<- factor(Drugfile2_labels_topgenes$Gene, levels=unique(Drugfile2_labels_topgenes$Gene))
                Drugfile2_labels_topgenes$value<-as.character(Drugfile2_labels_topgenes$value)
                
            }
            return(Drugfile2_labels_topgenes)
            
        } else {
            
            Drugfile1=Drugfile %>%filter(OR_cutoff >=input$ORselection_phenotype) %>%  filter(Gene %in% input$choosegene)   %>% 
                filter(Phecode_description %in% input$choosephenotype)%>% 
                select(Gene,OMIM,Coloc_phenotype,Eqtl_phenotype,Coloc_tissue,
                       HGMD,Genebass,
                       Clinvar,Open_Targets_Genetics_Portal,Tissue_specific,Eqtl_tissue,OE)
            
            Drugfile2<-reshape2::melt(Drugfile1, id.vars=c('Gene'))
            
            labels=Drugfile %>%filter(OR_cutoff >=input$ORselection_phenotype) %>%
                filter(Phecode_description %in% input$choosephenotype) %>% 
                select(Gene,OMIM_description,Coloc_phenotype_description,Eqtl_phenotype_description,Coloc_tissue_description,
                       HGMD_description,Genebass_description,Clinvar_description,
                       Open_Targets_Genetics_Portal_description,Tissue_specific, Eqtl_tissue_description,OE)
            
            colnames(labels)=colnames(Drugfile1)
            labels2<-reshape2::melt(labels, id.vars=c('Gene'),value.name = 'label')
            Drugfile2_labels=inner_join(Drugfile2,labels2, by = c("Gene", "variable") )
            #      Drugfile2_labels<-factor(Drugfile2_labels)
            Drugfile2_labels$Gene<- factor(Drugfile2_labels$Gene, levels=unique(Drugfile2_labels$Gene))
            Drugfile2_labels$value<-as.character(Drugfile2_labels$value)
            Drugfile2_labels=Drugfile2_labels %>% distinct()
            Drugfile2_labels$variable=gsub('_',' ', Drugfile2_labels$variable)
            Drugfile2_labels$label=gsub(";", "; ",     Drugfile2_labels$label)
            
            Drugfile2_labels$label=gsub("([a-z0-9]* [a-z0-9]* [a-z0-9]*) ", "\\1\n",     Drugfile2_labels$label)
            Drugfile2_labels_topgenes=Drugfile2_labels
            return(Drugfile2_labels_topgenes)
            
        }
    })  
    
    
    output$plot_evidence_phenotype <- plotly::renderPlotly({
        evid=datafile_react_evidence_phenotype() %>%  
            ggplot(aes(variable, Gene, fill= value, text=label)) + 
            geom_tile(colour = "white", lwd = 1.5,
                      linetype = 1) + theme(axis.text.x = element_text(angle = 45), axis.title.x=element_blank(),
                                            axis.title.y=element_blank()) +
            #  scale_y_discrete(limits = rev(levels(Drugfile2_labels_topgenes$Gene))) +
            scale_fill_manual(values=c("azure3","blueviolet")) +labs(fill='Evidence')
        heightplot=length(unique(evid$data$Phecode_Integer))
        if(heightplot<10)  {
            ggplotly(evid, tooltip="text") 
        }
        if(heightplot>=10 & heightplot<20)  {
            ggplotly(evid, tooltip="text",height=1500) 
        }else{
            
            
            ggplotly(evid, tooltip="text", height=2000) }
    })
    
    
}

shinyApp(ui, server)
