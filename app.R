## Tabs for genes/phenotypes
#filter on Gene 
# filter on phenotype
# use OR as slider to choose cutoff 

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

#setwd(dir = 'Documents/PhD/DrugGenetics/Dataset/Updated_dataset/Final_MI/Sept14/')
#library(ggplot2)

#1 Make Navbar and put tab for About and Score
#https://shiny.rstudio.com/gallery/navbar-example.html
#2. add link to rondo lab page 
#3/ change colors / change theme - https://splunktool.com/background-color-in-tabsetpanel-in-shiny
#https://www.w3schools.com/colors/colors_hex.asp
# color https://www.colourlovers.com/palette/3527826/Bootswatch_-Flatly
###2  merge tabs together for Gene and pheno 



## make tab for pheno and tab for Gene
# on Gene tab select Gene  and then show all phenotypes for 
#Drugfile<-fread('/Users/aineduffy/Documents/PhD/DrugGenetics/Dataset/Updated_dataset/Final_MI/Sept14/Allgenes_or4_prioritytable_miphenotypes_withdescription.txt.gz', data.table =F, quote="")

Drugfile<-fread('data/Allgenes_or4_prioritytable_miphenotypes_withdescription_sampled.txt', data.table =F, quote="")
names(Drugfile)=c('Gene','Phecode_Integer','Sider_Indication','Sider_drug','Open_Targets_Indication','Open_Targets_Clinical_Phase',
                  'Open_Targets_drug','OMIM','Coloc_phenotype','Eqtl_phenotype','Coloc_tissue','HGMD','Genebass','Clinvar',
                  'Open_Targets_Genetics_Portal','Tissue_specific','Eqtl_tissue','OE','Genescore_sum','Number_of_predictors','Order','Percent',
                  'Category','Phecode_description','Clinvar_description','Coloc_phenotype_description','Coloc_tissue_description',
                  'Eqtl_phenotype_description','Eqtl_tissue_description','Genebass_description','HGMD_description','OMIM_description','Open_Targets_Genetics_Portal_description','Druggable_gene','OR_cutoff','therapeutic')
#Drugfile<-Drugfile[c(1,2,24,34,36,5,6,3,23,20,19,35,22,8:18,4,7,25:33)] %>% distinct()

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
  #theme=shinythemes::shinytheme('flatly'),
 # shinythemes::themeSelector(),  # <--- Add this somewhere in the UI
  
  #   color = c("#F7FBFF", "#2171B5"),
  
  setBackgroundColor(
    color = c("#F5FFFA", "#E6E6FA"),
    
    gradient = "linear",#E6E6FA"
    direction = "bottom"
  ),
  
  #fillPage(
  
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  br(),
  #box(
  fluidRow(
    column(width=3),
    column(width=6,
  
  div(style = "text-align:center;

,
      padding-left: 4px; padding-right: 4px; padding-top: 4px;
      padding-bottom: 4px;padding:4px;
      dborder-color:black;background-color: white;",
          
          h1("Human genetics-driven priority index for drug indications", align='center'),
     
      br(),
      
          div(
            style = "text-align:center;",
            actionButton("Search_by_gene", "Search by Gene",class="btn btn-info"),
            
            actionButton("Search_by_phenotype", "Search by Indication",class="btn btn-info"),
            br(),
            br())
      ,

      br(), width=8,align="center"),
    )),
 br(),
 br(),
 br(),
 h4(paste('Prioritization tool for informing putative therepeutic targets for ',length(unique((Drugfile$Gene))), 'genes and ', length(unique((Drugfile$Phecode_Integer))), 'phenotypes.'), align='center'),
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
 
  #hr(), #add horizontal line
  # ),
  #fluidRow(
#  hr()


  #box(
   # div(style = "text-align:center; border-style: solid;
    #  #padding-left: 4px; padding-right: 4px; padding-top: 4px; padding-bottom: 4px;
     # dborder-color:black;background-color:#2C3E50 ;",class = "footer",
        
   #     actionButton("about", "About",class="btn-info"),
    #    tags$a(href='https://labs.icahn.mssm.edu/dolab/people/',
     #          tags$img(src="Mountsinailogo.jpeg", height='50', width='50'),
      #         ' Do Lab')),width = 12)
  


abouttitle='About'
Aboutinfo_page <- tabPanel(
  theme=shinythemes::shinytheme('flatly'),
  
  title='About',
#  fillPage(
    
    h3("Summary"),
    wellPanel(   
      p("This web application can be used to search for putative drug gene targets that have a PI sum that corresponds to an OR >4. We have provided evidence of genetic association and note that this should be a starting point for users to then conduct follow up analyses.")
    ),
      h3("Instructions"),
    wellPanel(   
     p("Users can use this website to explore associations by gene or by phenotype. All associations have a gene score sum which corresponds to an OR >= 4. This OR cut off can be refined to OR >= 5 or an OR >=6 using the OR cutoff slider."),
    p(span("When the Gene tab is selected users can search the PI for their gene of choice in the"), 
    strong("‘Search Gene Target’"),span("tab.  The user can choose to show all phenotype results or select a phenotype of interest in the"),
    strong("‘Choose Phenotypes:’"), span("dropdown bar.  When the phenotype tab is selected, users can select their phenotype of interest from the drop down"),
    strong("‘Select phenotype:’"), span("dropdown bar and these results can be further restricted by gene in the "),strong("‘Choose Gene:’ "), span('dropdown bar.'),
    br()),
    p(span("Once a gene or phenotype is searched on the left side panel, 
           this will generate a table of results in the", strong("Table"), span("tab. Users can click on the"),strong("Evidence"), span(" tab to generate a heatmap with the genetic predictors on the X axis and the associated gene/phenotypes on the Y axis. Genetic evidence for a predictor is colored in purple and the user can hover over the cell to get the phenotype description from the phenotype data source."))),
    br()),

    h3("Methods"),
    wellPanel(   ),

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
           'Do Lab'),
    

      

     ),
    
  
    br()
#  )
)
  

#  br(),
 # br(),
  #br(),
 # title = abouttitle
  #div(
   # img(src="Mountsinailogo.jpeg", height=100, width=100,
    #    style="margin:10px 10px"),
    #abouttitle
      
    #)
  #)
 
# titlePanel("Priority Index Score"),
 # tags$a(href='https://labs.icahn.mssm.edu/dolab/people/',
  #     #  tags$img(src="Mountsinailogo.jpeg", height='50', width='50'),
   #      '? Do Lab'),

  
#)
 #) 
  


Gene_page <- tabPanel(
  #shinythemes::themeSelector(),  # <--- Add this somewhere in the UI
  
  title='Gene',value="Gene",
  titlePanel("Priority Index Score"),
  theme=shinythemes::shinytheme('flatly'),
  
  
  #  tabsetPanel(id='Tab',
  #             tabPanel('Search by Gene', value='Gene',
  sidebarLayout(
    sidebarPanel(
      #tags$style(".well{background-color: #ABBAEA};"),
      
      # Add OR cutoff slider             
      sliderInput('ORselection', 'Use OR cutoff', min=4, max=6, value=5, step=1),
      # Add choices for DG and therpeutic target          
      # radioButtons('drug_genome',label='Druggable Gene', choices=c("Yes"="DRUGGABLE GENOME", "No"="Not druggable","Select either"="Either"),  selected="DRUGGABLE GENOME"),
      #radioButtons('therapeutic',label='Therapeutic Target', choices=c("Yes"=1, "No"=0, "Select either"="Either"),  selected=1),
      # Select Gene and then phenotypes          
      searchInput("selectgene","Search Gene Target:",placeholder = 'Search Gene',  btnSearch = icon("magnifying-glass"),btnReset = icon("xmark"),width = "100%"),
      
      #      selectizeInput("selectgene","Search Gene Target:", selected=NULL,choices=NULL,   options = list(create = FALSE, placeholder = 'Select Gene', maxOptions = 3)),
      #   selectizeInput("selectgene","Search Gene Target:", choices=NULL,   options = list(create = FALSE, placeholder = 'Search Gene',maxOptions = 1, 
      #      onDropdownOpen = I("function($dropdown) {if (!this.lastQuery.length) {this.close(); this.settings.openOnFocus = false;}}"),
      #     onType = I("function (str) {if (str === \"\") {this.close();}}"))),
      
      selectizeInput("selectphenotype","Choose Phenotypes:", choices=c(unique(druggene_pheno$Phecode_description)), selected='All', options = list(create = TRUE, placeholder = 'Phenotypes'))
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
#  tabPanel('Phenotypes',

Phenotype_page <- tabPanel(
  #shinythemes::themeSelector(),  # <--- Add this somewhere in the UI
  
  title='Phenotype', value="Phenotype",
  titlePanel("Priority Index Score"),
  #  theme=shinythemes::shinytheme('flatly'),
  
  
  # tabPanel('Search by Phenotype',value='pheno',
  sidebarLayout(
    sidebarPanel(
      
      # Add OR cutoff slider             
      sliderInput('ORselection_phenotype', 'Use OR cutoff', min=4, max=6, value=5, step=1),
      # Add choices for DG and therpeutic target          
      #  radioButtons('drug_genome',label='Druggable Gene', choices=c("Yes"="DRUGGABLE GENOME", "No"="Not druggable","Select either"="Either"),  selected="DRUGGABLE GENOME"),
      # radioButtons('therapeutic',label='Therapeutic Target', choices=c("Yes"=1, "No"=0, "Select either"="Either"),  selected=1),
      # Select Gene and then phenotypes          
      selectizeInput("choosephenotype","Select Phenotype:", choices=NULL,   options = list(create = TRUE, placeholder = 'Select Phenotype')),
      selectizeInput("choosegene","Choose Genes:", choices=NULL, selected='All', options = list(create = TRUE, placeholder = 'Genes'))
    ),
    
    
    mainPanel(
      tabsetPanel(
      
        tabPanel('Table',
               #box(div(style="dborder-color:black;background-color: white;",
               
                 DT::dataTableOutput("table_pheno")),
        #)),
        tabPanel('Evidence', 
                 plotly::plotlyOutput('plot_evidence_phenotype'))
        
      ))            
    
  ))



#  tabPanel('Phenotypes',
#          plotly::plotlyOutput('plot_pt')),
#)



ui <- fluidPage(
  navbarPage(
    id='tabset',
   title = 
    div(
    div(
    img(src="/PIs-3.png",
        height = "100px",width = "100px",style = "position: relative; margin:-35px -35px; display:right-align;")),''),
  position = 'fixed-top',
  header=tags$style(type="text/css", "body {padding-top: 70px;},
  
  "),
  #footer=tags$footer(title="Your footer here", align = "right", style = "
#position:absolute;
#bottom:0;
#width:100%;
#height:80px; /* Height of the footer */
#color: black;
#padding: 10px;
#background-color: white;
#z-index: 1000;"
 # ),
  
    theme=shinythemes::shinytheme('cosmo'),
  # shinythemes::themeSelector(),  # <--- Add this somewhere in the UI
  
  #theme = bs_theme(bootswatch = "cosmo"),
 # Create Right Side Logo/Image with Link       
 #tags$script(HTML("var header = $('.navbar > .container-fluid');
#header.append('<div style=\"float:right\"><a href=\"'https://labs.icahn.mssm.edu/dolab/people/'\"><img src=\"Logo.png\" alt=\"alt\" style=\"float:right;width:33px;height:41px;padding-top:10px;\"> </a></div>');
 #   console.log(header)")
 #),
 
 # titlePanel("Priority Index Score"),
 #tags$a(href='https://labs.icahn.mssm.edu/dolab/people/',
 #      tags$img(src="Mountsinailogo.jpeg", height='50', width='50'),
 #     'Do Lab')

  about_page,
  Gene_page, 
  Phenotype_page,
  Aboutinfo_page
))
#,

#tags$script(HTML("var header = $('.navbar > .container-fluid');
#header.append('<div style=\"float:right\"><a href=\"https://labs.icahn.mssm.edu/dolab/people/\"><img src=\"questionmark-white.png\" alt=\"alt\" style=\"float:right;width:35px;height:35px;padding-top:10px;\"> </a></div>');
 #   console.log(header)"),
  #          tags$script(HTML("var header = $('.navbar > .container-fluid');
 #header.append('<div style=\"float:right\"><h3>Do Lab</h3></div>');
#"))
      
            
#)
            
#)

#)



#  listchoic=unique(druggene_pheno$Phecode_description[druggene_pheno$Druggable_gene=='DRUGGABLE GENOME' & druggene_pheno$therapeutic==1])

server <- function(input, output,session) {
  
  #output$home_img <- renderImage({
  #list(src = "/Users/aineduffy/Documents/PhD/DrugGenetics/Dataset/Figures paper MI/TIS.png",
  #        list(src = "/Users/aineduffy/Documents/PhD/DrugGenetics/Dataset/TIS_searchpanel/PIs.png",
  #list(src = "/hpc/users/duffya02/MainIndication_website/TIS/TIS.png",
  #list(src = "TIS.png",
  
  #  width = "100%",
  # height = 330)
  
  #}, deleteFile = F)
  
  
  observeEvent(input$Search_by_gene, {
    updateTabsetPanel(session, "tabset" , selected="Gene")
  })
  
  observeEvent(input$Search_by_phenotype, {
    updateTabsetPanel(session, "tabset" , selected="Phenotype")
  })
  
  #observe({
   # or_gene_cutoff=unique(druggene1$Gene[druggene1$max>=input$ORselection])
   # updateTextInput(session, "selectgene", value=or_gene_cutoff, placeholder ='Search Gene' ) 
  #})
  
#  searchInput("selectgene","Search Gene Target:",placeholder = 'Search Gene',  btnSearch = icon("magnifying-glass"),btnReset = icon("xmark"),width = "100%"),
  
  ## panel Gene change drop down lists 
  
  #  listchoic_pheno=unique(druggene_pheno$Phecode_description[druggene_pheno$Gene=='FRRS1L' & druggene_pheno$max>=5 ])
  
  observe({
    listchoic_pheno=unique(druggene_pheno$Phecode_description[druggene_pheno$Gene==input$selectgene & druggene_pheno$max>=input$ORselection ])
    updateSelectizeInput(session, "selectphenotype", choices=listchoic_pheno,selected = NULL) 
    
  })
  
  #  or_gene_cutoff_pheno=unique(druggene_pheno$Phecode_description[druggene_pheno$max>=5])
  
  
  ## panel phenotype  change drop down lists 
  
  observe({
    or_gene_cutoff_pheno=unique(druggene_pheno$Phecode_description[druggene_pheno$max>=input$ORselection_phenotype])
    or_gene_cutoff_pheno1=or_gene_cutoff_pheno[or_gene_cutoff_pheno!='All']
    updateSelectizeInput(session, "choosephenotype", choices=or_gene_cutoff_pheno1) 
  })
  
  
  
  observe({
    listchoic_pheno=druggene_pheno$Gene[druggene_pheno$Phecode_description==input$choosephenotype & druggene_pheno$max>=input$ORselection_phenotype]
    listchoic_pheno1=c('All',listchoic_pheno)
    updateSelectizeInput(session, "choosegene", choices=listchoic_pheno1) 
    
  })
  
  
  datafile_react<- reactive({
    if ('All' %in% input$selectphenotype) {
      Drugfile1=Drugfile %>% filter(OR_cutoff >=input$ORselection)   %>% 
        filter(Gene %in% input$selectgene) %>% select(Gene, Phecode_Integer,Phecode_description,Open_Targets_Indication, Open_Targets_Clinical_Phase,Sider_Indication,Category, Number_of_predictors,Genescore_sum,OR_cutoff,Percent) %>% distinct()
      colnames(Drugfile1)=gsub('_',' ', colnames(Drugfile1))
      return(Drugfile1)
      
    } else {
      
      Drugfile1=Drugfile %>%filter(OR_cutoff >=input$ORselection) %>%
        filter(Gene %in% input$selectgene)  %>%  filter(Phecode_description %in% input$selectphenotype)  %>% select(Gene, Phecode_Integer,Phecode_description,Open_Targets_Indication, Open_Targets_Clinical_Phase,Sider_Indication,Category, Number_of_predictors,Genescore_sum,OR_cutoff, Percent)  %>% distinct()
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
        filter(Phecode_description %in% input$choosephenotype)%>% select( Phecode_Integer,Phecode_description,Gene,Open_Targets_Indication, Open_Targets_Clinical_Phase,Sider_Indication,Category, Number_of_predictors,Genescore_sum,OR_cutoff,Percent) %>% distinct()
      colnames(Drugfile1)=gsub('_',' ', colnames(Drugfile1))
      return(Drugfile1)
    } else {
      
      Drugfile1=Drugfile %>%filter(OR_cutoff >=input$ORselection_phenotype) %>%
        filter(Gene %in% input$choosegene)  %>%  filter(Phecode_description %in% input$choosephenotype) %>% select( Phecode_Integer,Phecode_description,Gene,Open_Targets_Indication, Open_Targets_Clinical_Phase,Sider_Indication,Category, Number_of_predictors,Genescore_sum,OR_cutoff,Percent) %>% distinct()
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
        filter(Gene %in% input$selectgene) %>% 
        select(Phecode_description,OMIM,Coloc_phenotype,Eqtl_phenotype,Coloc_tissue,
               HGMD,Genebass,
               Clinvar,Open_Targets_Genetics_Portal,Tissue_specific,Eqtl_tissue, OE)
      
      Drugfile2<-reshape2::melt(Drugfile1, id.vars=c('Phecode_description'))
      
      labels=Drugfile %>%filter(OR_cutoff >=input$ORselection) %>%
        filter(Gene %in% input$selectgene) %>% 
        select(Phecode_description,OMIM_description,Coloc_phenotype_description,Eqtl_phenotype_description,Coloc_tissue,
               HGMD_description,Genebass_description,Clinvar_description,
               Open_Targets_Genetics_Portal_description,Tissue_specific, Eqtl_tissue_description,OE)
      
      colnames(labels)=colnames(Drugfile1)
      labels2<-reshape2::melt(labels, id.vars=c('Phecode_description'),value.name = 'label')
      Drugfile2_labels=inner_join(Drugfile2,labels2, by = c("Phecode_description", "variable") )
      #      Drugfile2_labels$Phecode_description<- factor(Drugfile2_labels$Phecode_description, levels=unique(Drugfile2_labels$Phecode_description))
      Drugfile2_labels$value<-as.character(Drugfile2_labels$value)
      Drugfile2_labels=Drugfile2_labels%>% distinct() 
      Drugfile2_labels$variable=gsub('_',' ', Drugfile2_labels$variable)
      return(Drugfile2_labels)
    } else { 
      
      Drugfile1=Drugfile %>%filter(OR_cutoff >=input$ORselection) %>%
        filter(Gene %in% input$selectgene)  %>%  filter(Phecode_description %in% input$selectphenotype) %>% 
        select(Phecode_description,OMIM,Coloc_phenotype,Eqtl_phenotype,Coloc_tissue,
               HGMD,Genebass,
               Clinvar,Open_Targets_Genetics_Portal,Tissue_specific,Eqtl_tissue,OE)
      
      Drugfile2<-reshape2::melt(Drugfile1, id.vars=c('Phecode_description'))
      
      labels=Drugfile %>%filter(OR_cutoff >=input$ORselection) %>%
        filter(Gene %in% input$selectgene) %>% 
        select(Phecode_description,OMIM_description,Coloc_phenotype_description,Eqtl_phenotype_description,Coloc_tissue,
               HGMD_description,Genebass_description,Clinvar_description,
               Open_Targets_Genetics_Portal_description,Tissue_specific, Eqtl_tissue_description,OE)

      colnames(labels)=colnames(Drugfile1)
      labels2<-reshape2::melt(labels, id.vars=c('Phecode_description'),value.name = 'label')
      Drugfile2_labels=inner_join(Drugfile2,labels2, by = c("Phecode_description", "variable") )
      #      Drugfile2_labels$Phecode_description<- factor(Drugfile2_labels$Phecode_description, levels=unique(Drugfile2_labels$Phecode_description))
      Drugfile2_labels$value<-as.character(Drugfile2_labels$value)
       Drugfile2_labels=Drugfile2_labels%>% distinct() 
       Drugfile2_labels$variable=gsub('_',' ', Drugfile2_labels$variable)
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
    
    ggplotly(evid, tooltip="text") 
  })
  
  # Evidence phenotype
  
  ##Evidence tables
  
  datafile_react_evidence_phenotype<- reactive({
    
    
    if ('All' %in% input$choosegene) {
      
      Drugfile1=Drugfile %>%filter(OR_cutoff >=input$ORselection_phenotype) %>%
        filter(Phecode_description %in% input$choosephenotype)%>% 
        select(Gene,OMIM,Coloc_phenotype,Eqtl_phenotype,Coloc_tissue,
               HGMD,Genebass,
               Clinvar,Open_Targets_Genetics_Portal,Tissue_specific,Eqtl_tissue,OE)
      
      Drugfile2<-reshape2::melt(Drugfile1, id.vars=c('Gene'))
      
      labels=Drugfile %>%filter(OR_cutoff >=input$ORselection_phenotype) %>%
        filter(Phecode_description %in% input$choosephenotype) %>% 
        select(Gene,OMIM_description,Coloc_phenotype_description,Eqtl_phenotype_description,Coloc_tissue,
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
      return(Drugfile2_labels)
    } else {
      
      Drugfile1=Drugfile %>%filter(OR_cutoff >=input$ORselection_phenotype) %>%  filter(Gene %in% input$choosegene)   %>% 
        filter(Phecode_description %in% input$choosephenotype)%>% 
        select(Gene,OMIM,Coloc_phenotype,Eqtl_phenotype,Coloc_tissue,
               HGMD,Genebass,
               Clinvar,Open_Targets_Genetics_Portal,Tissue_specific,Eqtl_tissue,OE)
      
      Drugfile2<-reshape2::melt(Drugfile1, id.vars=c('Gene'))
      
      labels=Drugfile %>%filter(OR_cutoff >=input$ORselection_phenotype) %>%
        filter(Phecode_description %in% input$choosephenotype) %>% 
        select(Gene,OMIM_description,Coloc_phenotype_description,Eqtl_phenotype_description,Coloc_tissue,
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
      return(Drugfile2_labels)
      
    }
  })  
  
  
  output$plot_evidence_phenotype <- plotly::renderPlotly({
    evid=datafile_react_evidence_phenotype() %>%   
      ggplot(aes(variable, Gene, fill= value, text=label)) + 
      geom_tile(colour = "white", lwd = 1.5,
                linetype = 1) + theme(axis.text.x = element_text(angle = 45), axis.title.x=element_blank(),
                                      axis.title.y=element_blank()) +
      # scale_y_discrete(limits = rev(levels(Drugfile2_labels$Gene))) +
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
  
  
  
  # if ('All' %in% input$choosegene)  %>%  filter(Phecode_description %in% input$choosephenotype)
  
  
  
  
  
  
  
  
  
}

shinyApp(ui, server)


