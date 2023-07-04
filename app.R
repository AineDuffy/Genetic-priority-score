
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
Drugfile<-fread('data/Allgenes_prioritytable_miphenotypes_withdescription_updated.txt.gz', data.table =F, quote="")
druggene_pheno<-Drugfile %>% group_by(Gene, `Phecode description`) %>% mutate( max=max(gene_score_cutoff)) %>% distinct(Gene, `Phecode description`, max) 

AllPhenotypes <- sort(unique(Drugfile$`Phecode description`))
for (pheno in AllPhenotypes){
  PhenoData <- Drugfile %>% filter(`Phecode description` %in% pheno)
  assign(gsub(" ", "_", pheno), PhenoData)
}


# Define UI ----

about_page <- tabPanel(
    title = "Home",
    setBackgroundColor(
        color = c("#F5FFFA", "#E6E6FA"),
        gradient = "linear",#E6E6FA"
        direction = "bottom"
    ),
      tagList(
        lapply(1:6, function(i) br())),
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
    tagList(
        lapply(1:3, function(i) br())),
    h4(paste('A human genetics-guided priority tool for informing therapeutic target selection for',comma(length(unique((Drugfile$Gene)))), 'genes and', length(unique((Drugfile$`Phecode Integer`))), 'phenotypes.'), align='center'),
    h4(paste('Last update:',Sys.Date()), align='center'),
   tagList(
        lapply(1:4, function(i) br())),
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
     tagList(
        lapply(1:3, function(i) br())))


## About page 
source('Summary_about_page.R', local = TRUE)

theme_flatly <- shinytheme("flatly")
               
Gene_page <- tabPanel(

    title='Gene',value="Gene",
    titlePanel("Genetic Priority Score"),
    theme=theme_flatly,
    sidebarLayout(
        sidebarPanel(

            # Add OR cutoff slider             
            sliderInput('ORselection', 'Use score cutoff', min=0, max=2.1, value=1.5, step=0.3),
             # Select Gene and then phenotypes          
            searchInput("selectgene","Search Gene Target:",placeholder = 'Search Gene',  btnSearch = icon("search"),btnReset = icon("remove"),width = "100%"),
            helpText("Gene examples: ABCC8, JAK2, PCSK9"),
        ),
        mainPanel(
            tabsetPanel(
                tabPanel('Table',    
                         dataTableOutput("table")),
                tabPanel('Evidence', 
                         plotlyOutput('plot_evidence')),
                tabPanel('DOE Evidence', 
                         plotlyOutput('plot_evidence_doe'))                
            )                                    
        ))
)

Phenotype_page <- tabPanel(

    title='Phenotype', value="Phenotype",
    titlePanel("Genetic Priority Score"),
    sidebarLayout(
        sidebarPanel(
            # Add OR cutoff slider             
            sliderInput('scoreselection_phenotype', 'Use score cutoff',  min=0, max=2.1, value=1.5, step=0.3),
            selectizeInput("choosephenotype","Select Phenotype:", choices=AllPhenotypes,  multiple = FALSE,options = list(create = TRUE, placeholder = 'Select Phenotype'))
                                ),      
        mainPanel(
            tabsetPanel(
                tabPanel('Table',
                        dataTableOutput("table_pheno")),
                tabPanel('Evidence', 
                        plotlyOutput('plot_evidence_phenotype')),
                tabPanel('DOE Evidence', 
                         plotlyOutput('plot_evidence_phenotype_doe'))
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
        
        theme=shinytheme('cosmo'),
        about_page,
        Gene_page, 
        Phenotype_page,
        Aboutinfo_page
    ))

server <- function(input, output,session) {
      phenotypes <- reactiveValues()
  
  observeEvent(input$Search_by_gene, {
    updateTabsetPanel(session, "tabset", selected = "Gene")
  })
  
  observeEvent(input$Search_by_phenotype, {
    updateTabsetPanel(session, "tabset", selected = "Phenotype")
  })
  
  observe({
    or_gene_cutoff_pheno <- druggene_pheno %>%
      filter(max >= input$scoreselection_phenotype) %>%
      ungroup() %>%
      distinct(`Phecode description`) %>%
      arrange(`Phecode description`)
    or_gene_cutoff_pheno <- or_gene_cutoff_pheno$`Phecode description`
    # Store the current value of the selectizeInput
    current_value <- input$choosephenotype
    
    if (!is.null(current_value) && !(current_value %in% or_gene_cutoff_pheno)) {
      current_value <- NULL
    }
    # Update the choices and selected value of selectizeInput
    updateSelectizeInput(session, "choosephenotype", choices = or_gene_cutoff_pheno, selected = current_value)
    
    # Store the updated choices in reactiveValues for later use
    phenotypes$choices <- or_gene_cutoff_pheno
  })
  
    ## Gene Page
    
    datafile_react<- reactive({
        Drugfile1=Drugfile %>% filter(gene_score_cutoff >= input$scoreselection,
                                      Gene %in% input$selectgene) %>% select(Gene, `Phecode Integer`,`Phecode description`,`Open Targets Indication`, `Open Targets Clinical Phase`,`SIDER Indication`,Category, `Number of predictors`, `Number of DOE predictors`,`Genetic priority score`,`Genetic priority score - DOE`)   %>%
            arrange(desc(`Genetic priority score`))
        return(Drugfile1)
        
    })
    
    output$table <- renderDataTable({
        datafile_react()
    })
    
    output$table_pheno <- renderDataTable({
      #  datafile_react_pheno()
      Drugfile1 <- get(gsub(" ", "_", input$choosephenotype))  %>%
        filter(gene_score_cutoff >= input$scoreselection_phenotype) %>%
        distinct(Gene, `Phecode Integer`,`Phecode description`,`Open Targets Indication`, `Open Targets Clinical Phase`,`SIDER Indication`,Category, `Number of predictors`, `Number of DOE predictors`, `Genetic priority score`,`Genetic priority score - DOE`) %>% arrange(desc(`Genetic priority score`))
      return(Drugfile1)
    })
    
    ##Evidence tables
    
    datafile_react_evidence<- reactive({
        Drugfile1=Drugfile %>%filter(gene_score_cutoff >=input$scoreselection, Gene %in% input$selectgene) %>%
            arrange(desc(`Genetic priority score`)) %>%
            select(`Phecode description`, contains("phenotype predictor")) %>%
            rename_with(~ gsub(' phenotype predictor', '', .))
        
        Drugfile1<-reshape2::melt(Drugfile1, id.vars=c('Phecode description'))
        
        labels=Drugfile %>%filter(gene_score_cutoff >=input$scoreselection, Gene %in% input$selectgene) %>%
            arrange(desc(`Genetic priority score`))%>%
            select(contains("description")) %>%
            rename_with(~ gsub(' description', '', .)) %>% rename(`Phecode description` = Phecode)
        
        labels<-reshape2::melt(labels, id.vars=c('Phecode description'),value.name = 'label')
        Drugfile2_labels=inner_join(Drugfile1,labels, by = c("Phecode description", "variable") )
        Drugfile2_labels$value<-as.character(Drugfile2_labels$value)
        Drugfile2_labels=Drugfile2_labels%>% distinct()
        Drugfile2_labels$label=gsub(";", "; ",     Drugfile2_labels$label)
        Drugfile2_labels$label=gsub("([a-z0-9]* [a-z0-9]* [a-z0-9]*) ", "\\1\n",     Drugfile2_labels$label)
        return(Drugfile2_labels)
    })
    
    
    output$plot_evidence <- plotly::renderPlotly({
        evid=datafile_react_evidence() %>%
        ggplot(aes(x=variable, y=`Phecode description`,  fill= value, text=label)) +
            geom_tile(colour = "white", lwd = 1.5, linetype = 1) + 
            theme(axis.text.x = element_text(angle = 45), axis.title.x=element_blank(), axis.title.y=element_blank()) +
            scale_fill_manual(values=c("azure3","blueviolet")) +labs(fill='Evidence') +  scale_y_discrete(limits=rev)
        
        if(length(unique(evid$data$`Phecode description`))<7){
            heightplot=700
        } else {
            heightplot=length(unique(evid$data$`Phecode description`))*40
        }
        
        p <- ggplotly(evid, tooltip="text", height=heightplot) 
        p <- layout(p , xaxis = list(side = "top"))
        p
        
    })
    
    ### DOE Evidence     
    ##Evidence tables
    datafile_react_evidence_doe<- reactive({
        Drugfile1=Drugfile %>%filter(gene_score_cutoff >=input$scoreselection) %>%
            filter(Gene %in% input$selectgene )%>% arrange(desc(`Genetic priority score`)) %>%
            select(`Phecode description`, contains("phenotype predictor")) %>%
            rename_with(~ gsub(' phenotype predictor', '', .))
        
        Drugfile1<-reshape2::melt(Drugfile1, id.vars=c('Phecode description'))
        
        labels=Drugfile %>%filter(gene_score_cutoff >= input$scoreselection) %>%
            filter(Gene %in% input$selectgene ) %>% arrange(desc(`Genetic priority score`))%>%
            select(`Phecode description`, contains("DOE predicted")) %>%
            rename_with(~ gsub(' DOE predicted', '', .))
        
        labels<-reshape2::melt(labels, id.vars=c('Phecode description'),value.name = 'label')
        Drugfile2_labels=inner_join(Drugfile1,labels, by = c("Phecode description", "variable") )
        Drugfile2_labels$value<-as.character(Drugfile2_labels$value)
        #   Drugfile2_labels=Drugfile2_labels%>% distinct()
        Drugfile2_labels$label=gsub("([a-z0-9]* [a-z0-9]* [a-z0-9]*) ", "\\1\n",     Drugfile2_labels$label)
        return(Drugfile2_labels)
    })
    
    
    output$plot_evidence_doe <- plotly::renderPlotly({
        evid=datafile_react_evidence_doe() %>%
        ggplot( aes(x=variable, y=`Phecode description`, fill= label)) +
            geom_tile(colour = "white", lwd = 1.5, linetype = 1) + theme(axis.text.x = element_text(angle = 45), axis.title.x=element_blank(), axis.title.y=element_blank()) +
            scale_fill_manual(values=c("0"="azure3","GOF"= "dodgerblue1", "LOF"="violetred2")) +labs(fill='Evidence') +  scale_y_discrete(limits=rev)
        
        
        
        if(length(unique(evid$data$`Phecode description`))<7){
            heightplot=700
        } else {
            heightplot=length(unique(evid$data$`Phecode description`))*40
        }
        
        p<- ggplotly(evid, tooltip="text",height=heightplot)
        p <- layout(p , xaxis = list(side = "top"))
        p
    })
    
    ## Phenotype Page
    
    # Evidence phenotype
    
    ##Evidence tables
    
    datafile_react_evidence_phenotype<- reactive({
        
        Drugfile1=Drugfile %>%filter(gene_score_cutoff >=input$scoreselection_phenotype) %>%
            filter(`Phecode description` %in%  input$choosephenotype)%>%
            arrange(desc(`Genetic priority score`)) %>%
            select(Gene, `Genetic priority score`, contains("phenotype predictor")) %>%
            rename_with(~ gsub(' phenotype predictor', '', .))
        
        Drugfile1<-reshape2::melt(Drugfile1, id.vars=c('Gene', 'Genetic priority score'))
        
        labels=Drugfile %>%filter(gene_score_cutoff >=input$scoreselection_phenotype, `Phecode description` %in%  input$choosephenotype) %>%
            select( Gene, `Genetic priority score`, contains("description")) %>%
            rename_with(~ gsub(' description', '', .)) %>% rename(`Phecode description` = Phecode)
        
        labels<-reshape2::melt(labels, id.vars=c('Gene','Genetic priority score'),value.name = 'label')
        Drugfile2_labels=inner_join(Drugfile1,labels, by = c("Gene", "Genetic priority score", "variable") )
        Drugfile2_labels=Drugfile2_labels %>% distinct()
        Drugfile2_labels$label=gsub(";", "; ",     Drugfile2_labels$label)
        
        Drugfile2_labels$label=gsub("([a-z0-9]* [a-z0-9]* [a-z0-9]*) ", "\\1\n",     Drugfile2_labels$label)
        if( length(unique(Drugfile2_labels$Gene)) >100) {
            
            topgenepheno=head(Drugfile2_labels  %>% distinct(Gene,`Genetic priority score`),100)
            Drugfile2_labels_topgenes=inner_join(Drugfile2_labels,topgenepheno)
            Drugfile2_labels_topgenes=Drugfile2_labels_topgenes %>% select(-`Genetic priority score`)
            
        }else {
            
            Drugfile2_labels_topgenes=Drugfile2_labels %>% select(-`Genetic priority score`)
        }
        Drugfile2_labels_topgenes$Gene<- factor(Drugfile2_labels_topgenes$Gene, levels=unique(Drugfile2_labels_topgenes$Gene))
        Drugfile2_labels_topgenes$value<-as.character(Drugfile2_labels_topgenes$value)
        return(Drugfile2_labels_topgenes)
    })
    
    
    output$plot_evidence_phenotype <- plotly::renderPlotly({
        evid=datafile_react_evidence_phenotype() %>%
        ggplot(aes(variable, Gene, fill= value, text=label)) +
            geom_tile(colour = "white", lwd = 1.5, linetype = 1) +
            theme(axis.text.x = element_text(angle = 45), axis.title.x=element_blank(),axis.title.y=element_blank()) +
            scale_fill_manual(values=c("azure3","blueviolet")) +
            labs(fill='Evidence') +  scale_y_discrete(limits=rev)        

        if(length(unique(evid$data$Gene))<7){
            heightplot=700
        } else {
            heightplot=length(unique(evid$data$Gene))*40
        }  
        p<- ggplotly(evid, tooltip="text",height=heightplot)
        p <- layout(p , xaxis = list(side = "top"))
        p
        
    })
    
    datafile_react_evidence_phenotype_doe<- reactive({
        
        Drugfile1=Drugfile %>%filter(gene_score_cutoff >=input$scoreselection_phenotype, `Phecode description` %in%  input$choosephenotype) %>%
            arrange(desc(`Genetic priority score`)) %>%
            select(Gene, `Genetic priority score`, contains("phenotype predictor")) %>%
            rename_with(~ gsub(' phenotype predictor', '', .))
        
        Drugfile1<-reshape2::melt(Drugfile1, id.vars=c('Gene', 'Genetic priority score'))
        
        labels=Drugfile %>%filter(gene_score_cutoff >=input$scoreselection_phenotype, 
                                  `Phecode description` %in%  input$choosephenotype) %>%
            select( Gene, `Genetic priority score`, contains("DOE predicted")) %>%
            rename_with(~ gsub(' DOE predicted', '', .))
        
        labels<-reshape2::melt(labels, id.vars=c('Gene','Genetic priority score'),value.name = 'label')
        Drugfile2_labels=inner_join(Drugfile1,labels, by = c("Gene", "Genetic priority score", "variable") )
        
        if( length(unique(Drugfile2_labels$Gene)) >100) {
            
            topgenepheno=head(Drugfile2_labels  %>% distinct(Gene,`Genetic priority score`),100)
            Drugfile2_labels_topgenes=inner_join(Drugfile2_labels,topgenepheno)
            Drugfile2_labels_topgenes=Drugfile2_labels_topgenes %>% select(-`Genetic priority score`)
            
        }else {
            Drugfile2_labels_topgenes=Drugfile2_labels %>% select(-`Genetic priority score`)
        }
        
        Drugfile2_labels_topgenes$Gene<- factor(Drugfile2_labels_topgenes$Gene, levels=unique(Drugfile2_labels_topgenes$Gene))
        Drugfile2_labels_topgenes$value<-as.character(Drugfile2_labels_topgenes$value)
        
        return(Drugfile2_labels_topgenes)
    })
    
    output$plot_evidence_phenotype_doe <- plotly::renderPlotly({
        evid=datafile_react_evidence_phenotype_doe() %>%
         ggplot(aes(variable, Gene, fill= label, text=label)) +
            geom_tile(colour = "white", lwd = 1.5,
                      linetype = 1) + theme(axis.text.x = element_text(angle = 45), axis.title.x=element_blank(),  axis.title.y=element_blank()) +
            scale_fill_manual(values=c("0"="azure3","GOF"= "dodgerblue1", "LOF"="violetred2")) +labs(fill='Evidence') +
          scale_y_discrete(limits=rev)        
        if(length(unique(evid$data$Gene))<7){
            heightplot=700
        } else {
            heightplot=length(unique(evid$data$Gene))*40
        }
        
        p<-  ggplotly(evid, tooltip="text",height=heightplot)
        p <- layout(p , xaxis = list(side = "top"))
        p
        
    })
    
}

shinyApp(ui, server)
