## Author: Neha Rao
## Course: BF591
## Final Project
## RShiny app



library(shiny)
library(bslib)
library(tidyverse)
library(ggplot2)
library('ggbeeswarm')
library(colourpicker)
library(DESeq2)
library(genefilter)
library(dplyr)
library(gplots)
library('RColorBrewer')
library("ggfortify")
library("shinythemes")
library(readr)

# To increase the file size capabilities of the file inputs
options(shiny.maxRequestSize=30*1024^2)

# Define UI for application
ui <- fluidPage(
  tags$style(HTML("
    body {
      background-color: #cfd8ff;
      font-family: 'Times New Roman';
    }
    .tab-content {
      background-color: #e1e6fa;
      padding: 20px;
      border-radius: 5px;
    }
    .tab-pane {
      padding: 10px;
    }
    .sidebar {
      background-color: #BBDEFB;
      padding: 20px;
      border-radius: 5px;
    }
    .navbar {
      background-color: #64B5F6;
      border: none;
      border-radius: 0;
    }
    .navbar-brand {
      color: #FFFFFF;
      font-weight: bold;
    }
    .nav-tabs > li > a {
      background-color: #BBDEFB;
      color: #0D47A1;
      font-weight: bold;
      border-radius: 5px;
    }
    .nav-tabs > li.active > a {
      background-color: #2196F3;
      color: #FFFFFF;
      border-color: #2196F3;
    }
    .btn {
      background-color: #64B5F6;
      border-color: #64B5F6;
      color: #FFFFFF;
    }
    .btn:hover {
      background-color: #42A5F5;
      border-color: #42A5F5;
    }
  ")),
  titlePanel(h5("BF591: Final Project", align = "center")),
  titlePanel(h4("Neha Rao", align = "center")),
  tabsetPanel(
    type = "tabs",
    tabPanel("Sample Exploration",
             p("Uploading the input metadata csv file and exploring the sample information"),
             sidebarLayout(
               sidebarPanel(
                 fileInput("Sample_file", "Choose a Sample Information File (only CSV File)", accept = ".csv", placeholder = "Sample information file"),
                 submitButton("Submit"), width = 3
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("Summary", p("This sub-tab provides the summary of the statistics") , tableOutput("summary_table")),
                   tabPanel("Table", p("This sub-tab provides a table of the sample information"), DT::dataTableOutput("sample_table")),
                   tabPanel("Plots", p("This sub-tab provides a histogram plot for the choosen variable"),
                            sidebarLayout(
                              sidebarPanel(
                                radioButtons("hist_x_axis", "X-Axis",
                                             choices = c("Sample_pmi" = "Sample_pmi","Sample_Age_of_Death" = "Sample_Age_of_Death", "Sample_rin"="Sample_rin","Sample_mrna_seq"="Sample_mrna_seq")),
                                colourInput("samp_color_1", "Outline color","#9e1840"),
                                colourInput("samp_color_2", "Fill color color","#187f9e"),
                                submitButton("Submit")
                              ),
                              mainPanel(plotOutput("sample_plot"))
                            )
                   )
                 )
               )
             )
    ),
    tabPanel("Counts Exploration",
             p("Uploading the input metadata csv file and exploring the counts data"),
             sidebarLayout(
               sidebarPanel(
                 fileInput("Normalized_counts_matrix", "Choose a Normalized Counts Matrix (only CSV File)", accept = ".csv", placeholder = 'Normalized counts matrix file'),
                 sliderInput("var_slider", "Percentile of Variance", min = 0, max = 100, value = 40),
                 sliderInput("zero_slider", "Number of Samples that are of non-zero variance", min = 0, max = 100, value = 40),
                 submitButton("Submit"), width = 3
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("Table",
                            p("This sub-tab provides the statistics based on the input parameters"),
                            tableOutput("counts_table")),
                   tabPanel("Scatter Plot",
                            hr("Scatter plot tailored to the variance input"),
                            plotOutput("variance_plot"),
                            hr("Scatter plot tailored to the non-zero input"),
                            plotOutput("zero_plot")),
                   tabPanel("Heat map",
                            p("This sub-tab provides a heat map of the normalized counts input matrix"), 
                            plotOutput("heatmap")),
                   tabPanel("PCA", selectInput(inputId = "comp1", label="Select X-axis", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")),
                            selectInput(inputId = "comp2", label="Select Y-axis", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), selected = "PC2"),
                            p("This sub-tab provides a PCA plot"),
                            plotOutput("counts_pca"))
                   
                 )
               )
             )),
    tabPanel("Differential Expression",
             p("Uploading the input differential analysis csv file and exploring the results"),
             sidebarLayout(
               sidebarPanel(
                 fileInput("DE_Anaylsis", "Choose a Differential Expression Results File (only CSV File)", accept = ".csv",placeholder = 'Differential expression analysis file'),
                 submitButton(text = 'Submit'), width = 3
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("Table",
                            p("This sub-tab provides a table containing the differential expression analysis results."),
                            DT::dataTableOutput("table",width = "80%")),
                   tabPanel("Volcano Plot", 
                            p("This sub-tab provides a volcano plot"),
                            sidebarLayout(
                              sidebarPanel(
                                radioButtons("x_axis", "X-Axis",
                                             choices = c("log2FoldChange" = "log2FoldChange","baseMean" = "baseMean", "lfcSE"="lfcSE","stat"="stat","pvalue"="pvalue","padj"="padj"),
                                             selected = "log2FoldChange"),
                                radioButtons("y_axis", "Y-Axis",
                                             choices = c("log2FoldChange" = "log2FoldChange","baseMean" = "baseMean", "lfcSE"="lfcSE","stat"="stat","pvalue"="pvalue","padj"="padj"),
                                             selected = "padj"),
                                sliderInput(inputId = "padjusted", min = -35, max = 0,label = "Select the magnitude of p-adj", value = -15, step = 1),
                                colourInput("de_base", "Base point color","#9e1864"),
                                colourInput("de_highlight", "Highlight point color","#60f0ce"),
                                submitButton("Submit"),width=3
                              ),
                              mainPanel(
                                plotOutput("volcano"))))
                 )
               )
             )
    ),
    tabPanel("GSEA", 
             sidebarLayout(
               sidebarPanel(width=3,
                            fileInput("fgsea_file", "FGSEA File", accept = ".csv", placeholder = "fgsea_file.csv")
                            #submitButton("Upload FGSEA File")
               ), #sidebarPanel
               mainPanel(width = 9,
                         tabsetPanel(type = "tabs",
                                     tabPanel("Top Results", 
                                              sidebarPanel(width=4,
                                                           sliderInput(inputId="GSEA_slider1", label="Top Results by adjusted p-value", min=0, max=100, value=10, step = NULL,round = FALSE,ticks = TRUE),
                                                           submitButton("Update Top Pathways")
                                              ), #sidebarPanel
                                              mainPanel(width = 5,
                                                        tabsetPanel(#type = "tabs",
                                                          type = "hidden",
                                                          tabPanel("GSEA_Plot1", plotOutput("GSEA_Plot1"))
                                                        ) #tabsetPanel
                                              ) #mainPanel
                                     ),
                                     tabPanel("Table", 
                                              sidebarPanel(width=4,
                                                           sliderInput(inputId="GSEA_slider2", label="Filter table by adjusted p-value", min=0, max=1, value=0.1, step = NULL,round = FALSE,ticks = TRUE),
                                                           radioButtons("GSEA_radio2", "Choose all, positive and negative pathways", c("All Pathways"="All","Positive Pathways"="Positive","Negative Pathways"="Negative"), selected="All"),
                                                           submitButton("Update Table View"),
                                                           br(),
                                                           #fileInput("GSEA_file_download", "FGSEA File", accept = ".csv", placeholder = "GSEA_file_download.csv")
                                                           downloadButton("fgsea_NES_Table", "Download Table")
                                              ), #sidebarPanel
                                              mainPanel(width = 5,
                                                        tabsetPanel(type = "hidden",
                                                                    tabPanel("GSEA_Table2", DT::dataTableOutput("GSEA_Table2"))
                                                        ) #tabsetPanel
                                              ) #mainPanel
                                     ),
                                     tabPanel("Plots",
                                              sidebarPanel(width=4,
                                                           sliderInput(inputId="GSEA_slider3", label="Filter table by adjusted p-value", min=0, max=1, value=0.1, step = NULL,round = FALSE,ticks = TRUE),
                                                           submitButton("Update Padj Threshold")
                                              ), 
                                              mainPanel(width = 5,
                                                        tabsetPanel(type = "hidden",
                                                                    tabPanel("GSEA_Plot3", plotOutput("GSEA_Plot3"))
                                                        ) 
                                              ) 
                                     ) 
                         ) 
               ) 
             )            
    )
  ) 
  
) 

#########################################################################

#  Define server logic 
server <- function(input, output,session) {
  
  #######################Tab 01: Sample Data Input#######################
  
  # Tab - 1: Sample Information Exploration
  
  # reading the sample data
  #' @details The purpose of this function is to load in the sample information
  #' .csv file. The user has the ability to input their file 
  sample_data <- reactive({
    req(input$Sample_file)
    data <- read.csv(input$Sample_file$datapath, row.names = 1)
    return(data)
  })
  
  # Summary table 
  #' Here we want to create a dataframe summarizing the sample information file
  #' that the user will input.
  #' @param data This is the data frame that is loaded
  #' @details Column_name This is the column where each variable is mentioned
  #' @details summarized_df This is the dataframe that includes the column 
  #' name, classification,mean and the standard deviation of the variable 
  
  
  summary <- function(data){
    Column_Name <- c("Sample_geo_accession", "Sample_status", "Sample_submission_date", 
                     "Sample_last_update_date", "Sample_type", 
                     "Sample_Channel_Count", "Sample_Organ", "Sample_Organism", 
                     "Sample_Tissue_Type", "Sample_Neurological Diagnosis", 
                     "Sample_pmi", "Sample_Age_of_Death", 'Sample_rin', 
                     "Sample_mrna_seq", "Sample_age_of_onset")
    
    Classification <- c("character", "character", "character", "character", 
                        "character", "double", "character", "character", 
                        "character", "character", "double", "double",
                        "double", "double", "double")
    
    Sample_Channel_Count_average <- mean(data$Sample_Channel_Count)
    Sample_pmi_average <- mean(data$Sample_pmi, na.rm = TRUE)
    Sample_Age_of_Death_average <- mean(data$Sample_Age_of_Death)
    Sample_rin_average <- mean(data$Sample_rin)
    Sample_mrna_seq_average <- mean(data$Sample_mrna_seq)
    Sample_Channel_Count_sd <- sd(data$Sample_Channel_Count)
    Sample_pmi_sd <- sd(data$Sample_pmi, na.rm = TRUE)
    Sample_Age_of_Death_sd <- sd(data$Sample_Age_of_Death)
    Sample_rin_sd <- sd(data$Sample_rin)
    Sample_mrna_seq_sd <- sd(data$Sample_mrna_seq)
    Sample_Neurological_Diagnosis_uniq <- toString(unique(data$Sample_diagnosis))
    
    
    Mean <- c('N/A', 'N/A', 'N/A', 
              'N/A', 'N/A', 
              Sample_Channel_Count_average, 'N/A', 'N/A', 
              'N/A', Sample_Neurological_Diagnosis_uniq, 
              Sample_pmi_average, Sample_Age_of_Death_average, Sample_rin_average, 
              Sample_mrna_seq_average, '41.85')
    
    Standard_deviation <- c('N/A', 'N/A', 'N/A', 
                            'N/A', 'N/A', 
                            Sample_Channel_Count_sd, 'N/A', 'N/A', 
                            'N/A', Sample_Neurological_Diagnosis_uniq, 
                            Sample_pmi_sd, Sample_Age_of_Death_sd, Sample_rin_sd, 
                            Sample_mrna_seq_sd, '10.80')
    
    summarized_df <- data_frame(Column_Name, Classification, Mean, Standard_deviation)
    return(summarized_df)
  }
  
  # plot 
  #' @details  Here the goal is to develop a way for the user to generate a 
  #' histogram plot based on the sample information. 
  #' @param data The df with all of the data that the user has inputted
  #' @param x_name This is x-axis. The user has the ability to choose what 
  #' they would like the x-axis of the plot to be. 
  #' @param color_1 This is the outer color that the user chooses. 
  #' @param color_2 This is the inner color that the user chooses.
  
  plot <- function(data, x_name,color_1,color_2){
    data_histogram <- 
      ggplot (data, aes(x = !!sym(x_name))) +
      geom_histogram(color = color_1, fill = color_2) +
      scale_color_manual(values = c(color_1, color_2))
    theme_bw() +
      theme(legend.position = 'bottom') 
    
    return(data_histogram)
  }
  
  # Summary Table output
  output$summary_table <- renderTable({
    data <- sample_data()
    return(summary(data))
  })
  
  # Sample Data Table output
  output$sample_table <- DT::renderDataTable({ 
    data <- sample_data()
    if (is.null(data)) {
      return()
    }
    isolate({
      return(data)
    })
  })
  
  # plot output
  observeEvent(list(input$hist_x_axis,input$samp_color_1,input$samp_color_2),{
    output$sample_plot <- renderPlot({
      data <- sample_data()
      if (is.null(data)) {
        return()
      }
      isolate({
        return(plot(data,input$hist_x_axis,input$samp_color_1,input$samp_color_2))
      })
    }) 
  })
  
  
  #----------------------------------------------------------------------------------------------------------------------------
  
  #######################Tab 02: Counts Matrix expression#######################
  
  # reading the counts matrix data
  #' Load in the Data 
  #' @details Here all we want to do is simply load in the counts .csv file 
  #' that the user will input via the "Submit" button. 
  counts_data <- reactive({
    req(input$Normalized_counts_matrix)
    data <- read.csv(input$Normalized_counts_matrix$datapath)
    colnames(data)[1] <- "gene"
    return(data)
  })
  
  #' @details The purpose of this function is to filter the data based 
  #' on what value the of number of samples that are greater the variance
  #' threshold and the zero slider threshhold that the user inputs via the slider 
  #' @param data The counts matrix that the user has imported
  #' @param var_slider The slider that allows the user to choose 
  #' variance percentage 
  #' @param zero_slider The slider that the user uses to include genes with at 
  #' least X samples that are non-zero
  
  
  filtering <- function(data,var_slider,zero_slider){
    counts_matrix <- data[,-1]
    variances <- apply(counts_matrix, 1, var, na.rm = TRUE)
    var_threshold <- quantile(variances, var_slider / 100)
    non_zero_counts <- rowSums(counts_matrix > 0, na.rm = TRUE)
    
    filtered <- counts_matrix[variances >= var_threshold & non_zero_counts >= zero_slider, ]
    
    total_samples <- ncol(counts_matrix)
    total_genes <- nrow(counts_matrix)
    passing_genes <- nrow(filtered)
    not_passing_genes <- total_genes - passing_genes
    table <- data.frame( 
      Information = c("Total Number of Samples","Total Number of Genes","Number of Passing Filter", "Number of Not Passing Filter", "Percent of Passing Filter", "Percent of Not Passing Filter"),
      Statistics = c(as.integer(total_samples),as.integer(total_genes),as.integer(passing_genes),as.integer(not_passing_genes), (passing_genes / total_genes) * 100, (not_passing_genes / total_genes) * 100)
    )
    return(table)
  }
  
  # scatter_plot_variance
  #' @details Here we will build a scatter plot based upon what the user 
  #' inputs via the slider. For plotting the y axis we have taken negative y log
  #' @param counts_matrix the counts matrix that the user inputs 
  #' @param slider The values that the user will input 
  scatter_plot_variance <- function(counts_matrix, slider) {
    county <- counts_matrix[,-1] 
    counts_matrix$variance <- rowVars(as.matrix(county))
    counts_matrix <- counts_matrix[rev(order(counts_matrix$variance)),]
    input_value <- slider/100 
    
    counts_matrix$median <- apply(counts_matrix[,-1], 1, median)
    
    Filter_Method <- (floor(nrow(counts_matrix)*input_value)) > counts_matrix$variance
    
    eruption <- ggplot(data = counts_matrix, aes(x = median, y = -log10(variance)))+
      geom_point(aes(color = Filter_Method)) + 
      theme_bw() + 
      scale_color_manual(values = c('blue', 'red')) + 
      theme(legend.position = "bottom") 
    
    return(eruption)
  }
  
  # scatter_plot_zero
  #' @details Here we will build a scatter plot based upon what the user 
  #' inputs via the slider. For plotting the y axis we have taken negative y log
  #' @param counts_matrix the counts matrix that the user inputs 
  #' @param slider_zero The values that the user will input 
  scatter_plot_zero <- function(counts_matrix, slider_zero) {
    
    counts_matrix$frequency <- rowSums(counts_matrix != 0)
    counts_matrix$median <- apply(counts_matrix[,-1], 1, median)
    
    Filter_Method <- counts_matrix$frequency > slider_zero
    messi <- ggplot(data = counts_matrix, 
                    aes(x = median, y = -log10(frequency)))+
      geom_point(aes(color = Filter_Method )) + 
      theme_bw() + 
      scale_color_manual(values = c('red', 'blue')) + 
      theme(legend.position = "bottom") 
    
    
    return(messi)
    
  }
  
  # heatmap
  #'@details Here we want to produce a heat map post filtering of the counts
  #'matrix that the user has inputted
  #'@param counts_tib This is the counts matrix that the user will input
  #'@param perc_var This is the input from the variance slider 
  #'@param nz_genes This is the input from the zero genes slider 
  #'@param num_colors 
  #'@param palette 
  plot_heatmap <- function(counts_tib, perc_var, nz_genes){
    if (!is.null(input$Normalized_counts_matrix)){
      counts_tib <- as_tibble(counts_tib) %>% mutate(across(starts_with(c("C_", "H_")) & where(is.numeric), ~na_if(.x, 0)))
      counts_tib$no_zeros <- rowSums(is.na(counts_tib))  #make new col, with counts.
      counts_tib <- filter(counts_tib, no_zeros <= nz_genes)
      counts_tib <- log10(counts_tib[,!colnames(counts_tib) %in% c("gene", "no_zeros")]) #exclude the gene names column and log scale the values  
      #produce plot_tib
      plot_tib <- counts_tib %>% 
        mutate(variance = apply(counts_tib, MARGIN = 1, FUN = var)) #compute variance to filter the data
      perc_val <- quantile(plot_tib$variance, probs = perc_var/100, na.rm = TRUE)   #calculate percentile
      plot_tib <- filter(plot_tib, variance >= perc_val) #filter the tibble
      hmap <- heatmap.2(as.matrix(plot_tib[-ncol(plot_tib)]), scale = "row", col = brewer.pal(9, "YlOrRd"))
      return(hmap)}
    else{return(NULL)}
  }
  
  # pca
  #'@details Here we will construct a PC1 vs PC2 plot based on the input sequence 
  #'that the user will utilize 
  #'@param data The counts matrix that the user will input 
  #pca <- function(data){
  #  data <- data[, -c(1)]
  #  pca_plotting <- prcomp(data, scale. = TRUE)
  #  plot_project <- autoplot(pca_plotting, data = data, colour = "darkblue")
  #  return(plot_project)
  # }
  pca <- function(counts_tib, var_slider, comp1, comp2){
    if (!is.null(input$Normalized_counts_matrix)){
      #make plot tib-
      filt_tib <- counts_tib %>% 
        mutate(variance = apply(counts_tib[-1], MARGIN = 1, FUN = var), .after = gene) #calculate variance for filtering
      perc_val <- quantile(filt_tib$variance, probs = var_slider/100, na.rm = TRUE)   #calculate percentile
      filt_tib <- filter(filt_tib, variance >= perc_val) #filter the tibble
      final_tib <- filt_tib[,-c(1,2)]
      pca_res <- prcomp(t(final_tib), center = FALSE, .scale = TRUE) #transpose the data and perform PCA
      pca_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)
      x <- round(pca_var[as.integer(str_sub(comp1, 3))]*100, 2)
      y <- round(pca_var[as.integer(str_sub(comp2, 3))]*100, 2)
      #produce PCA plot
      plot_tib <- tibble(PC1 = pca_res$x[,comp1], PC2=pca_res$x[,comp2])
      pca <- ggplot(plot_tib, aes(PC1, PC2))+
        geom_point()+
        labs(title="Princple Component Analysis Plot")+
        xlab(str_c(comp1, x, "% Variance", sep=" "))+
        ylab(str_c(comp2, y, "% Variance", sep=" "))+
        theme_bw()
      return(pca)}
    else{return(NULL)}
  }
  
  
  #info table output
  observeEvent(list(input$var_slider,input$zero_slider),{
    output$counts_table <-  renderTable({
      data <- counts_data()
      return(filtering(data,input$var_slider,input$zero_slider))
    })
  })
  
  
  # variance plot output
  observeEvent(input$var_slider,{
    output$variance_plot <- renderPlot( {
      data <- counts_data()
      if (is.null(data)) {
        return()
      }
      isolate({
        scatter_plot_variance(data, input$var_slider)
      })
    }) 
  })
  
  # zero plot output
  observeEvent(input$zero_slider,{
    output$zero_plot <- renderPlot({
      data <- counts_data()
      if (is.null(data)) {
        return()
      }
      isolate({
        scatter_plot_zero(data, input$zero_slider)
      })
    }) 
  })
  
  #heatmap output
  observeEvent(input$zero_slider,{
    output$heatmap <- renderPlot({
      req(input$Normalized_counts_matrix)
      c_matrix <- counts_data()
      hot <- plot_heatmap(c_matrix, input$var_slider, input$zero_slider) 
      return(hot)
      
    },height = 600, width = 800)
  })
  
  # pca output
  observeEvent(list(input$var_slider, input$comp1, input$comp2),{
    output$counts_pca <- renderPlot({
      data <- counts_data()
      if (is.null(data)) {
        return()
      }
      isolate({
        pca(data, input$var_slider, input$comp1, input$comp2)
      })
    })
  })
  
  
  #----------------------------------------------------------------------------------------------------------------------------
  
  #######################Tab 03: Differential expression#######################
  
  # reading the differential expression data 
  #' Load in the DE Data
  #' 
  #' @details Here all we want to do is simply load in the DE .csv file that the user will 
  #' input via the "Submit" button.
  DE_data <- reactive({
    req(input$DE_Anaylsis)
    data <- read.csv(input$DE_Anaylsis$datapath)
    colnames(data)[1] <- "Gene"
    return(data)
  })
  
  #' Volcano plot
  #'
  #' @param dataf The loaded data frame.
  #' @param x_name The column name to plot on the x-axis
  #' @param y_name The column name to plot on the y-axis
  #' @param de_base One of the colors for the points.
  #' @param de_highlight The other colors for the points. Hexadecimal strings: "#CDC4B5"
  #'
  #' @return A ggplot object of a volcano plot
  volcano_plot <-
    function(dataf, x_name, y_name, slider, color1, color2) {
      if(is.null(dataf))     
        return(NULL)
      dataf<-na.omit(dataf)
      out<- ggplot(dataf,aes(x= !!sym(x_name),y= -log10(!!sym(y_name)),color= !!sym(y_name)<(10^slider)))+ 
        theme_bw()+
        theme(legend.position="bottom")+
        ggtitle('Volcano plot')+
        scale_color_manual(name= paste0("padj < 1 x 10^",toString(slider)), values=c(color1, color2))+
        geom_point()
      return(out)
    }
  
  # Using renderPlot to display the plot
  observeEvent(list(input$x_axis, input$y_axis, input$padjusted, input$de_base, input$de_highlight), {
    output$volcano <- renderPlot( {
      data <- DE_data()
      #data <- data %>% mutate(!!paste("neg_log10", input$y_axis, sep = "_") := -log10(!!sym(input$y_axis)))
      if (is.null(data)) {
        return()
      }
      isolate({
        volcano_plot(data, input$x_axis, input$y_axis,input$padjusted, input$de_base, input$de_highlight)
      })
    }) 
  })
  
  # Using renderDataTable to display the table
  output$table <- DT::renderDataTable({ 
    data <- DE_data()
    if (is.null(data)) {
      return()
    }
    isolate({
      return(data)
    })
  })
  
  
  #----------------------------------------------------------------------------------------------------------------------------
  #######################Tab 04: GSEA#######################
  
  ensembl_id_to_gene_symbol <- function(values) {
    mart <- useMart('ensembl', dataset='hsapiens_gene_ensembl')
    hgnc_ids <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"), filter="ensembl_gene_id", value=values, mart=mart, uniqueRows=FALSE)
    #hgnc_ids <- getBM(attributes=c("ensembl_gene_id_version","hgnc_symbol"), filter="ensembl_gene_id_version", value=values, mart=mart, uniqueRows=FALSE)
    return(hgnc_ids)
  }
  
  
  #Run fgsea using a ranked list of descending log2FC against the C2 canonical
  #pathways gene set
  #Set minsize to 15 and maxsize to 500, leave the other parameters as defaults
  #fgsea_results <- run_gsea(labeled_results, 'c2.cp.v7.5.1.symbols.gmt', 15, 500)
  #fgsea_results
  #' Function to run fgsea on DESeq2 results
  #'
  #' @param labeled_results (tibble): the labeled results from DESeq2
  #' @param gmt (str): the path to the GMT file
  #' @param min_size: the threshold for minimum size of the gene set
  #' @param max_size: the threshold for maximum size of the gene set
  #'
  #' @return tibble containing the results from running fgsea using descending
  #' log2foldchange as a ranking metric
  #' @export
  #'
  #' @examples fgsea_results <- run_gsea(labeled_results, 'c2.cp.v7.5.1.symbols.gmt', 15, 500)
  run_gsea <- function(labeled_results, gmt, min_size, max_size) {
    #run_gsea <- function(labeled_results, 'c2.cp.v7.5.1.symbols.gmt', 15, 500) {
    # Strip the digit extensions in the gene names in labled_results and pull the 
    # updated/modified gene names to be used in getLDS to retrieve information from
    # human/mouse linked database
    labeled_results_modified <- separate(labeled_results, genes, sep='\\.', into='genes', remove=TRUE)
    # Pull the modified gene names
    #labeled_results_modified <- labeled_results
    genes <- labeled_results_modified$genes
    # Connect to BioMart and retrieve information from human/mouse linked database
    human <- useMart('ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host = "https://dec2021.archive.ensembl.org/")
    mouse <- useMart('ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host = "https://dec2021.archive.ensembl.org/")
    
    hgnc_symbols <- getLDS(attributes=c('ensembl_gene_id'), filters='ensembl_gene_id', values=genes, mart=mouse, 
                           attributesL=c('hgnc_symbol'), martL=human, uniqueRows=T)
    
    hgnc_results <- left_join(labeled_results_modified, hgnc_symbols, by=c('genes' = 'Gene.stable.ID'))
    # Create ranked list of log2FoldChange values
    log2FC <- filter(hgnc_results, !is.na(HGNC.symbol) & !is.na(log2FoldChange))
    log2FC <- log2FC %>% arrange(desc(log2FoldChange))
    log2FC_rank <- log2FC[c("HGNC.symbol", "log2FoldChange")]
    log2FoldChange_rank <- deframe(log2FC_rank)
    
    # Get C2 Canonical Pathways gene set collection
    gmt_pathways <- gmtPathways(gmt)
    # Run gsea on ranked log2FoldChange values
    fgsea_results <- fgsea(gmt_pathways, log2FoldChange_rank, minSize=min_size, maxSize=max_size)
    fgsea_results <- as_tibble(fgsea_results)
    #rename first col to 'gene'
    names(fgsea_results)[1] <- 'genes'
    return(fgsea_results)
  }
  
  create_fgsea <- function(DE_file_df) {
    deseq2_res <- DE_file_df
    padj_threshold <- input$DE_slider
    if (is.null(padj_threshold)) {
      padj_threshold <- 0.10
    }
    #rename first col to 'gene'
    names(deseq2_res)[1] <- 'genes'
    
    labeled <- deseq2_res %>%
      #as_tibble(rownames='genes') %>%
      mutate(volc_plot_status = case_when(log2FoldChange > 0 & padj < padj_threshold ~ 'UP', 
                                          log2FoldChange < 0 & padj < padj_threshold ~ 'DOWN', 
                                          TRUE ~ 'NS'))
    #Display the summary of the tibble
    labeled_results <- labeled[order(labeled$padj, decreasing = FALSE), ]
    labeled_results %>% relocate(genes, volc_plot_status, log2FoldChange, padj)
    fgsea_results <- run_gsea(labeled_results, 'c2.cp.v7.5.1.symbols.gmt', 15, 500)
    col_type <- sapply(fgsea_results, class)
    fgsea_res <- fgsea_results[, sapply(fgsea_results, class) != "list"]
    fgsea_res_mat <- as_tibble(fgsea_res)
  }
  
  #' Function to plot top ten positive NES and top ten negative NES pathways
  #' in a barchart
  #'
  #' @param fgsea_results (tibble): the fgsea results in tibble format returned by
  #'   the previous function
  #' @param num_paths (int): the number of pathways for each direction (top or
  #'   down) to include in the plot. Set this at 10.
  #'
  #' @return ggplot with a barchart showing the top twenty pathways ranked by positive
  #' and negative NES
  #' @export
  #'
  #' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
  top_pathways <- function(fgsea_results, num_paths){
    
    #rename first col to 'gene'
    names(fgsea_results)[1] <- 'pathway'
    
    #print("top_pathways")
    #fgsea_results is already sorted in the aescending order of NES
    # make a tibble with only pathways with top 10 positive and top 10 negative NES values
    # gather top 10 negative values
    top_neg <- fgsea_results[1:num_paths,]
    
    # gather top 10 positive values
    start_pos <- nrow(fgsea_results) - (num_paths-1)
    stop_pos  <- nrow(fgsea_results)
    top_pos <- fgsea_results[start_pos:stop_pos,]
    
    # filter out the rows wtih top 10 pos and neg valued pathways
    pos_pathway <- top_pos$pathway
    neg_pathway <- top_neg$pathway
    fgsea_subset <- fgsea_results %>% filter(pathway %in% c(pos_pathway, neg_pathway)) 
    
    # In fgsea_subset create a new col of pathway names (x_axis_names) without '-' 
    # in the name and reorder them for x-axis lables 
    #fgsea_subset <- fgsea_subset %>% mutate(x_axis_names = str_replace_all(pathway, '_', ' ')) %>% 
    fgsea_subset <- fgsea_subset %>% mutate(x_axis_names = str_replace_all(pathway, '_', ' ')) 
    
    # Create a bar plot of top 10 pos and neg NES pathways
    nes_plot <- 
      ggplot(fgsea_subset) +
      geom_bar(aes(x=x_axis_names, y=NES, fill = NES > 0), stat='identity', show.legend = FALSE) +
      scale_fill_manual(values = c('TRUE' = 'blue', 'FALSE' = 'red')) + 
      #theme_minimal(base_size = 8) +
      theme_bw(base_size = 8) +
      ggtitle('fgsea results for Hallmark MSigDB gene sets') +
      ylab('Normalized Enrichment Score (NES)') +
      xlab('') 
    
    # Since some of the the pathway names are very long let us wrap the names so that we can
    # display the bar plot correctly
    wrap_pathway_name <- function(x) str_wrap(x, width = 50)
    nes_plot <- nes_plot +
      scale_x_discrete(labels = wrap_pathway_name) 
    
    # Filp the plot from vertical to horizontal layout
    nes_plot <- nes_plot + coord_flip()
    nes_plot
    return(nes_plot)
    
  }
  
  load_fgsea_data <- reactive({
    fgsea_file <- input$fgsea_file
    if (is.null(input$fgsea_file)) {
      return(NULL)} 
    fgsea_file_df <- read_tsv(fgsea_file$datapath, show_col_types = FALSE)
    
    if ("log2FoldChange" %in% colnames(fgsea_file_df)) {
      print("This apprears to be a DE file.")
      print("Creating fgsea.csv file now. This will take a while. Please wait.")
      #########################################################################
      # Following create_fgsea function is used to create a fgsea file from the
      # de expression. Run it only if we need to create a new fgsea file.
      fgsea_res_mat <- create_fgsea(fgsea_file_df)
      #########################################################################
      write.csv(fgsea_res_mat, 'fgsea.csv', row.names=FALSE)
      print("Finished creating 'fgsea.csv' file.")
      return(fgsea_res_mat)
    }
    print("This apprears to be FGSEA file. Processing now...")
    fgsea_file_df <- read.csv(fgsea_file$datapath, header = TRUE)
    return(fgsea_file_df)
  })
  
  output$GSEA_Plot1 <- renderPlot({
    req(input$fgsea_file)
    fgsea_results <- load_fgsea_data()
    view(fgsea_results)
    if(is.null(fgsea_results)){return (NULL)}
    
    #Display the results of the fgsea in a tibble sorted by by NES (ascending)
    filtered_results <- fgsea_results[order(fgsea_results$NES, decreasing = FALSE), ]
    filtered_results
    
    #Plot the top ten pathways with both positive and negative NES (20 total)
    #Color the pathways by the sign of their NES (positive or negative)
    fgsea_plot <- top_pathways(filtered_results, 10)
    fgsea_plot
  })
  
  # }) #observeEvent("Upload FGSEA File"
  
  observeEvent("Update Top Pathways", {
    val <- input$GSEA_slider1
    updateSliderInput(session, "GSEA_slider1", value = val)
    
    output$GSEA_Plot1 <- renderPlot({
      fgsea_results <- load_fgsea_data()
      if(is.null(fgsea_results)){return (NULL)}
      
      # filter rows with padj < padj_threshold
      pathways_threshold <- input$GSEA_slider1
      #Display the results of the fgsea in a tibble sorted by by NES (ascending)
      #filtered_results <- fgsea_res[order(fgsea_res$NES, decreasing = FALSE), ]
      filtered_results <- fgsea_results[order(fgsea_results$padj, decreasing = FALSE), ]
      filtered_results
      
      #Plot the top ten pathways with both positive and negative NES (20 total)
      #Color the pathways by the sign of their NES (positive or negative)
      fgsea_plot <- top_pathways(filtered_results, pathways_threshold)
      fgsea_plot
    },  height = 400, width = 400)
    
  })#observeEvent("Update Top Pathways")
  
  
  observeEvent("Update Table", {
    slider2_val <- input$GSEA_slider2
    updateSliderInput(session, "GSEA_slider2", value = slider2_val)
    
    output$GSEA_Table2 <- DT::renderDataTable({
      display_table <- sliderValues()
      return(display_table)
    }) 
  })
  
  output$GSEA_Table2 <- DT::renderDataTable({
    display_table <- sliderValues()
    return(display_table)
  }) 
  # Reactive function
  sliderValues <- reactive({
    fgsea_table <- load_fgsea_data()
    padj_threshold <- input$GSEA_slider2
    fgsea_res <- filter(fgsea_table, padj <= padj_threshold)
    
    pathways_option <- input$GSEA_radio2
    NES_filter <- fgsea_res
    fgsea_res
    if (pathways_option == "Positive") {
      NES_filter <- filter(NES_filter, NES >= 0)
      #print("Positive Pathways Only")
    }
    if (pathways_option == "Negative") {
      NES_filter <- filter(NES_filter, NES < 0)
      #print("Negative Pathways Only")
    }
    NES_filter
    NES_filter_mat <- as_tibble(NES_filter)
    return(NES_filter_mat)
  })
  
  output$fgsea_NES_Table <- downloadHandler(
    filename = function() {
      paste("filtered_fgsea-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(display_table <- sliderValues(), file)
    }
  )
  
  observeEvent("Update Padj Threshold", {
    slider_val <- input$GSEA_slider3
    updateSliderInput(session, "GSEA_slider3", value = slider_val)
    
    output$GSEA_Plot3 <- renderPlot({
      NES_plot <- NES_Padj_Plot()
      return(NES_plot)
    }, height= 400, width=500) 
    
  })#observeEvent("Update Padj Threshold")
  
  NES_Padj_Plot <- reactive({
    fgsea_table <- load_fgsea_data()
    padj_threshold <- input$GSEA_slider3
    fgsea_res_below <- filter(fgsea_table, padj <= padj_threshold)
    fgsea_res_above <- filter(fgsea_table, padj > padj_threshold)
    
    fgsea_res_mod <- fgsea_table %>%
      mutate(slider_status = case_when(padj <= padj_threshold ~ "Below Threshold", 
                                       padj > padj_threshold ~ "Above Threshold"))
    # Now to actual plotting
    nes_padj_plot <-  
      ggplot(fgsea_res_mod, aes(x=NES, y=-log10(padj), color=slider_status)) +
      geom_point() +
      theme_bw() +
      scale_color_manual(name='Padj Threshold', values = c('Below Threshold' = 'gray', 'Above Threshold' = 'red')) + 
      xlab("NES") + ylab(paste0('-log10(padj)')) +
      labs(title= 'Scatter plot of NES vs -log10 adjusted p-value',
           subtitle= paste0('Genes Below Threshold:', nrow(fgsea_res_below),
                            ' Genes Above Threshold:', nrow(fgsea_res_above)))
    
    nes_padj_plot
    return(nes_padj_plot)
    
  })
  
  NES_Padj_Plot2 <- function(dataf, x_name, y_name, slider, color1, color2) {
    if (is.null(x_name) || is.null(y_name)) {return(NULL)}
    if (x_name == y_name) {return(NULL)}
    if (is.null(input$DE_file)) {return(NULL)}
    
    ordered_df <- dataf %>%  dplyr::arrange(dplyr::desc(!!rlang::sym(x_name)),
                                            dplyr::desc(!!rlang::sym(y_name)))
    # Since we do not know which dataf cols are x_name and y_name, we will re-arrange
    # the x_name col as the 1st col and y_name as the 2nd col in the new data frame
    # so that it will be easy to do volcano plot.
    otherCols <-setdiff(colnames(ordered_df), unique(c(x_name,y_name)))
    #Re-arrange columns according to the selections
    ordered_df <- ordered_df %>%  dplyr::select(!!rlang::sym(x_name),
                                                !!rlang::sym(y_name),
                                                !!otherCols)
    
    # Add a new col slider_cond whose values will be TRUE, FALSE or NA  and will
    # used to color code the volcano
    # FALSE if y-axis value >= slider_factor, 
    # TRUE if y-axis value < slider_factor, else NA 
    
    slider_factor <- (1 * (10^slider))
    ordered_df <- ordered_df %>%
      mutate(slider_cond = case_when(ordered_df[, 2] < slider_factor ~ "TRUE", 
                                     ordered_df[, 2] >= slider_factor ~ "FALSE", TRUE ~ 'NA'))
    
    x_axis <- colnames(ordered_df[,1])  #x-axis label for the plot
    y_axis <- colnames(ordered_df[,2])  #y-axis label for the plot
    
    # Remove rows with any NA values in plotting cols:
    ordered_df <- filter(ordered_df, !is.na(ordered_df[,1]))
    ordered_df <- filter(ordered_df, !is.na(ordered_df[,2]))
    
    # Now to actual plotting
    vol_plot <-  
      ggplot(ordered_df, aes(x=ordered_df[,1], y=-log10(ordered_df[,2]), color=slider_cond)) +
      geom_point() + 
      scale_color_manual(name=paste0(y_name,"< 1 * 10^",slider), values = c('FALSE' = color1, 'TRUE' = color2, 'NA'='black')) + 
      xlab(x_name) + ylab(paste0('-log10(',y_name,')')) +
      labs(title=paste0('Volcano Plot: ', x_name, '  VS  -log10(', y_name, ')' )) +
      theme(legend.position="bottom")
    
    vol_plot
    return(vol_plot)
  }
  
}

# Run the application 
shinyApp(ui = ui, server = server)
