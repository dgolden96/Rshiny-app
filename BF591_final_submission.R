#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(GEOquery)
library(DT)
library(ggplot2)
library(shinythemes)
library(dashboard)
library(shinydashboard)
library(reshape)
library(colourpicker)
library(tidyverse)
library(ggbeeswarm)

# Define UI for application that draws a histogram
ui <- navbarPage("BF591 Final App",
                 theme = shinytheme("sandstone"),
                 tabPanel("Metadata Summary", 
                          sidebarLayout(
                              sidebarPanel(
                                  fileInput("meta", "Upload metadata in CSV format.", accept = ".csv"),
                                  sliderInput("colnum", "How many columns of the sample information table would you like to see?", 2, 69, 5),
                              ),
                              mainPanel(tabsetPanel(
                                  tabPanel("Summary Table", DT::dataTableOutput("summary_tb")),
                                  tabPanel("Sample Information Table", DT::dataTableOutput("sample_info_tb")),
                                  tabPanel("Example Metadata Plot", plotOutput("sample_plot"), plotOutput("sample_plot2"))
                              ))
                          )
                 ),
                 tabPanel("Counts Matrix Analysis",
                          sidebarLayout(
                              sidebarPanel(
                                  fileInput("counts", "Upload counts matrix in CSV format.", accept = ".csv"),
                                  sliderInput("percentile", "Select the minimum percentile of variance that a gene must meet to be included.", 0, 1.0, 0.5),
                                  sliderInput("non_na_minimum", "Select the minimum number of non-zero samples that a gene must have to be included.", 0, 69, 60),
                                  sliderInput("PC_num1", "Select which principal component to plot on the X-axis.", 1, 69, 1),
                                  sliderInput("PC_num2", "Select which principal component to plot on the Y-axis.", 1, 69, 2)
                              ),
                              mainPanel(tabsetPanel(
                                  tabPanel("Counts Summary", DT::dataTableOutput("filter_tb")),
                                  tabPanel("Diagnostic Plots", plotOutput("median_plot"), plotOutput("zero_plot")),
                                  tabPanel("Heatmap", "Warning: Heatmap may only load legibily under strict filtering criteria", plotOutput("heatmap")),
                                  tabPanel("PCA", plotOutput("PCA"))
                              ))
                          )
                 ),
                 tabPanel("Differential Gene Expression Analysis",
                          sidebarLayout(
                              sidebarPanel(
                                  fileInput("deseq", "Upload matrix of DESeq2 results in CSV format.", accept = ".csv"),
                                  radioButtons('variable1', 'Select the field to plot on the X-axis.', choices=c('baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')),
                                  radioButtons('variable2', 'Select the field to plot on the Y-axis.', choices=c('baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')),
                                  colourInput('color1', 'Choose the first color for your plot.', value = 'blue'),
                                  colourInput('color2', 'Choose the second color for your plot.', value = 'red'),
                                  sliderInput('slidervalue', 'Select a slider value.', -300, -1, -100),
                                  actionButton("volcano_button", "Volcano!", text="GO!")
                              ),
                              mainPanel(tabsetPanel(
                                  tabPanel("Results Table", DT::dataTableOutput("deseq_tb")),
                                  tabPanel("Volcano Plot", plotOutput('volc_plot'), DT::dataTableOutput("volc_table"))
                              ))
                          )
                 ),
                 tabPanel("Visualization and Individual Gene Expression",
                          sidebarLayout(
                              sidebarPanel(
                                  fileInput("meta_again", "Upload metadata in CSV format.", accept = ".csv"),
                                  fileInput("counts_again", "Upload counts matrix in CSV format.", accept = ".csv"),
                                  textInput("gene_select", "Choose a gene.", ""), 
                                  radioButtons("plot_select", "Choose a plot type.", choices=c("Bar plot", "Box plot", "Violin plot", "Beeswarm plot"))
                              ),
                              mainPanel(tabsetPanel(
                                  tabPanel("Categorical Plots", plotOutput("cat_plot")),
                                  tabPanel("Numerical Plot", plotOutput("num_plot"))
                              ))
                          )
                 )
)

# define server logic
server <- function(input, output) {
    
    # submit input button was screwing up reactivity of other plots
    observeEvent(input$volcano_button, {
        output$volc_plot <- renderPlot({
            req(input$deseq)
            volcano_plot(load_deseq(), input$variable1, input$variable2, input$slidervalue, input$color1, input$color2)
        })   
    })
    
    # increase max upload size
    options(shiny.maxRequestSize=30*1024^2)
    
    # the following are functions to load input files
    load_meta <- reactive({
        req(input$meta)
        metadata_df <- read.csv(input$meta$datapath, row.names=1)
        return(metadata_df)
    })
    
    load_counts <- reactive({
        req(input$counts)
        counts_df <- read.csv(input$counts$datapath, row.names=1) 
        return(counts_df)
    })
    
    load_deseq <- reactive({
        req(input$deseq)
        deseq_df <- read.csv(input$deseq$datapath, row.names=1,header=TRUE)
        return(deseq_df)
    })
    
    load_counts_again <- reactive({
        req(input$counts_again)
        counts_again_df <- read.csv(input$counts_again$datapath, row.names=1,header=TRUE)
        return(counts_again_df)
    })
    
    load_meta_again <- reactive({
        req(input$meta_again)
        meta_again_df <- read.csv(input$meta_again$datapath, row.names=1,header=TRUE)
        return(meta_again_df)
    })

    # generate summary table for counts summary tab
    # didn't realize until later that the numeric versions of RIN, PMI, and age fields are actually in the table already
    make_summary <- reactive ({
        req(input$meta)
        metadata_df <- load_meta()
        # section between these comments is for generating numerics from character fields containing numerical data
        RIN <- as.numeric(substring(metadata_df$characteristics_ch1.4, first=6))
        PMI <- as.numeric(substring(metadata_df$characteristics_ch1.2, first=6))
        diag <- substring(metadata_df$characteristics_ch1.1, first = 12)
        age <- as.numeric(substring(metadata_df$characteristics_ch1.3, first = 15))
        metadata_df$age_numeric <- age
        metadata_df$RIN_numeric <- RIN
        metadata_df$PMI_numeric <- PMI
        metadata_df$diag <- diag
        RIN<-RIN[!is.na(RIN)]
        PMI<-PMI[!is.na(PMI)]
        diag<-diag[!is.na(diag)]
        age<-age[!is.na(age)]
        mean_RIN <- format(round(mean(RIN), 2), nsmall = 2) 
        sd_RIN <- format(round(sd(RIN), 2), nsmall = 2) 
        mean_PMI <- format(round(mean(PMI), 2), nsmall = 2) 
        sd_PMI <- format(round(sd(PMI), 2), nsmall = 2) 
        mean_age <- format(round(mean(age), 2), nsmall = 2) 
        sd_age <- format(round(sd(age), 2), nsmall = 2)
        unique_diag <- unique(diag)
        # end of previous comment
        RIN_string <- paste(as.character(mean_RIN), as.character(sd_RIN), sep=" +/- ", collapse=NULL)
        PMI_string <- paste(as.character(mean_PMI), as.character(sd_PMI), sep=" +/- ", collapse=NULL)
        age_string <- paste(as.character(mean_age), as.character(sd_age), sep=" +/- ", collapse=NULL)
        distinct_diag <- paste(unique_diag[1], unique_diag[2], sep=", ", collapse=NULL)
        summary_df <- data.frame("Data Type" = c(class(diag), class(age), class(RIN), class(PMI)),
                                 "Mean (sd) or Distinct Values" = c(distinct_diag, age_string, RIN_string, PMI_string), 
                                 row.names = c("Donor Diagnosis", "Donor Age at Death", "RIN", "PMI"))
        names(summary_df) <- c("Data Type", "Mean (sd) or Distinct Values")
        return(summary_df)
    })
    
    # generate violin plot from metadata
    make_violin <- reactive ({
        req(input$meta)
        metadata_df <- load_meta()
        RIN <- as.numeric(substring(metadata_df$characteristics_ch1.4, first=6))
        PMI <- as.numeric(substring(metadata_df$characteristics_ch1.2, first=6))
        diag <- substring(metadata_df$characteristics_ch1.1, first = 12)
        age <- as.numeric(substring(metadata_df$characteristics_ch1.3, first = 15))
        metadata_df$age_numeric <- age
        metadata_df$RIN_numeric <- RIN
        metadata_df$PMI_numeric <- PMI
        metadata_df$diag <- diag
        plot <- ggplot(metadata_df, aes(x=diag, y=age_numeric)) +
            geom_violin() + 
            theme_gray() +
            geom_boxplot(width = 0.2) +
            scale_y_continuous(trans='log10') +
            xlab("Condition") +
            ylab("Age at Death") +
            ggtitle("Effect of Huntington's Disease on Life Expectancy") +
            theme(plot.title = element_text(hjust = 0.5)) 
        return(plot)
    })
    
    # generate density plot from metadata
    make_dplot <- reactive ({
        req(input$meta)
        metadata_df <- load_meta()
        plot <- ggplot(metadata_df, aes(x=age.of.death.ch1)) + 
            geom_density() +
            theme_gray() +
            xlab("Age at Death") +
            ylab("Count") +
            ggtitle("Density Plot of Donor Age at Death") +
            theme(plot.title = element_text(hjust = 0.5)) 
        return(plot)
            
    })
        
    # DT is a great package for Rshiny tables, lets you sort by field
    output$summary_tb <- DT::renderDataTable({
        req(input$meta)
        load_meta()
        make_summary()
    })
    
    output$sample_info_tb <- DT::renderDataTable({
        sample <- load_meta()
        sample <- dplyr::select(sample, 1:input$colnum)
        DT::datatable(sample, options = list(lengthMenu = c(5, 30, 50), pageLength = 5))
    })
    
    output$sample_plot <- renderPlot({
        req(input$meta)
        load_meta()
        make_violin()
    })
    
    output$sample_plot2 <- renderPlot({
        req(input$meta)
        load_meta()
        make_dplot()
    })

    # function that filters counts matrix by row-wise variance and row-wise zero counts
    filter_counts <- reactive ({
        req(input$counts)
        counts_df <- load_counts()
        counts_matrix <- as.matrix(counts_df)
        quants <- quantile(counts_matrix, probs = c(0.0, input$percentile))
        quants_numeric <- unname(quants)
        variance_to_beat <- quants_numeric[2]
        counts_df_copy <- counts_df
        counts_df$na_count <- apply(is.na(counts_df_copy), 1, sum)
        counts_df$variance <- apply(counts_df_copy, 1, var, na.rm=TRUE)
        filter_df <- dplyr::filter(counts_df, (ncol(counts_df) - na_count) >= input$non_na_minimum)
        filter_df <- dplyr::filter(filter_df, variance > variance_to_beat)
        return(filter_df)
    })
    
    # generate counts summary table using filter_counts and modify from there
    counts_summary <- reactive ({
        req(input$counts)
        counts_df <- load_counts()
        filter_df <- filter_counts()
        sample_number <- ncol(counts_df)
        gene_number <- nrow(counts_df)
        filter_gene_number <- nrow(filter_df)
        fail_gene_number <- gene_number - filter_gene_number
        filter_percent <- as.numeric((filter_gene_number / gene_number) * 100)
        filter_percent <- round(filter_percent, 2)
        fail_percent <- as.numeric(100 - filter_percent)
        fail_percent <- round(fail_percent, 2)
        counts_summary_df <- data.frame("Total Samples" = sample_number, "Total Genes" = gene_number, "Number of Genes Passing Filter" = filter_gene_number,
                                        "Percent of Genes Passing Filter" = filter_percent, "Number of Genes Failing Filter" = fail_gene_number,
                                        "Percent of Genes Failing Filter" = fail_percent, row.names = "Summary Values")
        colnames(counts_summary_df) <- c("Total Samples", "Total Genes", "Number Passing Filter", "Percent Passing Filter",
                                         "Number Failing Filter", "Percent Failing Filter")
        return(counts_summary_df)
    })
        
    output$filter_tb <- DT::renderDataTable({
        counts_summary()        
    })
    
    # generate our median count scatter plot for diagnostic plots tab
    make_diag_scatter_median <- reactive ({
        req(input$counts)
        counts_df <- load_counts()
        filter_df <- filter_counts()
        counts_df$variance <- apply(counts_df, 1, var, na.rm=TRUE)
        # need to subset dataframes because I added the variance field
        counts_df$row_median = apply(counts_df[,1:69], 1, median)
        filter_df$row_median = apply(filter_df[,1:69], 1, median)
        filterpass_genes <- row.names(filter_df)
        counts_df <- dplyr::mutate(counts_df, filterpass = (row.names(counts_df) %in% filterpass_genes))
        counts_df$filterpass <- factor(counts_df$filterpass, levels=c(TRUE, FALSE))
        plot <- ggplot(counts_df, aes(x=row_median, y=variance, col=as.factor(filterpass))) + 
            geom_point() +
            theme_gray() +
            scale_y_continuous(trans='log2') +
            xlab("Row-wise Median of Counts") +
            ylab("Row-wise Variance of Counts") +
            ggtitle("Median vs. Variance") +
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_colour_discrete("Passes Filter?") 
        return(plot)
    })
    
    output$median_plot <- renderPlot(
        make_diag_scatter_median()
    )
    
    # generate row-wise zero count scatter plot for diagnostic plots tab
    make_diag_scatter_zero <- reactive ({
        req(input$counts)
        counts_df <- load_counts()
        filter_df <- filter_counts()        
        filter_df$row_zero_count = rowSums(filter_df==0)
        counts_df$row_zero_count = rowSums(counts_df==0)
        counts_df$row_median = apply(counts_df[,1:69], 1, median)
        filter_df$row_median = apply(filter_df[,1:69], 1, median)
        filterpass_genes <- row.names(filter_df)
        counts_df <- dplyr::mutate(counts_df, filterpass = (row.names(counts_df) %in% filterpass_genes))
        counts_df$filterpass <- factor(counts_df$filterpass, levels=c(TRUE, FALSE))
        plot <- ggplot(counts_df, aes(x=row_median, y=row_zero_count, col=as.factor(filterpass))) + 
            geom_point() +
            theme_gray() +
            scale_y_continuous(trans='log2') +
            xlab("Row-wise Median of Counts") +
            ylab("Number of Zeroes per Row") +
            ggtitle("Median vs. Zeroes") +
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_colour_discrete("Passes Filter?") 
        return(plot) 
    })
    
    output$zero_plot <- renderPlot(
        make_diag_scatter_zero()
    )
    
    # generate counts matrix heatmap. GGplot wasn't working very well for this, so I just used the base R heatmap
    make_heatmap <- reactive ({
        req(input$counts)
        counts_df <- load_counts()
        filter_df <- filter_counts()
        filter_matrix <- as.matrix(filter_df[,1:69])
        matrix_melt <- melt(filter_matrix)
        plot <- heatmap(filter_matrix, Rowv = NA, Colv = NA, main = 'Heatmap of Filtered Counts Matrix')
        return(plot)
    })
    
    # tried to make heatmap bigger for better visualization, didn't really help
    output$heatmap <- renderPlot({
        make_heatmap()}, height = 400, width = 600
    )
    
    # function to run PCA on filtered counts matrix
    run_pca <- reactive ({
        req(input$counts)
        req(input$PC_num1)
        req(input$PC_num2)
        PC_1 <- paste0("PC", as.character(input$PC_num1))
        PC_2 <- paste0("PC", as.character(input$PC_num2))
        counts_df <- load_counts()
        filter_df <- filter_counts()
        filter_matrix <- as.matrix(filter_df[,1:69])
        pca <- prcomp(t(filter_matrix), scale = TRUE)
        sdev_vector <- pca$sdev
        PoV <- (sdev_vector**2)/sum(sdev_vector**2)
        pca_df <- data.frame(Sample=rownames(pca$x),
                             X_val=pca$x[,PC_1],
                             Y_val=pca$x[,PC_2])
        pca_plot <- ggplot(pca_df, aes(x=X_val, y=Y_val, label=Sample)) +
            geom_point() +
            theme_gray() +
            ggtitle("Principal Component Analysis of Expression Variance") +
            theme(plot.title = element_text(hjust = 0.5)) +
            xlab(paste0(PC_1, " - ", round(PoV[input$PC_num1] * 100),"% variance")) +
            ylab(paste0(PC_2, " - ", round(PoV[input$PC_num2] * 100),"% variance")) 
        return(pca_plot)
    })
    
    output$PCA <- renderPlot(
        run_pca()
    )
    
    output$deseq_tb <- DT::renderDataTable(
        return(load_deseq())
    )
    
    # function from assignment 7 to generate volcano plot for DESeq results data
    volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
        dataf_tib <- as_tibble(dataf)
        dataf_tib <- dplyr::mutate(dataf_tib, colorchoice = factor(case_when(!!sym(y_name)<1*10**slider ~ 'Lowest p-values', !!sym(y_name)>1*10**slider ~ 'Other p-values')))
        volc_plot <- ggplot(dataf_tib) +
            geom_point(aes(x=!!sym(x_name),y=-log10(!!sym(y_name)), color=colorchoice)) + 
            theme_gray() +
            ggtitle('Volcano Plot Output') +
            theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
            geom_hline(yintercept=(-slider)) +
            scale_color_manual(name="P-adj < 1e-100", values=c('Lowest p-values'=color1, 'Other p-values'=color2))
        return(volc_plot)
    }
        
    # generate table of filtered DESeq2 results data
    draw_volc_table <- function(dataf, slider) {
        df <- dplyr::filter(dataf, dataf$padj<1*10**slider)
        df$pvalue <- format(df$pvalue, digits=6, scientific=TRUE)
        df$padj <- format(df$padj, digits=6, scientific=TRUE)
        return(df)
    }
    
    output$volc_table <- DT::renderDataTable({
        req(input$deseq)
        volc_df <- load_deseq()
        draw_volc_table(volc_df, input$slidervalue)
    })
    
    # following 4 functions are for the plots in individual gene expression tab. They're mostly the same data just represented four ways.
    boxplot <- reactive ({
        counts <- load_counts_again()
        counts_trans <- as.data.frame(t(counts))
        metadata_df <- load_meta_again()
        diag <- substring(metadata_df$characteristics_ch1.1, first = 12)
        counts_trans$diag <- diag
        data <- dplyr::select(counts_trans, input$gene_select, diag)
        names(data) <- c('expr', 'diag')
        plot <- ggplot(data, aes(x=diag, y=expr)) +
            geom_boxplot() +
            theme_gray() +
            xlab("Condition") +
            ylab(paste0("Normalized Counts for ",input$gene_select)) +
            ggtitle("Individual Gene Expression Violin Plot") +
            theme(plot.title = element_text(hjust = 0.5)) 
        return(plot)
    })
    
    barplot <- reactive ({
        counts <- load_counts_again()
        counts_trans <- as.data.frame(t(counts))
        metadata_df <- load_meta_again()
        diag <- substring(metadata_df$characteristics_ch1.1, first = 12)
        counts_trans$diag <- diag
        data <- dplyr::select(counts_trans, input$gene_select, diag)
        names(data) <- c('expr', 'diag')
        plot <- ggplot(data, aes(x=diag, y=expr)) +
            geom_col() + 
            theme_gray() +
            xlab("Condition") +
            ylab(paste0("Normalized Counts for ",input$gene_select)) +
            ggtitle("Individual Gene Expression Violin Plot") +
            theme(plot.title = element_text(hjust = 0.5)) 
        return(plot)
    })
    
    violinplot <- reactive ({
        counts <- load_counts_again()
        counts_trans <- as.data.frame(t(counts))
        metadata_df <- load_meta_again()
        diag <- substring(metadata_df$characteristics_ch1.1, first = 12)
        counts_trans$diag <- diag
        data <- dplyr::select(counts_trans, input$gene_select, diag)
        names(data) <- c('expr', 'diag')
        plot <- ggplot(data, aes(x=diag, y=expr)) +
            geom_violin() +
            theme_gray() +
            xlab("Condition") +
            ylab(paste0("Normalized Counts for ",input$gene_select)) +
            ggtitle("Individual Gene Expression Violin Plot") +
            theme(plot.title = element_text(hjust = 0.5)) 
        return(plot)
    })
    
    beeswarmplot <- reactive ({
        counts <- load_counts_again()
        counts_trans <- as.data.frame(t(counts))
        metadata_df <- load_meta_again()
        diag <- substring(metadata_df$characteristics_ch1.1, first = 12)
        counts_trans$diag <- diag
        data <- dplyr::select(counts_trans, input$gene_select, diag)
        names(data) <- c('expr', 'diag')
        plot <- ggplot(data, aes(x=diag, y=expr)) +
            geom_beeswarm() +
            theme_gray() +
            xlab("Condition") +
            ylab(paste0("Normalized Counts for ",input$gene_select)) +
            ggtitle("Individual Gene Expression Violin Plot") +
            theme(plot.title = element_text(hjust = 0.5)) 
        return(plot)
    })
    
    # need to use some if else if else logic to select plot based on user input
    output$cat_plot <- renderPlot({
        req(input$counts_again)
        req(input$meta_again)
        req(input$gene_select)
        req(input$plot_select)
        if (input$plot_select == "Bar plot") {
            barplot()
        } else if (input$plot_select == "Box plot") {
            boxplot()
        } else if (input$plot_select == "Violin plot") {
            violinplot()
        } else {
            beeswarmplot()
        }
    })
    
    # bonus scatter plot!
    make_num_plot <- reactive ({
        counts <- load_counts_again()
        counts_trans <- as.data.frame(t(counts))
        metadata_df <- load_meta_again()
        age <- as.numeric(substring(metadata_df$characteristics_ch1.3, first = 15))
        counts_trans$age <- age
        data <- dplyr::select(counts_trans, input$gene_select, age)
        names(data) <- c('expr', 'age')
        plot <- ggplot(data, aes(x=age, y=expr)) +
            geom_point() +
            theme_gray() +
            xlab("Age at Death") +
            ylab(paste0("Normalized Counts for ",input$gene_select)) +
            ggtitle("Gene Expression vs. Age at Death") +
            theme(plot.title = element_text(hjust = 0.5)) 
        return(plot)
    })
    
    output$num_plot <- renderPlot({
        req(input$counts_again)
        req(input$meta_again)
        make_num_plot()
    })
}
# Run the application 
shinyApp(ui = ui, server = server)
