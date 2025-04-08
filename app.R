library(shiny)
library(edgeR)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pathview)
library(shinydashboard)
library(data.table)
library(readxl)
library(reshape2)
library(ggpmisc)
library(plotly)
library(dplyr)
library(readr)
library(DT)
library(shinyWidgets)
library(shinyjs)
library(httr)
library(jsonlite)
library(readxl)
library(stringr)
library(DBI)
library(dbplyr)
library(odbc)
library(dotenv)
library(RPostgres)
library(fgsea)
api_key_chatgpt <- Sys.getenv("CHAT_GPT_KEY")


###extra files
genenameskegg<- read.delim("kegg.txt")
genenames<- read.delim("gene_name_list.txt", header = F)
colnames(genenames)<- c("genes","name")
full_index<- read_excel("Full_GS115_index.xlsx")
colnames(full_index)[3]<- "genes"
genenames$genes<- paste0("lcl|",genenames$genes)

full_index$`Kegg ID`<- gsub("ppa:","",full_index$`Kegg ID`)

pathways<- read.delim("pichia_pathways.txt")
pathways$id<- gsub("ppa","", pathways$id)




call_openai_api <- function(question, model = "gpt-4o") {
  url <- "https://api.openai.com/v1/chat/completions"
  headers <- add_headers(
    "Authorization" = paste("Bearer", api_key_chatgpt),
    "Content-Type" = "application/json")
  
  system_msg <- "You are a specialized bioinformatics assistant with expertise in Pichia pastoris (Komagataella phaffii) biology. Analyze gene expression data to identify functional implications, affected pathways, and biological significance. Focus on metabolic impact, stress responses, and cellular adaptations. Provide well-structured, scientifically accurate interpretations. The stakes are evtremely high, several millios of dollars will be assigned to projects based on your analysis of the data"
  
  messages <- list(
    list(role = "system", content = system_msg),
    list(role = "user", content = question)
  )
  
  body <- list(
    model = model,
    messages = messages,
    temperature = 0.7)
  
  response <- POST(url, headers, body = toJSON(body, auto_unbox = TRUE))
  content <- content(response, "parsed")
  
  if(http_error(response)) {
    return(paste("Error:", content$error$message))
  }
  
  return(content$choices[[1]]$message$content)
}

#To maintain context awareness and carry the conversation

call_openai_api_with_history <- function(message_history, model = "gpt-4o") {
  url <- "https://api.openai.com/v1/chat/completions"
  headers <- add_headers(
    "Authorization" = paste("Bearer", api_key_chatgpt),
    "Content-Type" = "application/json")
  
  body <- list(
    model = model,
    messages = message_history,
    temperature = 0.7)
  
  response <- POST(url, headers, body = toJSON(body, auto_unbox = TRUE))
  content <- content(response, "parsed")
  
  if(http_error(response)) {
    return(paste("Error:", content$error$message))
  }
  
  return(content$choices[[1]]$message$content)
}

prepare_gene_rankings <- function(de_data) {
  gene_list <- de_data$logFC
  names(gene_list) <- de_data$genes
  gene_list <- sort(gene_list, decreasing = TRUE)
  return(gene_list)
}

prepare_pichia_pathways <- function() {
  all_pathways <- list()
  unique_pathways <- unique(pathways$id)
  
  for (pathway_id in unique_pathways) {
    pathway_name <- pathways$name[pathways$id == pathway_id]
    if (length(pathway_name) == 0) next
    
    pathway_genes <- genenameskegg$genes[genenameskegg$pathway == paste0("ppa", pathway_id)]
    
    if (length(pathway_genes) >= 5) {
      all_pathways[[as.character(pathway_name)]] <- pathway_genes
    }
  }
  
  return(all_pathways)
}

perform_gsea <- function(ranked_genes, pathways, nperm = 1000) {
  gsea_results <- fgsea(pathways = pathways,
                        stats = ranked_genes,
                        minSize = 5,
                        maxSize = 500,
                        nperm = nperm)
  
  gsea_results <- as.data.frame(gsea_results)
  gsea_results <- gsea_results[order(gsea_results$NES, decreasing = TRUE),]
  
  return(gsea_results)
}
gene_annotation_mapper <- function(gene_names, annotation_type, index_df = full_index) {
  if (annotation_type == "original") {
    return(gene_names)  # Return original names unchanged
  }
  
  # Define which column to use based on user selection
  target_column <- switch(annotation_type,
                          "gs115_mit" = "GS115 MIT Gene name",
                          "cbs7435" = "CBS7435 Gene name",
                          "kegg" = "Kegg Protein Name",
                          "gs115_ghent" = "GS115 Ghent Gene name",
                          "original" = NULL)
  
  if (is.null(target_column)) {
    return(gene_names)
  }
  
  # Try to match against multiple possible ID columns
  id_columns <- c("genes", "Kegg ID", "CBS7435 ID", "GS115 UGhent ID")
  
  # Create a function to look up a gene name
  map_gene <- function(gene) {
    # Try each possible ID column
    for (id_col in id_columns) {
      # Find matching rows
      matches <- index_df[index_df[[id_col]] == gene, ]
      if (nrow(matches) > 0 && !is.na(matches[[target_column]][1]) && 
          matches[[target_column]][1] != "") {
        return(matches[[target_column]][1])
      }
    }
    # Return original if no mapping found
    return(gene)
  }
  
  # Apply the mapping to each gene name
  mapped_names <- sapply(gene_names, map_gene)
  return(mapped_names)
}

dotenv::load_dot_env()
connect_to_snowflake <- function(warehouse, database, schema) {
  # Attempt to connect to Snowflake using provided credentials
  success <- tryCatch({
    myconn <- dbConnect(odbc::odbc(), Sys.getenv("SNOWFLAKE_USER"), role = 'shiny_app_role', PWD = Sys.getenv("SNOWFLAKE_PASSWORD"))
    
    # Set the warehouse, database, and schema to the specified values
    dbExecute(myconn, paste0("USE WAREHOUSE ", toupper(warehouse), ";"))
    dbExecute(myconn, paste0("USE DATABASE ", toupper(database), ";"))
    dbExecute(myconn, paste0("USE SCHEMA ", toupper(schema), ";"))
    
    # Assuming Snowflake session variables are available for retrieval like this
    current_warehouse <- dbGetQuery(myconn, "SELECT CURRENT_WAREHOUSE() AS WAREHOUSE;")
    current_database <- dbGetQuery(myconn, "SELECT CURRENT_DATABASE() AS DATABASE;")
    current_schema <- dbGetQuery(myconn, "SELECT CURRENT_SCHEMA() AS SCHEMA;")
    
    list(
      warehouse = current_warehouse$WAREHOUSE,
      database = current_database$DATABASE,
      schema = current_schema$SCHEMA
    )
  }, error = function(e) {
    message("Error connecting to Snowflake: ", e$message)
    NULL
  })
  
  if (!is.null(success)) {
    print("Connected to Snowflake successfully.")
  }
  
  return(myconn)
}

myconn<- connect_to_snowflake("EARLY_STAGE","NGS","GENOME_ANALYSIS")

#sample_names<- dbGetQuery(myconn, "SELECT DISTINCT EXPERIMENT_NAME FROM RNA_SEQ_READ_COUNTS")

fetch_sample_data <- function(selected_samples) {
  req(selected_samples)
  
  query <- sprintf(
    "SELECT * FROM RNA_SEQ_READ_COUNTS WHERE SAMPLE_NAME IN ('%s')",
    paste(selected_samples, collapse = "','")
  )
  
  df<- dbGetQuery(myconn, query)
  df$SAMPLE_NAME <- factor(df$SAMPLE_NAME, levels = selected_samples)
  df <- df[order(df$SAMPLE_NAME), ]
}



convert_to_dge_list <- function(df, group) {
  req(df)
  wide_df <- dcast(df, NAME ~ SAMPLE_NAME, value.var = "NUMREADS", fill = 0)
  
  mat <- as.matrix(wide_df[, -1, drop = FALSE])
  rownames(mat) <- wide_df$NAME
  colnames(mat) <- colnames(wide_df)[-1]
  
  dge <- DGEList(counts= mat,group = group)
  keep <- filterByExpr(dge)
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- normLibSizes(dge)
  return(dge)
}


#function used to remove selectedInput
selectizeOptions <- function () { 
  list(
    placeholder = 'Please select an option below',
    onInitialize = I('function() { this.setValue(""); }'),
    plugins = list("remove_button")
  )
}


fixname<- function(df){
  df<- plyr::join(df, genenames)
  df<- df %>%
    mutate(name = ifelse(name %in% NA, genes, name))
  #df<- df[!grepl("locus_tag",df$name),]
  df$name<- gsub( "\\[gene=","", df$name)
  df$name<- gsub( "\\]","", df$name)
  df$name<- gsub("\\[locus_tag=","",df$name)
  return(df)
}





within_strain_de_calculator <- function(y_obj, tp_ref, custom_name_index_df, custom_groups) {
  req(y_obj, tp_ref)
  
  df <- list()
  unique_groups <- unique(custom_groups)
  
  # Retrieve the correct numeric group for the selected reference
  tp_ref_group <- custom_name_index_df[custom_name_index_df$CustomName == tp_ref, "Group"]
  
  for (group in unique_groups) {
    if (group != tp_ref_group) {  # Compare correctly
      subset_idx <- which(custom_groups %in% c(tp_ref_group, group))
      
      if (length(unique(custom_groups[subset_idx])) < 2) {
        next  # Skip invalid comparisons
      }
      
      y_subset <- y_obj[, subset_idx]
      binary_group <- factor(ifelse(custom_groups[subset_idx] == tp_ref_group, "baseline", "current"))
      
      if (nlevels(binary_group) < 2) {
        next
      }
      
      design_subset <- model.matrix(~ binary_group)
      y_subset <- estimateDisp(y_subset, design_subset)
      fit <- glmFit(y_subset, design_subset)
      lrt <- glmLRT(fit)
      
      genes <- rownames(lrt$table)
      
      test_group_name <- custom_name_index_df[custom_name_index_df$Group == group, "CustomName"]
      ref_group_name <- custom_name_index_df[custom_name_index_df$Group == tp_ref_group, "CustomName"]
      
      comparison_name <- paste(test_group_name, "vs", ref_group_name)
      
      df[[comparison_name]] <- lrt$table %>%
        mutate(genes = genes) %>%
        mutate(cond = comparison_name) %>%
        mutate(FDR = p.adjust(PValue, method = "BH")) %>%
        select(genes, logFC, FDR, cond)
    }
  }
  
  if (length(df) == 0) {
    return(NULL)
  }
  
  df <- do.call(rbind, lapply(df, function(d) { rownames(d) <- NULL; d }))
  df <- fixname(df) 
  df$regulated <- "No change"
  df$regulated[df$logFC > 1 & df$FDR < 0.05] <- "Up-regulated"
  df$regulated[df$logFC < -1 & df$FDR < 0.05] <- "Down-regulated"
  
  df$cond <- factor(df$cond, levels = unique(df$cond))
  return(df)
}






kegg_map_generator<- function(df,pathway_id){
  df<- df %>%
    filter(!regulated %in% "No change")
  df_merged <- merge(x=df, y=genenameskegg, by.x="genes", by.y="GS115_MIT")[]
  gene_data<- df_merged$logFC  
  names(gene_data)<- df_merged$KEGG
  names(gene_data)<- gsub("ppa:","", names(gene_data))
  try({
    pathview::pathview(gene.data = gene_data, pathway.id = pathway_id, species = "ppa", 
                       gene.idtype = "KEGG", na.col = "transparent", low = "red", high = "green")
  }, silent = F)
  file.remove(paste0("ppa",pathway_id, ".xml"))
}


generate_volcano_tp <- function(df, logfc1, logfc2, annotation_type) {
  df$name <- gene_annotation_mapper(df$name, annotation_type)
  df$cond <- factor(df$cond, levels = unique(df$cond))
  
  # Filter data first
  df_filtered <- df %>%
    filter(abs(logFC) <= logfc2, abs(logFC) >= logfc1)
  
  # Assign colors manually outside ggplot aes
  df_filtered$color <- dplyr::case_when(
    df_filtered$regulated == "Up-regulated"   ~ "#7cae00",
    df_filtered$regulated == "Down-regulated" ~ "#f8766d",
    TRUE                                      ~ "#00bfc4"
  )
  
  tpvolcano <- ggplot(df_filtered,
                      aes(x = logFC, y = -log10(FDR),
                          text = paste0("Gene: ", name,
                                        "<br>logFC: ", round(logFC, 2),
                                        "<br>FDR: ", signif(FDR, 2)))) +
    geom_point(size = 0.5, color = df_filtered$color) +
    facet_wrap(~cond, scales = "free") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "gray40", linewidth = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40", linewidth = 0.5) +
    labs(x = "", y = "") +  # Remove ggplot labels
    theme_bw() +
    theme(
      text = element_text(face = "bold", size = 12),
      legend.position = "none"
    )
  
  ggplotly(tpvolcano, tooltip = "text") %>%
    layout(
      showlegend = FALSE,
      xaxis = list(title = "log₂ Fold Change"),
      yaxis = list(title = "-log₁₀ FDR")
    )
  
}


## For Download only:
generate_volcano_ggplot <- function(df, logfc1, logfc2, annotation_type) {
  df$name <- gene_annotation_mapper(df$name, annotation_type)
  df$cond <- factor(df$cond, levels = unique(df$cond))
  
  df_filtered <- df %>%
    filter(abs(logFC) <= logfc2, abs(logFC) >= logfc1)
  
  df_filtered$color <- dplyr::case_when(
    df_filtered$regulated == "Up-regulated"   ~ "#7cae00",
    df_filtered$regulated == "Down-regulated" ~ "#f8766d",
    TRUE                                      ~ "#00bfc4"
  )
  
  ggplot(df_filtered,
         aes(x = logFC, y = -log10(FDR))) +
    geom_point(size = 0.5, color = df_filtered$color) +
    facet_wrap(~cond, scales = "free") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "gray40", linewidth = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40", linewidth = 0.5) +
    labs(x = expression(Log[2]~"Fold Change"), y = expression("-log"[10]~"FDR")) +
    theme_bw() +
    theme(
      text = element_text(face = "bold", size = 12),
      legend.position = "none"
    )
}



# Define UI


jscode <- "shinyjs.refresh_page = function() { history.go(0); }"
ui <- dashboardPage(
  dashboardHeader(title = "Differential Expression"),
  dashboardSidebar(
    sidebarMenu( id = "sidebar_menu",
      menuItem("Experiment Selector", tabName = "experiment_selector", icon = icon("flask")),
      menuItem("TP Analysis", tabName = "tp", icon = icon("chart-bar")),
      menuItem("Full Index", tabName = "full_index", icon = icon("table"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(
        tabName = "experiment_selector",
        fluidRow(
          column(width = 4,
                 selectInput(
                   inputId = "experiment_selector",
                   label = "Select Experiment(s):",
                   choices = dbGetQuery(myconn, "SELECT DISTINCT EXPERIMENT_NAME FROM RNA_SEQ_READ_COUNTS")$EXPERIMENT_NAME,
                   multiple = TRUE
                 ),
                 uiOutput("all_sample_select"),
                 uiOutput("sample_selector_ui")
                 
          ),
          
          column(width = 8,
                 
                 uiOutput("custom_name_ui"),
                 br(), br(),
                 uiOutput("submit") 
          )
        ),
        
      ),
      tabItem(
        tabName = "tp",
        h2(""),
        fluidRow(
          column(
            width = 4,
            uiOutput("reference_selector"),
          ),
          column(
            width = 4,
            selectInput("annotation_type", 
                        "Select Gene Annotation Type:",
                        choices = c("Original Names" = "original",
                                    "GS115 MIT Gene Name" = "gs115_mit",
                                    "CBS7435 Gene Name" = "cbs7435",
                                    "Kegg Protein Name" = "kegg",
                                    "GS115 Ghent Gene Name" = "gs115_ghent"),
                        selected = "original"),
          ),
          column(
            width = 4,
            uiOutput("samples_to_use"),
          )
        ),
        fluidRow(
          column(
            width = 6,
            uiOutput("calculate_de"),
          ),
        
        ),
        fluidRow(
          tabBox(
            title = "",
            width = 16,
            height = "1000px",
            tabPanel(
              "Volcano plot",
              
              fluidRow(
                column(width = 4, downloadButton("download_volcano_tp", "Download Plot"))
              ),
              fluidRow(
                #column(width = 4, checkboxInput("tpshowlabels", "Label genes?", value = FALSE)),
                column(width = 4, sliderInput(inputId = "tp_volcano_logFC", label = "Adjust Fold change", min = 0, max = 20, value = c(0, 10)))
              ),
              plotlyOutput(
                "tp_volcano_plot",
                height = "600px", width = "1100px"
              )
            ),
            tabPanel(
              "Data table",
              fluidRow(
                column(
                  width = 2,
                  downloadButton("download_tp_de_table", "Download DE table"),
                ),
                column(
                  width = 3,
                  checkboxInput("tp_regulation", label = "Only show significant changes", value = FALSE)
                ),
                column(
                  width = 2,
                  sliderInput(
                    inputId = "tp_logFC",
                    label = "LogFC range",
                    min = -15,
                    max = 15,
                    value = c(-15, 15)
                  )
                ),
                column(
                  width = 4,
                  actionBttn("strain_chatgpt", "Ask Vera about this result")
                )
              ),
              wellPanel(
               
                uiOutput("gpt_response_container")
              ),
              DTOutput("tp_comparison")
            ),
            tabPanel(
              "Kegg pathway viewer",
              box(
                width = 12,
                column(
                  width = 4,
                  selectInput("tp_pathway", label = "Select Pathway", multiple = FALSE, choices = pathways$name)
                ),
                column(width = 4, uiOutput("drop_down_event_tp_kegg"))
              ),
              uiOutput("tp_pathway")
            ),
            tabPanel(
              "GSEA Analysis",
              fluidRow(
                column(width = 12,
                       box(
                         width = 12,
                         actionButton("run_gsea", "Run Gene Set Enrichment Analysis"),
                         downloadButton("download_gsea_results", "Download GSEA Results")
                       ),
                       box(
                         width = 12,
                         selectInput("gsea_contrast", "Select Contrast for Analysis:", choices = NULL)
                       )
                )
              ),
              fluidRow(
                column(width = 6,
                       plotOutput("gsea_enrichment_plot", height = "500px")
                ),
                column(width = 6, 
                       plotOutput("gsea_pathway_plot", height = "500px")
                )
              ),
              fluidRow(
                column(width = 12,
                       DT::dataTableOutput("gsea_results_table")
                )
              )
            )
          )
        )
      ),
      tabItem(
        tabName = "full_index",
        h2(""),
        fluidRow(DTOutput("full_index"))
      )
    )
  )
)


server <- function(input, output, session) {
  observeEvent(input$refresh, {
    js$refresh_page();
  })
  
  output$all_sample_select <- renderUI({
    req(input$experiment_selector)
    checkboxInput(
      inputId = "select_all_samples",
      label = "Select All Samples",
      value = FALSE
    )
  })
  
  output$sample_selector_ui <- renderUI({
    req(input$experiment_selector)
    selectizeInput(
      inputId = "sample_selector",
      label = "Select Sample(s):",
      choices = NULL,  
      multiple = TRUE
    )
  })
  
  observeEvent(input$experiment_selector, {
    req(input$experiment_selector)
    
    query <- sprintf(
      "SELECT DISTINCT SAMPLE_NAME FROM RNA_SEQ_READ_COUNTS WHERE EXPERIMENT_NAME IN ('%s')",
      paste(input$experiment_selector, collapse = "','")
    )
    
    sample_names <- dbGetQuery(myconn, query)
    
    updateSelectizeInput(
      session,
      inputId = "sample_selector",
      choices = as.vector(sample_names$SAMPLE_NAME),
      selected = NULL
    )
  }, ignoreNULL = FALSE)
  
  
  
  
  observeEvent(input$select_all_samples, {
    req(input$experiment_selector)
    query <- sprintf(
      "SELECT SAMPLE_NAME FROM RNA_SEQ_READ_COUNTS WHERE EXPERIMENT_NAME IN ('%s')",
      paste(input$experiment_selector, collapse = "','")
    )
    sample_names <- dbGetQuery(myconn, query)
    
    updateSelectizeInput(
      session,
      inputId = "sample_selector",
      selected = if (input$select_all_samples) sample_names$SAMPLE_NAME else character(0)
    )
  })
  
  
  output$custom_name_ui <- renderUI({
    req(input$sample_selector)
    tagList(
      helpText("Please ensure that the custom name is the same for replicates."),
      DT::dataTableOutput("custom_name_table")
    )
  })
  
  custom_names <- reactiveValues(data = data.frame(Sample = character(), CustomName = character(), stringsAsFactors = FALSE))
  
  observeEvent(input$sample_selector, {
    req(input$sample_selector)
    
    new_samples <- setdiff(input$sample_selector, custom_names$data$Sample)
    if (length(new_samples) > 0) {
      new_rows <- data.frame(Sample = new_samples, CustomName = rep("", length(new_samples)), stringsAsFactors = FALSE)
      custom_names$data <- rbind(custom_names$data, new_rows)
    }
    
    custom_names$data <- custom_names$data[custom_names$data$Sample %in% input$sample_selector, ]
  })
  
  observeEvent(input$custom_name_table_cell_edit, {
    info <- input$custom_name_table_cell_edit
    if (info$col == 1) {  
      custom_names$data[info$row, "CustomName"] <- info$value
    }
  })
  
  output$custom_name_table <- DT::renderDataTable({
    req(input$sample_selector)
    
    DT::datatable(
      custom_names$data,
      editable = list(target = "cell", disable = list(columns = 0)),
      rownames = FALSE,
      options = list(pageLength = 100)
    )
  })
  

  
  output$submit <- renderUI({
    req(length(get_grouping_from_custom_names()) > 0, !any(is.na(get_grouping_from_custom_names())))  
    actionButton("submit", "Submit")
  })
  
  
  
  custom_name_index <- reactive({
    req(custom_names$data)
    unique(data.frame(CustomName = custom_names$data$CustomName, Group = get_grouping_from_custom_names()))
  })
  
  observeEvent(input$submit, {
    sample_counts <- table(custom_names$data$CustomName)
    unique_samples <- length(sample_counts)
    min_replicates <- min(sample_counts)
    
    if (unique_samples < 2 || min_replicates < 2) {
      showModal(modalDialog(
        title = "Insufficient Replicates",
        "Duplicates are needed for DE calculation. Please update sample selection.",
        easyClose = TRUE
      ))
    }
  })
  
  observeEvent(input$submit, {
    sample_counts <- table(custom_names$data$CustomName)
    unique_samples <- length(sample_counts)
    min_replicates <- min(sample_counts)
    
    if (unique_samples < 2 || min_replicates < 2) {
      NULL
    }else{
      updateTabItems(session, "sidebar_menu", selected = "tp")  
      }
  })
  
 
  
  get_grouping_from_custom_names <- reactive({
    req(custom_names$data)
    unique_names <- unique(custom_names$data$CustomName[custom_names$data$CustomName != ""])
    group_mapping <- setNames(seq_along(unique_names), unique_names)
    group <- factor(as.numeric(group_mapping[custom_names$data$CustomName]))
    return(group)
  })
  
  observeEvent(input$submit, {
    output$reference_selector <- renderUI({
      selectInput("tp_ref", label = "Select Reference:", choices = unique(custom_names$data$CustomName))
      
    })
    output$calculate_de <- renderUI({
    actionButton("calculate_de", "Calculate")
    })
    
    output$samples_to_use <- renderUI({
      selectInput("samples_to_use", label = "Select Contrasts:", choices = setdiff(unique(custom_names$data$CustomName),input$tp_ref), multiple = TRUE)
    })
    
    
  })
  
 
  
  ## generate dge list
  y<- reactive({
    convert_to_dge_list(df = fetch_sample_data(input$sample_selector), group = get_grouping_from_custom_names() )
  })
 
  
  contrasts_to_keep<- reactive({
    req(input$samples_to_use)
    paste(input$samples_to_use, "vs", input$tp_ref)
  })
 
  de_results <- reactiveVal(NULL)
    
    observeEvent(input$calculate_de,{

    sample_counts <- table(custom_names$data$CustomName)
    unique_samples <- length(sample_counts)
    
    if (unique_samples < 2 || any(sample_counts < 2)) {
      return(NULL)  # Stop execution if validation fails
    }
    
    
    result<- within_strain_de_calculator(y(), input$tp_ref, custom_name_index(), get_grouping_from_custom_names())
    result<- result %>% filter(cond %in% contrasts_to_keep())
    de_results(result)
  })
  
  
  
  observeEvent(input$calculate_de, {
    
    output$tp_volcano_plot<- renderPlotly({
      if(is.null(de_results())){
        tpvolcano<- ""
      }
      else{
        df <- de_results()
        tpvolcano <- generate_volcano_tp(df, input$tp_volcano_logFC[1], input$tp_volcano_logFC[2], input$annotation_type)
      }
      tpvolcano
      }
      )
    
    
    tp_de_table <- reactive({
      df <- de_results()
      
      # Apply gene name mapping based on user selection
      df$name <- gene_annotation_mapper(df$name, input$annotation_type)
      
      df_filtered <- df %>% 
        select(name, cond, logFC, FDR, regulated) %>%
        filter(logFC > input$tp_logFC[1] & logFC < input$tp_logFC[2])
      
      if(input$tp_regulation == T){
        df_filtered <- df_filtered %>% filter(!regulated %in% "No change")
      }
      
      colnames(df_filtered) <- c("Gene", "Event", "logFC", "FDR", "Regulation")
      df_filtered$logFC <- round(df_filtered$logFC, 2)
      df_filtered$FDR <- sprintf("%.2e", df_filtered$FDR)
      
      return(df_filtered)
    })
    
    output$tp_comparison<- renderDataTable({
      tp_de_table()
    }) 
    ###KEGG 
    output$tp_pathway <- renderUI({
      df <- de_results()
      if(!is.null(input$event_tp_kegg)){
        requested_pathway<- pathways[pathways$name %in% input$tp_pathway,1]
        kegg_map_generator(df = df %>% filter(grepl(input$event_tp_kegg,df$cond)), requested_pathway)
        image_path<- paste0("ppa",requested_pathway,".pathview.png")
        file.remove(paste0("ppa",requested_pathway, ".png"))
        if((file.size(image_path) !=0) & (file.exists(image_path))){
          renderImage(list(src = image_path, width = "900px",  height = "800px"), deleteFile = T)
        }else{
          file.remove(paste0("ppa",requested_pathway, ".pathview.png"))  
          "Something went wrong. Try another pathway."
        }
      }
      else{ "Wait..."}
    })
    
    conversation_history <- reactiveVal(list())
    gene_context <- reactiveVal("")
    
    strain_gpt_response <- eventReactive(input$strain_chatgpt, {
      if(nrow(tp_de_table()) <= 50) {
        de_genes <- tp_de_table()
        
        withProgress(message = 'Analyzing genes...', {
          incProgress(0.1, detail = "Preparing gene data...")
          
          gene_info <- list()
          for(i in 1:nrow(de_genes)) {
            gene_id <- de_genes$Gene[i]
            
            gene_name <- genenames$name[genenames$genes == gene_id]
            if(length(gene_name) == 0) gene_name <- "Unknown"
            
            gene_pathways <- NA
            if(gene_id %in% genenameskegg$genes) {
              pathway_ids <- genenameskegg$pathway[genenameskegg$genes == gene_id]
              pathway_names <- pathways$name[pathways$id %in% pathway_ids]
              gene_pathways <- paste(pathway_names, collapse = ", ")
            }
            
            fold_change <- NA
            if("log2FoldChange" %in% colnames(de_genes)) {
              fold_change <- de_genes$log2FoldChange[i]
            }
            
            gene_info[[i]] <- paste0(
              "Gene: ", gene_id, 
              " (", gene_name, ")", 
              "\nRegulation: ", de_genes$Regulation[i],
              if(!is.na(fold_change)) paste0("\nLog2 Fold Change: ", round(fold_change, 2)) else "",
              if(!is.na(gene_pathways)) paste0("\nPathways: ", gene_pathways) else ""
            )
          }
          
          genes_context_text <- paste(unlist(gene_info), collapse = "\n\n")
          gene_context(genes_context_text)
          
          incProgress(0.2, detail = "Formulating question...")
          
          system_msg <- list(
            role = "system", 
            content = "You are a specialized bioinformatics assistant with expertise in Pichia pastoris (Komagataella phaffii) biology. Analyze gene expression data to identify functional implications, affected pathways, and biological significance. Focus on metabolic impact, stress responses, and cellular adaptations. Provide well-structured, scientifically accurate interpretations."
          )
          
          user_msg <- list(
            role = "user",
            content = paste(
              "I am analyzing differential gene expression between two strains of Pichia pastoris (Komagataella phaffii).",
              "Based on the differentially expressed genes below, what are the key functional implications, affected pathways, and potential biological significance?",
              "Please focus on metabolic impact, stress responses, and any strain-specific adaptations that might be occurring.\n\n",
              genes_context_text
            )
          )
          
          history <- list(system_msg, user_msg)
          conversation_history(history)
          
          incProgress(0.3, detail = "Sending to AI...")
          
          response <- call_openai_api_with_history(history, model = "gpt-4o")
          
          assistant_msg <- list(
            role = "assistant",
            content = response
          )
          
          conversation_history(c(history, list(assistant_msg)))
          
          incProgress(1.0, detail = "Analysis complete")
          return(response)
        })
      } else {
        return("Please select no more than 50 genes for AI analysis")
      }
    })
    
    followup_response <- eventReactive(input$submit_followup, {
      if(input$followup_question == "") {
        return("Please enter a follow-up question")
      }
      
      withProgress(message = 'Processing follow-up...', {
        history <- conversation_history()
        
        new_question <- list(
          role = "user",
          content = input$followup_question
        )
        
        updated_history <- c(history, list(new_question))
        
        response <- call_openai_api_with_history(updated_history, model = "gpt-4o")
        
        new_response <- list(
          role = "assistant",
          content = response
        )
        
        conversation_history(c(updated_history, list(new_response)))
        
        return(response)
      })
    })
    
    output$gpt_response_container <- renderUI({
      if(is.null(strain_gpt_response()) || strain_gpt_response() == "") {
        wellPanel(
          htmlOutput("strain_gpt_response")
        )
      } else {
        wellPanel(
          htmlOutput("strain_gpt_response"),
          hr(),
          textAreaInput("followup_question", "Ask a follow-up question about these genes:", ""),
          actionButton("submit_followup", "Ask Follow-up"),
          htmlOutput("followup_response"),
          hr(),
          downloadButton("download_conversation", "Download Conversation")
        )
      }
    })
    
    output$strain_gpt_response <- renderUI({
      HTML(markdown::markdownToHTML(text = strain_gpt_response(), fragment.only = TRUE))
    })
    
    output$followup_response <- renderUI({
      HTML(markdown::markdownToHTML(text = followup_response(), fragment.only = TRUE))
    })
    
    #DE TABLE DOWNLOAD LOCALLY
    output$download_tp_de_table <- downloadHandler(
      filename = function() { "DE_comparison.txt" },
      content = function(file) {
        write.table(as.data.frame(tp_de_table()),  file, quote= F, col.names = T, row.names = F, sep = "\t")
      })
    
    #DE TABLE DOWNLOAD ON ICA  
    observeEvent(input$download_tp_de_table_ica, {
      dir_tables <-paste('/Users/tbors', de_table_folder_name, sep="/")
      dir.create(dir_tables, showWarnings = FALSE)
      de_table_path <- paste(dir_tables, input$save_file, sep = "/")
      if (file.exists(de_table_path)){
        showModal(modalDialog("File already exists! Use another filename!", footer=NULL)) 
        Sys.sleep(2)
        removeModal()
      }else{
        showModal(modalDialog(paste("File will be saved to ICA. Check folder:",basename(dir_tables) ), footer=NULL)) 
        write.table(as.data.frame(tp_de_table()),  de_table_path, quote= F, col.names = T, row.names = F, sep = "\t")
        Sys.sleep(2)
        removeModal()
      }
    })
    
    # Add this to your server code
    output$download_conversation <- downloadHandler(
      filename = function() {
        paste("gene-analysis-conversation-", Sys.Date(), ".txt", sep="")
      },
      content = function(file) {
        history <- conversation_history()
        
        # Format the conversation for better readability
        conversation_text <- ""
        
        for (msg in history) {
          if (msg$role == "system") {
            # Skip system messages in the download
            next
          } else if (msg$role == "user") {
            conversation_text <- paste0(conversation_text, "Question: ", msg$content, "\n\n")
          } else if (msg$role == "assistant") {
            conversation_text <- paste0(conversation_text, "Answer: ", msg$content, "\n\n")
          }
        }
        
        # Add metadata
        header <- paste0(
          "Gene Expression Analysis Conversation\n",
          "Date: ", Sys.Date(), "\n",
          "Number of genes analyzed: ", nrow(tp_de_table()), "\n\n",
          "----------------\n\n"
        )
        
        full_text <- paste0(header, conversation_text)
        
        # Write to file
        writeLines(full_text, file)
      }
    )
    updateSelectInput(session, "gsea_contrast", 
                      choices = unique(de_results()$cond))
    
    gsea_results <- reactiveVal(NULL)
    
    observeEvent(input$run_gsea, {
      req(de_results(), input$gsea_contrast)
      
      withProgress(message = 'Running GSEA analysis...', {
        incProgress(0.2, detail = "Preparing gene rankings...")
        
        filtered_de <- de_results() %>% 
          filter(cond == input$gsea_contrast)
        
        ranked_genes <- prepare_gene_rankings(filtered_de)
        
        incProgress(0.5, detail = "Loading pathways...")
        
        pathway_database <- prepare_pichia_pathways()
        
        incProgress(0.7, detail = "Running enrichment analysis...")
        
        gsea_result <- perform_gsea(ranked_genes, pathway_database)
        
        gsea_results(gsea_result)
        
        incProgress(1.0, detail = "Analysis complete")
      })
    })
    
    output$gsea_results_table <- DT::renderDataTable({
      req(gsea_results())
      
      results <- gsea_results()
      
      results$pval <- format(results$pval, digits = 3, scientific = TRUE)
      results$padj <- format(results$padj, digits = 3, scientific = TRUE)
      
      display_data <- results[, c("pathway", "pval", "padj", "NES", "size")]
      colnames(display_data) <- c("Pathway", "P-value", "Adjusted P-value", "NES", "Size")
      
      DT::datatable(display_data, 
                    options = list(pageLength = 10, 
                                   autoWidth = TRUE,
                                   scrollX = TRUE),
                    rownames = FALSE) %>%
        DT::formatStyle('NES',
                        background = DT::styleColorBar(c(-max(abs(results$NES)), max(abs(results$NES))), 
                                                       'lightblue'),
                        backgroundSize = '98% 88%',
                        backgroundRepeat = 'no-repeat',
                        backgroundPosition = 'center')
    })
    
    output$gsea_enrichment_plot <- renderPlot({
      req(gsea_results())
      
      results <- gsea_results()
      
      top_pos <- head(results[results$NES > 0, ], 10)
      top_neg <- head(results[results$NES < 0, ], 10)
      top_pathways <- rbind(top_pos, top_neg)
      
      if(nrow(top_pathways) == 0) {
        plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(0, 0, "No significant pathways found", cex = 1.5)
        return()
      }
      
      ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = padj < 0.05)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("grey", "red"), name = "Significant") +
        coord_flip() +
        labs(title = "Top Enriched Pathways", x = "Pathway", y = "Normalized Enrichment Score") +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 10))
    })
    
    output$gsea_pathway_plot <- renderPlot({
      req(gsea_results(), input$gsea_results_table_rows_selected)
      
      results <- gsea_results()
      selected_row <- input$gsea_results_table_rows_selected
      
      if(length(selected_row) == 0) {
        plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(0, 0, "Select a pathway from the table to see enrichment plot", cex = 1.2)
        return()
      }
      
      selected_pathway <- results$pathway[selected_row]
      
      filtered_de <- de_results() %>% 
        filter(cond == input$gsea_contrast)
      
      ranked_genes <- prepare_gene_rankings(filtered_de)
      
      pathway_database <- prepare_pichia_pathways()
      pathway_genes <- pathway_database[[selected_pathway]]
      
      plotEnrichment(pathway_genes, ranked_genes) + 
        labs(title = selected_pathway) +
        theme_minimal()
    })
    
    output$download_gsea_results <- downloadHandler(
      filename = function() { 
        paste0("GSEA_results_", gsub(" ", "_", input$gsea_contrast), "_", Sys.Date(), ".csv") 
      },
      content = function(file) {
        write.csv(gsea_results(), file, row.names = FALSE)
      }
    )
    
    output$download_volcano_tp <- downloadHandler(
      filename = function() {"volcano_plot_tp.png"},
      content = function(file) {
        df <- de_results()
        if (is.null(df)) return(NULL)
        
        plot_to_save <- generate_volcano_ggplot(df, input$tp_volcano_logFC[1], input$tp_volcano_logFC[2], input$annotation_type)
        
        ggsave(file, plot = plot_to_save, width = 290, height = 265, units = "mm", device = "png")
      }
    )
    

  })
}

# Run the app
shinyApp(ui, server)