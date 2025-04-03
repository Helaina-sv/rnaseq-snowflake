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
library(ggplot2)
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

api_key_chatgpt <- Sys.getenv("CHAT_GPT_KEY")


###extra files
genenameskegg<- read.delim("kegg.txt")
genenames<- read.delim("gene_name_list.txt", header = F)
colnames(genenames)<- c("genes","name")
full_index<- read_excel("Full_GS115_index.xlsx")
colnames(full_index)[3]<- "genes"
genenames$genes<- paste0("lcl|",genenames$genes)

pathways<- read.delim("pichia_pathways.txt")
pathways$id<- gsub("ppa","", pathways$id)




call_openai_api <- function(question, model = "gpt-4o") {
  url <- "https://api.openai.com/v1/chat/completions"
  headers <- add_headers(
    "Authorization" = paste("Bearer", api_key_chatgpt),
    "Content-Type" = "application/json")
  
  system_msg <- "You are a specialized bioinformatics assistant with expertise in Pichia pastoris (Komagataella phaffii) biology. Analyze gene expression data to identify functional implications, affected pathways, and biological significance. Focus on metabolic impact, stress responses, and cellular adaptations. Provide well-structured, scientifically accurate interpretations."
  
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

dotenv::load_dot_env()
connect_to_snowflake <- function(warehouse, database, schema) {
  # Attempt to connect to Snowflake using provided credentials
  success <- tryCatch({
    myconn <- dbConnect(odbc::odbc(), Sys.getenv("SNOWFLAKE_USER"), role = 'ACCOUNTADMIN', PWD = Sys.getenv("SNOWFLAKE_PASSWORD"))
    
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


generate_volcano_tp <- function(df, logfc1, logfc2){
  df$cond<- factor(df$cond, levels = unique(df$cond))
  tpvolcano<- ggplot(df %>%
                       filter(abs(logFC)<= logfc2)%>%
                       filter(abs(logFC)>= logfc1),
                     aes(logFC, -log10(FDR),
                         label = name,
                         color = 
                           ifelse(regulated %in% "Up-regulated", "#7cae00",
                                  ifelse(regulated %in% "Down-regulated","#f8766d", "#00bfc4"
                                  )
                           ) ) ) +
    geom_point( size = 0.5) +
    scale_color_identity(guide = "none") +
    facet_wrap(~cond, scales = "free")+
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "gray40", linewidth = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40", linewidth = 0.5) +
    labs(x = expression(Log[2]~"Fold Change"), y = expression("-log"[10]~"FDR")) +
    theme_bw()+
    theme(text  = element_text( face = "bold", size = 12))
  
  tpvolcano
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
            width = 6,
            uiOutput("reference_selector"),
            
          ),
          column(
            width = 6,
            uiOutput("samples_to_use"),
          )
        ),
        fluidRow(
          column(
            width = 6,
            uiOutput("calculate_de"),
          )
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
                column(width = 4, checkboxInput("tpshowlabels", "Label genes?", value = FALSE)),
                column(width = 4, sliderInput(inputId = "tp_volcano_logFC", label = "Adjust Fold change", min = 0, max = 20, value = c(0, 10)))
              ),
              plotOutput(
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
    
    output$tp_volcano_plot<- renderPlot({
      if(is.null(de_results())){
        tpvolcano<- ""
      }
      else{
        df <- de_results()
        tpvolcano <- generate_volcano_tp(df, input$tp_volcano_logFC[1], input$tp_volcano_logFC[2])
        if (input$tpshowlabels == TRUE) {
          tpvolcano<- tpvolcano +
            geom_label_repel(aes(label = ifelse(abs(logFC)>1 & FDR < 0.05, name, NA)),
                             max.overlaps = 30)}
      }
      tpvolcano})
    
    
    #de table tab
    tp_de_table<- reactive({
      df<- de_results()
      df_filtered <- df %>% 
        select(name,cond,logFC,FDR, regulated)%>%
        filter(logFC > input$tp_logFC[1] & logFC < input$tp_logFC[2])
      if(input$tp_regulation == T){
        df_filtered <- df_filtered %>% filter(!regulated %in% "No change")}
      colnames(df_filtered)<- c("Gene","Event","logFC","FDR","Regulation")
      df_filtered$logFC <- round(df_filtered$logFC,2)
      df_filtered$FDR <- sprintf("%.2e",df_filtered$FDR)
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
    
    output$download_volcano_tp <- downloadHandler(
      filename = function() {"volcano_plot_tp.png"},content = function(file) {
        df <- de_results()
        if(is.null(df)){
          tpvolcano<- ""}
        else {
          tpvolcano <- generate_volcano_tp(df, input$tp_volcano_logFC[1], input$tp_volcano_logFC[2])
          if (input$tpshowlabels == TRUE) {
            tpvolcano<- tpvolcano+
              geom_label_repel(aes(label = ifelse(abs(logFC)>1 & FDR < 0.05, name, NA)),
                               max.overlaps = 30)
          }
        }
        ggsave(file, plot = tpvolcano, width = 290, height = 265, units = "mm", device = "png")
      }) 
  })
}

# Run the app
shinyApp(ui, server)