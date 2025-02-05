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




#function to use chatGPT API
call_openai_api <- function(question) {
  url <- "https://api.openai.com/v1/chat/completions"
  headers <- add_headers(
    "Authorization" = paste("Bearer", api_key_chatgpt),
    "Content-Type" = "application/json")
  body <- list(
    model = "gpt-3.5-turbo",
    messages = list(
      list(role = "user", content = question)),
    temperature = 0.7)
  response <- POST(url, headers, body = toJSON(body, auto_unbox = TRUE))
  content <- content(response, "parsed")
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
  
  dbGetQuery(myconn, query)
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
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Select samples", tabName = "upload"),
    menuItem("Group comparions", tabName = "tp"),
    menuItem("Annotations index", tabName = "full_index")
    
  )
)

jscode <- "shinyjs.refresh_page = function() { history.go(0); }"
ui <- dashboardPage(
  dashboardHeader(title = "Differential Expression"),
  dashboardSidebar(
    sidebarMenu(
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
        uiOutput("drop_down_tp_volcano"),
        br(), br(),
        uiOutput("calculate_de_ui") 
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
            uiOutput("volcano_samples_selector"),
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
                  actionBttn("strain_chatgpt", "Ask ChatGPT about this result")
                )
              ),
              wellPanel(
                textOutput("strain_gpt_response")
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
    if (info$col == 1) {  # Ensure updates only happen in the "CustomName" column
      custom_names$data[info$row, "CustomName"] <- info$value
    }
  })
  
  output$custom_name_table <- DT::renderDataTable({
    req(input$sample_selector)
    
    DT::datatable(
      custom_names$data,
      editable = list(target = "cell", disable = list(columns = 0)),
      rownames = FALSE,
      options = list(pageLength = 100)  # Set default rows per page to 100
    )
  })
  
  
  output$drop_down_tp_volcano <- renderUI({
    req(length(get_grouping_from_custom_names()) > 0, !any(is.na(get_grouping_from_custom_names())))
    
    selectInput("tp_ref",
                label = "Select reference:",
                choices = NULL,  
                selected = NULL
                )
  })
  
  output$calculate_de_ui <- renderUI({
    req(input$tp_ref)
    actionButton("calculate_de", "Calculate")
  })
  
  
  
  custom_name_index <- reactive({
    req(custom_names$data)
    unique(data.frame(CustomName = custom_names$data$CustomName, Group = get_grouping_from_custom_names()))
  })
  
  observeEvent(input$calculate_de, {
    sample_counts <- table(custom_names$data$CustomName)
    print(sample_counts)  # Debugging step to verify grouping
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
  
 
  
  observeEvent(custom_name_index(), {
    updateSelectInput(session, "tp_ref",
                      choices = unique(custom_names$data$CustomName),
                      selected = character(0))  # 🔹 Prevents automatic selection
  })
  
  get_grouping_from_custom_names <- reactive({
    req(custom_names$data)
    unique_names <- unique(custom_names$data$CustomName[custom_names$data$CustomName != ""])
    group_mapping <- setNames(seq_along(unique_names), unique_names)
    group <- factor(as.numeric(group_mapping[custom_names$data$CustomName]))
    return(group)
  })
  
  observeEvent(input$calculate_de,{
    output$reference_selector<- renderUI({
     # selectInput("volcano_reference", choices = c("cat","dog"), label = "Select Reference:" )
    })
  })
  
  ## generate dge list
  y<- reactive({
    convert_to_dge_list(df = fetch_sample_data(input$sample_selector), group = get_grouping_from_custom_names() )
  })
  
  
  
  de_results <- eventReactive(input$calculate_de, {
    req(custom_names$data)
    
    sample_counts <- table(custom_names$data$CustomName)
    unique_samples <- length(sample_counts)
    
    if (unique_samples < 2 || any(sample_counts < 2)) {
      return(NULL)  # Stop execution if validation fails
    }
    

      within_strain_de_calculator(y(), input$tp_ref, custom_name_index(), get_grouping_from_custom_names())
   
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
    
    
    #CHATGPT FUNCTION
    strain_gpt_response<- eventReactive(input$strain_chatgpt,{
     if(nrow(tp_de_table())<51){
        list_of_genes <- tp_de_table() %>% select(Gene, Regulation)
        list_of_genes<- paste(apply(list_of_genes, 1, paste, collapse = "\t"), collapse = "\n")
        question <- paste("I am running a differential transcriptomics expreriment between two strains of Pichia pastoris. What are the functional implications of this result?", list_of_genes)
        response <- as.vector(call_openai_api(question))
        return(response)
      }
      else{
        return("Please select no more than 50 genes for ChatGPT requests")
      }
    })
    output$strain_gpt_response<- renderText({
      strain_gpt_response()
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
