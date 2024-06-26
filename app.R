
# make sure these packages are installled
library("shiny")
library("shinyjs")
library("tidyverse")
library("mixtools")
library("dplyr")
library("ggplot2")
library("janitor")
library("refineR")
library("shinycssloaders")
library("viridisLite")
library("DT")
library("bestNormalize")
library("shinythemes")

# PLEASE 
addResourcePath("myfiles", "~/Documents/GitHub/lab_ref_tool/reference_range_GUItool/data_template")

`%notin%` <- negate(`%in%`)

Remove.Outliers <- function(x, na.rm = TRUE) 
{
  ## Find 25% and 75% Quantiles using inbuild function
  quant <- quantile(x, probs=c(.25, .75), na.rm = na.rm)
  
  ## Find Interquantile range and multiply it by 1.5 
  ## to derive factor for range calculation
  H <- 1.5 * IQR(x, na.rm = na.rm)
  
  y <- x
  
  ## fill the outlier elements with NA
  y[x < (quant[1] - H)] <- NA
  y[x > (quant[2] + H)] <- NA
  
  return(as.data.frame(y))
}


Hoff.QQ.plot <- function(x, alpha=0.05, from=NA, to=NA, xlim=range(x)) {
  
  x <- sort(x)
  qq.data <- as.data.frame(qqnorm(x, datax = TRUE, plot.it = FALSE))
  plot(y ~ x, 
       data = qq.data,
       type = "l",
       xlab = "Analyte Values",
       ylab = "Quantiles of the Normal Distibution",
       xlim = xlim)
  if (!is.na(from) & !is.na(to)) {
    linear <- subset(qq.data, x >= from & x <= to)
    lin.mod <- lm(y ~ x, data = linear)
    abline(lin.mod)
    RI <- (c(qnorm(alpha/2),qnorm(1-alpha/2))
           - lin.mod$coefficients[1])/lin.mod$coefficients[2]
    result <- c(RI[1], RI[2])
    names(result) <- c("lower", "upper")
    abline(h = c(qnorm(alpha/2), qnorm(1 - alpha/2)), col = "blue")
    abline(v = RI, col = "blue", lty = 2)
  }
  
}


Get.Hoff.Results <- function(x, alpha=0.05, 
                             from=NA, to=NA, xlim=range(x)) {
  x <- sort(x)
  
  qq.data <- as.data.frame(qqnorm(x, datax = TRUE, plot.it = FALSE))
  
  if (!is.na(from) & !is.na(to)) {
    linear <- subset(qq.data, x >= from & x <= to)
    lin.mod <- lm(y ~ x, data = linear)
    RI <- (c(qnorm(alpha/2),qnorm(1-alpha/2)) - 
             lin.mod$coefficients[1])/lin.mod$coefficients[2]
    result <- data.frame("lower limit" = RI[1],"upper limit"=RI[2])
  }
  return(result)
}

get_transform_func = function(x,
                              transform_type){
  if(transform_type == "yj"){
    transform_func = bestNormalize::yeojohnson(x,standardize = FALSE)
  }   
  if(transform_type == "none"){
    transform_func = bestNormalize::no_transform(x,standardize = FALSE)
  } 
  return(transform_func)
}

transform_data = function(x,
                          transform_function = NULL,
                          inverse = FALSE){
  return(
    predict(transform_function,newdata=x,inverse)
  )
}


# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("cerulean"),
  fluidRow(
    column(4,
           # Application title
           titlePanel("Reference Interval GUI Calculator")),
    column(4, fileInput(inputId = "data", 
                        label = "Import Data (.csv)", 
                        multiple = FALSE, 
                        accept = c(
                          "text/csv",
                          "text/comma-separated-values,text/plain",
                          ".csv"), 
                        buttonLabel = "Import File")),
    column(4, h3("What should my data look like?"),a("Download the data template",  href = "myfiles/data template.csv")
    )
  ),
  
  fluidRow(
    tabsetPanel(
      tabPanel("Data Import and View", fluid = TRUE,
               sidebarLayout(
                 sidebarPanel(
                   selectInput(inputId = "analyteSelect", 
                               label = "Select Analyte", 
                               choices = NULL), # Choices will be updated based on uploaded data
                   helpText("This is graphical user interface for calculating reference range intervals using a suite of methods. Import your data above and go through the tabs to start computing intervals. Check the data template in the top-right to make sure your data is in the right format.")
                 ),
                 mainPanel(
                   shinycssloaders::withSpinner(plotOutput(outputId = "analyteHistogram")),
                   shinycssloaders::withSpinner(DTOutput("analyteSummary"))
                 )
               )
      ),
      tabPanel("Mixture Modeling", fluid = TRUE,
               sidebarLayout(
                 sidebarPanel(selectInput(inputId = "Analyte_mix", 
                                          "Select Analyte", choices = NULL, selected = ""),
                              numericInput(inputId = "iterations", 
                                           label = "# of iterations", value = 10000),
                              numericInput(inputId = "population", 
                                           label = "Expected # of Populations", 
                                           min=2, max=4, value = 2),
                              numericInput(inputId = "age_lb_mix", 
                                           label = "age lowerbound ", value = 0),
                              numericInput(inputId = "age_ub_mix", 
                                           label = "age upperbound ", min=0, 
                                           max=150, value=100),
                              radioButtons("sex_mix", h3("Sex Selection"),
                                           choices = list("All" = 1, "Male" = 2,
                                                          "Female" = 3),selected = 1),
                              # Added radio buttons for transformation options
                              radioButtons(inputId = "transformation_mix", 
                                           label = "Pretransform Data",
                                           choices = list("None" = "none", 
                                                          "Yeo-johnson" = "yj"),
                                           selected = "none"),
                              actionButton(inputId = "go_mix", label = "Run"),
                              downloadButton("downloadDataMix", "Download"),
                              shiny::helpText("This tab uses Gaussian mixture modeling to compute reference range intervals implemented via the mixtools package in R. It assumes the data is generated by sampling different Gaussian-distributed population (e.g. a non-pathological population and pathological population).  Each population can have their own unique mean, variance, and frequency (lambda).  The lambda parameter refers to the estimated proportion of the population from that normal distribution.")
                 ),
                 mainPanel(
                   fluidRow(
                     shinycssloaders::withSpinner(plotOutput(outputId = "mixPlot")),
                     tableOutput(outputId = "ref.int.mix.table"),
                     shinycssloaders::withSpinner(plotOutput(outputId = "mixPlot_os"))
                   )
                 )
               )
      ),
      tabPanel("Hoffmann", fluid = TRUE,
               sidebarLayout(
                 sidebarPanel(selectInput(inputId="Analyte_hoff", "Select Analyte", choices = NULL, selected = ""),
                              
                              numericInput(inputId = "age_lb_hoff", label = "Age lowerbound ", value = 0),
                              numericInput(inputId = "age_ub_hoff", label = "Age upperbound ", min=0, max=150, value = 100),
                              radioButtons("sex_hoff", h3("Sex Selection"),
                                           choices = list("All" = 1, "Male" = 2,
                                                          "Female" = 3),selected = 1),
                              # i don't think transformations really make sen
                              # Added radio buttons for transformation options
                              # radioButtons(inputId = "transformation_hoff", label = "Transformation",
                              #              choices = list("None" = "none", "Yeo-johnson" = "yj"), selected = "none"),
                              numericInput(inputId = "linearLL", "Please select the range x-axis values between the black lines that is most linear (Best Estimate) \
                                           Lowerbound", value = -1e6),     
                              numericInput(inputId = "linearUL", "Upperbound",value= 1e6),
                              actionButton(inputId = "go_hoff", label = "Run"),
                              downloadButton("downloadDataHoffman", "Download"),
                              shiny::helpText("This tab uses the graphical Hoffman method to estimate refererence intervals. It requires user assement of the plot to generate reference range intervals",a("Read here for more information",href="https://academic.oup.com/ajcp/article/151/3/328/5208945"),"."),
                 ),
                 mainPanel(
                   fluidRow(
                     shinycssloaders::withSpinner(plotOutput(outputId = "hoff_plot")),
                     tableOutput(outputId = "hoff_table")
                   )
                 )
               )
      ),
      tabPanel("Nonparametric Quantiles", fluid = TRUE,
               sidebarLayout(
                 sidebarPanel(selectInput(inputId="Analyte_nonpara", "Select Analyte", choices = "", selected = ""),
                              
                              numericInput(inputId = "age_lb_nonpara", label = "Age lowerbound ", value = 0),
                              numericInput(inputId = "age_ub_nonpara", label = "Age upperbound ", min=0, max=150, value = 100),
                              radioButtons("sex_nonpara", h3("Sex Selection"),
                                           choices = list("All" = 1, "Male" = 2,
                                                          "Female" = 3),selected = 1),                              # Added radio buttons for transformation options
                              # radioButtons(inputId = "transformation_nonpara", label = "Transformation",
                              #              choices = list("None" = "none", "Yeo-johnson" = "yj"), selected = "none"),
                              actionButton(inputId = "go_nonpara", label = "Run"),
                              downloadButton("downloadDataNonpara", "Download"),
                              downloadButton("downloadDataNonparaTukey", "Download (outlier's removed)"),
                              helpText("For this tab, reference range is computed using the observed 97.5th and 2.5th quantiles in the data. This can be done with optional Tukey outlier removal. This quantile-based approach could be useful for data where direct approaches calculating reference intervals are applicable.")
                 ),
                 mainPanel(
                   fluidRow(
                     column(6,
                            shinycssloaders::withSpinner(plotOutput(outputId = "nonpara_plot")),
                            tableOutput(outputId = "nonpara_table")
                     ),
                     column(6,
                            shinycssloaders::withSpinner(plotOutput(outputId = "nonpara_plot_tukey")),
                            tableOutput(outputId = "nonpara_tukey_table")
                     )
                   )
                 )
               )
      ),
      tabPanel("RefineR", fluid = TRUE,
               sidebarLayout(
                 sidebarPanel(selectInput(inputId="Analyte_refineR", "Select Analyte", choices = "", selected = ""),
                              
                              numericInput(inputId = "age_lb_refineR", 
                                           label = "Age lowerbound ", value = 0),
                              numericInput(inputId = "age_ub_refineR", 
                                           label = "Age upperbound ", min=0, 
                                           max=150, value = 100),
                              numericInput(inputId = "RefineR_N_bootstrap", 
                                           label = "N bootstrap samples", min=0, 
                                           value = 0, step = 1),
                              selectInput("refineR_method", "Transformation", 
                                          choices = c("BoxCox", "modBoxCoxFast", 
                                                      "modBoxCox")),
                              radioButtons("sex_refiner", h3("Sex Selection"),
                                           choices = list("All" = 1, "Male" = 2,
                                                          "Female" = 3),
                                           selected = 1),
                              radioButtons(inputId = "plot_CI", 
                                           label = "Plot Confidence Interval", 
                                           choices = c("True","False")),
                              
                              actionButton(inputId = "go_refineR", 
                                           label = "Run"),
                              downloadButton("downloadDataRefineR", "Download"),
                              shiny::helpText("RefinerR works great for creating reference 
                                       intervals from indirect real world data and 
                                       it can accommodable non-normal skewed unimodal 
                                       distributions of non-pathological data. 
                                       It assumes the majority of made up of non-pathological test results.
                                       More details",
                                       a("here",href="https://www.nature.com/articles/s41598-021-95301-2")),
                 ),
                 mainPanel(
                   fluidRow(
                     column(6,
                            shinycssloaders::withSpinner(plotOutput(outputId = "refineR_plot")),
                            tableOutput(outputId = "refineR_table")
                     )
                   )
                 )
               )
      )
    )
  )
)



server <- function(input, output, session) {
  
  data<-eventReactive(input$data, {
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    inFile<-input$data
    if (is.null(inFile))
      return(NULL)
    
    read_csv(inFile$datapath) %>%
      clean_names()
  })
  
  
  # Observe the data and update dropdown choices
  observe({
    df <- data()
    if (is.null(df)) return()
    
    # Assuming 'sex', 'age', etc. are not analytes
    analytes <- setdiff(names(df), c("sex", "age", "dob", "sample_id", "measurement_type", "notes", "collection_date"))
    updateSelectInput(session, "analyteSelect", choices = analytes)
  })
  
  # Render histogram based on selected analyte
  output$analyteHistogram <- renderPlot({
    df <- data()
    selectedAnalyte <- input$analyteSelect
    if (is.null(df) || is.null(selectedAnalyte)) return()
    
    # Filter out NA values to avoid issues with histogram
    validData <- df[[selectedAnalyte]] %>% na.omit()
    
    # Generate histogram
   ggplot(data = data.frame(validData), aes(x = validData)) +
      geom_histogram(binwidth = (max(validData) - min(validData)) / 30) +
      labs(title = paste("Histogram of", selectedAnalyte), x = selectedAnalyte, y = "Frequency")+
      theme_minimal()
  })
  
  # Render summary statistics for the selected analyte with DT
  output$analyteSummary <- renderDT({
    df <- data()
    selectedAnalyte <- input$analyteSelect
    if (is.null(df) || is.null(selectedAnalyte)) return()
    
    # Direct summary statistics
    stats <- summary(df[[selectedAnalyte]] %>% na.omit())
    
    # Convert to data frame correctly
    summaryStats <- data.frame(
      Statistic = names(stats),
      Value = as.numeric(stats),
      stringsAsFactors = FALSE # Avoid factors to ensure DT handles it smoothly
    )
    
    # Pass to DT
    datatable(summaryStats, 
              options = list(pageLength = 6, searching = FALSE), 
              rownames = FALSE) %>%
      formatStyle(columns = c("Statistic"), 
                  fontWeight = 'bold', 
                  color = 'black', 
                  backgroundColor = 'lightgrey')
  }, server = FALSE)
  
  observeEvent(data(), {
    header <- colnames(data())[!colnames(data()) %in% c("sex","age","dob",
                                                        "sample_id","n","measurement_type",
                                                        "notes","collection_date")]
    
    # Can also set the label and select items
    updateSelectInput(session, "Analyte_mix",
                      label = paste("Select Analyte"),
                      choices = header,
                      selected = tail(header, 1)
    )
    
    # Can also set the label and select items
    updateSelectInput(session, "Analyte_hoff",
                      label = paste("Select Analyte"),
                      choices = header,
                      selected = tail(header, 1)
    )
    
    # Can also set the label and select items
    updateSelectInput(session, "Analyte_nonpara",
                      label = paste("Select Analyte"),
                      choices = header,
                      selected = tail(header, 1)
    )
    
    # Can also set the label and select items
    updateSelectInput(session, "Analyte_refineR",
                      label = paste("Select Analyte"),
                      choices = header,
                      selected = tail(header, 1)
    )
  })
  
  
  #mixtool specific information 
  analyte.mix<-eventReactive(input$go_mix, {analyte<-input$Analyte_mix})
  unit_mix<-eventReactive(input$go_mix, {unit<-input$unit_mix})
  population<-eventReactive(input$go_mix, {population<-input$population})
  
  
  sex_mix <- eventReactive(input$go_mix, {
    if(input$sex_mix==1){ sex_mix <- 'All'
    }
    if(input$sex_mix==2){ sex_mix <- 'Male'
    }
    if(input$sex_mix==3){ sex_mix <- 'Female'
    }
    sex_mix
  }
  )
  
  sex_hoff<- eventReactive(input$go_hoff, {
    if(input$sex_hoff==1){ sex_hoff <- 'All'
    }
    if(input$sex_hoff==2){ sex_hoff <- 'Male'
    }
    if(input$sex_hoff==3){ sex_hoff <- 'Female'
    }
    sex_hoff
  }
  )
  
  sex_nonpara <- eventReactive(input$go_nonpara, {
    if(input$sex_nonpara==1){ sex_nonpara <- 'All'
    }
    if(input$sex_nonpara==2){ sex_nonpara <- 'Male'
    }
    if(input$sex_nonpara==3){ sex_nonpara <- 'Female'
    }
    sex_nonpara
  }
  )
  
  sex_refiner <- eventReactive(input$go_refineR, {
    if(input$sex_refiner==1){ sex_refineR <- 'All'
    }
    if(input$sex_refiner==2){ sex_refineR <- 'Male'
    }
    if(input$sex_refiner==3){ sex_refineR <- 'Female'
    }
    sex_refineR
  }
  )
  
  age_lb_mix <- eventReactive(input$go_mix, {age_lb_mix<-input$age_lb_mix})
  age_ub_mix <- eventReactive(input$go_mix, {age_ub_mix<-input$age_ub_mix})
  
  data_clean_mix <- eventReactive(input$go_mix, {
    data_clean_mix_temp = data()[[analyte.mix()]]
    
    data_clean_mix_temp = data_clean_mix_temp[data()[["age"]] >= age_lb_mix()] 
    data_clean_mix_temp = data_clean_mix_temp[data()[["age"]] <= age_ub_mix()] 
    
    if(sex_mix() == 'Male'){
      data_clean_mix_temp = data_clean_mix_temp[grep(pattern = '^m',
                                                     x = data()[[which(tolower(colnames(data())) == 'sex')]],
                                                     ignore.case = TRUE)] 
    } else if (sex_mix() == 'Female'){
      data_clean_mix_temp = data_clean_mix_temp[grep(pattern = '^f',
                                                     x = data()[[which(tolower(colnames(data())) == 'sex')]],
                                                     ignore.case = TRUE)] 
    }
    
    data_clean_mix_temp = data_clean_mix_temp |>
      na.omit()
    print(data_clean_mix_temp)
    return(data_clean_mix_temp)
  })
  
  transformation_mix = eventReactive(input$go_mix,{
    print("Fetching transformation function...")
    transformation_func <- get_transform_func(data_clean_mix(), input$transformation_mix)
    print(transformation_func)
    return(transformation_func)
  })
  
  transformed_mix_data = eventReactive(transformation_mix(),{
    print("Transforming data...")
    transformed_mix_data1 <- transform_data(x = data_clean_mix(),
                                            transform_function = transformation_mix())
    print(transformed_mix_data1)
    transformed_mix_data1
  })
  
  # transformed_mix_data = eventReactive(transformation_mix(),{
  #   transformed_mix_data1 <- transform_data(x = data_clean_mix(),
  #                                         transform_function = transformation_mix())
  # })
  
  
  
  mixtool.fit<-eventReactive(input$go_mix, {
    mixtool.fit <- normalmixEM(x =transformed_mix_data(), 
                               lambda = NULL, 
                               mu = NULL, sigma = NULL, 
                               k = input$population,
                               mean.constr = NULL, 
                               sd.constr = NULL,
                               epsilon = 1e-08, 
                               maxit = input$iterations, 
                               maxrestarts=20,
                               verb = FALSE, 
                               fast=FALSE, ECM = FALSE,
                               arbmean = TRUE, arbvar = TRUE)
  })
  
  # send plot to outputs 
  
  output$mixPlot_os <- renderPlot({
    if( (is.null(mixtool.fit())) | (input$transformation_mix!="yj"))return(NULL)
    #req(transformed_data_mix())
    cols = viridisLite::viridis(n=population())
    
    
    p1 <- ggplot(data = as.data.frame(data_clean_mix(), 
                                      aes(x=data_clean_mix()))) +
      geom_histogram(mapping = aes(x = data_clean_mix(), y= ..density..), 
                     inherit.aes=FALSE,
                     fill = "white", color = "blue") + 
      xlab(paste("Analyte Value", sep="")) + 
      ylab("Density")+
      theme_minimal()
    
    max_lambda = max(mixtool.fit()$lambda)
    
    for (i in 1:population()) {
      #i<-2
      if(mixtool.fit()$lambda[i]==max_lambda){
        
        x.ulimit <- max(c(transform_data(qnorm(0.99,mixtool.fit()$mu[i],
                                         mixtool.fit()$sigma[i]),
                                   transform_function =
                                     transformation_mix(),inverse = TRUE),
                          max(data_clean_mix())*1.05))
        
        x.llimit <- min(c(transform_data(qnorm(0.01,mixtool.fit()$mu[i],
                                               mixtool.fit()$sigma[i]),
                                         transform_function =
                                           transformation_mix(),inverse = TRUE),
                          max(data_clean_mix())*.95))
        
        #find the 2.5th 97.5th percentile from the mixed model fit
        analyte.lln <- transform_data(qnorm(0.025,mixtool.fit()$mu[i],
                                            mixtool.fit()$sigma[i]),
                                      transform_function =
                                        transformation_mix(),inverse = TRUE)
        analyte.uln <- transform_data(qnorm(0.975,mixtool.fit()$mu[i],
                                            mixtool.fit()$sigma[i]),
                                      transform_function =
                                        transformation_mix(),inverse = TRUE)
        p1  =  p1 +
          geom_vline(xintercept = analyte.lln) +
          geom_text(aes(x = analyte.lln, y = 0, label="2.5 percentile"), size=4, angle=90, vjust=-0.4, hjust=-1.6) +
          geom_vline(xintercept = analyte.uln ) +
          geom_text(aes(x = analyte.uln, y = 0, label="97.5 percentile"), size=4, angle=90, vjust=1, hjust=-1.5) +
          ggtitle( paste0(analyte.mix(),' on original scale')) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_x_continuous(limits = c(x.llimit,x.ulimit))
        
        # ggsave(here::here("results/figure", "figure1.tiff"), width=4, height=4, dpi=600)
        return(p1)
      }
    }
  })
  
  output$mixPlot <- renderPlot({
    if(is.null(mixtool.fit()))return(NULL)
    #req(transformed_data_mix())
    cols = viridisLite::viridis(n=population())
    
    x.ulimit <- max((mixtool.fit()$mu)+10*(mixtool.fit()$sigma))
    x.llimit <- min((mixtool.fit()$mu)-10*(mixtool.fit()$sigma))
    
    p1 <- ggplot(data = as.data.frame(transformed_mix_data(), aes(x=transformed_mix_data()))) +
      geom_histogram(mapping = aes(x=transformed_mix_data(), y= ..density..), inherit.aes=FALSE,
                     fill = "white", color = "blue") + 
      xlab(paste("Analyte Value", sep="")) + 
      ylab("Density")+
      theme_minimal()
    
    max_lambda = max(mixtool.fit()$lambda)
    
    for (i in 1:population()) {
      #i<-2
      if(mixtool.fit()$lambda[i]==max_lambda){
        
        #find the 2.5th 97.5th percentile from the mixed model fit
        analyte.lln <- qnorm(0.025,mixtool.fit()$mu[i], mixtool.fit()$sigma[i]) 
        analyte.uln <- qnorm(0.975,mixtool.fit()$mu[i], mixtool.fit()$sigma[i])
        
        p1  =  p1 +
          stat_function(geom = "line",
                        fun = dnorm,
                        args = list(mean = mixtool.fit()$mu[i],
                                    sd = mixtool.fit()$sigma[i]),
                        color = cols[i], size=1.5) +
          geom_vline(xintercept = analyte.lln) +
          geom_text(aes(x = analyte.lln, y = 0, label="2.5 percentile"), size=4, angle=90, vjust=-0.4, hjust=-1.6) +
          geom_vline(xintercept = analyte.uln ) +
          geom_text(aes(x = analyte.uln, y = 0, label="97.5 percentile"), size=4, angle=90, vjust=1, hjust=-1.5) +
          ggtitle(analyte.mix()) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_x_continuous(limits = c(x.llimit,x.ulimit))
        
        # ggsave(here::here("results/figure", "figure1.tiff"), width=4, height=4, dpi=600)
        return(p1)
      }
    }
  })
  
  
  ref.int.mix <- eventReactive(input$go_mix, {
    my_df_mix <- data.frame()
    
    
    for (i in 1:population()) {
      datalist <- data.frame(i,
                             transform_data(x = mixtool.fit()$mu[i],
                                            transform_function = 
                                              transformation_mix(), 
                                            inverse = TRUE),
                             transform_data(x = qnorm(0.025, 
                                                      mixtool.fit()$mu[i], 
                                                      mixtool.fit()$sigma[i]),
                                            transform_function = 
                                              transformation_mix(), 
                                            inverse = TRUE), 
                             transform_data(x = qnorm(0.975,
                                                      mixtool.fit()$mu[i], 
                                                      mixtool.fit()$sigma[i]),
                                            transform_function = 
                                              transformation_mix(), 
                                            inverse = TRUE), 
                             mixtool.fit()$lambda[i])
      my_df_mix <- bind_rows(my_df_mix, datalist)
      
    }
    
    names(my_df_mix) <- c("Reference", "Mean", 
                          "Lower Limit", "Upper Limit",
                          "Lambda") 
    
    output$downloadDataMix <- downloadHandler(
      filename = function() {
        paste("mixmodel-data",
              "_sex-",sex_mix(),
              "_ages-",
              age_lb_mix(),"-to-",age_ub_mix(),"_",
              Sys.Date(), '.csv', sep='')
      },
      content = function(con) {
        write.csv(my_df_mix , con)
      }
    )
    return(my_df_mix)
  }
  )
  
  output$ref.int.mix.table <- renderTable(
    return(ref.int.mix())
  )
  
  
  observeEvent(input$go_mix, {show("mixPlot")})
  
  ## Hoffmann Method information
  analyte.hoff <- eventReactive(input$go_hoff, {analyte<-input$Analyte_hoff})
  unit.hoff <- eventReactive(input$go_hoff, {unit<-input$unit_hoff})
  
  # upperbound inputs
  LL.hoff<-eventReactive(input$go_hoff, {LL<-input$linearLL})
  UL.hoff<-eventReactive(input$go_hoff, {UL<-input$linearUL})
  
  # age inputs
  age_lb_hoff <- eventReactive(input$go_hoff, {age_lb_hoff<-input$age_lb_hoff})
  age_ub_hoff <- eventReactive(input$go_hoff, {age_ub_hoff<-input$age_ub_hoff})
  
  
  hoff.plot <- eventReactive(input$go_hoff, {
    
    print("update")
    my_df_hoff <- data.frame()
    data_clean_hoff = data()[[analyte.hoff()]] 
    data_clean_hoff = data_clean_hoff[data()[["age"]] >= age_lb_hoff()] 
    data_clean_hoff = data_clean_hoff[data()[["age"]] <= age_ub_hoff()] 
    data_clean_hoff = data_clean_hoff |> na.omit()
    
    if(sex_hoff() == 'Male'){
      data_clean_hoff = data_clean_hoff[grep(pattern='^m',
                                             x=data()[[which(tolower(colnames(data()))=='sex')]],ignore.case=TRUE)] 
    } else if (sex_hoff() == 'Female'){
      data_clean_hoff = data_clean_hoff[grep(pattern='^f',
                                             x=data()[[which(tolower(colnames(data()))=='sex')]],ignore.case=TRUE)] 
    }
    hoff.limits <- c(LL.hoff(),UL.hoff())
    
    transformation_hoff = get_transform_func(data_clean_hoff, input$transformation_mix)
    
    transformed_hoff_data = transform_data(x = data_clean_hoff,
                                              transform_function = transformation_hoff)
  
    
    my_df_hoff <- Get.Hoff.Results(transformed_hoff_data,from=hoff.limits[1],
                                   to=hoff.limits[2])
    
    output$hoff_table<-renderTable(
      return(my_df_hoff))
    
    hoff.plot1 <- Hoff.QQ.plot(transformed_hoff_data,from=hoff.limits[1],
                               to=hoff.limits[2])
    
    output$downloadDataHoffman <- downloadHandler(
      filename = function() {
        paste('Hoffman-data',
              "_sex-",sex_hoff(),
              "_ages-",
              age_lb_hoff(),"-to-",age_ub_hoff(),"_",
              Sys.Date(), '.csv', sep='')
      },
      content = function(con) {
        write.csv(my_df_hoff, con)
      }
    )
    
    
    return(hoff.plot1)
    
  })
  
  
  # send plot to outputs 
  output$hoff_plot <- renderPlot({hoff.plot()})
  
  # NON-Parametric Area
  analyte_nonpara<-eventReactive(input$go_nonpara, {analyte_nonpara<-input$Analyte_nonpara})
  unit_nonpara<-eventReactive(input$go_nonpara, {unit_nonpara<-input$unit_nonpara})
  age_lb_nonpara <- eventReactive(input$go_nonpara, {age_lb_nonpara <- input$age_lb_nonpara})
  age_ub_nonpara <- eventReactive(input$go_nonpara, {age_ub_nonpara <- input$age_ub_nonpara})
  
  ############
  # outlier removed section
  ############
  data.tukey <- eventReactive(input$go_nonpara, {
    print(head(data()[[analyte_nonpara()]]))
    
    data.clean = data()[[analyte_nonpara()]]
    data.clean = data.clean[data()[["age"]] >= age_lb_nonpara()] 
    data.clean = data.clean[data()[["age"]] <= age_ub_nonpara()] 
    
    data.clean= data.clean |> na.omit()
    if(sex_nonpara() == 'Male'){
      data.clean= data.clean[grep(pattern='^m',
                                  x=data()[[which(tolower(colnames(data()))=='sex')]],ignore.case=TRUE)] 
    } else if (sex_nonpara() == 'Female'){
      data.clean= data.clean[grep(pattern='^f',
                                  x=data()[[which(tolower(colnames(data()))=='sex')]],ignore.case=TRUE)] 
    }
    
    tukey<-Remove.Outliers(data.clean)
    
    return(tukey)
    
  })
  
  # create reference range plot 
  nonpara.plot.tukey<-eventReactive(input$go_nonpara, {
    
    my_df_nonpara <- data.frame()
    nonpara_ref_tukey <- quantile(data.tukey(), probs = c(0.025,0.975), na.rm = TRUE)
    
    my_df_nonpara_tukey <- rbind(my_df_nonpara, nonpara_ref_tukey)
    names(my_df_nonpara_tukey) <- c("Lower Limit (2.5%)", "Upper Limit (97.5%)")
    output$nonpara_tukey_table <- renderTable(
      return(my_df_nonpara_tukey))
    
    # num_bin <- as.integer(1+2.5*log10(count(data.tukey())))
    # print(paste("this is the bin ", num_bin, sep=""))
    
    nonpara.plot.tukey <-
      ggplot(as.data.frame(data.tukey()), aes(sample=data.tukey()))+
      geom_histogram(mapping=aes(x=array(as.numeric(unlist(data.tukey()))),
                                 y= ..density..), 
                     inherit.aes=FALSE,
                     fill = "white", color = "blue") +
      xlab(paste("Analyte Values", sep="")) +
      ylab("Density")+
      ggtitle(paste(analyte_nonpara(), " with Tukey outliers removed", sep="")) + 
      theme(plot.title = element_text(hjust = 0.5))+
      geom_vline(xintercept = nonpara_ref_tukey, color = "black")+
      theme_minimal()
    
    output$downloadDataNonparaTukey <- downloadHandler(
      filename = function() {
        paste('nonpara-tukey-data',"_sex-",sex_nonpara(),'-', Sys.Date(), '.csv', sep='')
      },
      content = function(con) {
        write.csv(my_df_nonpara_tukey, con)
      }
    )
    return(nonpara.plot.tukey)
  })
  
  # send plot to outputs
  output$nonpara_plot_tukey <- renderPlot({nonpara.plot.tukey()})
  
  # without outlier removeal
  nonpara.plot <- eventReactive(input$go_nonpara, {
    my_df_nonpara <- data.frame()
    
    data_clean_nonpara = data()[[analyte_nonpara()]]
    data_clean_nonpara = data_clean_nonpara[data()[["age"]] >= age_lb_nonpara()] 
    data_clean_nonpara = data_clean_nonpara[data()[["age"]] <= age_ub_nonpara()] 
    data_clean_nonpara = data_clean_nonpara |> na.omit()
    
    data_clean_nonpara = data_clean_nonpara|> na.omit()
    if(sex_nonpara() == 'Male'){
      data_clean_nonpara= data_clean_nonpara[grep(pattern='^m',
                                                  x=data()[[which(tolower(colnames(data()))=='sex')]],ignore.case=TRUE)] 
    } else if (sex_nonpara() == 'Female'){
      data_clean_nonpara= data_clean_nonpara[grep(pattern='^f',
                                                  x=data()[[which(tolower(colnames(data()))=='sex')]],ignore.case=TRUE)] 
    }
    
    nonpara_limits <- quantile(data_clean_nonpara, 
                               probs = c(0.025,0.975), 
                               na.rm = TRUE)
    
    my_df_nonpara <- rbind(my_df_nonpara, 
                           nonpara_limits)
    
    names(my_df_nonpara) <- c("Lower Limit (2.5%)", 
                              "Upper Limit (97.5%)")
    
    output$nonpara_table<- renderTable(
      return(my_df_nonpara))
    
    data.clean <- data.frame(data_clean_nonpara)
    
    # create plot object to return 
    nonpara_plot <-
      ggplot(data.clean, aes(sample=data.clean))+
      geom_histogram(mapping=aes(x=array(data.clean[,1]), 
                                 y= ..density..), 
                     inherit.aes=FALSE,
                     fill = "white", 
                     color = "blue") +
      xlab(paste("Analyte Values", sep="")) +
      ylab("Density")+
      ggtitle(analyte_nonpara()) + 
      theme(plot.title = element_text(hjust = 0.5))+
      geom_vline(xintercept = nonpara_limits, color= "black")+
      theme_minimal()
    
    output$downloadDataNonpara <- downloadHandler(
      filename = function() {
        paste('nonpara-data',
              "_sex-",sex_nonpara(),
              "_ages-",
              age_lb_nonpara(),"-to-",age_ub_nonpara(),"_",
              Sys.Date(), '.csv', sep='')
      },
      content = function(con) {
        write.csv(my_df_nonpara, con)
      }
    )
    return(nonpara_plot)
  })
  
  
  output$nonpara_plot <- renderPlot({nonpara.plot()})
  
  analyte_refineR <- eventReactive(input$go_refineR, {analyte_refineR <- input$Analyte_refineR})
  unit_refineR <- eventReactive(input$go_refineR, {unit_refineR <- input$unit_refineR})
  age_lb_refineR <- eventReactive(input$go_refineR, {age_lb_refineR <- input$age_lb_refineR})
  age_ub_refineR <- eventReactive(input$go_refineR, {age_ub_refineR <- input$age_ub_refineR})
  n_bs_refineR <- eventReactive(input$go_refineR, {age_ub_refineR <- input$RefineR_N_bootstrap})
  
  ### refineR 
  plot.refineR <- eventReactive(input$go_refineR, {
    
    data.clean = data()[[analyte_refineR()]]
    data.clean = data.clean[data()[["age"]] >= age_lb_refineR()] 
    data.clean = data.clean[data()[["age"]] <= age_ub_refineR()] 
    data.clean = data.clean |> na.omit()
    
    if(sex_refiner() == 'Male'){
      data.clean= data.clean[grep(pattern='^m',
                                  x=data()[[which(tolower(colnames(data()))=='sex')]],ignore.case=TRUE)] 
    } else if (sex_refiner() == 'Female'){
      data.clean= data.clean[grep(pattern='^f',
                                  x=data()[[which(tolower(colnames(data()))=='sex')]],ignore.case=TRUE)] 
    }
    
    RI <- findRI(Data =  data.clean, NBootstrap = floor(n_bs_refineR()))
    plot_refineR <- plot(RI,
                         showCI = (input$plot_CI == "True"))
    
    output$downloadDataRefineR <- downloadHandler(
      filename = function() {
        paste("refineR-data",
              "_sex-",sex_refiner(),
              "_ages-",
              age_lb_refineR(),"-to-",age_ub_refineR(),"_",
              Sys.Date(), '.csv', sep='')
      },
      content = function(con) {
        write.csv(getRI(RI), con)
      }
    )
    
    output$refineR_table <- renderTable(getRI(RI))
    
    return(plot_refineR)
  })
  
  output$refineR_plot <- renderPlot({ plot.refineR() })
  
}

# Run the application 
shinyApp(ui = ui, server = server)