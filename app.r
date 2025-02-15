# library(shiny)
# library(shinydashboard)
# library(tidyverse)
# library(ggplot2)
# library(lubridate)
# library(DT)

# # Define UI
# ui <- dashboardPage(
#   dashboardHeader(title = "VR Puzzle Analysis Dashboard"),
#   dashboardSidebar(
#     sidebarMenu(
#       menuItem("Raw Data", tabName = "raw", icon = icon("table")),
#       menuItem("Demographics", tabName = "demo", icon = icon("users")),
#       menuItem("Duration Analysis", tabName = "duration", icon = icon("clock")),
#       menuItem("Condition Comparison", tabName = "conditions", icon = icon("balance-scale")),
#       menuItem("Controller Metrics", tabName = "controllers", icon = icon("gamepad")),
#       menuItem("Rotational Analysis", tabName = "rotation", icon = icon("sync"))
#     ),
#     selectInput("participant", "Select Participant:", choices = NULL),
#     sliderInput("ageRange", "Age Range:",
#                 min = 19, max = 32, value = c(19, 32))
#   ),
  
#   dashboardBody(
#     tabItems(
#       # Demographic Tab
#       tabItem(tabName = "demo",
#               fluidRow(
#                 box(plotOutput("demoPlot1"), width = 6),
#                 box(plotOutput("demoPlot2"), width = 6)
#               ),
#               fluidRow(
#                 valueBoxOutput("totalParticipants"),
#                 valueBoxOutput("avgAge"),
#                 valueBoxOutput("maleFemaleRatio")
#               )
#       ),
      
#       # Duration Analysis Tab
#       tabItem(tabName = "duration",
#               fluidRow(
#                 box(plotOutput("durationDist"), width = 12)
#               ),
#               fluidRow(
#                 box(plotOutput("conditionDuration"), width = 12)
#               )
#       ),
      
#       # Controller Metrics Tab
#       tabItem(tabName = "controllers",
#               tabBox(width = 12,
#                      tabPanel("Left Controller",
#                               plotOutput("leftSpeedPlot"),
#                               plotOutput("leftRotationPlot")),
#                      tabPanel("Right Controller",
#                               plotOutput("rightSpeedPlot"),
#                               plotOutput("rightRotationPlot"))
#               )
#       ),
      
#       # Condition Comparison Tab
#       tabItem(tabName = "conditions",
#               fluidRow(
#                 box(plotOutput("conditionComparison"), width = 12)
#               )
#       ),
      
#       # Rotational Analysis Tab
#       tabItem(tabName = "rotation",
#               fluidRow(
#                 box(plotOutput("rotationAnalysis"), width = 12)
#               )
#       ),
      
#       # Raw Data Tab
#       tabItem(tabName = "raw",
#               fluidRow(
#                 box(DT::dataTableOutput("rawDataTable"), width = 12)
#               )
#       )
#     )
#   )
# )

# # Server Logic
# server <- function(input, output, session) {
  
#   # Set your main directory (adjust the path as needed)
#   main_directory <- "C:\\Users\\tusha\\Downloads\\PuzzleProjectBIDS\\PuzzleProjectBIDS"
  
#   # Load participants data and filter out anomalies and by age range
#   participants_data <- reactive({
#     data <- read_tsv(file.path(main_directory, "participants.tsv"), 
#                      show_col_types = FALSE) %>% 
#       filter(Anomaly == "No",
#              between(Age, input$ageRange[1], input$ageRange[2]))
#     updateSelectInput(session, "participant", 
#                       choices = unique(data$participant_id))
#     data
#   })
  
#   # Load Main Log data (e.g., _events.tsv files)
#   main_log <- reactive({
#     log_files <- list.files(main_directory, pattern = "\\_events.tsv$", 
#                             recursive = TRUE, full.names = TRUE)
#     map_df(log_files, ~ {
#       subject_id <- str_extract(.x, "sub-\\d+")
#       read_tsv(.x, col_types = cols()) %>% 
#         mutate(participant_id = subject_id)
#     })
#   })
  
#   # Reactive loading of controller data for the selected participant
#   controller_data <- reactive({
#     req(input$participant)
#     list(
#       left = read_controller_data("Left"),
#       right = read_controller_data("Right")
#     )
#   })
  
#   # Function to read controller movement files
#   read_controller_data <- function(type) {
#     pattern <- paste0(type, "ControllerMovement_motion.tsv$")
#     list.files(main_directory, pattern = pattern, 
#                recursive = TRUE, full.names = TRUE) %>% 
#       map_df(~ {
#         subject_id <- str_extract(.x, "sub-\\d+")
#         read_tsv(.x, col_types = cols()) %>% 
#           mutate(participant_id = subject_id)
#       }) %>% 
#       filter(participant_id == input$participant)
#   }
  
#   ## ----------------- Demographics -----------------
#   # Plot: Participants by Race and Gender
#   output$demoPlot1 <- renderPlot({
#     participants_data() %>%
#       group_by(Race, Gender) %>%
#       summarise(Count = n(), .groups = "drop") %>%
#       ggplot(aes(x = Race, y = Count, fill = Gender)) +
#       geom_bar(stat = "identity", position = "stack") +
#       theme_minimal() +
#       labs(title = "Participants by Race and Gender", x = "Race", y = "Count")
#   })
  
#   # Plot: Age Distribution
#   output$demoPlot2 <- renderPlot({
#     ggplot(participants_data(), aes(x = Age)) +
#       geom_histogram(fill = "skyblue", bins = 20) +
#       theme_minimal() +
#       labs(title = "Age Distribution", x = "Age", y = "Count")
#   })
  
#   # Value Box: Total Participants
#   output$totalParticipants <- renderValueBox({
#     valueBox(nrow(participants_data()), "Total Participants", icon = icon("users"))
#   })
  
#   # Value Box: Average Age
#   output$avgAge <- renderValueBox({
#     valueBox(round(mean(participants_data()$Age)), "Average Age", icon = icon("calculator"))
#   })
  
#   # Value Box: Gender Distribution (Males and Females)
#   output$maleFemaleRatio <- renderValueBox({
#     data <- participants_data()
#     num_males <- sum(data$Gender == "Male", na.rm = TRUE)
#     num_females <- sum(data$Gender == "Female", na.rm = TRUE)
#     valueBox(paste("M:", num_males, "F:", num_females), "Gender Distribution", icon = icon("venus-mars"))
#   })
  
#   ## ----------------- Duration Analysis -----------------
#   # Plot: Distribution of Experiment Durations
#   output$durationDist <- renderPlot({
#     ggplot(participants_data(), aes(x = duration)) +
#       geom_histogram(binwidth = 55, fill = "skyblue", color = "black") +
#       labs(title = "Distribution of Experiment Durations", 
#            x = "Duration (s)", y = "Count") +
#       theme_minimal()
#   })
  
#   # Plot: Average Duration by Condition (using Main Log data)
#   output$conditionDuration <- renderPlot({
#     req(input$participant)
#     main_log() %>%
#       filter(participant_id == input$participant) %>%
#       group_by(trial_type) %>%
#       summarise(avg_duration = mean(duration, na.rm = TRUE)) %>%
#       ggplot(aes(x = trial_type, y = avg_duration, fill = trial_type)) +
#       geom_col() +
#       theme_minimal() +
#       labs(title = "Average Duration by Condition", 
#            x = "Condition", y = "Average Duration (s)")
#   })
  
#   ## ----------------- Condition Comparison -----------------
#   # Plot: Boxplot comparing duration across conditions
#   output$conditionComparison <- renderPlot({
#     req(input$participant)
#     main_log() %>%
#       filter(participant_id == input$participant) %>%
#       ggplot(aes(x = trial_type, y = duration, fill = trial_type)) +
#       geom_boxplot() +
#       theme_minimal() +
#       labs(title = "Duration Comparison by Condition", 
#            x = "Condition", y = "Duration (s)")
#   })
  
#   ## ----------------- Controller Metrics -----------------
#   # Left Controller Speed Plot
#   output$leftSpeedPlot <- renderPlot({
#     req(controller_data()$left)
#     calculate_speed_metrics(controller_data()$left) %>%
#       ggplot(aes(x = time, y = speed)) +
#       geom_line(color = "blue") +
#       labs(title = "Left Controller Speed Over Time", x = "Time", y = "Speed")
#   })
  
#   # Right Controller Speed Plot
#   output$rightSpeedPlot <- renderPlot({
#     req(controller_data()$right)
#     calculate_speed_metrics(controller_data()$right) %>%
#       ggplot(aes(x = time, y = speed)) +
#       geom_line(color = "green") +
#       labs(title = "Right Controller Speed Over Time", x = "Time", y = "Speed")
#   })
  
#   # Function to calculate linear speed (using positional differences)
#   calculate_speed_metrics <- function(data) {
#     data %>%
#       mutate(
#         delta_t = time - lag(time),
#         speed = sqrt((pos_x - lag(pos_x))^2 +
#                      (pos_y - lag(pos_y))^2 +
#                      (pos_z - lag(pos_z))^2) / delta_t
#       ) %>%
#       filter(!is.na(speed))
#   }
  
#   # Left Controller Angular Rotation Plot
#   output$leftRotationPlot <- renderPlot({
#     req(controller_data()$left)
#     calculate_rotation_metrics(controller_data()$left) %>%
#       ggplot(aes(x = time, y = angular_speed)) +
#       geom_line(color = "blue") +
#       labs(title = "Left Controller Angular Speed Over Time", x = "Time", y = "Angular Speed")
#   })
  
#   # Right Controller Angular Rotation Plot
#   output$rightRotationPlot <- renderPlot({
#     req(controller_data()$right)
#     calculate_rotation_metrics(controller_data()$right) %>%
#       ggplot(aes(x = time, y = angular_speed)) +
#       geom_line(color = "green") +
#       labs(title = "Right Controller Angular Speed Over Time", x = "Time", y = "Angular Speed")
#   })
  
#   # Function to calculate (simplified) angular speed using rotation differences
#   calculate_rotation_metrics <- function(data) {
#     data %>%
#       mutate(
#         delta_t = time - lag(time),
#         angular_change = sqrt((rot_x - lag(rot_x))^2 +
#                               (rot_y - lag(rot_y))^2 +
#                               (rot_z - lag(rot_z))^2),
#         angular_speed = angular_change / delta_t
#       ) %>%
#       filter(!is.na(angular_speed))
#   }
  
#   ## ----------------- Rotational Analysis -----------------
#   # Plot: Distribution of angular speeds combining both controllers
#   output$rotationAnalysis <- renderPlot({
#     req(controller_data())
#     left_data <- calculate_rotation_metrics(controller_data()$left) %>% mutate(controller = "Left")
#     right_data <- calculate_rotation_metrics(controller_data()$right) %>% mutate(controller = "Right")
#     combined <- bind_rows(left_data, right_data)
#     ggplot(combined, aes(x = angular_speed, fill = controller)) +
#       geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
#       labs(title = "Distribution of Angular Speed", x = "Angular Speed", y = "Count") +
#       theme_minimal()
#   })
  
#   ## ----------------- Raw Data Tab -----------------
#   output$rawDataTable <- DT::renderDataTable({
#     DT::datatable(participants_data(), options = list(pageLength = 10))
#   })
# }

# # Run the application
# shinyApp(ui, server)

library(shiny)
library(shinydashboard)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(DT)
library(stringr)
library(purrr)

# Define UI with several tabs to reflect the R Markdown outputs
ui <- dashboardPage(
  dashboardHeader(title = "VR Puzzle Analysis Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Raw Data", tabName = "raw", icon = icon("table")),
      menuItem("Demographics", tabName = "demo", icon = icon("users")),
      menuItem("Duration Analysis", tabName = "duration", icon = icon("clock")),
      menuItem("Condition Comparison", tabName = "conditions", icon = icon("balance-scale")),
      menuItem("Controller Metrics", tabName = "controllers", icon = icon("gamepad")),
      menuItem("Rotational Analysis", tabName = "rotation", icon = icon("sync"))
    ),
    # In some plots we filter by age and also by participant
    selectInput("participant", "Select Participant:", choices = NULL),
    sliderInput("ageRange", "Age Range:",
                min = 19, max = 32, value = c(19, 32))
  ),
  
  dashboardBody(
    tabItems(
      # Demographics Tab
      tabItem(tabName = "demo",
              fluidRow(
                box(plotOutput("demoPlot1"), width = 6),
                box(plotOutput("demoPlot2"), width = 6)
              )
      ),
      
      # Duration Analysis Tab
      tabItem(tabName = "duration",
              fluidRow(
                box(plotOutput("durationDist"), width = 12)
              )
      ),
      
      # Condition Comparison Tab – using a tabBox with three sub-tabs
      tabItem(tabName = "conditions",
              tabBox(width = 12,
                     tabPanel("One-step Movement", plotOutput("conditionOneStepPlot")),
                     tabPanel("Two-step Movement", plotOutput("conditionTwoStepPlot")),
                     tabPanel("Overall Conditions", plotOutput("overallConditionPlot"))
              )
      ),
      
      # Controller Metrics Tab – overall controller speed plots
      tabItem(tabName = "controllers",
              tabBox(width = 12,
                     tabPanel("Left Controller Speed", plotOutput("leftControllerSpeedPlot")),
                     tabPanel("Right Controller Speed", plotOutput("rightControllerSpeedPlot"))
              )
      ),
      
      # Rotational Analysis Tab – condition‐specific rotational speed plots for left/right controllers
      tabItem(tabName = "rotation",
              tabBox(width = 12,
                     tabPanel("Left Controller Rotation",
                              fluidRow(
                                box(plotOutput("leftRotationOneStepPlot"), width = 6),
                                box(plotOutput("leftRotationTwoStepPlot"), width = 6)
                              )),
                     tabPanel("Right Controller Rotation",
                              fluidRow(
                                box(plotOutput("rightRotationOneStepPlot"), width = 6),
                                box(plotOutput("rightRotationTwoStepPlot"), width = 6)
                              ))
              )
      ),
      
      # Raw Data Tab
      tabItem(tabName = "raw",
              fluidRow(
                box(DT::dataTableOutput("rawDataTable"), width = 12)
              )
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  
  # Set your main directory (adjust the path as needed; note the PDF used "D:/IIT/PuzzleProject/restructuredDataset")
  main_directory <- "C:\\Users\\tusha\\Downloads\\PuzzleProjectBIDS\\PuzzleProjectBIDS"
  
  # Load all participants (do not filter yet)
  all_participants <- reactive({
    read_tsv(file.path(main_directory, "participants.tsv"), show_col_types = FALSE)
  })
  
  # Separate anomaly participants and non-anomalies
  anomaly_data <- reactive({
    all_participants() %>% filter(Anomaly == "Yes")
  })
  
  participants_data <- reactive({
    all_participants() %>% 
      filter(Anomaly == "No") %>%
      filter(between(Age, input$ageRange[1], input$ageRange[2]))
  })
  
  # Update participant select input based on non-anomalous participants
  observe({
    updateSelectInput(session, "participant", choices = unique(participants_data()$participant_id))
  })
  
  # Load Main Log data (combining all _events.tsv files)
  MainLog <- reactive({
    log_files <- list.files(main_directory, pattern = "\\_events.tsv$", recursive = TRUE, full.names = TRUE)
    map_df(log_files, ~ {
      subject_id <- str_extract(.x, "sub-\\w+")
      read_tsv(.x, col_types = cols()) %>% 
        mutate(participant_id = subject_id)
    })
  })
  
  # Load motion data for Left and Right Controllers
  LeftControllerData <- reactive({
    files <- list.files(main_directory, pattern = "LeftControllerMovement_motion.tsv$", recursive = TRUE, full.names = TRUE)
    if(length(files)==0) return(tibble())
    map_df(files, function(f) {
      subject_id <- str_extract(f, "sub-\\w+")
      read_tsv(f, col_types = cols()) %>% mutate(participant_id = subject_id)
    })
  })
  
  RightControllerData <- reactive({
    files <- list.files(main_directory, pattern = "RightControllerMovement_motion.tsv$", recursive = TRUE, full.names = TRUE)
    if(length(files)==0) return(tibble())
    map_df(files, function(f) {
      subject_id <- str_extract(f, "sub-\\w+")
      read_tsv(f, col_types = cols()) %>% mutate(participant_id = subject_id)
    })
  })
  
  ## Define condition-specific logs (for conditions 0,1,2,3)
  condition_0_log <- reactive({
    MainLog() %>% filter(trial_type == "Condition: 0") %>%
      mutate(participant_number = as.numeric(str_extract(participant_id, "\\d+"))) %>%
      arrange(participant_number)
  })
  condition_1_log <- reactive({
    MainLog() %>% filter(trial_type == "Condition: 1") %>%
      mutate(participant_number = as.numeric(str_extract(participant_id, "\\d+"))) %>%
      arrange(participant_number)
  })
  condition_2_log <- reactive({
    MainLog() %>% filter(trial_type == "Condition: 2") %>%
      mutate(participant_number = as.numeric(str_extract(participant_id, "\\d+"))) %>%
      arrange(participant_number)
  })
  condition_3_log <- reactive({
    MainLog() %>% filter(trial_type == "Condition: 3") %>%
      mutate(participant_number = as.numeric(str_extract(participant_id, "\\d+"))) %>%
      arrange(participant_number)
  })
  
  # For the duration plot, compute min, max, avg durations
  min_duration <- reactive({ min(participants_data()$duration, na.rm = TRUE) })
  max_duration <- reactive({ max(participants_data()$duration, na.rm = TRUE) })
  avg_duration <- reactive({ mean(participants_data()$duration, na.rm = TRUE) })
  
  ## ---- Functions from the R Markdown file ----
  
  # Function to calculate speed in a given time interval (for condition-specific analysis)
  calculate_speed_in_interval <- function(motion_data, start_time, end_time) {
    interval_data <- motion_data %>%
      filter(time >= start_time, time <= end_time) %>%
      mutate(
        delta_t = time - lag(time),
        delta_x = pos_x - lag(pos_x),
        delta_y = pos_y - lag(pos_y),
        delta_z = pos_z - lag(pos_z),
        speed = sqrt(delta_x^2 + delta_y^2 + delta_z^2) / delta_t
      ) %>%
      filter(!is.na(speed) & delta_t > 0)
    round(mean(interval_data$speed, na.rm = TRUE), 4)
  }
  
  # Function to normalize quaternions (from rotation data in degrees)
  normalize_quaternion <- function(x, y, z) {
    x <- x * pi / 180
    y <- y * pi / 180
    z <- z * pi / 180
    qx = sin(x / 2) * cos(y / 2) * cos(z / 2) - cos(x / 2) * sin(y / 2) * sin(z / 2)
    qy = cos(x / 2) * sin(y / 2) * cos(z / 2) + sin(x / 2) * cos(y / 2) * sin(z / 2)
    qz = cos(x / 2) * cos(y / 2) * sin(z / 2) - sin(x / 2) * sin(y / 2) * cos(z / 2)
    qw = cos(x / 2) * cos(y / 2) * cos(z / 2) + sin(x / 2) * sin(y / 2) * sin(z / 2)
    tibble(w = qw, x = qx, y = qy, z = qz)
  }
  
  # Function to calculate angular speed given two quaternion values
  calculate_angular_speed <- function(q1, q2) {
    q1_inv <- tibble(w = q1$w, x = -q1$x, y = -q1$y, z = -q1$z)
    delta_w <- q1_inv$w * q2$w - q1_inv$x * q2$x - q1_inv$y * q2$y - q1_inv$z * q2$z
    delta_w <- pmin(pmax(delta_w, -1), 1)
    theta <- 2 * acos(delta_w)
    angular_speed <- theta * (180 / pi)
    angular_speed
  }
  
  # Function to calculate average rotational speed in a time interval
  calculate_rotational_speed_in_interval <- function(motion_data, start_time, end_time) {
    interval_data <- motion_data %>%
      filter(time >= start_time & time <= end_time) %>%
      mutate(
        delta_t = time - lag(time),
        lag_rot_x = lag(rot_x),
        lag_rot_y = lag(rot_y),
        lag_rot_z = lag(rot_z)
      ) %>%
      filter(!is.na(delta_t) & delta_t > 0)
    q1 <- normalize_quaternion(interval_data$lag_rot_x, interval_data$lag_rot_y, interval_data$lag_rot_z)
    q2 <- normalize_quaternion(interval_data$rot_x, interval_data$rot_y, interval_data$rot_z)
    rotational_speeds <- calculate_angular_speed(q1, q2) / interval_data$delta_t
    avg_rotational_speed <- mean(rotational_speeds, na.rm = TRUE)
    avg_rotational_speed
  }
  
  ## ---------------- Demographics ----------------
  # Plot 1: Participants by Race and Gender (using the color scheme from the PDF)
  output$demoPlot1 <- renderPlot({
    summary_data <- participants_data() %>%
      group_by(Race, Gender) %>%
      summarise(Count = n(), .groups = "drop")
    ggplot(summary_data, aes(x = Race, y = Count, fill = Gender)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = c("Male" = "skyblue", "Female" = "pink", "Non-binary" = "yellow")) +
      theme_minimal() +
      labs(title = "Participants by Race and Gender",
           x = "Race",
           y = "Number of Participants",
           fill = "Gender") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Plot 2: Age Distribution Histogram
  output$demoPlot2 <- renderPlot({
    ggplot(participants_data(), aes(x = Age)) +
      geom_histogram(fill = "lightgreen", bins = 10, color = "black") +
      theme_minimal() +
      labs(title = "Age Distribution",
           x = "Age",
           y = "Count")
  })
  
  ## ---------------- Duration Analysis ----------------
  output$durationDist <- renderPlot({
    # Convert duration from seconds to minutes for plotting and annotation
    ggplot(participants_data(), aes(x = duration/60)) +
      geom_histogram(binwidth = 55/60, fill = "skyblue", color = "black") +
      labs(title = "Distribution of Experiment Durations in Minutes",
           subtitle = paste(nrow(participants_data()), "Participants"),
           x = "Duration (Minutes)", y = "Count") +
      geom_vline(xintercept = min_duration()/60, linetype = "dashed") +
      geom_vline(xintercept = max_duration()/60, linetype = "dashed") +
      geom_vline(xintercept = avg_duration()/60, linetype = "solid") +
      annotate("text", x = min_duration()/60, y = 5, label = paste("Min:", round(min_duration()/60, 2)), angle = 90, vjust = -0.5) +
      annotate("text", x = max_duration()/60, y = 5, label = paste("Max:", round(max_duration()/60, 2)), angle = 90, vjust = -0.5) +
      annotate("text", x = avg_duration()/60, y = 5, label = paste("Avg:", round(avg_duration()/60, 2)), angle = 90, vjust = -0.5) +
      theme_minimal()
  })
  
  ## ---------------- Condition Comparison ----------------
  # One-step Movement (Conditions 0 & 2)
  output$conditionOneStepPlot <- renderPlot({
    combined_data <- bind_rows(
      condition_0_log() %>% mutate(condition = "Condition 0"),
      condition_2_log() %>% mutate(condition = "Condition 2")
    ) %>% select(participant_number, duration, condition)
    combined_data$condition <- factor(combined_data$condition, levels = c("Condition 0", "Condition 2"))
    line_colors <- combined_data %>%
      pivot_wider(names_from = condition, values_from = duration) %>%
      mutate(color_flag = ifelse(`Condition 0` > `Condition 2`, "highlight", "normal")) %>%
      select(participant_number, color_flag)
    combined_data <- left_join(combined_data, line_colors, by = "participant_number")
    ggplot(combined_data, aes(x = condition, y = duration, group = participant_number)) +
      geom_line(aes(color = color_flag)) +
      geom_point(aes(color = color_flag), size = 2) +
      stat_summary(fun = mean, geom = "point", aes(group = 1), size = 5, color = "red") +
      scale_color_manual(values = c("highlight" = "red", "normal" = "black")) +
      scale_x_discrete(labels = c("Condition 0" = "One-step without clutter",
                                  "Condition 2" = "One-step with clutter")) +
      labs(title = "Time to Solve Puzzle (One-step Movement)",
           x = "Condition",
           y = "Time (seconds)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  })
  
  # Two-step Movement (Conditions 1 & 3)
  output$conditionTwoStepPlot <- renderPlot({
    combined_data <- bind_rows(
      condition_1_log() %>% mutate(condition = "Condition 1"),
      condition_3_log() %>% mutate(condition = "Condition 3")
    ) %>% select(participant_number, duration, condition)
    combined_data$condition <- factor(combined_data$condition, levels = c("Condition 1", "Condition 3"))
    # Identify outliers from anomaly data
    outliers <- str_extract(anomaly_data()$participant_id, "\\d+")
    line_colors <- combined_data %>%
      pivot_wider(names_from = condition, values_from = duration) %>%
      mutate(color_flag = ifelse(participant_number %in% as.numeric(outliers), "highlight", "normal")) %>%
      select(participant_number, color_flag)
    combined_data <- left_join(combined_data, line_colors, by = "participant_number")
    ggplot(combined_data, aes(x = condition, y = duration, group = participant_number)) +
      geom_line(aes(color = color_flag)) +
      geom_point(size = 2) +
      stat_summary(fun = mean, geom = "point", aes(group = 1), size = 5, color = "red") +
      scale_color_manual(values = c("highlight" = "red", "normal" = "black")) +
      scale_x_discrete(labels = c("Condition 1" = "Two-step without clutter",
                                  "Condition 3" = "Two-step with clutter")) +
      labs(title = "Time to Solve Puzzle (Two-step Movement)",
           x = "Condition",
           y = "Time (seconds)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  })
  
  # Overall Conditions (combining all conditions by clutter and step type)
  output$overallConditionPlot <- renderPlot({
    allConditions <- MainLog() %>%
      mutate(
        clutter_group = case_when(
          trial_type %in% c("Condition: 0", "Condition: 1") ~ "Open",
          trial_type %in% c("Condition: 2", "Condition: 3") ~ "Cluttered"
        ),
        step_type = case_when(
          trial_type %in% c("Condition: 0", "Condition: 2") ~ "One Step",
          trial_type %in% c("Condition: 1", "Condition: 3") ~ "Two Steps"
        )
      )
    avg_duration_conditions <- allConditions %>%
      group_by(clutter_group, step_type) %>%
      summarise(avg_duration = mean(duration), .groups = "drop")
    slope_data <- avg_duration_conditions %>%
      pivot_wider(names_from = clutter_group, values_from = avg_duration) %>%
      rename(WithoutClutter = Open, WithClutter = Cluttered)
    ggplot() +
      geom_point(data = allConditions, aes(x = clutter_group, y = duration, shape = clutter_group, color = step_type), 
                 alpha = 0.2, size = 1.5) +
      geom_point(data = avg_duration_conditions, aes(x = clutter_group, y = avg_duration, shape = clutter_group, color = step_type), 
                 size = 3) +
      geom_segment(data = slope_data, aes(x = 1, xend = 2, y = WithoutClutter, yend = WithClutter, color = step_type), 
                   size = 1) +
      scale_shape_manual(values = c("Open" = 17, "Cluttered" = 16)) +
      scale_color_manual(values = c("One Step" = "red", "Two Steps" = "blue")) +
      scale_x_discrete(limits = c("Open", "Cluttered")) +
      labs(title = "Time to Solve Puzzle in Different Conditions",
           x = "",
           y = "Time to Solve (seconds)",
           shape = "Desk",
           color = "Steps") +
      theme_minimal()
  })
  
  ## ---------------- Controller Metrics ----------------
  # Overall average left controller speed by participant
  output$leftControllerSpeedPlot <- renderPlot({
    motion_data_left <- LeftControllerData() %>%
      group_by(participant_id) %>%
      mutate(
        delta_t = time - lag(time),
        delta_x = pos_x - lag(pos_x),
        delta_y = pos_y - lag(pos_y),
        delta_z = pos_z - lag(pos_z),
        speed = sqrt(delta_x^2 + delta_y^2 + delta_z^2) / delta_t,
        participant_number = as.numeric(str_extract(participant_id, "\\d+"))
      ) %>%
      filter(!is.na(speed) & delta_t > 0)
    average_speed_left <- motion_data_left %>%
      group_by(participant_number) %>%
      summarise(avg_speed = mean(speed, na.rm = TRUE))
    ggplot(average_speed_left, aes(x = participant_number, y = avg_speed)) +
      geom_point(size = 3, color = "black") +
      labs(title = "Average Left Controller Speed by Participant",
           x = "Participant ID",
           y = "Average Speed (m/s)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Overall average right controller speed by participant
  output$rightControllerSpeedPlot <- renderPlot({
    motion_data_right <- RightControllerData() %>%
      group_by(participant_id) %>%
      mutate(
        delta_t = time - lag(time),
        delta_x = pos_x - lag(pos_x),
        delta_y = pos_y - lag(pos_y),
        delta_z = pos_z - lag(pos_z),
        speed = sqrt(delta_x^2 + delta_y^2 + delta_z^2) / delta_t,
        participant_number = as.numeric(str_extract(participant_id, "\\d+"))
      ) %>%
      filter(!is.na(speed) & delta_t > 0)
    average_speed_right <- motion_data_right %>%
      group_by(participant_number) %>%
      summarise(avg_speed = mean(speed, na.rm = TRUE))
    ggplot(average_speed_right, aes(x = participant_number, y = avg_speed)) +
      geom_point(size = 3, color = "black") +
      labs(title = "Average Right Controller Speed by Participant",
           x = "Participant ID",
           y = "Average Speed (m/s)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  ## ---------------- Rotational Analysis ----------------
  # Left Controller Rotation – One-step (Conditions 0 vs 2)
  output$leftRotationOneStepPlot <- renderPlot({
    participants_list <- unique(LeftControllerData()$participant_id)
    result <- map_dfr(participants_list, function(participant) {
      participant_log <- MainLog() %>% 
        filter(participant_id == participant, trial_type %in% c("Condition: 0", "Condition: 2"))
      start_times_condition_0 <- participant_log$onset[participant_log$trial_type == "Condition: 0"]
      end_times_condition_0 <- participant_log$end_timestamp[participant_log$trial_type == "Condition: 0"]
      start_times_condition_2 <- participant_log$onset[participant_log$trial_type == "Condition: 2"]
      end_times_condition_2 <- participant_log$end_timestamp[participant_log$trial_type == "Condition: 2"]
      participant_motion_data <- LeftControllerData() %>% filter(participant_id == participant)
      avg_speed_condition_0 <- calculate_rotational_speed_in_interval(participant_motion_data, start_times_condition_0, end_times_condition_0)
      avg_speed_condition_2 <- calculate_rotational_speed_in_interval(participant_motion_data, start_times_condition_2, end_times_condition_2)
      tibble(
        participant_id = rep(as.numeric(str_extract(participant, "\\d+")), 2),
        avg_speed = c(avg_speed_condition_0, avg_speed_condition_2),
        condition = rep(c("Condition 0", "Condition 2"), 1)
      )
    })
    ggplot(result, aes(x = condition, y = avg_speed, group = participant_id)) +
      geom_line() +
      geom_point(size = 2) +
      stat_summary(fun = mean, geom = "point", aes(group = 1), size = 5, color = "red") +
      scale_x_discrete(labels = c("Condition 0" = "One-step without clutter",
                                  "Condition 2" = "One-step with clutter")) +
      labs(title = "Left Controller Rotational Speed (One-step)",
           x = "Condition",
           y = "Speed (°/s)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  })
  
  # Left Controller Rotation – Two-step (Conditions 1 vs 3)
  output$leftRotationTwoStepPlot <- renderPlot({
    participants_list <- unique(LeftControllerData()$participant_id)
    result <- map_dfr(participants_list, function(participant) {
      participant_log <- MainLog() %>% 
        filter(participant_id == participant, trial_type %in% c("Condition: 1", "Condition: 3"))
      start_times_condition_1 <- participant_log$onset[participant_log$trial_type == "Condition: 1"]
      end_times_condition_1 <- participant_log$end_timestamp[participant_log$trial_type == "Condition: 1"]
      start_times_condition_3 <- participant_log$onset[participant_log$trial_type == "Condition: 3"]
      end_times_condition_3 <- participant_log$end_timestamp[participant_log$trial_type == "Condition: 3"]
      participant_motion_data <- LeftControllerData() %>% filter(participant_id == participant)
      avg_speed_condition_1 <- calculate_rotational_speed_in_interval(participant_motion_data, start_times_condition_1, end_times_condition_1)
      avg_speed_condition_3 <- calculate_rotational_speed_in_interval(participant_motion_data, start_times_condition_3, end_times_condition_3)
      tibble(
        participant_id = rep(as.numeric(str_extract(participant, "\\d+")), 2),
        avg_speed = c(avg_speed_condition_1, avg_speed_condition_3),
        condition = rep(c("Condition 1", "Condition 3"), 1)
      )
    })
    ggplot(result, aes(x = condition, y = avg_speed, group = participant_id)) +
      geom_line() +
      geom_point(size = 2) +
      stat_summary(fun = mean, geom = "point", aes(group = 1), size = 5, color = "red") +
      scale_x_discrete(labels = c("Condition 1" = "Two-step without clutter",
                                  "Condition 3" = "Two-step with clutter")) +
      labs(title = "Left Controller Rotational Speed (Two-step)",
           x = "Condition",
           y = "Speed (°/s)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  })
  
  # Right Controller Rotation – One-step (Conditions 0 vs 2)
  output$rightRotationOneStepPlot <- renderPlot({
    participants_list <- unique(RightControllerData()$participant_id)
    result <- map_dfr(participants_list, function(participant) {
      participant_log <- MainLog() %>% 
        filter(participant_id == participant, trial_type %in% c("Condition: 0", "Condition: 2"))
      start_times_condition_0 <- participant_log$onset[participant_log$trial_type == "Condition: 0"]
      end_times_condition_0 <- participant_log$end_timestamp[participant_log$trial_type == "Condition: 0"]
      start_times_condition_2 <- participant_log$onset[participant_log$trial_type == "Condition: 2"]
      end_times_condition_2 <- participant_log$end_timestamp[participant_log$trial_type == "Condition: 2"]
      participant_motion_data <- RightControllerData() %>% filter(participant_id == participant)
      avg_speed_condition_0 <- calculate_rotational_speed_in_interval(participant_motion_data, start_times_condition_0, end_times_condition_0)
      avg_speed_condition_2 <- calculate_rotational_speed_in_interval(participant_motion_data, start_times_condition_2, end_times_condition_2)
      tibble(
        participant_id = rep(as.numeric(str_extract(participant, "\\d+")), 2),
        avg_speed = c(avg_speed_condition_0, avg_speed_condition_2),
        condition = rep(c("Condition 0", "Condition 2"), 1)
      )
    })
    ggplot(result, aes(x = condition, y = avg_speed, group = participant_id)) +
      geom_line() +
      geom_point(size = 2) +
      stat_summary(fun = mean, geom = "point", aes(group = 1), size = 5, color = "red") +
      scale_x_discrete(labels = c("Condition 0" = "One-step without clutter",
                                  "Condition 2" = "One-step with clutter")) +
      labs(title = "Right Controller Rotational Speed (One-step)",
           x = "Condition",
           y = "Speed (°/s)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  })
  
  # Right Controller Rotation – Two-step (Conditions 1 vs 3)
  output$rightRotationTwoStepPlot <- renderPlot({
    participants_list <- unique(RightControllerData()$participant_id)
    result <- map_dfr(participants_list, function(participant) {
      participant_log <- MainLog() %>% 
        filter(participant_id == participant, trial_type %in% c("Condition: 1", "Condition: 3"))
      start_times_condition_1 <- participant_log$onset[participant_log$trial_type == "Condition: 1"]
      end_times_condition_1 <- participant_log$end_timestamp[participant_log$trial_type == "Condition: 1"]
      start_times_condition_3 <- participant_log$onset[participant_log$trial_type == "Condition: 3"]
      end_times_condition_3 <- participant_log$end_timestamp[participant_log$trial_type == "Condition: 3"]
      participant_motion_data <- RightControllerData() %>% filter(participant_id == participant)
      avg_speed_condition_1 <- calculate_rotational_speed_in_interval(participant_motion_data, start_times_condition_1, end_times_condition_1)
      avg_speed_condition_3 <- calculate_rotational_speed_in_interval(participant_motion_data, start_times_condition_3, end_times_condition_3)
      tibble(
        participant_id = rep(as.numeric(str_extract(participant, "\\d+")), 2),
        avg_speed = c(avg_speed_condition_1, avg_speed_condition_3),
        condition = rep(c("Condition 1", "Condition 3"), 1)
      )
    })
    ggplot(result, aes(x = condition, y = avg_speed, group = participant_id)) +
      geom_line() +
      geom_point(size = 2) +
      stat_summary(fun = mean, geom = "point", aes(group = 1), size = 5, color = "red") +
      scale_x_discrete(labels = c("Condition 1" = "Two-step without clutter",
                                  "Condition 3" = "Two-step with clutter")) +
      labs(title = "Right Controller Rotational Speed (Two-step)",
           x = "Condition",
           y = "Speed (°/s)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  })
  
  ## ---------------- Raw Data ----------------
  output$rawDataTable <- DT::renderDataTable({
    DT::datatable(participants_data(), options = list(pageLength = 10))
  })
  
}

# Run the application
shinyApp(ui, server)

