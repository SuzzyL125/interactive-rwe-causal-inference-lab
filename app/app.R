library(shiny)
library(ggplot2)
root <- if (file.exists("R/simulate_cohort.R")) "." else ".."
source(file.path(root,"R/simulate_cohort.R")); source(file.path(root,"R/causal_methods.R")); source(file.path(root,"R/plots.R"))

ui <- fluidPage(
  titlePanel("Interactive RWE Causal Inference Lab"),
  sidebarLayout(
    sidebarPanel(
      numericInput("n","Cohort size",10000,1000,100000,1000),
      sliderInput("txprev","Drug A prevalence",.1,.9,.4,.05),
      sliderInput("age","Age → treatment selection",0,1.5,.35,.05),
      sliderInput("comorb","Comorbidity → treatment selection",0,1.5,.55,.05),
      sliderInput("severity","Severity → treatment selection",0,2,.8,.05),
      sliderInput("unmeasured","Unmeasured confounding",0,1.5,0,.05),
      sliderInput("positivity","Positivity stress",.5,3,1,.1),
      actionButton("run","Run scenario",class="btn-primary"),
      helpText("Synthetic educational data; not clinical evidence.")
    ),
    mainPanel(tabsetPanel(
      tabPanel("Cohort",uiOutput("cards"),tableOutput("table1")),
      tabPanel("Balance & overlap",uiOutput("warning"),plotOutput("love"),plotOutput("overlap")),
      tabPanel("Treatment effect",tableOutput("effects"),plotOutput("forest")),
      tabPanel("What if the method is wrong?",uiOutput("bias_text"),plotOutput("bias")),
      tabPanel("Methods",h3("Analysis choices"),p("Crude, regression standardization, propensity-score matching, stabilized IPTW, and doubly robust AIPW are compared on the risk-difference scale."),h3("Core lesson"),p("Good measured balance does not eliminate bias from unmeasured confounding or severe positivity violations."))
    ))
  )
)

server <- function(input,output,session) {
  scenario <- eventReactive(input$run,{
    d <- simulate_rwe_cohort(input$n,input$txprev,input$age,input$comorb,input$severity,input$unmeasured,input$positivity)
    a <- estimate_effects(d); list(d=d,a=a,b=balance_table(d,a$ps,a$weights))
  },ignoreNULL=FALSE)
  output$cards <- renderUI({x<-scenario();d<-x$d;fluidRow(column(4,wellPanel(h3(format(nrow(d),big.mark=",")),p("Patients"))),column(4,wellPanel(h3(sum(d$treatment)),p("Drug A"))),column(4,wellPanel(h3(sum(d$outcome)),p("Outcome events"))))})
  output$table1 <- renderTable({d<-scenario()$d;data.frame(Group=c("Drug A","Drug B"),N=c(sum(d$treatment),sum(!d$treatment)),Mean_age=c(mean(d$age[d$treatment==1]),mean(d$age[d$treatment==0])),Comorbidity=c(mean(d$comorbidity[d$treatment==1]),mean(d$comorbidity[d$treatment==0])),Outcome_percent=100*c(mean(d$outcome[d$treatment==1]),mean(d$outcome[d$treatment==0])))},digits=1)
  output$warning <- renderUI({
    x <- scenario(); poor <- mean(x$a$ps < .05 | x$a$ps > .95) > .05
    if (poor) {
      div(class="alert alert-warning",strong("Limited positivity detected.")," More than 5% of estimated scores are near 0 or 1.")
    } else {
      div(class="alert alert-success","Propensity-score overlap is acceptable under this diagnostic.")
    }
  })
  output$love <- renderPlot(plot_love(scenario()$b)); output$overlap <- renderPlot(plot_overlap(scenario()$d,scenario()$a$ps))
  output$effects <- renderTable({x<-scenario()$a$effects;x[,c("method","estimate","lower","upper")]},digits=3)
  output$forest <- renderPlot(plot_forest(scenario()$a$effects,scenario()$a$true_ate))
  output$bias_text <- renderUI({x<-scenario();div(class="alert alert-info",sprintf("True simulated ATE: %.3f. Increasing unmeasured confounding shows why propensity scores cannot adjust variables that were not measured.",x$a$true_ate))})
  output$bias <- renderPlot(plot_bias(scenario()$a$effects,scenario()$a$true_ate))
}
shinyApp(ui,server)
