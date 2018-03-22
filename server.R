setwd("/Users/hhsieh/Documents/ANTABIS/RASp/RAS species list/Three Bigs")
data <- read.csv("data_m.csv", sep = ",", header = T, row.names = NULL)
#data <- read.csv("ThreeBigs_n_2.0.csv", sep = ",", header = T, row.names = NULL)
#data_s <- subset(data, Match.type != "")
#data_s <- subset(data_s, Kingdom !="" & Kingdom != "P. Micheli ex Haller")
#library(stringr)

#Authority_year <- function(Authority_accepted) {
#  regexp <- "[[:digit:]]+"
#  return(str_extract(Authority_accepted, regexp))
#}

#Authority <- as.character(data_s$Authority_accepted)
#year <- unlist(lapply(Authority, Authority_year))
#data_m <- data.frame(data_s, year)
#data_m <- subset(data_m, year!= "")

library(data.table)
#data_m <- data.frame(data_s, year)
#data_m <- subset(data_m, year!= "")[,c(16,19:24,44)]
#colnames(data_m) <- c("AphiaIDs", "Kingdoms", "Phyla", "Classes", "Orders", "Families", "Genera", "year")
#write.csv(data_m, "data_m.csv", row.names = FALSE)

library(shiny)

shinyServer(function(input, output) {
  viewData <- reactive({
    df <- subset(data_m, Kingdoms == input$taxa | Phyla == input$taxa | Classes == input$taxa | Orders == input$taxa
                 | Families == input$taxa | Genera == input$taxa)
  })

  #taxaaccum <- reactive({
  #  df <- subset(data_m, Kingdoms == input$taxa | Phyla == input$taxa | Classes == input$taxa | Orders == input$taxa
  #               | Families == input$taxa | Genera == input$taxa)
  #  dt = as.data.table(unique(df))
  #  setkey(dt, "year")
  #  if (input$rank == "Phylum" | input$rank == "phylum") {
  #    dt[, id := as.numeric(factor(Phyla, levels = unique(Phyla)))]
  #  } else if (input$rank == "Class" | input$rank == "class") {
  #    dt[, id := as.numeric(factor(Classes, levels = unique(Classes)))]
  #  } else if (input$rank == "Order" | input$rank == "order") {
  #    dt[, id := as.numeric(factor(Orders, levels = unique(Orders)))]
  #  } else if (input$rank == "Family" | input$rank == "family") {
  #    dt[, id := as.numeric(factor(Families, levels = unique(Families)))]
  #  } else if (input$rank == "Genus" | input$rank == "genus") {
  #    dt[, id := as.numeric(factor(Genera, levels = unique(Genera)))]
  #  } else if (input$rank == "Species" | input$rank == "species") {
  #    dt[, id := as.numeric(factor(AphiaIDs, levels = unique(AphiaIDs)))]
  #  }

  #  dt.out <- dt[J(unique(year)), mult = "last"]#[, Phylum := NULL]
  #  dt.out[, id := cummax(id)]
  #  numtaxa <- cummax(as.numeric(factor(dt$id)))
  #  taxa_dt <- aggregate(numtaxa, list(year = dt$year), max)
  #  colnames(taxa_dt) <- c("year", "taxa count")
  #  taxa_dt <- taxa_dt
  #})


  modelfit <- reactive({
    df <- subset(data_m, Kingdoms == input$taxa | Phyla == input$taxa | Classes == input$taxa | Orders == input$taxa
                 | Families == input$taxa | Genera == input$taxa)
    dt = as.data.table(unique(df))
    setkey(dt, "year")
    if (input$rank == "Phylum" | input$rank == "phylum") {
      dt[, id := as.numeric(factor(Phyla, levels = unique(Phyla)))]
    } else if (input$rank == "Class" | input$rank == "class") {
      dt[, id := as.numeric(factor(Classes, levels = unique(Classes)))]
    } else if (input$rank == "Order" | input$rank == "order") {
      dt[, id := as.numeric(factor(Orders, levels = unique(Orders)))]
    } else if (input$rank == "Family" | input$rank == "family") {
      dt[, id := as.numeric(factor(Families, levels = unique(Families)))]
    } else if (input$rank == "Genus" | input$rank == "genus") {
      dt[, id := as.numeric(factor(Genera, levels = unique(Genera)))]
    } else if (input$rank == "Species" | input$rank == "species") {
      dt[, id := as.numeric(factor(AphiaIDs, levels = unique(AphiaIDs)))]
    }

    dt.out <- dt[J(unique(year)), mult = "last"]#[, Phylum := NULL]
    dt.out[, id := cummax(id)]
    numtaxa <- cummax(as.numeric(factor(dt$id)))
    taxa_dt <- aggregate(numtaxa, list(year = dt$year), max)
    colnames(taxa_dt) <- c("year", "taxa count")
    N_obs <- taxa_dt$'taxa count'
    times <- as.numeric(taxa_dt$year)

    SS<-getInitial(N_obs~SSlogis(times,alpha,xmid,scale),data=data.frame(N_obs=N_obs,times=times))
    K_start <- SS["alpha"]
    R_start <- 1/SS["scale"]
    N0_start <- SS["alpha"]/(exp(SS["xmid"]/SS["scale"])) + 1

    log_formula<-formula(N_obs ~ K * N0 * exp(R * times) / (K + N0 * (exp(R * times) - 1)))
    m<-nls(log_formula,start = list(K = K_start, R = R_start, N0 = N0_start))

    #corr_coef <- cor(N_obs,predict(m))
    lines(times,predict(m),col="red",lty=2,lwd=2)
    n = length(times)

    ## add model predictions
    K = summary(m)$coefficient[1]
    R = summary(m)$coefficient[2]
    N0 = summary(m)$coefficient[3]

    ## add variances - first, find standard errors
    K_se = summary(m)$coefficients[4]
    R_se = summary(m)$coefficients[5]
    N0_se = summary(m)$coefficients[6]

    ## compute standard deviations
    K_sd = K_se * sqrt(n)
    R_sd = R_se * sqrt(n)
    N0_sd = N0_se * sqrt(n)

    # compute upper bounds of model prediction
    UP = (K + K_sd) * (N0 + N0_sd) * exp((R + R_sd)*times)/((K + K_sd)+(N0 + N0_sd)*(exp((R + R_sd)*times)-1))
    lines(times, UP, col = 'red', lty = "dashed")
    LW = (K - K_sd) * (N0 - N0_sd) * exp((R - R_sd)*times)/((K - K_sd)+(N0 - N0_sd)*(exp((R - R_sd)*times)-1))
    lines(times, LW, col ='red', lty = 'dashed')
    taxa_dt <- taxa_dt
  })

  ranklable <- reactive({
      if (input$rank == "Phylum") {
      paste("phyla")
    } else if (input$rank == "Class") {
      paste("classes")
    } else if (input$rank == "Order") {
      paste("orders")
    } else if (input$rank == "Family") {
      paste("families")
    } else if (input$rank == "Genus") {
      paste("genera")
    } else if (input$rank == "Species") {
      paste("species")
    }
})


  output$dataview <- renderTable({
    head(viewData(), n = 6)
  }, caption = "Brief Data View", caption.placement = getOption("xtable.caption.placement", "top"), cex = 5)

  output$taxacurve <- renderPlot({

    if (input$model == FALSE) {
      df <- subset(data_m, Kingdoms == input$taxa | Phyla == input$taxa | Classes == input$taxa | Orders == input$taxa
                   | Families == input$taxa | Genera == input$taxa)
      dt = as.data.table(unique(df))
      setkey(dt, "year")
      if (input$rank == "Phylum" | input$rank == "phylum") {
        dt[, id := as.numeric(factor(Phyla, levels = unique(Phyla)))]
        ranklabel = "phyla"
      } else if (input$rank == "Class" | input$rank == "class") {
        dt[, id := as.numeric(factor(Classes, levels = unique(Classes)))]
        ranklabel = "class"
      } else if (input$rank == "Order" | input$rank == "order") {
        dt[, id := as.numeric(factor(Orders, levels = unique(Orders)))]
        ranklabel = "order"
      } else if (input$rank == "Family" | input$rank == "family") {
        dt[, id := as.numeric(factor(Families, levels = unique(Families)))]
        ranklabel = "family"
      } else if (input$rank == "Genus" | input$rank == "genus") {
        dt[, id := as.numeric(factor(Genera, levels = unique(Genera)))]
        ranklabel = "genus"
      } else if (input$rank == "Species" | input$rank == "species") {
        dt[, id := as.numeric(factor(AphiaIDs, levels = unique(AphiaIDs)))]
        ranklabel = "species"
      }

      dt.out <- dt[J(unique(year)), mult = "last"]#[, Phylum := NULL]
      dt.out[, id := cummax(id)]
      numtaxa <- cummax(as.numeric(factor(dt$id)))
      taxa_dt <- aggregate(numtaxa, list(year = dt$year), max)
      colnames(taxa_dt) <- c("year", "taxacount")
      taxa_dt <- taxa_dt
      #plot(taxa_dt, xlab = "Year", ylab = paste("Number of", tolower(input$rank), sep = " "), main = input$taxa, ylim = c(0, max(taxa_dt$"taxa count")*1.35))
      minx <- min(as.vector(taxa_dt$year))
      maxx <- max(as.vector(taxa_dt$year))
      ylab = paste("Number of", ranklabel, sep = " ")
      p <- ggplot(taxa_dt, aes(x = year, y = taxacount, colour = "#FF9999", group = 1)) + geom_point()
      p <- p + labs(x = "Year", y = ylab) + ggtitle(input$taxa) + scale_x_discrete(breaks = c(seq(minx, maxx, 25))) + theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), axis.text.y = element_text(angle = 60, hjust = 1), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
      p


      #plot(taxaaccum(), xlab = "Year", ylab = paste("Number of", tolower(input$rank), sep = " "), main = input$taxa, ylim = c(0, max(taxaaccum()$"taxa count")*1.35))
    } else if(input$model == TRUE) {
      df <- subset(data_m, Kingdoms == input$taxa | Phyla == input$taxa | Classes == input$taxa | Orders == input$taxa
                   | Families == input$taxa | Genera == input$taxa)
      dt = as.data.table(unique(df))
      setkey(dt, "year")
      if (input$rank == "Phylum" | input$rank == "phylum") {
        dt[, id := as.numeric(factor(Phyla, levels = unique(Phyla)))]
        ranklabel = "phyla"

      } else if (input$rank == "Class" | input$rank == "class") {
        dt[, id := as.numeric(factor(Classes, levels = unique(Classes)))]
        ranklabel = "class"
      } else if (input$rank == "Order" | input$rank == "order") {
        dt[, id := as.numeric(factor(Orders, levels = unique(Orders)))]
        ranklabel = "order"
      } else if (input$rank == "Family" | input$rank == "family") {
        dt[, id := as.numeric(factor(Families, levels = unique(Families)))]
        ranklabel = "family"
      } else if (input$rank == "Genus" | input$rank == "genus") {
        dt[, id := as.numeric(factor(Genera, levels = unique(Genera)))]
        ranklabel = "genus"
      } else if (input$rank == "Species" | input$rank == "species") {
        dt[, id := as.numeric(factor(AphiaIDs, levels = unique(AphiaIDs)))]
        ranklabel = "species"
      }

      dt.out <- dt[J(unique(year)), mult = "last"]#[, Phylum := NULL]
      dt.out[, id := cummax(id)]
      numtaxa <- cummax(as.numeric(factor(dt$id)))
      taxa_dt <- aggregate(numtaxa, list(year = dt$year), max)
      colnames(taxa_dt) <- c("year", "taxacount")
      #taxa_dt <- taxa_dt
      #plot(taxa_dt, xlab = "Year", ylab = paste("Number of", tolower(input$rank), sep = " "), main = input$taxa, ylim = c(0, max(taxa_dt$"taxa count")*1.35))
      #
      minx <- min(as.vector(taxa_dt$year))
      maxx <- max(as.vector(taxa_dt$year))
      ylab = paste("Number of", ranklabel, sep = " ")

      N_obs <- taxa_dt$'taxacount'
      times <- as.numeric(taxa_dt$year)

      SS<-stats::getInitial(N_obs~SSlogis(times,alpha,xmid,scale),data=data.frame(N_obs=N_obs,times=times))
      K_start <- SS["alpha"]
      R_start <- 1/SS["scale"]
      N0_start <- SS["alpha"]/(exp(SS["xmid"]/SS["scale"])) + 1

      log_formula<-formula(N_obs ~ K * N0 * exp(R * times) / (K + N0 * (exp(R * times) - 1)))
      m<-nls(log_formula,start = list(K = K_start, R = R_start, N0 = N0_start))

      corr_coef <- cor(N_obs,predict(m))
      lines(times,predict(m),col="red",lty=2,lwd=2)
      n = length(times)

      ## add model predictions
      K = summary(m)$coefficient[1]
      R = summary(m)$coefficient[2]
      N0 = summary(m)$coefficient[3]

      ## add variances - first, find standard errors
      K_se = summary(m)$coefficients[4]
      R_se = summary(m)$coefficients[5]
      N0_se = summary(m)$coefficients[6]

      ## compute standard deviations
      K_sd = K_se * sqrt(n)
      R_sd = R_se * sqrt(n)
      N0_sd = N0_se * sqrt(n)

      # compute upper bounds of model prediction
      UP = (K + K_sd) * (N0 + N0_sd) * exp((R + R_sd)*times)/((K + K_sd)+(N0 + N0_sd)*(exp((R + R_sd)*times)-1))
      #lines(times, UP, col = 'red', lty = "dashed")
      LW = (K - K_sd) * (N0 - N0_sd) * exp((R - R_sd)*times)/((K - K_sd)+(N0 - N0_sd)*(exp((R - R_sd)*times)-1))
      #lines(times, LW, col ='red', lty = 'dashed')
      #
      p <- ggplot(taxa_dt, aes(x = year, y = taxacount, colour = "#FF9999", group = 1)) + geom_point()
      p <- p + labs(x = "Year", y = ylab) + ggtitle(input$taxa) + scale_x_discrete(breaks = c(seq(minx, maxx, 25))) + theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), axis.text.y = element_text(angle = 60, hjust = 1), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
      p <- p + geom_line()
      p <- p + geom_ribbon(aes(ymin = LW, ymax = UP), linetype = 2, alpha = 0.1)
      #p <- plotly::ggplotly(p)
      p
    }
  })

})




