library(shiny)
library(ggplot2)
library(dplyr)
library(rootSolve)
library(shinythemes)

# Full Green-Ampt Model

# Newton method to find the root of the Implicit version of the Green & Ampt model
newtnp = function(tstart, tend, tp, tpp, rfi, amtinf, soil.para, tol = 0.00001, ...){
  
  # Soil parameters
  Sav = soil.para$Sav
  Ks  = soil.para$Ksat
  M = soil.para$ths - soil.para$thi
  
  tnp = tend
  
  bftry = (tnp-tstart)*rfi + amtinf
  
  arglog = 1.0 + bftry/(M*Sav)
  hnp = bftry - M*Sav*log(arglog) - Ks*(tnp - tp + tpp)
  dhnp= rfi - rfi*(1.0/arglog) - Ks;
  tnp = tnp - hnp/dhnp
  
  error = abs(hnp)
  itermax = 300
  iter = 0
  
  while(error > tol & iter < itermax){
    
    iter = iter + 1
    dhnp = rfi - rfi*(1/arglog) - Ks
    tnp = tnp - hnp/dhnp
    
    if(tnp < 0){
      tnp = 0
    }
    
    bftry = (tnp - tstart)*rfi + amtinf
    arglog = 1 + bftry/(M*Sav)
    hnp = bftry - M*Sav*log(arglog) - Ks*(tnp - tp + tpp)
    error = abs(hnp)
    
    # Print for convergence testing
    #print(c(iter, error, tnp))
    
    if(is.na(hnp)) {
      break;
      
    } else if(error < tol |iter == itermax) {
      break;
      
    }
    
  }
  
  tnpnew = tnp
  return(tnpnew)
}

# Chu's iterative solution
sschu<- function (Fpi, Time, tp, tpp, soil.para, tol = 0.00001, ...){
  
  M = soil.para$ths - soil.para$thi
  Sav = soil.para$Sav
  Ks = soil.para$Ksat
  
  Fpt = Ks*Time + M*Sav*log(1+ Fpi/(M*Sav))
  iter<- 0
  itermax <- 300
  error = abs(Fpt  -Fpi)
  
  while(error > tol & (iter< itermax)) {
    Fpi = Fpt
    Fpt = Ks*Time + M*Sav*log(1+ Fpi/(M*Sav))
    iter<- iter+1
    error = abs(Fpi- Fpt)
    # print(c(iter, error, Fpt))
  }
  if(error > tol & (iter == itermax)){
    stop("maximum iteration reached!")
  }
  return(Fpt)
}

# ------------------------------------------------------------------------------
# Start of the Green & Ampt infiltration model
# ------------------------------------------------------------------------------

# Inputs
# Time - Numeric vector indicating the time rainfall longs
# rainfall - Numeric vector indicating the volume of rainfall over time (hr)
# soil.para - Soil parameters: ths = saturated water content, thi = initial water content
#   Sav = soil suction (mm), and Ks = Saturated hydraulic conductivity (mm/hr)
# smax - Maximum soil depression storage (mm)
# tstep - Simulation's time step (hr)

green_ampt_full <- function(Time, rainfall, soil.para, smax = 0.5, tstep, tol = 0.1e-5, ...){
  
  ths = soil.para$ths
  thi = soil.para$thi
  Sav = soil.para$Sav
  Ks  = soil.para$Ksat
  
  M = ths - thi
  
  # delt = Time[2] - Time[1]
  
  # Time to poinding (tp) and new, compressed time, Tnew
  
  # Tnew = Time - tp + t0
  
  Inf.capacity<- function (Fp){
    A = Ks
    B = Ks*Sav*M
    fp<- A + B/Fp
    return(fp)
  }
  
  time.new = function(Time, tp, t0) (Time - tp + t0)
  
  fdata<- data.frame(Time = Time, rainfall = rainfall,
                     frate = rep(0, length(Time)),
                     fp = rep(0, length(Time)),
                     Fcum = rep(0, length(Time)),
                     stor = rep(0, length(Time)),
                     Fp = rep(0, length(Time)),
                     ro = rep(0, length(Time)))
  
  tpf = rep(0, length(Time))
  tppf = rep(0, length(Time))
  
  delt = Time[2] - Time[1]
  Fp = M*Sav/(rainfall[1]/Ks - 1)
  tp = Fp/rainfall[1]
  Fp_old = 0
  fdata$Fp[1] = max(Fp,0)
  fp_old = Inf.capacity(Fp)
  npp = "npp"
  
  if(tp < 0 | tp > Time[1]) { 
    
    # Means R < Ks No ponding durign this period
    Fcum = Time[1]*rainfall[1]
    frate = min(rainfall[1], Inf.capacity(Fcum))
    rexcess = 0
    dF = 0
    fdata$stor[1] = 0
    fdata$ro[1] = 0
    fdata$frate[1] = frate
    fdata$Fcum[1] = Fcum
    
    endPond = Time[1]
    new_rain = rainfall[2]
    endPond = Time[1]
    tpf[1] = npp
    tppf[1] = npp
    
  } else if(tp < Time[1] & tp > 0){
    
    # tp = Ks*M*Sav/((rainfall[1]-Ks))
    t0 = Fp/Ks - M*(Sav/Ks)*log(1+ Fp/(M*Sav))
    
    Tnew = time.new(Time = Time[1], tp = tp, t0 = t0)
    Fcum = sschu(Fpi = Fp, Time = Tnew, soil.para, tol = tol)
    fp = Inf.capacity(Fp = Fcum)
    delt_new = Time[1] - tp
    
    frate = min(fp, new_rain)
    
    dP = delt_new*rainfall[1]
    dF = Fcum - Fp
    rexcess = dP - dF
    new_rain = rainfall[2]
    
    # dS = rexcess - smax
    fdata$frate[1] = frate
    fdata$Fcum[1] = Fcum
    fdata$stor[1] = min(smax, rexcess)
    fdata$ro[1] = max(0, rexcess - smax)
    fdata$Fp[1] = Fp
    fdata$fp[1] = fp
    
    # fp[1] = Inf.capacity(Fp)
    tppf[1] = t0
    tpf[1] = tp
    
    
  }
  
  
  for(i in 2: (nrow(fdata))){
    # for(i in 2:7){
    
    tstart = Time[i-1]
    tend = Time[i]
    
    tt = Time[i]
    
    frate_old = fdata$frate[i-1]
    Fcum_old = fdata$Fcum[i-1]
    stor_old = fdata$stor[i-1]
    scap_old = max(0, smax - stor_old)
    delt = tt - Time[i-1]
    new_rain = rainfall[i]
    Fp_old = fdata$Fp[i-1]
    fp_old = fdata$fp[i-1]
    
    
    if(stor_old == 0) {
      
      if(fp_old > 0 & new_rain > fp_old){
        tps = 0
        tp = tps + tstart
        Fpt = M*Sav*Ks/(fp_old - Ks)
        fp_new = Inf.capacity(Fpt)
      } else if(new_rain > Ks){
        Fpt = M*Sav/(new_rain/Ks - 1)
        tps = Fpt/new_rain
        tp = tps + tstart
        fp_new = Inf.capacity(Fpt)
      }else {
        Fpt = 0
        tp = 9999
        tps = 9999
        # fp_new = Inf.capacity(Fcum_old)
      }
      
      t0 = (Fpt - M*Sav*log(1.0 + Fpt/(M*Sav)))/Ks
      t0 = ifelse(is.na(t0), 0, t0)
      
      # if(new_rain < fp_old & fp_old > 0){
      #          t0 = tstart + t0
      # }
      
      tppf[i] = ifelse(t0 == 0, npp, t0)
      tpf[i] = ifelse(tp == 9999, npp, tp)
      
      
      if(tp > tend){
        
        Fcum_new = Fcum_old + delt*new_rain
        fp_new = Inf.capacity(Fcum_new)
        frate_new = min(fp_new, new_rain)
        rexcess_new = 0
        dF_new = 0
        fdata$stor[i] = 0
        fdata$ro[i] = 0
        fdata$frate[i] = frate_new
        fdata$Fcum[i] = Fcum_new
        fdata$fp[i] = fp_new
        endPond = Time[i]
        fdata$fp[i] = fp_new
        tpf[i] = npp
        
      } else {
        
        if(new_rain > fp_old) {
          endPond = tend + delt
        } else {
          
          new_rain2 = (tend - tp)*new_rain/delt
          
          amt_inf = Fcum_old + smax
          
          endPond = newtnp(tstart = tstart, tend = tend,
                           tp = tp, tpp = t0, rfi = new_rain2,
                           amtinf = amt_inf, soil.para = soil.para)
          
        }
        
        if(endPond >= tend | new_rain > fp_old){
          
          Tnew = time.new(Time = tend, tp = (tp), t0 = t0)
          Fcum_old = ifelse(i == 2, Fpt, Fcum_old)
          Fcum_new = sschu(Fpi = Fcum_old, Time = Tnew,
                           soil.para = soil.para, tol = tol)
          
          fp_new = Inf.capacity(Fcum_new)
          frate_new = min(fp_new, new_rain)
          
          dF_new = Fcum_new - Fcum_old
          
          deltp = tend - tp
          dP_new = deltp*new_rain
          
          rexcess_new = dP_new - dF_new
          
          stor_new = min(smax, stor_old + rexcess_new)
          stor_new = max(0, stor_new)
          
          ro_new = ifelse(rexcess_new <= scap_old, 0,
                          (rexcess_new - scap_old))
          
          fdata$Fcum[i] = Fcum_new
          fdata$frate[i] = frate_new
          fdata$stor[i] = stor_new
          fdata$ro[i] = ro_new
          fdata$Fp[i] = Fp_old
          fdata$fp[i] = fp_new
          tpf[i] = tp
          tppf[i] = t0
          
        } else {
          
          time_ndnp = tend - endPond
          Tnew = time.new(Time = endPond, tp = tp, t0 = t0)
          
          Fcum_new_pnd = sschu(Fpi = Fcum_old, Time = Tnew,
                               soil.para = soil.para, tol = tol)
          frate_new_pnd = Inf.capacity(Fcum_new_pnd)
          
          # nonponded portion
          Fcum_new_npnd = (tend-endPond)*new_rain
          Fcum_new = Fcum_new_pnd + Fcum_new_npnd
          
          frate_new = Inf.capacity(Fcum_new)
          
          dF_new = Fcum_new - Fcum_old
          dP_new = delt*new_rain
          # rexcess_new = dP_new - Fcum_new
          
          rexcess_new = dP_new - dF_new
          stor_new = min(smax, stor_old + rexcess_new)
          stor_new = max(0, stor_new)
          
          ro_new = ifelse(rexcess_new <= scap_old, 0,
                          (rexcess_new - scap_old))
          
          fdata$Fcum[i] = Fcum_new
          fdata$frate[i] = frate_new
          fdata$stor[i] = stor_new
          fdata$ro[i] = ro_new
          fdata$Fp[i] = Fp_old
          fdata$fp[i] = frate_new
          tpf[i] = tp
          tppf[i] = t0
        }
        
        
      }
      
    } else if(stor_old > 0){
      
      itp = as.numeric(tpf[i-1])
      itpp = as.numeric(tppf[i-1])
      
      if(new_rain > fp_old) {
        endPond = tend + delt
        
      } else {
        
        amt_inf = Fcum_old + stor_old
        
        endPond = newtnp(tstart = tstart, tend = tend,
                         tp = itp, tpp = itpp, rfi = new_rain,
                         amtinf = amt_inf, soil.para = soil.para)
        
      }
      
      if(endPond >= tend | new_rain > fp_old){
        
        Tnew = time.new(Time = tend, tp = itp, t0 = itpp)
        
        Fcum_new = sschu(Fpi = Fcum_old, Time = Tnew,
                         soil.para = soil.para, tol = tol)
        
        fp_new = Inf.capacity(Fcum_new)
        frate_new = min(fp_new, new_rain)
        dF_new = Fcum_new - Fcum_old
        
        
        dP_new = delt*new_rain
        rexcess_new = dP_new - dF_new
        
        stor_new = min(smax, stor_old + rexcess_new)
        stor_new = max(0, stor_new)
        
        ro_new = ifelse(rexcess_new <= scap_old, 0,
                        (rexcess_new - scap_old))
        
        fdata$Fcum[i] = Fcum_new
        fdata$frate[i] = frate_new
        fdata$stor[i] = stor_new
        fdata$ro[i] = ro_new
        fdata$Fp[i] = Fp_old
        fdata$fp[i] = fp_new
        tpf[i] = itp
        tppf[i] = itpp
        
      } else {
        
        time_ndnp = tend - endPond
        Tnew = time.new(Time = endPond, tp = itp, t0 = itpp)
        
        Fcum_new_pnd = sschu(Fpi = Fcum_old, Time = Tnew,
                             soil.para = soil.para, tol = tol)
        
        fp_new_pnd = Inf.capacity(Fcum_new_pnd)
        
        frate_new_pnd = min(fp_new_pnd, new_rain)
        
        # nonponded portion
        Fcum_new_npnd = (tend-endPond)*rainfall[i]
        
        frate_new_npnd = min(Inf.capacity(Fcum_new_npnd), new_rain)
        frate_new = mean(frate_new_pnd, frate_new_npnd)
        
        Fcum_new = Fcum_new_pnd + Fcum_new_npnd
        
        dF_new = Fcum_new - Fcum_old
        dP_new = delt*new_rain
        
        
        rexcess_new = dP_new - dF_new
        
        stor_new = min(smax, stor_old + rexcess_new)
        stor_new = max(0, stor_new)
        
        ro_new = ifelse(rexcess_new <= scap_old, 0,
                        (rexcess_new - scap_old))
        
        fdata$Fcum[i] = Fcum_new
        fdata$frate[i] = frate_new
        fdata$stor[i] = stor_new
        fdata$ro[i] = ro_new
        fdata$Fp[i] = Fp_old
        fdata$fp[i] = fp_new_pnd
        tppf[i] = t0
        tpf[i] = tp
      }
      
      
      
    }
    
  }
  
  # Need to infiltrate remaining surface storage
  
  rem_stor = fdata$stor[nrow(fdata)]
  
  if(rem_stor > 0){
    
    stor_old = rem_stor
    tstart = fdata$Time[nrow(fdata)]
    tend = 1000
    Fcum_old = fdata$Fcum[nrow(fdata)]
    
    amt_inf = Fcum_old + stor_old
    new_rain = 0.0
    endPond = newtnp(tstart = tstart, tend = tend,
                     tp = itp, tpp = itpp, rfi = new_rain,
                     amtinf = amt_inf, soil.para = soil.para)
    
    
    endPond_i = endPond %% tstep
    endPond = ifelse(endPond_i > 0, endPond + tstep - endPond_i, endPond)
    
    add_times = seq(tstart + tstep, (endPond), by = tstep)
    
    aa = nrow(fdata) + 1
    bb = nrow(fdata) + length(add_times)
    
    fdata_new = as.data.frame(matrix(NA, length(add_times), ncol(fdata)))
    names(fdata_new) = names(fdata)
    fdata_new$Time = add_times
    fdata_new$rainfall = 0
    fdata_new$ro = 0
    
    nno = length(Time)
    
    ii = 1
    
    last_time = fdata$Time[nrow(fdata)]
    tstart_new = last_time
    ai = nno + ii
    
    while(stor_old > 0 & tstart_new < (endPond - tstep)){
      
      fdata[ai, ] = rep(NA, ncol(fdata))
      tstart_new = last_time + tstep*ii
      fdata$Time[ai] = tstart_new
      
      Tnew = time.new(Time = tstart_new, tp = itp , t0 = 0)
      
      Fcum_new = sschu(Fpi = Fcum_old, Time = Tnew,
                       soil.para = soil.para, tol = tol)
      
      
      dF_new = Fcum_new - Fcum_old
      
      fp_new = Inf.capacity(Fcum_new)
      dP_new = 0
      rexcess_new = dP_new - dF_new
      
      stor_new = min(smax, stor_old + rexcess_new)
      stor_new = max(0, stor_new)
      
      ro_new = 0
      fdata$rainfall[ai] = 0
      fdata$Fcum[ai] = Fcum_new
      fdata$frate[ai] = frate_new
      fdata$stor[ai] = stor_new
      fdata$ro[ai] = ro_new
      fdata$Fp[ai] = Fp_old
      fdata$fp[ai] = frate_new
      tpf[ai] = itp
      tppf[ai] = itpp
      stor_old = stor_new
      Fcum_old = Fcum_new
      
      ii = ii + 1
      ai = nno + ii
      
    }
    
    fdata[ai, ] =    rep(0, ncol(fdata))
    fdata$Time[ai] = max(fdata$Time) + tstep
    fdata[ai, "Fcum"] = Fcum_new + stor_new
    fdata[ai, "frate"] = Inf.capacity(Fcum_new + stor_new)
    
  }
  
  fdata$cum_ro = cumsum(fdata$ro)
  fdata$cum_rain = cumsum(fdata$rainfall*tstep)
  
  tpf = rep(tpf[length(tpf)], 100)
  tppf = rep(tppf[length(tppf)], 100)
  tpf = tpf[1:nrow(fdata)]
  tppf = tppf[1:nrow(fdata)]
  
  
  options(warn = -1)
  fdata = apply(X = fdata, MARGIN = 2, FUN = function(x) round(x, digits = 3))
  fdata = as.data.frame(fdata)
  fdata$tp = ifelse (tpf == "npp", "npp", round(as.numeric(tpf), digits = 3))
  fdata$tpp = ifelse (tppf == "npp", "npp", round(as.numeric(tppf), digits = 3))
  options(warn = 0)
  
  tot_rain = max(fdata$cum_rain)
  tot_ro = max(fdata$cum_ro)
  tot_inf = max(fdata$Fcum)
  peak_frate = max(fdata$frate)
  
  balance  = tot_rain - (tot_ro + tot_inf)
  
  wat_bal = c(total_rainfall = tot_rain,
              total_inf = tot_inf,
              tot_runoff = tot_ro,
              peak_inf_rate = peak_frate,
              inf_plot_ro = tot_ro + tot_inf,
              balance = tot_rain - tot_inf - tot_ro)
  
  wat_bal = round(wat_bal, digits = 3)
  
  return(list(infil_data = fdata, water_balance = wat_bal))
}

# ------------------------------------------------------------------------------
# End of the Green & Ampt infiltration model
# ------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# This code runs the graphical interface of the Green-Ampt method
# -------------------------------------------------------------------------------

ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("Green-Ampt Infiltration Model"),
  
  # Subtitle Panel (Corrected)
  tags$div(
    style = "background-color: #f0f0f0; padding: 10px; border-bottom: 1px solid #ccc; text-align: center;",
    p(
      "Developed by: Josept D. Revuelta-Acosta",
      br(),
      "Email: jrevuelta@uv.mx",
      br(),
      "University: Universidad Veracruzana - Course: Hydrology",
      style = "font-style: italic;"
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("rainfall_duration", "Rainfall Duration", value = 7, min = 1),
      textAreaInput("rainfall_intensity", "Rainfall Intensity (cm/time step) - comma separated values", value = ("0, 3.0, 0.1, 0.1, 1.0, 0.4, 0.6"), rows = 5),
      numericInput("dt", "Time Step (hours)", value = 1, min = 0.1),
      numericInput("Ks", "Saturated Hydraulic Conductivity (cm/hr)", value = 0.044),
      numericInput("psi_f", "Wetting Front Suction Head (cm)", value = 22.4),
      sliderInput("theta_i", "Initial Water Content (fraction)", min = 0, max = 1, value = 0.25, step = 0.01),
      sliderInput("theta_s", "Saturation Water Content (fraction)", min = 0, max = 1, value = 0.499, step = 0.01),
      actionButton("run_model", "Run Model")
    ),
    
    mainPanel(
      fluidRow(
        column(6,
               plotOutput("rainfall_plot", height = "300px"),
               plotOutput("infiltration_plot", height = "300px")
        ),
        column(6,
               plotOutput("capacity_plot", height = "300px"),
               plotOutput("cumulative_plot", height = "300px")
        )
      ),
      verbatimTextOutput("results_table")
    )
  )
)

server <- function(input, output) {
  observeEvent(input$run_model, {
    # Rainfall input handling (Corrected)
    rainfall_input_text <- input$rainfall_intensity
    rainfall_input <- as.numeric(trimws(unlist(strsplit(rainfall_input_text, ","))))
    
    if (length(rainfall_input)!= input$rainfall_duration) {
      showNotification("Rainfall data points must match duration", type = "error")
      return()
    }
    
    Time <- seq(from = 0, to = length(rainfall_input)-1, by = 1)*input$dt
    soil_para = list(ths = input$theta_s, thi = input$theta_i, Ksat = input$Ks, Sav = input$psi_f)
    
    # Overall results
    results.over <- green_ampt_full(Time = Time, rainfall = rainfall_input, soil.para = soil_para, tstep = input$dt)
    
    # Result for plots
    results <- data.frame(time = results.over$infil_data$Time, rainfall = results.over$infil_data$rainfall,
                                 infiltration = results.over$infil_data$frate,
                                 runoff = results.over$infil_data$ro, cumulative_infiltration = results.over$infil_data$Fcum,
                                 infiltration_capacity = results.over$infil_data$fp)
    
    create_plot <- function(data, x_var, y_var, title, color, y_label) {
      ggplot(data, aes(x = time, y =.data[[y_var]])) +
        geom_line(color = color, linewidth = 1.2) +
        labs(title = title, x = "Time (hours)", y = y_label) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10)
        )
    }
    
    output$rainfall_plot <- renderPlot({
      ggplot(results, aes(x = time, y = rainfall)) +
        geom_col(fill = "skyblue", color = "black") +
        labs(title = "Rainfall", x = "Time (hours)", y = "Rainfall (mm/time step)") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    })
    
    output$infiltration_plot <- renderPlot({
      ggplot(results, aes(x = time)) +
        geom_line(aes(y = infiltration, color = "Infiltration"), linewidth = 1.2) +
        geom_line(aes(y = runoff, color = "Runoff"), linewidth = 1.2) +
        geom_col(aes(y = rainfall, fill = "Rainfall"), width = 0.2, alpha = 0.5) +
        labs(title = "Infiltration, Runoff, and Rainfall", x = "Time (hours)", y = "Amount (mm/time step)") +
        scale_color_manual(values = c("Infiltration" = "green", "Runoff" = "red")) +
        scale_fill_manual(values = c("Rainfall" = "skyblue")) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = c(0.8, 0.8),
          legend.title = element_blank()
        )
    })
    
    output$capacity_plot <- renderPlot({
      create_plot(results, "time", "infiltration_capacity", "Infiltration Capacity", "darkgreen", "Infiltration Capacity (mm/hr)")
    })
    
    output$cumulative_plot <- renderPlot({
      create_plot(results, "time", "cumulative_infiltration", "Cumulative Infiltration", "purple", "Cumulative Infiltration (mm)")
    })
    
    output$results_table <- renderPrint({
      print(results.over)
    })
  })
}

shinyApp(ui = ui, server = server)