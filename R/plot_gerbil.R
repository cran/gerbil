plot_gerbil <- function(gerb, type = "Univariate", vars = NULL, file = NULL, sep = FALSE, imp = NULL,
                        obs.col = "blue4", imp.col = "brown2", obs.lty = 1, imp.lty = 2, obs.lwd = NULL, 
                        imp.lwd = NULL,  height = NULL, width = NULL, legend.spot = "topright", legend.text = legend.text, 
                        mfrow = c(3, 2), trace.type = "Mean", log = NULL, pch = c(1, 3), partial = "imputed", ...)
{
  
  if (tolower(type) == "univariate" | type == 1) {
    plot_gerbil.uni(gerb = gerb, file = file, vars = vars, imp = imp, obs.col = obs.col, imp.col = imp.col, sep = sep,
                    obs.lty = obs.lty, imp.lty = imp.lty, use.log = log, legend.text = legend.text, legend.spot = legend.spot,
                    height = height, width = width, mfrow = mfrow, obs.lwd = obs.lwd, imp.lwd = imp.lwd, ...)
  } else if (tolower(type) == "bivariate" | tolower(type) == "multivariate" | type == 2) {
    #if (length(imp) > 1) {
    #  warning("imp can not have length greater than 1 for bivariate plotting. Will use only the first element.")
    #}
    plot_gerbil.mult(gerb = gerb, file = file, vars = vars, imp = imp, obs.col = obs.col, imp.col = imp.col, 
                     lty = c(obs.lty, imp.lty), lwd = c(obs.lwd, imp.lwd), sep = sep, legend.text = legend.text, 
                     legend.spot = legend.spot, height = height, width = width, mfrow = mfrow, pch = pch, use.log = log, partial = partial, ...)
  } else if (tolower(type) == "ts" | tolower(type) == "time.series" | type == 3) {
    plot_gerbil.ts(gerb = gerb, file = file, vars = vars, imps = imp, obs.col = obs.col, imp.col = imp.col, sep = sep, legend.text = legend.text, 
                   legend.spot = legend.spot, height = height, width = width, mfrow = mfrow, type = trace.type, 
                   obs.lty = obs.lty, imp.lty = imp.lty, obs.lwd = obs.lwd, imp.lwd = imp.lwd, ...)
  } else {
    stop("Input 'type' not recognized.")
  }
  
}

plot_gerbil.ts <- function (gerb, file = NULL, vars = NULL, obs.col = "blue4", imp.col = "brown2", obs.lty = 2, imp.lty = 1,
                            imps = NULL, sep = FALSE, height = NULL, width = NULL,
                            type = "Mean", plot.obs = TRUE, obs.lwd = NULL, imp.lwd = NULL, mfrow = c(3, 2), legend.text = NULL, 
                            legend.spot = "topright", main = NULL, xlab = NULL, ylab = NULL, bty = NULL, ...) {
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  plot.legend <- TRUE
  if (is.logical(legend.text)) {
    if(!legend) {
      plot.legend <- FALSE
    } else {
      legend.text <- NULL
    }
  } else if (length(legend.text) == 1) {
    if (legend.text == "n") {
      plot.legend <- FALSE
    } else if (legend.text == "y") {
      legend.text <- NULL
    }
  }

  if (length(gerb$chainSeq) == 0) {
    stop("The gerbil object was created with trace = FALSE. Cannot create trace plots.")
  }

  if (type == 1 | tolower(type) == "mean") {
    gerb.tmp <- gerb$chainSeq$means.mis
    gerb.tmp.obs <- gerb$chainSeq$means.obs
  } else if (type == 2 | tolower(type) == "variance") {
    if (is.element("vars.mis", names(gerb$chainSeq))) {
      gerb.tmp <- gerb$chainSeq$vars.mis
      gerb.tmp.obs <- gerb$chainSeq$vars.obs
    } else {
      stop("Input 'gerb' was produced with 'calc.var = FALSE'.")
    }
  } else if (type == 3 | tolower(type) == "latent.mean") {
    gerb.tmp <- gerb$chainSeq$means.mis.latent
    gerb.tmp.obs <- NULL
    if (plot.obs) {
      warning("Will not plot means of the latent process for observed values.")
      plot.obs <- FALSE
    }
  } else if (type == 4 | tolower(type) == "latent.variance") {
    if (is.element("vars.mis.latent", names(gerb$chainSeq))) {
      gerb.tmp <- gerb$chainSeq$vars.mis.latent
      gerb.tmp.obs <- NULL
      if (plot.obs) {
        warning("Will not plot variances of the latent process for observed values.")
        plot.obs <- FALSE
      }
    } else {
      stop("Input 'gerb' was produced with 'calc.var = FALSE'.")
    }
  } else {
    stop ("Invalid value of input 'type'.")
  }
  
  if (length(vars) == 0) {
    #vars <- dimnames(gerb.tmp)[[1]]
    vars <- colnames(gerb$missing.latent)[colSums(gerb$missing.latent == 1 | gerb$missing.latent == 2) > 0]
  }
  
  if (length(file) == 0 | length(vars) == 1) {
    sep <- TRUE
    supp.suff <- TRUE
  } else {
    supp.suff <- FALSE
  }
  
  if (sep) {
    mfrow <- c(1, 1)
  }
  
  if(length(height) == 0) {
    if (!sep) {
      height <- 11
    } else {
      height <- 6.5
    }
  }
  
  if(length(width) == 0) {
    if (!sep) {
      width <- 8.5
    } else {
      width <- 6.5
    }
  }
  
  if (length(file) == 0) {
    graphics::par(mfrow = c(1, 1), ask = FALSE)
  } else {
    file.type <- "pdf"
    if (substr(file, nchar(file) - 3, nchar(file)) ==
        ".pdf") {
      file <- substr(file, 1, nchar(file) - 4)
    } else if (substr(file, nchar(file) - 3, nchar(file)) ==
               ".png") {
      file <- substr(file, 1, nchar(file) - 4)
      file.type <- "png"
    }
    if (!sep) {
      if (file.type == "pdf") {
        file <- paste(file, ".pdf", sep = "")
        grDevices::pdf(file = file, width = width, height = height)
      } else if (file.type == "png") {
        file <- paste(file, ".png", sep = "")
        grDevices::png(file = file, width = width, height = height,
                       units = "in", res = 500)
      }
      graphics::par(mfrow = mfrow, ask = FALSE)
    }
  }
  
  latent.vars <- dimnames(gerb.tmp)[[1]]
  orig.vars <- dimnames(gerb[[1]][[1]])[[2]]

  nams.out <- gerb$nams.out

  vars.new <- vars.old <- NULL

  for (i in 1:length(vars)) {
    if (is.element(vars[i], latent.vars)) {
      tmp <- vars[i]
      tmp1 <- nams.out[names(nams.out) == vars[i]]
    } else if (is.element(vars[i], orig.vars)) {
      tmp <- names(nams.out)[nams.out == vars[i]]
      tmp1 <- rep(vars[i], length(tmp))
    } else {
      message(paste0("Variable name ", vars[i], " is invalid and will not be used for plotting."), "\n", appendLF = FALSE)
      tmp <- tmp1 <- NULL
    }

    if (length(vars.new) == 0 & length(tmp) > 0) {
      vars.new <- tmp
      vars.old <- tmp1
    } else if (length(tmp) > 0) {
      vars.new[(length(vars.new) + 1):(length(vars.new) + length(tmp))] <- tmp
      vars.old[(length(vars.old) + 1):(length(vars.old) + length(tmp1))] <- tmp1
    }
  }

  if (length(vars.new) == 0) {
    stop("No valid variable names were given.")
  }

  vars <- vars.new
  
  if (is.logical(imps)) {
    if (imps) {
      imps <- as.numeric(dimnames(gerb.tmp)[[3]])
    } else {
      imps <- 1
    }
  }
  if(length(imps) == 0) {
    imps <- 1
  }
  
  if(length(obs.lwd) == 0) {
    if (!sep & length(file) != 0) {
      obs.lwd <- 1
    } else {
      obs.lwd <- 2
    }
  }

  if(length(imp.lwd) == 0) {
    if (!sep & length(file) != 0) {
      imp.lwd <- 1
    } else {
      imp.lwd <- 2
    }
  }
  
  if (length(imps) == 1) {
    ltys <- imp.lty
    lwds <- imp.lwd
  } else {
    if(length(imp.col) == 1) {
      imp.col <- rep(imp.col, length(imps))
      imp.col <- imp.col[1:length(imps)]
    }
    if(length(imp.lwd) == 1) {
      imp.lwd <- rep(imp.lwd, length(imps))
      imp.lwd <- imp.lwd[1:length(imps)]
    }
    if(length(imp.lty) == 1) {
      imp.lty <- c(imp.lty, 3:(length(imps) + 1))
    } else {
      imp.lty <- rep(imp.lty, length(imps))
      imp.lty <- imp.lty[1:length(imps)]
    }
    ltys <- imp.lty
    lwds <- imp.lwd
  }
  
  for (i in 1:length(vars)) {
    
    n.obs <- sum(gerb[[2]][,vars.old[i]] == 0)
    n.mis <- sum(gerb[[2]][,vars.old[i]] != 0)
    
    cont <- TRUE
    if (n.mis == 0) {
      cont <- FALSE
      #warning(paste0("Will not produce univariate plot for ", var, ": No imputed values for the variable."))
      message(paste0("Will not produce trace plots for ", vars[i], ": No imputed values for the variable."), "\n", appendLF = FALSE)
    } else {
      cont <- TRUE
    }
    
    if (cont) {
      
      if (plot.obs) {
        z <- gerb.tmp.obs[vars[i]]
      } else {
        z <- NULL
      }
      ylim <- c(min(c(gerb.tmp[vars[i], , imps], z)),max(c(gerb.tmp[vars[i], , imps], z)))
      
      y <- gerb.tmp[vars[i], , imps[1]]
      x <- as.numeric(dimnames(gerb.tmp)[[2]])
      
      if (sep & length(file) > 0) {
        if (!supp.suff) {
          suff <- paste0("_", vars[i], "_time_series")
        } else {
          suff <- ""
        }
        if (file.type == "pdf") {
          grDevices::pdf(file = paste(file, suff, ".pdf", sep = ""), width = width, height = height)
        } else if (file.type == "png") {
          grDevices::png(file = paste(file, suff, ".png", sep = ""),
                         width = width, height = height, units = "in", res = 500)
        }
      }
      
      if (length(main) == 0) {
        main1 <- vars[i]
      } else {
        main1 <- main 
      }
      if (length(ylab) == 0) {
        ylab1 <- vars[i]
      } else {
        ylab1 <- ylab 
      }
      if (length(xlab) == 0) {
        xlab1 <- "Iteration"
      } else {
        xlab1 <- xlab 
      }
      plot(x, y, type = "l", lty = ltys[1], col = imp.col[1], main = main1, ylab = ylab1, 
           xlab = xlab1, ylim = ylim, lwd = lwds[1], panel.first = grid(), ...)
      
      if (length(imps) > 1) {
        for(j in 2:length(imps)) {
          y <- gerb.tmp[vars[i], , imps[j]]
          lines(x, y, lty = ltys[j], col = imp.col[j], lwd = lwds[j], ...)
        }
      }
      
      leg <- paste0("Imputation ",imps)

      cols <- imp.col
      
      if (plot.obs) {
        #obs.col <- rgb(0.3,0.2,0.4,0.6)
        abline(h = z, col = obs.col, lty = obs.lty, lwd = obs.lwd)
        leg <- c(leg, "Observed")
        cols <- c(cols, obs.col)
        ltys <- c(ltys, obs.lty)
        lwds <- c(lwds, obs.lwd)
      }

      if (length(bty) == 0) {
        bty <- "n"
      }

      if (plot.legend) {
        if (length(legend.text) > 0) {
          leg[1:length(legend.text)] <- legend.text
        }
        legend(legend.spot, legend = leg, bty = bty, col = cols, lty = ltys, lwd = lwds, ...)
      }      

      if (sep & length(file) > 0) {
        grDevices::dev.off()
      }
    }
  }
  if (!sep & length(file) > 0) {
    grDevices::dev.off()
  }
  
}


plot_gerbil.uni <- function (gerb, file = NULL, vars = NULL, imp = 1, obs.col = "blue4", imp.col = "brown2", sep = FALSE, obs.lty = 1, imp.lty = 2,
                             use.log = NULL, legend.text = NULL, legend.spot = "bottomleft", height = NULL, width = NULL, mfrow = c(3, 2), obs.lwd = NULL, imp.lwd = NULL, 
                             main = NULL, xlab = NULL, ylab = NULL, bty = NULL, cex = NULL, ...) {

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  plot.legend <- TRUE
  if (is.logical(legend.text)) {
    if(!legend) {
      plot.legend <- FALSE
    } else {
      legend.text <- NULL
    }
  } else if (length(legend.text) == 1) {
    if (legend.text == "n") {
      plot.legend <- FALSE
    } else if (legend.text == "y") {
      legend.text <- NULL
    }
  }

  if(length(obs.lwd) == 0) {
    if (!sep & length(file) > 0) {
      obs.lwd <- 1
    } else {
      obs.lwd <- 2
    }
  }

  if(length(imp.lwd) == 0) {
    if (!sep & length(file) > 0) {
      imp.lwd <- 1
    } else {
      imp.lwd <- 2
    }
  }
  
  if (is.logical(imp)) {
    if (imp) {
      imp <- as.numeric(dimnames(gerb)[[3]])
    } else {
      imp <- 1
    }
  }
  if(length(imp) == 0) {
    imp <- 1
  }

  #imp <- imp[1]

  if(length(imp.col) == 1) {
    imp.col <- rep(imp.col, length(imp))
    imp.col <- imp.col[1:length(imp)]
  }
  if(length(imp.lwd) == 1) {
    imp.lwd <- rep(imp.lwd, length(imp))
    imp.lwd <- imp.lwd[1:length(imp)]
  }
  if(length(imp.lty) == 1) {
    imp.lty <- c(imp.lty, 3:(length(imp) + 1))
  } else {
    imp.lty <- rep(imp.lty, length(imp))
    imp.lty <- imp.lty[1:length(imp)]
  }

  col <- c(obs.col, imp.col)
  lty <- c(obs.lty, imp.lty)
  lwd <- c(obs.lwd, imp.lwd)
  
  if (length(vars) == 0) {
    #vars <- rownames(gerb[[3]])
    #vars <- rownames(gerb$summary)[gerb$summary[, 4] > 0]
    vars <- colnames(gerb$missing)[colSums(gerb$missing == 1 | gerb$missing == 2) > 0]
  } else if (is.numeric(vars)) {
    vars <- rownames(gerb[[3]])[vars]
  }
  
  if (length(file) == 0 | length(vars) == 1) {
    sep <- TRUE
    supp.suff <- TRUE
  } else {
    supp.suff <- FALSE
  }
  
  if (sep) {
    mfrow <- c(1, 1)
  }
  
  if(length(height) == 0) {
    if (!sep) {
      height <- 11
    } else {
      height <- 11
    }
  }
  
  if(length(width) == 0) {
    if (!sep) {
      width <- 8.5
    } else {
      width <- 8.5
    }
  }
  
  if (length(file) == 0) {
    graphics::par(mfrow = c(1, 1), ask = FALSE)
  } else {
    file.type <- "pdf"
    if (substr(file, nchar(file) - 3, nchar(file)) ==
        ".pdf") {
      file <- substr(file, 1, nchar(file) - 4)
    } else if (substr(file, nchar(file) - 3, nchar(file)) ==
               ".png") {
      file <- substr(file, 1, nchar(file) - 4)
      file.type <- "png"
    }
    if (!sep) {
      if (file.type == "pdf") {
        file <- paste(file, ".pdf", sep = "")
        grDevices::pdf(file = file, width = width, height = height)
      } else if (file.type == "png") {
        file <- paste(file, ".png", sep = "")
        grDevices::png(file = file, width = width, height = height,
                       units = "in", res = 500)
      }
      graphics::par(mfrow = mfrow, ask = FALSE)
    }
  }
  
  if (is.logical(use.log)) {
    if(use.log) {
      use.log <- vars
    }
  }
 
  for(i in 1:length(vars)) {
    
    #type <- as.character(gerb[[3]][c(i), "Variable.Type"])
    #var <- rownames(gerb[[3]])[c(i)]
    
    type <- as.character(gerb[[3]][vars[i], "Variable.Type"])
    var <- vars[i]
    
    if (length(imp) > 1) {
      dat <- matrix(NA, NROW(gerb[[1]][[imp[1]]]), length(imp))
      for(j in 1:length(imp)) {
        dat[, j] <- gerb[[1]][[imp[j]]][, var]
        if (length(gerb$ineligibles) > 0) {
          dat[gerb$ineligibles[, var], j] <- NA
        }
      }
    } else {
      dat <- gerb[[1]][[imp]][, var]
    }
    obs <- gerb[[2]][, var]
    
    if(is.element(var, use.log)) {
      if(sum(dat <= 0, na.rm = TRUE) > 0) {
        dat[dat <= 0] <- NA
        #warning(paste0("Non-positive values detected for ", var, ": Cases will be removed prior to log transformation."))
        message(paste0("Non-positive values detected for ", var, ": Cases will be removed prior to log transformation."), "\n", appendLF = FALSE)
      }
      dat <- log(dat)
      var.nam <- paste0("log(",var,")")
    } else {
      var.nam <- var
    }

    cens <- gerb[[2]][, var] == 2
    if(sum(cens, na.rm = TRUE) > 0) {
      #obs <- obs[!cens]
      #dat <- dat[!cens]
      if(length(imp) > 1) {
        dat[cens, ] <- NA
      } else {
        dat[cens] <- NA
      }
      #warning(paste0("Censored values detected when plotting ", var, ": Censored cases will be removed."))
      message(paste0("Censored values detected when plotting ", var, ": Censored cases will be removed."), "\n", appendLF = FALSE)
    }
    
    if (sum(obs == 1, na.rm = TRUE) == 0) {
      cont <- FALSE
      #warning(paste0("Will not produce univariate plot for ", var, ": No imputed values for the variable."))
      message(paste0("Will not produce univariate plot for ", var, ": No imputed values for the variable."), "\n", appendLF = FALSE)
    } else {
      cont <- TRUE
    }
    
    #leg <- c("Observed", "Imputed")
    if (length(imp) == 1) { 
      leg <- c("Observed", "Imputed")
    } else {
      leg <- c("Observed", paste0("Imputed: ", imp))
    }
    ns <- table(obs[!is.na(as.matrix(dat)[, 1])])
    #leg <- paste0(leg, ": n = ", ns[c("0", "1")])
    
    type <- as.character(gerb[[3]][var, "Variable.Type"])
     
    if ((type == "binary" | type == "categorical" | type == "ordinal") & cont) {
      
      if (length(file) == 0) {
        par(mfrow = c(1, 1))
      }
      
      tab <- table(as.matrix(dat)[, 1], obs)
      colnames(tab)[colnames(tab) == "1"] <- as.character(imp[1])
      if (NCOL(dat) > 1) {
        tab.tmp <- tab
        tab <- matrix(NA, NROW(tab), 1 + NCOL(dat))
        rownames(tab) <- rownames(tab.tmp)
        colnames(tab) <- c("0", as.character(imp))
        tab[, colnames(tab.tmp)] <- tab.tmp
        for (j in 2:length(imp)) {
          tab[, as.character(imp[j])] <- table(dat[, j], obs)[, "1"]
        }
      }
    
      ests <- 100 * t(tab)/colSums(tab)
      rownames(ests) <- leg
      #ests <- 100 * rbind(
      #  Observed = tab[, "0"]/sum(tab[, "0"]),
      #  Imputed = tab[, "1"]/sum(tab[, "1"])
      #)
      #if (ests[1] <= ests[length(ests)]) {
      #  loc <- "topleft"
      #} else {
      #  loc <- "topright"
      #}
      loc <- legend.spot
      
      if (sep & length(file) > 0) {
        if (!supp.suff) {
          suff <- paste0("_", var)
        } else {
          suff <- ""
        }
        if (file.type == "pdf") {
          grDevices::pdf(file = paste(file, suff, ".pdf", sep = ""), width = width, height = height)
        } else if (file.type == "png") {
          grDevices::png(file = paste(file, suff, ".png", sep = ""),
                         width = width, height = height, units = "in", res = 500)
        }
      }
      
      if (length(main) == 0) {
        main1 <- var.nam
      } else {
        main1 <- main 
      }
      if (length(ylab) == 0) {
        ylab1 <- "Frequency (%)"
      } else {
        ylab1 <- ylab 
      }
      if (length(xlab) == 0) {
        xlab1 <- ""
      } else {
        xlab1 <- xlab 
      }
      if (length(bty) == 0) {
        bty1 <- "n"
      } else {
        bty1 <- bty 
      }

      if (plot.legend) {
        if (length(legend.text) > 0) {
          leg[1:length(legend.text)] <- legend.text
        }
        barplot(ests, beside = T, legend.text = T, col = col, args.legend = list(legend = leg, bty = bty1, x = loc), main = main1, ylab = ylab1, xlab = xlab1, ...)
      } else {
        barplot(ests, beside = T, col = col, main = main1, ylab = ylab1, xlab = xlab1, ...)
      }
      #grid()
      abline(h = 0)
      
      if (sep & length(file) > 0) {
        grDevices::dev.off()
      }
      
    } else if (type == "semicont" & cont) {
      
      mass <- gerb$mass.final[var]
      
      if (length(file) == 0) {
        #Edited by Max from (c(1, 2))
        par(mfrow = c(2, 1))
      }

      #pos <- as.numeric(dat != mass)
      #tab <- table(pos, obs)
      #ests <- 100 * rbind(
      #  Observed = tab[, "0"]/sum(tab[, "0"]),
      #  Imputed = tab[, "1"]/sum(tab[, "1"])
      #)
      #colnames(ests) = paste0(c("Equals ","Does not equal "),mass)
      #if (ests[1] <= ests[length(ests)]) {
      #  loc <- "topleft"
      #} else {
      #  loc <- "topright"
      #}

      pos <- as.numeric(as.matrix(dat)[, 1] != mass)
      tab <- table(pos, obs)
      colnames(tab)[colnames(tab) == "1"] <- as.character(imp[1])
      if (NCOL(dat) > 1) {
        tab.tmp <- tab
        tab <- matrix(NA, NROW(tab), 1 + NCOL(dat))
        rownames(tab) <- rownames(tab.tmp)
        colnames(tab) <- c("0", as.character(imp))
        tab[, colnames(tab.tmp)] <- tab.tmp
        for (j in 2:length(imp)) {
          pos <- as.numeric(dat[, j] != mass)
          tab[, as.character(imp[j])] <- table(pos, obs)[, "1"]
        }
      }
    
      ests <- 100 * t(tab)/colSums(tab)
      rownames(ests) <- leg
      colnames(ests) <- paste0(c("Equals ","Does not equal "),mass)
      
      loc <- legend.spot

      if (sep & length(file) > 0) {
        if (!supp.suff) {
          suff <- paste0("_", var)
        } else {
          suff <- ""
        }
        if (file.type == "pdf") {
          grDevices::pdf(file = paste(file, suff, "_binary.pdf", sep = ""), width = width, height = height)
        } else if (file.type == "png") {
          grDevices::png(file = paste(file, suff, "_binary.png", sep = ""),
                         width = width, height = height, units = "in", res = 500)
        }
      }
      
      if (length(main) == 0) {
        main1 <- paste0(var.nam, ": Binary")
      } else {
        main1 <- main
      }
      if (length(ylab) == 0) {
        ylab1 <- "Frequency (%)"
      } else {
        ylab1 <- ylab
      }
      if (length(xlab) == 0) {
        xlab1 <- ""
      } else {
        xlab1 <- xlab 
      }
      if (length(bty) == 0) {
        bty1 <- "n"
      } else {
        bty1 <- bty 
      }
      if (length(cex) == 0) {
        #cex1 <- .8
        cex1 <- 1
      } else {
        cex1 <- cex 
      }

      if (plot.legend) {
        if (length(legend.text) > 0) {
          leg[1:length(legend.text)] <- legend.text
        }
        barplot(ests, beside=T, legend.text = T, col = col, 
                args.legend = list(legend = leg, bty = bty1, x = loc, cex = cex1), main = main1, ylab = ylab1, xlab = xlab1, ...)
      } else {
        barplot(ests, beside=T, col = col, main = main1, ylab = ylab1, xlab = xlab1, ...)
      }
      #grid()
      abline(h=0)
      
      if (sep & length(file) > 0) {
        grDevices::dev.off()
      }
      
      dat[dat == mass] <- NA
      dat <- as.matrix(dat)

      d1 <- dat[obs == 0, 1]
      d1 <- d1[!is.na(d1)]
      d1 <- density(d1)

      if(is.element(var, use.log)) {
        d1x <- exp(d1$x)
        log.plot <- "x"        
      } else {
        d1x <- d1$x
        log.plot <- ""
      }      
      d1y <- d1$y

      xlim <- c(min(d1x), max(d1x))
      ylim <- c(min(d1y), max(d1y))

      d2y <- d2x <- list()

      for(j in 1:length(imp)) {

        d2 <- dat[obs == 1, j]
        d2 <- d2[!is.na(d2)]
        d2 <- density(d2)

        if(is.element(var, use.log)) {
          d2x[[j]] <- exp(d2$x)       
        } else {
          d2x[[j]] <- d2$x
        }      
        d2y[[j]] <- d2$y

        xlim.tmp <- c(min(d2x[[j]]), max(d2x[[j]]))
        ylim.tmp <- c(min(d2y[[j]]), max(d2y[[j]]))

        xlim <- c(min(xlim[1], xlim.tmp[1]), max(xlim[2], xlim.tmp[2]))
        ylim <- c(min(ylim[1], ylim.tmp[1]), max(ylim[2], ylim.tmp[2]))

      }

      #leg <- c("Observed","Imputed")
      #ns <- c(length(d1),length(d2))
      #leg <- paste0(leg,": n = ",ns)   
      
      if (sep & length(file) > 0) {
        if (!supp.suff) {
          suff <- paste0("_", var)
        } else {
          suff <- ""
        }
        if (file.type == "pdf") {
          grDevices::pdf(file = paste(file, suff, "_continuous.pdf", sep = ""), width = width, height = height)
        } else if (file.type == "png") {
          grDevices::png(file = paste(file, suff, "_continuous.png", sep = ""),
                         width = width, height = height, units = "in", res = 500)
        }
      }
      
      if (length(main) == 0) {
        main1 <- paste0(var.nam, ": Continuous")
      } else {
        main1 <- main
      }
      if (length(ylab) == 0) {
        ylab1 <- "Density"
      } else {
        ylab1 <- ylab
      }
      if (length(xlab) == 0) {
        xlab1 <- ""
      } else {
        xlab1 <- xlab 
      }

      plot(d1x, d1y, xlim = xlim, ylim = ylim, lty = lty[1], col = col[1], type = "l", xlab = xlab1, ylab = ylab1, 
           main = main1, lwd = lwd[1], log = log.plot, panel.first = grid(), ...)
      for(j in 1:length(imp)) {
        lines(d2x[[j]], d2y[[j]], xlim = xlim, ylim = ylim, lty = lty[j + 1], col = col[j + 1], type = "l", lwd = lwd[j + 1])
      }
      
      #Edited by Max below, to include the cex option
      if (plot.legend) {
        if (length(legend.text) > 0) {
          leg[1:length(legend.text)] <- legend.text
        }
        legend(legend.spot, legend = leg, bty = bty1, col = col, lty = lty, lwd = lwd, cex = cex1, ...)
      }
      
      if (sep & length(file) > 0) {
        grDevices::dev.off()
      }
      
    } else if (cont) {
      
      if (length(file) == 0) {
        par(mfrow = c(1, 1))
      }
     
      dat <- as.matrix(dat)

      d1 <- dat[obs == 0, 1]
      d1 <- d1[!is.na(d1)]
      d1 <- density(d1)

      if(is.element(var, use.log)) {
        d1x <- exp(d1$x)
        log.plot <- "x"        
      } else {
        d1x <- d1$x
        log.plot <- ""
      }      
      d1y <- d1$y

      xlim <- c(min(d1x), max(d1x))
      ylim <- c(min(d1y), max(d1y))

      d2y <- d2x <- list()

      for(j in 1:length(imp)) {

        d2 <- dat[obs == 1, j]
        d2 <- d2[!is.na(d2)]
        d2 <- density(d2)

        if(is.element(var, use.log)) {
          d2x[[j]] <- exp(d2$x)       
        } else {
          d2x[[j]] <- d2$x
        }      
        d2y[[j]] <- d2$y

        xlim.tmp <- c(min(d2x[[j]]), max(d2x[[j]]))
        ylim.tmp <- c(min(d2y[[j]]), max(d2y[[j]]))

        xlim <- c(min(xlim[1], xlim.tmp[1]), max(xlim[2], xlim.tmp[2]))
        ylim <- c(min(ylim[1], ylim.tmp[1]), max(ylim[2], ylim.tmp[2]))

      }

      #lty <- c(imp.lty, obs.lty)
      #lty <- c(obs.lty, ltys)
      #col <- c(1, 2)
      
      if (sep & length(file) > 0) {
        if (!supp.suff) {
          suff <- paste0("_", var)
        } else {
          suff <- ""
        }
        if (file.type == "pdf") {
          grDevices::pdf(file = paste(file, suff, ".pdf", sep = ""), width = width, height = height)
        } else if (file.type == "png") {
          grDevices::png(file = paste(file, suff, ".png", sep = ""),
                         width = width, height = height, units = "in", res = 500)
        }
      }

      if (length(main) == 0) {
        main1 <- var.nam
      } else {
        main1 <- main 
      }
      if (length(ylab) == 0) {
        ylab1 <- "Density"
      } else {
        ylab1 <- ylab 
      } 
      if (length(bty) == 0) {
        bty1 <- "n"
      } else {
        bty1 <- bty 
      }
      if (length(cex) == 0) {
        cex1 <- .8
      } else {
        cex1 <- cex 
      }
      if (length(xlab) == 0) {
        xlab1 <- ""
      } else {
        xlab1 <- xlab 
      }

      plot(d1x, d1y, xlim = xlim, ylim = ylim, lty = lty[1], col = col[1], type = "l", xlab = xlab1, 
           ylab = ylab1, main = main1, lwd = lwd[1], log = log.plot, panel.first = grid(), ...)
      for(j in 1:length(imp)) {
        lines(d2x[[j]], d2y[[j]], xlim = xlim, ylim = ylim, lty = lty[j + 1], col = col[j + 1], type = "l", lwd = lwd[j + 1])
      }
      if (plot.legend) {
        if (length(legend.text) > 0) {
          leg[1:length(legend.text)] <- legend.text
        }
        legend(legend.spot, legend = leg, bty = bty1, col = col, lty = lty, lwd = lwd, ...)
      }
      
      if (sep & length(file) > 0) {
        grDevices::dev.off()
      }
    }
  }
  if (!sep & length(file) > 0) {
    grDevices::dev.off()
  }
}


plot_gerbil.mult <- function(gerb, file = NULL, vars = NULL, imp = 1, obs.col = "blue4", imp.col = "brown2", sep = FALSE, legend.text = NULL, 
                             legend.spot = "bottomleft", height = NULL, width = NULL, mfrow = c(3, 2), pch = c(1, 3), use.log = NULL, 
                             lwd = NULL, main = NULL, xlab = NULL, ylab = NULL, bty = NULL, cex = NULL, lty = NULL, xaxt = "s", yaxt = "s", 
                             partial = "imputed", ...) {
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  plot.legend <- TRUE
  if (is.logical(legend.text)) {
    if(!legend) {
      plot.legend <- FALSE
    } else {
      legend.text <- NULL
    }
  } else if (length(legend.text) == 1) {
    if (legend.text == "n") {
      plot.legend <- FALSE
    } else if (legend.text == "y") {
      legend.text <- NULL
    }
  }
  
  if (is.logical(imp)) {
    if (imp) {
      imp <- as.numeric(dimnames(gerb)[[3]])
    } else {
      imp <- 1
    }
  }
  if(length(imp) == 0) {
    imp <- 1
  }
  
  if(length(imp.col) == 1) {
    imp.col <- rep(imp.col, length(imp))
    imp.col <- imp.col[1:length(imp)]
  }

  col <- c(obs.col, imp.col)

  obs.pch <- pch[1]
  if (length(pch) == 1) {  
    imp.pch <- pch[1]
  } else {
    imp.pch <- pch[-1]
  }

  if(length(imp.pch) == 1) {
    imp.pch <- rep(imp.pch, length(imp))
    imp.pch <- imp.pch[1:length(imp)]
  }

  pch <- c(obs.pch, imp.pch)  

  obs.lwd <- lwd[1]
  if (length(lwd) == 1) {  
    imp.lwd <- lwd[1]
  } else {
    imp.lwd <- lwd[-1]
  }

  if(length(imp.lwd) == 1) {
    imp.lwd <- rep(imp.lwd, length(imp))
    imp.lwd <- imp.lwd[1:length(imp)]
  }

  lwd <- c(obs.lwd, imp.lwd)

  obs.lty <- lty[1]
  if (length(lty) == 1) {  
    imp.lty <- lty[1]
  } else {
    imp.lty <- lty[-1]
  }

  if(length(imp.lty) == 1) {
    imp.lty <- rep(imp.lty, length(imp))
    imp.lty <- imp.lty[1:length(imp)]
  }

  lty <- c(obs.lty, imp.lty)

  #imp <- imp[1]
  
  if (is.logical(vars)) {
      if (vars) {
          vars <- colnames(gerb$missing)[colSums(gerb$missing == 1 | gerb$missing == 2) > 0]
      }
      else {
          vars <- colnames(gerb$missing)
      }
  }
  else if (length(vars) == 0) {
      vars <- colnames(gerb$missing)
  }
  
  if (length(file) == 0 | length(vars) == 2) {
    sep <- TRUE
    supp.suff <- TRUE
  } else {
    supp.suff <- FALSE
  }
  
  if (sep) {
    mfrow <- c(1, 1)
  }
  
  if(length(vars) <= 1) {
    stop("Need at least 2 variables for bivariate plotting.")
  }
  
  if(length(height) == 0) {
    if (!sep) {
      height <- 11
    } else {
      height <- 6.5
    }
  }
  
  if(length(width) == 0) {
    if (!sep) {
      width <- 8.5
    } else {
      width <- 6.5
    }
  }
  
  if (is.logical(use.log)) {
    if(use.log) {
      use.log <- vars
    }
  }
  
  if (length(file) == 0) {
    #graphics::par(mfrow = c(1, 1), ask = FALSE)
  } else {
    file.type <- "pdf"
    if (substr(file, nchar(file) - 3, nchar(file)) ==
        ".pdf") {
      file <- substr(file, 1, nchar(file) - 4)
    } else if (substr(file, nchar(file) - 3, nchar(file)) ==
               ".png") {
      file <- substr(file, 1, nchar(file) - 4)
      file.type <- "png"
    }
    if (!sep) {
      if (file.type == "pdf") {
        file <- paste(file, ".pdf", sep = "")
        grDevices::pdf(file = file, width = width, height = height)
      } else if (file.type == "png") {
        file <- paste(file, ".png", sep = "")
        grDevices::png(file = file, width = width, height = height,
                       units = "in", res = 500)
      }
      graphics::par(mfrow = mfrow, ask = FALSE)
    }
  }

  if (length(imp) == 1) { 
    leg <- c("Observed", "Imputed")
    leg <- c("Obs.", "Imp.")
  } else {
    leg <- c("Observed", paste0("Imputed ", imp))
  }
  #ns <- table(obs[!is.na(as.matrix(dat)[, 1])])

  if (length(legend.text) > 0) {
    leg[1:length(legend.text)] <- legend.text
  }
  
  plot.num <- 0
  neg.vars <- NA

  for(j in 2:length(vars)) {
    
    for(i in 1:(j-1)) {
      
      types <- gerb[[3]][c(vars[i], vars[j]), "Variable.Type"]
      #vars1 <- rownames(gerb[[3]])[c(i, j)]
      vars1 <- c(vars[i], vars[j])
      
      is.box <- types == "categorical" | types == "ordinal" | types == "binary"
      
      if(!is.box[1] & is.box[2]) {
        var1 <- vars1[2]
        var2 <- vars1[1]
      } else if (is.box[1] & is.box[2]) {

        dat1 <- gerb[[1]][[imp[1]]][, vars1[1]]
        dat2 <- gerb[[1]][[imp[1]]][, vars1[2]]

        if (length(gerb$ineligibles) > 0) {
          dat1[gerb$ineligibles[, vars1[1]]] <- NA
          dat2[gerb$ineligibles[, vars1[2]]] <- NA
        }

        #n1 <- length(table(gerb[[1]][[imp]][, vars1[1]]))
        #n2 <- length(table(gerb[[1]][[imp]][, vars1[2]]))
        n1 <- length(table(dat1))
        n2 <- length(table(dat2))

        if(n2 > n1) {
          var2 <- vars1[2]
          var1 <- vars1[1]
        } else {
          var2 <- vars1[1]
          var1 <- vars1[2]
        }
      } else {
        var1 <- vars1[1]
        var2 <- vars1[2]
      }
      
      if (partial == "imputed") {
        miss <- gerb[[2]][, var1] == 1 | gerb[[2]][, var2] == 1
      } else {
        miss <- gerb[[2]][, var1] == 1 & gerb[[2]][, var2] == 1
      }

      if (length(imp) > 1) {
        dat1 <- dat2 <- matrix(NA, NROW(gerb[[1]][[imp[1]]]), length(imp))
        for(k in 1:length(imp)) {
          c.d1 <- class(gerb[[1]][[imp[k]]][, var1])[1]
          c.d2 <- class(gerb[[1]][[imp[k]]][, var2])[1]
          if (c.d1 == "character" | c.d1 == "factor" | c.d1 == "ordered") {
            d1.tmp <- as.character(gerb[[1]][[imp[k]]][, var1])
          } else {
            d1.tmp <- gerb[[1]][[imp[k]]][, var1]
          }
          if (c.d2 == "character" | c.d2 == "factor" | c.d2 == "ordered") {
            d2.tmp <- as.character(gerb[[1]][[imp[k]]][, var2])
          } else {
            d2.tmp <- gerb[[1]][[imp[k]]][, var2]
          }
          dat1[, k] <- d1.tmp
          dat2[, k] <- d2.tmp
          if (length(gerb$ineligibles) > 0) {
            dat1[gerb$ineligibles[, var1], k] <- NA
            dat2[gerb$ineligibles[, var2], k] <- NA
          }
        }
      } else {
        c.d1 <- class(gerb[[1]][[imp]][, var1])[1]
        c.d2 <- class(gerb[[1]][[imp]][, var2])[1]
        if (c.d1 == "character" | c.d1 == "factor" | c.d1 == "ordered") {
          d1.tmp <- as.character(gerb[[1]][[imp]][, var1])
        } else {
          d1.tmp <- gerb[[1]][[imp]][, var1]
        }
        if (c.d2 == "character" | c.d2 == "factor" | c.d2 == "ordered") {
          d2.tmp <- as.character(gerb[[1]][[imp]][, var2])
        } else {
          d2.tmp <- gerb[[1]][[imp]][, var2]
        }
        dat1 <- d1.tmp
        dat2 <- d2.tmp

        if (length(gerb$ineligibles) > 0) {
          dat1[gerb$ineligibles[, var1]] <- NA
          dat2[gerb$ineligibles[, var2]] <- NA
        }
      }

      cens <- gerb[[2]][, var1] == 2 | gerb[[2]][, var2] == 2
      if(sum(cens, na.rm = TRUE) > 0) {
        #miss <- miss[!cens]
        #dat1 <- dat1[!cens]
        #dat2 <- dat2[!cens]
        if(length(imp) > 1) {
          dat1[cens, ] <- NA
          dat2[cens, ] <- NA
        } else {
          dat1[cens] <- NA
          dat2[cens] <- NA
        }
        #warning(paste0("Censored values detected when plotting ", var1, " and ", var2, ": Censored cases will be removed."))
        message(paste0("Censored values detected when plotting ", var1, " and ", var2, ": Censored cases will be removed."), "\n", appendLF = FALSE)
      }
      
      if (sum(miss, na.rm = TRUE) == 0) {
        cont <- FALSE
        #warning(paste0("Will not produce bivariate plot between ", var1, " and ", var2, ": No imputed values for either variable."))
        message(paste0("Will not produce bivariate plot between ", var1, " and ", var2, ": No imputed values for either variable."), "\n", appendLF = FALSE)
      } else {
        cont <- TRUE
      }
      if(is.element(var1, use.log)) {
        if (sum(dat1 <= 0, na.rm = TRUE) > 0) {
          dat1[dat1 <= 0] <- NA
          if (!is.element(var1, neg.vars)) {
            #warning(paste0("Non-positive values detected for ", var1, ": Cases will be removed prior to log transformation."))
            message(paste0("Non-positive values detected for ", var1, ": Cases will be removed prior to log transformation."), "\n", appendLF = FALSE)
            neg.vars[(length(neg.vars) + 1)] <- var1 
          }
        }
        dat1 <- log10(dat1)
        var.nam1 <- paste0("log(",var1,")")
      } else {
        var.nam1 <- var1
      }
      
      if(is.element(var2, use.log)) {
        if (sum(dat2 <= 0, na.rm = TRUE) > 0) {
          dat2[dat2 <= 0] <- NA
          if (!is.element(var2, neg.vars)) {
            #warning(paste0("Non-positive values detected for ", var2, ": Cases will be removed prior to log transformation."))
            message(paste0("Non-positive values detected for ", var2, ": Cases will be removed prior to log transformation."), "\n", appendLF = FALSE)
            neg.vars[(length(neg.vars) + 1)] <- var2 
          }
        }
        dat2 <- log10(dat2)
        var.nam2 <- paste0("log(",var2,")")
      } else {
        var.nam2 <- var2
      }
      
      if (sum(is.box) == 2) {
        #warning(paste0("Will not produce bivariate plot between ", var1, " and ", var2, ": Neither variable is continuous."))
        #message(paste0("Will not produce bivariate plot between ", var1, " and ", var2, ": Neither variable is continuous."), "\n", appendLF = FALSE)
        #cont <- FALSE
      }
      
      #leg <- c("Observed", "Imputed")
      #leg <- paste0(leg, ": n = ", ns[c("0", "1")])
      
      if (!supp.suff) {
        suff <- paste0("_", var2, "_vs_", var1)
      } else {
        suff <- ""
      }
      
      if (sep & length(file) > 0 & cont) {
        if (file.type == "pdf") {
          grDevices::pdf(file = paste(file, suff, ".pdf", sep = ""),
                         width = width, height = height)
        } else if (file.type == "png") {
          grDevices::png(file = paste(file, suff, ".png", sep = ""),
                         width = width, height = height, units = "in", res = 500)
        }
      }

      if (sum(is.box) == 0) {
        num.plots <- length(imp)
      } else  {
        num.plots <- 1
      }
      plot.num <- plot.num + cont * num.plots
      
      if(sum(is.box) == 0 & cont) {
        
        ylim <- c(min(dat2, na.rm = TRUE), max(dat2, na.rm = TRUE))
        xlim <- c(min(dat1, na.rm = TRUE), max(dat1, na.rm = TRUE))
        
        #pch <- c("o","x")
        #col <- c(2,1)

        if (length(ylab) == 0) {
          #ylab1 <- var.nam2
          ylab1 <- var2
        } else {
          ylab1 <- ylab 
        }
        if (length(bty) == 0) {
          bty1 <- "n"
        } else {
          bty1 <- bty 
        }
        if (length(xlab) == 0) {
          #xlab1 <- var.nam1
          xlab1 <- var1
        } else {
          xlab1 <- xlab 
        }

        dat1 <- as.matrix(dat1)
        dat2 <- as.matrix(dat2)

        if(is.element(var1, use.log)) {
          x1 <- exp(dat1[!miss, ])
          x2 <- exp(dat1[miss, ])
          log.plot <- "x"
        } else {
          x1 <- dat1[!miss, ]
          x2 <- dat1[miss, ]
          log.plot <- ""
        }    
        
        if(is.element(var2, use.log)) {
          y1 <- exp(dat2[!miss, ])
          y2 <- exp(dat2[miss, ])
          log.plot <- paste0(log.plot, "y")
        } else {
          y1 <- dat2[!miss, ]
          y2 <- dat2[miss, ]
        }    

        if (length(main) == 0) {
          main1 <- paste0(var.nam2," vs. ", var.nam1)
        } else {
          main1 <- main 
        }

        if (NCOL(dat1) > 1) {
          message(paste0("Both variables ", var1, " and ", var2, " are continuous: Separate plots will be made for each element of 'imp'."), "\n", appendLF = FALSE)
        }

        for(k in 1:NCOL(dat1)) {
          if (NCOL(dat1) > 1) {
            main2 <- paste0(leg[k + 1], ": ", main1)
          } else {
            main2 <- main1
          }
          #plot(dat1[!miss], dat2[!miss], pch = pch[1], col = col[1], xlim = xlim, ylim = ylim, main = main2, xlab = xlab1, ylab = ylab1, ...)
          plot(dat1[!miss, k], dat2[!miss, k], pch = pch[1], col = col[1], xlim = xlim, ylim = ylim, main = main2, xlab = xlab1, ylab = ylab1, xaxt = "n", yaxt = "n", ...)
          ats1 <- par("xaxp")
          ats2 <- par("yaxp")
          ats1 <- seq(ats1[1], ats1[2], length.out = ats1[3] + 1)
          if(is.element(var1, use.log)) {
            labs1 <- 10^ats1
          } else {
            labs1 <- ats1
          }
          if (is.element(xaxt, c("s", "l", "t"))) {
            axis(side = 1, at = ats1, labels = labs1)
          }
          ats2 <- seq(ats2[1], ats2[2], length.out = ats2[3] + 1)
          if(is.element(var2, use.log)) {
            labs2 <- 10^ats2
          } else {
            labs2 <- ats2
          }
          if (is.element(yaxt, c("s", "l", "t"))) {
            axis(side = 2, at = ats2, labels = labs2)
          }
          abline(lm(dat2[!miss, k] ~ dat1[!miss, k]), lty = lty[1], col = col[1], lwd = lwd[1])
          points(dat1[miss, k], dat2[miss, k], pch = pch[k + 1], col = col[k + 1])
          abline(lm(dat2[miss, k] ~ dat1[miss, k]), lty = lty[k + 1], col = col[k + 1], lwd = lwd[k + 1])
          if (plot.legend) {
            legend(legend.spot, legend = leg[c(1, k + 1)], bty = bty1, col = col[c(1, k + 1)], pch = pch[c(1, k + 1)], ...)
          }
        }
        
      } else if (sum(is.box) == 1 & cont) {

        tmp.df <- data.frame(Y = as.matrix(dat2)[, 1], X = as.matrix(dat1)[, 1], R = miss)
        if (NCOL(dat1) > 1) {
          n.m <- sum(miss)
          tmp.df1 <- matrix(NA, n.m * (NCOL(dat1) - 1), 3)
          colnames(tmp.df1) <- colnames(tmp.df)
          for (k in 2:NCOL(dat1)) {
            tmp.df1[(((k - 1) - 1) * n.m + 1):((k - 1) * n.m), ] <- cbind(Y = dat2[miss, k], X = dat1[miss, k], R = k)
          }
          tmp.df[, "R"] <- as.numeric(tmp.df[, "R"])
          tmp.df1 <- data.frame(tmp.df1)
          tmp.df <- rbind(tmp.df, tmp.df1)
        }
        tmp.df <- tmp.df[!is.na(tmp.df$Y) & !is.na(tmp.df$X),]
        tmp.df$Y <- as.numeric(as.character(tmp.df$Y))
        tab <- table(tmp.df$X)
        p <- length(tab)
        p1 <- NCOL(dat1) + 1

        sep1 <- 1
        at <- 1:(p1 * p) + c(t(matrix(rep(sep1 * (0:(p - 1)), p1), p, p1)))
        here <- matrix(at, p1, p)
        here <- colMeans(here)

        if (length(main) == 0) {
          main1 <- paste0(var.nam2, " vs. ", var.nam1)
        } else {
          main1 <- main 
        }
        if (length(ylab) == 0) {
          #ylab1 <- var.nam2
          ylab1 <- var2
        } else {
          ylab1 <- ylab 
        } 
        if (length(bty) == 0) {
          bty1 <- "n"
        } else {
          bty1 <- bty 
        }
        if (length(xlab) == 0) {
          xlab1 <- var.nam1
        } else {
          xlab1 <- xlab 
        }

        boxplot(Y ~ R + X, data = tmp.df, at = at, col = col, xaxt = "n", yaxt = "n", ylab = ylab1, xlab = xlab1, main = main1, ...)
        ats2 <- par("yaxp")
        ats2 <- seq(ats2[1], ats2[2], length.out = ats2[3] + 1)

        if(is.element(var2, use.log)) {
          labs2 <- 10^ats2
        } else {
          labs2 <- ats2
        }
        if (is.element(yaxt, c("s", "l", "t"))) {
          axis(side = 2, at = ats2, labels = labs2)
        }
        if (is.element(xaxt, c("s", "l", "t"))) {
          axis(side = 1, at = here, labels = names(tab), tick = TRUE)
        }
        if (plot.legend) {
          legend(legend.spot, legend = leg, fill = col, bty = bty1, ...)
        }
      } else if (sum(is.box) == 2 & cont) {
        
        #miss <- gerb[[2]][, var1] == 1 | gerb[[2]][, var2] == 1
        #dat1 <- gerb[[1]][[imp]][, var1]
        #dat2 <- gerb[[1]][[imp]][, var2]
        
        tab <- table(as.matrix(dat1)[, 1], as.matrix(dat2)[, 1], miss)
        dimnames(tab)[[3]] <- c("0","1")
        if (NCOL(dat1) > 1) {
          p1 <- NCOL(dat1)
          tab.tmp <- array(NA, c(dim(tab)[1], dim(tab)[2], p1 + 1))
          dimnames(tab.tmp) <- list(dimnames(tab)[[1]], dimnames(tab)[[2]], 0:p1)
          tab.tmp[, , 1:2] <- tab
          for(k in 2:p1) {
            tab.tmp1 <- table(dat1[, k], dat2[, k], miss)
            tab.tmp[, , k + 1] <- tab.tmp1[, , "TRUE"]        
          }
          tab <- tab.tmp
        }

        sums.tab <- apply(tab, 3, sum)
        for (k in 1:dim(tab)[3]) {
          tab[, , k] <- tab[, , k]/sums.tab[k]       
        }
        tab <- 100 * tab
        
        nams1 <- dimnames(tab)[[1]]
        nams2 <- dimnames(tab)[[2]]
        nams3 <- dimnames(tab)[[3]]
        
        if(sum(is.na(as.numeric(nams1))) == 0) {
          nams1 <- as.numeric(nams1)
        }
        if(sum(is.na(as.numeric(nams2))) == 0) {
          nams2 <- as.numeric(nams2)
        }
        
        x1.tmp <- factor(rep(nams1, length(nams2)))
        x2.tmp <- factor(c(t(matrix(nams2, length(nams2), length(nams1)))))
        x3.tmp <- matrix(nams3, length(nams3), length(x1.tmp))
        x3.tmp <- c(t(x3.tmp))
        x1.tmp <- rep(x1.tmp, length(nams3))
        x2.tmp <- rep(x2.tmp, length(nams3))
        tmp <- data.frame(x1 = x1.tmp, x2 = x2.tmp, x3 = x3.tmp)
        #tmp <- data.frame(
        #  x1 = factor(rep(nams1, length(nams2))),
        #  x2 = factor(c(t(matrix(nams2, length(nams2), length(nams1)))))
        #)
        #tmp <- rbind(data.frame(tmp, x3 = nams3[1]), data.frame(tmp, x3 = nams3[2]))
        
        tmp1 <- rep(NA,NROW(tmp))
        for(ind in 1:NROW(tmp)) {
          tmp1[ind] <- tab[tmp[ind, 1], tmp[ind, 2], tmp[ind, 3]]
        }
        
        tmp[, 3] <- as.character(tmp[, 3])
        #tmp[tmp[, 3] == "FALSE", 3] <- "Observed"
        #tmp[tmp[, 3] == "TRUE", 3] <- "Imputed"
        tmp[, 3] <- leg[match(tmp[, 3], nams3)]
        colnames(tmp) <- c(var1, var2, "Status")
        
        tmp2 <- data.frame(tmp, Frequency = tmp1)
        tmp2$Status <- factor(tmp2$Status, levels = leg)     

        # require(lattice)
        
        #pch = c(1, 3)
        #pch <- c("o","x")
        #cols=c("blue4", "brown2")
        
        form <- paste0(var2, "~", "Frequency", "|", var1)
        form <- as.formula(form)
        
        if (length(file) == 0) {
          size.tmp <- dev.size()
          dev.off()
          dev.new(height = size.tmp[2], width = size.tmp[1])
        }
        
        if(!sep) {
          plot.new()
          num.tmp <- plot.num %% prod(mfrow)
          if(num.tmp == 0) {
            num.tmp <- prod(mfrow)
          }
          n.y <- ceiling(num.tmp/mfrow[2])
          n.x <- num.tmp - mfrow[2]*(n.y - 1)
          if (length(cex) == 0) {
            cex2 <- 0.7
            cex1 <- 0.8
          } else {
            cex2 <- cex
            cex1 <- cex
          }
        } else {
          if (length(cex) == 0) {
            cex2 <- 1
            cex1 <- 1
          } else {
            cex2 <- cex
            cex1 <- cex
          }
          n.x <- n.y <- 1
        }

        if (length(main) == 0) {
          main1 <- paste0(var.nam2," vs. ", var.nam1)
        } else {
          main1 <- main 
        }
        if (length(bty) == 0) {
          bty1 <- "n"
        } else {
          bty1 <- bty 
        }
        if (length(xlab) == 0) {
          xlab1 <- "Frequency (%)"
        } else {
          xlab1 <- xlab 
        }
        if (length(ylab) == 0) {
          ylab1 <- ""
        } else {
          ylab1 <- ylab 
        }
        if (plot.legend) {
          Status <- tmp2$Status
          plot(
            lattice::dotplot(form,
                             data = tmp2,
                             horizontal = TRUE,
                             #panel = panel.superpose,
                             group = Status, pch = pch, col = col,
                             scales = list(cex = cex2),
                             par.strip.text = list(cex = cex2),
                             key = list(space = "right",
                                        #transparent = TRUE,
                                        points = list(pch = pch, col = col),
                                        text = list(leg, cex = cex2),
                                        cex.title = cex2, cex = cex2),
                             main = list(main1, cex = cex1), xlab = list(xlab1, cex = cex2), ylab = list(ylab1, cex = cex2), cex = cex2),
            split = c(n.x, n.y, mfrow[2], mfrow[1]), cex = cex2,
            newpage = FALSE
          )
        } else {
          Status <- tmp2$Status
          plot(
            lattice::dotplot(form,
                             data = tmp2,
                             horizontal = TRUE,
                             #panel = panel.superpose,
                             group = Status, pch = pch, col = col,
                             scales = list(cex = cex2),
                             par.strip.text = list(cex = cex2),
                             main = list(main1, cex = cex1), xlab = list(xlab1, cex = cex2), ylab = list(ylab1, cex = cex2), cex = cex2),
            split = c(n.x, n.y, mfrow[2], mfrow[1]), cex = cex2,
            newpage = FALSE
          )
        }
      }
      if (sep & length(file) > 0 & cont) {
        grDevices::dev.off()
      }
    }
  }
  if (!sep & length(file) > 0) {
    grDevices::dev.off()
  }
}


#' @title Plotting for gerbil objects
#'
#' @description
#' Using a \code{gerbil} object as an input, this function gives
#' diagnostic plots for selected variables
#'
#' @details
#' Three types of plots may be produced:
#' 1) Univariate (produced by setting \code{type = 1}): Compares the marginal distribution of observed and imputed values of a given variable.  Density plots are produced for continuous variables, and bar plots are given for binary, categorical, and ordinal variables.  For semi-continuous variables, two plots are constructed: a) a bar plot for the binary portion of the variable and 2) a density plot for the continuous portion.
#' 2) Bivariate (produced by setting \code{type = 2}): Compares the bivariate distributions of observed and imputed values of two variables.  Scatter plots are produced if both variables are continuous or semi-continuous, box plots are produced if one variable is continuous or semi-continuous and the other is not, and a lattice plot is produced if neither variable is continuous or semi-continuous. For bivariate plots, imputed observations are those that have one or more of the values of the pair missing within the original dataset.
#' 3) Trace lines (produced by setting \code{type = 3}): Plots a pre-specified parameter across iterations of MCMC in order to examine convergence for a given variable.  Parameters that may be plotted include means (\code{trace.type = 1}) and variances (\code{trace.type = 2}).
#'
#' Multiple plots may be created, as determined by the variable names listed in the parameter \code{y}. For univariate and trace plots, one plot is created for
#' each variable listed in \code{y}. For bivariate plotting, one plot is created for each combination of two elements within the vector \code{y} (as such, \code{y} must have a length of at least two in this case).
#' For trace plotting, elements of \code{y} should correspond to column names in the dataset that has been expanded to include binary indicators for categorical and semi-continuous variables.
#' If multiple plots are to be created, it is recommended to specify a file for output using the parameter \code{file}, in which case separate
#' files will be created for each plot (if \code{sep = TRUE}) or all plots will be written to the same file (if \code{sep = FALSE}).
#'
#' The only required input is a parameter \code{x} which is a \code{gerbil} object.
#'
#' @param x A \code{gerbil} object containing the imputed data.
#' @param y A vector listing the column names of the imputed data for which plots should be created. See details. By default, \code{y} contains all columns of the data that required imputation.
#' @param type A scalar used to specify the type of plots that will be created.  Options include univariate (marginal) plots (\code{type = 1}), bivariate plots (\code{type = 2}), and trace plots (\code{type = 3}). See details. Defaults to \code{type = 1}.
#' @param imp A scalar or vector indicating which of the multiply imputed datasets should be used for plotting.  Defaults to \code{imp = 1}. Setting \code{imp = TRUE} will include all imputed datasets.
#' @param col The color used for plotting -- should be a vector of length equal to \code{imp + 1}. The first element references plotting of observed data, and remaining elements reference plotting of imputed data. 
#' @param lty The line type used for plotting imputed values with trace lines or density plots -- should be a vector of length equal to \code{imp + 1}. The first element references plotting of observed data, and remaining elements reference plotting of imputed data.
#' @param lwd The line width used for density and trace line plotting -- should be a vector of length equal to \code{imp + 1}. The first element references plotting of observed data, and remaining elements reference plotting of imputed data. 
#' @param pch A length-2 vector that indicates the plotting symbol to be used for imputed and observed values in scatter and lattice plots.
#' @param log A character vector that includes names of variables of which a log transformation is to be taken prior to plotting.
#' @param legend A character or expression vector to appear in the legend. If \code{FALSE} or \code{'n'}, no legend is created. Defaults to \code{c("Observed", "Imputed", ...)}.
#' @param legend.loc The location of the legend in the plots. 
#' @param mfrow The layout of plots across a single page when there are to be multiple plots per page (as is the case when \code{file} is non-\code{NULL} and \code{sep = FALSE}).
#' @param trace.type The type of trace plot to be created (only valid when \code{type = 3}).  See details.  Defaults to \code{trace.type = 1}.
#' @param file A character string giving the name of file that will be created in the home directory containing plots. The name should have a \code{.pdf} or \code{.png} extension. If \code{NULL} (the default), no file is created. 
#' @param sep If \code{sep = TRUE}, separate plots will be generated for each outcome.  Applicable only if plots are saved to file (\code{plot.file} is \code{non-NULL}). To change display of plots produced as output, use \code{\link[graphics]{par}}.
#' @param height The height of the graphics region (in inches) when a pdf is created.
#' @param width The width of the graphics region (in inches) when a pdf is created.
#' @param partial Indicates how partially imputed pairs are handled in bivariate plotting. If \code{'imputed'}, cases with at least one missing variable in a pair are considered imputed. Otherwise (\code{partial = 'observed'}), only cases with both variables in the pair missing are considered imputed.  
#' @param ... Arguments to be passed to methods, such as \code{plot}.
#' 
#' @return No returned value, but instead plots are generated in the workspace or written to a specified directory. 
#'
#' @examples
#' \donttest{
#' #Load the India Human Development Survey-II dataset
#' data(ihd_mcar) 
#' 
#' # Create a gerbil object
#' 
#' imps.gerbil <- gerbil(ihd_mcar, m = 1, ords = "education_level", semi = "farm_labour_days", 
#'        bincat = "job_field")
#'
#' # Univariate plotting of all variables to a file
#' plot(imps.gerbil, type = 1, file = file.path(tempdir(), "gerbil_univariate.pdf"))
#' 
#' # Bivariate plotting of all variables to a file
#' plot(imps.gerbil, type = 2, file = file.path(tempdir(), "gerbil_bivariate.pdf"))
#' 
#' # Trace plotting of all variables to a file
#' plot(imps.gerbil, type = 3, file = file.path(tempdir(), "gerbil_ts.pdf"))
#' 
#' # Univariate plotting of one variable (not to a file)
#' plot(imps.gerbil, type = 1, y = "job_field")
#' 
#' # Bivariate plotting of one pair of variables (not to a file)
#' plot(imps.gerbil, type = 2, y = c("job_field", "income"))
#' 
#' # Bivariate plotting of one pair of variables (not to a file) with income logged
#' plot(imps.gerbil, type = 2, y = c("job_field", "income"), log = "income")
#' }
#' 
#' @export
#'
#'
#' @importFrom stats as.formula density
#' @importFrom graphics abline axis barplot boxplot grid legend lines par plot plot.new points
#' @importFrom grDevices dev.new dev.off dev.size
#' 
plot.gerbil <- function(x, y = NULL, type = "Univariate", imp = 1, col = NULL, lty = NULL, lwd = NULL, pch = NULL, 
                        log = NULL, legend = NULL, legend.loc = "topright", mfrow = c(3, 2), trace.type = "Mean", 
                        file = NULL, sep = FALSE, height = NULL, width = NULL, partial = "imputed", ...)
{

  if(length(file) == 0) {
    #if ((tolower(type) == "univariate" | type == 1) & length(y) != 1) {
    #  warning("y should have length 1 for univariate plotting.")
    #} else 
    if ((tolower(type) == "bivariate" | tolower(type) == "multivariate" | type == 2) & length(y) < 2) {
      warning("y should have length of at least 2 for bivariate plotting.")
    } 
    #else if ((tolower(type) == "ts" | tolower(type) == "time.series" | type == 3) & length(y) != 1) {
    #  warning("y should have length 1 for trace plotting.")
    #} 
  }

  if (length(pch) == 0) {
    pch <- c(1, 3)
  } else if (length(pch) == 1) {
    pch <- c(pch, pch)
  }
  if (length(col) == 0) {
    col <- c("blue4", "brown2")
  } else if (length(col) == 1) {
    col <- rep(col, 2)
  }
  if (length(lty) == 0) {
    lty <- c(2, 1)
  } else if (length(lty) == 1) {
    lty <- rep(lty, 2)
  }

  obs.col <- col[1]
  imp.col <- col[-1]
  
  obs.lty <- lty[1]
  imp.lty <- lty[-1]

  if (length(lwd) <= 1) {
    imp.lwd <- obs.lwd <- lwd
  } else {
    obs.lwd <- lwd[1]
    imp.lwd <- lwd[-1]
  }

  plot_gerbil(gerb = x, vars = y, type = type, 
              imp = imp, trace.type = trace.type, log = log, 
              obs.col = obs.col, imp.col = imp.col, legend.text = legend, 
              legend.spot = legend.loc, obs.lty = obs.lty, imp.lty = imp.lty, 
              obs.lwd = obs.lwd, imp.lwd = imp.lwd, pch = pch, 
              file = file, sep = sep, mfrow = mfrow, height = height, 
              width = width, partial = partial, ...)

}