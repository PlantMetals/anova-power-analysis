#--------------------------------------------------------------------------
#
#             FINDING THE REQUIRED SAMPLE SIZE AND ESTIMATING
#                       THE POWER OF A ONE-WAY ANOVA
#
#                               Version 1.00
#                              September 2026
#
#                          by Dr. Filip Poscic
#                     https://github.com/PlantMetals 
#
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# 1. Information about the script version and the author
# --------------------------------------------------------------------------

author <- "Dr. Filip Poscic"
link <- "https://github.com/PlantMetals"
version <- "1.00"
date <- "September 2026"

# --------------------------------------------------------------------------
# 2. Helper function for correct input
# --------------------------------------------------------------------------

get_numeric <- function(prompt, condition = NULL, error_msg = "Invalid input. Please try again.") {
  repeat {
    val <- as.numeric(readline(prompt))
    if (length(val) == 0 || is.na(val)) {
      cat(error_msg, "\n")
      next
    }
    if (!is.null(condition) && !condition(val)) {
      cat(error_msg, "\n")
      next
    }
    return(val)
  }
}

# --------------------------------------------------------------------------
# 3. Table formatting functions
# --------------------------------------------------------------------------

table_width <- 77
inside_width <- table_width - 4

format_blank_line <- function() {
  sprintf("│ %-*s │", inside_width, "")
}

format_separator <- function() {
  paste0("├", strrep("─", inside_width + 2), "┤")
}

format_title_line <- function(title) {
  sprintf("│ %-*s │", inside_width, title)
}

format_table_line <- function(label, value) {
  left  <- paste0(label, ":")
  right <- as.character(value)
  sprintf("│ %-*s %15s │", inside_width - 16, left, right)
}

# --------------------------------------------------------------------------
# 4. Function to print the summary table
# --------------------------------------------------------------------------

print_summary_table <- function(mode, n, final_result, cohens_f,
                                v1, a, delta, ms_within, alpha,
                                version, date, target_power = NULL,
                                method = "delta_ms", user_cohens_f = NULL) {
  
  cat("┌──[ ANOVA POWER ANALYSIS ]────────────────────[ Version ", version, " - ", date, " ]──┐\n", sep = "")
  cat(format_blank_line(), "\n", sep = "")
  cat(format_title_line("USER INPUTS"), "\n", sep = "")
  
  cat(format_table_line("Number of groups (a)", a), "\n", sep = "")
  
  if (method == "delta_ms") {
    cat(format_table_line("delta (min. difference between extreme means)", delta), "\n", sep = "")
    cat(format_table_line("MS within", ms_within), "\n", sep = "")
    cat(format_table_line("Cohen's f (computed from delta/MS)", round(cohens_f, 8)), "\n", sep = "")
  } else { # cohens_f
    cat(format_table_line("Cohen's f (user input)", user_cohens_f), "\n", sep = "")
    # No delta/MS lines
  }
  
  cat(format_table_line("alpha (Type I error)", alpha), "\n", sep = "")
  
  if (mode == "required") {
    cat(format_table_line("Target power (1 - beta)", target_power), "\n", sep = "")
  } else {
    cat(format_table_line("Fixed n per group", n), "\n", sep = "")
  }
  
  cat(format_separator(), "\n", sep = "")
  cat(format_blank_line(), "\n", sep = "")
  
  # ---- Output summary ----
  if (mode == "required") {
    cat(format_title_line("OUTPUT SUMMARY (required sample size)"), "\n", sep = "")
    cat(format_table_line("Required independent n per group", n), "\n", sep = "")
    cat(format_table_line("Total replicates (n * a)", n * a), "\n", sep = "")
  } else {
    cat(format_title_line("OUTPUT SUMMARY (fixed sample size)"), "\n", sep = "")
    cat(format_table_line("Independent replicates per group", n), "\n", sep = "")
  }
  
  cat(format_separator(), "\n", sep = "")
  cat(format_blank_line(), "\n", sep = "")
  
  # ---- Power details ----
  cat(format_title_line("POWER DETAILS"), "\n", sep = "")
  cat(format_table_line("v1 = a - 1", v1), "\n", sep = "")
  cat(format_table_line("v2 = a * (n - 1)", a * (n - 1)), "\n", sep = "")
  cat(format_table_line("Achieved power (1 - beta)", round(final_result$power, 8)), "\n", sep = "")
  cat(format_table_line("Non-centrality parameter phi", round(final_result$phi, 8)), "\n", sep = "")
  cat(format_table_line("Critical F", round(final_result$fcrit, 8)), "\n", sep = "")
  
  cat("└", strrep("─", inside_width + 2), "┘\n", sep = "")
}

# --------------------------------------------------------------------------
# 5. Power calculation functions
# --------------------------------------------------------------------------

# Using delta and MS_within
compute_power <- function(n, a, delta, ms_within, alpha, v1) {
  v2 <- a * (n - 1)
  phi <- sqrt((n * delta^2) / (2 * a * ms_within))
  lambda <- a * phi^2
  fcrit <- qf(1 - alpha, v1, v2)
  power <- 1 - pf(fcrit, v1, v2, lambda)
  list(phi = phi, lambda = lambda, power = power, v1 = v1, v2 = v2, fcrit = fcrit)
}

# Using Cohen's f directly
compute_power_cohens <- function(n, a, cohens_f, alpha, v1) {
  v2 <- a * (n - 1)
  phi <- sqrt(n) * cohens_f
  lambda <- a * phi^2
  fcrit <- qf(1 - alpha, v1, v2)
  power <- 1 - pf(fcrit, v1, v2, lambda)
  list(phi = phi, lambda = lambda, power = power, v1 = v1, v2 = v2, fcrit = fcrit)
}

# --------------------------------------------------------------------------
# 6. Required sample size search functions (with progress indicator)
# --------------------------------------------------------------------------

find_required_n <- function(a, delta, ms_within, alpha, target_power, v1, initial_max_n = 100) {
  max_n <- initial_max_n
  n <- 2 # start of the search
  repeat {
    cat("  Searching n from", n, "to", max_n, "...\n")   # progress indicator
    while (n <= max_n) {
      result <- compute_power(n, a, delta, ms_within, alpha, v1)
      if (!is.na(result$power) && result$power >= target_power) return(n)
      n <- n + 1
    }
    # If we reach here, no n in [previous_range] met the target so we double max_n and continue the outer loop
    max_n <- max_n * 2
    if (max_n > 100000) return(NA)
  }
}

find_required_n_cohens <- function(a, cohens_f, alpha, target_power, v1, initial_max_n = 100) {
  max_n <- initial_max_n
  n <- 2 # start of the search
  repeat {
    cat("  Searching n from", n, "to", max_n, "...\n")   # progress indicator
    while (n <= max_n) {
      result <- compute_power_cohens(n, a, cohens_f, alpha, v1)
      if (!is.na(result$power) && result$power >= target_power) return(n)
      n <- n + 1
    }
    # If we reach here, no n in [previous_range] met the target so we double max_n and continue the outer loop
    max_n <- max_n * 2
    if (max_n > 100000) return(NA)
  }
}

# --------------------------------------------------------------------------
# 7. Function to draw the final power plot (used for display and saving)
# --------------------------------------------------------------------------

draw_final_power_plot <- function(v1, target_power, chart_choice,
                                  phi_grid, final_result, required_n,
                                  x_start, y_floor_trans, 
                                  y_transformed_ticks, y_ticks,
                                  a, method, cohens_f_input) {
  
  if (chart_choice == 1) {
    # Pearson & Hartley inverted log scale
    plot(NULL, NULL,
         xlim = c(1.0, 4.5), ylim = c(y_floor_trans, max(y_transformed_ticks)),
         yaxt = "n",
         xlab = expression(paste("Noncentrality index (", phi, ")")),
         ylab = "Probability of correctly rejecting the null hypothesis (Power)",
         main = paste0("Pearson & Hartley (1951) power curves (v1 = ", v1, ")"),
         xaxs = "i", yaxs = "i")
    
    abline(v = seq(1.0, 4.5, by = 0.5), col = "gray95")
    abline(h = y_transformed_ticks, col = "gray95")
    axis(2, at = y_transformed_ticks, labels = sprintf("%.2f", y_ticks),
         las = 2, cex.axis = 0.7, hadj = 0.9)
    
    target_y_trans <- -log10(1 - target_power)
    segments(x0 = 1.0, y0 = target_y_trans, x1 = 4.5, y1 = target_y_trans,
             col = "red", lty = 2, lwd = 1.5)
  } else {
    # Standard linear grid (log y-axis)
    plot(NULL, NULL,
         xlim = c(0.0, 4.5), ylim = c(0.10, 1), log = "y",
         xlab = expression(paste("Noncentrality index (", phi, ")")),
         ylab = "Probability of correctly rejecting the null hypothesis (Power)",
         main = paste0("Type-II power curves (v1 = ", v1, ")"),
         xaxs = "i")
    
    abline(v = seq(0.0, 4.5, by = 0.5), col = "gray95")
    abline(h = axTicks(2), col = "gray95")
    segments(x0 = 0.0, y0 = target_power, x1 = 4.5, y1 = target_power,
             col = "red", lty = 2, lwd = 1.5)
  }
  
  # ---- Draw the final (optimal) curve ----
  # Compute the power curve for the final n (using final_result$v2)
  # We need to recompute the curve for all phi values using the same v2
  # The curve is based on the final n's degrees of freedom.
  lambda_grid <- a * phi_grid^2
  fcrit <- qf(1 - alpha, v1, final_result$v2)
  curve_power <- 1 - pf(fcrit, v1, final_result$v2, lambda_grid)
  
  if (chart_choice == 1) {
    curve_plot <- -log10(1 - pmin(pmax(curve_power, 0.05), 0.99))
    point_y <- -log10(1 - final_result$power)
    y_floor <- y_floor_trans
  } else {
    curve_plot <- pmin(pmax(curve_power, 0.10), 1)
    point_y <- final_result$power
    y_floor <- 0.10
  }
  
  # Draw final curve in bold blue
  lines(phi_grid, curve_plot, col = "#1f77b4", lwd = 3)
  
  # Final annotations
  segments(x0 = final_result$phi, y0 = y_floor, x1 = final_result$phi, y1 = point_y,
           col = "darkgreen", lty = 4, lwd = 1.8)
  segments(x0 = x_start, y0 = point_y, x1 = final_result$phi, y1 = point_y,
           col = "darkgreen", lty = 4, lwd = 1.8)
  points(final_result$phi, point_y, col = "darkgreen", pch = 19, cex = 1.5)
  
  text_y_pos <- if (chart_choice == 1) (point_y - 0.15) else (point_y * 0.85)
  text(x = final_result$phi + 0.05, y = text_y_pos,
       labels = paste0("Optimal:\nn = ", required_n, ", v2 = ", final_result$v2,
                       "\nφ = ", round(final_result$phi, 4),
                       "\nPower = ", round(final_result$power, 4)),
       adj = 0, col = "darkgreen", font = 2, cex = 0.8)
  
  legend("bottomright",
         legend = c("Minimum required sample curve", "Target power threshold"),
         col = c("#1f77b4", "red"),
         lty = c(1, 2),
         lwd = c(3, 1.5),
         bty = "n", cex = 0.8)
}

# --------------------------------------------------------------------------
# 8. Batch processing function
# --------------------------------------------------------------------------

run_batch <- function() {
  cat("\n--- BATCH MODE ---\n")
  cat("Input CSV file must have columns:\n")
  cat("  mode, a, alpha, method, delta, ms_within, cohens_f, target_power, fixed_n\n")
  cat("  Please refer to the instructions\n\n")
  
  input_file <- readline(prompt = "Enter path and file name (with .csv) to input CSV file: ")
  if (!file.exists(input_file)) {
    cat("File not found.\n")
    return()
  }
  
  output_file <- readline(prompt = "Enter path and file name (with .csv) for output CSV file: ")
  if (output_file == "") {
    cat("No output file specified. Aborting.\n")
    return()
  }
  
  # Read CSV with stringsAsFactors = FALSE to avoid factors
  batch_data <- tryCatch({
    read.csv(input_file, stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("Error reading CSV:", e$message, "\n")
    return(NULL)
  })
  if (is.null(batch_data)) return()
  
  # Required columns
  required_cols <- c("mode", "a", "alpha", "method", "delta", "ms_within", 
                     "cohens_f", "target_power", "fixed_n")
  missing <- setdiff(required_cols, names(batch_data))
  if (length(missing) > 0) {
    cat("Missing columns:", paste(missing, collapse=", "), "\n")
    return()
  }
  
  # Define output columns (input columns + computed ones)
  out_cols <- c(required_cols, "n", "v1", "v2", "cohens_f_used", "phi", 
                "lambda", "fcrit", "power", "status")
  
  # Pre-allocate a list to store results (more efficient than rbind in loop)
  results_list <- vector("list", nrow(batch_data))
  
  cat("Processing", nrow(batch_data), "rows...\n")
  
  for (i in 1:nrow(batch_data)) {
    row <- batch_data[i, ]
    status <- "OK"
    n_out <- NA
    v1 <- NA
    v2 <- NA
    cohens_f_used <- NA
    phi <- NA
    lambda <- NA
    fcrit <- NA
    power <- NA
    
    # --- Validation ---
    mode <- trimws(tolower(as.character(row$mode)))
    if (!mode %in% c("required", "fixed")) {
      status <- paste("Invalid mode:", row$mode)
      # Build row with NA for computed fields
      results_list[[i]] <- list(
        mode = row$mode, a = row$a, alpha = row$alpha, method = row$method,
        delta = row$delta, ms_within = row$ms_within, cohens_f = row$cohens_f,
        target_power = row$target_power, fixed_n = row$fixed_n,
        n = NA, v1 = NA, v2 = NA, cohens_f_used = NA,
        phi = NA, lambda = NA, fcrit = NA, power = NA, status = status
      )
      next
    }
    
    a <- as.numeric(row$a)
    if (is.na(a) || a < 2 || a != round(a)) {
      status <- "a must be integer >= 2"
      results_list[[i]] <- list(
        mode = row$mode, a = row$a, alpha = row$alpha, method = row$method,
        delta = row$delta, ms_within = row$ms_within, cohens_f = row$cohens_f,
        target_power = row$target_power, fixed_n = row$fixed_n,
        n = NA, v1 = NA, v2 = NA, cohens_f_used = NA,
        phi = NA, lambda = NA, fcrit = NA, power = NA, status = status
      )
      next
    }
    a <- as.integer(a)
    
    alpha <- as.numeric(row$alpha)
    if (is.na(alpha) || alpha <= 0 || alpha >= 1) {
      status <- "alpha must be between 0 and 1"
      results_list[[i]] <- list(
        mode = row$mode, a = row$a, alpha = row$alpha, method = row$method,
        delta = row$delta, ms_within = row$ms_within, cohens_f = row$cohens_f,
        target_power = row$target_power, fixed_n = row$fixed_n,
        n = NA, v1 = NA, v2 = NA, cohens_f_used = NA,
        phi = NA, lambda = NA, fcrit = NA, power = NA, status = status
      )
      next
    }
    
    method <- trimws(tolower(as.character(row$method)))
    if (!method %in% c("delta_ms", "cohens_f")) {
      status <- paste("Invalid method:", row$method)
      results_list[[i]] <- list(
        mode = row$mode, a = row$a, alpha = row$alpha, method = row$method,
        delta = row$delta, ms_within = row$ms_within, cohens_f = row$cohens_f,
        target_power = row$target_power, fixed_n = row$fixed_n,
        n = NA, v1 = NA, v2 = NA, cohens_f_used = NA,
        phi = NA, lambda = NA, fcrit = NA, power = NA, status = status
      )
      next
    }
    
    if (method == "delta_ms") {
      delta <- as.numeric(row$delta)
      ms_within <- as.numeric(row$ms_within)
      if (is.na(delta) || delta <= 0) {
        status <- "delta must be positive"
        results_list[[i]] <- list(
          mode = row$mode, a = row$a, alpha = row$alpha, method = row$method,
          delta = row$delta, ms_within = row$ms_within, cohens_f = row$cohens_f,
          target_power = row$target_power, fixed_n = row$fixed_n,
          n = NA, v1 = NA, v2 = NA, cohens_f_used = NA,
          phi = NA, lambda = NA, fcrit = NA, power = NA, status = status
        )
        next
      }
      if (is.na(ms_within) || ms_within <= 0) {
        status <- "ms_within must be positive"
        results_list[[i]] <- list(
          mode = row$mode, a = row$a, alpha = row$alpha, method = row$method,
          delta = row$delta, ms_within = row$ms_within, cohens_f = row$cohens_f,
          target_power = row$target_power, fixed_n = row$fixed_n,
          n = NA, v1 = NA, v2 = NA, cohens_f_used = NA,
          phi = NA, lambda = NA, fcrit = NA, power = NA, status = status
        )
        next
      }
      cohens_f_used <- sqrt(delta^2 / (2 * a * ms_within))
    } else { # cohens_f
      cohens_f_input <- as.numeric(row$cohens_f)
      if (is.na(cohens_f_input) || cohens_f_input <= 0) {
        status <- "cohens_f must be positive"
        results_list[[i]] <- list(
          mode = row$mode, a = row$a, alpha = row$alpha, method = row$method,
          delta = row$delta, ms_within = row$ms_within, cohens_f = row$cohens_f,
          target_power = row$target_power, fixed_n = row$fixed_n,
          n = NA, v1 = NA, v2 = NA, cohens_f_used = NA,
          phi = NA, lambda = NA, fcrit = NA, power = NA, status = status
        )
        next
      }
      cohens_f_used <- cohens_f_input
      delta <- NA
      ms_within <- NA
    }
    
    v1 <- a - 1
    
    # Mode-specific computations
    if (mode == "required") {
      target_power <- as.numeric(row$target_power)
      if (is.na(target_power) || target_power <= 0 || target_power >= 1) {
        status <- "target_power must be between 0 and 1"
        results_list[[i]] <- list(
          mode = row$mode, a = row$a, alpha = row$alpha, method = row$method,
          delta = row$delta, ms_within = row$ms_within, cohens_f = row$cohens_f,
          target_power = row$target_power, fixed_n = row$fixed_n,
          n = NA, v1 = NA, v2 = NA, cohens_f_used = NA,
          phi = NA, lambda = NA, fcrit = NA, power = NA, status = status
        )
        next
      }
      if (method == "delta_ms") {
        n_out <- find_required_n(a, delta, ms_within, alpha, target_power, v1)
      } else {
        n_out <- find_required_n_cohens(a, cohens_f_used, alpha, target_power, v1)
      }
      if (is.na(n_out)) {
        status <- "Required n not found (max_n exceeded)"
        results_list[[i]] <- list(
          mode = row$mode, a = row$a, alpha = row$alpha, method = row$method,
          delta = row$delta, ms_within = row$ms_within, cohens_f = row$cohens_f,
          target_power = row$target_power, fixed_n = row$fixed_n,
          n = NA, v1 = NA, v2 = NA, cohens_f_used = NA,
          phi = NA, lambda = NA, fcrit = NA, power = NA, status = status
        )
        next
      }
      if (method == "delta_ms") {
        res <- compute_power(n_out, a, delta, ms_within, alpha, v1)
      } else {
        res <- compute_power_cohens(n_out, a, cohens_f_used, alpha, v1)
      }
      phi <- res$phi
      lambda <- res$lambda
      power <- res$power
      fcrit <- res$fcrit
      v2 <- res$v2
      # fixed_n is not used
    } else { # fixed
      fixed_n <- as.numeric(row$fixed_n)
      if (is.na(fixed_n) || fixed_n < 2 || fixed_n != round(fixed_n)) {
        status <- "fixed_n must be integer >= 2"
        results_list[[i]] <- list(
          mode = row$mode, a = row$a, alpha = row$alpha, method = row$method,
          delta = row$delta, ms_within = row$ms_within, cohens_f = row$cohens_f,
          target_power = row$target_power, fixed_n = row$fixed_n,
          n = NA, v1 = NA, v2 = NA, cohens_f_used = NA,
          phi = NA, lambda = NA, fcrit = NA, power = NA, status = status
        )
        next
      }
      fixed_n <- as.integer(fixed_n)
      n_out <- fixed_n
      if (method == "delta_ms") {
        res <- compute_power(fixed_n, a, delta, ms_within, alpha, v1)
      } else {
        res <- compute_power_cohens(fixed_n, a, cohens_f_used, alpha, v1)
      }
      phi <- res$phi
      lambda <- res$lambda
      power <- res$power
      fcrit <- res$fcrit
      v2 <- res$v2
      # target_power not used
    }
    
    # Build the successful row
    results_list[[i]] <- list(
      mode = row$mode, a = row$a, alpha = row$alpha, method = row$method,
      delta = row$delta, ms_within = row$ms_within, cohens_f = row$cohens_f,
      target_power = row$target_power, fixed_n = row$fixed_n,
      n = n_out, v1 = v1, v2 = v2, cohens_f_used = cohens_f_used,
      phi = phi, lambda = lambda, fcrit = fcrit, power = power, status = status
    )
  }
  
  # Convert list of lists to data frame
  results <- do.call(rbind, lapply(results_list, as.data.frame, stringsAsFactors = FALSE))
  
  # Write output
  tryCatch({
    write.csv(results, file = output_file, row.names = FALSE)
    cat("Batch processing completed. Results written to", output_file, "\n")
  }, error = function(e) {
    cat("Error writing output file:", e$message, "\n")
  })
  
  # Show summary
  cat("\nSummary:\n")
  cat("  Total rows processed:", nrow(results), "\n")
  cat("  Successful:", sum(results$status == "OK"), "\n")
  cat("  Errors:   ", sum(results$status != "OK"), "\n")
  if (sum(results$status != "OK") > 0) {
    cat("  Error details:\n")
    print(table(results$status[results$status != "OK"]))
  }
}

# --------------------------------------------------------------------------
# 9. Main execution loop
# --------------------------------------------------------------------------

repeat {
  cat(
    "┌──[ ANOVA POWER ANALYSIS ]────────────────────[ Version ", version, " - ", date, " ]──┐\n",
    "│                                                                           │\n",
    "│   Finding the required sample size and estimating the power on an ANOVA   │\n",
    "│                                                                           │\n",
    "│                                                     by Dr. Filip Poscic   │\n",
    "├────────────────────────────────────────[ ", link, " ]─┤\n",
    "│                                                                           │\n",
    "│ Select your choice (1, 2, 3, 4, or 5) and press ENTER                     │\n",
    "│                                                                           │\n",
    "│  1. Determine required n per group to achieve target power                │\n",
    "│  2. Calculate statistical power for a fixed sample size                   │\n",
    "│  3. Run batch analysis from CSV file (no graphs)                          │\n",
    "│  4. About the statistical methods, instructions and references            │\n",
    "│  5. Exit the program                                                      │\n",
    "└───────────────────────────────────────────────────────────────────────────┘\n",
    sep = ""
  )
  
  task_choice <- get_numeric(
    prompt = "",
    condition = function(x) x %in% c(1, 2, 3, 4, 5),
    error_msg = "Please enter 1, 2, 3, 4, or 5."
  )
  
  if (task_choice == 5) {
    cat("\nExiting program. Goodbye!\n")
    break
  }
  
  if (task_choice == 4) {
    cat(
      "┌──[ ANOVA POWER ANALYSIS ]────────────────────────[ Version", version, "-", date, "]──┐\n",
      "│                                                                           │\n",
      "│ INTRODUCTION TO POWER ANALYSIS                                            │\n",
      "│                                                                           │\n",
      "│  What is it? Power analysis is a statistical planning tool that           │\n",
      "│  determines the minimum sample size needed to have a good chance          │\n",
      "│  ('power') of detecting a true difference between groups, if one exists.  │\n",
      "│                                                                           │\n",
      "│  Why do it?                                                               │\n",
      "│   • Too small n means high risk of missing a real effect (false           │\n",
      "│     negative). Your experiment would be wasted and the results            │\n",
      "│     uninformative. This problem is, unfortunately, overrepresented in     │\n",
      "│     virtually every branch of modern-day science.                         │\n",
      "│   • Too large n means waste of resources, time, or (in animal/clinical    │\n",
      "│     studies) excessive use of subjects.                                   │\n",
      "│                                                                           │\n",
      "│  It ensures your results are reliable and defensible. A low-power study   │\n",
      "│  that finds 'no significant difference' is inconclusive – it does NOT     │\n",
      "│  prove that no effect exists, only that your study was possibly too small │\n",
      "│  to detect it. Power analysis protects you from this trap. In any case it │\n",
      "│  is good practice to take advantage also of, for example, equivalence     │\n",
      "│  testing (not part of this script) to confidently conclude there is no    │\n",
      "│  meaningful effect.                                                       │\n",
      "├───────────────────────────────────────────────────────────────────────────┤\n",
      "│                                                                           │\n",
      "│ STATISTICAL VOCABULARY                                                    │\n",
      "│                                                                           │\n",
      "│  n       : Sample size per group (number of independent replicates).      │\n",
      "│  alpha   : Significance level (Type I error rate). Probability of         │\n",
      "│            falsely rejecting H_0 when it is true (e.g., 0.05).            │\n",
      "│  beta    : Type II error rate. Probability of failing to reject a         │\n",
      "│            false H_0.                                                     │\n",
      "│  power   : Power = 1 - beta. Probability of correctly rejecting H_0       │\n",
      "│            (detecting a true effect). Target is often 0.80 or 0.90.       │\n",
      "│  v1      : Numerator degrees of freedom = a - 1 (between groups).         │\n",
      "│            Together with v2, they define the shape of the F-distribution. │\n",
      "│  v2      : Denominator degrees of freedom = a * (n - 1) (within groups).  │\n",
      "│            Larger v2 increases power for a given effect size.             │\n",
      "│  phi     : Non-centrality parameter [sqrt(lambda/a)]. Measures the        │\n",
      "│            overall effect magnitude, combining sample size and effect     │\n",
      "│            size.                                                          │\n",
      "│  lambda  : Non-centrality parameter = a * phi^2. Used in the non-central  │\n",
      "│            F-distribution to compute power.                               │\n",
      "│  F_crit  : Critical F-value from the central F-distribution at alpha.     │\n",
      "│            The test rejects H_0 if observed F > F_crit.                   │\n",
      "│                                                                           │\n",
      "│  Graphs – what they show and how to read them:                            │\n",
      "│   The graphs display power (y-axis) as a function of non-centrality phi   │\n",
      "│   (x-axis). As phi increases (larger sample size or stronger effect),     │\n",
      "│   power rises. The bold blue curve represents the power for your          │\n",
      "│   calculated sample size (n). The red dashed line marks your target       │\n",
      "│   power. The green annotation shows the exact phi and power at the        │\n",
      "│   optimal n. In the animation mode, grey curves show how power changes    │\n",
      "│   as n increases step by step, helping you visualise the search process.  │\n",
      "│                                                                           │\n",
      "│  Two graph styles:                                                        │\n",
      "│   1. Pearson & Hartley (1951) – a classic transformation that uses        │\n",
      "│      -log10(1-power) on the y-axis. This linearises the curves and        │\n",
      "│      was traditionally used to read power from printed charts.            │\n",
      "│   2. Standard log-scale – power is plotted on a logarithmic axis from     │\n",
      "│      0.10 to 1.00. This makes high-power values (e.g., 0.90 vs 0.95)      │\n",
      "│      easier to distinguish.                                               │\n",
      "│                                                                           │\n",
      "│  Cohen's f approximation                                                  │\n",
      "│   This program estimates Cohen's f (a measure of effect size) from the    │\n",
      "│   user-specified largest expected difference between treatment means (δ), │\n",
      "│   number of groups in one-way ANOVA (a), and the within-group mean square │\n",
      "│   (MS within) (also known as the error term in ANOVA):                    │\n",
      "│                                                                           │\n",
      "│           Cohen's f approximation = √[δ² / (2 · a · MS-within)]           │\n",
      "│                                                                           │\n",
      "│   This is the general approximate relationship given by Sokal & Rohlf     │\n",
      "│   (2012, Box 9.12). It does not assume a particular pattern of treatment  │\n",
      "│   means. When the means are in fact equally spaced, this approximation    │\n",
      "│   yields a conservative (larger) estimate of δ (i.e., it overestimates    │\n",
      "│   the true detectable difference and thus leads to a slightly larger      │\n",
      "│   required sample size than the exact equal-spacing formula).             │\n",
      "│                                                                           │\n",
      "│   The exact definition of Cohen's f is based on the pooled standard       │\n",
      "│   deviation (σ²) of all treatment means (μᵢ) around the grand mean (μ):   │\n",
      "│                                                                           │\n",
      "│                Cohen's f exact = √{Σ[pᵢ·(μᵢ − μ)²] / σ²}                  │\n",
      "│                                                                           │\n",
      "│   where pᵢ = nᵢ/N, nᵢ is the sample size of group i, and N is the total   │\n",
      "│   sample size across all groups. This is equal to:                        │\n",
      "│                                                                           │\n",
      "│                Cohen's f exact = √(SS-among / SS-within)                  │\n",
      "│                                                                           │\n",
      "│   Thus, Cohen's f approximation used in its indirect calculation, by      │\n",
      "│   providing δ and MS-witihin, is a conservative rule of thumb for         │\n",
      "│   planning experiments when the exact pattern of means is unknown.        │\n",
      "├───────────────────────────────────────────────────────────────────────────┤\n",
      "│                                                                           │\n",
      "│ INSTRUCTIONS FOR BATCH MODE                                               │\n",
      "│                                                                           │\n",
      "│  Following is the input file preparation explanation. To use batch mode   │\n",
      "│  (option 3), prepare a CSV file with these columns:                       │\n",
      "│                                                                           │\n",
      "│                mode, a, alpha, method, delta, ms_within,                  │\n",
      "│                cohens_f, target_power, fixed_n                            │\n",
      "│                                                                           │\n",
      "│   • mode:         'required' or 'fixed'                                   │\n",
      "│                (required means you search for minimum n at desired power, │\n",
      "│                 fixed means n is already fixed and you search for power)  │\n",
      "│   • a:            integer >= 2                                            │\n",
      "│   • alpha:        0 < alpha < 1                                           │\n",
      "│   • method:       'delta_ms' or 'cohens_f'                                │\n",
      "│   • delta:        positive (if method=delta_ms), else NA or empty (,,)    │\n",
      "│   • ms_within:    positive (if method=delta_ms), else NA or empty (,,)    │\n",
      "│   • cohens_f:     positive (if method=cohens_f), else NA or empty (,,)    │\n",
      "│   • target_power: 0 < power < 1 (required only if mode='required')        │\n",
      "│   • fixed_n:      integer >= 2 (required only if mode='fixed')            │\n",
      "│                                                                           │\n",
      "│  Example rows:                                                            │\n",
      "│                 required,3,0.05,delta_ms,2,1.5,,0.8,                      │\n",
      "│                 fixed,4,0.01,cohens_f,,,0.35,,15                          │\n",
      "│                                                                           │\n",
      "│  The output CSV will include all input columns plus computed results:     │\n",
      "│  n, v1, v2, cohens_f_used, phi, lambda, fcrit, power, and a status column │\n",
      "│  indicating success or error for each row. If error check your inputs.    │\n",
      "├───────────────────────────────────────────────────────────────────────────┤\n",
      "│                                                                           │\n",
      "│ REFERENCES                                                                │\n",
      "│                                                                           │\n",
      "│  Pearson, E. S. & Hartley H. O. Charts of the power function for analysis │\n",
      "│    of variance tests, derived from the non-central F-distribution.        │\n",
      "│    Biometrika 38, 112-130 (1951).                                         │\n",
      "│                                                                           │\n",
      "│  Sokal, R. R. & Rohlf, F. J. Biometry: The principles and practice of     │\n",
      "│   statistics in biological research, 4th ed. W.H. Freeman and Company,    │\n",
      "│   New York, NY, USA (2012).                                               │\n",
      "├───────────────────────────────────────────────────────────────────────────┤\n",
      "│                                                                           │\n",
      "│ ABOUT THE AUTHOR                                                          │\n",
      "│                                                                           │\n",
      "│  I am a plant scientist and educator with a focus on abiotic stress       │\n",
      "│  responses, phytoremediation, biometry, and sustainable agriculture.      │\n",
      "│  My research investigates how plants adapt to metal toxicity, nutrient    │\n",
      "│  deficiencies, drought, salinity, and climate extremes – using these      │\n",
      "│  insights to develop eco-friendly strategies for soil health and crop     │\n",
      "│  resilience. The ultimate goal is to support food security and            │\n",
      "│  environmental restoration.                                               │\n",
      "│                                                                           │\n",
      "│  This software is completely free and open source. You can find it on     │\n",
      "│  my GitHub channel: https://github.com/PlantMetals                        │\n",
      "│  For questions or feedback, please reach me at:                           │\n",
      "│  Filip.Poscic@gmail.com  or  BlueSky @filipposcic.bsky.social             │\n",
      "└───────────────────────────────────────────────────────────────────────────┘\n",
      sep = ""
    )
    next
  }
  
  if (task_choice == 3) {
    run_batch()
    cat("\nPress ENTER to return to the main menu...")
    readline()
    next
  }
  
  # ---- Common parameters: a and alpha ----
  a <- get_numeric(
    prompt = " Insert the number of separate treatments (a) you are planning to have:\n",
    condition = function(x) x >= 2 && x %% 1 == 0,
    error_msg = "a must be an integer >= 2."
  )
  
  alpha <- get_numeric(
    prompt = " Insert the Type-I error rate threshold α (significance level, e.g., 0.05):\n",
    condition = function(x) x > 0 && x < 1,
    error_msg = "Alpha must be between 0 and 1."
  )
  
  # ---- Ask for effect size method ----
  cat(" Do you want to specify Cohen's f directly instead of delta and MS within? (y/n):\n")
  use_f <- tolower(trimws(readline(prompt = "")))
  while (!use_f %in% c("y", "n")) {
    cat("Please enter 'y' or 'n'.\n")
    use_f <- tolower(trimws(readline(prompt = "")))
  }
  
  method <- if (use_f == "y") "cohens_f" else "delta_ms"
  
  # ---- Read effect size accordingly ----
  if (method == "cohens_f") {
    cohens_f_input <- get_numeric(
      prompt = " Insert the value of Cohen's f (effect size):\n",
      condition = function(x) x > 0,
      error_msg = "Cohen's f must be positive."
    )
    delta <- NA
    ms_within <- NA
    cohens_f_used <- cohens_f_input
  } else {
    delta <- get_numeric(
      prompt = " Insert the expected minimum difference (δ) between the most extreme means:\n",
      condition = function(x) x > 0,
      error_msg = "δ must be positive."
    )
    ms_within <- get_numeric(
      prompt = " Insert the estimated MS within (error) from preliminary experiments:\n",
      condition = function(x) x > 0,
      error_msg = "MS within must be positive."
    )
    cohens_f_used <- sqrt(delta^2 / (2 * a * ms_within))
  }
  
  v1 <- a - 1
  
  # ---- Task-specific ----
  if (task_choice == 1) {
    target_power <- get_numeric(
      prompt = " Insert the target statistical power (P) (e.g., 0.80 for 80%):\n",
      condition = function(x) x > 0 && x < 1,
      error_msg = "Power must be between 0 and 1."
    )
    
    # ---- Calculate required n based on method ----
    if (method == "delta_ms") {
      required_n <- find_required_n(a = a, delta = delta, ms_within = ms_within, alpha = alpha, target_power = target_power, v1 = v1)
      final_result <- compute_power(n = required_n, a = a, delta = delta, ms_within = ms_within, alpha = alpha, v1 = v1)
    } else {
      required_n <- find_required_n_cohens(a = a, cohens_f = cohens_f_input, alpha = alpha, target_power = target_power, v1 = v1)
      final_result <- compute_power_cohens(n = required_n, a = a, cohens_f = cohens_f_input, alpha = alpha, v1 = v1)
    }
    
    if (is.na(required_n)) {
      cat("No valid sample size found. Please adjust parameters.\n")
      next
    }
    
    # ---- Plot options (only for task 1) ----
    want_plot <- "n"
    chart_choice <- 1
    animation_mode <- "y"   # default to full animation if plot is chosen
    
    cat(" Generate graphical output? (y/n):\n")
    repeat {
      want_plot <- tolower(trimws(readline(prompt = "")))
      if (want_plot %in% c("y", "n")) break
      cat("Please enter 'y' or 'n'.\n")
    }
    
    if (want_plot == "y") {
      cat(" Do you want the full animation (y) or just the final plot (n)? (y/n):\n")
      repeat {
        animation_mode <- tolower(trimws(readline(prompt = "")))
        if (animation_mode %in% c("y", "n")) break
        cat("Please enter 'y' or 'n'.\n")
      }
      
      cat(
        "┌───────────────────────────────────────────────────────────────────────────┐\n",
        "│ Select your preferred chart format (1 or 2) and press ENTER               │\n",
        "│                                                                           │\n",
        "│  1. Classic Pearson & Hartley (1951) chart layout                         │\n",
        "│  2. Power chart in a standard Log10 grid layout                           │\n",
        "└───────────────────────────────────────────────────────────────────────────┘\n",
        sep = ""
      )
      chart_choice <- get_numeric(
        prompt = "",
        condition = function(x) x %in% c(1, 2),
        error_msg = "Please enter 1 or 2."
      )
    }
    
    # ---- Graphical output ----
    if (want_plot == "y") {
      
      # ---- Base plot setup ----
      x_start <- if (chart_choice == 2) 0.0 else 1.0
      phi_grid <- seq(x_start, 4.5, length.out = 500)
      
      if (chart_choice == 1) {
        y_ticks <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.92, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99)
        y_transformed_ticks <- -log10(1 - y_ticks)
        y_floor_trans <- -log10(1 - 0.05)
      } else {
        y_ticks <- NULL   # not used in standard layout
        y_transformed_ticks <- NULL
        y_floor_trans <- NULL
      }
      
      # ---- Draw the final plot (function defined above) ----
      # We'll call draw_final_power_plot to create the plot on the current device.
      # For animation, we need to draw the base plot first and then add curves iteratively.
      
      if (animation_mode == "n") {
        # Just draw the final plot directly (no animation)
        draw_final_power_plot(v1, target_power, chart_choice,
                              phi_grid, final_result, required_n,
                              x_start, y_floor_trans, 
                              y_transformed_ticks, y_ticks,
                              a, method, cohens_f_input)
      } else {
        # ---- Full animation ----
        if (chart_choice == 1) {
          plot(NULL, NULL,
               xlim = c(1.0, 4.5), ylim = c(y_floor_trans, max(y_transformed_ticks)),
               yaxt = "n",
               xlab = expression(paste("Noncentrality index (", phi, ")")),
               ylab = "Probability of correctly rejecting the null hypothesis (Power)",
               main = paste0("Pearson & Hartley (1951) power curves (v1 = ", v1, ")"),
               xaxs = "i", yaxs = "i")
          
          abline(v = seq(1.0, 4.5, by = 0.5), col = "gray95")
          abline(h = y_transformed_ticks, col = "gray95")
          axis(2, at = y_transformed_ticks, labels = sprintf("%.2f", y_ticks),
               las = 2, cex.axis = 0.7, hadj = 0.9)
          
          target_y_trans <- -log10(1 - target_power)
          segments(x0 = 1.0, y0 = target_y_trans, x1 = 4.5, y1 = target_y_trans,
                   col = "red", lty = 2, lwd = 1.5)
        } else {
          plot(NULL, NULL,
               xlim = c(0.0, 4.5), ylim = c(0.10, 1), log = "y",
               xlab = expression(paste("Noncentrality index (", phi, ")")),
               ylab = "Probability of correctly rejecting the null hypothesis (Power)",
               main = paste0("Type-II power curves (v1 = ", v1, ")"),
               xaxs = "i")
          
          abline(v = seq(0.0, 4.5, by = 0.5), col = "gray95")
          abline(h = axTicks(2), col = "gray95")
          segments(x0 = 0.0, y0 = target_power, x1 = 4.5, y1 = target_power,
                   col = "red", lty = 2, lwd = 1.5)
        }
        
        # Animation loop
        animation_speed <- 0.15
        for (n in 2:required_n) {
          if (method == "delta_ms") {
            current <- compute_power(n, a, delta, ms_within, alpha, v1)
          } else {
            current <- compute_power_cohens(n, a, cohens_f_input, alpha, v1)
          }
          
          # Vectorized power curve for this n
          lambda_grid <- a * phi_grid^2
          fcrit <- qf(1 - alpha, v1, current$v2)
          curve_power <- 1 - pf(fcrit, v1, current$v2, lambda_grid)
          
          if (chart_choice == 1) {
            curve_plot <- -log10(1 - pmin(pmax(curve_power, 0.05), 0.99))
            point_y <- -log10(1 - current$power)
            y_floor <- y_floor_trans
          } else {
            curve_plot <- pmin(pmax(curve_power, 0.10), 1)
            point_y <- current$power
            y_floor <- 0.10
          }
          
          if (n < required_n) {
            lines(phi_grid, curve_plot, col = "gray75", lwd = 0.8)
            plot_snapshot <- recordPlot()
            segments(x0 = current$phi, y0 = y_floor, x1 = current$phi, y1 = point_y,
                     col = "gray60", lty = 3, lwd = 1.2)
            segments(x0 = x_start, y0 = point_y, x1 = current$phi, y1 = point_y,
                     col = "gray60", lty = 3, lwd = 1.2)
            points(current$phi, point_y, col = "gray40", pch = 1, cex = 0.8)
            text(x = current$phi + 0.05, y = point_y,
                 labels = paste0("n=", n, "\nv2=", current$v2),
                 col = "gray50", cex = 0.65, adj = 0)
            Sys.sleep(animation_speed)
            replayPlot(plot_snapshot)
          } else {
            # Final curve: we replicate the final drawing part here.
            lines(phi_grid, curve_plot, col = "#1f77b4", lwd = 3)
            segments(x0 = current$phi, y0 = y_floor, x1 = current$phi, y1 = point_y,
                     col = "darkgreen", lty = 4, lwd = 1.8)
            segments(x0 = x_start, y0 = point_y, x1 = current$phi, y1 = point_y,
                     col = "darkgreen", lty = 4, lwd = 1.8)
            points(current$phi, point_y, col = "darkgreen", pch = 19, cex = 1.5)
            
            text_y_pos <- if (chart_choice == 1) (point_y - 0.15) else (point_y * 0.85)
            text(x = current$phi + 0.05, y = text_y_pos,
                 labels = paste0("Optimal:\nn = ", n, ", v2 = ", current$v2,
                                 "\nφ = ", round(current$phi, 4),
                                 "\nPower = ", round(current$power, 4)),
                 adj = 0, col = "darkgreen", font = 2, cex = 0.8)
            
            legend("bottomright",
                   legend = c("Sub-optimal paths", "Minimum required sample curve", "Target power threshold"),
                   col = c("gray75", "#1f77b4", "red"),
                   lty = c(1, 1, 2),
                   lwd = c(0.8, 3, 1.5),
                   bty = "n", cex = 0.8)
          }
        }
      } # end of animation mode
      
      # ---- Save plot (if requested) ----
      # open new device, draw final plot, close
      Sys.sleep(0.5)
      cat("Save the final plot? (y/n):\n")
      save_plot <- tolower(trimws(readline(prompt = "")))
      
      if (save_plot == "y") {
        cat("Enter file name (without extension) or press Enter for default 'Power_plot'\n")
        fname <- trimws(readline(prompt = ""))
        if (fname == "") fname <- "Power_plot"
        
        cat(
          "┌───────────────────────────────────────────────────────────────────────────┐\n",
          "│ Choose format (1 or 2) and press ENTER                                    │\n",
          "│                                                                           │\n",
          "│  1. PNG                                                                   │\n",
          "│  2. PDF                                                                   │\n",
          "└───────────────────────────────────────────────────────────────────────────┘\n",
          sep = ""
        )
        fmt <- get_numeric(
          prompt = "",
          condition = function(x) x %in% c(1, 2),
          error_msg = "Please enter 1 or 2."
        )
        
        # open device explicitly, draw the final plot, then close
        if (fmt == 1) {
          png(filename = paste0(fname, ".png"), width = 800, height = 600)
          draw_final_power_plot(v1, target_power, chart_choice,
                                phi_grid, final_result, required_n,
                                x_start, y_floor_trans, 
                                y_transformed_ticks, y_ticks,
                                a, method, cohens_f_input)
          dev.off()
          cat("Plot saved as", paste0(fname, ".png"), "\n")
        } else {
          pdf(file = paste0(fname, ".pdf"), width = 8, height = 6)
          draw_final_power_plot(v1, target_power, chart_choice,
                                phi_grid, final_result, required_n,
                                x_start, y_floor_trans, 
                                y_transformed_ticks, y_ticks,
                                a, method, cohens_f_input)
          dev.off()
          cat("Plot saved as", paste0(fname, ".pdf"), "\n")
        }
      }
    } # end if (want_plot == "y")
    
    # ---- Print summary table ----
    print_summary_table(
      mode = "required",
      n = required_n,
      final_result = final_result,
      cohens_f = cohens_f_used,
      v1 = v1,
      a = a,
      delta = delta,
      ms_within = ms_within,
      alpha = alpha,
      version = version,
      date = date,
      target_power = target_power,
      method = method,
      user_cohens_f = if (method == "cohens_f") cohens_f_input else NULL
    )
    
  } else {  # task_choice == 2: fixed n
    fixed_n <- get_numeric(
      prompt = " Insert the fixed number of independent replicates per group (n):\n",
      condition = function(x) x >= 2 && x %% 1 == 0,
      error_msg = "n must be an integer >= 2."
    )
    
    if (method == "delta_ms") {
      final_result <- compute_power(n = fixed_n, a = a, delta = delta,
                                    ms_within = ms_within, alpha = alpha, v1 = v1)
    } else {
      final_result <- compute_power_cohens(n = fixed_n, a = a, cohens_f = cohens_f_input,
                                           alpha = alpha, v1 = v1)
    }
    
    print_summary_table(
      mode = "fixed",
      n = fixed_n,
      final_result = final_result,
      cohens_f = cohens_f_used,
      v1 = v1,
      a = a,
      delta = delta,
      ms_within = ms_within,
      alpha = alpha,
      version = version,
      date = date,
      target_power = NULL,
      method = method,
      user_cohens_f = if (method == "cohens_f") cohens_f_input else NULL
    )
  }
  
  # ---- Pause before looping back ----
  cat("\nPress ENTER to return to the main menu...")
  readline()
}