# Build a small DigiBog_Hydro configuration for testing the daily mean.
#
# 5 x 5 grid, 10 timesteps per day, 4 days. Small enough to run in about a
# second, but the water table must genuinely move within each day, otherwise
# a daily mean would equal the end-of-day snapshot and the test would prove
# nothing. Rainfall therefore differs from day to day.
#
# Files are written to tests/small_config/.

out_dir <- "tests/small_config"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

x_extent <- 5
y_extent <- 5
n_cells  <- x_extent * y_extent

daily_timesteps <- 10
n_days   <- 4
t_extent <- daily_timesteps * n_days

# --- column status -----------------------------------------------------
# Corners off, rest of the perimeter Dirichlet, interior active. Written in
# x-outer, y-inner order to match how the model reads it.
status <- matrix("on", nrow = x_extent, ncol = y_extent)
status[1, ] <- "diri"
status[x_extent, ] <- "diri"
status[, 1] <- "diri"
status[, y_extent] <- "diri"
status[1, 1] <- "off"
status[1, y_extent] <- "off"
status[x_extent, 1] <- "off"
status[x_extent, y_extent] <- "off"

# t() then as.vector gives x-outer, y-inner ordering
status_vec <- as.vector(t(status))

# --- per-column fields -------------------------------------------------
no_layers  <- rep(1, n_cells)      # model adds a ponding layer itself
thickness  <- rep(200.0, n_cells)  # cm
k_layer    <- rep(0.001, n_cells)  # cm/s
s_layer    <- rep(0.3, n_cells)    # drainable porosity
wt_bc      <- rep(100.0, n_cells)  # cm above base
baltitude  <- rep(0.0, n_cells)    # flat substrate
aet_status <- rep(0, n_cells)      # no gully columns

# --- climate -----------------------------------------------------------
# Rates in cm/s, held constant through each day. Varied between days so the
# water table rises on some days and falls on others.
rainfall <- c(2.4e-05, 0.0, 8.0e-06, 1.5e-05)
pet      <- c(3.7e-07, 3.7e-07, 3.7e-07, 3.7e-07)

# --- parameters (13 lines, output_stat last) ---------------------------
params <- c(t_extent,          # runtime in timesteps
            daily_timesteps,   # output interval in timesteps
            daily_timesteps,   # timesteps per day
            x_extent,
            y_extent,
            3,                 # z extent
            n_cells,           # steady columns
            100,               # spatial step, cm
            0,                 # ponding depth, cm
            0,                 # steady threshold, cm
            50,                # n_aet
            100,               # aet extinction depth, cm
            1)                 # output_stat

writeLines(as.character(params), file.path(out_dir, "hydro_parameters.txt"))
writeLines(status_vec, file.path(out_dir, "hydro_column_status.txt"))
writeLines(as.character(no_layers), file.path(out_dir, "hydro_no_layers.txt"))
writeLines(format(thickness, scientific = FALSE),
           file.path(out_dir, "hydro_thickness.txt"))
writeLines(format(k_layer, scientific = FALSE),
           file.path(out_dir, "hydro_k.txt"))
writeLines(format(s_layer, scientific = FALSE),
           file.path(out_dir, "hydro_s.txt"))
writeLines(format(wt_bc, scientific = FALSE),
           file.path(out_dir, "hydro_wt_bc_input.txt"))
writeLines(format(baltitude, scientific = FALSE),
           file.path(out_dir, "hydro_baltitude.txt"))
writeLines(as.character(aet_status),
           file.path(out_dir, "hydro_aet_status.txt"))
writeLines(formatC(rainfall, format = "e", digits = 12),
           file.path(out_dir, "hydro_rainfall.txt"))
writeLines(formatC(pet, format = "e", digits = 12),
           file.path(out_dir, "hydro_pet.txt"))

cat("wrote small config to", out_dir, "\n")
cat("  grid          :", x_extent, "x", y_extent, "=", n_cells, "cells\n")
cat("  active cells  :", sum(status_vec == "on"), "\n")
cat("  timesteps/day :", daily_timesteps, "\n")
cat("  days          :", n_days, " (t_extent =", t_extent, ")\n")
