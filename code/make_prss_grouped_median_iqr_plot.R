suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(deldir)
})

set.seed(42)

root <- getwd()
out_file <- file.path(root, "data", "processed", "prss_area_vs_time_median_rate_iqr_by_pressure_group_capped.png")

tracking <- read_parquet(file.path(root, "data", "processed", "bdb-tracking-all.parquet")) %>%
  filter(!is.na(nflId), !is.na(x), !is.na(y))
players <- read_csv(file.path(root, "data", "bdb", "players.csv"), show_col_types = FALSE) %>%
  select(nflId, officialPosition)
pff <- read_csv(file.path(root, "data", "bdb", "pffScoutingData.csv"), show_col_types = FALSE)

tracking <- tracking %>% left_join(players, by = "nflId")

fallback_positions <- c(
  "QB", "C", "G", "T", "TE", "FB", "RB", "WR",
  "DE", "DT", "NT", "ILB", "OLB", "MLB", "LB",
  "CB", "DB", "FS", "SS"
)
back_cap_yards <- 3
side_cap_yards <- 3

poly_area <- function(x, y) {
  if (length(x) < 3 || length(y) < 3) return(NA_real_)
  x2 <- c(x, x[1]); y2 <- c(y, y[1])
  abs(sum(x2[-1] * y2[-length(y2)] - x2[-length(x2)] * y2[-1])) / 2
}

clip_polygon <- function(x, y, inside_fn, intersect_fn) {
  n <- length(x)
  if (n < 3) return(list(x = numeric(0), y = numeric(0)))
  out_x <- numeric(0); out_y <- numeric(0)
  for (i in seq_len(n)) {
    j <- if (i == 1) n else i - 1
    x1 <- x[j]; y1 <- y[j]; x2 <- x[i]; y2 <- y[i]
    in1 <- inside_fn(x1, y1); in2 <- inside_fn(x2, y2)
    if (in2) {
      if (!in1) {
        p <- intersect_fn(x1, y1, x2, y2)
        out_x <- c(out_x, p[1]); out_y <- c(out_y, p[2])
      }
      out_x <- c(out_x, x2); out_y <- c(out_y, y2)
    } else if (in1) {
      p <- intersect_fn(x1, y1, x2, y2)
      out_x <- c(out_x, p[1]); out_y <- c(out_y, p[2])
    }
  }
  list(x = out_x, y = out_y)
}

clip_xmin <- function(x, y, xmin) {
  clip_polygon(x, y,
    inside_fn = function(px, py) px >= xmin,
    intersect_fn = function(x1, y1, x2, y2) {
      dx <- x2 - x1
      t <- if (abs(dx) < 1e-12) 0 else (xmin - x1) / dx
      c(xmin, y1 + t * (y2 - y1))
    }
  )
}

clip_xmax <- function(x, y, xmax) {
  clip_polygon(x, y,
    inside_fn = function(px, py) px <= xmax,
    intersect_fn = function(x1, y1, x2, y2) {
      dx <- x2 - x1
      t <- if (abs(dx) < 1e-12) 0 else (xmax - x1) / dx
      c(xmax, y1 + t * (y2 - y1))
    }
  )
}

clip_ymin <- function(x, y, ymin) {
  clip_polygon(x, y,
    inside_fn = function(px, py) py >= ymin,
    intersect_fn = function(x1, y1, x2, y2) {
      dy <- y2 - y1
      t <- if (abs(dy) < 1e-12) 0 else (ymin - y1) / dy
      c(x1 + t * (x2 - x1), ymin)
    }
  )
}

clip_ymax <- function(x, y, ymax) {
  clip_polygon(x, y,
    inside_fn = function(px, py) py <= ymax,
    intersect_fn = function(x1, y1, x2, y2) {
      dy <- y2 - y1
      t <- if (abs(dy) < 1e-12) 0 else (ymax - y1) / dy
      c(x1 + t * (x2 - x1), ymax)
    }
  )
}

make_back_semicircle <- function(cx, cy, r, play_dir, n_pts = 60, far = 200) {
  dir    <- tolower(as.character(play_dir))
  if (dir == "right") {
    thetas <- seq(3*pi/2, pi/2, length.out = n_pts)
    arc_x  <- cx + r * cos(thetas)
    arc_y  <- cy + r * sin(thetas)
    list(x = c(arc_x, cx + far, cx + far),
         y = c(arc_y, cy + r,   cy - r))
  } else {
    thetas <- seq(pi/2, -pi/2, length.out = n_pts)
    arc_x  <- cx + r * cos(thetas)
    arc_y  <- cy + r * sin(thetas)
    list(x = c(arc_x, cx - far, cx - far),
         y = c(arc_y, cy - r,   cy + r))
  }
}

clip_to_convex_polygon <- function(sx, sy, cx, cy) {
  out <- list(x = sx, y = sy)
  n   <- length(cx)
  for (i in seq_len(n)) {
    j  <- if (i == n) 1L else i + 1L
    ax <- cx[i]; ay <- cy[i]
    bx <- cx[j]; by <- cy[j]
    ex <- bx - ax; ey <- by - ay
    out <- clip_polygon(out$x, out$y,
      inside_fn    = function(px, py) (px - ax) * ey - (py - ay) * ex >= 0,
      intersect_fn = function(x1, y1, x2, y2) {
        dx <- x2 - x1; dy <- y2 - y1
        denom <- dx * ey - dy * ex
        if (abs(denom) < 1e-12) return(c(x1, y1))
        t <- ((ax - x1) * ey - (ay - y1) * ex) / denom
        c(x1 + t * dx, y1 + t * dy)
      }
    )
    if (length(out$x) < 3) break
  }
  out
}

clip_qb_polygon <- function(x, y, qb_x, qb_y, play_dir, back_cap, n_pts = 60) {
  if (!is.finite(qb_x) || !is.finite(qb_y) || !is.finite(back_cap) || back_cap <= 0)
    return(list(x = x, y = y))
  semi <- make_back_semicircle(qb_x, qb_y, back_cap, play_dir, n_pts)
  clip_to_convex_polygon(x, y, semi$x, semi$y)
}

pressure_labels <- pff %>%
  filter(pff_role == "Pass Rush") %>%
  group_by(gameId, playId) %>%
  summarise(
    any_sack = any(pff_sack == 1, na.rm = TRUE),
    any_hurry_or_hit = any(pff_hurry == 1 | pff_hit == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pressure_group = case_when(
      any_sack ~ "Sack",
      any_hurry_or_hit ~ "Pressure (Hit/Hurry)",
      TRUE ~ "Neither"
    ),
    pressure_group = factor(pressure_group, levels = c("Sack", "Pressure (Hit/Hurry)", "Neither"))
  )

pff_pocket_ids <- pff %>%
  filter(pff_role %in% c("Pass Block", "Pass Rush"), !is.na(nflId)) %>%
  distinct(gameId, playId, nflId) %>%
  group_by(gameId, playId) %>%
  summarise(pocket_ids = list(unique(nflId)), .groups = "drop")

valid_plays <- tracking %>%
  group_by(gameId, playId) %>%
  summarise(
    has_snap = any(event == "ball_snap", na.rm = TRUE),
    has_release = any(event %in% c("pass_forward", "qb_sack", "qb_spike"), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(has_snap, has_release) %>%
  inner_join(pressure_labels, by = c("gameId", "playId")) %>%
  left_join(pff_pocket_ids, by = c("gameId", "playId"))

sample_per_group <- 300L
sampled_plays <- valid_plays %>%
  group_split(pressure_group) %>%
  map_dfr(function(df) df %>% slice_sample(n = min(sample_per_group, nrow(df))))

get_area_vs_time <- function(game_id, play_id, pocket_ids = NULL) {
  play_data <- tracking %>% filter(gameId == game_id, playId == play_id)
  snap_frame <- suppressWarnings(min(play_data$frameId[play_data$event == "ball_snap"], na.rm = TRUE))
  release_frame <- suppressWarnings(min(play_data$frameId[play_data$event %in% c("pass_forward", "qb_sack", "qb_spike")], na.rm = TRUE))
  if (!is.finite(snap_frame) || !is.finite(release_frame) || release_frame < snap_frame) return(NULL)

  cleaned_ids <- unique(pocket_ids[!is.na(pocket_ids)])
  has_pff_ids <- length(cleaned_ids) > 0

  play_data <- if (has_pff_ids) {
    play_data %>%
      filter(frameId >= snap_frame, frameId <= release_frame,
             (officialPosition == "QB" | nflId %in% cleaned_ids),
             !is.na(x), !is.na(y))
  } else {
    play_data %>%
      filter(frameId >= snap_frame, frameId <= release_frame,
             officialPosition %in% fallback_positions,
             !is.na(x), !is.na(y))
  }

  frames <- sort(unique(play_data$frameId))

  map_dfr(frames, function(f) {
    frame_data <- play_data %>%
      filter(frameId == f) %>%
      distinct(x, y, .keep_all = TRUE) %>%
      mutate(ptNum = row_number())

    if (nrow(frame_data) < 3) return(NULL)

    tryCatch({
      qb_seed <- frame_data %>% filter(officialPosition == "QB") %>% slice(1)
      if (nrow(qb_seed) == 0) return(NULL)

      vor <- deldir(frame_data$x, frame_data$y, rw = c(0, 120, 0, 53.3))
      tiles <- tile.list(vor)
      idx <- qb_seed$ptNum[[1]]
      if (!is.finite(idx) || idx < 1 || idx > length(tiles)) return(NULL)

      tile <- tiles[[idx]]
      capped <- clip_qb_polygon(
        tile$x, tile$y,
        qb_x = qb_seed$x[[1]],
        qb_y = qb_seed$y[[1]],
        play_dir = qb_seed$playDirection[[1]],
        back_cap = back_cap_yards
      )
      area <- poly_area(capped$x, capped$y)
      if (!is.finite(area)) return(NULL)

      tibble(
        gameId = game_id,
        playId = play_id,
        frames_since_snap = f - snap_frame,
        qb_area = area
      )
    }, error = function(e) NULL)
  })
}

area_vs_time <- imap_dfr(seq_len(nrow(sampled_plays)), function(i, ...) {
  if (i %% 100 == 0) message("Progress: ", i, "/", nrow(sampled_plays))
  get_area_vs_time(
    sampled_plays$gameId[i],
    sampled_plays$playId[i],
    pocket_ids = sampled_plays$pocket_ids[[i]]
  )
}) %>%
  left_join(sampled_plays %>% select(gameId, playId, pressure_group), by = c("gameId", "playId"))

rate_summary <- area_vs_time %>%
  mutate(second_bin = round((frames_since_snap / 10), 1)) %>%
  arrange(gameId, playId, frames_since_snap) %>%
  group_by(gameId, playId, pressure_group) %>%
  mutate(area_change = qb_area - lag(qb_area)) %>%
  ungroup() %>%
  filter(!is.na(area_change)) %>%
  group_by(pressure_group, second_bin) %>%
  summarise(
    n = n(),
    median_rate = median(area_change, na.rm = TRUE),
    q25 = quantile(area_change, 0.25, na.rm = TRUE),
    q75 = quantile(area_change, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n >= 20)

palette_vals <- c(
  "Sack" = "#b22222",
  "Pressure (Hit/Hurry)" = "#d97706",
  "Neither" = "#1f77b4"
)

p <- ggplot(rate_summary, aes(x = second_bin, y = median_rate, color = pressure_group, fill = pressure_group)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.20, color = NA) +
  geom_line(linewidth = 1.15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray45") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgreen") +
  annotate("text", x = 0.05, y = Inf, label = "snap", vjust = 1.4, color = "darkgreen", size = 3.5) +
  scale_color_manual(values = palette_vals) +
  scale_fill_manual(values = palette_vals) +
  labs(
    title = "Median Frame-to-Frame Pocket Area Change by Pressure Outcome",
    subtitle = sprintf("Line = median | Band = IQR | Semicircle cap (r=%.1f yd)", back_cap_yards),
    x = "Seconds since snap",
    y = "Change in QB area (sq yards / frame)",
    color = "Play Outcome",
    fill = "Play Outcome"
  ) +
  theme_minimal(base_size = 13)

ggsave(out_file, p, width = 10.5, height = 5.8, dpi = 180)

cat("Saved:", out_file, "\n")
cat("Sampled plays by group:\n")
print(sampled_plays %>% count(pressure_group))
cat("Rate summary rows by group:\n")
print(rate_summary %>% count(pressure_group))
