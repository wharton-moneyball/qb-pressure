suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(deldir)
})

set.seed(42)

# =============================================================================
# QB Cell Measurement — Dynamic Back Cap
#
#
#       dynamic_back_cap = qb_speed (yards/sec) * t_horizon (seconds)
#
#
#   The cap is clamped to [MIN_BACK_CAP, MAX_BACK_CAP] to handle noise and
#   the edge case where the pocket has already broken down (Lucian's concern).
#
# Parameters to tune:
#   T_HORIZON_SEC   — how many seconds ahead we project QB movement
#   MIN_BACK_CAP    — floor (yards), prevents cell from vanishing
#   MAX_BACK_CAP    — ceiling (yards), 9 yd covers full 7-step drop per PFF
#   SIDE_CAP_YARDS  — lateral cap, still fixed (unchanged from prior work)
# =============================================================================

root <- getwd()
out_file <- file.path(
  root, "data", "processed",
  "qb_cell_dynamic_back_cap_median_rate_iqr.png"
)

# --- Tunable parameters ------------------------------------------------------
T_HORIZON_SEC  <- 0.5   # seconds of forward-looking speed window
MIN_BACK_CAP   <- 1.5   # yards — floor
MAX_BACK_CAP   <- 9.0   # yards — ceiling (full 7-step drop)
SIDE_CAP_YARDS <- 3.0   # yards — lateral cap (unchanged)
SAMPLE_N       <- 500L  # number of plays to sample
FPS            <- 10    # tracking data frames per second
# -----------------------------------------------------------------------------

# --- Load data ----------------------------------------------------------------
tracking <- read_parquet(
  file.path(root, "data", "processed", "bdb-tracking-all.parquet")
) %>%
  filter(!is.na(nflId), !is.na(x), !is.na(y))

players <- read_csv(
  file.path(root, "data", "bdb", "players.csv"),
  show_col_types = FALSE
) %>%
  select(nflId, officialPosition)

pff <- read_csv(
  file.path(root, "data", "bdb", "pffScoutingData.csv"),
  show_col_types = FALSE
)

tracking <- tracking %>% left_join(players, by = "nflId")

# --- Fallback position & event filter ------------------------------------------------
fallback_positions <- c(
  "QB", "C", "G", "T", "TE", "FB", "RB", "WR",
  "DE", "DT", "NT", "ILB", "OLB", "MLB", "LB",
  "CB", "DB", "FS", "SS"
)

SNAP_EVENTS    <- c("ball_snap", "autoevent_ballsnap")

RELEASE_EVENTS <- c(
  "pass_forward", "autoevent_passforward", "autoevent_passinterrupted",
  "qb_sack", "qb_strip_sack", "qb_spike", "run"
)
# =============================================================================
# Helper: polygon area (shoelace formula)
# =============================================================================
poly_area <- function(x, y) {
  if (length(x) < 3 || length(y) < 3) return(NA_real_)
  x2 <- c(x, x[1])
  y2 <- c(y, y[1])
  abs(sum(x2[-1] * y2[-length(y2)] - x2[-length(x2)] * y2[-1])) / 2
}

# =============================================================================
# Sutherland-Hodgman polygon clipping (one half-plane at a time)
# =============================================================================
clip_polygon <- function(x, y, inside_fn, intersect_fn) {
  n <- length(x)
  if (n < 3) return(list(x = numeric(0), y = numeric(0)))
  out_x <- numeric(0)
  out_y <- numeric(0)
  for (i in seq_len(n)) {
    j  <- if (i == 1) n else i - 1
    x1 <- x[j]; y1 <- y[j]
    x2 <- x[i]; y2 <- y[i]
    in1 <- inside_fn(x1, y1)
    in2 <- inside_fn(x2, y2)
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
    inside_fn    = function(px, py) px >= xmin,
    intersect_fn = function(x1, y1, x2, y2) {
      dx <- x2 - x1
      t  <- if (abs(dx) < 1e-12) 0 else (xmin - x1) / dx
      c(xmin, y1 + t * (y2 - y1))
    }
  )
}

clip_xmax <- function(x, y, xmax) {
  clip_polygon(x, y,
    inside_fn    = function(px, py) px <= xmax,
    intersect_fn = function(x1, y1, x2, y2) {
      dx <- x2 - x1
      t  <- if (abs(dx) < 1e-12) 0 else (xmax - x1) / dx
      c(xmax, y1 + t * (y2 - y1))
    }
  )
}

clip_ymin <- function(x, y, ymin) {
  clip_polygon(x, y,
    inside_fn    = function(px, py) py >= ymin,
    intersect_fn = function(x1, y1, x2, y2) {
      dy <- y2 - y1
      t  <- if (abs(dy) < 1e-12) 0 else (ymin - y1) / dy
      c(x1 + t * (x2 - x1), ymin)
    }
  )
}

clip_ymax <- function(x, y, ymax) {
  clip_polygon(x, y,
    inside_fn    = function(px, py) py <= ymax,
    intersect_fn = function(x1, y1, x2, y2) {
      dy <- y2 - y1
      t  <- if (abs(dy) < 1e-12) 0 else (ymax - y1) / dy
      c(x1 + t * (x2 - x1), ymax)
    }
  )
}

# =============================================================================
# Core clip function — now accepts a dynamic back_cap per call
# =============================================================================
make_back_semicircle <- function(cx, cy, r, play_dir, n_pts = 60, far = 200) {
  dir    <- tolower(as.character(play_dir))
  if (dir == "right") {
    # arc: bottom (cx, cy-r) → left (cx-r, cy) → top (cx, cy+r), clockwise
    thetas <- seq(3*pi/2, pi/2, length.out = n_pts)
    arc_x  <- cx + r * cos(thetas)
    arc_y  <- cy + r * sin(thetas)
    # extend top and bottom forward (positive x) to leave forward space open
    list(x = c(arc_x, cx + far, cx + far),
         y = c(arc_y, cy + r,   cy - r))
  } else {
    # arc: top (cx, cy+r) → right (cx+r, cy) → bottom (cx, cy-r), clockwise
    thetas <- seq(pi/2, -pi/2, length.out = n_pts)
    arc_x  <- cx + r * cos(thetas)
    arc_y  <- cy + r * sin(thetas)
    # extend top and bottom forward (negative x) to leave forward space open
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

# =============================================================================
# Compute per-frame QB speed for all plays
#
# Speed is computed as the Euclidean distance between consecutive frames,
# converted to yards/sec (multiply by FPS).
# We then compute a rolling forward mean over T_HORIZON_SEC worth of frames
# so that back_cap at frame f reflects where the QB is heading next.
# =============================================================================
t_horizon_frames <- round(T_HORIZON_SEC * FPS)  # e.g. 0.5s * 10fps = 5 frames

play_windows <- tracking %>%
  group_by(gameId, playId) %>%
  summarise(
    snap_frame    = suppressWarnings(min(frameId[event %in% SNAP_EVENTS],    na.rm = TRUE)),
    release_frame = suppressWarnings(min(frameId[event %in% RELEASE_EVENTS], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  filter(is.finite(snap_frame), is.finite(release_frame), release_frame >= snap_frame)


qb_speeds <- tracking %>%
  filter(officialPosition == "QB") %>%
  arrange(gameId, playId, frameId) %>%
  group_by(gameId, playId, nflId) %>%
  inner_join(play_windows, by = c("gameId", "playId")) %>%
  filter(frameId >= snap_frame, frameId <= release_frame) %>%
  mutate(
    speed_forward = slider::slide_dbl(
      s, 
      mean,
      .before = 0, .after = t_horizon_frames - 1, .complete = FALSE
    )
  ) %>%
  ungroup() %>%
  mutate(
    # Clamp to [MIN_BACK_CAP, MAX_BACK_CAP]
    dynamic_back_cap = pmax(MIN_BACK_CAP, pmin(speed_forward * T_HORIZON_SEC, MAX_BACK_CAP)),
    # Replace NA (e.g. last frames where lead() is missing) with MIN_BACK_CAP
    dynamic_back_cap = if_else(is.na(dynamic_back_cap), MIN_BACK_CAP, dynamic_back_cap)
  ) %>%
  select(gameId, playId, nflId, frameId, speed = s, speed_forward, dynamic_back_cap)

# =============================================================================
# PFF pocket player IDs (Pass Blockers + Pass Rushers)
# =============================================================================
pff_pocket_ids <- pff %>%
  filter(pff_role %in% c("Pass Block", "Pass Rush"), !is.na(nflId)) %>%
  distinct(gameId, playId, nflId) %>%
  group_by(gameId, playId) %>%
  summarise(pocket_ids = list(unique(nflId)), .groups = "drop")

# =============================================================================
# Sample valid plays
# =============================================================================
valid_plays <- tracking %>%
  group_by(gameId, playId) %>%
  summarise(
    has_snap    = any(event %in% SNAP_EVENTS,    na.rm = TRUE),
    has_release = any(event %in% RELEASE_EVENTS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(has_snap, has_release) %>%
  left_join(pff_pocket_ids, by = c("gameId", "playId"))

sampled_plays <- valid_plays %>%
  slice_sample(n = min(SAMPLE_N, nrow(valid_plays)))

cat("Sampled plays:", nrow(sampled_plays), "\n")

# =============================================================================
# Main computation: frame-by-frame QB cell area with dynamic back cap
# =============================================================================
get_area_vs_time <- function(game_id, play_id, pocket_ids = NULL) {

  play_data <- tracking %>% filter(gameId == game_id, playId == play_id)

  snap_frame    <- suppressWarnings(
    min(play_data$frameId[play_data$event %in% SNAP_EVENTS], na.rm = TRUE)
  )
  release_frame <- suppressWarnings(
    min(play_data$frameId[play_data$event %in% RELEASE_EVENTS], na.rm = TRUE)
  )

  if (!is.finite(snap_frame) || !is.finite(release_frame)) return(NULL)
  if (release_frame < snap_frame) return(NULL)

  # Speed lookup for this play's QB
  speed_lookup <- qb_speeds %>%
    filter(gameId == game_id, playId == play_id) %>%
    select(frameId, dynamic_back_cap)

  cleaned_ids <- unique(pocket_ids[!is.na(pocket_ids)])
  has_pff_ids <- length(cleaned_ids) > 0

  play_data <- if (has_pff_ids) {
    play_data %>%
      filter(
        frameId >= snap_frame, frameId <= release_frame,
        (officialPosition == "QB" | nflId %in% cleaned_ids),
        !is.na(x), !is.na(y)
      )
  } else {
    play_data %>%
      filter(
        frameId >= snap_frame, frameId <= release_frame,
        officialPosition %in% fallback_positions,
        !is.na(x), !is.na(y)
      )
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

      # Look up the dynamic back cap for this frame; fall back to MIN if missing
      back_cap <- speed_lookup %>%
        filter(frameId == f) %>%
        pull(dynamic_back_cap)
      back_cap <- if (length(back_cap) == 0 || !is.finite(back_cap[[1]])) {
        MIN_BACK_CAP
      } else {
        back_cap[[1]]
      }

      vor   <- deldir(frame_data$x, frame_data$y, rw = c(0, 120, 0, 53.3))
      tiles <- tile.list(vor)
      idx   <- qb_seed$ptNum[[1]]
      if (!is.finite(idx) || idx < 1 || idx > length(tiles)) return(NULL)

      tile   <- tiles[[idx]]
      capped <- clip_qb_polygon(
        tile$x, tile$y,
        qb_x     = qb_seed$x[[1]],
        qb_y     = qb_seed$y[[1]],
        play_dir = qb_seed$playDirection[[1]],
        back_cap = back_cap
      )

      area <- poly_area(capped$x, capped$y)
      if (!is.finite(area)) return(NULL)

      tibble(
        gameId            = game_id,
        playId            = play_id,
        frameId           = f,
        frames_since_snap = f - snap_frame,
        qb_area           = area,
        back_cap_used     = back_cap   # keep for diagnostics
      )
    }, error = function(e) NULL)
  })
}

# --- Run ---------------------------------------------------------------------
area_vs_time <- imap_dfr(seq_len(nrow(sampled_plays)), function(i, ...) {
  if (i %% 100 == 0) message("Progress: ", i, "/", nrow(sampled_plays))
  get_area_vs_time(
    sampled_plays$gameId[i],
    sampled_plays$playId[i],
    pocket_ids = sampled_plays$pocket_ids[[i]]
  )
})

cat("Total frame-play rows:", nrow(area_vs_time), "\n")

# =============================================================================
# Summarise: median rate of change + IQR per 0.1-second bin
# =============================================================================
rate_summary <- area_vs_time %>%
  mutate(second_bin = round(frames_since_snap / FPS, 1)) %>%
  arrange(gameId, playId, frames_since_snap) %>%
  group_by(gameId, playId) %>%
  mutate(area_change = qb_area - lag(qb_area)) %>%
  ungroup() %>%
  filter(!is.na(area_change)) %>%
  group_by(second_bin) %>%
  summarise(
    n           = n(),
    median_rate = median(area_change, na.rm = TRUE),
    q25         = quantile(area_change, 0.25, na.rm = TRUE),
    q75         = quantile(area_change, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n >= 20)

# =============================================================================
# Diagnostic: distribution of dynamic_back_cap values
# =============================================================================
cat("\nDynamic back cap summary (yards):\n")
print(summary(area_vs_time$back_cap_used))

# =============================================================================
# Plot
# =============================================================================
p <- ggplot(rate_summary, aes(x = second_bin)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "tomato", alpha = 0.25) +
  geom_line(aes(y = median_rate), color = "tomato", linewidth = 1.2) +
  geom_hline(yintercept = 0,  linetype = "dashed", color = "gray45") +
  geom_vline(xintercept = 0,  linetype = "dashed", color = "darkgreen") +
  annotate(
    "text", x = 0.05, y = Inf, label = "snap",
    vjust = 1.4, color = "darkgreen", size = 3.5
  ) +
  labs(
    title    = "Median Frame-to-Frame QB Cell Area Change (Dynamic Back Cap)",
    subtitle = sprintf(
      "Back cap = speed × %.1fs  |  clamped [%.1f, %.1f] yds  |  Semicircle cap",
      T_HORIZON_SEC, MIN_BACK_CAP, MAX_BACK_CAP
    ),
    x = "Seconds since snap",
    y = "Change in QB cell area (sq yards / frame)"
  ) +
  theme_minimal(base_size = 13)

ggsave(out_file, p, width = 10.5, height = 5.8, dpi = 180)
cat("Saved:", out_file, "\n")
cat("Rate summary rows:", nrow(rate_summary), "\n")

summary(area_vs_time$back_cap_used)
# =============================================================================
yards_to_meters <- 0.9144

players_full <- read_csv(
  file.path(root, "data", "bdb", "players.csv"),
  show_col_types = FALSE
)

qb_max_speed <- tracking %>%
  filter(officialPosition == "QB", !is.na(s)) %>%
  left_join(players_full %>% select(nflId, displayName), by = "nflId") %>%
  group_by(nflId, displayName) %>%
  summarise(
    max_speed_yds_per_sec = max(s, na.rm = TRUE),
    max_speed_m_per_sec   = max(s, na.rm = TRUE) * yards_to_meters,
    .groups = "drop"
  ) %>%
  arrange(desc(max_speed_m_per_sec))

write_csv(qb_max_speed, file.path(root, "data", "processed", "qb_max_speed.csv"))
