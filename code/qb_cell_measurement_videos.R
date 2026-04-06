suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(deldir)
})

# =============================================================================
# make_qb_cell_videos.R
#
# Renders 5 play videos showing the dynamic QB cell (Voronoi cell with
# dynamic back cap). Selects plays where dynamic_back_cap varies most
# across frames (max - min), to best illustrate the JP dynamic cap idea.
#
# Requires: area_vs_time, qb_speeds, sampled_plays, tracking, players
# from qb_cell_measurement.R to already be in the environment.
# =============================================================================

# --- Check ffmpeg ------------------------------------------------------------
if (Sys.which("ffmpeg") == "") {
  stop("ffmpeg not found. Install via:\n  Mac: brew install ffmpeg\n  Linux: sudo apt install ffmpeg")
}

set.seed(42)

root_dir    <- getwd()
out_dir     <- file.path(root_dir, "data", "processed", "qb_cell_videos")
frames_root <- file.path(out_dir, "frames")
dir.create(out_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(frames_root, recursive = TRUE, showWarnings = FALSE)

# --- Reuse event definitions from qb_cell_measurement.R ---------------------
# (SNAP_EVENTS, RELEASE_EVENTS, SIDE_CAP_YARDS, MIN_BACK_CAP must be defined)

fallback_positions <- c(
  "QB", "C", "G", "T", "TE", "FB", "RB", "WR",
  "DE", "DT", "NT", "ILB", "OLB", "MLB", "LB",
  "CB", "DB", "FS", "SS"
)

# =============================================================================
# Select 5 plays where dynamic_back_cap varies most (max - min)
# =============================================================================
play_back_cap_range <- area_vs_time %>%
  group_by(gameId, playId) %>%
  summarise(
    back_cap_range = max(back_cap_used, na.rm = TRUE) - min(back_cap_used, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(back_cap_range))

# Get snap/release frame boundaries for valid play filtering
play_meta <- tracking %>%
  group_by(gameId, playId) %>%
  summarise(
    snap_frame    = suppressWarnings(min(frameId[event %in% SNAP_EVENTS],    na.rm = TRUE)),
    release_frame = suppressWarnings(min(frameId[event %in% RELEASE_EVENTS], na.rm = TRUE)),
    has_qb        = any(officialPosition == "QB", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    snap_frame    = ifelse(is.infinite(snap_frame),    NA_real_, snap_frame),
    release_frame = ifelse(is.infinite(release_frame), NA_real_, release_frame),
    frame_span    = release_frame - snap_frame + 1,
    valid_play    = !is.na(snap_frame) & !is.na(release_frame) &
      release_frame >= snap_frame & has_qb & frame_span >= 10
  )

# Get QB display name for each play
qb_names <- read_csv(
  file.path(root_dir, "data", "bdb", "players.csv"),
  show_col_types = FALSE
) %>%
  filter(officialPosition == "QB") %>%
  distinct(nflId, displayName) %>%
  inner_join(
    tracking %>% distinct(gameId, playId, nflId),
    by = "nflId"
  ) %>%
  group_by(gameId, playId) %>%
  slice(1) %>%
  ungroup() %>%
  select(gameId, playId, qb = displayName)

selected_plays <- play_back_cap_range %>%
  inner_join(play_meta %>% filter(valid_play), by = c("gameId", "playId")) %>%
  left_join(sampled_plays %>% select(gameId, playId, pocket_ids), by = c("gameId", "playId")) %>%
  left_join(qb_names, by = c("gameId", "playId")) %>%
  distinct(gameId, playId, .keep_all = TRUE) %>%
  slice_head(n = 5)

if (nrow(selected_plays) < 5) {
  warning(sprintf("Only found %d valid plays to render.", nrow(selected_plays)))
}

message("Selected plays:")
print(selected_plays %>% select(gameId, playId, qb, back_cap_range))

# =============================================================================
# Helper
# =============================================================================
safe_shift_label <- function(x) {
  if (is.na(x) || !is.finite(x)) return("NA")
  sprintf("%+.2f", x)
}

# =============================================================================
# Render function
# =============================================================================
render_play <- function(play_row, play_rank) {
  gid            <- play_row$gameId[[1]]
  pid            <- play_row$playId[[1]]
  qb_name        <- play_row$qb[[1]]
  back_cap_range <- play_row$back_cap_range[[1]]
  snap_frame     <- play_row$snap_frame[[1]]
  release_frame  <- play_row$release_frame[[1]]
  pocket_ids     <- unique(play_row$pocket_ids[[1]])
  pocket_ids     <- pocket_ids[!is.na(pocket_ids)]
  has_pff_ids    <- length(pocket_ids) > 0
  
  play_data <- tracking %>%
    filter(gameId == gid, playId == pid,
           frameId >= snap_frame, frameId <= release_frame)
  
  frames <- sort(unique(play_data$frameId))
  if (length(frames) < 2) return(NULL)
  
  # Determine QB team for offense/defense coloring
  qb_team <- play_data %>%
    filter(officialPosition == "QB", !is.na(team)) %>%
    slice(1) %>%
    pull(team)
  if (length(qb_team) == 0) qb_team <- NA_character_
  
  # Pre-fetch back_cap values for this play from area_vs_time
  play_back_caps <- area_vs_time %>%
    filter(gameId == gid, playId == pid) %>%
    select(frameId, back_cap_used)
  
  play_slug      <- sprintf("play_%02d_g%s_p%s", play_rank, gid, pid)
  play_frame_dir <- file.path(frames_root, play_slug)
  dir.create(play_frame_dir, recursive = TRUE, showWarnings = FALSE)
  
  prev_area <- NA_real_
  
  for (i in seq_along(frames)) {
    f <- frames[[i]]
    
    frame_all <- play_data %>%
      filter(frameId == f) %>%
      mutate(side = ifelse(!is.na(qb_team) & team == qb_team, "Offense", "Defense"))
    
    frame_calc <- if (has_pff_ids) {
      frame_all %>%
        filter(officialPosition == "QB" | nflId %in% pocket_ids) %>%
        distinct(x, y, .keep_all = TRUE) %>%
        mutate(ptNum = row_number())
    } else {
      frame_all %>%
        filter(officialPosition %in% fallback_positions) %>%
        distinct(x, y, .keep_all = TRUE) %>%
        mutate(ptNum = row_number())
    }
    
    # Fallback if too few players
    if (nrow(frame_calc) < 3) {
      frame_calc <- frame_all %>%
        filter(officialPosition %in% fallback_positions) %>%
        distinct(x, y, .keep_all = TRUE) %>%
        mutate(ptNum = row_number())
    }
    
    qb_area   <- NA_real_
    shift     <- NA_real_
    qb_tile   <- NULL
    this_back_cap <- MIN_BACK_CAP
    
    if (nrow(frame_calc) >= 3) {
      vor <- tryCatch(
        deldir(frame_calc$x, frame_calc$y, rw = c(0, 120, 0, 53.3)),
        error = function(e) NULL
      )
      
      if (!is.null(vor)) {
        qb_row <- frame_calc %>% filter(officialPosition == "QB") %>% slice(1)
        
        if (nrow(qb_row) > 0) {
          tiles <- tile.list(vor)
          idx   <- qb_row$ptNum[[1]]
          
          if (is.finite(idx) && idx >= 1 && idx <= length(tiles)) {
            tile <- tiles[[idx]]
            
            # Look up dynamic back cap for this frame
            bc <- play_back_caps %>% filter(frameId == f) %>% pull(back_cap_used)
            this_back_cap <- if (length(bc) == 0 || !is.finite(bc[[1]])) MIN_BACK_CAP else bc[[1]]
            
            capped <- clip_qb_polygon(
              tile$x, tile$y,
              qb_x     = qb_row$x[[1]],
              qb_y     = qb_row$y[[1]],
              play_dir = qb_row$playDirection[[1]],
              back_cap = this_back_cap
            )
            
            qb_area <- poly_area(capped$x, capped$y)
            if (is.finite(prev_area) && is.finite(qb_area)) shift <- qb_area - prev_area
            if (is.finite(qb_area)) prev_area <- qb_area
            if (length(capped$x) >= 3) qb_tile <- list(x = capped$x, y = capped$y)
          }
        }
      }
    }
    
    # --- Draw frame ----------------------------------------------------------
    frame_file <- file.path(play_frame_dir, sprintf("frame_%04d.png", i))
    png(frame_file, width = 1600, height = 900, res = 150)
    par(mar = c(3.8, 3.8, 3.8, 1.2), xaxs = "i", yaxs = "i")
    
    plot(NA,
         xlim = c(0, 120), ylim = c(0, 53.3), asp = 1,
         xlab = "Field X (yards)", ylab = "Field Y (yards)",
         main = sprintf("Play %d/5  |  gameId %s, playId %s  |  QB: %s",
                        play_rank, gid, pid, qb_name)
    )
    
    # Field
    rect(0, 0, 120, 53.3, col = "#1f5d2d", border = "#f0f0f0", lwd = 1.4)
    abline(v = seq(10, 110, by = 10), col = adjustcolor("white", alpha.f = 0.25), lty = 1)
    abline(h = c(0, 53.3), col = "white", lwd = 1.2)
    
    # QB cell: filled region + explicit boundary line
    if (!is.null(qb_tile)) {
      polygon(qb_tile$x, qb_tile$y,
              col    = adjustcolor("gold", alpha.f = 0.30),
              border = NA
      )
      lines(
        c(qb_tile$x, qb_tile$x[1]),
        c(qb_tile$y, qb_tile$y[1]),
        col = "gold", lwd = 2.5
      )
    }
    
    # Players
    defense <- frame_all %>% filter(side == "Defense")
    offense <- frame_all %>% filter(side == "Offense")
    qb_now  <- frame_all %>% filter(officialPosition == "QB") %>% slice(1)
    
    if (nrow(defense) > 0)
      points(defense$x, defense$y, pch = 21, cex = 1.05, bg = "#d95f5f", col = "black")
    if (nrow(offense) > 0)
      points(offense$x, offense$y, pch = 21, cex = 1.05, bg = "#5ab4e6", col = "black")
    if (nrow(qb_now) > 0)
      points(qb_now$x, qb_now$y, pch = 23, cex = 1.8, bg = "gold", col = "black", lwd = 1.2)
    
    # Info box
    rect(1.0, 44.0, 72.0, 53.0, col = adjustcolor("black", alpha.f = 0.52), border = NA)
    text(2.0, 51.9,
         sprintf("Frame %d (t+%d)  |  snap=%d  release=%d",
                 f, f - snap_frame, snap_frame, release_frame),
         adj = c(0, 1), col = "white", cex = 0.9)
    text(2.0, 50.1,
         sprintf("QB cell area: %s sq yds",
                 ifelse(is.finite(qb_area), sprintf("%.2f", qb_area), "NA")),
         adj = c(0, 1), col = "white", cex = 0.9)
    text(2.0, 48.3,
         sprintf("Area change (dA): %s sq yds/frame", safe_shift_label(shift)),
         adj = c(0, 1), col = "white", cex = 0.9)
    text(2.0, 46.5,
         sprintf("Dynamic back cap: %.2f yds  |  play range: %.2f yds",
                 this_back_cap, back_cap_range),
         adj = c(0, 1), col = "#ffd966", cex = 0.9)
    text(2.0, 44.7,
         sprintf("Semicircle cap (r=back_cap)  |  T_horizon: %.1fs",
                 T_HORIZON_SEC),
         adj = c(0, 1), col = "#ffd966", cex = 0.9)
    
    legend("bottomleft",
           legend = c("Offense", "Defense", "QB", "QB cell (dynamic back cap)"),
           pch    = c(21, 21, 23, 15),
           pt.bg  = c("#5ab4e6", "#d95f5f", "gold", adjustcolor("gold", alpha.f = 0.3)),
           col    = c("black", "black", "black", "gold"),
           bty = "n", cex = 0.85
    )
    
    dev.off()
  }
  
  # --- Encode video with ffmpeg ----------------------------------------------
  video_path  <- file.path(out_dir, paste0(play_slug, ".mp4"))
  ffmpeg_args <- c(
    "-y",
    "-framerate", "10",
    "-i", file.path(play_frame_dir, "frame_%04d.png"),
    "-c:v", "libx264",
    "-pix_fmt", "yuv420p",
    "-movflags", "+faststart",
    video_path
  )
  
  status <- system2("ffmpeg", args = ffmpeg_args)
  if (!identical(status, 0L)) {
    warning(sprintf("ffmpeg failed for %s", play_slug))
    return(NULL)
  }
  
  unlink(play_frame_dir, recursive = TRUE, force = TRUE)
  
  tibble(
    play_rank      = play_rank,
    gameId         = gid,
    playId         = pid,
    qb             = qb_name,
    back_cap_range = back_cap_range,
    video_path     = video_path
  )
}

# =============================================================================
# Run
# =============================================================================
video_rows <- map_dfr(seq_len(nrow(selected_plays)), function(i) {
  message("Rendering play ", i, "/", nrow(selected_plays))
  render_play(selected_plays[i, ], i)
})

if (nrow(video_rows) == 0) stop("No videos were rendered.")

# --- Combine into single reel ------------------------------------------------
concat_file <- file.path(out_dir, "concat_list.txt")
writeLines(
  sprintf("file '%s'", normalizePath(video_rows$video_path, winslash = "/")),
  con = concat_file
)

combined_video <- file.path(out_dir, "qb_cell_5_play_reel.mp4")
system2("ffmpeg", args = c(
  "-y", "-f", "concat", "-safe", "0",
  "-i", concat_file,
  "-c:v", "libx264", "-pix_fmt", "yuv420p", "-movflags", "+faststart",
  combined_video
))

write_csv(video_rows, file.path(out_dir, "video_manifest.csv"))

message("Done. Combined reel: ", combined_video)
print(video_rows)

