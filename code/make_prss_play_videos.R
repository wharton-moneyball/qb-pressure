suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(deldir)
})

# Render 5 PRSS play videos with highlighted QB pocket region (Voronoi cell)
# and frame-level pocket area shift annotation.

set.seed(42)

root_dir <- getwd()
out_dir <- file.path(root_dir, "data", "processed", "prss_play_videos")
frames_root <- file.path(out_dir, "frames")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(frames_root, recursive = TRUE, showWarnings = FALSE)

tracking_path <- file.path(root_dir, "data", "processed", "bdb-tracking-all.parquet")
players_path <- file.path(root_dir, "data", "bdb", "players.csv")
pff_path <- file.path(root_dir, "data", "bdb", "pffScoutingData.csv")
summary_path <- file.path(root_dir, "data", "processed", "all_summary_k5.rds")

tracking <- read_parquet(tracking_path) %>%
  filter(!is.na(nflId), !is.na(x), !is.na(y))

players <- read_csv(players_path, show_col_types = FALSE) %>%
  select(nflId, officialPosition, displayName)
pff <- read_csv(pff_path, show_col_types = FALSE)

tracking <- tracking %>%
  left_join(players, by = "nflId")

all_summary_k5 <- readRDS(summary_path)

release_events <- c("pass_forward", "qb_sack", "qb_spike")
fallback_positions <- c(
  "QB", "C", "G", "T", "TE", "FB", "RB", "WR",
  "DE", "DT", "NT", "ILB", "OLB", "MLB", "LB",
  "CB", "DB", "FS", "SS"
)
back_cap_yards <- 3
side_cap_yards <- 3

poly_area <- function(x, y) {
  if (length(x) < 3 || length(y) < 3) return(NA_real_)
  x2 <- c(x, x[1])
  y2 <- c(y, y[1])
  abs(sum(x2[-1] * y2[-length(y2)] - x2[-length(x2)] * y2[-1])) / 2
}

clip_polygon <- function(x, y, inside_fn, intersect_fn) {
  n <- length(x)
  if (n < 3) return(list(x = numeric(0), y = numeric(0)))

  out_x <- numeric(0)
  out_y <- numeric(0)

  for (i in seq_len(n)) {
    j <- if (i == 1) n else i - 1
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
  clip_polygon(
    x, y,
    inside_fn = function(px, py) px >= xmin,
    intersect_fn = function(x1, y1, x2, y2) {
      dx <- x2 - x1
      t <- if (abs(dx) < 1e-12) 0 else (xmin - x1) / dx
      c(xmin, y1 + t * (y2 - y1))
    }
  )
}

clip_xmax <- function(x, y, xmax) {
  clip_polygon(
    x, y,
    inside_fn = function(px, py) px <= xmax,
    intersect_fn = function(x1, y1, x2, y2) {
      dx <- x2 - x1
      t <- if (abs(dx) < 1e-12) 0 else (xmax - x1) / dx
      c(xmax, y1 + t * (y2 - y1))
    }
  )
}

clip_ymin <- function(x, y, ymin) {
  clip_polygon(
    x, y,
    inside_fn = function(px, py) py >= ymin,
    intersect_fn = function(x1, y1, x2, y2) {
      dy <- y2 - y1
      t <- if (abs(dy) < 1e-12) 0 else (ymin - y1) / dy
      c(x1 + t * (x2 - x1), ymin)
    }
  )
}

clip_ymax <- function(x, y, ymax) {
  clip_polygon(
    x, y,
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

pff_pocket_ids <- pff %>%
  filter(pff_role %in% c("Pass Block", "Pass Rush"), !is.na(nflId)) %>%
  distinct(gameId, playId, nflId) %>%
  group_by(gameId, playId) %>%
  summarise(pocket_ids = list(unique(nflId)), .groups = "drop")

play_meta <- tracking %>%
  group_by(gameId, playId) %>%
  summarise(
    snap_frame = min(frameId[event == "ball_snap"], na.rm = TRUE),
    release_frame = min(frameId[event %in% release_events], na.rm = TRUE),
    has_qb = any(officialPosition == "QB", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    snap_frame = ifelse(is.infinite(snap_frame), NA_real_, snap_frame),
    release_frame = ifelse(is.infinite(release_frame), NA_real_, release_frame),
    frame_span = release_frame - snap_frame + 1,
    valid_play = !is.na(snap_frame) & !is.na(release_frame) & release_frame >= snap_frame & has_qb & frame_span >= 25
  )

selected_plays <- all_summary_k5 %>%
  filter(!is.na(max_neg_change), !is.na(qb)) %>%
  arrange(max_neg_change) %>%
  inner_join(play_meta %>% filter(valid_play), by = c("gameId", "playId")) %>%
  left_join(pff_pocket_ids, by = c("gameId", "playId")) %>%
  select(gameId, playId, qb, max_neg_change, snap_frame, release_frame, pocket_ids) %>%
  distinct(gameId, playId, .keep_all = TRUE) %>%
  slice_head(n = 5)

if (nrow(selected_plays) < 5) {
  stop("Could not find 5 valid plays to render.")
}

message("Selected plays:")
print(selected_plays)

safe_shift_label <- function(x) {
  if (is.na(x) || !is.finite(x)) return("NA")
  sprintf("%+.2f", x)
}

render_play <- function(play_row, play_rank) {
  gid <- play_row$gameId[[1]]
  pid <- play_row$playId[[1]]
  qb_name <- play_row$qb[[1]]
  max_prss <- play_row$max_neg_change[[1]]
  snap_frame <- play_row$snap_frame[[1]]
  release_frame <- play_row$release_frame[[1]]
  pocket_ids <- unique(play_row$pocket_ids[[1]])
  pocket_ids <- pocket_ids[!is.na(pocket_ids)]
  has_pff_ids <- length(pocket_ids) > 0

  play_data <- tracking %>%
    filter(gameId == gid, playId == pid, frameId >= snap_frame, frameId <= release_frame)

  frames <- sort(unique(play_data$frameId))
  if (length(frames) < 2) return(NULL)

  qb_team <- play_data %>%
    filter(frameId == snap_frame, officialPosition == "QB", !is.na(team)) %>%
    slice(1) %>%
    pull(team)

  if (length(qb_team) == 0) {
    qb_team <- play_data %>%
      filter(officialPosition == "QB", !is.na(team)) %>%
      slice(1) %>%
      pull(team)
  }

  if (length(qb_team) == 0) qb_team <- NA_character_

  play_slug <- sprintf("play_%02d_g%s_p%s", play_rank, gid, pid)
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

    if (nrow(frame_calc) < 3) {
      frame_calc <- frame_all %>%
        filter(officialPosition %in% fallback_positions) %>%
        distinct(x, y, .keep_all = TRUE) %>%
        mutate(ptNum = row_number())
    }

    qb_area <- NA_real_
    shift <- NA_real_
    qb_tile <- NULL

    if (nrow(frame_calc) >= 3) {
      vor <- tryCatch(
        deldir(frame_calc$x, frame_calc$y, rw = c(0, 120, 0, 53.3)),
        error = function(e) NULL
      )

      if (!is.null(vor)) {
        qb_row <- frame_calc %>%
          filter(officialPosition == "QB") %>%
          slice(1)

        if (nrow(qb_row) > 0) {
          tiles <- tile.list(vor)
          idx <- qb_row$ptNum[[1]]
          if (is.finite(idx) && idx >= 1 && idx <= length(tiles)) {
            tile <- tiles[[idx]]
            capped <- clip_qb_polygon(
              tile$x, tile$y,
              qb_x = qb_row$x[[1]],
              qb_y = qb_row$y[[1]],
              play_dir = qb_row$playDirection[[1]],
              back_cap = back_cap_yards
            )

            qb_area <- poly_area(capped$x, capped$y)
            if (is.finite(prev_area) && is.finite(qb_area)) {
              shift <- qb_area - prev_area
            }
            if (is.finite(qb_area)) {
              prev_area <- qb_area
            }

            if (length(capped$x) >= 3) {
              qb_tile <- list(x = capped$x, y = capped$y)
            }
          }
        }
      }
    }

    frame_file <- file.path(play_frame_dir, sprintf("frame_%04d.png", i))

    png(frame_file, width = 1600, height = 900, res = 150)
    par(mar = c(3.8, 3.8, 3.8, 1.2), xaxs = "i", yaxs = "i")

    plot(NA,
      xlim = c(0, 120), ylim = c(0, 53.3), asp = 1,
      xlab = "Field X (yards)", ylab = "Field Y (yards)",
      main = sprintf("Play %d/5  |  gameId %s, playId %s  |  QB: %s", play_rank, gid, pid, qb_name)
    )

    rect(0, 0, 120, 53.3, col = "#1f5d2d", border = "#f0f0f0", lwd = 1.4)
    abline(v = seq(10, 110, by = 10), col = adjustcolor("white", alpha.f = 0.25), lty = 1)
    abline(h = c(0, 53.3), col = "white", lwd = 1.2)

    if (!is.null(qb_tile)) {
      polygon(qb_tile$x, qb_tile$y,
        col = adjustcolor("gold", alpha.f = 0.30),
        border = "gold", lwd = 2
      )
    }

    defense <- frame_all %>% filter(side == "Defense")
    offense <- frame_all %>% filter(side == "Offense")
    qb_now <- frame_all %>% filter(officialPosition == "QB") %>% slice(1)

    if (nrow(defense) > 0) {
      points(defense$x, defense$y, pch = 21, cex = 1.05, bg = "#d95f5f", col = "black")
    }
    if (nrow(offense) > 0) {
      points(offense$x, offense$y, pch = 21, cex = 1.05, bg = "#5ab4e6", col = "black")
    }
    if (nrow(qb_now) > 0) {
      points(qb_now$x, qb_now$y, pch = 23, cex = 1.8, bg = "gold", col = "black", lwd = 1.2)
    }

    # Info box
    rect(1.0, 46.5, 66.0, 53.0, col = adjustcolor("black", alpha.f = 0.52), border = NA)
    text(2.0, 51.9, sprintf("Frame %d (relative %d)  |  Event window: snap=%d to release=%d", f, f - snap_frame, snap_frame, release_frame), adj = c(0, 1), col = "white", cex = 0.9)
    text(2.0, 50.1, sprintf("Pocket area (QB Voronoi): %s", ifelse(is.finite(qb_area), sprintf("%.2f", qb_area), "NA")), adj = c(0, 1), col = "white", cex = 0.9)
    text(2.0, 48.3, sprintf("Pocket area shift (dA): %s sq yd/frame", safe_shift_label(shift)), adj = c(0, 1), col = "white", cex = 0.95)
    text(2.0, 46.9, sprintf("Play PRSS (k=5): %.2f  |  Semicircle cap: r=%.1f yds", max_prss, back_cap_yards), adj = c(0, 1), col = "#ffd966", cex = 0.9)

    legend("bottomleft",
      legend = c("Offense", "Defense", "QB", "Pocket region (QB Voronoi cell)"),
      pch = c(21, 21, 23, 15),
      pt.bg = c("#5ab4e6", "#d95f5f", "gold", adjustcolor("gold", alpha.f = 0.3)),
      col = c("black", "black", "black", "gold"),
      bty = "n", cex = 0.85
    )

    dev.off()
  }

  video_path <- file.path(out_dir, paste0(play_slug, ".mp4"))
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

  # Optional cleanup of frame PNGs to save space
  unlink(play_frame_dir, recursive = TRUE, force = TRUE)

  tibble(
    play_rank = play_rank,
    gameId = gid,
    playId = pid,
    qb = qb_name,
    max_neg_change = max_prss,
    video_path = video_path
  )
}

video_rows <- map_dfr(seq_len(nrow(selected_plays)), function(i) {
  render_play(selected_plays[i, ], i)
})

if (nrow(video_rows) == 0) {
  stop("No videos were rendered.")
}

concat_file <- file.path(out_dir, "concat_list.txt")
writeLines(sprintf("file '%s'", normalizePath(video_rows$video_path, winslash = "/")), con = concat_file)

combined_video <- file.path(out_dir, "prss_5_play_reel.mp4")
combine_args <- c(
  "-y",
  "-f", "concat",
  "-safe", "0",
  "-i", concat_file,
  "-c:v", "libx264",
  "-pix_fmt", "yuv420p",
  "-movflags", "+faststart",
  combined_video
)

system2("ffmpeg", args = combine_args)

manifest_path <- file.path(out_dir, "video_manifest.csv")
write_csv(video_rows, manifest_path)

message("Rendered videos:")
print(video_rows)
message("Combined reel: ", combined_video)
message("Manifest: ", manifest_path)
