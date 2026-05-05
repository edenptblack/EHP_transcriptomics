# R/logging.R

start_log = function(script_name, log_dir = 'logs') {
  options(readr.show_progress = FALSE)
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  log_file = file.path(log_dir, sprintf('%s_%s.log', script_name,
                                        format(Sys.time(), '%Y%m%d_%H%M%S')))
  
  # Direct output
  log_con = file(log_file, open = 'wt')
  sink(log_con, append = TRUE, split = TRUE)
  sink(log_con, append = TRUE, type = 'message')
  
  cat(sprintf('%s started at %s\n', script_name, Sys.time()))
  cat(sprintf('R version: %s\n', R.version.string))
  
  invisible(log_con)
}

end_log = function(log_con) {
  cat(sprintf('\nFinished at %s\n', Sys.time()))
  sink(type = 'message')
  sink()
  close(log_con)
}

quiet_ggsave = function(...) {
  invisible(capture.output(ggsave(...), type = 'message'))
}