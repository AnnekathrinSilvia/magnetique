format_url <- function(url, target) {
  if (is.na(target)) {
    return("")
  }
  sprintf("<a target='_blank' href='%s%s'>%s</a>", url, target, target)
}


