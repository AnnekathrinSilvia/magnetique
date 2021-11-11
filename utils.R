format_url <- function(url, target) {
  if (is.na(target)) {
    return("")
  }
  sprintf("<a target='_blank' href='%s%s'>%s</a>", url, target, target)
}

team_df <- function() {
  data.frame(
    name = c(
      "Annekathrin Ludt",
      "Christoph Dieterich",
      "Enio Gjerga",
      "Etienne Boileau",
      "Federico Marini",
      "Thiago Britto-Borges")
    # email
    # twitter
    # orcid etc
  )
}

