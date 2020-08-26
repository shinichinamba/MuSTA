get_group_name <- function() system("id -gn", intern = TRUE)
get_user_name <- function() system("id -un", intern = TRUE)
get_group_id <- function() system("id -gr", intern = TRUE)
get_user_id <- function() system("id -ur", intern = TRUE)