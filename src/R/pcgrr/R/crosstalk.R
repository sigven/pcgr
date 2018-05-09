

jqueryLib_ct <- function() {
  htmltools::htmlDependency(
    "jquery", "1.11.3",
    system.file("lib/jquery", package = "crosstalk"),
    script = "jquery.min.js"
  )
}

ionrangesliderLibs_ct <- function() {
  list(
    jqueryLib_ct(),
    htmltools::htmlDependency("ionrangeslider", "2.1.2",
                   system.file("lib/ionrangeslider", package = "crosstalk"),
                   script = "js/ion.rangeSlider.min.js",
                   # ion.rangeSlider also needs normalize.css, which is already included in
                   # Bootstrap.
                   stylesheet = c("css/ion.rangeSlider.css",
                                  "css/ion.rangeSlider.skinShiny.css")
    ),
    htmltools::htmlDependency("strftime", "0.9.2",
                   system.file("lib/strftime", package = "crosstalk"),
                   script = "strftime-min.js"
    )
  )
}

hasDecimals_ct <- function(value) {
  truncatedValue <- round(value)
  return (!identical(value, truncatedValue))
}

controlLabel_ct <- function(controlName, label) {
  if (is.null(label)) {
    NULL
  } else {
    shiny::tags$label(class = "control-label", `for` = controlName, label)
  }
}


# Given a vector or list, drop all the NULL items in it
dropNulls_ct <- function(x) {
  x[!vapply(x, is.null, FUN.VALUE=logical(1))]
}

# Format a number without sci notation, and keep as many digits as possible (do
# we really need to go beyond 15 digits?)
formatNoSci_ct <- function(x) {
  if (is.null(x)) return(NULL)
  format(x, scientific = FALSE, digits = 15)
}

#' Range filter control
#'
#' Creates a slider widget that lets users filter observations based on a range
#' of values.
#'
#' @param id An HTML element ID; must be unique within the web page
#' @param label A human-readable label
#' @param sharedData \code{SharedData} object with the data to filter
#' @param column A one-sided formula whose values will be used for this slider.
#'   The column must be of type \code{\link{Date}}, \code{\link{POSIXt}}, or
#'   numeric.
#' @param step Specifies the interval between each selectable value on the
#'   slider (if \code{NULL}, a heuristic is used to determine the step size). If
#'   the values are dates, \code{step} is in days; if the values are times
#'   (POSIXt), \code{step} is in seconds.
#' @param round \code{TRUE} to round all values to the nearest integer;
#'   \code{FALSE} if no rounding is desired; or an integer to round to that
#'   number of digits (for example, 1 will round to the nearest 10, and -2 will
#'   round to the nearest .01). Any rounding will be applied after snapping to
#'   the nearest step.
#' @param ticks \code{FALSE} to hide tick marks, \code{TRUE} to show them
#'   according to some simple heuristics.
#' @param animate \code{TRUE} to show simple animation controls with default
#'   settings; \code{FALSE} not to; or a custom settings list, such as those
#'   created using \code{\link{animationOptions}}.
#' @param width The width of the slider control (see
#'   \code{\link[htmltools]{validateCssUnit}} for valid formats)
#' @param sep Separator between thousands places in numbers.
#' @param pre A prefix string to put in front of the value.
#' @param post A suffix string to put after the value.
#' @param dragRange This option is used only if it is a range slider (with two
#'   values). If \code{TRUE} (the default), the range can be dragged. In other
#'   words, the min and max can be dragged together. If \code{FALSE}, the range
#'   cannot be dragged.
#' @param timeFormat Only used if the values are Date or POSIXt objects. A time
#'   format string, to be passed to the Javascript strftime library. See
#'   \url{https://github.com/samsonjs/strftime} for more details. The allowed
#'   format specifications are very similar, but not identical, to those for R's
#'   \code{\link{strftime}} function. For Dates, the default is \code{"\%F"}
#'   (like \code{"2015-07-01"}), and for POSIXt, the default is \code{"\%F \%T"}
#'   (like \code{"2015-07-01 15:32:10"}).
#' @param timezone Only used if the values are POSIXt objects. A string
#'   specifying the time zone offset for the displayed times, in the format
#'   \code{"+HHMM"} or \code{"-HHMM"}. If \code{NULL} (the default), times will
#'   be displayed in the browser's time zone. The value \code{"+0000"} will
#'   result in UTC time.
#'
#' @examples
#' ## Only run examples in interactive R sessions
#' if (interactive()) {
#'
#' sd <- SharedData$new(mtcars)
#' filter_slider("mpg", "Miles per gallon", sd, "mpg")
#'
#' }
#' @export
filter_slider_allelic_fraction <- function(id, label, sharedData, column, step = NULL,
                                           round = FALSE, ticks = TRUE, animate = FALSE, width = NULL, sep = ",",
                                           pre = NULL, post = NULL, timeFormat = NULL,
                                           timezone = NULL, dragRange = TRUE)
{
  # TODO: Check that this works well with factors
  # TODO: Handle empty data frame, NA/NaN/Inf/-Inf values

  if (is.character(column)) {
    column <- lazyeval::f_new(as.symbol(column))
  }

  df <- sharedData$data(withKey = TRUE)
  col <- lazyeval::f_eval(column, df)
  values <- na.omit(col)
  min <- 0
  max <- 1
  value <- range(values)

  ord <- order(col)
  options <- list(
    values = col[ord],
    keys = df$key_[ord],
    group = sharedData$groupName()
  )

  # If step is NULL, use heuristic to set the step size.
  findStepSize <- function(min, max, step) {
    if (!is.null(step)) return(step)

    range <- max - min
    # If short range or decimals, use continuous decimal with ~100 points
    if (range < 2 || hasDecimals_ct(min) || hasDecimals_ct(max)) {
      step <- pretty(c(min, max), n = 100)
      step[2] - step[1]
    } else {
      1
    }
  }

  if (inherits(min, "Date")) {
    if (!inherits(max, "Date") || !inherits(value, "Date"))
      stop("`min`, `max`, and `value must all be Date or non-Date objects")
    dataType <- "date"

    if (is.null(timeFormat))
      timeFormat <- "%F"

  } else if (inherits(min, "POSIXt")) {
    if (!inherits(max, "POSIXt") || !inherits(value, "POSIXt"))
      stop("`min`, `max`, and `value must all be POSIXt or non-POSIXt objects")
    dataType <- "datetime"

    if (is.null(timeFormat))
      timeFormat <- "%F %T"

  } else {
    dataType <- "number"
  }

  step <- findStepSize(min, max, step)
  # Avoid ugliness from floating point errors, e.g.
  # findStepSize(min(quakes$mag), max(quakes$mag), NULL)
  # was returning 0.01999999999999957 instead of 0.2
  step <- signif(step, 14)

  if (dataType %in% c("date", "datetime")) {
    # For Dates, this conversion uses midnight on that date in UTC
    to_ms <- function(x) 1000 * as.numeric(as.POSIXct(x))

    # Convert values to milliseconds since epoch (this is the value JS uses)
    # Find step size in ms
    step  <- to_ms(max) - to_ms(max - step)
    min   <- to_ms(min)
    max   <- to_ms(max)
    value <- to_ms(value)
  }

  range <- max - min

  # Try to get a sane number of tick marks
  if (ticks) {
    n_steps <- range / step

    # Make sure there are <= 10 steps.
    # n_ticks can be a noninteger, which is good when the range is not an
    # integer multiple of the step size, e.g., min=1, max=10, step=4
    scale_factor <- ceiling(n_steps / 10)
    n_ticks <- n_steps / scale_factor

  } else {
    n_ticks <- NULL
  }

  sliderProps <- dropNulls_ct(list(
    `data-type` = if (length(value) > 1) "double",
    `data-min` = formatNoSci_ct(min),
    `data-max` = formatNoSci_ct(max),
    `data-from` = formatNoSci_ct(value[1]),
    `data-to` = if (length(value) > 1) formatNoSci_ct(value[2]),
    `data-step` = formatNoSci_ct(step),
    `data-grid` = ticks,
    `data-grid-num` = n_ticks,
    `data-grid-snap` = FALSE,
    `data-prettify-separator` = sep,
    `data-prefix` = pre,
    `data-postfix` = post,
    `data-keyboard` = TRUE,
    `data-keyboard-step` = step / (max - min) * 100,
    `data-drag-interval` = dragRange,
    # The following are ignored by the ion.rangeSlider, but are used by Shiny.
    `data-data-type` = dataType,
    `data-time-format` = timeFormat,
    `data-timezone` = timezone
  ))

  # Replace any TRUE and FALSE with "true" and "false"
  sliderProps <- lapply(sliderProps, function(x) {
    if (identical(x, TRUE)) "true"
    else if (identical(x, FALSE)) "false"
    else x
  })

  sliderTag <- htmltools::div(
    class = "form-group crosstalk-input",
    class = "crosstalk-input-slider js-range-slider",
    id = id,

    style = if (!is.null(width)) paste0("width: ", htmltools::validateCssUnit(width), ";"),
    if (!is.null(label)) controlLabel_ct(id, label),
    do.call(shiny::tags$input, sliderProps),
    shiny::tags$script(type = "application/json",
                `data-for` = id,
                jsonlite::toJSON(options, dataframe = "columns", pretty = TRUE)
    )
  )

  # Add animation buttons
  if (identical(animate, TRUE))
    animate <- shiny::animationOptions()

  if (!is.null(animate) && !identical(animate, FALSE)) {
    if (is.null(animate$playButton))
      animate$playButton <- shiny::icon('play', lib = 'glyphicon')
    if (is.null(animate$pauseButton))
      animate$pauseButton <- shiny::icon('pause', lib = 'glyphicon')

    sliderTag <- htmltools::tagAppendChild(
      sliderTag,
      shiny::tags$div(class='slider-animate-container',
               shiny::tags$a(href='#',
                      class='slider-animate-button',
                      'data-target-id'=id,
                      'data-interval'=animate$interval,
                      'data-loop'=animate$loop,
                      span(class = 'play', animate$playButton),
                      span(class = 'pause', animate$pauseButton)
               )
      )
    )
  }

  htmltools::browsable(htmltools::attachDependencies(
    sliderTag,
    c(ionrangesliderLibs_ct(), crosstalk::crosstalkLibs())
  ))
}
