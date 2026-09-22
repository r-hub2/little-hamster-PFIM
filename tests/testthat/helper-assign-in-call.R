# Named first argument on expect_warning / expect_error is not assignment
# (testthat 3.2.3: Can't specify `...` without `pattern`). Use <- inside the call.
# Replacement prop(x, "a") = v is not this pattern.

.pfimExpectNamedAssignRe = paste0(
  "expect_(warning|error|message|condition|no_warning|silent)",
  "\\s*\\(\\s*\\{?\\s*[A-Za-z.][A-Za-z0-9.]*\\s*="
)

.pfimTestFilesWithNamedAssign = function( files ) {
  hits = vapply( files, function( f ) {
    txt = paste( readLines( f, warn = FALSE ), collapse = "\n" )
    grepl( .pfimExpectNamedAssignRe, txt, perl = TRUE )
  }, logical( 1L ) )
  files[ hits ]
}
