# Tests for amRml's error codes and the functions that build and throw errors.

# The error codes

test_that(".RULES has unique rule IDs and known stages", {
  stages <- c("input", "orb", "profile", "feasibility", "eligibility", "matrix", "fit")

  expect_false(anyDuplicated(.RULES$rule_id) > 0)
  expect_true(all(.RULES$stage %in% stages))
})

# Building an error

test_that(".amrError() builds a classed error with its fields, without throwing", {
  expect_silent(cnd <- .amrError(
    "orb_moved",
    c("Headline.", x = "detail"),
    observed = list(recorded = "/a"),
    threshold = list(expected = "/b")
  ))

  expect_s3_class(cnd, c("amrml_orb_moved", "amrml_error", "rlang_error"))
  expect_equal(cnd$rule_id, "orb_moved")
  expect_equal(cnd$stage, "orb")
  expect_equal(cnd$observed, list(recorded = "/a"))
  expect_equal(cnd$threshold, list(expected = "/b"))
})

test_that("an unknown rule_id is an internal error", {
  expect_error(.amrError("no_such_rule", "x"), class = "amrml_internal_error")
  expect_error(.amrError(c("orb_moved", "orb_moved"), "x"), class = "amrml_internal_error")
})

# Throwing an error

test_that(".amrAbort() throws the error and blames the function the user called", {
  user_facing <- function() helper(call = rlang::current_env())
  helper <- function(call) .amrAbort("orb_moved", c("Headline.", x = "detail"), call = call)

  cnd <- rlang::catch_cnd(user_facing(), classes = "error")
  expect_s3_class(cnd, c("amrml_orb_moved", "amrml_error"))
  expect_equal(rlang::call_name(cnd$call), "user_facing")
  expect_match(conditionMessage(cnd), "detail")
})
