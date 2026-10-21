# Tests for planMatrices() and its `modes` choice.

test_that("modes keeps the plan and its members to the modes asked for", {
  eligibility <- eligibleScopes(profileORB(readORB(fxOrb())))
  planned <- unique(planMatrices(eligibility)$matrices$mode_id)
  other <- setdiff(.MODES$mode_id, planned)[[1]]

  kept <- planMatrices(eligibility, modes = planned)
  dropped <- planMatrices(eligibility, modes = other)

  expect_equal(kept$matrices, planMatrices(eligibility)$matrices)
  expect_equal(nrow(dropped$matrices), 0L)
  expect_equal(nrow(dropped$eligibility$members), 0L)
})

test_that("unknown modes are rejected", {
  eligibility <- eligibleScopes(profileORB(readORB(fxOrb())))

  expect_error(planMatrices(eligibility, modes = "nope"), class = "amrml_invalid_argument")
})
