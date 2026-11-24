<!-- _modified from [EmbeddedArtistry](https://embeddedartistry.com/blog/2017/08/04/a-github-pull-request-template-for-your-projects/)_
_referenced with modifications from [pycytominer](https://github.com/cytomining/pycytominer/blob/master/.github/PULL_REQUEST_TEMPLATE.md)_ -->

# Description

<!--
Thank you for your contribution to amRml!

Please _succinctly_ summarize your proposed change: what motivated it, anything
special you had to do, and related issues (`#<number>` links them).
-->

## Pipeline stage

<!-- Which stage of the amRml pipeline this changes, e.g. "reading the ORB",
"eligibility gate", "matrix build", "fit", "consolidation". -->

## v1.0 requirement served

<!-- Which requirement from the amRml v1.0 design this addresses. -->

## Contract

<!-- What this stage takes in, what it produces (files, columns), and the errors
it raises, with the rule_id of each. -->

## Differences from amRml

<!-- Anything that behaves differently from the current amRml, and why. -->

## Decisions for collaborators

<!-- Open questions this PR needs answered, each phrased with options. -->

## How to try it

```r
# A few lines a reviewer can run.
```

## What kind of change(s) are included?

- [ ] Feature (adds or updates new capabilities)
- [ ] Bug fix (fixes an issue).
- [ ] Enhancement (adds functionality).
- [ ] Breaking change (these changes would cause existing functionality to not work as expected).

# Checklist

Please ensure that all boxes are checked before indicating that this pull request is ready for review.

- [ ] I have read and followed the [CONTRIBUTING.md](CONTRIBUTING.md) guidelines.
- [ ] I have searched for existing content to ensure this is not a duplicate. I have deduplicated my codebase.
- [ ] I have performed a self-review of these additions (including spelling, grammar, and related).
- [ ] I have implemented basic code linting, styling, and indentation (e.g., with lintr, styler).
- [ ] I have added comments to my code to help provide understanding.
- [ ] I have added a test which covers the code changes found within this PR.
- [ ] Every new error uses a `rule_id` declared in `.RULES`, and its test asserts the class, not the message.
- [ ] I have added necessary data (or sample data) to test run these scripts.
- [ ] **Reviewer assignment**: Tag a relevant team member to review and approve the changes.
