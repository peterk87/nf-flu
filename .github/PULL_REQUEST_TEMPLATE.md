<!--
# CFIA-NCFAD/nf-flu pull request

Many thanks for contributing to CFIA-NCFAD/nf-flu!

Please fill in the appropriate checklist below (delete whatever is not relevant).
These are the most common things requested on pull requests (PRs).

Remember that PRs should be made against the dev branch, unless you're preparing a pipeline release.

Learn more about contributing: [CONTRIBUTING.md](https://github.com/CFIA-NCFAD/nf-flu/tree/master/.github/CONTRIBUTING.md)
-->

## PR checklist

- [ ] This comment contains a description of changes (with reason).
- [ ] If you've fixed a bug or added code that should be tested, add tests!
  - [ ] If you've added a new tool - ensure that you've added the version info to a `versions.yml` and added it to the `ch_versions` channel in the workflow you added the tool to.
  - [ ] If you've added a new tool - have you followed the pipeline conventions in the [contribution docs](https://github.com/CFIA-NCFAD/nf-flu/tree/master/.github/CONTRIBUTING.md)
  - [ ] If necessary, also make a PR on [the CFIA-NCFAD/nf-test-datasets repo](https://github.com/CFIA-NCFAD/nf-test-datasets/pull/new)
- [ ] Ensure the test suite passes (`nextflow run . -profile test_{illumina,nanopore},docker`).
- [ ] Usage Documentation in `docs/usage.md` is updated.
- [ ] Output Documentation in `docs/output.md` is updated.
- [ ] `CHANGELOG.md` is updated.
- [ ] `README.md` is updated (including new tool citations and authors/contributors).
