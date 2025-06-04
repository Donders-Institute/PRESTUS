# Contributing to PRESTUS

Thank you for your interest in contributing to the open-source development of PRESTUS!

We welcome all types of contributions, including code submissions, documentation, tests, and feedback.

PRESTUS emerged from the need for tailored simulation solutions in the realm of transcranial ultrasound stimulation (TUS) research. As TUS research advances, new equipment and experimental paradigm requirements continually arise, often outpacing the capabilities of existing simulation tools. PRESTUS was developed to address these challenges by providing a flexible set of functions for ultrasonic simulations. This adaptability allows researchers to efficiently prototype, test, and refine simulation workflows that are responsive to the unique and evolving demands of their studies.

Your contributions help maintain and usefully expand this evolving tool for the non-invasive brain stimulation (NIBS) community.

## Getting Started
- Fork the repository and clone it locally.
- Install dependencies as described in [README.md](README.md).
- Please review open issues and discussions to avoid duplicating work and to align with ongoing development.

## Coding Standards

PRESTUS is a research tool, and we strive for transparency, reproducibility, and clarity—even if the workflows are evolving and at times built for specific use cases in scientific work. Please consider the following guiding principles:

- **Motivation:** When proposing changes, open an issue or comment on an existing one to describe the scientific or technical motivation behind your contribution.
- **Documentation:** Clearly document any new features, changes, or known limitations directly in the code and in your pull request description. If your contribution addresses messy or ambiguous data, describe your approach and any assumptions made.
- **Reproducibility:** Strive to make your changes as reproducible as possible. If you use external data, scripts, or parameters, provide instructions or references so others can replicate your results.
- **Transparency:** If your solution involves workarounds, non-standard methods, or "messy" scientific practices, explain these choices in your code comments and pull request. This helps future users and contributors understand the context and rationale.
- **Code Quality:** Write clear, descriptive commit messages. Limit the impact of your updates to the target features/fixes, and minimize unintended effects on other workflows or users.
- **Backward Compatibility:** Remain backward-compatible where possible. If not feasible, document the necessary changes and impacts in your code, issue, and pull request.
- **Data Handling:** If your contribution involves new datasets or changes to data processing, include metadata and a description of data provenance, quality, and any preprocessing steps.

## How to Submit Changes
- Create a new branch for your feature or fix.
- Ensure (to the best of your ability) that your features work as expected and the pipeline remains functional. Unit tests have not yet been implemented, so please verify your changes manually and describe your testing process.
- Submit a pull request to the `development` branch with a clear description. Before this pull request can be accepted, it will need to be reviewed by another contributor.
- [Optional] Attach a log or summary to the pull request, stating that the "pipeline has finished successfully" or describing any issues encountered.
- [Optional] Add yourself to the contributor list in the readme (you will automatically be recognized by GitHub's contribution tracker). Follow [these instructions](https://allcontributors.org/docs/en/bot/usage).

## Code of Conduct
Please read our [Code of Conduct](CODE_OF_CONDUCT.md) before contributing.

## Contact
If you have questions, please [open an issue](https://github.com/Donders-Institute/PRESTUS/issues).
