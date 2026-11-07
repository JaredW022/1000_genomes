# Acknolwedgement of Agentic AI use

Codex was used in the production of this analysis, aiding in the writing and testing of code, organization, and final review. I chose the methods of analysis, selection and presentation of figures, set up measures regarding privacy of data. All proposed changes were reviewed prior to implementaiton.

Below is a more detailed listing of AI use. Prompts are paraphrased for clarity. 

## Notable prompts

### Project setup and scope

Initial setup included defining what my the agent should do. Quality and clarity of output was addressed though setting of limitations and expectations for the Agent.

| Prompt | Outcome |
| --- | --- |
| Explain proposed changes and wait for approval before editing code. | I defined rules for code review, notice on large compute jobs, brought awareness to data privacy, for andcommits. |

### Part 1 · cis-eQTL analysis

To identify cis-eQTL within the +/- 500kb window, a workflow plan was drafted. From this, I selected a model to fit the associations, and selected example plots which communicate the outputs.

| Prompt | Outcome |
| --- | --- |
| Draft the initial cis-eQTL workflow. | Codex drafted R code to match samples, define cis windows, and test SNP–gene pairs. |
| Use `lm()` and report the strongest association per gene. | I chose ordinary linear regression, with the metric of lowest raw p-value per tested gene. Codex fitted this model. |
| Plot the IL17RA association and the surrounding region. | Codex aided in the syntax for creating the selected boxplot and association regression graph. |
| Can readers reproduce the figures by knitting the notebook? | After my approval, Codex fixed a error in the reproduction step of our analysis plots. |

### Part 2 · Height polygenic risk score

Similarly, a workflow was devised for calculating polygenic scores for height.

| Prompt | Outcome |
| --- | --- |
| Walk through the scoring workflow before writing code. | We discussed GWAS filtering, position matching, and allele alignment. Codex explained how LD clumping selects lead variants and how the score weights effect-allele counts. |
| Set up the notebook and make a simple histogram. | Codex helped plot the sample against the reference distribution and checked the rendered figure and PNG export. |

### Reproduction instructions

Utilized Codex to ensure that all required tools to replicate our results were included. 

| Prompt| Outcome |
| --- | --- |
| Consolidate package setup | Codex grouped dependencies and downloads. I decided to direct users to each sections respective Rmd file for replication of the graphs. |
