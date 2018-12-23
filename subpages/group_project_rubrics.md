---
title: "Group project rubrics"
output:
  html_document:
    includes:
      before_body: ../include/nav.html
      after_body: ../include/nothing.html
    toc: true
---

This page contains more specific requirements for each project deliverable as well as their points breakdown.

## Proposal (5 pts)
A one-page proposal as `project_proposal.md`. Please make sure that you have incorporated the feedback and comments by the professors/TAs on your initial proposal. Your project proposal includes:

**Motivation and background work (1 pt)**

- State the background and rationale for the study.
- What is the overall hypothesis and objectives for the project?
- Why is it important to answer this question? Identify knowledge gaps.

**Division of labour (1 pt)**

- Who is going to do what? State assignment of tasks and projected contributions for each group members.
- Please provide a table of group members with their background, degree, affiliations and job assignments.

**Dataset (1 pt)**

- What kind of data are you working with? What is the general description and characteristics of the data?
- State the technology used to generate the data.
- If your data is small enough, create a read only subdirectory called `data` and keep your data there (see example project repo figure).
- Add a `README.md` file in the subdirectory to describe your data in more detail, e.g., how many rows, how many columns, experimental design etc.

**Aims and methodology (2 pt)**

- What are the specific questions/aims that you will address in order to achieve the final goals of the project?
- For each of the aims, what computational/statistical approaches did you use? (Note: your methodology might change as you are progressing, exploring your data and reading through the literature)
- Provide references for any sections as needed.

Please provide a link to your `project_proposal.md` file in your group repo `README.md`. This is a group-level deliverable so you'd submit it by opening an issue in your group repository, with the issue named "Project proposal of [group name]", and provide the revision SHA and tag the TAs and instructors who are mentors of your group by <span style="color: red">_Feb. 13th, 2019_</span>.

## Progress report (7 pts)

Your group progress report entitled `progress_report.md`.

**What has changed based on the final proposal (1 pt.)**

- Did your dataset change? If so, why?
- Have you decided to do a different analysis than what was mentioned in your proposal? If so, Why?
- Are there any changes in task assignments of group members?

**What is the progress of the analyses (4 pts.)**

- Since your initial proposal, you should have decided more concretely on what methods to use for each step of your analyses, and employed some of those methods.
- Briefly and concisely explain your methodology and progress for the aims you have investigated so far. Which parts were modified and which parts remained the same?
- What R packages or other tools are you using for your analyses? You do not need to provide your scripts in your report.
- Provide the links to any markdown reports within your repo to refer to the relevant analysis. 
- Provide references.

**Results (2 pts.)**
- What are your primary results? Were you able to answer your hypothesis? Did you have any positive results? If no, postulate a discussion as to why that may be. Provide plots and/or tables to present your results.
- List some challenges that you encountered? How will you address them?

Please provide a link to your `progress_report.md` file in your group repo `README.md`. This is a group-level deliverable so you'd submit it by opening an issue in your group repo can call it "Progress report of [group name]"", and provide the revision SHA and tag the TAs and instructors who are mentors of your project by the due date. 

## GitHub repository (20 pts)

This should contain the materials (or associated live links) an instructor would need to evaluate your work and that a group member would need to reproduce/reuse/extend the work. This is a group-level deliverable.

**GitHub organization (10 pts.)**

- Create subdirectories [e.g. src (source code), data, results, docs, etc.] for different parts of your analyses. Provide a `README.md` file for each subdirectory to summarize its contents and analyses.
- The `README.md` file in the root directory of your the group repo should contain the summary of your whole work. This file should summarize the project content, your analysis, and main results. Add links as needed when you are elaborating on your analysis or results. Your project instructors should be able to go through all of your work just by following your `README.md` file. You can include a flow chart if you find it useful. Keep it brief and concise.

Example of a group repository: 

<img src="images/example_group_repo.png" style="width: 70%; height: 70%"/>

**Analysis and results (10 pts)**

- Think of it as a technical report of your presentation. Add brief comments and text throughout explaining the aim and conclusions of each part. You do not need to talk about all the analysis you have done but make sure you are addressing the important files, such as “inputs” (e.g. the dataset(s) you analyzed) and “outputs” (e.g. plain text files of your key statistical results or PDFs of your main figures)
- Make sure your code is readable and could be run by instructors. So, it is necessary to provide valid addresses for the files you use in your scripts. Add comments for different chunks of your code and document it.

You do not need to open an issue and tag TAs for this. Just make sure that your last commit is by the date of the first presentation


## Poster (20 pts)

Your poster as `poster.pdf` should be in your root group repo. For all the things you need for your poster, you can create a subdirectory poster and keep your stuff there.

**Poster session**

- We will have a poster session on Wed. Apr 3th, 9:30 AM - 12 pm in the Atrium of the ESB. This is a group-level deliverable so all group members should be present for the poster presentation.
- You will print your poster in 48’‘(W)x36’’(H) which is the standard conference poster size. Apologies for the cost! If you really need help for that, sometimes supervisors will help especially if your poster could be presented in other academic poster sessions.

**Poster evaluation**

Your poster is evaluated in 5 different categories. You can find the rubric for poster evaluation [here](https://github.com/STAT540-UBC/STAT540-UBC.github.io/blob/master/homework/posterRubric2012.pdf).

Please provide a link to your poster.pdf file in your group repo `README.md`. This is a group-level deliverable so one of your group member should submit it by linking it in an issue with name "Poster of [group name]", and provide the revision SHA and tag the TAs and instructors who are mentors of your project by <span style="color: red">_Apr 4th, 2019_</span>.

## Individual report (8 pts)

A one-page individual report as `individual_report.md` should be in your individual private repository. It includes:
**A concise summary of contributions of each group member (4 pts)**

- Describes the tasks and contributions of each of your group members
- Do you think the task assignments were fairly assigned and appropriate given each member’s background and skills? If no, how would you change it?

**Your specific contributions and comments (4 pts)**

- Explain what are your contributions to the project?
- What worked well and what did not? What was the most challenging or rewarding moment during your group project?
- How did you find having members of different backgrounds for a scientific project? 
- What have you learned from this group project? 
- Any other comments on how the group project could have been structured differently (i.e. delivery requirements, assessment") 

Please provide a link to your `individual_report.md` file in your private repository `README.md`. This is not a group-level deliverable so every student should submit it by linking it in an issue with name Individual report of [your private repository name], provide the revision SHA and tag the TAs by the due date 

Example of a personal repository

<img src="images/example_personal_repo.png" style="width: 70%; height: 70%" align="middle"/>
