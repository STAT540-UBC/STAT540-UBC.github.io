---
title: "Assignment Tips"
output:
  html_document:
    includes:
      before_body: ../include/nav.html
      after_body: ../include/nothing.html
    toc: true
---

This is for tips on the assignment, but also useful for seminars and any future project you might have. 

### When working on the assignment  

*Start early*. Even if you are already fluent with all the seminars materials, it'd still take time to answer all the questions. 

When you get stuck or when you run into an error, ask yourself these questions: 

- Am I in the right working directory ? 
- Is this material covered in one of the seminars?
- Can I google the information to find how to do this? 
- Is there an R package that can more efficiently do what I'm attempting now? 
- Am I using the right parameters for this function?  (hint, type `?function_name()` where function_name is the name of the function, to check if you're using it right) 

Your Rmarkdown can't be knitted? Is it because you didn't define all the variables inside your Rmarkdown? (your code chunk might run in Rstudio if all the variables are all defined in the Rstudio environment, but it wouldn't run as the Rmarkdown is being knitted if the variables aren't defined in the Rmarkdown document.)

Well-structured dataframes and reusable functions can often lighten your work load. 

Overall presentation and mechanics refer to the fluency, neatness, easiness to read. For example:

- Using headings/subheadings to distinguishes different sections/questions
- Explain what you're doing to show your understanding, ie sandwiching your code and result with some explanations and interpretation. We don't want to see just a graph or some R outputs standing alone in a question. 
- Use inline R code whenever you refers to the value of a variable in a block of text. 
- Hide useless messages or warnings (though make sure the warnings are harmless) using code chunk options. But don't hide any code that generates graphs or are critical steps in your analysis, obviously. 
- Comment on your code so that everyone, including yourself, can easily follow through your steps.
    - In R, comment follows the number sign. For example `# What I did here`  
- There are a few occasions where, instead of just printing an object with R, you could format the info in an nice table. 
    - The `kable()` function from `knitr` package.
    - Also look into the packages `xtable`, `pander` for making pretty HTML tables.

You might find [cheatsheets from Rstudio](https://www.rstudio.com/resources/cheatsheets/) useful, in terms of graphing, making an awesome R markdown, etc.

### Make it easy for people to access your work 

Reduce the friction for TAs and profs to get the hard-working source code and commentary (the R markdown) __and__ the markdown that has all your code and outputs. 

To create the markdown file from Rmarkdown, set the output of the Rmarkdown to "github_document".

```
---
title: "Homework assingment"
author: "Santina Lin"
date: "February 7, 2017"
output: github_document
---

```
When you submit your homework, knit your Rmarkdown into markdown. Commit and push all of the following:
- The Rmarkdown file
- The markdown file that's created
- The folder that's created. It contains your figures.  


### Make it easy for others to run your code

  * In exactly one, very early R chunk, load any necessary packages, so your dependencies are obvious.
  * In exactly one, very early R chunk, import anything coming from an external file. This will make it easy for someone to see which data files are required, edit to reflect their locals paths if necessary, etc. There are situations where you might not keep data in the repo itself.
  * Pretend you are someone else. Clone a fresh copy of your own repo from GitHub, fire up a new RStudio session and try to knit your R markdown file. Does it "just work"? It should!
  