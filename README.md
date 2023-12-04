# Reproducible research: version control and R

### **Questions 1, 2, and 3:** 

https://github.com/ashsilvv/logistic_growth/blob/069146d06eeb7300843f47e5bf6cd5f44c999d30/README.md 

### **Question 4** 
**A) Execute the code to produce the paths of two random walks. What do you observe?**

The "random_walk" function is called twice to generate two sets of random walk data (data1 and data2), and then these datasets are used to create two separate plots (plot1 and plot2). The plots are arranged side by side (with the "grid.arrange" function) so they can be compared.
  
In the "random_walk" function, the initial position of each walk is set at (0,0) at time 1. Then, a random walk is generated using a "for" loop. This iterates a new step of a specific size (0.25) at a random angle from the previous step, at each new timestep. The source of the randomness of each walk is the angle generation: i.e., how the direction of each step is determined. The "runif(1, min = 0, max = 2*pi)" function is used to create a random angle between 0 and 2π for each step. Then, then "x" and "y" coordinates of each new step are calculated using this information in addition to the coordinates of the previous step. 

This "for" loop is repeated for "n_steps" number of times to generate the random walk. Then, when you call "data1 <- random_walk(500)", it sets n_steps to 500: telling the "random_walk" function to generate a dataset for a random walk with 500 steps. This is repeated for each of the two datasets so that both walks are the same length. 

In each of the two plots arranged side by side, we see a path connecting a series of points. Each point corresponds to a position in the two-dimensional space (x, y). Between each point is 1 step, and each walk has 500 steps total. Both paths begin at time step 1 at the coordinates (0,0). The paths have a colour gradient to indicate the timeframe that each step occurred in, where the earliest steps are dark blue and the latest steps are light blue. This helps to understand the progression of the random walk over time. Each step is the same length (0.25) from the previous point, and each step moves forward at a random angle from the previous step, as I previously described. This introduces variability in the direction of movement, so each plot looks different. Having the two walks visualised next to each other helps to illustrate the inherent randomness in the path's trajectory. 

**B) Investigate the term random seeds. What is a random seed and how does it work?**

Random seeds are used to create variables that take on random values, in a way which ensures that the results are reproducible. By setting a seed, you change the underlying state of the random number generator such that it yields the same sequence of pseudo-random numbers every time across different runs (as long as you use the same seed value). You need to set the seed before each call to a random number generation function in order for the numbers to come from the same sequence. Random seeds are extremely important for the reproducibility of code because they ensure that anyone running your code will obtain identical outputs.

The numbers or sequences produced from random seeds are called "pseudo-random" because computers aren't capable of generating truly random numbers, so they use they use algorithms to produce numbers or sequences which appear random. Also, the numeric seed value used to set the seed doesn't actually have any impact on the sequence produced: it's only there so that the sequence can be easily reproduced. 

In R, random seeds can be created using the "set.seed()" function. For example, you can set a random seed with the code "set.seed(69)", then ask R to generate random numbers (e.g., with the "runif()" function). Each time that you set the same seed before asking for random numbers to be generated, the numbers produced will be the same, because they come from the same sequence. 

**C) & D) Edit the script to make a reproducible simulation of Brownian motion. Commit the file and push it to your forked reproducible-research homework repo. Show the edit you made to the code in the comparison view.**

See below. Edits were made directly to the "random_walk.R" code in the "question-4-code" file.
![comparison](https://github.com/ashsilvv/reproducible-research_homework/assets/150149935/5717d64a-c689-487b-9329-9d5660d2998a)


### Question 5

**A) Import the data for double-stranded DNA (dsDNA) viruses taken from the Supplementary Materials of the original paper into Posit Cloud (the csv file is in the question-5-data folder). How many rows and columns does the table have?**

There are 33 rows and 13 columns in this dataset. 

**B) What transformation can you use to fit a linear model to the data? Apply the transformation.**

You can apply a log transformation to both virion volume (y) and genome length (x). This helps to meet one assumption of linear models: that the relationship between X and Y is linear. If this relationship is linear, then data points appear to fall along a straight line in a scatterplot. You can see how I applied a log transformation and checked the linearity of the relationship between X and Y below - see within the "aes()" function:
```
virus <- read.csv("Cui_etal2014.csv")
# reading in the data

install.packages("ggplot2")
library(ggplot2)
#installing packages

###### plotting the relationship:

ggplot(virus, aes(x = log(`Genome.length..kb.`), y = log(`Virion.volume..nm.nm.nm.`))) +
  # feeding in the data and telling R which are X and Y variables
  # also applying a log transformation to both of these!
  geom_point() +
  # telling R we want a scatterplot
  labs(x = "log [Genome length (kb)]",
       y = "log [Virion volume (nm3)]") 
  # making axes labels 
```
![linear](https://github.com/ashsilvv/reproducible-research_homework/assets/150149935/78d7085e-fdc9-4be4-bea3-db82664412e7)


**C) Find the exponent (α) and scaling factor (β) of the allometric law for dsDNA viruses and write the p-values from the model you obtained, are they statistically significant? Compare the values you found to those shown in Table 2 of the paper, did you find the same values?**

The exponent (α) is equal to the regression coefficient, or the slope, of the linear model. In my linear model, this was **1.5152**. This slope is significantly different from 0 (p = 6.44e-10), indicating a statistically significant linear relationship between log(genome length) and log(virion volume).  

The scaling factor (β) is equal to the exponent of the intercept of the linear model (because it needs to be back-transformed from the log scale). In my linear model, this was **1181.8071** (AKA e^7.0748). The intercept was also significantly different from 0 (p = 2.28e-10), which predicts that virion volume (y) has a certain statistically significant baseline value even when genome length (x) is zero, or that virion volume (y) has some statistically significant inherent value which is not accounted for by genome length (x). 

The adjusted R-squared value for this linear model is 0.7042. This indicates that about 70% of the variance in virion volume can be explained by genome length, which is a high percentage. 

The value for the exponent (α) in my linear model is the same as the value found in the paper for dsDNA viruses, when rounded to 2 decimal places (1.52). Likewise, the scaling factor (β) derived from my linear model is also the same as the value found in the paper for dsDNA viruses, when rounded to the nearest whole number (1,182). 

**D) Write the code to reproduce the figure shown below.**

 ```
virus <- read.csv("Cui_etal2014.csv")
head(virus)
#reading in the data and having a look at it

library(ggplot2)
#loading in the necessary packages

##### plotting the relationship between log genome length and log virion volume

ggplot(virus, aes(x = log(`Genome.length..kb.`), y = log(`Virion.volume..nm.nm.nm.`))) +
  # feeding in the data and setting the X and Y variables
  geom_point() +
  # telling R we want a scatterplot
  geom_smooth(method="lm") +
  # fitting a regression line with standard error margins
  labs(x = "log [Genome length (kb)]",
       y = "log [Virion volume (nm3)]") +
  #making axes labels 
  theme_light() + 
# setting a theme
  theme(
  axis.title.x = element_text(size=9, face = "bold"),
  axis.title.y = element_text(size=9, face = "bold"))
  # changing axis label font size and making it bold
 ```

**E) What is the estimated volume of a 300 kb dsDNA virus?**

I used the terms for α and β taken from the summary of my linear model, rather than the rounded terms from the paper, to achieve a more precise estimate. I calculated the volume of a 300kb virus using two different approaches - one using V=βL^α, and one using y=mx+c. I got the same answer with each approach (within 7 nm3 of each other).

V = βL^α

β = 1181.807

L = 300

α = 1.5152

V = 1181.807 * 300 ^ 1.5152 = **6,697,006 nm3**

ALTERNATIVELY: 

y = mx + c 

x = log(300) = 5.703782

m = 1.5152

c = log(1181.807) = 7.0748

y = 1.5152 * 5.703782 + 7.0748 = 15.71717

y = e^15.71717 = **6,696,999 nm3**

### Bonus question

**A) Explain the difference between reproducibility and replicability in scientific research.**

Answer

**B) How can git and GitHub be used to enhance the reproducibility and replicability of your work?**

Answer

**C) What limitations do git and GitHub have?**

Answer

## Instruction

The homework for this Computer skills practical is divided into 5 questions for a total of 100 points (plus an optional bonus question worth 10 extra points). First, fork this repo and make sure your fork is made **Public** for marking. Answers should be added to the # INSERT ANSWERS HERE # section above in the **README.md** file of your forked repository.

Questions 1, 2 and 3 should be answered in the **README.md** file of the `logistic_growth` repo that you forked during the practical. To answer those questions here, simply include a link to your logistic_growth repo.

**Submission**: Please submit a single **PDF** file with your candidate number (and no other identifying information), and a link to your fork of the `reproducible-research_homework` repo with the completed answers. All answers should be on the `main` branch.

## Assignment questions 

1) (**10 points**) Annotate the **README.md** file in your `logistic_growth` repo with more detailed information about the analysis. Add a section on the results and include the estimates for $N_0$, $r$ and $K$ (mention which *.csv file you used).
   
2) (**10 points**) Use your estimates of $N_0$ and $r$ to calculate the population size at $t$ = 4980 min, assuming that the population grows exponentially. How does it compare to the population size predicted under logistic growth? 

3) (**20 points**) Add an R script to your repository that makes a graph comparing the exponential and logistic growth curves (using the same parameter estimates you found). Upload this graph to your repo and include it in the **README.md** file so it can be viewed in the repo homepage.
   
4) (**30 points**) Sometimes we are interested in modelling a process that involves randomness. A good example is Brownian motion. We will explore how to simulate a random process in a way that it is reproducible:

   - A script for simulating a random_walk is provided in the `question-4-code` folder of this repo. Execute the code to produce the paths of two random walks. What do you observe? (10 points)
   - Investigate the term **random seeds**. What is a random seed and how does it work? (5 points)
   - Edit the script to make a reproducible simulation of Brownian motion. Commit the file and push it to your forked `reproducible-research_homework` repo. (10 points)
   - Go to your commit history and click on the latest commit. Show the edit you made to the code in the comparison view (add this image to the **README.md** of the fork). (5 points)

5) (**30 points**) In 2014, Cui, Schlub and Holmes published an article in the *Journal of Virology* (doi: https://doi.org/10.1128/jvi.00362-14) showing that the size of viral particles, more specifically their volume, could be predicted from their genome size (length). They found that this relationship can be modelled using an allometric equation of the form **$`V = \beta L^{\alpha}`$**, where $`V`$ is the virion volume in nm<sup>3</sup> and $`L`$ is the genome length in nucleotides.

   - Import the data for double-stranded DNA (dsDNA) viruses taken from the Supplementary Materials of the original paper into Posit Cloud (the csv file is in the `question-5-data` folder). How many rows and columns does the table have? (3 points)
   - What transformation can you use to fit a linear model to the data? Apply the transformation. (3 points)
   - Find the exponent ($\alpha$) and scaling factor ($\beta$) of the allometric law for dsDNA viruses and write the p-values from the model you obtained, are they statistically significant? Compare the values you found to those shown in **Table 2** of the paper, did you find the same values? (10 points)
   - Write the code to reproduce the figure shown below. (10 points)

  <p align="center">
     <img src="https://github.com/josegabrielnb/reproducible-research_homework/blob/main/question-5-data/allometric_scaling.png" width="600" height="500">
  </p>

  - What is the estimated volume of a 300 kb dsDNA virus? (4 points)

**Bonus** (**10 points**) Explain the difference between reproducibility and replicability in scientific research. How can git and GitHub be used to enhance the reproducibility and replicability of your work? what limitations do they have? (e.g. check the platform [protocols.io](https://www.protocols.io/)).
