* changes for gradient model matter a lot more than in the original model
* smoother fitness landscape
* easier to optimize
* basically change rate matters way more with gradient
* probably not an effect of mutation rate

Moving forward:
* tells us that the way we're changing the environment with gradient mimics large change rates
* change rate doesn't mean the same thing across both models
* can't answer original question
* gradient model does make a difference with respect to change rate
* to answer the original question, need to change the axes
Todo: compare graphs we have but pretending the stuff in the low change rate gradient model to high change rate non-gradient model
* genes spread out more with environmental change whether gradient or non-gradient
* increase in mutation rate at static environment decreases number of coding sites significantly
* qualitatively, results are the same, just noisier

* Soooo turns out that mutation rate was .01 for environmental change runs, not .001 like for mutation rate runs

* TODO: fix .01 to .001 error, find correct mutation rate - in what file??
* TODO: fix change rate issue to be smaller, change every certain number of updates, or probabalistic
* TODO: probabilisitc if < 1, and then that's - Scratch that
* TODO: see old model at higher change rates, see new model at lower change rates
* for lower change rates in gradient model: change 1 bit every 2, every 5, every 10 ,every 100, every 1000 generations, etc.

* how many updates should we wait between changes: 10,000, 5,000, 1,000, 500, 100, 50, 25, 10, 5, 2 - only need to run with gradient model

* TODO: which model do we want to use?
* tempted to use the gradient data, especially if it's cleaner data
* gradient and non-gradient are polar opposites, if they have same behavior then don't have to worry what we pick, probably pick gradient cause easier to explain

* root of problem is that in gradient model theres 128 bits in model =/= 5,000 something bits in non-gradient model
* could do is inject neutrality into gradient model's environment to get number of bits closer to 5,000
* but goal is not to compare