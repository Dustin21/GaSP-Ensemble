# Create GaSP job file

gaspfit <- function(seed = 100, xdescrip = "NULL", ydescrip = "NULL", regmod = "constmod.mtx",
corfn = "PowerExponential", RandomError = FALSE, 
   tries = 10,
   CrossValidate = FALSE, Predict = TRUE, Validate = FALSE, ytrue = "ytrue.mtx", Visualize = FALSE, platform = "unix", m) 
{ 
   fn = paste(paste("fit", m, sep = ""), "gsp", sep = ".")
    
   cat(file = fn, "RandomNumberSeed = ", seed, "\n")

   cat(file = fn, append = TRUE, "#Input Files:\n")
   if (xdescrip != "NULL")
      cat(file = fn, append = TRUE, "XDescription <", xdescrip, "\n")
   cat(file = fn, append = TRUE, "X <", paste0("xc", m, ".mtx"), "\n")
   if (ydescrip != "NULL")
      cat(file = fn, append = TRUE, "YDescription <", ydescrip, "\n")
   cat(file = fn, append = TRUE, "Y <", paste0("yc", m, ".mtx"), "\n")

   cat(file = fn, append = TRUE, "RegressionModel <", regmod, "\n")
   if (corfn == "Gauss")
   {
      cat(file = fn, append = TRUE, "CorrelationFamily = PowerExponential\n")
      cat(file = fn, append = TRUE, "Alpha.Max = 0\n")
      cat(file = fn, append = TRUE, "Alpha.Min = 0\n")
   }
   else if (corfn == "Matern2")	
   { 
      cat(file = fn, append = TRUE, "CorrelationFamily = Matern\n")
  	cat(file = fn, append = TRUE, "Derivatives.Min = 2\n")
      cat(file = fn, append = TRUE, "Derivatives.Max = 2\n")	
   }
   else if (corfn == "Matern")	
      cat(file = fn, append = TRUE, "CorrelationFamily = Matern\n")
   else if (corfn == "Exponential")
   {
      cat(file = fn, append = TRUE, "CorrelationFamily = PowerExponential\n")
      cat(file = fn, append = TRUE, "Alpha.Max = 0\n")
      cat(file = fn, append = TRUE, "Alpha.Min = 0\n")
   }
   else if (corfn == "PowerExponential")	
      cat(file = fn, append = TRUE, "CorrelationFamily = PowerExponential\n")
 
   if (RandomError)
      cat(file = fn, append = TRUE, "RandomError = Yes\n")
   else
      cat(file = fn, append = TRUE, "RandomError = No\n")

   cat(file = fn, append = TRUE, "# MLE\n")
   cat(file = fn, append = TRUE, "CriticalLogLikelihoodDifference = 0\n")
   cat(file = fn, append = TRUE, "Tries  = ", tries, "\n", sep = "")
   cat(file = fn, append = TRUE, "Fit\n")
   cat(file = fn, append = TRUE, "#RegressionModel > beta.mtx\n")
   cat(file = fn, append = TRUE, "#StochasticProcessModel > corpar.mtx\n")
 
   if (Visualize)
   {
      cat(file = fn, append = TRUE, "MainEffectPercentage = 0\n")
      cat(file = fn, append = TRUE, "InteractionEffectPercentage = 5\n")
      cat(file = fn, append = TRUE, "Visualize\n")
      cat(file = fn, append = TRUE, "ANOVAPercentages > anova.mtx\n")	
      cat(file = fn, append = TRUE,  "MainEffects > maineff.mtx\n")
      cat(file = fn, append = TRUE,  "JointEffects > jointeff.mtx\n")    
   }
 
   if (CrossValidate)
   {
      cat(file = fn, append = TRUE, "CrossValidate\n")
      cat(file = fn, append = TRUE, "#CrossValidations > cv.mtx \n")
   }
   
   if (Validate)
   { 
      cat(file = fn, append = TRUE, "YTrue <", ytrue, "\n") 
   }		
   if (Predict)
   { 
      cat(file = fn, append = TRUE, "XPrediction < ", "xp.mtx\n")
     cat(file = fn, append = TRUE, "Predict\n");
     cat(file = fn, append = TRUE, "YPrediction > ", paste(paste("ypr", m, sep=""), "mtx\n", sep="."))  }			
cat(file = fn, append = TRUE, "# Output files:\n")
cat(file = fn, append = TRUE, "#YDescription > sm.mtx \n")
cat(file = fn, append = TRUE, "Quit\n")
cat()
}		
