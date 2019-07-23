Exercises for LSN, 2019

## IMPORTANT! ## 

## PLEASE, DO NOT RUN THE NOTEBOOKS IF NOT NECESSARY. ##

## Some of the exercises have in them a "os.system" command, which executes directly ##
## C++ programms. While very useful to me while I was doing them, running everything ##
## could take sometimes more then 15min. ##

I decided to operate this way for two main reasons:
  • I could have everything under control in one enviroment.
  • For very computationally demanding programs, it's a way to show that I "didn't cheat": as a final test, I alwyas run 
    the full notebook one last time, which means I cannot "copy" the graphs from someone else. 
  • I find it pretty cool.


## COMMENTS ##

  • Ex_4: in some of the graphs in here, the y axis cuts a small part of the data. Due to the size of the exercise and
    the very small nature of the problem, I didn't want to manually adjust the axix where the problem accured.
    
  • Ex_8: the exercise hasn't been done in the most efficient way. I decided to run two "for" loops in the notebook instead
    of the C++ program just for learning purposes. While I know that running the lopps inside the main program would have taken
    a few seconds instead of minutes, I wanted to test the capabilities (and limits) of Python.
    
  • Ex_10: in order to have uncorrelated rundom numbers between the different parallelized runs, I decided to built a random
    number generator function in the "classes.h" file, which modifies the seeds according to the rank.
    While I cannot say that it's 100% perfect, I hope that it's a good alternative to what was suggested.


## DISCLAIMER ##  
I runed every single notebook one last time before uploading it to the repository. 
While I can confirm that everything worked perfectly on my local system (MacBook Pro mid-2014),
I cannot confirm that it will not have problems on another systems, e.g. the results are different
or it just doesn't run.

## ERRORS ##
LSN_exercise_07 has a quite big mistake in the calculation for the function $g(r)$.

©Leonardo Alchieri, 2019
