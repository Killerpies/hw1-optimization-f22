## CPSC 471/510 HW 1 - Gradient-based Optimization Demo



Run the python code from PyCharm with working directory set to top-level project, or
from command prompt with working directory set to top-level project:

`PYTHONPATH=. python src/optimization.py`

Requires `python3`, `numpy`, and `matplotlib`

This runs a script that will display a 3D surface *z=f(x, y)* as both 3D surface and
contour plot.

See lecture notes for math details

This is one of the few times I'm going to ask you to actually take a derivative,
and write some code to implement that derivative directly.

If you run the script now, it will create a random surface with constant term (0th term)
and at least one additional term using a product of `cosine` and `sine` terms.

Again see lecture notes for math details.
This version also has some *barrier* function terms to keep things in bounds.
I have provided the derivative components for those.

Each time you run the script, it will loads a set of terms, and a starting location from files in the `data` folder.

Use these to start, but if you define a new file name, it will create a random
set up.

After completing the code to apply the derivative, run the code several times with
 different parameters for the alpha and momentum terms, and both starting points given.
Watch https://www.youtube.com/watch?v=k8fTYJPd3_I for details on the momentum terms.


  > NOTE:  For initial testing of gradient calculation, you may want to start with alpha = 0.05 and mom=0.0

The new code you need to write is in `optimization.py`:
1) On Line 408, replace 'first.last.yy' with your user name (email)
2) In the `def calc_gradient_term(position, terms_list):` starting on Line 160,
  you need to correct the lines 178 and 179 to properly calculate the gradient terms.


After writing and testing the correct gradient code, run the script with different `alpha` values of with 0.9, 0.5, 0.05, 0.005, and `mom` of 0.0, 0.8 with the `two_terms.txt` file,
or with either of the other files and max terms set to 2.  
Note the impact on values and number of steps required.

This should require at least 16 runs of the script after getting the gradient correct.
(4 `alpha` *  2 `mom` * 2 `start_points`).

Choose a reasonable `alpha` and `mom` combination, set max terms = 4, and run
with the other term files, and both starting points.  Another 6 runs.

Then try a couple of random start positions in with the 3 term file.
This should require about 22 runs in total.

Each run saves an image of the contour plot in the `plots` folder.  
Choose four of the more interesting plots, and upload to Scholar with
your written solutions to HW1 in one submission; only the latest submission is graded.


* * *

### Submission

You're not done yet!
Make sure your modified code is pushed to gitlab.
Make sure images are submitted with HW 1 to Scholar.
