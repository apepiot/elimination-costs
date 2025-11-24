#
# Copyright (c) 2001-2011 by Ziena Optimization LLC
# All Rights Reserved.                                 
#
# Example nonlinear optimization problem formulated as an AMPL model.
#

# Define variables and enforce that they be non-negative.
var x{j in 1..3} >= 0;
var y;
# Objective function to be minimized.
minimize obj:

         1000 - x[1]^2 - 2*x[2]^2 - x[3]^2 - x[1]*x[2] - x[1]*x[3];

let y := 8*x[1] + 14*x[2];

# Equality constraint.
s.t. c1: y + 7*x[3] - 56 = 0;

# Inequality constraint.
s.t. c2: x[1]^2 + x[2]^2 + x[3]^2 -25 >= 0;

data;

# Define initial point.
let x[1] := 2;
let x[2] := 2;
let x[3] := 2;
