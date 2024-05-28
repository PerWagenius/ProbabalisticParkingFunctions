import numpy as np
import matplotlib.pyplot as plt
import random

import itertools
from sympy import *
import seaborn as sns

def is_parking(n,k,p,s, boundary = "HO-HO"):
    """
    Function that runs one simulation to see whether the sequence s of length n is a parking function. Parameters are:
        
        n = number of car/spots. Spots are counted starting from 1 up to n.
        k = Number of coin flips.
        p = Bernoulli parameter of coin flips.
        s = The preference vector, a list of n numbers between 1 and n
        boundary = The boundary conditions. Options are
            "HO-HO" = "Hold on" conditions on both endpoints (cars don't leave on a bad flip)
            "LDC-LDC" = "Leave and don't comeback" conditions on both endpoints (cars leave on a bad flip)
            "HO-LDC" = HO condition on the left and LDC condition on the right
            "LDC-HO" = LDC condition on the left and HO condition on the right
            "circle" = Boundary condition is equivalent to going mod number of spots
        Boundary conditions are assumed to be "Hold on"

        TODO: Add option for other boundary conditions
    """
    if min(s) not in range(1,n + 1) or max(s) not in range(1, n + 1):
        return False
    parked_spots = [] # keeps track of spots where a car has parked
    for i in range(0,n): # each part of this loop is the trajectory of a car
        if s[i] not in parked_spots: # if spot is available, the car can park
            parked_spots += [s[i]]
        else:
            current_spot = s[i] # if not, we need to keep track of the spot where the car is at each step
            for j in range(k): # each car has k coin flips
                flip = random.uniform(0,1) # flip is the "coin flip" variable: flip >= p means tails, flip < p means heads 
                if flip >= p and current_spot != 1: # if coin flip is tails and the car is not at the left boundary, the car goes back one spot
                    current_spot = current_spot - 1
                elif flip >= p and current_spot == 1: # case where coin flip is tails and the car is exactly at the left boundary
                    if boundary == "HO-HO" or boundary == "HO-LDC": # if left boundary condition is HO, car stays put
                        current_spot = current_spot
                    elif boundary == "LDC-HO" or boundary == "LDC-LDC": # if left boundary condition is LDC, it is not a parking function
                        return False
                    elif boundary == "circle": # if boundary condition is circle, car goes to the right end spot
                        current_spot = n
                elif flip < p and current_spot < n: # if coin flip is heads and the car is not at the right boundary, the car goes forward one spot
                        current_spot = current_spot + 1
                elif flip < p and current_spot == n: # case where coin flip is heads and the car is at the right boundary
                    if boundary == "HO-HO" or boundary == "LDC-HO":
                        current_spot = current_spot # if right boundary is HO, car stays put
                    elif boundary == "HO-LDC" or boundary == "LDC-LDC": # if right boundary is LDC, it is not a parking function
                        return False 
                    elif boundary == "circle": #if boundary condition is circle, car goes to the left end spot
                        current_spot = 1
                
                if current_spot in parked_spots: # after moving the car, we need to check whether the current spot is available or not. 
                                                 #If it isn't, we need to flip the coin again
                    continue
                else: # if it is available, park the car. We add the current spot to our list of parked spots
                      # and may go to the next car by breaking the current (j) loop 
                    parked_spots += [current_spot]
                    break
    if len(parked_spots) == n: # we have a parking function if and only if every car parked, so our list of spots with parked cars should have size 
                               # equal to the number of cars
        return True
    else: # otherwise it is not a parking function
        return False

def sample_parking(n,k,p,s, sample_size, boundary):
    """
        Function that runs "sample_size" simulations in order to estimate the probability that a given sequence is a parking function. Parameters are:

        n = number of car/spots. Spots are counted starting from 1 up to n.
        k = Number of coin flips.
        p = Bernoulli parameter of coin flips.
        s = The preference vector, a list of n numbers between 1 and n.
        sample_size = number of times the simulation will run

        Boundary conditions are assumed to be "Hold on"

        TODO: Add option for other boundary conditions
    """
    c = 0 # We need a variable that is going to count the number of times we got an outcome of parking function.
    for i in range(0,sample_size): # this loop simply runs the "one time" simulation "sample_size" times, over and over
        if is_parking(n, k, p, s, boundary): # if the outcome is that "s" is a parking function, we add one to our counter
            c += 1 
        else: # otherwise, run it again
            continue
    return c/sample_size # return number of "true" outcomes divided by total number of simulations

def kpi_curve(n,p,s, sample_size, sample_times, krange, boundary):
    """
        Function that computes the kpi curve, with x coordinates being the number k of coin flips,
        and y coordinate being the probability that s is a parking function. Parameters are:

        n = number of car/spots. Spots are counted starting from 1 up to n.
        k = Number of coin flips.
        p = Bernoulli parameter of coin flips.
        s = The preference vector, a list of n numbers between 1 and n.
        sample_size = number of times the simulation will run
        sample_times = number of separate times, the simulation will run. In total, simulation runs sample_size * sample_times times
        krange = domain of x axis for number of coin flips

        TODO: Add option for other boundary conditions
    """
    probs = [] # list of y coordinates of the curve. Starts empty and we add values to it one at a time
    for k in krange: # for every value in the krange, we run it once. Note that by definition this curve is only defined on integer points (number of coin flips).
                     # We can still plot a "continuous" curve in the end by interpolating points.
        pk = 0 # pk is the y coordinate. Since we are running the (sample_size) simulation sample_times times, we need to take an average in the end
        for j in range(0,sample_times): # add up every probability for each sample_time result
            pk += sample_parking(n,k,p,s, sample_size, boundary)
        pk = pk/sample_times # take the average
        probs += [pk] # add y coordinate to the x coordinate k
    return probs

def make_grid(n):
    temp = [np.arange(1 , n + 1) for i in range(n)]
    grid = np.meshgrid(*temp)
    return np.reshape(grid, (n, n ** n)).T

def nonzero_probability_counting(n, k, p = 0.5, threshold = 0, sample_size = 100, boundary = "HO-HO"):
    """
    Function that counts the number of preference vectors that are parking functions with nonzero probability given the number of cars and coin flips.
    Parameters are:

    n = number of cars/spots. Spots are counted starting from 1 up to n.
    k = number of coin flips.
    p = Bernoulli parameter.
    threshold = Threshold parameter for counting number of Parking Functions with probability > threshold.
    sample_size = number of times the simulation will run.

    TODO: Add version to count probability > p instead of just > 0.
    """
    grid = make_grid(n)
    c = 0
    for s in grid:
        pk = sample_parking(n, k, p, s, sample_size, boundary)
        if pk > threshold:
            c += 1
        if pk <= threshold:
            print(s)
    return c

#Per's Code
def is_parking_function_recursive(sequence, p, n, attempts_left, k, occupied=None, bc="circle"):
    """
    Function that finds the probability of a given sequence of parking

    sequence = staring preference vector
    p = probability of moving right. Will take a value i.e. 1/2 or a symbolic variable
    n = number of cars and spots
    attempts_left = Put in k, this is used in the recursive function
    k = number of coin flips.
    occupied = list of occupied spots, use [False] * n or leave it empty
    bc = boundary conditions; options are circle, sticky, and leave and don't come back
    """
    if occupied==None:
        occupied=[False] * n
    # Base case: if sequence is empty, all cars have successfully parked
    if len(sequence) == 0:
        return 1

    # Base case: if all spots are occupied, all cars have successfully parked
    if all(occupied):
        return 1  
    
    total_probability = 0  # Initialize the total probability of successful parking for this sequence
    pref = sequence[0] - 1  # Convert the preferred spot to 0-based indexing (the first spot is index 0 etc...)
    
    if occupied[pref]:  # If the desired spot is already occupied
        if attempts_left == 0:
            return 0  # No attempts left, so return 0 as parking failed
        else:
            if bc == "circle":
                # Move left or right in the circle
                if sequence[0] == 1:
                    new_sequence_left = (n,) + sequence[1:]  # Wrap around to the end of the parking lot
                else:
                    new_sequence_left = (sequence[0] - 1,) + sequence[1:]
                if sequence[0] == n:
                    new_sequence_right = (1,) + sequence[1:]
                else:
                    new_sequence_right = (sequence[0] + 1,) + sequence[1:]
                # Recursively calculate the probability of successful parking for both left and right movements
                total_probability += p * is_parking_function_recursive(new_sequence_right, p, n, attempts_left - 1, k, tuple(occupied), bc)
                total_probability += (1 - p) * is_parking_function_recursive(new_sequence_left, p, n, attempts_left - 1, k, tuple(occupied), bc)
            elif bc == "sticky":
                # Move left or right with sticky boundary condition
                if sequence[0] == 1:
                    new_sequence_left = sequence  # Stay in the same spot
                else:
                    new_sequence_left = (sequence[0] - 1,) + sequence[1:]
                if sequence[0] == n:
                    new_sequence_right = sequence  # Stay in the same spot
                else:
                    new_sequence_right = (sequence[0] + 1,) + sequence[1:]
                # Recursively calculate the probability of successful parking for both left and right movements
                total_probability += p * is_parking_function_recursive(new_sequence_right, p, n, attempts_left - 1, k, tuple(occupied), bc)
                total_probability += (1 - p) * is_parking_function_recursive(new_sequence_left, p, n, attempts_left - 1, k, tuple(occupied), bc)
            elif bc == "leave and don't come back":
                # Move left or right with "leave and don't come back" boundary condition
                if sequence[0] == 1:
                    new_sequence_left = None  # Moving left from leftmost spot gives 0 probability
                else:
                    new_sequence_left = (sequence[0] - 1,) + sequence[1:]
                if sequence[0] == n:
                    new_sequence_right = None  # Moving right from rightmost spot gives 0 probability
                else:
                    new_sequence_right = (sequence[0] + 1,) + sequence[1:]
                # Recursively calculate the probability of successful parking for both left and right movements
                if new_sequence_right is not None:
                    total_probability += p * is_parking_function_recursive(new_sequence_right, p, n, attempts_left - 1, k, tuple(occupied), bc)
                if new_sequence_left is not None:
                    total_probability += (1 - p) * is_parking_function_recursive(new_sequence_left, p, n, attempts_left - 1, k, tuple(occupied), bc)
    else:  # If the preferred spot is free
        new_occupied = list(occupied)  # Make a copy of the occupied list to update it
        new_occupied[pref] = True  # Mark the preferred spot as occupied
        remaining_sequence = sequence[1:]  # Remove the first car from the sequence
        
        # Recursively calculate the probability of successful parking for the remaining sequence
        total_probability += is_parking_function_recursive(tuple(remaining_sequence), p, n, k, k, tuple(new_occupied), bc)
        
    return total_probability  # Return the total probability of successful parking for this sequence

