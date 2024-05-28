import numpy as np
import matplotlib.pyplot as plt
import random


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
