#This is some code that Steve made to work on our problem.

p=var('p')

# Given a parking function this calculates maximum
# amount a car has to move to park
def Lcount(pref):
    n = len(pref)
    spots = set(range(n))
    L = 0
    for i in range(n):
        x = pref[i]
        y = min(set(range(x,n)).intersection(spots))
        spots.discard(y)
        L = max([L,y-x])
    return L

# Over all parking functions of [n] this counts the number
# of parking function where each car moves at most k
def F(n,k):
    count = 0
    PF = [0]*n
    while True:
        for C in Permutations(PF):
            count += Lcount(C) <= k
        for i in range (n):
            if PF[i]<n-1-i:
                PF[i] += 1
                for j in range(i):
                    PF[j] = PF[i]
                break
        else:
            break
    return count

# Sume the probabilities of [n]^n of being a parking function
def probF(n,k):
    totalprob = 0
    F = [0]*n
    while True:
        totalprob += q
        for i in range(n):
            if F[i] == n-1:
                F[i] = 0
            else:
                F[i] += 1
                break
        else:
            break
    if not type(totalprob)==Integer:
        totalprob = totalprob.expand()
    return totalprob

# Given function F computes expression for probability
# that it will work as parking under model. For example
#      probFrecurse([3, 3, 2, 0, 2], set(range(5)), 2)
# returns (in some form)
#      2*p^3 - 5*p^2 + 3*p

def probFrecurse(F,open_spots,k):
    if len(open_spots) == 0:
        return 1
    i = len(F) - len(open_spots) # which car is parking
    if F[i] in open_spots:
        new_open_spots = copy(open_spots)
        new_open_spots.discard(F[i])
        return probFrecurse(F,new_open_spots,k)
    else:
        below = set(range(F[i]-k,F[i]))
        below = below.intersection(open_spots)
        above = set(range(F[i],F[i]+k+1))
        above = above.intersection(open_spots)
        if above and below:
            new_open_spots_aa = copy(open_spots)
            new_open_spots_bb = copy(open_spots)
            new_open_spots_aa.discard(max(below))
            new_open_spots_bb.discard(min(above))
            return (1-p)*probFrecurse(F,new_open_spots_aa,k)+p*probFrecurse(F,new_open_spots_bb,k)
        if above and not below:
            new_open_spots_bb = copy(open_spots)
            new_open_spots_bb.discard(min(above))
            return p*probFrecurse(F,new_open_spots_bb,k)
        if below and not above:
            new_open_spots_aa = copy(open_spots)
            new_open_spots_aa.discard(max(below))
            return (1-p)*probFrecurse(F,new_open_spots_aa,k)
        if not above and not below:
            return 0

for n in [2..8]:
    for k in [0..n]:
        print n, k, ":", F(n,k), probF(n,k)