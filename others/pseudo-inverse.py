from sympy import *



mat = [

    [0,          -1,  0        ],
    [-1/sqrt(2),  0,  1/sqrt(2)],
    [-1/sqrt(2),  0, -1/sqrt(2)]

]


M = Matrix(mat)



# I think since this is all symbolic it should be OK to use the text-book formulas taught 
# in a linear algebra class (e.g. see the list of special cases in the Wikipedia article on the 
# Moore–Penrose pseudoinverse). For numerical evaluation pinv uses the singular value decomposition 
# (svd) instead.

# You have linearly independent rows (full row rank), so you can use the formula for a 'right' inverse:
N = M.H * (M * M.H) ** -1
print(N)


# For full column rank, replace M with M.H, transpose the result, and simplify 
# to get the following formula for the 'left' inverse:
N = (M.H * M) ** -1 * M.H
print(N)


#########

print(M.pinv())





