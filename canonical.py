from sympy import *
# Reference: "A Gentle Introduction to Optimization" by B.Guenin, J.Könemann and L.Tunçel
# All definitons and ideas are generated from the above book

def canonical_form(A, b, c, k, basis):
    init_printing()

    # Correct the basis index
    basis = [(x-1) for x in basis]

    # Condition 1: A_B is an identity matrix
    # Idea: multiply inverse of A_B on both sides of the equation
    A_B = A[:, basis]
    # Update A and b for the new constraint matrix
    A_ = A_B.inv()*A
    b_ = A_B.inv()*b

    # Condition 2: C_B = 0
    # Idea: 1)  Since y.T*b - y.T*A*x = 0 for arbitrary vector y,
    #           we have z(x):= c.T*x + k = y.T*b - (c.T - y.T*A) * x + k
    #       2)  Note that C_B = (c.T - y.T*A) = 0, then c_B.T = y.T*A_B

    # Calculate the y in step 2)
    y = A_B.T.inv()*c[:,basis].T
    # Update c and k for the new objective function
    c = c - y.T*A
    k = k + (y.T*b)[0]

    # Printing the results
    x = MatrixSymbol('x', c.cols, 1)
    print("\nThe canonical form of the LP is: \n")
    print("maximize")
    pprint(c*x + Matrix([k]))
    print("\nsubject to")
    pprint(A_ * x + b_)

    return A_, b_, c, k

def main():
    print("""Please input the following row by row:
1. number of rows of the constraints matrix A
2. each row of A seperated by space
3. b in a single line seperated by space
4. c in a single line seperated by space
5. basis in a single line seperated by space\n""")
    m = int(input())

    # Input A row by row
    A = []
    for i in range(m):
        row = [S(int(x)) for x in input().split()]
        A.append(row)
    A = Matrix(A)

    # Input vector b
    b = [S(int(x)) for x in input().split()]
    b = Matrix(b)

    # Input vector c
    c = [S(int(x)) for x in input().split()]
    c = Matrix([c])

    # Define constant k
    k = 0

    # Input the new basis
    basis = [int(x) for x in input().split()]

    # Call the function
    canonical_form(A, b, c, k, basis)


if __name__ == "__main__":
    main()