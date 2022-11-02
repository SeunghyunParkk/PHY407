
# Define Hermite polynomial
def H(n,x):
    H = []
    # H_0(x) = 1
    for i in range(0,1):
        H.append(1)
    # H_1(x) = 2*x
    for i in range(1,2):
        H.append(2*x)
    # Calculate H_n+1 using for loop
    for i in range(2,n+1):
        H.append(2*x*H[i-1] - 2*n*H[i-2])
    return H[n]

print(H(3,1))

def HER(x,n):
    if n==0:
        return 1.0 + 0.0*x
    elif n==1:
        return 2.0*x
    else:
        return 2.0*x*HER(x,n-1) -2.0*(n-1)*HER(x,n-2)
print(HER(4,1))


def hermite(n,x):

    if n==0:

        return 1

    elif n==1:

        return 2*x

    else:

        return 2*x*hermite(n-1,x)-2*(n-1)*hermite(n-2,x)



