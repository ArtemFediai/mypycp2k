
def fun(a=0,
        b=0,
        c=0,
        d=0):
    return a + b + c + d



def fun1(a):
    return(fun(a=a,
               b=2,
               c=2,
               d=2))


print(fun(1,1,1,1))


print(fun1(a=2))

print('hello')

a = 1
b = 2
c = 3

def f3(x):
    return a + b + c + x

print(f3(1))