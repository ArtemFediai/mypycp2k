
def my_print(x, silent=False, *args, **kwargs):
    if not silent:
        print(x, *args, **kwargs)
    else:
        pass

my_print('mo', silent=True)