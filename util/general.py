def timeit(fn):
    """
    measures execution time. used as decorator
    """
    from time import perf_counter
    from functools import wraps

    @wraps(fn)
    def inner(*args, **kwargs):
        start = perf_counter()
        result = fn(*args, **kwargs)
        end = perf_counter()
        elapsed = end - start

        args_ = [str(a) for a in args]
        kwargs_ = ['{0}={1}'.format(k, v) for (k, v) in kwargs.items()]
        all_args = args_ + kwargs_
        args_str = ','.join(all_args)
        # print('{0}({1}) took {2:.6f}s to run.'.format(fn.__name__,
        #                                               args_str,
        #                                               elapsed))
        print('{0} took {1:.6f}s to run.'.format(fn.__name__, elapsed))
        return result
    return inner


def parse_mixed_list(input_str):
    """
    '2-5,7,15-17,12' --> [2, 3, 4, 5, 7, 12, 15, 16, 17]
    str --> list
    """
    single_frames = set()

    for element in input_str.split(','):
        parts = [int(x) for x in element.split('-')]
        if len(parts) == 1:
            single_frames.add(parts[0])
        else:
            for frame in range(min(parts), max(parts) + 1):
                single_frames.add(frame)

    return list(single_frames)