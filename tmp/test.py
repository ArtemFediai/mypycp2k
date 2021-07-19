class Artem:
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    @classmethod
    def from_dict(cls, **kwargs):
        return cls(**kwargs)

    def __repr__(self):
        return_thing = ''
        for this_key, this_value in self.__dict__.items():
            return_thing += this_key + ' -> ' + this_value + '\n'
        return return_thing
        # return f'{self.__dict__}'


artem = Artem(height='180', weight='88')


print(artem)

artem_as_dict = artem.__dict__

artem_2 = Artem.from_dict(**artem_as_dict)

print(artem_2)